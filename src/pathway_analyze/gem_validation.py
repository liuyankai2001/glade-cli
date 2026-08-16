from __future__ import annotations

import json
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.util.solver import linear_reaction_coefficients

from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir
from src.pathway_analyze.target_id import validate_target_compound_id


KEGG_REST_BASE_URL = "https://rest.kegg.jp"
HTTP_TIMEOUT = 30
HTTP_RETRIES = 3
REQUEST_SLEEP_SECONDS = 0.2

DEFAULT_GROWTH_FRACTION = 0.1
DEFAULT_FLUX_THRESHOLD = 1e-8
DEFAULT_FVA_FRACTION = 0.99
DEFAULT_PFBA_FRACTION = 1.0
DEFAULT_PROCESSES = 1
DEFAULT_COMPARTMENT = "c"
DEFAULT_REACTION_UPPER_BOUND = 1000.0

DEFAULT_MODE = "per-solution"
DEFAULT_COFACTOR_MODE = "strict_l1"

KEGG_ID_PATTERN = re.compile(r"^[CDG]\d{5}$")
KEGG_METABOLITE_ANNOTATION_KEYS = ("kegg.compound", "kegg.drug", "kegg.glycan")
REACTION_ANNOTATION_KEY = "kegg.reaction"

GENERIC_COFACTOR_IDS = {
    "C00028",  # Acceptor
    "C00030",  # Reduced acceptor
    "C00125",  # Cytochrome c
    "C00126",  # Reduced cytochrome c
    "C00138",  # Reduced ferredoxin
    "C00139",  # Oxidized ferredoxin
    "C00268",  # Dihydrobiopterin
    "C00272",  # Tetrahydrobiopterin
    "C00342",  # Thioredoxin
    "C00343",  # Oxidized thioredoxin
    "C03024",  # Reduced NADPH---hemoprotein reductase
    "C03161",  # Oxidized NADPH---hemoprotein reductase
    "C14818",  # Reduced ferredoxin [iron-sulfur cluster]
}


def validation_depth_output_dir(base_gap_dir: str | Path, depth: int) -> Path:
    """返回指定 gap depth 对应的 GEM 验证输出目录。"""

    return gap_depth_output_dir(base_gap_dir, depth) / "gem_validation"


@dataclass(frozen=True)
class ReactionRecord:
    reaction_id: str
    name: str
    equation: str
    left_stoichiometry: Tuple[Tuple[str, float], ...]
    right_stoichiometry: Tuple[Tuple[str, float], ...]


@dataclass(frozen=True)
class CompoundRecord:
    compound_id: str
    name: str


@dataclass(frozen=True)
class AddedReaction:
    model_reaction_id: str
    kegg_reaction_id: str
    direction: str
    added_to_model: bool
    source_status: str
    source_solution_ids: Tuple[int, ...]
    source_step_indexes: Tuple[str, ...]
    produced_compound_ids: Tuple[str, ...]
    equation: str
    created_metabolite_ids: Tuple[str, ...]
    generic_compound_ids: Tuple[str, ...]


@dataclass(frozen=True)
class PreparedModel:
    model: cobra.Model
    target_demand_id: str
    target_metabolite_id: str
    biomass_reaction_id: str
    route_reactions: Tuple[AddedReaction, ...]
    opened_generic_compound_ids: Tuple[str, ...]
    created_metabolite_ids: Tuple[str, ...]
    cofactor_mode: str
    issues: Tuple[str, ...]


@dataclass(frozen=True)
class ValidationResult:
    summary_row: Dict[str, Any]
    flux_rows: Tuple[Dict[str, Any], ...]


def safe_mkdir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def normalize_annotation_values(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if item is not None]
    return [str(value)]


def normalize_kegg_id(token: str) -> str:
    token = str(token).strip()
    if ":" in token:
        token = token.split(":", 1)[1]
    return token


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value: Any) -> int:
    return int(float(value))


def sanitize_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")


def parse_kegg_flatfile(text: str) -> Dict[str, List[str]]:
    fields: Dict[str, List[str]] = {}
    current_key: str | None = None
    for raw_line in text.splitlines():
        if raw_line.strip() == "///":
            break
        key = raw_line[:12].strip()
        value = raw_line[12:].rstrip()
        if key:
            current_key = key
            fields.setdefault(key, []).append(value.strip())
        elif current_key is not None:
            fields[current_key].append(value.strip())
    return fields


def first_or_empty(values: Sequence[str]) -> str:
    return values[0] if values else ""


def split_equation(equation: str) -> Tuple[str, str]:
    for arrow in ("<=>", "<=", "=>", "="):
        if arrow in equation:
            left, right = equation.split(arrow, 1)
            return left.strip(), right.strip()
    return equation.strip(), ""


def parse_equation_side(side: str) -> Tuple[Tuple[str, float], ...]:
    rows: List[Tuple[str, float]] = []
    if not side:
        return tuple()
    for raw_part in side.split("+"):
        part = raw_part.strip()
        if not part:
            continue
        tokens = part.split()
        amount = 1.0
        compound_id = ""
        if len(tokens) == 1:
            compound_id = normalize_kegg_id(tokens[0])
        elif re.fullmatch(r"\d+(?:\.\d+)?", tokens[0]):
            amount = float(tokens[0])
            compound_id = normalize_kegg_id(tokens[1])
        else:
            compound_id = normalize_kegg_id(tokens[-1])
        if KEGG_ID_PATTERN.match(compound_id):
            rows.append((compound_id, amount))
    return tuple(rows)


def directional_stoichiometry(
    reaction: ReactionRecord,
    direction: str,
) -> Tuple[Tuple[Tuple[str, float], ...], Tuple[Tuple[str, float], ...]]:
    normalized = str(direction).strip().lower()
    if normalized in {"reverse", "right_to_left", "backward", "rtl"}:
        return reaction.right_stoichiometry, reaction.left_stoichiometry
    return reaction.left_stoichiometry, reaction.right_stoichiometry


class KeggRestClient:
    def __init__(
        self,
        cache_dir: str | Path,
        timeout: int = HTTP_TIMEOUT,
        request_sleep_seconds: float = REQUEST_SLEEP_SECONDS,
    ) -> None:
        self.cache_dir = safe_mkdir(cache_dir).resolve()
        self.timeout = timeout
        self.request_sleep_seconds = max(0.0, float(request_sleep_seconds))
        self._text_cache: Dict[str, str] = {}
        self._compound_name_cache: Dict[str, str] = {}
        self._compound_record_cache: Dict[str, CompoundRecord] = {}
        self._reaction_cache: Dict[str, ReactionRecord] = {}

    def _cache_path(self, namespace: str, key: str, suffix: str = ".txt") -> Path:
        return safe_mkdir(self.cache_dir / namespace) / f"{key}{suffix}"

    def _fetch_text(self, url: str, cache_key: str) -> str:
        if cache_key in self._text_cache:
            return self._text_cache[cache_key]
        namespace, raw_key = cache_key.split(":", 1)
        cache_path = self._cache_path(namespace, raw_key)
        if cache_path.exists():
            text = cache_path.read_text(encoding="utf-8")
            self._text_cache[cache_key] = text
            return text

        last_error: Exception | None = None
        for attempt in range(1, HTTP_RETRIES + 1):
            try:
                request = Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urlopen(request, timeout=self.timeout) as response:
                    text = response.read().decode("utf-8")
                self._text_cache[cache_key] = text
                cache_path.write_text(text, encoding="utf-8")
                if self.request_sleep_seconds > 0:
                    time.sleep(self.request_sleep_seconds)
                return text
            except (HTTPError, URLError, TimeoutError, ConnectionError) as exc:
                last_error = exc
                time.sleep(attempt)
        raise RuntimeError(f"Failed to fetch KEGG URL after retries: {url}") from last_error

    def get_compound_record(self, compound_id: str) -> CompoundRecord:
        compound_id = normalize_kegg_id(compound_id)
        if compound_id in self._compound_record_cache:
            return self._compound_record_cache[compound_id]
        text = self._fetch_text(
            f"{KEGG_REST_BASE_URL}/get/cpd:{compound_id}",
            f"compound:{compound_id}",
        )
        fields = parse_kegg_flatfile(text)
        name = first_or_empty(fields.get("NAME", []))
        if ";" in name:
            name = name.split(";", 1)[0]
        record = CompoundRecord(compound_id=compound_id, name=name)
        self._compound_record_cache[compound_id] = record
        self._compound_name_cache[compound_id] = name
        return record

    def get_compound_name(self, compound_id: str) -> str:
        compound_id = normalize_kegg_id(compound_id)
        if compound_id not in self._compound_name_cache:
            self._compound_name_cache[compound_id] = self.get_compound_record(compound_id).name
        return self._compound_name_cache[compound_id]

    def get_reaction(self, reaction_id: str) -> ReactionRecord:
        reaction_id = normalize_kegg_id(reaction_id)
        if reaction_id in self._reaction_cache:
            return self._reaction_cache[reaction_id]
        text = self._fetch_text(
            f"{KEGG_REST_BASE_URL}/get/rn:{reaction_id}",
            f"reaction:{reaction_id}",
        )
        fields = parse_kegg_flatfile(text)
        equation = " ".join(fields.get("EQUATION", []))
        left_side, right_side = split_equation(equation)
        record = ReactionRecord(
            reaction_id=reaction_id,
            name=first_or_empty(fields.get("NAME", [])),
            equation=equation,
            left_stoichiometry=parse_equation_side(left_side),
            right_stoichiometry=parse_equation_side(right_side),
        )
        self._reaction_cache[reaction_id] = record
        return record


def build_kegg_metabolite_index(model: cobra.Model) -> Dict[str, List[cobra.Metabolite]]:
    index: Dict[str, List[cobra.Metabolite]] = {}
    for metabolite in model.metabolites:
        annotations = getattr(metabolite, "annotation", {}) or {}
        for key in KEGG_METABOLITE_ANNOTATION_KEYS:
            for raw_id in normalize_annotation_values(annotations.get(key)):
                compound_id = normalize_kegg_id(raw_id)
                if KEGG_ID_PATTERN.match(compound_id):
                    index.setdefault(compound_id, []).append(metabolite)
    return index


def build_kegg_reaction_index(model: cobra.Model) -> Dict[str, List[cobra.Reaction]]:
    index: Dict[str, List[cobra.Reaction]] = {}
    for reaction in model.reactions:
        annotations = getattr(reaction, "annotation", {}) or {}
        for raw_id in normalize_annotation_values(annotations.get(REACTION_ANNOTATION_KEY)):
            reaction_id = normalize_kegg_id(raw_id)
            if reaction_id.startswith("R"):
                index.setdefault(reaction_id, []).append(reaction)
    return index


def choose_metabolite(
    candidates: Sequence[cobra.Metabolite],
    preferred_compartment: str,
) -> cobra.Metabolite | None:
    if not candidates:
        return None
    for metabolite in candidates:
        if metabolite.compartment == preferred_compartment:
            return metabolite
    for metabolite in candidates:
        if metabolite.compartment == "c":
            return metabolite
    return candidates[0]


def get_or_create_metabolite(
    model: cobra.Model,
    compound_id: str,
    metabolite_index: Dict[str, List[cobra.Metabolite]],
    client: KeggRestClient,
    preferred_compartment: str,
) -> Tuple[cobra.Metabolite, bool]:
    compound_id = normalize_kegg_id(compound_id)
    existing = choose_metabolite(metabolite_index.get(compound_id, []), preferred_compartment)
    if existing is not None:
        return existing, False

    metabolite_id = sanitize_id(f"kegg_{compound_id}_{preferred_compartment}")
    if metabolite_id in model.metabolites:
        metabolite = model.metabolites.get_by_id(metabolite_id)
        metabolite_index.setdefault(compound_id, []).append(metabolite)
        return metabolite, False

    try:
        name = client.get_compound_name(compound_id)
    except Exception:
        name = compound_id

    metabolite = cobra.Metabolite(
        id=metabolite_id,
        name=name or compound_id,
        compartment=preferred_compartment,
    )
    metabolite.annotation["kegg.compound"] = compound_id
    model.add_metabolites([metabolite])
    metabolite_index.setdefault(compound_id, []).append(metabolite)
    return metabolite, True


def directional_net_stoichiometry(reaction: ReactionRecord, direction: str) -> Dict[str, float]:
    consumed_stoichiometry, produced_stoichiometry = directional_stoichiometry(reaction, direction)
    net: Dict[str, float] = {}
    for compound_id, amount in consumed_stoichiometry:
        net[compound_id] = net.get(compound_id, 0.0) - float(amount)
    for compound_id, amount in produced_stoichiometry:
        net[compound_id] = net.get(compound_id, 0.0) + float(amount)
    return {compound_id: coefficient for compound_id, coefficient in net.items() if abs(coefficient) > 1e-12}


def add_kegg_reaction(
    model: cobra.Model,
    *,
    route_scope: str,
    reaction_id: str,
    direction: str,
    client: KeggRestClient,
    metabolite_index: Dict[str, List[cobra.Metabolite]],
    preferred_compartment: str,
    reaction_upper_bound: float,
) -> Tuple[cobra.Reaction | None, Tuple[str, ...], Tuple[str, ...], str | None]:
    reaction = client.get_reaction(reaction_id)
    net_stoichiometry = directional_net_stoichiometry(reaction, direction)
    if not net_stoichiometry:
        return None, tuple(), tuple(), f"{reaction_id} {direction} has zero net stoichiometry"

    cobra_reaction = cobra.Reaction(
        id=sanitize_id(f"HET_{route_scope}_{reaction_id}_{direction}"),
        name=reaction.name or reaction_id,
        lower_bound=0.0,
        upper_bound=reaction_upper_bound,
    )
    cobra_reaction.annotation["kegg.reaction"] = reaction_id
    cobra_reaction.notes["kegg_direction"] = direction
    cobra_reaction.notes["kegg_equation"] = reaction.equation

    created_metabolites: List[str] = []
    generic_compounds: List[str] = []
    cobra_stoichiometry: Dict[cobra.Metabolite, float] = {}
    for compound_id, coefficient in net_stoichiometry.items():
        metabolite, was_created = get_or_create_metabolite(
            model=model,
            compound_id=compound_id,
            metabolite_index=metabolite_index,
            client=client,
            preferred_compartment=preferred_compartment,
        )
        cobra_stoichiometry[metabolite] = cobra_stoichiometry.get(metabolite, 0.0) + coefficient
        if was_created:
            created_metabolites.append(metabolite.id)
        if compound_id in GENERIC_COFACTOR_IDS:
            generic_compounds.append(compound_id)

    cobra_reaction.add_metabolites(cobra_stoichiometry)
    model.add_reactions([cobra_reaction])
    return cobra_reaction, tuple(sorted(set(created_metabolites))), tuple(sorted(set(generic_compounds))), None


def add_target_demand(
    model: cobra.Model,
    target_compound: str,
    client: KeggRestClient,
    metabolite_index: Dict[str, List[cobra.Metabolite]],
    preferred_compartment: str,
    upper_bound: float,
) -> Tuple[str, str, bool]:
    target_metabolite, was_created = get_or_create_metabolite(
        model=model,
        compound_id=target_compound,
        metabolite_index=metabolite_index,
        client=client,
        preferred_compartment=preferred_compartment,
    )
    demand_id = sanitize_id(f"DM_TARGET_{target_metabolite.id}")
    if demand_id in model.reactions:
        demand = model.reactions.get_by_id(demand_id)
        demand.lower_bound = 0.0
        demand.upper_bound = upper_bound
    else:
        demand = model.add_boundary(
            target_metabolite,
            type="demand",
            reaction_id=demand_id,
            lb=0.0,
            ub=upper_bound,
        )
    return demand.id, target_metabolite.id, was_created


def open_generic_cofactor_sinks(
    model: cobra.Model,
    compound_ids: Iterable[str],
    client: KeggRestClient,
    metabolite_index: Dict[str, List[cobra.Metabolite]],
    preferred_compartment: str,
    bound: float,
) -> Tuple[str, ...]:
    opened: List[str] = []
    for compound_id in sorted(set(compound_ids) & GENERIC_COFACTOR_IDS):
        metabolite, _ = get_or_create_metabolite(
            model=model,
            compound_id=compound_id,
            metabolite_index=metabolite_index,
            client=client,
            preferred_compartment=preferred_compartment,
        )
        sink_id = sanitize_id(f"SK_GENERIC_{metabolite.id}")
        if sink_id in model.reactions:
            sink = model.reactions.get_by_id(sink_id)
            sink.lower_bound = -bound
            sink.upper_bound = bound
        else:
            model.add_boundary(
                metabolite,
                type="sink",
                reaction_id=sink_id,
                lb=-bound,
                ub=bound,
            )
        opened.append(compound_id)
    return tuple(opened)


def get_primary_objective_reaction_id(model: cobra.Model) -> str:
    coefficients = linear_reaction_coefficients(model)
    if not coefficients:
        raise ValueError("Model has no linear reaction objective.")
    reaction, _ = max(coefficients.items(), key=lambda item: abs(item[1]))
    return reaction.id


def apply_growth_floor(model: cobra.Model, biomass_reaction_id: str, required_growth: float) -> None:
    biomass = model.reactions.get_by_id(biomass_reaction_id)
    biomass.lower_bound = max(biomass.lower_bound, required_growth)


def solution_rows_for_scope(all_steps_df: pd.DataFrame, solution_ids: Sequence[int]) -> pd.DataFrame:
    if all_steps_df.empty:
        return all_steps_df.copy()
    return (
        all_steps_df[all_steps_df["solution_id"].astype(int).isin(solution_ids)]
        .copy()
        .sort_values(["solution_id", "step_index"])
    )


def append_reaction_metadata(
    metadata_by_model_id: Dict[str, Dict[str, Any]],
    *,
    model_reaction_id: str,
    kegg_reaction_id: str,
    direction: str,
    added_to_model: bool,
    source_status: str,
    solution_id: int,
    step_index: int,
    produced_compound_id: str,
    equation: str,
    created_metabolite_ids: Sequence[str],
    generic_compound_ids: Sequence[str],
) -> None:
    metadata = metadata_by_model_id.setdefault(
        model_reaction_id,
        {
            "model_reaction_id": model_reaction_id,
            "kegg_reaction_id": kegg_reaction_id,
            "direction": direction,
            "added_to_model": added_to_model,
            "source_status": source_status,
            "source_solution_ids": set(),
            "source_step_indexes": set(),
            "produced_compound_ids": set(),
            "equation": equation,
            "created_metabolite_ids": set(),
            "generic_compound_ids": set(),
        },
    )
    metadata["source_solution_ids"].add(solution_id)
    metadata["source_step_indexes"].add(f"{solution_id}:{step_index}")
    if produced_compound_id:
        metadata["produced_compound_ids"].add(produced_compound_id)
    metadata["created_metabolite_ids"].update(created_metabolite_ids)
    metadata["generic_compound_ids"].update(generic_compound_ids)


def finalize_reaction_metadata(metadata_by_model_id: Dict[str, Dict[str, Any]]) -> Tuple[AddedReaction, ...]:
    rows: List[AddedReaction] = []
    for metadata in metadata_by_model_id.values():
        rows.append(
            AddedReaction(
                model_reaction_id=metadata["model_reaction_id"],
                kegg_reaction_id=metadata["kegg_reaction_id"],
                direction=metadata["direction"],
                added_to_model=bool(metadata["added_to_model"]),
                source_status=metadata["source_status"],
                source_solution_ids=tuple(sorted(metadata["source_solution_ids"])),
                source_step_indexes=tuple(sorted(metadata["source_step_indexes"])),
                produced_compound_ids=tuple(sorted(metadata["produced_compound_ids"])),
                equation=metadata["equation"],
                created_metabolite_ids=tuple(sorted(metadata["created_metabolite_ids"])),
                generic_compound_ids=tuple(sorted(metadata["generic_compound_ids"])),
            )
        )
    return tuple(sorted(rows, key=lambda row: (row.kegg_reaction_id, row.direction, row.model_reaction_id)))


def prepare_model_for_scope(
    base_model: cobra.Model,
    step_rows: pd.DataFrame,
    *,
    target_compound: str,
    client: KeggRestClient,
    route_scope: str,
    biomass_reaction_id: str,
    required_growth: float,
    preferred_compartment: str,
    reaction_upper_bound: float,
    cofactor_mode: str,
) -> PreparedModel:
    model = base_model.copy()
    metabolite_index = build_kegg_metabolite_index(model)
    reaction_index = build_kegg_reaction_index(model)
    metadata_by_model_id: Dict[str, Dict[str, Any]] = {}
    added_reaction_by_key: Dict[Tuple[str, str], cobra.Reaction] = {}
    issues: List[str] = []
    all_created_metabolites: set[str] = set()
    generic_compounds_in_route: set[str] = set()

    for row in step_rows.to_dict("records"):
        solution_id = safe_int(row["solution_id"])
        step_index = safe_int(row["step_index"])
        source_status = str(row.get("status", "heterologous"))
        reaction_id = str(row["reaction_id"])
        direction = str(row["direction"])
        produced_compound_id = str(row.get("produced_compound_id", ""))
        equation = str(row.get("equation", ""))

        if source_status == "endogenous":
            endogenous_reactions = reaction_index.get(reaction_id, [])
            if not endogenous_reactions:
                issues.append(f"{route_scope}: endogenous {reaction_id} was not found in model annotations")
                continue
            for endogenous_reaction in endogenous_reactions:
                append_reaction_metadata(
                    metadata_by_model_id,
                    model_reaction_id=endogenous_reaction.id,
                    kegg_reaction_id=reaction_id,
                    direction=direction,
                    added_to_model=False,
                    source_status=source_status,
                    solution_id=solution_id,
                    step_index=step_index,
                    produced_compound_id=produced_compound_id,
                    equation=equation,
                    created_metabolite_ids=tuple(),
                    generic_compound_ids=tuple(),
                )
            continue

        key = (reaction_id, direction)
        created_metabolites: Tuple[str, ...] = tuple()
        generic_compounds: Tuple[str, ...] = tuple()
        if key in added_reaction_by_key:
            cobra_reaction = added_reaction_by_key[key]
        else:
            cobra_reaction, created_metabolites, generic_compounds, issue = add_kegg_reaction(
                model,
                route_scope=route_scope,
                reaction_id=reaction_id,
                direction=direction,
                client=client,
                metabolite_index=metabolite_index,
                preferred_compartment=preferred_compartment,
                reaction_upper_bound=reaction_upper_bound,
            )
            if issue:
                issues.append(f"{route_scope}: {issue}")
                continue
            if cobra_reaction is None:
                continue
            added_reaction_by_key[key] = cobra_reaction
            all_created_metabolites.update(created_metabolites)
            generic_compounds_in_route.update(generic_compounds)

        append_reaction_metadata(
            metadata_by_model_id,
            model_reaction_id=cobra_reaction.id,
            kegg_reaction_id=reaction_id,
            direction=direction,
            added_to_model=True,
            source_status=source_status,
            solution_id=solution_id,
            step_index=step_index,
            produced_compound_id=produced_compound_id,
            equation=equation,
            created_metabolite_ids=created_metabolites,
            generic_compound_ids=generic_compounds,
        )

    if cofactor_mode == "relaxed":
        opened_generic_compounds = open_generic_cofactor_sinks(
            model=model,
            compound_ids=generic_compounds_in_route,
            client=client,
            metabolite_index=metabolite_index,
            preferred_compartment=preferred_compartment,
            bound=reaction_upper_bound,
        )
    elif cofactor_mode == "strict_l1":
        opened_generic_compounds = tuple()
        if generic_compounds_in_route:
            issues.append(
                f"{route_scope}: strict_l1 cofactor mode left generic carriers closed: "
                + ";".join(sorted(generic_compounds_in_route))
            )
    else:
        raise ValueError("cofactor_mode must be relaxed or strict_l1")

    target_demand_id, target_metabolite_id, target_created = add_target_demand(
        model=model,
        target_compound=target_compound,
        client=client,
        metabolite_index=metabolite_index,
        preferred_compartment=preferred_compartment,
        upper_bound=reaction_upper_bound,
    )
    if target_created:
        all_created_metabolites.add(target_metabolite_id)

    apply_growth_floor(model, biomass_reaction_id, required_growth)
    model.objective = model.reactions.get_by_id(target_demand_id)
    model.objective_direction = "max"

    return PreparedModel(
        model=model,
        target_demand_id=target_demand_id,
        target_metabolite_id=target_metabolite_id,
        biomass_reaction_id=biomass_reaction_id,
        route_reactions=finalize_reaction_metadata(metadata_by_model_id),
        opened_generic_compound_ids=opened_generic_compounds,
        created_metabolite_ids=tuple(sorted(all_created_metabolites)),
        cofactor_mode=cofactor_mode,
        issues=tuple(issues),
    )


def flux_value(fluxes: pd.Series | None, reaction_id: str) -> float:
    if fluxes is None or reaction_id not in fluxes.index:
        return float("nan")
    return float(fluxes.loc[reaction_id])


def fva_bounds(fva_df: pd.DataFrame | None, reaction_id: str) -> Tuple[float, float]:
    if fva_df is None or reaction_id not in fva_df.index:
        return float("nan"), float("nan")
    return float(fva_df.loc[reaction_id, "minimum"]), float(fva_df.loc[reaction_id, "maximum"])


def can_carry_flux(minimum: float, maximum: float, threshold: float) -> bool:
    return max(abs(safe_float(minimum, 0.0)), abs(safe_float(maximum, 0.0))) > threshold


def is_required_flux(minimum: float, maximum: float, threshold: float) -> bool:
    return minimum > threshold and maximum > threshold


def classify_validation(
    *,
    mode: str,
    fba_product_flux: float,
    route_reactions: Sequence[AddedReaction],
    fva_df: pd.DataFrame | None,
    flux_threshold: float,
) -> str:
    if not math.isfinite(fba_product_flux) or fba_product_flux <= flux_threshold:
        return "FAIL_NO_PRODUCT_FLUX"
    if mode == "pooled":
        return "PASS_POOLED_PRODUCT_FLUX"
    if fva_df is None:
        return "PASS_PRODUCT_FLUX_FVA_NOT_RUN"

    blocked = []
    for reaction in route_reactions:
        if not reaction.added_to_model:
            continue
        minimum, maximum = fva_bounds(fva_df, reaction.model_reaction_id)
        if not can_carry_flux(minimum, maximum, flux_threshold):
            blocked.append(reaction.kegg_reaction_id)
    if blocked:
        return "WARN_PRODUCT_BUT_ROUTE_STEP_BLOCKED"

    optional = []
    for reaction in route_reactions:
        if not reaction.added_to_model:
            continue
        minimum, maximum = fva_bounds(fva_df, reaction.model_reaction_id)
        if not is_required_flux(minimum, maximum, flux_threshold):
            optional.append(reaction.kegg_reaction_id)
    if optional:
        return "PASS_ROUTE_CAN_CARRY_FLUX_OPTIONAL"
    return "PASS_ROUTE_REQUIRED_AT_TARGET_FLUX"


def validate_prepared_model(
    prepared: PreparedModel,
    *,
    route_scope: str,
    validation_mode: str,
    solution_ids: Sequence[int],
    baseline_growth: float,
    required_growth: float,
    fva_fraction: float,
    pfba_fraction: float,
    flux_threshold: float,
    run_fva: bool,
    processes: int,
) -> ValidationResult:
    fba_solution = prepared.model.optimize()
    fba_status = fba_solution.status
    fba_product_flux = safe_float(fba_solution.objective_value, float("nan"))
    fba_fluxes = fba_solution.fluxes if fba_status == "optimal" else None

    pfba_status = "not_run"
    pfba_product_flux = float("nan")
    pfba_growth_flux = float("nan")
    pfba_total_abs_flux = float("nan")
    pfba_fluxes: pd.Series | None = None
    if fba_status == "optimal" and fba_product_flux > flux_threshold:
        try:
            pfba_solution = pfba(prepared.model, fraction_of_optimum=pfba_fraction)
            pfba_status = pfba_solution.status
            if pfba_status == "optimal":
                pfba_fluxes = pfba_solution.fluxes
                pfba_product_flux = flux_value(pfba_fluxes, prepared.target_demand_id)
                pfba_growth_flux = flux_value(pfba_fluxes, prepared.biomass_reaction_id)
                pfba_total_abs_flux = float(pfba_fluxes.abs().sum())
        except Exception as exc:
            pfba_status = f"error: {exc}"

    fva_df: pd.DataFrame | None = None
    fva_status = "not_run"
    fva_reaction_ids = [
        prepared.target_demand_id,
        prepared.biomass_reaction_id,
        *[reaction.model_reaction_id for reaction in prepared.route_reactions if reaction.model_reaction_id],
    ]
    fva_reaction_ids = list(dict.fromkeys(fva_reaction_ids))
    if run_fva and fba_status == "optimal" and fba_product_flux > flux_threshold:
        try:
            fva_df = flux_variability_analysis(
                prepared.model,
                reaction_list=fva_reaction_ids,
                fraction_of_optimum=fva_fraction,
                processes=processes,
            )
            fva_status = "optimal"
        except Exception as exc:
            fva_status = f"error: {exc}"

    flux_rows: List[Dict[str, Any]] = []
    active_route_reactions: List[str] = []
    pfba_active_route_reactions: List[str] = []
    fva_capable_route_reactions: List[str] = []
    required_route_reactions: List[str] = []
    blocked_route_reactions: List[str] = []

    for reaction in prepared.route_reactions:
        fva_min, fva_max = fva_bounds(fva_df, reaction.model_reaction_id)
        fba_flux = flux_value(fba_fluxes, reaction.model_reaction_id)
        pfba_flux = flux_value(pfba_fluxes, reaction.model_reaction_id)
        can_carry = can_carry_flux(fva_min, fva_max, flux_threshold) if fva_df is not None else False
        required = is_required_flux(fva_min, fva_max, flux_threshold) if fva_df is not None else False
        active_pfba = math.isfinite(pfba_flux) and abs(pfba_flux) > flux_threshold

        if active_pfba or can_carry:
            active_route_reactions.append(reaction.kegg_reaction_id)
        if active_pfba:
            pfba_active_route_reactions.append(reaction.kegg_reaction_id)
        if can_carry:
            fva_capable_route_reactions.append(reaction.kegg_reaction_id)
        if required:
            required_route_reactions.append(reaction.kegg_reaction_id)
        if fva_df is not None and reaction.added_to_model and not can_carry:
            blocked_route_reactions.append(reaction.kegg_reaction_id)

        flux_rows.append(
            {
                "route_scope": route_scope,
                "validation_mode": validation_mode,
                "source_solution_ids": ";".join(str(item) for item in reaction.source_solution_ids),
                "source_step_indexes": ";".join(reaction.source_step_indexes),
                "kegg_reaction_id": reaction.kegg_reaction_id,
                "direction": reaction.direction,
                "model_reaction_id": reaction.model_reaction_id,
                "added_to_model": reaction.added_to_model,
                "source_status": reaction.source_status,
                "produced_compound_ids": ";".join(reaction.produced_compound_ids),
                "equation": reaction.equation,
                "fba_flux": fba_flux,
                "pfba_flux": pfba_flux,
                "fva_minimum": fva_min,
                "fva_maximum": fva_max,
                "fva_can_carry_flux": can_carry,
                "fva_required_at_target_flux": required,
                "pfba_active": active_pfba,
                "created_metabolite_ids": ";".join(reaction.created_metabolite_ids),
                "generic_compound_ids": ";".join(reaction.generic_compound_ids),
            }
        )

    target_fva_min, target_fva_max = fva_bounds(fva_df, prepared.target_demand_id)
    biomass_fba_flux = flux_value(fba_fluxes, prepared.biomass_reaction_id)
    validation_status = classify_validation(
        mode=validation_mode,
        fba_product_flux=fba_product_flux,
        route_reactions=prepared.route_reactions,
        fva_df=fva_df,
        flux_threshold=flux_threshold,
    )

    summary_row = {
        "route_scope": route_scope,
        "validation_mode": validation_mode,
        "solution_ids": ";".join(str(item) for item in solution_ids),
        "validation_status": validation_status,
        "baseline_growth": baseline_growth,
        "required_growth": required_growth,
        "target_demand_reaction_id": prepared.target_demand_id,
        "target_metabolite_id": prepared.target_metabolite_id,
        "biomass_reaction_id": prepared.biomass_reaction_id,
        "fba_status": fba_status,
        "fba_product_flux": fba_product_flux,
        "fba_growth_flux": biomass_fba_flux,
        "pfba_status": pfba_status,
        "pfba_product_flux": pfba_product_flux,
        "pfba_growth_flux": pfba_growth_flux,
        "pfba_total_abs_flux": pfba_total_abs_flux,
        "fva_status": fva_status,
        "fva_fraction_of_product_optimum": fva_fraction,
        "target_fva_minimum": target_fva_min,
        "target_fva_maximum": target_fva_max,
        "route_reaction_count": len(prepared.route_reactions),
        "active_route_reaction_count": len(set(active_route_reactions)),
        "pfba_active_route_reaction_count": len(set(pfba_active_route_reactions)),
        "fva_capable_route_reaction_count": len(set(fva_capable_route_reactions)),
        "required_route_reaction_count": len(set(required_route_reactions)),
        "blocked_route_reaction_count": len(set(blocked_route_reactions)),
        "active_route_reaction_ids": ";".join(sorted(set(active_route_reactions))),
        "pfba_active_route_reaction_ids": ";".join(sorted(set(pfba_active_route_reactions))),
        "fva_capable_route_reaction_ids": ";".join(sorted(set(fva_capable_route_reactions))),
        "required_route_reaction_ids": ";".join(sorted(set(required_route_reactions))),
        "blocked_route_reaction_ids": ";".join(sorted(set(blocked_route_reactions))),
        "created_metabolite_ids": ";".join(prepared.created_metabolite_ids),
        "cofactor_mode": prepared.cofactor_mode,
        "cofactor_relaxed": prepared.cofactor_mode == "relaxed",
        "opened_generic_compound_ids": ";".join(prepared.opened_generic_compound_ids),
        "issues": " | ".join(prepared.issues),
    }
    return ValidationResult(summary_row=summary_row, flux_rows=tuple(flux_rows))


def parse_solution_ids(
    raw_value: Any,
    available_solution_ids: Sequence[int],
) -> Tuple[int, ...]:
    if raw_value is None or raw_value == "":
        return tuple(sorted(set(available_solution_ids)))

    if isinstance(raw_value, str):
        raw_parts: Iterable[Any] = re.split(r"[,;\s]+", raw_value.strip())
    elif isinstance(raw_value, int):
        raw_parts = (raw_value,)
    else:
        raw_parts = raw_value

    selected = list(dict.fromkeys(
        int(part)
        for part in raw_parts
        if str(part).strip()
    ))
    available = set(available_solution_ids)
    missing = sorted(set(selected) - available)
    if missing:
        raise ValueError(f"Requested solution IDs not found in gap output: {missing}")
    return tuple(selected)


def load_medium(medium_path: str | Path) -> Dict[str, float]:
    with Path(medium_path).open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    return {str(key): float(value) for key, value in raw.items()}


def validate_gap_output(
    *,
    target_compound: str,
    gap_dir: Path,
    model_path: Path,
    medium_path: Path,
    cache_dir: Path,
    output_dir: Path,
    mode: str,
    solution_ids: Sequence[int] | None = None,
    growth_fraction: float = DEFAULT_GROWTH_FRACTION,
    flux_threshold: float = DEFAULT_FLUX_THRESHOLD,
    fva_fraction: float = DEFAULT_FVA_FRACTION,
    pfba_fraction: float = DEFAULT_PFBA_FRACTION,
    run_fva: bool = True,
    processes: int = DEFAULT_PROCESSES,
    preferred_compartment: str = DEFAULT_COMPARTMENT,
    reaction_upper_bound: float = DEFAULT_REACTION_UPPER_BOUND,
    cofactor_mode: str = DEFAULT_COFACTOR_MODE,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if mode not in {"per-solution", "pooled", "both"}:
        raise ValueError("mode must be per-solution, pooled, or both")
    if cofactor_mode == "strict":
        cofactor_mode = "strict_l1"
    if cofactor_mode not in {"relaxed", "strict_l1"}:
        raise ValueError("cofactor_mode must be relaxed or strict_l1")

    steps_path = gap_dir / "all_solution_steps.csv"
    if not steps_path.exists():
        raise FileNotFoundError(f"Missing gap step table: {steps_path}")

    all_steps_df = pd.read_csv(steps_path)
    if all_steps_df.empty:
        raise ValueError(f"No solution steps found in {steps_path}")

    available_solution_ids = tuple(sorted(int(value) for value in all_steps_df["solution_id"].unique()))
    selected_solution_ids = tuple(solution_ids or available_solution_ids)

    client = KeggRestClient(cache_dir)
    base_model = cobra.io.load_json_model(str(model_path))
    base_model.medium = load_medium(medium_path)

    baseline_growth = safe_float(base_model.slim_optimize(), float("nan"))
    if not math.isfinite(baseline_growth) or baseline_growth <= flux_threshold:
        raise ValueError("Base model has no meaningful growth under the configured medium.")
    required_growth = growth_fraction * baseline_growth
    biomass_reaction_id = get_primary_objective_reaction_id(base_model)

    validation_modes = ("per-solution", "pooled") if mode == "both" else (mode,)
    summary_rows: List[Dict[str, Any]] = []
    flux_rows: List[Dict[str, Any]] = []

    for validation_mode in validation_modes:
        scopes: List[Tuple[str, Tuple[int, ...], pd.DataFrame]] = []
        if validation_mode == "per-solution":
            for solution_id in selected_solution_ids:
                rows = solution_rows_for_scope(all_steps_df, [solution_id])
                scopes.append((f"solution_{solution_id}", (solution_id,), rows))
        elif validation_mode == "pooled":
            rows = solution_rows_for_scope(all_steps_df, selected_solution_ids)
            scope_name = "pooled_" + "_".join(str(item) for item in selected_solution_ids)
            scopes.append((scope_name, tuple(selected_solution_ids), rows))

        for route_scope, scope_solution_ids, step_rows in scopes:
            prepared = prepare_model_for_scope(
                base_model,
                step_rows,
                target_compound=target_compound,
                client=client,
                route_scope=route_scope,
                biomass_reaction_id=biomass_reaction_id,
                required_growth=required_growth,
                preferred_compartment=preferred_compartment,
                reaction_upper_bound=reaction_upper_bound,
                cofactor_mode=cofactor_mode,
            )
            result = validate_prepared_model(
                prepared,
                route_scope=route_scope,
                validation_mode=validation_mode,
                solution_ids=scope_solution_ids,
                baseline_growth=baseline_growth,
                required_growth=required_growth,
                fva_fraction=fva_fraction,
                pfba_fraction=pfba_fraction,
                flux_threshold=flux_threshold,
                run_fva=run_fva,
                processes=processes,
            )
            summary_rows.append(result.summary_row)
            flux_rows.extend(result.flux_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_df = pd.DataFrame(summary_rows)
    flux_df = pd.DataFrame(flux_rows)
    summary_df.to_csv(output_dir / "gem_validation_summary.csv", index=False, encoding="utf-8-sig")
    flux_df.to_csv(output_dir / "gem_validation_route_fluxes.csv", index=False, encoding="utf-8-sig")
    return summary_df, flux_df


def gem_validate(config: Any) -> dict[str, Any]:
    """使用 ``RunConfig`` 验证 gap 候选路径的 GEM 通量可行性。"""

    target_compound = validate_target_compound_id(config.target_name)
    model_path = Path(config.model_path).expanduser().resolve()
    medium_path = Path(config.medium_path).expanduser().resolve()
    gap_root_dir = Path(config.gap_output_path).expanduser().resolve()
    expansion_depth = int(getattr(config, "depth", 0))
    gap_dir = gap_depth_output_dir(gap_root_dir, expansion_depth)
    cache_dir = Path(config.cache_dir).expanduser().resolve() / "kegg"
    output_dir = validation_depth_output_dir(gap_root_dir, expansion_depth)

    requested_mode = str(
        getattr(config, "validation_mode", DEFAULT_MODE)
    ).strip().lower()
    mode = "per-solution" if requested_mode == "per" else requested_mode
    solution_ids = getattr(
        config,
        "solutions",
        getattr(config, "validation_solution_ids", None),
    )
    cofactor_mode = str(
        getattr(config, "validation_cofactor_mode", DEFAULT_COFACTOR_MODE)
    ).strip().lower()
    skip_fva = bool(getattr(config, "validation_skip_fva", False))

    if mode not in {"per-solution", "pooled", "both"}:
        raise ValueError("validation_mode must be per, pooled, or both")
    if cofactor_mode == "strict":
        cofactor_mode = "strict_l1"
    if cofactor_mode not in {"strict_l1", "relaxed"}:
        raise ValueError("validation_cofactor_mode must be strict_l1 or relaxed")
    if not model_path.is_file():
        raise FileNotFoundError(f"GEM model file not found: {model_path}")
    if model_path.suffix.lower() != ".json":
        raise ValueError(f"Only JSON GEM models are supported: {model_path}")
    if not medium_path.is_file():
        raise FileNotFoundError(f"Medium file not found: {medium_path}")

    steps_path = gap_dir / "all_solution_steps.csv"
    if not steps_path.is_file():
        raise FileNotFoundError(
            "Gap analysis output not found. Run the gap command first: "
            f"{steps_path}"
        )
    try:
        all_steps_df = pd.read_csv(steps_path)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"Gap analysis produced no candidate solutions: {steps_path}") from exc
    if all_steps_df.empty:
        raise ValueError(f"Gap analysis produced no candidate solutions: {steps_path}")
    if "solution_id" not in all_steps_df.columns:
        raise ValueError(f"Missing solution_id column in gap output: {steps_path}")

    available_solution_ids = tuple(
        sorted(int(value) for value in all_steps_df["solution_id"].unique())
    )
    selected_solution_ids = parse_solution_ids(
        solution_ids,
        available_solution_ids,
    )

    print(
        f"[INFO] validating {len(selected_solution_ids)} gap solution(s) "
        f"for {target_compound} at depth {expansion_depth}"
    )
    print(f"[INFO] gap solution directory: {gap_dir}")
    print(f"[INFO] validation output directory: {output_dir}")
    summary_df, flux_df = validate_gap_output(
        target_compound=target_compound,
        gap_dir=gap_dir,
        model_path=model_path,
        medium_path=medium_path,
        cache_dir=cache_dir,
        output_dir=output_dir,
        mode=mode,
        solution_ids=selected_solution_ids,
        growth_fraction=DEFAULT_GROWTH_FRACTION,
        flux_threshold=DEFAULT_FLUX_THRESHOLD,
        fva_fraction=DEFAULT_FVA_FRACTION,
        pfba_fraction=DEFAULT_PFBA_FRACTION,
        run_fva=not skip_fva,
        processes=DEFAULT_PROCESSES,
        preferred_compartment=DEFAULT_COMPARTMENT,
        reaction_upper_bound=DEFAULT_REACTION_UPPER_BOUND,
        cofactor_mode=cofactor_mode,
    )

    summary_csv = output_dir / "gem_validation_summary.csv"
    fluxes_csv = output_dir / "gem_validation_route_fluxes.csv"
    statuses = (
        [str(value) for value in summary_df["validation_status"].tolist()]
        if not summary_df.empty and "validation_status" in summary_df.columns
        else []
    )
    return {
        "ok": True,
        "target_compound": target_compound,
        "expansion_depth": expansion_depth,
        "gap_dir": str(gap_dir.resolve()),
        "validation_mode": mode,
        "solution_ids": list(selected_solution_ids),
        "cofactor_mode": cofactor_mode,
        "fva_skipped": skip_fva,
        "validation_dir": str(output_dir.resolve()),
        "validation_summary_csv": str(summary_csv.resolve()),
        "validation_route_fluxes_csv": str(fluxes_csv.resolve()),
        "validation_row_count": int(len(summary_df)),
        "route_flux_row_count": int(len(flux_df)),
        "validation_statuses": statuses,
    }


def run_validation(config: Any) -> dict[str, Any]:
    """命令行入口；JSON 读取和 ``RunConfig`` 构造由 ``main.py`` 负责。"""

    result = gem_validate(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result
