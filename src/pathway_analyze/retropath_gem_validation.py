"""Strict GEM validation for fully reconstructed RetroPath candidate DAGs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba

from src.pathway_analyze.gem_validation import (
    DEFAULT_COMPARTMENT,
    DEFAULT_FLUX_THRESHOLD,
    DEFAULT_FVA_FRACTION,
    DEFAULT_GROWTH_FRACTION,
    DEFAULT_REACTION_UPPER_BOUND,
    KeggRestClient,
    add_kegg_reaction,
    apply_growth_floor,
    build_kegg_metabolite_index,
    build_kegg_reaction_index,
    get_primary_objective_reaction_id,
    load_medium,
    normalize_annotation_values,
    safe_float,
    sanitize_id,
)
from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEP_COLUMNS,
    CANDIDATE_STEPS_FILE_NAME,
    REJECTED_ROUTE_COLUMNS,
    REJECTED_ROUTES_FILE_NAME,
)
from src.pathway_analyze.retropath_mnxref import MnxrefIndex
from src.pathway_analyze.retropath_pipeline import (
    PIPELINE_RESULT_FILE_NAME,
    RETROPATH_PIPELINE_SCHEMA,
)
from src.pathway_analyze.retropath_stoichiometry import (
    CompletedReactionHypothesis,
    CompoundProperty,
    ReconstructionRejection,
    StepReconstructionResult,
    enumerate_candidate_hypotheses,
    load_p2_compound_properties,
    parse_formula,
    reconstruct_retropath_step,
)
from src.pathway_analyze.target_id import validate_target_compound_id


RETROPATH_GEM_VALIDATION_SCHEMA = "retropath_gem_validation.v1"
VALIDATION_MANIFEST_FILE_NAME = "validation_manifest.json"
STOICHIOMETRY_HYPOTHESES_FILE_NAME = "stoichiometry_hypotheses.csv"
STOICHIOMETRY_TERMS_FILE_NAME = "stoichiometry_terms.csv"
VALIDATION_SUMMARY_FILE_NAME = "gem_validation_summary.csv"
VALIDATION_FLUX_FILE_NAME = "gem_validation_route_fluxes.csv"
REJECTED_HYPOTHESES_FILE_NAME = "rejected_hypotheses.csv"

DEFAULT_MAX_STEP_HYPOTHESES = 8
DEFAULT_MAX_CANDIDATE_COMBINATIONS = 32
# Keep the forced-route flux safely above common LP feasibility tolerances.
DEFAULT_ROUTE_MIN_FLUX = 1e-4

HYPOTHESIS_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "step_id",
    "hypothesis_id",
    "rule_id",
    "source_mnxr_id",
    "source_reference",
    "source_equation",
    "source_orientation",
    "evidence_grade",
    "balance_status",
    "cofactor_reconstruction_status",
    "recovered_compound_ids",
)
TERM_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "step_id",
    "hypothesis_id",
    "side",
    "coefficient",
    "role",
    "compound_id",
    "source_mnxm_id",
    "name",
    "formula",
    "charge",
    "inchi",
    "inchikey",
    "smiles",
    "xrefs",
)
SUMMARY_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "combination_id",
    "candidate_status",
    "validation_status",
    "stoichiometry_hypothesis_ids",
    "baseline_growth",
    "required_growth",
    "fba_status",
    "target_flux",
    "growth_flux",
    "pfba_status",
    "pfba_target_flux",
    "fva_status",
    "route_step_count",
    "active_route_step_count",
    "blocked_route_step_ids",
    "created_metabolite_ids",
    "issues",
    "combination_truncated",
)
FLUX_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "combination_id",
    "step_index",
    "step_id",
    "step_source",
    "model_reaction_id",
    "direction_sign",
    "fba_flux",
    "directed_fba_flux",
    "pfba_flux",
    "directed_pfba_flux",
    "fva_minimum",
    "fva_maximum",
    "required_minimum_flux",
    "active",
)
REJECTION_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "step_id",
    "rule_id",
    "source_mnxr_id",
    "reason_code",
    "reason_detail",
)


@dataclass(frozen=True)
class _CandidateInputs:
    target_compound: str
    depth: int
    retropath_dir: Path
    pipeline_path: Path
    pipeline: Mapping[str, Any]
    routes: tuple[Mapping[str, str], ...]
    steps: tuple[Mapping[str, str], ...]
    rejected_routes: tuple[Mapping[str, str], ...]
    compound_mapping_path: Path


@dataclass(frozen=True)
class _RouteReaction:
    step_index: int
    step_id: str
    step_source: str
    model_reaction_id: str
    direction_sign: float


@dataclass(frozen=True)
class _ModelValidation:
    row: Mapping[str, Any]
    flux_rows: tuple[Mapping[str, Any], ...]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _atomic_write_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(text)
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> str:
    text = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        indent=2,
        allow_nan=False,
    ) + "\n"
    return _atomic_write_text(path, text)


def _render_csv(columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> str:
    import io

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=tuple(columns),
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def _as_int(value: Any, field_name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be an integer")
    try:
        normalized = int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer") from exc
    if normalized < minimum:
        raise ValueError(f"{field_name} must be >= {minimum}")
    return normalized


def _read_csv(path: Path, columns: Sequence[str]) -> tuple[Mapping[str, str], ...]:
    if not path.is_file():
        raise FileNotFoundError(f"RetroPath candidate artifact not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise ValueError(f"RetroPath candidate artifact schema mismatch: {path}")
        return tuple(dict(row) for row in reader)


def _verified_artifact(
    retropath_dir: Path,
    artifacts: Mapping[str, Any],
    key: str,
    file_name: str,
) -> Path:
    record = artifacts.get(key)
    if not isinstance(record, Mapping):
        raise ValueError(f"pipeline artifacts.{key} is missing")
    path = (retropath_dir / file_name).resolve()
    if Path(str(record.get("path") or "")).expanduser().resolve() != path:
        raise ValueError(f"pipeline artifacts.{key}.path is stale")
    if not path.is_file() or _sha256_file(path) != record.get("sha256"):
        raise ValueError(f"pipeline artifacts.{key} checksum mismatch")
    return path


def _load_candidate_inputs(config: Any) -> _CandidateInputs:
    target = validate_target_compound_id(config.target_name)
    depth = _as_int(getattr(config, "depth", 0), "depth")
    retropath_dir = (
        gap_depth_output_dir(Path(config.gap_output_path), depth) / "retropath"
    ).expanduser().resolve()
    pipeline_path = retropath_dir / PIPELINE_RESULT_FILE_NAME
    if not pipeline_path.is_file():
        raise FileNotFoundError(
            f"RetroPath pipeline result not found: {pipeline_path}"
        )
    try:
        pipeline = json.loads(pipeline_path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid RetroPath pipeline result: {pipeline_path}") from exc
    if (
        not isinstance(pipeline, dict)
        or pipeline.get("schema_version") != RETROPATH_PIPELINE_SCHEMA
        or pipeline.get("ok") is not True
        or pipeline.get("status") != "retropath_candidates_found"
    ):
        raise ValueError("RetroPath pipeline did not produce candidate routes")
    if pipeline.get("target_compound") != target:
        raise ValueError("RetroPath pipeline target does not match current input")
    if _as_int(pipeline.get("expansion_depth"), "expansion_depth") != depth:
        raise ValueError("RetroPath pipeline depth does not match the requested depth")
    artifacts = pipeline.get("artifacts")
    if not isinstance(artifacts, Mapping):
        raise ValueError("RetroPath pipeline artifacts are missing")
    route_path = _verified_artifact(
        retropath_dir,
        artifacts,
        "candidate_routes",
        CANDIDATE_ROUTES_FILE_NAME,
    )
    step_path = _verified_artifact(
        retropath_dir,
        artifacts,
        "candidate_steps",
        CANDIDATE_STEPS_FILE_NAME,
    )
    rejected_path = _verified_artifact(
        retropath_dir,
        artifacts,
        "rejected_routes",
        REJECTED_ROUTES_FILE_NAME,
    )
    mapping_record = artifacts.get("compound_mapping")
    if not isinstance(mapping_record, Mapping):
        raise ValueError("pipeline compound_mapping artifact is missing")
    mapping_path = Path(str(mapping_record.get("path") or "")).expanduser().resolve()
    expected_mapping = (retropath_dir / "input" / "compound_mapping.csv").resolve()
    if mapping_path != expected_mapping:
        raise ValueError("pipeline compound_mapping path is stale")
    if not mapping_path.is_file() or _sha256_file(mapping_path) != mapping_record.get(
        "sha256"
    ):
        raise ValueError("pipeline compound_mapping checksum mismatch")
    routes = _read_csv(route_path, CANDIDATE_ROUTE_COLUMNS)
    steps = _read_csv(step_path, CANDIDATE_STEP_COLUMNS)
    rejected = _read_csv(rejected_path, REJECTED_ROUTE_COLUMNS)
    if len(routes) != _as_int(pipeline.get("candidate_count"), "candidate_count"):
        raise ValueError("candidate count does not match pipeline result")
    return _CandidateInputs(
        target_compound=target,
        depth=depth,
        retropath_dir=retropath_dir,
        pipeline_path=pipeline_path,
        pipeline=pipeline,
        routes=routes,
        steps=steps,
        rejected_routes=rejected,
        compound_mapping_path=mapping_path,
    )


def _selected_candidate_ranks(
    raw: Any,
    available: Sequence[int],
) -> tuple[int, ...]:
    if raw is None:
        raise ValueError("--retropath-candidates was not requested")
    if raw == []:
        return tuple(sorted(available))
    selected = tuple(dict.fromkeys(_as_int(value, "candidate rank", minimum=1) for value in raw))
    missing = sorted(set(selected) - set(available))
    if missing:
        raise ValueError(f"RetroPath candidate ranks not found: {missing}")
    return selected


def _annotation_values(value: Any) -> set[str]:
    return {
        str(item).strip().upper()
        for item in normalize_annotation_values(value)
        if str(item).strip()
    }


def _formula_matches(metabolite: cobra.Metabolite, prop: CompoundProperty) -> bool:
    try:
        model_formula = parse_formula(str(metabolite.formula or ""))
        expected_formula = parse_formula(prop.formula)
    except ValueError:
        return False
    charge = metabolite.charge
    return model_formula == expected_formula and charge is not None and int(charge) == prop.charge


class _MetaboliteResolver:
    def __init__(self, model: cobra.Model, compartment: str) -> None:
        self.model = model
        self.compartment = compartment
        self.by_mnxm: dict[str, list[cobra.Metabolite]] = defaultdict(list)
        self.by_inchikey: dict[str, list[cobra.Metabolite]] = defaultdict(list)
        self.by_kegg: dict[str, list[cobra.Metabolite]] = defaultdict(list)
        self.by_chebi: dict[str, list[cobra.Metabolite]] = defaultdict(list)
        self.created: set[str] = set()
        for metabolite in model.metabolites:
            annotation = metabolite.annotation or {}
            for value in _annotation_values(annotation.get("metanetx.chemical")):
                self.by_mnxm[value].append(metabolite)
            for value in _annotation_values(annotation.get("inchi_key")):
                self.by_inchikey[value].append(metabolite)
            for value in _annotation_values(annotation.get("kegg.compound")):
                self.by_kegg[value.removeprefix("CPD:")].append(metabolite)
            for value in _annotation_values(annotation.get("chebi")):
                self.by_chebi[value.removeprefix("CHEBI:")].append(metabolite)

    def _choose(
        self,
        values: Iterable[cobra.Metabolite],
        prop: CompoundProperty,
    ) -> cobra.Metabolite | None:
        unique = {item.id: item for item in values}.values()
        compatible = [item for item in unique if _formula_matches(item, prop)]
        compatible.sort(
            key=lambda item: (
                0 if item.compartment == self.compartment else 1,
                item.id,
            )
        )
        return compatible[0] if compatible else None

    def resolve(self, prop: CompoundProperty) -> cobra.Metabolite:
        candidates: list[cobra.Metabolite] = []
        if prop.source_mnxm_id:
            candidates.extend(self.by_mnxm.get(prop.source_mnxm_id.upper(), []))
        if prop.inchikey:
            candidates.extend(self.by_inchikey.get(prop.inchikey.upper(), []))
        if re.fullmatch(r"C\d{5}", prop.compound_id):
            candidates.extend(self.by_kegg.get(prop.compound_id.upper(), []))
        for xref in prop.xrefs:
            prefix, separator, identifier = str(xref).partition(":")
            if not separator:
                continue
            if prefix.lower() == "kegg":
                candidates.extend(self.by_kegg.get(identifier.upper(), []))
            elif prefix.lower() == "chebi":
                candidates.extend(self.by_chebi.get(identifier.upper(), []))
        selected = self._choose(candidates, prop)
        if selected is not None:
            return selected
        metabolite_id = sanitize_id(
            f"p8_{prop.compound_id}_{self.compartment}"
        )
        if metabolite_id in self.model.metabolites:
            existing = self.model.metabolites.get_by_id(metabolite_id)
            if not _formula_matches(existing, prop):
                raise ValueError(f"created metabolite identity conflict: {metabolite_id}")
            return existing
        metabolite = cobra.Metabolite(
            id=metabolite_id,
            name=prop.name or prop.compound_id,
            formula=prop.formula,
            charge=prop.charge,
            compartment=self.compartment,
        )
        if prop.source_mnxm_id:
            metabolite.annotation["metanetx.chemical"] = prop.source_mnxm_id
        if prop.inchikey:
            metabolite.annotation["inchi_key"] = prop.inchikey
        if re.fullmatch(r"C\d{5}", prop.compound_id):
            metabolite.annotation["kegg.compound"] = prop.compound_id
        self.model.add_metabolites([metabolite])
        self.created.add(metabolite.id)
        return metabolite


def _combination_id(
    candidate_id: str,
    hypotheses: Sequence[CompletedReactionHypothesis],
) -> str:
    return "RP2GEM:" + _canonical_sha256(
        {
            "candidate_id": candidate_id,
            "hypothesis_ids": sorted(item.hypothesis_id for item in hypotheses),
            "schema_version": RETROPATH_GEM_VALIDATION_SCHEMA,
        }
    )


def _add_rp2_reaction(
    model: cobra.Model,
    resolver: _MetaboliteResolver,
    hypothesis: CompletedReactionHypothesis,
    minimum_flux: float,
    upper_bound: float,
) -> cobra.Reaction:
    reaction = cobra.Reaction(
        id=sanitize_id(hypothesis.hypothesis_id),
        name=hypothesis.hypothesis_id,
        lower_bound=minimum_flux,
        upper_bound=upper_bound,
    )
    stoichiometry: dict[cobra.Metabolite, float] = defaultdict(float)
    for term in hypothesis.terms:
        metabolite = resolver.resolve(term.compound)
        sign = -1.0 if term.side == "left" else 1.0
        stoichiometry[metabolite] += sign * term.coefficient
    reaction.add_metabolites(dict(stoichiometry))
    reaction.annotation["metanetx.reaction"] = hypothesis.source_mnxr_id
    reaction.notes["rp2_stoichiometry_hypothesis"] = hypothesis.hypothesis_id
    reaction.notes["rr02_rule_id"] = hypothesis.rule_id
    model.add_reactions([reaction])
    imbalance = reaction.check_mass_balance()
    if imbalance:
        raise ValueError(
            f"model-mapped RP2 reaction is unbalanced: {hypothesis.step_id}: {imbalance}"
        )
    return reaction


def _force_endogenous_reaction(
    reactions: Sequence[cobra.Reaction],
    direction: str,
    minimum_flux: float,
) -> tuple[cobra.Reaction, float]:
    ordered = sorted(reactions, key=lambda item: item.id)
    if direction == "right_to_left":
        compatible = [item for item in ordered if item.lower_bound < 0]
        if not compatible:
            raise ValueError("endogenous reaction cannot carry right-to-left flux")
        reaction = compatible[0]
        reaction.upper_bound = min(reaction.upper_bound, -minimum_flux)
        return reaction, -1.0
    compatible = [item for item in ordered if item.upper_bound > 0]
    if not compatible:
        raise ValueError("endogenous reaction cannot carry left-to-right flux")
    reaction = compatible[0]
    reaction.lower_bound = max(reaction.lower_bound, minimum_flux)
    return reaction, 1.0


def _property_for_target(
    target: str,
    known: Mapping[str, CompoundProperty],
    hypotheses: Sequence[CompletedReactionHypothesis],
) -> CompoundProperty:
    if target in known:
        return known[target]
    for hypothesis in hypotheses:
        for term in hypothesis.terms:
            if term.compound.compound_id == target:
                return term.compound
    raise ValueError(f"target structure property is unavailable: {target}")


def _validate_combination(
    *,
    base_model: cobra.Model,
    target: str,
    candidate_rank: int,
    candidate_id: str,
    steps: Sequence[Mapping[str, str]],
    hypotheses: Sequence[CompletedReactionHypothesis],
    known_properties: Mapping[str, CompoundProperty],
    kegg_client: KeggRestClient,
    baseline_growth: float,
    required_growth: float,
    biomass_reaction_id: str,
    minimum_flux: float,
    run_fva: bool,
    combination_truncated: bool,
) -> _ModelValidation:
    model = base_model.copy()
    resolver = _MetaboliteResolver(model, DEFAULT_COMPARTMENT)
    reaction_index = build_kegg_reaction_index(model)
    metabolite_index = build_kegg_metabolite_index(model)
    hypothesis_by_step = {item.step_id: item for item in hypotheses}
    route_reactions: list[_RouteReaction] = []
    issues: list[str] = []
    combination_id = _combination_id(candidate_id, hypotheses)
    try:
        for row in sorted(
            steps,
            key=lambda item: _as_int(item.get("step_index"), "step_index", minimum=1),
        ):
            step_index = _as_int(row.get("step_index"), "step_index", minimum=1)
            step_id = str(row.get("step_id") or "").strip()
            source = str(row.get("step_source") or "").strip()
            if source == "retropath":
                hypothesis = hypothesis_by_step.get(step_id)
                if hypothesis is None:
                    raise ValueError(f"missing RP2 hypothesis for step {step_id}")
                reaction = _add_rp2_reaction(
                    model,
                    resolver,
                    hypothesis,
                    minimum_flux,
                    DEFAULT_REACTION_UPPER_BOUND,
                )
                route_reactions.append(
                    _RouteReaction(step_index, step_id, source, reaction.id, 1.0)
                )
                continue
            if source != "kegg_expansion":
                raise ValueError(f"unsupported candidate step source: {source}")
            reaction_ids = [
                value.strip()
                for value in str(row.get("reaction_option_ids") or "").split(";")
                if value.strip()
            ]
            if len(reaction_ids) != 1 or not re.fullmatch(r"R\d{5}", reaction_ids[0]):
                raise ValueError(f"KEGG step has no unique Rxxxxx reaction: {step_id}")
            reaction_id = reaction_ids[0]
            direction = str(row.get("direction") or "").strip()
            status = str(row.get("status") or "").strip()
            if status == "endogenous":
                reaction, sign = _force_endogenous_reaction(
                    reaction_index.get(reaction_id, []),
                    direction,
                    minimum_flux,
                )
            else:
                reaction, _, _, issue = add_kegg_reaction(
                    model,
                    route_scope=sanitize_id(combination_id),
                    reaction_id=reaction_id,
                    direction=direction,
                    client=kegg_client,
                    metabolite_index=metabolite_index,
                    preferred_compartment=DEFAULT_COMPARTMENT,
                    reaction_upper_bound=DEFAULT_REACTION_UPPER_BOUND,
                )
                if issue or reaction is None:
                    raise ValueError(issue or f"cannot add KEGG reaction {reaction_id}")
                reaction.lower_bound = max(reaction.lower_bound, minimum_flux)
                sign = 1.0
            route_reactions.append(
                _RouteReaction(step_index, step_id, source, reaction.id, sign)
            )

        target_property = _property_for_target(target, known_properties, hypotheses)
        target_metabolite = resolver.resolve(target_property)
        demand_id = sanitize_id(f"DM_P8_{target_metabolite.id}")
        if demand_id in model.reactions:
            demand = model.reactions.get_by_id(demand_id)
            demand.lower_bound = minimum_flux
            demand.upper_bound = DEFAULT_REACTION_UPPER_BOUND
        else:
            demand = model.add_boundary(
                target_metabolite,
                type="demand",
                reaction_id=demand_id,
                lb=minimum_flux,
                ub=DEFAULT_REACTION_UPPER_BOUND,
            )
        apply_growth_floor(model, biomass_reaction_id, required_growth)
        model.objective = demand
        model.objective_direction = "max"
        solution = model.optimize()
        fba_status = solution.status
        fluxes = solution.fluxes if fba_status == "optimal" else None
        target_flux = (
            float(fluxes.get(demand.id, float("nan")))
            if fluxes is not None
            else float("nan")
        )
        growth_flux = (
            float(fluxes.get(biomass_reaction_id, float("nan")))
            if fluxes is not None
            else float("nan")
        )
        pfba_status = "not_run"
        pfba_fluxes = None
        pfba_target_flux = float("nan")
        if fba_status == "optimal" and target_flux >= minimum_flux:
            try:
                pfba_solution = pfba(model, fraction_of_optimum=1.0)
                pfba_status = pfba_solution.status
                if pfba_status == "optimal":
                    pfba_fluxes = pfba_solution.fluxes
                    pfba_target_flux = float(pfba_fluxes.get(demand.id, float("nan")))
            except Exception as exc:
                pfba_status = f"error: {exc}"
        fva_status = "not_run"
        fva_df = None
        if run_fva and fba_status == "optimal" and target_flux >= minimum_flux:
            try:
                fva_df = flux_variability_analysis(
                    model,
                    reaction_list=[
                        demand.id,
                        biomass_reaction_id,
                        *[item.model_reaction_id for item in route_reactions],
                    ],
                    fraction_of_optimum=DEFAULT_FVA_FRACTION,
                    processes=1,
                )
                fva_status = "optimal"
            except Exception as exc:
                fva_status = f"error: {exc}"

        flux_rows: list[dict[str, Any]] = []
        blocked: list[str] = []
        for item in route_reactions:
            fba_flux = (
                float(fluxes.get(item.model_reaction_id, float("nan")))
                if fluxes is not None
                else float("nan")
            )
            pfba_flux = (
                float(pfba_fluxes.get(item.model_reaction_id, float("nan")))
                if pfba_fluxes is not None
                else float("nan")
            )
            directed_fba = item.direction_sign * fba_flux
            directed_pfba = item.direction_sign * pfba_flux
            active = math.isfinite(directed_fba) and directed_fba >= minimum_flux * 0.999
            if not active:
                blocked.append(item.step_id)
            fva_min = float("nan")
            fva_max = float("nan")
            if fva_df is not None and item.model_reaction_id in fva_df.index:
                fva_min = float(fva_df.loc[item.model_reaction_id, "minimum"])
                fva_max = float(fva_df.loc[item.model_reaction_id, "maximum"])
            flux_rows.append(
                {
                    "candidate_rank": candidate_rank,
                    "candidate_id": candidate_id,
                    "combination_id": combination_id,
                    "step_index": item.step_index,
                    "step_id": item.step_id,
                    "step_source": item.step_source,
                    "model_reaction_id": item.model_reaction_id,
                    "direction_sign": item.direction_sign,
                    "fba_flux": fba_flux,
                    "directed_fba_flux": directed_fba,
                    "pfba_flux": pfba_flux,
                    "directed_pfba_flux": directed_pfba,
                    "fva_minimum": fva_min,
                    "fva_maximum": fva_max,
                    "required_minimum_flux": minimum_flux,
                    "active": str(active).lower(),
                }
            )
        if fba_status != "optimal":
            validation_status = "FAIL_GEM_INFEASIBLE"
        elif target_flux < minimum_flux or growth_flux + 1e-12 < required_growth:
            validation_status = "FAIL_TARGET_OR_GROWTH_FLUX"
        elif blocked:
            validation_status = "FAIL_ROUTE_STEP_NO_FLUX"
        else:
            validation_status = "PASS_STRICT_ROUTE_FLUX"
        row = {
            "candidate_rank": candidate_rank,
            "candidate_id": candidate_id,
            "combination_id": combination_id,
            "candidate_status": "",
            "validation_status": validation_status,
            "stoichiometry_hypothesis_ids": ";".join(
                item.hypothesis_id for item in hypotheses
            ),
            "baseline_growth": baseline_growth,
            "required_growth": required_growth,
            "fba_status": fba_status,
            "target_flux": target_flux,
            "growth_flux": growth_flux,
            "pfba_status": pfba_status,
            "pfba_target_flux": pfba_target_flux,
            "fva_status": fva_status,
            "route_step_count": len(route_reactions),
            "active_route_step_count": len(route_reactions) - len(blocked),
            "blocked_route_step_ids": ";".join(blocked),
            "created_metabolite_ids": ";".join(sorted(resolver.created)),
            "issues": " | ".join(issues),
            "combination_truncated": str(combination_truncated).lower(),
        }
        return _ModelValidation(row, tuple(flux_rows))
    except Exception as exc:
        return _ModelValidation(
            {
                "candidate_rank": candidate_rank,
                "candidate_id": candidate_id,
                "combination_id": combination_id,
                "candidate_status": "",
                "validation_status": "FAIL_MODEL_PREPARATION",
                "stoichiometry_hypothesis_ids": ";".join(
                    item.hypothesis_id for item in hypotheses
                ),
                "baseline_growth": baseline_growth,
                "required_growth": required_growth,
                "fba_status": "not_run",
                "target_flux": float("nan"),
                "growth_flux": float("nan"),
                "pfba_status": "not_run",
                "pfba_target_flux": float("nan"),
                "fva_status": "not_run",
                "route_step_count": len(steps),
                "active_route_step_count": 0,
                "blocked_route_step_ids": ";".join(
                    str(row.get("step_id") or "") for row in steps
                ),
                "created_metabolite_ids": ";".join(sorted(resolver.created)),
                "issues": str(exc),
                "combination_truncated": str(combination_truncated).lower(),
            },
            tuple(),
        )


def _hypothesis_rows(
    rank: int,
    hypothesis: CompletedReactionHypothesis,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    summary = {
        "candidate_rank": rank,
        "candidate_id": hypothesis.candidate_id,
        "step_id": hypothesis.step_id,
        "hypothesis_id": hypothesis.hypothesis_id,
        "rule_id": hypothesis.rule_id,
        "source_mnxr_id": hypothesis.source_mnxr_id,
        "source_reference": hypothesis.source_reference,
        "source_equation": hypothesis.source_equation,
        "source_orientation": hypothesis.source_orientation,
        "evidence_grade": hypothesis.evidence_grade,
        "balance_status": hypothesis.balance_status,
        "cofactor_reconstruction_status": (
            hypothesis.cofactor_reconstruction_status
        ),
        "recovered_compound_ids": ";".join(hypothesis.recovered_compound_ids),
    }
    terms = [
        {
            "candidate_rank": rank,
            "candidate_id": hypothesis.candidate_id,
            "step_id": hypothesis.step_id,
            "hypothesis_id": hypothesis.hypothesis_id,
            "side": term.side,
            "coefficient": term.coefficient,
            "role": term.role,
            "compound_id": term.compound.compound_id,
            "source_mnxm_id": term.compound.source_mnxm_id,
            "name": term.compound.name,
            "formula": term.compound.formula,
            "charge": term.compound.charge,
            "inchi": term.compound.inchi,
            "inchikey": term.compound.inchikey,
            "smiles": term.compound.smiles,
            "xrefs": ";".join(term.compound.xrefs),
        }
        for term in hypothesis.terms
    ]
    return summary, terms


def validate_retropath_candidates(config: Any) -> dict[str, Any]:
    """Reconstruct and strictly validate selected P5 candidates."""

    if str(getattr(config, "validation_mode", "per")).strip().lower() not in {
        "per",
        "per-solution",
    }:
        raise ValueError("RetroPath validation supports only per-candidate mode")
    if str(getattr(config, "validation_cofactor_mode", "strict")).strip().lower() not in {
        "strict",
        "strict_l1",
    }:
        raise ValueError("RetroPath validation supports only strict cofactor mode")
    inputs = _load_candidate_inputs(config)
    available_ranks = [
        _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1)
        for row in inputs.routes
    ]
    selected_ranks = _selected_candidate_ranks(
        getattr(config, "retropath_candidates", None),
        available_ranks,
    )
    model_path = Path(config.model_path).expanduser().resolve()
    medium_path = Path(config.medium_path).expanduser().resolve()
    if not model_path.is_file() or model_path.suffix.lower() != ".json":
        raise FileNotFoundError(f"JSON GEM model not found: {model_path}")
    if not medium_path.is_file():
        raise FileNotFoundError(f"medium file not found: {medium_path}")
    base_model = cobra.io.load_json_model(str(model_path))
    base_model.medium = load_medium(medium_path)
    baseline_growth = safe_float(base_model.slim_optimize(), float("nan"))
    if not math.isfinite(baseline_growth) or baseline_growth <= DEFAULT_FLUX_THRESHOLD:
        raise ValueError("base model has no meaningful growth under the configured medium")
    required_growth = DEFAULT_GROWTH_FRACTION * baseline_growth
    biomass_reaction_id = get_primary_objective_reaction_id(base_model)
    known_properties = load_p2_compound_properties(inputs.compound_mapping_path)
    output_dir = inputs.retropath_dir / "gem_validation"
    kegg_client = KeggRestClient(Path(config.cache_dir).resolve() / "kegg")
    run_fva = not bool(getattr(config, "validation_skip_fva", False))

    hypothesis_rows: list[dict[str, Any]] = []
    term_rows: list[dict[str, Any]] = []
    rejection_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    flux_rows: list[dict[str, Any]] = []
    candidate_statuses: dict[int, str] = {}

    rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    mnxref_dir = Path(config.data_dir) / "retropath" / "mnxref" / "3.0"
    with MnxrefIndex(mnxref_dir, rules_path) as mnxref:
        for rank in selected_ranks:
            route = next(
                row
                for row in inputs.routes
                if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1)
                == rank
            )
            candidate_id = str(route.get("candidate_id") or "").strip()
            steps = tuple(
                sorted(
                    (
                        row
                        for row in inputs.steps
                        if str(row.get("candidate_id") or "").strip()
                        == candidate_id
                    ),
                    key=lambda item: _as_int(
                        item.get("step_index"),
                        "step_index",
                        minimum=1,
                    ),
                )
            )
            rp_steps = [
                row
                for row in steps
                if str(row.get("step_source") or "").strip() == "retropath"
            ]
            reconstructions: list[StepReconstructionResult] = []
            for step in rp_steps:
                try:
                    result = reconstruct_retropath_step(
                        step,
                        mnxref,
                        known_properties,
                        max_hypotheses=DEFAULT_MAX_STEP_HYPOTHESES,
                    )
                except Exception as exc:
                    result = StepReconstructionResult(
                        candidate_id=candidate_id,
                        step_id=str(step.get("step_id") or ""),
                        status="incomplete",
                        hypotheses=tuple(),
                        rejections=(
                            ReconstructionRejection(
                                candidate_id,
                                str(step.get("step_id") or ""),
                                str(step.get("rule_ids") or ""),
                                "",
                                "reconstruction_error",
                                str(exc),
                            ),
                        ),
                        truncated=False,
                    )
                reconstructions.append(result)
                for rejection in result.rejections:
                    rejection_rows.append(
                        {
                            "candidate_rank": rank,
                            **rejection.to_dict(),
                        }
                    )
                for hypothesis in result.hypotheses:
                    hypothesis_row, hypothesis_terms = _hypothesis_rows(
                        rank,
                        hypothesis,
                    )
                    hypothesis_rows.append(hypothesis_row)
                    term_rows.extend(hypothesis_terms)
            combinations, combination_truncated = enumerate_candidate_hypotheses(
                reconstructions,
                max_combinations=DEFAULT_MAX_CANDIDATE_COMBINATIONS,
            )
            if not combinations:
                candidate_statuses[rank] = "FAIL_NO_COMPLETE_STOICHIOMETRY"
                summary_rows.append(
                    {
                        "candidate_rank": rank,
                        "candidate_id": candidate_id,
                        "combination_id": "",
                        "candidate_status": candidate_statuses[rank],
                        "validation_status": "FAIL_STOICHIOMETRY",
                        "stoichiometry_hypothesis_ids": "",
                        "baseline_growth": baseline_growth,
                        "required_growth": required_growth,
                        "fba_status": "not_run",
                        "target_flux": float("nan"),
                        "growth_flux": float("nan"),
                        "pfba_status": "not_run",
                        "pfba_target_flux": float("nan"),
                        "fva_status": "not_run",
                        "route_step_count": len(steps),
                        "active_route_step_count": 0,
                        "blocked_route_step_ids": ";".join(
                            str(item.get("step_id") or "") for item in steps
                        ),
                        "created_metabolite_ids": "",
                        "issues": "one or more RP2 steps lack a complete balanced template",
                        "combination_truncated": "false",
                    }
                )
                continue
            validations: list[_ModelValidation] = []
            for combination in combinations:
                validation = _validate_combination(
                    base_model=base_model,
                    target=inputs.target_compound,
                    candidate_rank=rank,
                    candidate_id=candidate_id,
                    steps=steps,
                    hypotheses=combination,
                    known_properties=known_properties,
                    kegg_client=kegg_client,
                    baseline_growth=baseline_growth,
                    required_growth=required_growth,
                    biomass_reaction_id=biomass_reaction_id,
                    minimum_flux=DEFAULT_ROUTE_MIN_FLUX,
                    run_fva=run_fva,
                    combination_truncated=combination_truncated,
                )
                validations.append(validation)
            passed = any(
                item.row["validation_status"] == "PASS_STRICT_ROUTE_FLUX"
                for item in validations
            )
            candidate_statuses[rank] = (
                "PASS_STRICT_HYPOTHESIS_EXISTS"
                if passed
                else "FAIL_STRICT_GEM"
            )
            for validation in validations:
                row = dict(validation.row)
                row["candidate_status"] = candidate_statuses[rank]
                summary_rows.append(row)
                flux_rows.extend(validation.flux_rows)

        index_manifest = mnxref.manifest

    output_dir.mkdir(parents=True, exist_ok=True)
    artifact_texts = {
        STOICHIOMETRY_HYPOTHESES_FILE_NAME: _render_csv(
            HYPOTHESIS_COLUMNS,
            hypothesis_rows,
        ),
        STOICHIOMETRY_TERMS_FILE_NAME: _render_csv(TERM_COLUMNS, term_rows),
        VALIDATION_SUMMARY_FILE_NAME: _render_csv(SUMMARY_COLUMNS, summary_rows),
        VALIDATION_FLUX_FILE_NAME: _render_csv(FLUX_COLUMNS, flux_rows),
        REJECTED_HYPOTHESES_FILE_NAME: _render_csv(
            REJECTION_COLUMNS,
            rejection_rows,
        ),
    }
    artifact_records: dict[str, dict[str, Any]] = {}
    for name, text in artifact_texts.items():
        path = output_dir / name
        artifact_records[name] = {
            "path": str(path.resolve()),
            "sha256": _atomic_write_text(path, text),
        }
    manifest_path = output_dir / VALIDATION_MANIFEST_FILE_NAME
    manifest: dict[str, Any] = {
        "schema_version": RETROPATH_GEM_VALIDATION_SCHEMA,
        "created_at": datetime.now(UTC).isoformat(),
        "target_compound": inputs.target_compound,
        "expansion_depth": inputs.depth,
        "selected_candidate_ranks": list(selected_ranks),
        "candidate_statuses": {
            str(key): value for key, value in sorted(candidate_statuses.items())
        },
        "parameters": {
            "growth_fraction": DEFAULT_GROWTH_FRACTION,
            "route_minimum_flux": DEFAULT_ROUTE_MIN_FLUX,
            "strict_cofactor_mode": True,
            "max_step_hypotheses": DEFAULT_MAX_STEP_HYPOTHESES,
            "max_candidate_combinations": DEFAULT_MAX_CANDIDATE_COMBINATIONS,
            "fva_run": run_fva,
        },
        "inputs": {
            "pipeline_result": {
                "path": str(inputs.pipeline_path),
                "sha256": _sha256_file(inputs.pipeline_path),
            },
            "candidate_routes_sha256": inputs.pipeline["artifacts"][
                "candidate_routes"
            ]["sha256"],
            "candidate_steps_sha256": inputs.pipeline["artifacts"][
                "candidate_steps"
            ]["sha256"],
            "model_path": str(model_path),
            "model_sha256": _sha256_file(model_path),
            "medium_path": str(medium_path),
            "medium_sha256": _sha256_file(medium_path),
            "rr02_path": str(rules_path),
            "rr02_sha256": _sha256_file(rules_path),
            "mnxref_index_path": index_manifest["index_path"],
            "mnxref_index_sha256": index_manifest["index_sha256"],
        },
        "artifacts": artifact_records,
        "formal_promotion_allowed": any(
            row.get("validation_status") == "PASS_STRICT_ROUTE_FLUX"
            for row in summary_rows
        ),
    }
    manifest_sha256 = _atomic_write_json(manifest_path, manifest)
    from src.pathway_analyze.retropath_promotion import (
        materialize_retropath_solutions,
    )

    promotion = materialize_retropath_solutions(
        config,
        validation_manifest_path=manifest_path,
    )
    return {
        "ok": True,
        "search_engine": "retropath",
        "target_compound": inputs.target_compound,
        "expansion_depth": inputs.depth,
        "selected_candidate_ranks": list(selected_ranks),
        "candidate_statuses": manifest["candidate_statuses"],
        "validation_dir": str(output_dir.resolve()),
        "validation_manifest": str(manifest_path.resolve()),
        "validation_manifest_sha256": manifest_sha256,
        "summary_row_count": len(summary_rows),
        "flux_row_count": len(flux_rows),
        "formal_promotion_allowed": manifest["formal_promotion_allowed"],
        "formal_solution_ids": promotion["formal_solution_ids"],
        "solution_mappings": promotion["solution_mappings"],
        "promotion_manifest": promotion["promotion_manifest"],
        "promotion_manifest_sha256": promotion["promotion_manifest_sha256"],
    }


def run_retropath_validation(config: Any) -> dict[str, Any]:
    result = validate_retropath_candidates(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "DEFAULT_MAX_CANDIDATE_COMBINATIONS",
    "DEFAULT_ROUTE_MIN_FLUX",
    "RETROPATH_GEM_VALIDATION_SCHEMA",
    "run_retropath_validation",
    "validate_retropath_candidates",
]
