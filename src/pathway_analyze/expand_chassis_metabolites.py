"""按 KEGG 反应把底盘可提供化合物集合做同步分层扩展。

``A0`` 来自 GEM/FBA 生成的 ``producible_kegg_compounds.csv``。第 ``n``
层扩展只允许使用冻结的 ``A(n-1)`` 中全部底物，并把首次生成的化合物记为
``Fn``；累计集合 ``An = A0 ∪ F1 ... ∪ Fn``。每个新化合物都会保留生成它
的定向 KEGG 反应，供 gap 搜索命中扩展锚点后恢复完整反应路线。
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple

import pandas as pd

from src.pathway_analyze.kegg_gap_analyze import (
    ENDOGENOUS_DIRECTION_COMPARTMENTS,
    KeggRestClient,
    classify_reaction_resolution,
    directional_stoichiometry,
    load_endogenous_direction_index,
    make_reaction_option,
    option_is_endogenous,
)


EXPANSION_RULE_VERSION = "chassis_forward_expansion.v1"
COMPOUND_ID_PATTERN = re.compile(r"^C\d{5}$")
FRONTIER_FILE_TEMPLATE = "chassis_frontier_depth_{depth}.csv"
EXPANDED_FILE_TEMPLATE = "chassis_expanded_reachable_depth_{depth}.csv"
MANIFEST_FILE_NAME = "chassis_expansion_manifest.json"
DEFAULT_EXPANSION_PROGRESS_INTERVAL = 250

# 这些 ID 表示未明确指定的电子载体或大分子伙伴。即使它们碰巧出现在 A0，
# 也不能据此宣称一个具体 KEGG 反应已经具备全部可实现底物。
GENERIC_CARRIER_IDS = frozenset(
    {
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
)

EXPANSION_COLUMNS = (
    "source",
    "depth",
    "met_id",
    "met_name",
    "compartment",
    "kegg_id",
    "bridge_reaction_id",
    "bridge_direction",
    "bridge_substrate_ids",
    "bridge_product_ids",
    "bridge_is_endogenous",
    "thermo_direction",
    "oxygen_required",
    "nadph_burden",
    "sam_burden",
    "coa_burden",
    "electron_risk_level",
    "electron_risk_score",
    "enzyme_ecs",
    "ko_ids",
)


@dataclass(frozen=True)
class ForwardExpansionWitness:
    """一个扩展化合物及其最小深度定向 KEGG 生成反应。"""

    product_compound: str
    depth: int
    reaction_id: str
    direction: str
    substrate_compounds: Tuple[str, ...]
    product_compounds: Tuple[str, ...]
    is_endogenous: bool
    thermo_direction: str
    oxygen_required: bool
    nadph_burden: float
    sam_burden: float
    coa_burden: float
    electron_risk_level: str
    electron_risk_score: int
    enzyme_ecs: Tuple[str, ...]
    ko_ids: Tuple[str, ...]

    @property
    def signature(self) -> Tuple[str, str, str, Tuple[str, ...]]:
        return (
            self.product_compound,
            self.reaction_id,
            self.direction,
            self.substrate_compounds,
        )


@dataclass(frozen=True)
class ExpansionBundle:
    """供 gap 搜索读取的累计集合、深度和反应证据索引。"""

    depth: int
    base_compounds: frozenset[str]
    reachable_compounds: frozenset[str]
    depth_by_compound: Mapping[str, int]
    witnesses_by_product: Mapping[str, Tuple[ForwardExpansionWitness, ...]]
    expanded_file: Path
    manifest: Mapping[str, Any]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _frontier_path(output_dir: Path, depth: int) -> Path:
    return output_dir / FRONTIER_FILE_TEMPLATE.format(depth=depth)


def _expanded_path(output_dir: Path, depth: int) -> Path:
    return output_dir / EXPANDED_FILE_TEMPLATE.format(depth=depth)


def _manifest_path(output_dir: Path) -> Path:
    return output_dir / MANIFEST_FILE_NAME


def _normalize_id_tuple(values: Iterable[Any]) -> Tuple[str, ...]:
    return tuple(sorted({str(value).strip() for value in values if str(value).strip()}))


def _split_ids(value: Any) -> Tuple[str, ...]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return tuple()
    return _normalize_id_tuple(str(value).split(";"))


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _as_int(value: Any, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _witness_sort_key(witness: ForwardExpansionWitness) -> Tuple[Any, ...]:
    thermo_penalty = {"favored": 0, "neutral": 1}.get(
        witness.thermo_direction,
        2,
    )
    return (
        thermo_penalty,
        witness.electron_risk_score,
        0 if witness.is_endogenous else 1,
        witness.sam_burden,
        witness.nadph_burden,
        witness.coa_burden,
        witness.reaction_id,
        witness.direction,
        witness.substrate_compounds,
    )


def _witness_to_row(witness: ForwardExpansionWitness) -> Dict[str, Any]:
    return {
        "source": "kegg_expansion",
        "depth": witness.depth,
        "met_id": "",
        "met_name": "",
        "compartment": "",
        "kegg_id": witness.product_compound,
        "bridge_reaction_id": witness.reaction_id,
        "bridge_direction": witness.direction,
        "bridge_substrate_ids": ";".join(witness.substrate_compounds),
        "bridge_product_ids": ";".join(witness.product_compounds),
        "bridge_is_endogenous": witness.is_endogenous,
        "thermo_direction": witness.thermo_direction,
        "oxygen_required": witness.oxygen_required,
        "nadph_burden": witness.nadph_burden,
        "sam_burden": witness.sam_burden,
        "coa_burden": witness.coa_burden,
        "electron_risk_level": witness.electron_risk_level,
        "electron_risk_score": witness.electron_risk_score,
        "enzyme_ecs": ";".join(witness.enzyme_ecs),
        "ko_ids": ";".join(witness.ko_ids),
    }


def _row_to_witness(row: Mapping[str, Any]) -> ForwardExpansionWitness:
    return ForwardExpansionWitness(
        product_compound=str(row.get("kegg_id", "")).strip(),
        depth=_as_int(row.get("depth")),
        reaction_id=str(row.get("bridge_reaction_id", "")).strip(),
        direction=str(row.get("bridge_direction", "")).strip(),
        substrate_compounds=_split_ids(row.get("bridge_substrate_ids")),
        product_compounds=_split_ids(row.get("bridge_product_ids")),
        is_endogenous=_as_bool(row.get("bridge_is_endogenous")),
        thermo_direction=str(row.get("thermo_direction", "neutral")).strip(),
        oxygen_required=_as_bool(row.get("oxygen_required")),
        nadph_burden=_as_float(row.get("nadph_burden")),
        sam_burden=_as_float(row.get("sam_burden")),
        coa_burden=_as_float(row.get("coa_burden")),
        electron_risk_level=str(row.get("electron_risk_level", "none")).strip(),
        electron_risk_score=_as_int(row.get("electron_risk_score")),
        enzyme_ecs=_split_ids(row.get("enzyme_ecs")),
        ko_ids=_split_ids(row.get("ko_ids")),
    )


def _load_base_dataframe(base_path: Path) -> tuple[pd.DataFrame, set[str]]:
    if not base_path.is_file():
        raise FileNotFoundError(
            "Chassis A0 file not found. Run the chassis command first: "
            f"{base_path}"
        )
    df = pd.read_csv(base_path)
    if "kegg_id" not in df.columns:
        raise ValueError(f"Missing kegg_id column in chassis A0 file: {base_path}")
    df = df.copy()
    df["kegg_id"] = df["kegg_id"].astype(str).str.strip()
    base_compounds = {
        compound_id
        for compound_id in df["kegg_id"]
        if COMPOUND_ID_PATTERN.fullmatch(compound_id)
    }
    if not base_compounds:
        raise ValueError(f"Chassis A0 contains no Cxxxxx KEGG compounds: {base_path}")
    return df, base_compounds


def _base_rows(base_df: pd.DataFrame) -> list[Dict[str, Any]]:
    rows: list[Dict[str, Any]] = []
    for raw_row in base_df.to_dict(orient="records"):
        row = {column: "" for column in EXPANSION_COLUMNS}
        for column in ("source", "met_id", "met_name", "compartment", "kegg_id"):
            if column in raw_row:
                row[column] = raw_row[column]
        row["source"] = "producible"
        row["depth"] = 0
        rows.append(row)
    return rows


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=EXPANSION_COLUMNS).to_csv(
        path,
        index=False,
        encoding="utf-8-sig",
    )


def _read_manifest(path: Path) -> Dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _manifest_matches(
    manifest: Mapping[str, Any] | None,
    *,
    a0_sha256: str,
) -> bool:
    return bool(
        manifest
        and manifest.get("algorithm_version") == EXPANSION_RULE_VERSION
        and manifest.get("a0_sha256") == a0_sha256
        and manifest.get("direction_policy") == "bidirectional_reject_disfavored"
        and manifest.get("compound_policy") == "Cxxxxx_strict_all_substrates"
    )


def _clean_generated_outputs(output_dir: Path) -> None:
    for pattern in (
        "chassis_frontier_depth_*.csv",
        "chassis_expanded_reachable_depth_*.csv",
        MANIFEST_FILE_NAME,
    ):
        for path in output_dir.glob(pattern):
            if path.is_file():
                path.unlink()


def _load_witnesses_from_dataframe(
    df: pd.DataFrame,
) -> Dict[str, Tuple[ForwardExpansionWitness, ...]]:
    mutable: Dict[str, list[ForwardExpansionWitness]] = defaultdict(list)
    if "depth" not in df.columns:
        return {}
    for row in df.to_dict(orient="records"):
        if _as_int(row.get("depth")) <= 0:
            continue
        witness = _row_to_witness(row)
        if (
            not COMPOUND_ID_PATTERN.fullmatch(witness.product_compound)
            or not witness.reaction_id
            or witness.direction not in {"left_to_right", "right_to_left"}
        ):
            continue
        mutable[witness.product_compound].append(witness)
    return {
        product: tuple(sorted({item.signature: item for item in items}.values(), key=_witness_sort_key))
        for product, items in mutable.items()
    }


def load_expansion_bundle(
    *,
    base_path: Path,
    output_dir: Path,
    depth: int,
) -> ExpansionBundle:
    """读取并校验指定深度的累计扩展集合。"""

    if depth < 0:
        raise ValueError("depth must be greater than or equal to 0")
    base_df, base_compounds = _load_base_dataframe(base_path)
    del base_df
    a0_sha256 = _sha256_file(base_path)
    if depth == 0:
        return ExpansionBundle(
            depth=0,
            base_compounds=frozenset(base_compounds),
            reachable_compounds=frozenset(base_compounds),
            depth_by_compound={compound_id: 0 for compound_id in base_compounds},
            witnesses_by_product={},
            expanded_file=base_path,
            manifest={},
        )

    manifest_path = _manifest_path(output_dir)
    manifest = _read_manifest(manifest_path)
    if not _manifest_matches(manifest, a0_sha256=a0_sha256):
        raise ValueError(
            "Requested chassis expansion is missing or stale. "
            f"Run the expand command for depth {depth} first."
        )
    if _as_int(manifest.get("max_depth")) < depth:
        raise ValueError(
            f"Expansion depth {depth} has not been generated. "
            f"Run the expand command for depth {depth} first."
        )

    expanded_path = _expanded_path(output_dir, depth)
    if not expanded_path.is_file():
        raise FileNotFoundError(
            f"Missing cumulative expansion file for depth {depth}: {expanded_path}"
        )
    df = pd.read_csv(expanded_path)
    required = {
        "source",
        "depth",
        "kegg_id",
        "bridge_reaction_id",
        "bridge_direction",
        "bridge_substrate_ids",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            f"Expansion file is missing required columns {missing}: {expanded_path}"
        )

    depth_by_compound = {compound_id: 0 for compound_id in base_compounds}
    for row in df.to_dict(orient="records"):
        compound_id = str(row.get("kegg_id", "")).strip()
        compound_depth = _as_int(row.get("depth"))
        if not COMPOUND_ID_PATTERN.fullmatch(compound_id):
            continue
        previous = depth_by_compound.get(compound_id)
        if previous is None or compound_depth < previous:
            depth_by_compound[compound_id] = compound_depth
    witnesses = _load_witnesses_from_dataframe(df)
    return ExpansionBundle(
        depth=depth,
        base_compounds=frozenset(base_compounds),
        reachable_compounds=frozenset(depth_by_compound),
        depth_by_compound=depth_by_compound,
        witnesses_by_product=witnesses,
        expanded_file=expanded_path,
        manifest=manifest or {},
    )


def _reaction_witnesses(
    *,
    reaction: Any,
    depth: int,
    available_compounds: set[str],
    previous_frontier: set[str],
    endogenous_reactions: set[str],
    endogenous_direction_index: Mapping[str, frozenset[str]],
    rejection_counts: Counter[str] | None = None,
) -> Tuple[ForwardExpansionWitness, ...]:
    resolution = classify_reaction_resolution(reaction)
    if resolution.is_incomplete or resolution.is_multistep or resolution.hard_blocker:
        if rejection_counts is not None:
            rejection_counts[f"resolution:{resolution.reason}"] += 1
        return tuple()

    witnesses: list[ForwardExpansionWitness] = []
    for direction in ("left_to_right", "right_to_left"):
        consumed, produced = directional_stoichiometry(reaction, direction)
        substrates = _normalize_id_tuple(compound_id for compound_id, _ in consumed)
        products = _normalize_id_tuple(compound_id for compound_id, _ in produced)
        if not substrates or not products:
            if rejection_counts is not None:
                rejection_counts["empty_substrates_or_products"] += 1
            continue
        if not all(COMPOUND_ID_PATTERN.fullmatch(item) for item in (*substrates, *products)):
            if rejection_counts is not None:
                rejection_counts["non_Cxxxxx_compound"] += 1
            continue
        if GENERIC_CARRIER_IDS.intersection((*substrates, *products)):
            if rejection_counts is not None:
                rejection_counts["generic_carrier"] += 1
            continue
        substrate_set = set(substrates)
        if not substrate_set.issubset(available_compounds):
            if rejection_counts is not None:
                rejection_counts["missing_required_substrate"] += 1
            continue
        if not substrate_set.intersection(previous_frontier):
            if rejection_counts is not None:
                rejection_counts["not_driven_by_previous_frontier"] += 1
            continue

        for product_id in products:
            if product_id in available_compounds or product_id in substrate_set:
                if rejection_counts is not None:
                    rejection_counts["product_already_reachable"] += 1
                continue
            option = make_reaction_option(
                compound_id=product_id,
                reaction=reaction,
                direction=direction,
                ignored_common_compounds=set(),
            )
            if option.screening.thermo_direction == "disfavored":
                if rejection_counts is not None:
                    rejection_counts["thermodynamically_disfavored"] += 1
                continue
            witnesses.append(
                ForwardExpansionWitness(
                    product_compound=product_id,
                    depth=depth,
                    reaction_id=reaction.reaction_id,
                    direction=direction,
                    substrate_compounds=substrates,
                    product_compounds=products,
                    is_endogenous=option_is_endogenous(
                        option,
                        endogenous_reactions,
                        endogenous_direction_index,
                    ),
                    thermo_direction=option.screening.thermo_direction,
                    oxygen_required=option.screening.oxygen_required,
                    nadph_burden=option.screening.nadph_burden,
                    sam_burden=option.screening.sam_burden,
                    coa_burden=option.screening.coa_burden,
                    electron_risk_level=option.electron_requirement.risk_level,
                    electron_risk_score=option.electron_requirement.risk_score,
                    enzyme_ecs=tuple(reaction.enzyme_ecs),
                    ko_ids=tuple(reaction.ko_ids),
                )
            )
    return tuple(witnesses)


def ensure_expansion_depth(
    *,
    base_path: Path,
    output_dir: Path,
    model_path: Path,
    cache_dir: Path,
    requested_depth: int,
    client: KeggRestClient | None = None,
    progress_interval: int = DEFAULT_EXPANSION_PROGRESS_INTERVAL,
) -> ExpansionBundle:
    """生成或复用直到 ``requested_depth`` 的同步 KEGG 扩展。"""

    expansion_started = time.perf_counter()
    if requested_depth < 1:
        raise ValueError("expand depth must be a positive integer")
    if progress_interval < 1:
        raise ValueError("expand progress interval must be a positive integer")
    base_df, base_compounds = _load_base_dataframe(base_path)
    valid_a0_rows = sum(
        bool(COMPOUND_ID_PATTERN.fullmatch(str(value).strip()))
        for value in base_df["kegg_id"]
    )
    print(
        "[INFO] chassis A0 loaded: "
        f"rows={len(base_df)}, valid_Cxxxxx_rows={valid_a0_rows}, "
        f"unique_compounds={len(base_compounds)}, "
        f"duplicate_rows={max(0, valid_a0_rows - len(base_compounds))}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    a0_sha256 = _sha256_file(base_path)
    manifest_path = _manifest_path(output_dir)
    manifest = _read_manifest(manifest_path)

    manifest_matches = _manifest_matches(manifest, a0_sha256=a0_sha256)
    if not manifest_matches:
        if manifest is None:
            print("[INFO] no reusable expansion manifest found; starting from A0")
        else:
            print(
                "[WARNING] expansion manifest is stale or incompatible; "
                "rebuilding generated expansion files"
            )
        _clean_generated_outputs(output_dir)
        manifest = None

    existing_depth = _as_int((manifest or {}).get("max_depth"))
    if existing_depth >= requested_depth:
        try:
            cached_bundle = load_expansion_bundle(
                base_path=base_path,
                output_dir=output_dir,
                depth=requested_depth,
            )
            print(
                "[INFO] reusing cached chassis expansion: "
                f"requested_depth={requested_depth}, cached_max_depth={existing_depth}, "
                f"reachable_compounds={len(cached_bundle.reachable_compounds)}, "
                f"elapsed={time.perf_counter() - expansion_started:.1f}s"
            )
            return cached_bundle
        except (FileNotFoundError, ValueError):
            # manifest 存在但派生 CSV 缺失或损坏时，expand 命令负责自愈重建；
            # gap 读取端仍会直接报错，避免搜索阶段隐式联网。
            _clean_generated_outputs(output_dir)
            manifest = None
            existing_depth = 0
            print(
                "[WARNING] cached expansion files are missing or invalid; "
                "rebuilding from A0"
            )

    if client is None:
        client = KeggRestClient(cache_dir)
    print("[INFO] loading endogenous GEM reaction directions")
    endogenous_direction_index = load_endogenous_direction_index(
        model_path,
        allowed_compartments=ENDOGENOUS_DIRECTION_COMPARTMENTS,
    )
    endogenous_reactions = set(endogenous_direction_index)
    endogenous_direction_count = sum(
        len(directions) for directions in endogenous_direction_index.values()
    )
    print(
        "[INFO] endogenous reaction index loaded: "
        f"reactions={len(endogenous_reactions)}, "
        f"directions={endogenous_direction_count}"
    )
    print("[INFO] loading KEGG compound-reaction index")
    compound_reaction_index = client.get_compound_reaction_index()
    compound_reaction_links = sum(
        len(reaction_ids) for reaction_ids in compound_reaction_index.values()
    )
    print(
        "[INFO] KEGG compound-reaction index loaded: "
        f"compounds={len(compound_reaction_index)}, "
        f"compound_reaction_links={compound_reaction_links}"
    )

    if existing_depth > 0:
        try:
            existing_bundle = load_expansion_bundle(
                base_path=base_path,
                output_dir=output_dir,
                depth=existing_depth,
            )
        except (FileNotFoundError, ValueError):
            _clean_generated_outputs(output_dir)
            existing_depth = 0
            manifest = None
            existing_bundle = None
    else:
        existing_bundle = None

    if existing_bundle is not None:
        available = set(existing_bundle.reachable_compounds)
        previous_frontier = {
            compound_id
            for compound_id, compound_depth in existing_bundle.depth_by_compound.items()
            if compound_depth == existing_depth
        }
        all_witnesses = [
            witness
            for witnesses in existing_bundle.witnesses_by_product.values()
            for witness in witnesses
        ]
        layers_meta: Dict[str, Any] = dict((manifest or {}).get("layers", {}))
        print(
            "[INFO] resuming chassis expansion: "
            f"existing_depth={existing_depth}, "
            f"reachable_compounds={len(available)}, "
            f"previous_frontier={len(previous_frontier)}"
        )
    else:
        available = set(base_compounds)
        previous_frontier = set(base_compounds)
        all_witnesses = []
        layers_meta = {}

    for depth in range(existing_depth + 1, requested_depth + 1):
        layer_started = time.perf_counter()
        input_frontier_count = len(previous_frontier)
        candidate_reaction_ids = sorted(
            {
                reaction_id
                for compound_id in previous_frontier
                for reaction_id in compound_reaction_index.get(compound_id, tuple())
            }
        )
        print(
            f"[INFO] chassis expansion depth {depth} started: "
            f"input_frontier={input_frontier_count}, "
            f"available_compounds={len(available)}, "
            f"candidate_reactions={len(candidate_reaction_ids)}"
        )
        prefetch_started = time.perf_counter()
        for start in range(0, len(candidate_reaction_ids), progress_interval):
            end = min(start + progress_interval, len(candidate_reaction_ids))
            client.prefetch_reactions(candidate_reaction_ids[start:end])
            prefetch_progress_elapsed = time.perf_counter() - prefetch_started
            prefetch_rate = (
                end / prefetch_progress_elapsed
                if prefetch_progress_elapsed > 0
                else 0.0
            )
            prefetch_remaining = len(candidate_reaction_ids) - end
            prefetch_eta = (
                prefetch_remaining / prefetch_rate if prefetch_rate > 0 else 0.0
            )
            print(
                f"[INFO] depth {depth} KEGG reaction prefetch: "
                f"{end}/{len(candidate_reaction_ids)} "
                f"({100.0 * end / len(candidate_reaction_ids):.1f}%), "
                f"elapsed={prefetch_progress_elapsed:.1f}s, "
                f"ETA={prefetch_eta:.1f}s"
            )
        if not candidate_reaction_ids:
            print(f"[INFO] depth {depth} KEGG reaction prefetch: no candidates")
        prefetch_elapsed = time.perf_counter() - prefetch_started
        layer_witnesses: list[ForwardExpansionWitness] = []
        rejection_counts: Counter[str] = Counter()
        loaded_reaction_count = 0
        missing_reaction_count = 0
        for reaction_id in candidate_reaction_ids:
            reaction = client.try_get_reaction(reaction_id)
            if reaction is None:
                missing_reaction_count += 1
                continue
            loaded_reaction_count += 1
            layer_witnesses.extend(
                _reaction_witnesses(
                    reaction=reaction,
                    depth=depth,
                    available_compounds=available,
                    previous_frontier=previous_frontier,
                    endogenous_reactions=endogenous_reactions,
                    endogenous_direction_index=endogenous_direction_index,
                    rejection_counts=rejection_counts,
                )
            )

        # 一个产物仅登记在首次可达的层；同层的多个生成证据全部保留。
        unique_witnesses = {
            witness.signature: witness
            for witness in layer_witnesses
            if witness.product_compound not in available
        }
        ordered_layer = sorted(unique_witnesses.values(), key=lambda item: (
            item.product_compound,
            *_witness_sort_key(item),
        ))
        new_frontier = {witness.product_compound for witness in ordered_layer}
        eligible_reaction_count = len(
            {witness.reaction_id for witness in ordered_layer}
        )
        all_witnesses.extend(ordered_layer)
        available.update(new_frontier)

        frontier_path = _frontier_path(output_dir, depth)
        expanded_path = _expanded_path(output_dir, depth)
        _write_rows(frontier_path, [_witness_to_row(item) for item in ordered_layer])
        cumulative_rows = [
            *_base_rows(base_df),
            *[
                _witness_to_row(item)
                for item in sorted(
                    all_witnesses,
                    key=lambda witness: (
                        witness.depth,
                        witness.product_compound,
                        *_witness_sort_key(witness),
                    ),
                )
            ],
        ]
        _write_rows(expanded_path, cumulative_rows)
        layers_meta[str(depth)] = {
            "frontier_file": frontier_path.name,
            "expanded_file": expanded_path.name,
            "input_frontier_count": input_frontier_count,
            "candidate_reaction_count": len(candidate_reaction_ids),
            "loaded_reaction_count": loaded_reaction_count,
            "missing_reaction_count": missing_reaction_count,
            "eligible_reaction_count": eligible_reaction_count,
            "raw_witness_count": len(layer_witnesses),
            "frontier_compound_count": len(new_frontier),
            "witness_count": len(ordered_layer),
            "cumulative_compound_count": len(available),
            "rejection_counts": dict(sorted(rejection_counts.items())),
            "prefetch_elapsed_seconds": round(prefetch_elapsed, 3),
            "layer_elapsed_seconds": round(
                time.perf_counter() - layer_started,
                3,
            ),
        }
        print(
            f"[INFO] chassis expansion depth {depth} completed: "
            f"loaded_reactions={loaded_reaction_count}, "
            f"missing_reactions={missing_reaction_count}, "
            f"eligible_reactions={eligible_reaction_count}, "
            f"new_compounds={len(new_frontier)}, "
            f"reaction_witnesses={len(ordered_layer)}, "
            f"cumulative_compounds={len(available)}, "
            f"prefetch={prefetch_elapsed:.1f}s, "
            f"elapsed={time.perf_counter() - layer_started:.1f}s"
        )
        if rejection_counts:
            rejection_summary = ", ".join(
                f"{reason}={count}"
                for reason, count in sorted(rejection_counts.items())
            )
            print(
                f"[INFO] chassis expansion depth {depth} filters: "
                f"{rejection_summary}"
            )
        if not new_frontier and depth < requested_depth:
            print(
                f"[WARNING] chassis expansion depth {depth} produced an empty "
                "frontier; subsequent depths will remain empty"
            )
        print(f"[INFO] frontier file written to: {frontier_path}")
        print(f"[INFO] cumulative expansion file written to: {expanded_path}")
        previous_frontier = new_frontier

    manifest_payload = {
        "algorithm_version": EXPANSION_RULE_VERSION,
        "a0_sha256": a0_sha256,
        "a0_file": str(base_path.resolve()),
        "max_depth": requested_depth,
        "direction_policy": "bidirectional_reject_disfavored",
        "compound_policy": "Cxxxxx_strict_all_substrates",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "layers": layers_meta,
    }
    manifest_path.write_text(
        json.dumps(manifest_payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"[INFO] expansion manifest written to: {manifest_path}")
    print(
        "[INFO] chassis expansion generation completed: "
        f"requested_depth={requested_depth}, "
        f"reachable_compounds={len(available)}, "
        f"elapsed={time.perf_counter() - expansion_started:.1f}s"
    )
    return load_expansion_bundle(
        base_path=base_path,
        output_dir=output_dir,
        depth=requested_depth,
    )


def expand_chassis_metabolites(config: Any) -> dict[str, Any]:
    """使用 ``config.depth`` 生成指定累计深度的底盘 KEGG 扩展集合。"""

    analysis_started = time.perf_counter()
    depth = int(getattr(config, "depth", 0))
    base_path = Path(config.chassis_producible_csv).expanduser().resolve()
    output_dir = Path(config.chassis_output_path).expanduser().resolve()
    model_path = Path(config.model_path).expanduser().resolve()
    cache_dir = Path(config.cache_dir).expanduser().resolve() / "kegg"
    progress_interval = int(
        getattr(
            config,
            "expand_progress_interval",
            DEFAULT_EXPANSION_PROGRESS_INTERVAL,
        )
    )
    if depth < 1:
        raise ValueError("expand command requires depth >= 1")
    if progress_interval < 1:
        raise ValueError("expand_progress_interval must be greater than or equal to 1")
    if not model_path.is_file():
        raise FileNotFoundError(f"GEM model file not found: {model_path}")

    print("[INFO] starting chassis KEGG expansion")
    print(f"[INFO] chassis A0 file: {base_path}")
    print(f"[INFO] GEM model: {model_path}")
    print(f"[INFO] KEGG cache directory: {cache_dir}")
    print(f"[INFO] expansion output directory: {output_dir}")
    print(f"[INFO] requested expansion depth: {depth}")
    print(f"[INFO] KEGG prefetch progress interval: {progress_interval}")
    bundle = ensure_expansion_depth(
        base_path=base_path,
        output_dir=output_dir,
        model_path=model_path,
        cache_dir=cache_dir,
        requested_depth=depth,
        progress_interval=progress_interval,
    )
    frontier_path = _frontier_path(output_dir, depth)
    manifest_path = _manifest_path(output_dir)
    elapsed = time.perf_counter() - analysis_started
    frontier_compound_count = sum(
        1 for value in bundle.depth_by_compound.values() if value == depth
    )
    expansion_compound_count = (
        len(bundle.reachable_compounds) - len(bundle.base_compounds)
    )
    print(
        "[INFO] chassis KEGG expansion completed: "
        f"A0={len(bundle.base_compounds)}, "
        f"added={expansion_compound_count}, "
        f"frontier_depth_{depth}={frontier_compound_count}, "
        f"reachable={len(bundle.reachable_compounds)}, "
        f"elapsed={elapsed:.1f}s"
    )
    print(f"[INFO] requested frontier file: {frontier_path.resolve()}")
    print(f"[INFO] requested cumulative file: {bundle.expanded_file.resolve()}")
    print(f"[INFO] expansion manifest file: {manifest_path.resolve()}")
    return {
        "ok": True,
        "depth": depth,
        "progress_interval": progress_interval,
        "base_compound_count": len(bundle.base_compounds),
        "expanded_compound_count": expansion_compound_count,
        "reachable_compound_count": len(bundle.reachable_compounds),
        "frontier_compound_count": frontier_compound_count,
        "elapsed_seconds": elapsed,
        "frontier_file": str(frontier_path.resolve()),
        "expanded_reachable_file": str(bundle.expanded_file.resolve()),
        "manifest_file": str(manifest_path.resolve()),
    }


def run_expand(config: Any) -> dict[str, Any]:
    """命令入口；参数解析和 ``config.depth`` 赋值由 ``main.py`` 负责。"""

    result = expand_chassis_metabolites(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result
