"""Write one user-selected main-enzyme set to ``design_manifest.json``.

This module is deliberately a commit boundary.  The main-protein selection
package produces ranked artifacts, while this writer validates that those
artifacts still describe the route currently stored in the design manifest
before recording the user's choice.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from filelock import FileLock, Timeout

from src.main_protein_selection.common import get_solution_steps, heterologous_requirements
from src.main_protein_selection.genome_ec_filter import (
    UNSPECIFIED_GENOME, read_current_genome_selection, validate_genome_ec_filter,
)
from src.main_protein_selection.models import (
    SUPPORTED_MAIN_ENZYME_SETS_SCHEMAS,
    parse_main_enzyme_sets,
    MainEnzymeSet,
    MainEnzymeSetsResult,
)
from src.pathway_analyze.target_id import validate_target_compound_id
from src.write_manifest.store import (
    read_design_manifest,
    manifest_update_lock,
    update_design_manifest,
)


MAIN_ENZYME_MANIFEST_SCHEMA_VERSION = "main_enzyme_manifest_selection.v6"

MAIN_ENZYME_DOWNSTREAM_SECTIONS = (
    "protein_selection",
    "auxiliary_protein_selection",
    "enzyme_system_selection",
    "electron_inference",
    "electron_balance",
    "cds_selection",
    "expression_box_selection",
    "expression_cassette_assembly",
    "expression_parts_draft",
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
    "final_design_report",
)

_SUPPORTED_DIRECTIONS = {
    "supported",
    "verified",
    "compatible",
    "forward",
    "reversible",
}
_SUPPORTED_SPECIFICITIES = {"exact", "supported"}


def _selection_dir(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "main_protein_selection"
    )


def _sets_path(config: Any) -> Path:
    return _selection_dir(config) / "main_enzyme_sets.json"


def _candidate_selection_path(config: Any) -> Path:
    return _selection_dir(config) / "main_enzyme_selection.json"


def _read_sets(path: Path) -> MainEnzymeSetsResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶组合，请先运行：main-enzyme-sets -i <输入文件>"
        )
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"主酶组合结果必须是 JSON 对象，请重新运行 main-enzyme-sets：{path}")
    if payload.get("schema_version") not in SUPPORTED_MAIN_ENZYME_SETS_SCHEMAS:
        raise ValueError(
            "主酶组合使用旧版结果格式，请重新运行 main-enzyme-sets"
        )
    try:
        return parse_main_enzyme_sets(payload)
    except ValueError as exc:
        raise ValueError(f"主酶组合结果格式无效：{path}") from exc


def _positive_set_id(config: Any) -> int:
    raw_value = getattr(config, "main_enzyme_set", None)
    if raw_value is None:
        raise ValueError(
            "未指定主酶组合，请使用 write --main-enzyme-set N"
        )
    try:
        set_id = int(raw_value)
    except (TypeError, ValueError) as exc:
        raise ValueError("主酶组合编号必须是正整数") from exc
    if set_id < 1:
        raise ValueError("主酶组合编号必须是正整数")
    return set_id


def _current_solution_context(
    manifest: Mapping[str, Any],
) -> tuple[int, int, list[dict[str, Any]]]:
    solution_id, steps = get_solution_steps(dict(manifest))
    solution = manifest.get("solution")
    if not isinstance(solution, Mapping):
        raise ValueError(
            "manifest 尚未写入路线，请先运行 write --solution N"
        )
    try:
        expansion_depth = int(solution.get("expansion_depth") or 0)
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest 中的路线扩展深度无效") from exc
    return (
        solution_id,
        expansion_depth,
        steps,
    )


def _select_set(
    result: MainEnzymeSetsResult,
    set_id: int,
) -> MainEnzymeSet:
    selected = next(
        (item for item in result.sets if item.set_id == set_id),
        None,
    )
    if selected is not None:
        return selected

    available = [item.set_id for item in result.sets]
    if not available:
        reasons = "；".join(result.blocking_reasons)
        suffix = f"；原因：{reasons}" if reasons else ""
        raise ValueError(f"当前没有可写入的主酶组合{suffix}")
    raise ValueError(
        f"不存在主酶组合 {set_id}；可用组合编号：{available}"
    )


def _review_items(
    result: MainEnzymeSetsResult,
    selected: MainEnzymeSet,
) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []

    fit_steps = [
        item.step_index
        for item in selected.step_assignments
        if item.reaction_fit_status != "verified"
    ]
    if fit_steps:
        items.append(
            {
                "review_type": "reaction_fit",
                "status": "pending",
                "step_indexes": sorted(fit_steps),
                "reason": "部分主酶与目标反应的匹配仍需复核",
            }
        )

    direction_steps = [
        item.step_index
        for item in selected.step_assignments
        if item.direction_verdict.lower() not in _SUPPORTED_DIRECTIONS
    ]
    if direction_steps:
        items.append(
            {
                "review_type": "reaction_direction",
                "status": "pending",
                "step_indexes": sorted(direction_steps),
                "reason": "部分主酶缺少充分的生理反应方向证据",
            }
        )

    specificity_steps = [
        item.step_index
        for item in selected.step_assignments
        if item.specificity_status.lower() not in _SUPPORTED_SPECIFICITIES
    ]
    if specificity_steps:
        items.append(
            {
                "review_type": "reaction_specificity",
                "status": "pending",
                "step_indexes": sorted(specificity_steps),
                "reason": "部分主酶的底物或产物特异性仍需复核",
            }
        )


    if not result.search_complete:
        items.append(
            {
                "review_type": "combination_search",
                "status": "pending",
                "step_indexes": [],
                "reason": (
                    "主酶组合搜索未完成，当前选择可能不是全局最优组合"
                ),
            }
        )

    if selected.status == "review_required" and not items:
        items.append(
            {
                "review_type": "general",
                "status": "pending",
                "step_indexes": selected.review_required_step_indexes,
                "reason": "该主酶组合仍包含需要人工复核的证据",
            }
        )
    return items


def _manifest_payload(
    result: MainEnzymeSetsResult,
    selected: MainEnzymeSet,
    solution_prediction: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    unresolved_reviews = _review_items(result, selected)
    if solution_prediction and solution_prediction.get("source") == "retropath":
        unresolved_reviews.append({
            "review_type": "predicted_route",
            "status": "pending",
            "step_indexes": sorted(
                assignment.step_index
                for assignment in selected.step_assignments
                if assignment.reaction_id.startswith("RP2STOICH:")
            ),
            "reason": "RetroPath预测路线仍需实验验证",
            "candidate_id": solution_prediction.get("candidate_id"),
            "combination_id": solution_prediction.get("combination_id"),
            "promotion_id": solution_prediction.get("promotion_id"),
        })
    selection_status = (
        "user_selected_pending_review"
        if unresolved_reviews
        else "user_selected"
    )
    warnings = list(
        dict.fromkeys([*result.warnings, *selected.warnings])
    )

    protein_payloads: list[dict[str, Any]] = []
    for protein in selected.proteins:
        payload = {
            "accession": protein.accession,
            "protein_name": protein.protein_name,
            "organism_name": protein.organism_name,
            "reviewed": protein.reviewed,
            "cofactors": protein.cofactors,
            "transmembrane": protein.transmembrane,
            "capable_step_indexes": protein.capable_step_indexes,
            "assigned_step_indexes": protein.assigned_step_indexes,
        }
        protein_payloads.append(payload)

    assignment_payloads = [
        assignment.model_dump(mode="json")
        for assignment in selected.step_assignments
    ]

    manifest_payload = {
        "schema_version": MAIN_ENZYME_MANIFEST_SCHEMA_VERSION,
        "selection_status": selection_status,
        "selected_set_id": selected.set_id,
        "selected_solution_id": result.selected_solution_id,
        "expansion_depth": result.expansion_depth,
        "chassis_key": result.chassis_key,
        "source": {
            "artifact": (
                "main_protein_selection/main_enzyme_sets.json"
            ),
            "sets_schema_version": result.schema_version,
            "algorithm_version": result.algorithm_version,
        },
        "coverage": {
            "complete": selected.coverage_complete,
            "required_step_indexes": result.required_step_indexes,
            "covered_step_indexes": selected.covered_step_indexes,
        },
        "proteins": protein_payloads,
        "step_assignments": assignment_payloads,
        "evaluation": {
            "set_status": selected.status,
            "search_complete": result.search_complete,
            "metrics": selected.metrics.model_dump(mode="json"),
            "reasons": selected.reasons,
        },
        "unresolved_reviews": unresolved_reviews,
        "warnings": warnings,
    }
    if result.genome_ec_filter is not None:
        genome_filter = result.genome_ec_filter
        manifest_payload["genome_ec_filter"] = genome_filter.model_dump(mode="json")
        manifest_payload["coverage"].update({
            "required_step_indexes": genome_filter.required_step_indexes,
            "external_required_step_indexes": result.required_step_indexes,
            "genome_supported_step_indexes": genome_filter.supported_step_indexes,
            "covered_step_indexes": sorted(set(selected.covered_step_indexes) | set(genome_filter.supported_step_indexes)),
        })
    return manifest_payload


def _step_ranges(step_indexes: list[int]) -> str:
    if not step_indexes:
        return ""
    ordered = sorted(set(step_indexes))
    ranges: list[str] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(
            str(start) if start == previous else f"{start}-{previous}"
        )
        start = previous = value
    ranges.append(
        str(start) if start == previous else f"{start}-{previous}"
    )
    return ",".join(ranges)


def write_main_enzyme_set(config: Any) -> dict[str, Any]:
    """Validate and commit one ranked main-enzyme set to the manifest."""

    directory = _selection_dir(config)
    if not directory.is_dir():
        raise FileNotFoundError("未找到主酶组合，请先运行 main-enzyme-sets")
    try:
        with FileLock(str(directory / ".main_enzyme_calculation.lock"), timeout=0):
            with manifest_update_lock(config.manifest_output_path):
                return _write_main_enzyme_set(config)
    except Timeout as exc:
        raise ValueError("主酶计算正在运行，请稍后重试") from exc


def _write_main_enzyme_set(config: Any) -> dict[str, Any]:

    set_id = _positive_set_id(config)
    target_compound = validate_target_compound_id(config.target_name)
    manifest_path = Path(config.manifest_output_path).expanduser()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if not recorded_target:
        raise ValueError(
            "manifest 尚未写入路线，请先运行 write --solution N"
        )
    if recorded_target != target_compound:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target}，"
            f"与当前输入目标 {target_compound} 不一致"
        )

    solution_id, depth, steps = _current_solution_context(manifest)

    result = _read_sets(_sets_path(config))
    if result.expansion_depth != depth or result.selected_solution_id != solution_id:
        raise ValueError("Selected route changed; rerun main-enzyme and main-enzyme-sets")
    genome_path = getattr(config, "genome_file", UNSPECIFIED_GENOME)
    preselected = getattr(result, "preselected_reactions", {})
    fixed_records, fixed_paths = {}, {}
    if preselected:
        from src.main_protein_selection.preselected_sets import validate_preselected_snapshot
        fixed_records = validate_preselected_snapshot(result, manifest, genome_path)
        for group in preselected.values():
            for protein in group.proteins:
                fixed_paths.setdefault(protein.accession, protein.sequence_file)
        if result.status == "already_selected":
            raise ValueError("主酶已选齐，无需写入组合")
    genome_filter = validate_genome_ec_filter(
        result.genome_ec_filter, heterologous_requirements(steps), result.chassis_key, genome_path,
    )
    if genome_filter is not None:
        if genome_filter.no_external_enzymes_required and not preselected:
            raise ValueError("无需外源酶和 CDS 构建，无需选择主酶组合")
    selected = _select_set(result, set_id)
    required = {int(step["step_index"]): str(step.get("reaction_id") or "") for step in heterologous_requirements(steps)}
    external = set(genome_filter.external_required_step_indexes) if genome_filter is not None else set(required)
    fixed_steps = {index for group in preselected.values() for index in group.step_indexes}
    if set(result.required_step_indexes) != external | fixed_steps or any(
        assignment.step_index not in required or assignment.reaction_id != required[assignment.step_index]
        for assignment in selected.step_assignments
    ):
        raise ValueError("Selected route steps or reactions changed; rerun main-enzyme and main-enzyme-sets")
    payload = _manifest_payload(
        result,
        selected,
        (
            manifest.get("solution", {}).get("prediction")
            if isinstance(manifest.get("solution"), Mapping)
            and isinstance(manifest.get("solution", {}).get("prediction"), Mapping)
            else None
        ),
    )

    from src.protein_selection.download import download_protein_sequence
    from src.protein_selection.sequence_store import entry_sequence, read_cached_pair
    from src.write_manifest.main_enzyme_format import selected_sequences, sequence_file_value
    from src.write_manifest.reaction_selection import (
        compact_reaction_selection, from_step_selection, is_reaction_selection, merge_manual_selection,
    )
    sequences = selected_sequences(_selection_dir(config), result, selected)
    paths = {}
    for accession, sequence in sequences.items():
        if accession in fixed_records:
            if fixed_records[accession]["sequence"] != sequence:
                raise ValueError(f"自动候选与固定蛋白序列冲突：{accession}")
            paths[accession] = fixed_paths[accession]
            continue
        downloaded = download_protein_sequence(config, uniprot=accession)
        fasta = Path(downloaded["fasta_path"])
        metadata = read_cached_pair(fasta.parent, accession)
        if downloaded["accession"] != accession or metadata is None or entry_sequence(metadata["uniprot_entry"]) != sequence:
            raise ValueError(f"本地蛋白序列与当前选中候选序列不一致：{accession}；请核对缓存并重新运行 main-enzyme")
        paths[accession] = sequence_file_value(fasta)
    payload["review_required_step_indexes"] = sorted({
        step for review in payload["unresolved_reviews"] for step in review["step_indexes"]
    })
    for assignment in payload["step_assignments"]:
        assignment["sequence_file"] = paths[assignment["accession"]]
    current_selection = manifest.get("main_enzyme_selection")
    context_payload = payload
    payload = from_step_selection(payload, steps)
    if preselected:
        if set(payload["reaction_assignments"]) & set(preselected):
            raise ValueError("自动组合试图覆盖固定反应，请重新生成组合")
        groups = {reaction: current_selection["reaction_assignments"][reaction] for reaction in preselected}
        payload["reaction_assignments"] = {**groups, **payload["reaction_assignments"]}
        payload["review_required_step_indexes"] = sorted(
            (set(payload["review_required_step_indexes"]) - fixed_steps)
            | (set(current_selection["review_required_step_indexes"]) & fixed_steps))
        payload = merge_manual_selection(payload, None, manifest)
    else:
        payload = merge_manual_selection(payload, current_selection, manifest)
    unchanged = (
        is_reaction_selection(current_selection)
        and compact_reaction_selection(current_selection) == payload
        and not current_selection.get("_compact_details_error")
    )
    if not unchanged:
        from src.main_protein_selection.lifecycle import retained_uploaded_draft
        same_selection = _same_selection(current_selection, payload, steps, paths, context_payload)
        update_design_manifest(
            manifest_path,
            target_compound_id=target_compound,
            sections={"main_enzyme_selection": payload, **(retained_uploaded_draft(manifest) if not same_selection else {})},
            discard_sections=MAIN_ENZYME_DOWNSTREAM_SECTIONS if not same_selection else (),
        )

    proteins = {}
    for group in payload["reaction_assignments"].values():
        for protein in group["proteins"]:
            proteins.setdefault(protein["accession"], set()).update(group["step_indexes"])
    ordered = sorted(proteins, key=lambda accession: (min(proteins[accession]), accession))
    pending_review = bool(payload["review_required_step_indexes"])
    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "路径编号": solution_id,
        "主酶组合编号": selected.set_id,
        "主酶数量": len(proteins),
        "覆盖步骤": sorted({index for indexes in proteins.values() for index in indexes}),
        "主酶分配": [
            (
                f"{accession}: Step "
                f"{_step_ranges(sorted(proteins[accession]))}"
            )
            for accession in ordered
        ],
        "选择状态": (
            "已选择，仍有待复核项"
            if pending_review
            else "已选择"
        ),
        "待复核类型": [
            item["review_type"]
            for item in context_payload["unresolved_reviews"]
        ],
        "清单是否更新": not unchanged,
        "清单文件": str(manifest_path.resolve()),
    }


def _same_selection(current: Any, payload: Mapping[str, Any], steps: list[dict], paths: dict, context: dict) -> bool:
    """A storage-only conversion must not invalidate existing CDSs."""
    if not isinstance(current, Mapping) or current.get("schema_version") != MAIN_ENZYME_MANIFEST_SCHEMA_VERSION:
        return False
    if current.get("_compact_details_error"):
        return False
    if current.get("selected_solution_id") != payload["selected_solution_id"]:
        return False
    if any(current.get(key) != context.get(key) for key in ("chassis_key", "genome_ec_filter")):
        return False
    from src.write_manifest.reaction_selection import from_step_selection, is_reaction_selection
    from src.write_manifest.main_enzyme_format import resolve_sequence_file
    previous = current if is_reaction_selection(current) else from_step_selection(current, steps, paths=paths)
    def identity(value):
        return {reaction: (tuple(group["step_indexes"]), {
            protein["accession"]: resolve_sequence_file(protein["sequence_file"]) for protein in group["proteins"]})
            for reaction, group in value["reaction_assignments"].items()}
    return identity(previous) == identity(payload)


def run_write_main_enzyme_set(config: Any) -> dict[str, Any]:
    """CLI entry point for ``write --main-enzyme-set N``."""

    result = write_main_enzyme_set(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "MAIN_ENZYME_MANIFEST_SCHEMA_VERSION",
    "run_write_main_enzyme_set",
    "write_main_enzyme_set",
]
