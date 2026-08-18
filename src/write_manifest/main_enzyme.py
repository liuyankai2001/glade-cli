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

from src.main_protein_selection import (
    candidate_pool_fingerprint_from_rows,
    main_enzyme_selection_fingerprint,
)
from src.main_protein_selection.build_main_enzyme_sets import (
    shortlist_decision_fingerprint,
)
from src.main_protein_selection.common import (
    get_solution_steps,
    heterologous_requirements,
    read_csv,
)
from src.main_protein_selection.models import (
    MainEnzymeSelectionResult,
    MainEnzymeSet,
    MainEnzymeSetsResult,
)
from src.main_protein_selection.provenance import (
    solution_fingerprint,
    stable_json_hash,
)
from src.pathway_analyze.target_id import validate_target_compound_id
from src.write_manifest.store import (
    read_design_manifest,
    update_design_manifest,
)


MAIN_ENZYME_MANIFEST_SCHEMA_VERSION = "main_enzyme_manifest_selection.v1"

MAIN_ENZYME_DOWNSTREAM_SECTIONS = (
    "protein_selection",
    "auxiliary_protein_selection",
    "enzyme_system_selection",
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


def _candidate_csv_path(config: Any) -> Path:
    return _selection_dir(config) / "step_main_enzyme_candidates.csv"


def _read_sets(path: Path) -> MainEnzymeSetsResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶组合，请先运行：main-enzyme-sets -i <输入文件>"
        )
    try:
        return MainEnzymeSetsResult.model_validate_json(
            path.read_text(encoding="utf-8-sig")
        )
    except ValueError as exc:
        raise ValueError(f"主酶组合结果格式无效：{path}") from exc


def _read_candidate_selection(path: Path) -> MainEnzymeSelectionResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选，请重新运行：main-enzyme -i <输入文件>"
        )
    try:
        return MainEnzymeSelectionResult.model_validate_json(
            path.read_text(encoding="utf-8-sig")
        )
    except ValueError as exc:
        raise ValueError(f"主酶候选结果格式无效：{path}") from exc


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
) -> tuple[int, int, list[dict[str, Any]], str]:
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
        solution_fingerprint(solution_id, steps),
    )


def _validate_route_source(
    result: MainEnzymeSetsResult,
    selection: MainEnzymeSelectionResult,
    *,
    solution_id: int,
    expansion_depth: int,
    current_solution_fingerprint: str,
) -> None:
    if result.selected_solution_id != solution_id:
        raise ValueError(
            "主酶组合对应的路径与 manifest 当前路径不一致，"
            "请重新运行 main-enzyme-sets"
        )
    if result.expansion_depth != expansion_depth:
        raise ValueError(
            "主酶组合对应的扩展深度已经变化，请重新运行 "
            "main-enzyme-sets"
        )
    if result.solution_fingerprint != current_solution_fingerprint:
        raise ValueError(
            "manifest 中的路线已经发生变化，请重新运行 "
            "main-enzyme-sets"
        )

    if selection.selected_solution_id != solution_id:
        raise ValueError(
            "主酶候选对应的路径已经变化，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )
    if selection.expansion_depth != expansion_depth:
        raise ValueError(
            "主酶候选对应的扩展深度已经变化，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )
    if selection.solution_fingerprint != current_solution_fingerprint:
        raise ValueError(
            "主酶候选对应的路线已经变化，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )
    if selection.chassis_key != result.chassis_key:
        raise ValueError(
            "主酶候选和主酶组合使用了不同底盘，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )


def _validate_candidate_source(
    result: MainEnzymeSetsResult,
    selection: MainEnzymeSelectionResult,
    candidate_rows: list[dict[str, Any]],
    *,
    current_solution_fingerprint: str,
    required_steps: list[dict[str, Any]],
    electron_inference: Mapping[str, Any],
) -> None:
    current_detail_fingerprint = shortlist_decision_fingerprint(
        candidate_rows
    )
    if not selection.shortlist_decision_fingerprint:
        raise ValueError(
            "主酶候选缺少组合校验指纹，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )
    if (
        selection.shortlist_decision_fingerprint
        != current_detail_fingerprint
    ):
        raise ValueError(
            "主酶候选详情已经变化，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )

    stored_selection_fingerprint = result.source_artifacts.get(
        "candidate_selection_fingerprint"
    )
    current_selection_fingerprint = main_enzyme_selection_fingerprint(
        selection
    )
    if (
        not stored_selection_fingerprint
        or stored_selection_fingerprint != current_selection_fingerprint
    ):
        raise ValueError(
            "主酶候选结果已经变化，请重新运行 main-enzyme-sets"
        )

    current_pool_fingerprint = candidate_pool_fingerprint_from_rows(
        solution_fingerprint_value=current_solution_fingerprint,
        chassis_key=selection.chassis_key,
        required_steps=required_steps,
        candidate_rows=candidate_rows,
        electron_inference=electron_inference,
    )
    if result.candidate_pool_fingerprint != current_pool_fingerprint:
        raise ValueError(
            "路线、候选或电子系统上下文已经变化，请重新运行 "
            "main-enzyme-sets"
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


def _validate_selected_set(
    result: MainEnzymeSetsResult,
    selected: MainEnzymeSet,
    required_steps: list[dict[str, Any]],
) -> None:
    if not selected.coverage_complete or selected.uncovered_step_indexes:
        raise ValueError("主酶组合没有完整覆盖全部必需步骤，不能写入")

    required_by_index = {
        int(step["step_index"]): str(step.get("reaction_id") or "")
        for step in required_steps
    }
    required_indexes = sorted(required_by_index)
    if result.required_step_indexes != required_indexes:
        raise ValueError(
            "主酶组合的必需步骤与 manifest 当前路线不一致，请重新运行 "
            "main-enzyme-sets"
        )
    if selected.covered_step_indexes != required_indexes:
        raise ValueError("主酶组合没有精确覆盖 manifest 中的必需步骤")

    for assignment in selected.step_assignments:
        expected_reaction = required_by_index.get(assignment.step_index)
        if expected_reaction != assignment.reaction_id:
            raise ValueError(
                f"Step {assignment.step_index} 的反应已发生变化："
                f"组合记录为 {assignment.reaction_id}，"
                f"manifest 当前为 {expected_reaction or '不存在'}"
            )

    expected_fingerprint = stable_json_hash(
        {
            "candidate_pool_fingerprint": result.candidate_pool_fingerprint,
            "accessions": sorted(
                protein.accession for protein in selected.proteins
            ),
            "assignments": [
                {
                    "step_index": assignment.step_index,
                    "accession": assignment.accession,
                    "candidate_rank": assignment.candidate_rank,
                }
                for assignment in selected.step_assignments
            ],
        }
    )
    if selected.set_fingerprint != expected_fingerprint:
        raise ValueError(
            "主酶组合指纹校验失败，结果可能已损坏；请重新运行 "
            "main-enzyme-sets"
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

    if selected.metrics.electron_reassessment_status == "review_required":
        items.append(
            {
                "review_type": "electron_assessment",
                "status": "deferred",
                "step_indexes": selected.review_required_step_indexes,
                "reason": selected.electron_assessment,
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
) -> dict[str, Any]:
    unresolved_reviews = _review_items(result, selected)
    selection_status = (
        "user_selected_pending_review"
        if unresolved_reviews
        else "user_selected"
    )
    warnings = list(
        dict.fromkeys([*result.warnings, *selected.warnings])
    )

    return {
        "schema_version": MAIN_ENZYME_MANIFEST_SCHEMA_VERSION,
        "selection_status": selection_status,
        "selected_set_id": selected.set_id,
        "selected_set_fingerprint": selected.set_fingerprint,
        "selected_solution_id": result.selected_solution_id,
        "expansion_depth": result.expansion_depth,
        "chassis_key": result.chassis_key,
        "source": {
            "artifact": (
                "main_protein_selection/main_enzyme_sets.json"
            ),
            "sets_schema_version": result.schema_version,
            "algorithm_version": result.algorithm_version,
            "solution_fingerprint": result.solution_fingerprint,
            "candidate_pool_fingerprint": (
                result.candidate_pool_fingerprint
            ),
            "input_fingerprint": result.input_fingerprint,
            "candidate_selection_fingerprint": (
                result.source_artifacts.get(
                    "candidate_selection_fingerprint", ""
                )
            ),
            "candidate_detail_fingerprint": (
                result.source_artifacts.get(
                    "candidate_detail_fingerprint", ""
                )
            ),
        },
        "coverage": {
            "complete": selected.coverage_complete,
            "required_step_indexes": result.required_step_indexes,
            "covered_step_indexes": selected.covered_step_indexes,
        },
        "proteins": [
            {
                "accession": protein.accession,
                "protein_name": protein.protein_name,
                "organism_name": protein.organism_name,
                "reviewed": protein.reviewed,
                "sequence_sha256": protein.sequence_sha256,
                "cofactors": protein.cofactors,
                "capable_step_indexes": protein.capable_step_indexes,
                "assigned_step_indexes": protein.assigned_step_indexes,
            }
            for protein in selected.proteins
        ],
        "step_assignments": [
            assignment.model_dump(mode="json")
            for assignment in selected.step_assignments
        ],
        "evaluation": {
            "set_status": selected.status,
            "search_complete": result.search_complete,
            "metrics": selected.metrics.model_dump(mode="json"),
            "electron_assessment": selected.electron_assessment,
            "reasons": selected.reasons,
        },
        "unresolved_reviews": unresolved_reviews,
        "warnings": warnings,
    }


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

    (
        solution_id,
        expansion_depth,
        steps,
        current_solution_fingerprint,
    ) = _current_solution_context(manifest)

    result = _read_sets(_sets_path(config))
    selection = _read_candidate_selection(
        _candidate_selection_path(config)
    )
    _validate_route_source(
        result,
        selection,
        solution_id=solution_id,
        expansion_depth=expansion_depth,
        current_solution_fingerprint=current_solution_fingerprint,
    )

    candidate_path = _candidate_csv_path(config)
    if not candidate_path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选详情，请重新运行："
            "main-enzyme -i <输入文件>"
        )
    candidate_rows = read_csv(candidate_path)
    required_steps = heterologous_requirements(steps)
    electron_inference = manifest.get("electron_inference")
    if not isinstance(electron_inference, Mapping):
        electron_inference = {}
    _validate_candidate_source(
        result,
        selection,
        candidate_rows,
        current_solution_fingerprint=current_solution_fingerprint,
        required_steps=required_steps,
        electron_inference=electron_inference,
    )

    selected = _select_set(result, set_id)
    _validate_selected_set(result, selected, required_steps)
    payload = _manifest_payload(result, selected)

    current_selection = manifest.get("main_enzyme_selection")
    unchanged = (
        isinstance(current_selection, Mapping)
        and current_selection.get("schema_version")
        == MAIN_ENZYME_MANIFEST_SCHEMA_VERSION
        and current_selection.get("selected_set_fingerprint")
        == selected.set_fingerprint
    )
    if unchanged:
        updated_manifest = manifest
    else:
        updated_manifest = update_design_manifest(
            manifest_path,
            target_compound_id=target_compound,
            sections={"main_enzyme_selection": payload},
            discard_sections=MAIN_ENZYME_DOWNSTREAM_SECTIONS,
        )

    proteins = sorted(
        selected.proteins,
        key=lambda item: (
            min(item.assigned_step_indexes),
            item.accession,
        ),
    )
    pending_review = bool(payload["unresolved_reviews"])
    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "路径编号": solution_id,
        "主酶组合编号": selected.set_id,
        "主酶数量": selected.protein_count,
        "覆盖步骤": selected.covered_step_indexes,
        "主酶分配": [
            (
                f"{protein.accession}: Step "
                f"{_step_ranges(protein.assigned_step_indexes)}"
            )
            for protein in proteins
        ],
        "选择状态": (
            "已选择，仍有待复核项"
            if pending_review
            else "已选择"
        ),
        "待复核类型": [
            item["review_type"]
            for item in payload["unresolved_reviews"]
        ],
        "清单是否更新": not unchanged,
        "清单文件": str(manifest_path.resolve()),
        "清单版本": updated_manifest["revision"],
    }


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
