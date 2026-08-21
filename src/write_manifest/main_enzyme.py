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
    MAIN_ENZYME_SELECTION_SCHEMA_VERSION,
    MAIN_ENZYME_SETS_SCHEMA_VERSION,
    MainEnzymeSelectionResult,
    MainEnzymeSet,
    MainEnzymeSetsResult,
)
from src.main_protein_selection.literature_activity.storage import (
    artifact_fingerprint as literature_artifact_fingerprint,
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


MAIN_ENZYME_MANIFEST_SCHEMA_VERSION = "main_enzyme_manifest_selection.v2"
LITERATURE_SCHEMA_VERSION = "literature_activity_evidence.v1"
LITERATURE_RETRIEVAL_STRATEGY = "literature_experimental_activity"

MAIN_ENZYME_DOWNSTREAM_SECTIONS = (
    "protein_selection",
    "auxiliary_protein_selection",
    "enzyme_system_selection",
    "cds_selection",
    "expression_box_selection",
    "expression_cassette_assembly",
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
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


def _literature_artifact_path(config: Any) -> Path:
    return _selection_dir(config) / "literature_activity_evidence.json"


def _read_literature_artifact(
    config: Any,
) -> tuple[Path, dict[str, Any]] | None:
    path = _literature_artifact_path(config)
    if not path.is_file():
        return None
    try:
        artifact = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"文献酶活证据文件不是有效 JSON：{path}") from exc
    if not isinstance(artifact, dict):
        raise ValueError(f"文献酶活证据文件根节点必须是对象：{path}")
    if artifact.get("schema_version") != LITERATURE_SCHEMA_VERSION:
        raise ValueError(
            "不支持的文献酶活证据格式："
            f"{artifact.get('schema_version')!r}"
        )
    evidence = artifact.get("evidence")
    if not isinstance(evidence, list) or any(
        not isinstance(item, dict) for item in evidence
    ):
        raise ValueError("文献酶活证据文件的 evidence 必须是对象数组")

    stored_fingerprint = str(
        artifact.get("artifact_fingerprint") or ""
    ).strip()
    computed_fingerprint = literature_artifact_fingerprint(artifact)
    if not stored_fingerprint or stored_fingerprint != computed_fingerprint:
        raise ValueError(
            "文献酶活证据文件指纹校验失败，请重新运行 main-enzyme"
        )
    return path, artifact


def _publication_ids(record: Mapping[str, Any]) -> list[str]:
    values: list[str] = []
    for field, prefix in (
        ("doi", "DOI"),
        ("pmid", "PubMed"),
        ("pmcid", "PMC"),
    ):
        value = str(record.get(field) or "").strip()
        if not value:
            continue
        if value.upper().startswith(f"{prefix}:".upper()):
            values.append(value)
        else:
            values.append(f"{prefix}:{value}")
    return values


def _selected_literature_provenance(
    config: Any,
    selected: MainEnzymeSet,
    selection: MainEnzymeSelectionResult,
) -> dict[str, Any] | None:
    selected_candidates: list[tuple[Any, Any]] = []
    for assignment in selected.step_assignments:
        matches = [
            candidate
            for candidate in selection.candidates_by_step.get(
                assignment.step_index,
                [],
            )
            if candidate.accession == assignment.accession
            and candidate.candidate_rank == assignment.candidate_rank
            and candidate.reaction_id == assignment.reaction_id
        ]
        if len(matches) != 1:
            raise ValueError(
                f"主酶候选结果中 Step {assignment.step_index}、"
                f"{assignment.accession}、排名 {assignment.candidate_rank} "
                "的记录不唯一，请重新运行 main-enzyme"
            )
        candidate = matches[0]
        if LITERATURE_RETRIEVAL_STRATEGY in candidate.retrieval_strategies:
            selected_candidates.append((assignment, candidate))
    if not selected_candidates:
        return None

    loaded = _read_literature_artifact(config)
    if loaded is None:
        raise FileNotFoundError(
            "选中的主酶使用了文献非标准酶活证据，但缺少 "
            f"{_literature_artifact_path(config)}；请重新运行 main-enzyme"
        )
    _, artifact = loaded
    evidence_records = list(artifact.get("evidence") or [])
    evidence_by_id = {
        str(item.get("evidence_id") or "").strip(): item
        for item in evidence_records
        if str(item.get("evidence_id") or "").strip()
    }

    references: list[dict[str, Any]] = []
    for assignment, candidate in selected_candidates:
        expected_artifact_reason = (
            "literature_artifact_sha256:"
            f"{artifact['artifact_fingerprint']}"
        )
        if expected_artifact_reason not in candidate.reasons:
            raise ValueError(
                "文献酶活证据文件与主酶候选绑定的指纹不一致；"
                "请重新运行 main-enzyme 和 main-enzyme-sets"
            )
        evidence_ids = list(candidate.retrieval_query_ids)
        if not evidence_ids:
            raise ValueError(
                f"文献酶活证据文件中缺少 Step {assignment.step_index} "
                f"候选 {assignment.accession} 的关联证据"
            )
        missing_evidence_ids = sorted(
            set(evidence_ids) - set(evidence_by_id)
        )
        if missing_evidence_ids:
            raise ValueError(
                "文献酶活证据文件缺少主酶候选引用的证据："
                f"{missing_evidence_ids}；请重新运行 main-enzyme"
            )

        assignment_levels: list[str] = []
        for evidence_id in sorted(set(evidence_ids)):
            record = evidence_by_id[evidence_id]
            identity = (
                int(record.get("step_index") or 0),
                str(record.get("reaction_id") or "").strip().upper(),
                str(record.get("resolved_accession") or "").strip().upper(),
            )
            expected = (
                assignment.step_index,
                assignment.reaction_id.upper(),
                assignment.accession.upper(),
            )
            if identity != expected:
                raise ValueError(
                    f"文献证据 {evidence_id} 与选中主酶的身份不一致"
                )
            level = str(record.get("evidence_level") or "").strip()
            fit_status = str(record.get("fit_status") or "").strip()
            review_status = str(
                record.get("review_status") or ""
            ).strip()
            if level not in {"A", "B"} or fit_status != "verified_with_risk":
                raise ValueError(
                    f"文献证据 {evidence_id} 不满足进入主酶组合的条件"
                )
            if review_status != "pending":
                raise ValueError(
                    f"文献证据 {evidence_id} 的审核状态不允许写入"
                )
            assignment_levels.append(level)

            publications = _publication_ids(record)
            references.append({
                "evidence_id": evidence_id,
                "step_index": assignment.step_index,
                "reaction_id": assignment.reaction_id,
                "accession": assignment.accession,
                "evidence_level": level,
                "assay_type": str(record.get("assay_type") or ""),
                "fit_status": fit_status,
                "review_status": review_status,
                "doi": str(record.get("doi") or ""),
                "pmid": str(record.get("pmid") or ""),
                "pmcid": str(record.get("pmcid") or ""),
                "publication_ids": publications,
                "limitations": list(record.get("limitations") or []),
            })

        highest_level = "A" if "A" in assignment_levels else "B"
        expected_confidence = f"literature_grade_{highest_level.lower()}"
        if candidate.reaction_confidence != expected_confidence:
            raise ValueError(
                f"Step {assignment.step_index} 的最高文献证据等级"
                "与主酶候选不一致；请重新运行 main-enzyme 和 "
                "main-enzyme-sets"
            )

    return {
        "artifact": (
            "main_protein_selection/literature_activity_evidence.json"
        ),
        "artifact_schema_version": artifact["schema_version"],
        "algorithm_version": artifact.get("algorithm_version"),
        "artifact_fingerprint": artifact["artifact_fingerprint"],
        "review_status": "pending",
        "evidence": references,
    }


def _read_sets(path: Path) -> MainEnzymeSetsResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶组合，请先运行：main-enzyme-sets -i <输入文件>"
        )
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if payload.get("schema_version") != MAIN_ENZYME_SETS_SCHEMA_VERSION:
        raise ValueError(
            "主酶组合使用旧版辅助角色格式，请重新运行 main-enzyme-sets"
        )
    try:
        return MainEnzymeSetsResult.model_validate(payload)
    except ValueError as exc:
        raise ValueError(f"主酶组合结果格式无效：{path}") from exc


def _read_candidate_selection(path: Path) -> MainEnzymeSelectionResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选，请重新运行：main-enzyme -i <输入文件>"
        )
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if payload.get("schema_version") != MAIN_ENZYME_SELECTION_SCHEMA_VERSION:
        raise ValueError(
            "主酶候选使用旧版辅助角色格式，请重新运行 main-enzyme"
        )
    try:
        return MainEnzymeSelectionResult.model_validate(payload)
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
    literature_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    unresolved_reviews = _review_items(result, selected)
    literature_references = (
        list(literature_provenance.get("evidence") or [])
        if literature_provenance
        else []
    )
    pending_literature_steps = sorted({
        int(item["step_index"])
        for item in literature_references
        if item.get("review_status") == "pending"
    })
    if pending_literature_steps:
        unresolved_reviews.append({
            "review_type": "literature_activity",
            "status": "pending",
            "step_indexes": pending_literature_steps,
            "reason": "选中的非标准酶活来自论文抽取，仍需人工复核",
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
            "sequence_sha256": protein.sequence_sha256,
            "cofactors": protein.cofactors,
            "capable_step_indexes": protein.capable_step_indexes,
            "assigned_step_indexes": protein.assigned_step_indexes,
            "enzyme_system_types": protein.enzyme_system_types,
            "auxiliary_requirements": [
                item.model_dump(mode="json")
                for item in protein.auxiliary_requirements
            ],
        }
        evidence_ids = sorted({
            str(item["evidence_id"])
            for item in literature_references
            if item.get("accession") == protein.accession
        })
        if evidence_ids:
            payload["literature_evidence_ids"] = evidence_ids
        protein_payloads.append(payload)

    assignment_payloads: list[dict[str, Any]] = []
    for assignment in selected.step_assignments:
        payload = assignment.model_dump(mode="json")
        related = [
            item
            for item in literature_references
            if item.get("step_index") == assignment.step_index
            and item.get("accession") == assignment.accession
        ]
        if related:
            payload.update({
                "literature_evidence_ids": sorted({
                    str(item["evidence_id"]) for item in related
                }),
                "literature_publication_ids": sorted({
                    publication_id
                    for item in related
                    for publication_id in item.get("publication_ids", [])
                }),
                "literature_review_status": "pending",
            })
        assignment_payloads.append(payload)

    manifest_payload = {
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
        "proteins": protein_payloads,
        "step_assignments": assignment_payloads,
        "auxiliary_requirement_status": (
            selected.auxiliary_requirement_status
        ),
        "auxiliary_requirements": [
            item.model_dump(mode="json")
            for item in selected.auxiliary_requirements
        ],
        "evaluation": {
            "set_status": selected.status,
            "search_complete": result.search_complete,
            "metrics": selected.metrics.model_dump(mode="json"),
            "electron_assessment": selected.electron_assessment,
            "auxiliary_requirement_status": (
                selected.auxiliary_requirement_status
            ),
            "reasons": selected.reasons,
        },
        "unresolved_reviews": unresolved_reviews,
        "warnings": warnings,
    }
    if literature_provenance is not None:
        manifest_payload["literature_activity"] = literature_provenance
        manifest_payload["source"].update({
            "literature_activity_artifact": literature_provenance["artifact"],
            "literature_activity_artifact_fingerprint": (
                literature_provenance["artifact_fingerprint"]
            ),
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
    literature_provenance = _selected_literature_provenance(
        config,
        selected,
        selection,
    )
    payload = _manifest_payload(
        result,
        selected,
        literature_provenance,
    )

    current_selection = manifest.get("main_enzyme_selection")
    unchanged = (
        isinstance(current_selection, Mapping)
        and current_selection.get("schema_version")
        == MAIN_ENZYME_MANIFEST_SCHEMA_VERSION
        and current_selection.get("selected_set_fingerprint")
        == selected.set_fingerprint
        and current_selection.get("literature_activity")
        == payload.get("literature_activity")
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
    pending_auxiliary = [
        item
        for item in payload["auxiliary_requirements"]
        if item["selection_status"] == "pending_user_selection"
    ]
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
        "辅助蛋白需求状态": payload["auxiliary_requirement_status"],
        "待用户选择的辅助角色": sorted({
            item["role"] for item in pending_auxiliary
        }),
        "待选择辅助角色数量": len(pending_auxiliary),
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
