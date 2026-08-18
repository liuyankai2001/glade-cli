"""Compact read-only view of main-enzyme candidates grouped by route step."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.main_protein_selection.common import get_solution_steps, read_manifest
from src.main_protein_selection.models import MainEnzymeSelectionResult
from src.main_protein_selection.provenance import solution_fingerprint


FIT_STATUS_NAMES = {
    "verified": "已验证",
    "verified_with_risk": "已验证，但需要复核",
}

DIRECTION_NAMES = {
    "supported": "支持",
    "unknown": "待确认",
    "contradicted": "矛盾",
}

CONFIDENCE_NAMES = {
    "high": "高",
    "medium": "中",
    "low": "低",
}


def _selection_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "main_protein_selection"
        / "main_enzyme_selection.json"
    )


def _candidate_csv_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "main_protein_selection"
        / "step_main_enzyme_candidates.csv"
    )


def _read_selection(path: Path) -> MainEnzymeSelectionResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选，请先运行：main-enzyme -i <输入文件>"
        )
    try:
        return MainEnzymeSelectionResult.model_validate_json(
            path.read_text(encoding="utf-8-sig")
        )
    except ValueError as exc:
        raise ValueError(f"主酶候选结果格式无效：{path}") from exc


def _load_current_selection_context(
    config: Any,
) -> tuple[
    MainEnzymeSelectionResult,
    int,
    list[dict[str, Any]],
    dict[int, dict[str, Any]],
]:
    selection = _read_selection(_selection_path(config))
    manifest = read_manifest(config.manifest_output_path)
    solution_id, steps = get_solution_steps(manifest)
    current_fingerprint = solution_fingerprint(solution_id, steps)

    if selection.selected_solution_id != solution_id:
        raise ValueError(
            "主酶候选对应的路线与 manifest 当前路线不一致，"
            "请重新运行 main-enzyme"
        )
    if selection.solution_fingerprint != current_fingerprint:
        raise ValueError(
            "manifest 中的路线已经发生变化，请重新运行 main-enzyme"
        )

    steps_by_index = {
        int(step.get("step_index") or 0): step
        for step in steps
        if int(step.get("step_index") or 0) > 0
    }
    return selection, solution_id, steps, steps_by_index


def get_main_enzyme_candidates_info(config: Any) -> dict[str, Any]:
    """读取当前已选路线的主酶候选，并按正向步骤编号分组。"""

    raw_step = getattr(config, "step", None)
    try:
        selected_step = int(raw_step) if raw_step is not None else None
    except (TypeError, ValueError) as exc:
        raise ValueError(f"step 必须是正整数，当前值：{raw_step!r}") from exc
    if selected_step is not None and selected_step < 1:
        raise ValueError("step 必须是正整数")

    selection, solution_id, _, steps_by_index = _load_current_selection_context(
        config
    )
    if selected_step is not None and selected_step not in steps_by_index:
        raise ValueError(
            f"路线中不存在 Step {selected_step}；"
            f"可用步骤：{sorted(steps_by_index)}"
        )

    if selected_step is None:
        step_indexes = sorted(selection.candidates_by_step)
    else:
        step_indexes = [selected_step]

    grouped_candidates: list[dict[str, Any]] = []
    for step_index in step_indexes:
        step = steps_by_index.get(step_index)
        if step is None:
            raise ValueError(
                f"主酶候选引用了 manifest 中不存在的步骤：{step_index}"
            )
        candidates = sorted(
            selection.candidates_by_step.get(step_index, []),
            key=lambda candidate: candidate.candidate_rank,
        )
        candidate_note = ""
        if not candidates:
            if str(step.get("status") or "").strip() == "endogenous":
                candidate_note = "该步骤为内源反应，无需选择主酶"
            else:
                candidate_note = "该步骤没有可用的主酶候选"
        grouped_candidates.append({
            "步骤编号": step_index,
            "反应ID": step.get("reaction_id"),
            "候选数量": len(candidates),
            "说明": candidate_note,
            "候选酶": [
                {
                    "排名": candidate.candidate_rank,
                    "UniProt": candidate.accession,
                    "蛋白名称": candidate.protein_name,
                    "来源物种": candidate.organism_name,
                    "EC编号": candidate.ec_number,
                    "Reviewed": candidate.reviewed,
                    "蛋白评分": candidate.protein_score,
                    "反应匹配": FIT_STATUS_NAMES.get(
                        candidate.reaction_fit_status,
                        candidate.reaction_fit_status,
                    ),
                    "反应匹配分数": candidate.reaction_fit_score,
                    "方向判断": DIRECTION_NAMES.get(
                        candidate.direction_verdict,
                        candidate.direction_verdict,
                    ),
                    "方向置信度": CONFIDENCE_NAMES.get(
                        candidate.direction_confidence,
                        candidate.direction_confidence,
                    ),
                }
                for candidate in candidates
            ],
        })

    def filter_step_indexes(values: list[int]) -> list[int]:
        if selected_step is None:
            return values
        return [value for value in values if value == selected_step]

    result: dict[str, Any] = {
        "运行成功": True,
        "目标化合物": str(config.target_name),
        "路径编号": solution_id,
        "候选生成状态": selection.status,
    }
    if selected_step is not None:
        result["查看步骤"] = selected_step
    result.update({
        "未覆盖步骤": filter_step_indexes(selection.uncovered_step_indexes),
        "方向待确认步骤": filter_step_indexes(
            selection.direction_risk_step_indexes
        ),
        "步骤候选": grouped_candidates,
    })
    return result


def _as_int(value: Any) -> int | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _as_float(value: Any) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _as_optional_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if not text:
        return None
    return text in {"1", "true", "yes"}


def _split_values(value: Any) -> list[str]:
    return [
        item.strip()
        for item in str(value or "").split(";")
        if item.strip()
    ]


def _split_evidence(value: Any) -> list[str]:
    return [
        item.strip()
        for item in str(value or "").split(" | ")
        if item.strip()
    ]


def _json_object(value: Any) -> dict[str, Any]:
    text = str(value or "").strip()
    if not text:
        return {}
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"候选详情中的 JSON 字段无效：{text!r}") from exc
    if not isinstance(parsed, dict):
        raise ValueError(f"候选详情中的 JSON 字段必须是对象：{text!r}")
    return parsed


def _read_candidate_detail_row(
    path: Path,
    *,
    solution_id: int,
    step_index: int,
    candidate_rank: int,
) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选详情，请重新运行：main-enzyme -i <输入文件>"
        )
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        matches = [
            dict(row)
            for row in csv.DictReader(handle)
            if _as_int(row.get("solution_id")) == solution_id
            and _as_int(row.get("step_index")) == step_index
            and _as_int(row.get("candidate_rank")) == candidate_rank
        ]
    if not matches:
        raise ValueError(
            f"候选详情 CSV 中缺少 Step {step_index} 排名 {candidate_rank}；"
            "请重新运行 main-enzyme"
        )
    if len(matches) > 1:
        raise ValueError(
            f"候选详情 CSV 中 Step {step_index} 排名 {candidate_rank} 不唯一；"
            "请重新运行 main-enzyme"
        )
    return matches[0]


def get_main_enzyme_candidate_info(config: Any) -> dict[str, Any]:
    """读取指定正向步骤中某一排名主酶候选的完整信息。"""

    raw_rank = getattr(config, "main_enzyme_candidate", None)
    raw_step = getattr(config, "step", None)
    if raw_rank is None:
        raise ValueError("未指定候选排名，请使用 --main-enzyme-candidate N")
    if raw_step is None:
        raise ValueError("--main-enzyme-candidate 必须与 --step 一起使用")
    try:
        candidate_rank = int(raw_rank)
        step_index = int(raw_step)
    except (TypeError, ValueError) as exc:
        raise ValueError("候选排名和 step 必须是正整数") from exc
    if candidate_rank < 1:
        raise ValueError("主酶候选排名必须是正整数")
    if step_index < 1:
        raise ValueError("step 必须是正整数")

    selection, solution_id, _, steps_by_index = _load_current_selection_context(
        config
    )
    step = steps_by_index.get(step_index)
    if step is None:
        raise ValueError(
            f"路线中不存在 Step {step_index}；"
            f"可用步骤：{sorted(steps_by_index)}"
        )

    step_candidates = sorted(
        selection.candidates_by_step.get(step_index, []),
        key=lambda item: item.candidate_rank,
    )
    if not step_candidates:
        if str(step.get("status") or "").strip() == "endogenous":
            raise ValueError(f"Step {step_index} 为内源反应，无需选择主酶")
        raise ValueError(f"Step {step_index} 没有可用的主酶候选")

    candidate = next(
        (
            item
            for item in step_candidates
            if item.candidate_rank == candidate_rank
        ),
        None,
    )
    if candidate is None:
        available_ranks = [item.candidate_rank for item in step_candidates]
        raise ValueError(
            f"Step {step_index} 不存在排名 {candidate_rank} 的候选；"
            f"可用排名：{available_ranks}"
        )

    detail = _read_candidate_detail_row(
        _candidate_csv_path(config),
        solution_id=solution_id,
        step_index=step_index,
        candidate_rank=candidate_rank,
    )
    expected_values = {
        "accession": candidate.accession,
        "reaction_id": str(step.get("reaction_id") or ""),
    }
    for field, expected in expected_values.items():
        actual = str(detail.get(field) or "").strip()
        if actual.upper() != expected.strip().upper():
            raise ValueError(
                f"主酶候选 JSON 与 CSV 的 {field} 不一致，"
                "请重新运行 main-enzyme"
            )
    csv_sequence_sha256 = str(detail.get("sequence_sha256") or "").strip()
    if (
        candidate.sequence_sha256
        and csv_sequence_sha256
        and candidate.sequence_sha256.lower() != csv_sequence_sha256.lower()
    ):
        raise ValueError(
            "主酶候选 JSON 与 CSV 的序列摘要不一致，"
            "请重新运行 main-enzyme"
        )

    coverage: list[dict[str, Any]] = []
    for covered_step_index in sorted(selection.candidates_by_step):
        covered_step = steps_by_index.get(covered_step_index, {})
        for covered_candidate in selection.candidates_by_step[covered_step_index]:
            if covered_candidate.accession == candidate.accession:
                coverage.append({
                    "步骤编号": covered_step_index,
                    "反应ID": covered_step.get("reaction_id"),
                    "候选排名": covered_candidate.candidate_rank,
                })
                break

    return {
        "运行成功": True,
        "目标化合物": str(config.target_name),
        "路径编号": solution_id,
        "步骤编号": step_index,
        "反应ID": step.get("reaction_id"),
        "候选排名": candidate_rank,
        "基本信息": {
            "UniProt": candidate.accession,
            "Entry名称": detail.get("entry_name"),
            "蛋白名称": candidate.protein_name,
            "基因名称": _split_values(detail.get("gene_names")),
            "别名": _split_values(detail.get("aliases")),
            "来源物种": candidate.organism_name,
            "物种ID": _as_int(detail.get("organism_id")),
            "Reviewed": candidate.reviewed,
            "蛋白存在证据": detail.get("protein_existence"),
        },
        "路线覆盖": {
            "当前反应名称": detail.get("reaction_name"),
            "当前生成化合物ID": detail.get("produced_compound_id"),
            "当前生成化合物名称": detail.get("produced_compound_name"),
            "覆盖步骤": [item["步骤编号"] for item in coverage],
            "覆盖反应": [item["反应ID"] for item in coverage],
            "各步骤排名": coverage,
        },
        "评分与筛选": {
            "综合评分": candidate.protein_score,
            "分项评分": _json_object(detail.get("score_breakdown")),
            "评估排名": _as_int(detail.get("evaluation_rank")),
            "选择状态": detail.get("selection_status"),
            "反应匹配": FIT_STATUS_NAMES.get(
                candidate.reaction_fit_status,
                candidate.reaction_fit_status,
            ),
            "反应匹配分数": candidate.reaction_fit_score,
            "反应匹配规则": _split_values(
                detail.get("reaction_fit_rule_ids")
            ),
            "反应匹配证据": _split_evidence(
                detail.get("reaction_fit_evidence")
            ),
            "推荐理由": _split_values(detail.get("reasons")),
            "警告": _split_evidence(detail.get("warnings")),
        },
        "催化信息": {
            "EC编号": _split_values(
                detail.get("ec_numbers") or detail.get("ec_number")
            ),
            "催化活性": detail.get("catalytic_activities"),
            "辅因子": _split_values(detail.get("cofactors")),
            "Rhea编号": _split_values(detail.get("rhea_ids")),
            "匹配Rhea编号": _split_values(detail.get("matched_rhea_ids")),
            "匹配KO编号": _split_values(detail.get("matched_ko_ids")),
            "KEGG基因编号": _split_values(detail.get("kegg_gene_ids")),
        },
        "方向证据": {
            "方向判断": DIRECTION_NAMES.get(
                candidate.direction_verdict,
                candidate.direction_verdict,
            ),
            "方向置信度": CONFIDENCE_NAMES.get(
                candidate.direction_confidence,
                candidate.direction_confidence,
            ),
            "证据级别": detail.get("direction_evidence_level"),
            "证据来源": _split_values(
                detail.get("direction_evidence_source_ids")
            ),
            "证据说明": _split_evidence(detail.get("direction_evidence")),
            "所需Rhea方向编号": _split_values(
                detail.get("required_rhea_direction_ids")
            ),
        },
        "检索证据": {
            "检索策略": _split_values(detail.get("retrieval_strategy")),
            "检索查询ID": _split_values(detail.get("retrieval_query_id")),
            "反应置信度": detail.get("reaction_confidence"),
            "Selenzyme": {
                "排名": _as_int(detail.get("selenzyme_rank")),
                "评分": _as_float(detail.get("selenzyme_score")),
                "反应相似度": _as_float(
                    detail.get("selenzyme_reaction_similarity")
                ),
                "风险状态": detail.get("selenzyme_risk_status"),
                "匹配反应": detail.get("selenzyme_matched_reaction_id"),
                "分类距离": _as_float(
                    detail.get("selenzyme_taxonomic_distance")
                ),
                "使用方向": detail.get("selenzyme_direction_used"),
                "方向优选": _as_optional_bool(
                    detail.get("selenzyme_direction_preferred")
                ),
            },
        },
        "表达与结构": {
            "亚细胞定位": _split_values(detail.get("subcellular_locations")),
            "亚基信息": detail.get("subunit"),
            "功能说明": detail.get("function_comments"),
            "翻译后修饰": detail.get("ptm_comments"),
            "特征注释": detail.get("feature_annotations"),
            "结构域": _split_values(detail.get("domain_ids")),
            "关键词": _split_values(detail.get("keywords")),
        },
        "系统组成信息": {
            "角色": detail.get("role"),
            "候选角色": detail.get("candidate_role"),
            "组分类型": detail.get("component_type"),
            "系统锚点UniProt": detail.get("system_anchor_accession"),
            "ComplexPortal编号": detail.get("complex_portal_ac"),
            "ComplexPortal名称": detail.get("complex_portal_name"),
            "复合物匹配依据": detail.get("complex_match_basis"),
            "复合物证据类型": detail.get("complex_evidence_type"),
            "组分计量": _json_object(detail.get("component_stoichiometry")),
        },
        "文献与数据库": {
            "文献": _split_values(detail.get("publication_ids")),
            "交叉引用": _split_values(detail.get("cross_references")),
        },
        "序列信息": {
            "长度": candidate.length,
            "序列版本": candidate.sequence_version,
            "SHA256": candidate.sequence_sha256,
            "序列": candidate.sequence,
        },
    }


def run_main_enzyme_candidate_info(config: Any) -> dict[str, Any]:
    """CLI entry point for one ranked main-enzyme candidate detail."""

    result = get_main_enzyme_candidate_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def run_main_enzyme_candidates_info(config: Any) -> dict[str, Any]:
    """CLI-compatible entry point for ``info --main-enzyme-candidates``."""

    result = get_main_enzyme_candidates_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "get_main_enzyme_candidate_info",
    "get_main_enzyme_candidates_info",
    "run_main_enzyme_candidate_info",
    "run_main_enzyme_candidates_info",
]
