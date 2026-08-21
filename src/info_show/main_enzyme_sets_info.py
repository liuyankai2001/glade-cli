"""Chinese read-only views for ranked main-enzyme combinations."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Mapping

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
    read_manifest,
)
from src.main_protein_selection.models import (
    MAIN_ENZYME_SELECTION_SCHEMA_VERSION,
    MAIN_ENZYME_SETS_SCHEMA_VERSION,
    MainEnzymeSelectionResult,
    MainEnzymeSet,
    MainEnzymeSetsResult,
)
from src.main_protein_selection.provenance import solution_fingerprint
from src.info_show.main_enzyme_candidates_info import (
    CONFIDENCE_NAMES,
    DIRECTION_NAMES,
    FIT_STATUS_NAMES,
)


SET_STATUS_NAMES = {
    "complete": "可直接选择",
    "review_required": "需要复核",
}

RESULT_STATUS_NAMES = {
    "complete": "已完成",
    "review_required": "已完成，但首选组合需要复核",
    "infeasible": "没有完整组合",
    "truncated": "搜索未完成",
    "stale_input": "输入已经过期",
    "source_unavailable": "候选来源不可用",
}

SPECIFICITY_NAMES = {
    "exact": "精确匹配",
    "supported": "有证据支持",
    "unknown": "待确认",
    "conflict": "冲突",
}

CARRIER_STATUS_NAMES = {
    "not_required": "无需额外检查",
    "single_multistep_enzyme": "同一多步酶，需按整体反应复核",
    "compatibility_review_required": "电子载体兼容性待确认",
    "external_regeneration_required": "需要外部电子再生",
}

ELECTRON_STATUS_NAMES = {
    "not_required": "无需重新评估",
    "auxiliary_role_identified": "已识别辅助角色，待用户选择具体蛋白",
    "review_required": "需要重新评估",
}

AUXILIARY_ROLE_NAMES = {
    "p450_reductase": "P450还原酶（CPR）",
    "ferredoxin": "铁氧还蛋白（ferredoxin）",
    "ferredoxin_reductase": "铁氧还蛋白还原酶",
    "thioredoxin": "硫氧还蛋白（thioredoxin）",
    "thioredoxin_reductase": "硫氧还蛋白还原酶",
    "generic_electron_transfer_partner": "通用电子转移伙伴",
    "oxygenase_electron_partner": "加氧酶电子转移伙伴",
}

AUXILIARY_STATUS_NAMES = {
    "not_required": "无需辅助蛋白",
    "pending_user_selection": "待用户选择",
    "integrated": "已融合在主酶中",
    "mixed": "部分已融合，部分待选择",
    "integrated_in_main_enzyme": "已融合在主酶中",
}


def _auxiliary_requirement_info(requirement: Any) -> dict[str, Any]:
    return {
        "辅助角色": AUXILIARY_ROLE_NAMES.get(
            requirement.role, requirement.role
        ),
        "必要性": "必需" if requirement.necessity == "required" else "可能需要",
        "置信度": {
            "high": "高",
            "medium": "中",
            "low": "低",
        }.get(requirement.confidence, requirement.confidence),
        "选择状态": AUXILIARY_STATUS_NAMES.get(
            requirement.selection_status,
            requirement.selection_status,
        ),
        "相关步骤": requirement.step_indexes,
        "支持的主酶": requirement.main_enzyme_accessions,
        "载体ID": requirement.carrier_ids,
        "判定依据": requirement.evidence,
    }

_MESSAGE_TRANSLATIONS = {
    "covers every required heterologous step": (
        "覆盖全部需要主酶的异源反应步骤"
    ),
    "No complete combination was found in the selected Top-N shortlist.": (
        "当前 Top-N 候选中没有能够覆盖全部步骤的主酶组合"
    ),
    "The selected multi-step enzyme replaces decomposed carrier steps; "
    "reassess its whole-reaction electron acceptor and cofactor reoxidation "
    "before downstream design.": (
        "所选多步酶替代了拆分的电子载体步骤；进入后续设计前，"
        "需要按整体反应重新确认电子受体和辅因子再氧化方式"
    ),
    "The selected set still requires an external electron regeneration "
    "mechanism.": "该组合仍需要外部电子再生机制",
    "Electron-carrier compatibility remains unresolved across the selected "
    "proteins.": "所选蛋白之间的电子载体兼容性仍待确认",
    "No route-level electron-system reassessment is required.": (
        "无需重新评估路线层面的电子系统"
    ),
    "Additional lower-ranked combinations require review.": (
        "部分较低排名的组合仍需复核"
    ),
    "Main-enzyme candidate retrieval was unavailable; rerun main-enzyme.": (
        "主酶候选检索不可用，请重新运行 main-enzyme"
    ),
}


def _sets_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "main_protein_selection"
        / "main_enzyme_sets.json"
    )


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


def _read_selection(path: Path) -> MainEnzymeSelectionResult:
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


def _translate_message(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    translated = _MESSAGE_TRANSLATIONS.get(text)
    if translated:
        return translated

    match = re.fullmatch(r"uses (\d+) distinct protein\(s\)", text)
    if match:
        return f"使用 {match.group(1)} 个不同的主酶蛋白"

    match = re.fullmatch(
        r"([^:]+): one annotation names the whole-reaction start and end",
        text,
    )
    if match:
        return f"{match.group(1)}：单条注释明确给出了整体反应的起点和终点"

    match = re.fullmatch(
        r"([^:]+): annotations name the route start, end, and intermediates",
        text,
    )
    if match:
        return f"{match.group(1)}：注释明确给出了路线起点、终点和中间体"

    match = re.fullmatch(
        r"([^:]+): annotations explicitly name every compressed Step",
        text,
    )
    if match:
        return f"{match.group(1)}：注释明确覆盖了每个被合并的反应步骤"

    match = re.fullmatch(
        r"Search stopped at max_search_nodes=(\d+); reported sets may be "
        r"incomplete\.",
        text,
    )
    if match:
        return (
            f"搜索达到节点上限 {match.group(1)}，当前组合列表可能不完整"
        )

    match = re.match(
        r"No constructable, non-conflicting Top-N candidate for Step (.+)",
        text,
    )
    if match:
        return f"以下步骤没有可构建且无冲突的 Top-N 候选：{match.group(1)}"

    replacements = (
        ("rerun main-enzyme", "请重新运行 main-enzyme"),
        ("rerun main-enzyme-sets", "请重新运行 main-enzyme-sets"),
        ("selected route fingerprint changed", "所选路线指纹已经变化"),
        ("selected solution ID changed", "所选路径编号已经变化"),
        ("selected expansion depth changed", "所选扩展深度已经变化"),
    )
    for source, target in replacements:
        text = text.replace(source, target)
    return text


def _translated_messages(values: list[str]) -> list[str]:
    return [
        translated
        for value in values
        if (translated := _translate_message(value))
    ]


def _load_current_context(
    config: Any,
) -> tuple[MainEnzymeSetsResult, int, list[dict[str, Any]]]:
    result = _read_sets(_sets_path(config))
    manifest = read_manifest(config.manifest_output_path)
    solution_id, steps = get_solution_steps(manifest)
    current_solution_fingerprint = solution_fingerprint(solution_id, steps)
    solution = manifest.get("solution")
    expansion_depth = int(
        solution.get("expansion_depth") or 0
        if isinstance(solution, Mapping)
        else 0
    )

    if result.selected_solution_id != solution_id:
        raise ValueError(
            "主酶组合对应的路径与 manifest 当前路径不一致，"
            "请重新运行 main-enzyme-sets"
        )
    if result.expansion_depth != expansion_depth:
        raise ValueError(
            "主酶组合对应的扩展深度已经变化，请重新运行 main-enzyme-sets"
        )
    if result.solution_fingerprint != current_solution_fingerprint:
        raise ValueError(
            "manifest 中的路线已经发生变化，请重新运行 main-enzyme-sets"
        )

    # Failed runs have no selectable sets.  Their blocking reasons are still
    # useful to the user and do not require the unavailable source artifacts.
    if not result.sets:
        return result, solution_id, steps

    selection = _read_selection(_selection_path(config))
    if selection.selected_solution_id != solution_id:
        raise ValueError(
            "主酶候选对应的路径已经变化，请重新运行 main-enzyme"
        )
    if selection.solution_fingerprint != current_solution_fingerprint:
        raise ValueError(
            "主酶候选对应的路线已经变化，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )
    if selection.chassis_key != result.chassis_key:
        raise ValueError(
            "主酶候选的底盘配置已经变化，请依次重新运行 "
            "main-enzyme 和 main-enzyme-sets"
        )

    candidate_csv_path = _candidate_csv_path(config)
    if not candidate_csv_path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选详情，请重新运行：main-enzyme -i <输入文件>"
        )
    candidate_rows = read_csv(candidate_csv_path)
    current_detail_fingerprint = shortlist_decision_fingerprint(candidate_rows)
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

    electron_inference = manifest.get("electron_inference")
    if not isinstance(electron_inference, Mapping):
        electron_inference = {}
    current_pool_fingerprint = candidate_pool_fingerprint_from_rows(
        solution_fingerprint_value=current_solution_fingerprint,
        chassis_key=selection.chassis_key,
        required_steps=heterologous_requirements(steps),
        candidate_rows=candidate_rows,
        electron_inference=electron_inference,
    )
    if result.candidate_pool_fingerprint != current_pool_fingerprint:
        raise ValueError(
            "路线、候选或电子系统上下文已经变化，请重新运行 "
            "main-enzyme-sets"
        )
    return result, solution_id, steps


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
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def _protein_assignments(enzyme_set: MainEnzymeSet) -> list[str]:
    proteins = sorted(
        enzyme_set.proteins,
        key=lambda item: (
            min(item.assigned_step_indexes),
            item.accession,
        ),
    )
    return [
        f"{protein.accession}（Step {_step_ranges(protein.assigned_step_indexes)}）"
        for protein in proteins
    ]


def _compact_set(enzyme_set: MainEnzymeSet) -> dict[str, Any]:
    pending_roles = sorted({
        AUXILIARY_ROLE_NAMES.get(item.role, item.role)
        for item in enzyme_set.auxiliary_requirements
        if item.selection_status == "pending_user_selection"
    })
    return {
        "组合编号": enzyme_set.set_id,
        "状态": SET_STATUS_NAMES.get(enzyme_set.status, enzyme_set.status),
        "蛋白数量": enzyme_set.protein_count,
        "主酶分工": _protein_assignments(enzyme_set),
        "需要复核步骤": enzyme_set.review_required_step_indexes,
        "风险计数": {
            "反应匹配待复核": (
                enzyme_set.metrics.reaction_fit_risk_count
            ),
            "方向待确认": enzyme_set.metrics.direction_risk_count,
            "特异性待确认": enzyme_set.metrics.specificity_risk_count,
        },
        "电子系统": CARRIER_STATUS_NAMES.get(
            enzyme_set.metrics.carrier_compatibility_status,
            enzyme_set.metrics.carrier_compatibility_status,
        ),
        "辅助蛋白需求": AUXILIARY_STATUS_NAMES.get(
            enzyme_set.auxiliary_requirement_status,
            enzyme_set.auxiliary_requirement_status,
        ),
        "待选择辅助角色": pending_roles,
    }


def get_main_enzyme_sets_info(config: Any) -> dict[str, Any]:
    """Return a compact ranked list of all generated main-enzyme sets."""

    result, solution_id, _ = _load_current_context(config)
    return {
        "运行成功": result.ok,
        "目标化合物": str(config.target_name),
        "路径编号": solution_id,
        "组合生成状态": RESULT_STATUS_NAMES.get(
            result.status, result.status
        ),
        "搜索是否完整": result.search_complete,
        "最少主酶数量": result.minimum_protein_count,
        "组合数量": len(result.sets),
        "推荐组合编号": result.sets[0].set_id if result.sets else None,
        "主酶组合": [_compact_set(item) for item in result.sets],
        "阻塞原因": _translated_messages(result.blocking_reasons),
        "警告": _translated_messages(result.warnings),
    }


def _selected_set(config: Any) -> tuple[MainEnzymeSetsResult, MainEnzymeSet]:
    raw_set_id = getattr(config, "main_enzyme_set", None)
    if raw_set_id is None:
        raise ValueError("未指定组合编号，请使用 --main-enzyme-set N")
    try:
        set_id = int(raw_set_id)
    except (TypeError, ValueError) as exc:
        raise ValueError("主酶组合编号必须是正整数") from exc
    if set_id < 1:
        raise ValueError("主酶组合编号必须是正整数")

    result, _, _ = _load_current_context(config)
    selected = next(
        (item for item in result.sets if item.set_id == set_id),
        None,
    )
    if selected is None:
        available = [item.set_id for item in result.sets]
        if not available:
            reasons = _translated_messages(result.blocking_reasons)
            suffix = f"；原因：{'；'.join(reasons)}" if reasons else ""
            raise ValueError(f"当前没有可查看的主酶组合{suffix}")
        raise ValueError(
            f"不存在主酶组合 {set_id}；可用组合编号：{available}"
        )
    return result, selected


def get_main_enzyme_set_info(config: Any) -> dict[str, Any]:
    """Return full details for one ranked main-enzyme set."""

    result, enzyme_set = _selected_set(config)
    members = sorted(
        enzyme_set.proteins,
        key=lambda item: (
            min(item.assigned_step_indexes),
            item.accession,
        ),
    )
    metrics = enzyme_set.metrics
    return {
        "运行成功": True,
        "目标化合物": str(config.target_name),
        "路径编号": result.selected_solution_id,
        "组合编号": enzyme_set.set_id,
        "组合状态": SET_STATUS_NAMES.get(
            enzyme_set.status, enzyme_set.status
        ),
        "完整覆盖": enzyme_set.coverage_complete,
        "蛋白数量": enzyme_set.protein_count,
        "需要复核步骤": enzyme_set.review_required_step_indexes,
        "主酶成员": [
            {
                "UniProt": protein.accession,
                "蛋白名称": protein.protein_name,
                "来源物种": protein.organism_name,
                "Reviewed": protein.reviewed,
                "辅因子": protein.cofactors,
                "酶系统类型": protein.enzyme_system_types,
                "辅助蛋白需求": [
                    _auxiliary_requirement_info(item)
                    for item in protein.auxiliary_requirements
                ],
                "可催化步骤": protein.capable_step_indexes,
                "实际负责步骤": protein.assigned_step_indexes,
            }
            for protein in members
        ],
        "步骤分配": [
            {
                "步骤编号": assignment.step_index,
                "反应ID": assignment.reaction_id,
                "UniProt": assignment.accession,
                "候选排名": assignment.candidate_rank,
                "蛋白评分": assignment.protein_score,
                "底盘适配评分": assignment.host_fit_score,
                "反应匹配": FIT_STATUS_NAMES.get(
                    assignment.reaction_fit_status,
                    assignment.reaction_fit_status,
                ),
                "反应匹配分数": assignment.reaction_fit_score,
                "方向判断": DIRECTION_NAMES.get(
                    assignment.direction_verdict,
                    assignment.direction_verdict,
                ),
                "方向置信度": CONFIDENCE_NAMES.get(
                    assignment.direction_confidence,
                    assignment.direction_confidence,
                ),
                "产物特异性": SPECIFICITY_NAMES.get(
                    assignment.specificity_status,
                    assignment.specificity_status,
                ),
            }
            for assignment in enzyme_set.step_assignments
        ],
        "组合评价": {
            "来源物种数量": metrics.organism_count,
            "最低反应匹配分数": metrics.min_reaction_fit_score,
            "平均反应匹配分数": metrics.mean_reaction_fit_score,
            "最低蛋白评分": metrics.min_protein_score,
            "平均蛋白评分": metrics.mean_protein_score,
            "最低底盘适配评分": metrics.min_host_fit_score,
            "平均底盘适配评分": metrics.mean_host_fit_score,
            "Reviewed比例": metrics.reviewed_fraction,
            "反应匹配待复核数": metrics.reaction_fit_risk_count,
            "方向待确认数": metrics.direction_risk_count,
            "低方向置信度数": metrics.low_direction_confidence_count,
            "特异性待确认数": metrics.specificity_risk_count,
            "精确特异性匹配数": metrics.exact_specificity_count,
            "辅助蛋白需求状态": AUXILIARY_STATUS_NAMES.get(
                enzyme_set.auxiliary_requirement_status,
                enzyme_set.auxiliary_requirement_status,
            ),
            "待选择辅助角色数": metrics.pending_auxiliary_role_count,
        },
        "电子系统": {
            "载体兼容性": CARRIER_STATUS_NAMES.get(
                metrics.carrier_compatibility_status,
                metrics.carrier_compatibility_status,
            ),
            "是否需要重新评估": ELECTRON_STATUS_NAMES.get(
                metrics.electron_reassessment_status,
                metrics.electron_reassessment_status,
            ),
            "说明": _translate_message(enzyme_set.electron_assessment),
        },
        "辅助蛋白角色": [
            _auxiliary_requirement_info(item)
            for item in enzyme_set.auxiliary_requirements
        ],
        "推荐理由": _translated_messages(enzyme_set.reasons),
        "警告": _translated_messages(enzyme_set.warnings),
        "组合指纹": enzyme_set.set_fingerprint,
    }


def run_main_enzyme_sets_info(config: Any) -> dict[str, Any]:
    """CLI-compatible entry point for ``info --main-enzyme-sets``."""

    result = get_main_enzyme_sets_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def run_main_enzyme_set_info(config: Any) -> dict[str, Any]:
    """CLI-compatible entry point for ``info --main-enzyme-set N``."""

    result = get_main_enzyme_set_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "get_main_enzyme_set_info",
    "get_main_enzyme_sets_info",
    "run_main_enzyme_set_info",
    "run_main_enzyme_sets_info",
]
