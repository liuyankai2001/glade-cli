"""Chinese read-only views for ranked main-enzyme combinations."""

from __future__ import annotations

import json
from src.main_protein_selection.transmembrane import display as membrane_display
import re
from pathlib import Path
from typing import Any

from src.main_protein_selection.common import (
    get_solution_steps,
    heterologous_requirements,
    read_manifest,
)
from src.main_protein_selection.genome_ec_filter import validate_genome_ec_filter
from src.main_protein_selection.models import (
    SUPPORTED_MAIN_ENZYME_SETS_SCHEMAS,
    parse_main_enzyme_sets,
    MainEnzymeSet,
    MainEnzymeSetsResult,
)
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
    "already_selected": "主酶已选齐，无需生成或写入组合",
}

SPECIFICITY_NAMES = {
    "exact": "精确匹配",
    "supported": "有证据支持",
    "unknown": "待确认",
    "conflict": "冲突",
}


_MESSAGE_TRANSLATIONS = {
    "includes every preselected protein": "包含人工涉及反应下的全部固定蛋白",
    "covers every required heterologous step": (
        "覆盖全部需要主酶的异源反应步骤"
    ),
    "No complete combination was found in the selected Top-N shortlist.": (
        "当前 Top-N 候选中没有能够覆盖全部步骤的主酶组合"
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


def _read_sets(path: Path) -> MainEnzymeSetsResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶组合，请先运行：main-enzyme-sets -i <输入文件>"
        )
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError("主酶组合结果必须是 JSON 对象，请重新运行 main-enzyme-sets")
    if payload.get("schema_version") not in SUPPORTED_MAIN_ENZYME_SETS_SCHEMAS:
        raise ValueError(
            "主酶组合使用旧版结果格式，请重新运行 main-enzyme-sets"
        )
    try:
        return parse_main_enzyme_sets(payload)
    except ValueError as exc:
        raise ValueError(f"主酶组合结果格式无效：{path}") from exc


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
    if getattr(result, "preselected_reactions", None):
        from src.main_protein_selection.preselected_sets import validate_preselected_snapshot
        validate_preselected_snapshot(result, manifest, getattr(config, "genome_file", None))
    if result.genome_ec_filter is not None:
        validate_genome_ec_filter(result.genome_ec_filter, heterologous_requirements(steps),
                                 result.chassis_key, getattr(config, "genome_file", None))
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
    return {
        "组合编号": enzyme_set.set_id,
        "状态": SET_STATUS_NAMES.get(enzyme_set.status, enzyme_set.status),
        "蛋白数量": enzyme_set.protein_count,
        "主酶分工": _protein_assignments(enzyme_set),
        "跨膜信息": {protein.accession: membrane_display(protein.transmembrane)
                     for protein in enzyme_set.proteins},
        "需要复核步骤": enzyme_set.review_required_step_indexes,
        "风险计数": {
            "反应匹配待复核": (
                enzyme_set.metrics.reaction_fit_risk_count
            ),
            "方向待确认": enzyme_set.metrics.direction_risk_count,
            "特异性待确认": enzyme_set.metrics.specificity_risk_count,
        },
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
        ) if not (result.genome_ec_filter and result.genome_ec_filter.no_external_enzymes_required
                  and not getattr(result, "preselected_reactions", None))
        else "无需外源酶和 CDS 构建",
        "基因组支持步骤": result.genome_ec_filter.supported_step_indexes if result.genome_ec_filter else [],
        "需外源酶步骤": result.required_step_indexes,
        "搜索是否完整": result.search_complete,
        "最少主酶数量": result.minimum_protein_count,
        "组合数量": len(result.sets),
        "推荐组合编号": result.sets[0].set_id if result.sets else None,
        "主酶组合": [_compact_set(item) for item in result.sets],
        **_fixed_details(result),
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
    if result.status == "already_selected":
        raise ValueError("主酶已选齐，无需选择或查看组合")
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
        **_fixed_details(result),
        "主酶成员": [
            {
                ("蛋白ID" if protein.reviewed is None else "UniProt"): protein.accession,
                "蛋白名称": protein.protein_name or "未知",
                "跨膜信息": membrane_display(protein.transmembrane),
                "来源物种": protein.organism_name or "未知",
                "Reviewed": protein.reviewed if protein.reviewed is not None else "未知",
                "辅因子": protein.cofactors if protein.reviewed is not None else "未知",
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
        },
        "推荐理由": _translated_messages(enzyme_set.reasons),
        "警告": _translated_messages(enzyme_set.warnings),
    }


def _fixed_details(result):
    fixed = getattr(result, "preselected_reactions", {})
    if not fixed:
        return {}
    return {"固定反应": [{"反应ID": reaction, "步骤编号": group.step_indexes,
                         "已选蛋白": [protein.accession for protein in group.proteins],
                         "跨膜信息": {protein.accession: membrane_display(protein.transmembrane)
                                      for protein in group.proteins}}
                        for reaction, group in fixed.items()],
            "评分范围": "反应、蛋白、底盘评分及风险计数仅评估自动选择部分；数量、物种及Reviewed比例按完整方案统计"}


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
