"""Read-only compact view of completed gap-search results."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.kegg_gap_analyze import (
    ELECTRON_INFERENCE_VERSION,
    gap_depth_output_dir,
)
from src.pathway_analyze.target_status import (
    NO_PATHWAY_FOUND_STATUS,
    ROUTES_FOUND_STATUS,
    TARGET_ALREADY_AVAILABLE_STATUS,
    target_already_available_message,
)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到 gap 分析结果，请先运行 gap 命令："
            f"{path}"
        )
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"gap 运行配置不是有效 JSON：{path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"gap 运行配置根节点必须是对象：{path}")
    return value


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到 gap 分析结果，请先运行 gap 命令："
            f"{path}"
        )
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _read_optional_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    return _read_csv(path)


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes"}


def _split(value: Any) -> list[str]:
    return [
        item.strip()
        for item in str(value or "").split(";")
        if item.strip()
    ]


def _risk_name(value: Any) -> str:
    return {
        "none": "无",
        "low": "低",
        "medium": "中",
        "high": "高",
    }.get(str(value or "").strip().lower(), str(value or "").strip())


def _balance_name(value: Any) -> str:
    return {
        "not_applicable": "不涉及电子载体",
        "internally_balanced": "路线内已平衡",
        "unbalanced": "存在净消耗",
        "unresolved": "载体配对待确认",
    }.get(str(value or "").strip(), str(value or "").strip())


def _electron_status_name(value: Any) -> str:
    return {
        "not_required": "不需要额外电子系统",
        "internally_balanced": "电子载体已在路线内平衡",
        "internally_balanced_carrier_check_required": (
            "电子已内部平衡，但真实载体兼容性待确认"
        ),
        "carrier_check_required": "真实载体或电子伙伴待确认",
        "external_regeneration_required": "需要外部电子载体再生系统",
    }.get(str(value or "").strip(), str(value or "").strip())


def _solution_source(row: dict[str, str]) -> str:
    return (
        "RetroPath"
        if str(row.get("solution_source") or "").strip().lower() == "retropath"
        else "KEGG"
    )


def get_gap_info(config: Any) -> dict[str, Any]:
    """Return one compact view of all KEGG and RetroPath routes at a depth."""

    depth = int(getattr(config, "depth", 0))
    if depth < 0:
        raise ValueError("depth 必须大于等于 0")
    gap_dir = gap_depth_output_dir(config.gap_output_path, depth).resolve()
    run_config_path = gap_dir / "run_config.json"
    retropath_pipeline_path = gap_dir / "retropath" / "pipeline_result.json"
    if not run_config_path.is_file() and not retropath_pipeline_path.is_file():
        raise FileNotFoundError(
            "未找到 gap 分析结果，请先运行 gap 命令；已检查："
            f"{run_config_path}；{retropath_pipeline_path}"
        )
    run_config = _read_json(run_config_path) if run_config_path.is_file() else {}
    retropath_pipeline = (
        _read_json(retropath_pipeline_path)
        if retropath_pipeline_path.is_file()
        else {}
    )
    target_id = str(getattr(config, "target_name", "") or "").strip().upper()
    kegg_status = str(run_config.get("status") or "").strip()
    retropath_status = str(retropath_pipeline.get("status") or "").strip()
    search_sources = [
        source
        for source, present in (
            ("KEGG", bool(run_config)),
            ("RetroPath", bool(retropath_pipeline)),
        )
        if present
    ]

    if not run_config and retropath_pipeline.get("ok") is False:
        detail = str(retropath_pipeline.get("detail") or "").strip()
        return {
            "运行成功": False,
            "运行状态": "RetroPath 运行失败",
            "原始状态代码": retropath_status,
            "提示": detail,
            "目标化合物ID": target_id,
            "Gap深度": depth,
            "搜索来源": search_sources,
            "候选路线数": 0,
            "警告": [detail] if detail else [],
        }

    target_already_available = (
        kegg_status == TARGET_ALREADY_AVAILABLE_STATUS
        or retropath_status
        in {
            "retropath_target_already_reachable",
            "retropath_source_in_sink",
        }
    )
    if target_already_available:
        message = str(
            run_config.get("message")
            or retropath_pipeline.get("message")
            or ""
        ).strip()
        reachable_count = _as_int(run_config.get("reachable_compound_count"))
        if not run_config:
            input_summary = retropath_pipeline.get("input_summary")
            if isinstance(input_summary, dict):
                reachable_count = _as_int(
                    input_summary.get("reachable_compound_count")
                )
        return {
            "运行成功": True,
            "运行状态": "目标化合物已在底盘细胞中",
            "原始状态代码": kegg_status or retropath_status,
            "提示": message or target_already_available_message(target_id),
            "目标化合物ID": target_id,
            "Gap深度": depth,
            "搜索来源": search_sources,
            "需要新增合成路径": False,
            "可达化合物数": reachable_count,
            "警告": [],
        }

    solutions_path = gap_dir / "solutions.csv"
    solutions = _read_optional_csv(solutions_path)
    expected_solutions = (
        kegg_status == ROUTES_FOUND_STATUS
        or _as_int(retropath_pipeline.get("candidate_count")) > 0
    )
    if expected_solutions and not solutions_path.is_file():
        raise FileNotFoundError(
            "gap 运行记录显示已找到路线，但缺少共用路线文件："
            f"{solutions_path}"
        )
    if retropath_pipeline.get("ok") is False:
        solutions = [
            row for row in solutions if _solution_source(row) != "RetroPath"
        ]

    kegg_rejected = _read_optional_csv(
        gap_dir / "rejected_reaction_routes.csv"
    )
    retropath_rejected = _read_optional_csv(
        gap_dir / "retropath" / "rejected_routes.csv"
    )
    required_electron_columns = {
        "electron_balance_status",
        "electron_system_status",
        "requires_external_electron_regeneration",
        "requires_carrier_compatibility_check",
    }
    missing_electron_columns = (
        required_electron_columns.difference(solutions[0]) if solutions else set()
    )
    kegg_solutions = [
        row for row in solutions if _solution_source(row) == "KEGG"
    ]
    retropath_solutions = [
        row for row in solutions if _solution_source(row) == "RetroPath"
    ]
    if missing_electron_columns or (
        kegg_solutions
        and run_config.get("electron_inference_version")
        != ELECTRON_INFERENCE_VERSION
    ):
        raise ValueError(
            "gap 结果使用旧版电子推断，请重新运行："
            f"gap -i <输入文件> -d {depth}"
        )

    warnings: list[str] = []
    if run_config:
        recorded_depth = _as_int(run_config.get("expansion_depth"), depth)
        if recorded_depth != depth:
            warnings.append(
                "KEGG gap 运行记录的深度与请求深度不一致"
                f"（{recorded_depth} != {depth}）"
            )
        recorded_target = str(run_config.get("target") or "").strip().upper()
        if recorded_target and recorded_target != target_id:
            warnings.append(
                "KEGG gap 运行记录的目标与当前目标不一致"
                f"（{recorded_target} != {target_id}）"
            )
    if retropath_pipeline:
        recorded_depth = _as_int(
            retropath_pipeline.get("expansion_depth"), depth
        )
        if recorded_depth != depth:
            warnings.append(
                "RetroPath 运行记录的深度与请求深度不一致"
                f"（{recorded_depth} != {depth}）"
            )
        recorded_target = str(
            retropath_pipeline.get("target_compound") or ""
        ).strip().upper()
        if recorded_target and recorded_target != target_id:
            warnings.append(
                "RetroPath 运行记录的目标与当前目标不一致"
                f"（{recorded_target} != {target_id}）"
            )
        if retropath_pipeline.get("ok") is False:
            detail = str(retropath_pipeline.get("detail") or "").strip()
            warnings.append(
                "RetroPath 本次运行失败"
                + (f"：{detail}" if detail else "")
            )
    if any(_as_bool(row.get("prediction_review_required")) for row in solutions):
        warnings.append("部分 RetroPath 预测路线尚需计量、GEM 和人工验证")

    route_summaries = [
        {
            "路线编号": _as_int(row.get("solution_id")),
            "路线来源": (
                "RetroPath预测"
                if _solution_source(row) == "RetroPath"
                else "KEGG"
            ),
            "总步骤数": _as_int(row.get("total_steps")),
            "异源步骤数": _as_int(row.get("heterologous_steps")),
            "可达起始化合物": _split(row.get("reachable_anchor_compounds")),
            "需要氧气步骤数": _as_int(row.get("oxygen_required_steps")),
            "热力学不利步骤数": _as_int(row.get("thermo_disfavored_steps")),
            "最大电子风险": _risk_name(row.get("max_electron_risk_level")),
            "电子载体平衡": _balance_name(row.get("electron_balance_status")),
            "电子系统结论": _electron_status_name(
                row.get("electron_system_status")
            ),
            "需要额外电子再生系统": _as_bool(
                row.get("requires_external_electron_regeneration")
            ),
            "需要确认真实载体兼容性": _as_bool(
                row.get("requires_carrier_compatibility_check")
            ),
            "可以推荐": _as_bool(row.get("eligible_for_recommendation")),
            "需要预测结果复核": _as_bool(
                row.get("prediction_review_required")
            ),
        }
        for row in solutions
    ]

    messages = []
    for value in (
        run_config.get("message"),
        retropath_pipeline.get("message"),
    ):
        message = str(value or "").strip()
        if message and message not in messages:
            messages.append(message)
    reachable_count = _as_int(run_config.get("reachable_compound_count"))
    if not run_config:
        input_summary = retropath_pipeline.get("input_summary")
        if isinstance(input_summary, dict):
            reachable_count = _as_int(
                input_summary.get("reachable_compound_count")
            )
    aggregate_status = ROUTES_FOUND_STATUS if solutions else NO_PATHWAY_FOUND_STATUS

    return {
        "运行成功": True,
        "运行状态": (
            "已找到候选合成路径"
            if solutions
            else "未找到候选合成路径"
        ),
        "原始状态代码": aggregate_status,
        "KEGG状态代码": kegg_status,
        "RetroPath状态代码": retropath_status,
        "提示": "；".join(messages),
        "目标化合物ID": target_id,
        "Gap深度": depth,
        "搜索来源": search_sources,
        "需要新增合成路径": True,
        "可达化合物数": reachable_count,
        "候选路线数": len(solutions),
        "KEGG候选路线数": len(kegg_solutions),
        "RetroPath候选路线数": len(retropath_solutions),
        "可推荐路线数": sum(
            _as_bool(row.get("eligible_for_recommendation"))
            for row in solutions
        ),
        "需要预测复核路线数": sum(
            _as_bool(row.get("prediction_review_required"))
            for row in solutions
        ),
        "拒绝路线数": len(kegg_rejected) + len(retropath_rejected),
        "KEGG拒绝路线数": len(kegg_rejected),
        "RetroPath拒绝路线数": len(retropath_rejected),
        "搜索模式": str(run_config.get("search_mode_used") or ""),
        "使用回退搜索": _as_bool(run_config.get("did_fallback")),
        "使用RetroPath缓存": (
            _as_bool(retropath_pipeline.get("cache_hit"))
            if retropath_pipeline
            else False
        ),
        "路线摘要": route_summaries,
        "警告": warnings,
    }


def run_gap_info(config: Any) -> dict[str, Any]:
    result = get_gap_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["get_gap_info", "run_gap_info"]
