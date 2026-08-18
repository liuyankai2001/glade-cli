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


def get_gap_info(config: Any) -> dict[str, Any]:
    """Return a compact Chinese summary for one gap-search depth."""

    depth = int(getattr(config, "depth", 0))
    if depth < 0:
        raise ValueError("depth 必须大于等于 0")
    gap_dir = gap_depth_output_dir(config.gap_output_path, depth).resolve()
    run_config = _read_json(gap_dir / "run_config.json")
    solutions = _read_csv(gap_dir / "solutions.csv")
    rejected = _read_csv(gap_dir / "rejected_reaction_routes.csv")
    required_electron_columns = {
        "electron_balance_status",
        "electron_system_status",
        "requires_external_electron_regeneration",
        "requires_carrier_compatibility_check",
    }
    missing_electron_columns = (
        required_electron_columns.difference(solutions[0]) if solutions else set()
    )
    if (
        run_config.get("electron_inference_version")
        != ELECTRON_INFERENCE_VERSION
        or missing_electron_columns
    ):
        raise ValueError(
            "gap 结果使用旧版电子推断，请重新运行："
            f"gap -i <输入文件> -d {depth}"
        )

    target_id = str(getattr(config, "target_name", "") or "").strip().upper()
    warnings: list[str] = []
    recorded_depth = _as_int(run_config.get("expansion_depth"), depth)
    if recorded_depth != depth:
        warnings.append(
            f"gap 运行记录的深度与请求深度不一致（{recorded_depth} != {depth}）"
        )
    recorded_target = str(run_config.get("target") or "").strip().upper()
    if recorded_target and recorded_target != target_id:
        warnings.append(
            f"gap 运行记录的目标与当前目标不一致（{recorded_target} != {target_id}）"
        )

    route_summaries = [
        {
            "路线编号": _as_int(row.get("solution_id")),
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
        }
        for row in solutions
    ]

    return {
        "运行成功": True,
        "目标化合物ID": target_id,
        "Gap深度": depth,
        "可达化合物数": _as_int(run_config.get("reachable_compound_count")),
        "候选路线数": len(solutions),
        "可推荐路线数": sum(
            _as_bool(row.get("eligible_for_recommendation"))
            for row in solutions
        ),
        "拒绝路线数": len(rejected),
        "搜索模式": str(run_config.get("search_mode_used") or ""),
        "使用回退搜索": _as_bool(run_config.get("did_fallback")),
        "路线摘要": route_summaries,
        "警告": warnings,
    }


def run_gap_info(config: Any) -> dict[str, Any]:
    result = get_gap_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["get_gap_info", "run_gap_info"]
