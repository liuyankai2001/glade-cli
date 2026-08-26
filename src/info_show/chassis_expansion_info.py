"""Read-only view of one completed chassis-expansion depth."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.expand_chassis_metabolites import (
    EXPANSION_RULE_VERSION,
)


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到底盘扩展结果，请先运行 expand 命令："
            f"{path}"
        )
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"底盘扩展 manifest 不是有效 JSON：{path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"底盘扩展 manifest 根节点必须是对象：{path}")
    return value


def _read_reachable_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"未找到指定深度的累计可达化合物表：{path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"depth", "kegg_id"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                f"累计可达化合物表必须包含 depth 和 kegg_id 列：{path}"
            )
        return [dict(row) for row in reader]


def _as_int(value: Any, field: str) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} 不是有效整数：{value!r}") from exc


def _as_text_list(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    return [
        item.strip()
        for item in str(value or "").split(";")
        if item.strip()
    ]


def get_chassis_expansion_info(config: Any, depth: int) -> dict[str, Any]:
    """Return a compact Chinese summary for one existing expansion depth."""

    requested_depth = int(depth)
    if requested_depth < 1:
        raise ValueError("扩展深度必须大于等于 1；原始底盘请使用 --chassis")

    output_dir = Path(config.chassis_output_path).expanduser().resolve()
    manifest_path = output_dir / "chassis_expansion_manifest.json"
    manifest = _read_manifest(manifest_path)
    if manifest.get("algorithm_version") != EXPANSION_RULE_VERSION:
        raise ValueError(
            "底盘扩展结果使用旧版扩展策略，请重新运行 "
            f"expand -d {requested_depth}"
        )
    max_depth = _as_int(manifest.get("max_depth", 0), "max_depth")
    if requested_depth > max_depth:
        raise ValueError(
            f"尚未生成深度 {requested_depth} 的底盘扩展结果；"
            f"当前最大深度为 {max_depth}，请先运行 expand -d {requested_depth}"
        )

    layers = manifest.get("layers")
    if not isinstance(layers, dict):
        raise ValueError(f"底盘扩展 manifest 缺少 layers：{manifest_path}")
    layer = layers.get(str(requested_depth))
    if not isinstance(layer, dict):
        raise ValueError(
            f"底盘扩展 manifest 缺少深度 {requested_depth} 的记录：{manifest_path}"
        )

    expanded_name = Path(str(layer.get("expanded_file") or "")).name
    if not expanded_name:
        raise ValueError(f"深度 {requested_depth} 缺少 expanded_file 记录")
    expanded_path = output_dir / expanded_name
    rows = _read_reachable_rows(expanded_path)

    compounds_by_depth: dict[int, set[str]] = {}
    first_depth_by_compound: dict[str, int] = {}
    for row in rows:
        kegg_id = str(row.get("kegg_id") or "").strip().upper()
        if not kegg_id:
            continue
        row_depth = _as_int(row.get("depth", 0), "depth")
        compounds_by_depth.setdefault(row_depth, set()).add(kegg_id)
        first_depth_by_compound[kegg_id] = min(
            row_depth,
            first_depth_by_compound.get(kegg_id, row_depth),
        )

    observed_base_count = len(compounds_by_depth.get(0, set()))
    observed_frontier_count = len(compounds_by_depth.get(requested_depth, set()))
    observed_cumulative_count = len(first_depth_by_compound)
    first_layer = layers.get("1")
    base_count = _as_int(
        first_layer.get("input_frontier_count", observed_base_count)
        if isinstance(first_layer, dict)
        else observed_base_count,
        "input_frontier_count",
    )
    frontier_count = _as_int(
        layer.get("frontier_compound_count", observed_frontier_count),
        "frontier_compound_count",
    )
    cumulative_count = _as_int(
        layer.get("cumulative_compound_count", observed_cumulative_count),
        "cumulative_compound_count",
    )
    target_id = str(getattr(config, "target_name", "") or "").strip().upper()
    target_first_depth = first_depth_by_compound.get(target_id)

    warnings: list[str] = []
    missing_reactions = _as_int(
        layer.get("missing_reaction_count", 0),
        "missing_reaction_count",
    )
    if missing_reactions:
        warnings.append(f"本层有 {missing_reactions} 个 KEGG 反应未能加载")
    carrier_review_count = _as_int(
        layer.get("carrier_review_required_compound_count", 0),
        "carrier_review_required_compound_count",
    )
    if carrier_review_count:
        warnings.append(
            f"本层有 {carrier_review_count} 个新增化合物依赖电子载体，"
            "需要辅助系统与载体兼容性复核"
        )
    return {
        "运行成功": True,
        "目标化合物ID": target_id,
        "查看扩展深度": requested_depth,
        "目标在该深度内可达": (
            target_first_depth is not None
            and target_first_depth <= requested_depth
        ),
        "目标首次可达深度": target_first_depth,
        "初始可达化合物数": base_count,
        "本层新增化合物数": frontier_count,
        "累计可达化合物数": cumulative_count,
        "本层候选反应数": _as_int(
            layer.get("candidate_reaction_count", 0),
            "candidate_reaction_count",
        ),
        "本层可用反应数": _as_int(
            layer.get("eligible_reaction_count", 0),
            "eligible_reaction_count",
        ),
        "本层缺失反应数": missing_reactions,
        "本层载体支持反应数": _as_int(
            layer.get("carrier_supported_reaction_count", 0),
            "carrier_supported_reaction_count",
        ),
        "本层载体支持化合物数": _as_int(
            layer.get("carrier_supported_compound_count", 0),
            "carrier_supported_compound_count",
        ),
        "本层载体复核化合物数": carrier_review_count,
        "本层所需辅助角色": _as_text_list(
            layer.get("required_auxiliary_roles")
        ),
        "本层电子载体": _as_text_list(layer.get("electron_carrier_ids")),
        "警告": warnings,
    }


def run_chassis_expansion_info(config: Any, depth: int) -> dict[str, Any]:
    result = get_chassis_expansion_info(config, depth)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["get_chassis_expansion_info", "run_chassis_expansion_info"]
