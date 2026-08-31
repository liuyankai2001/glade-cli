from __future__ import annotations

import json
from pathlib import Path
from typing import Any


TARGET_ALREADY_AVAILABLE_STATUS = "target_already_available_in_chassis"
ROUTES_FOUND_STATUS = "routes_found"
NO_PATHWAY_FOUND_STATUS = "no_pathway_found"


def target_already_available_message(target_compound: str) -> str:
    return (
        f"目标化合物 {target_compound} 已在底盘细胞可生成集合中，"
        "无需新增合成路径。"
    )


def read_target_already_available_status(
    gap_dir: str | Path,
    target_compound: str,
) -> dict[str, Any] | None:
    path = Path(gap_dir) / "run_config.json"
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"无法读取 gap 运行状态：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"gap 运行状态根节点必须是对象：{path}")
    if payload.get("status") != TARGET_ALREADY_AVAILABLE_STATUS:
        return None
    recorded_target = str(payload.get("target") or "").strip().upper()
    expected_target = str(target_compound or "").strip().upper()
    if recorded_target and recorded_target != expected_target:
        raise ValueError(
            f"gap 运行状态目标不一致：{recorded_target} != {expected_target}"
        )
    return payload


__all__ = [
    "NO_PATHWAY_FOUND_STATUS",
    "ROUTES_FOUND_STATUS",
    "TARGET_ALREADY_AVAILABLE_STATUS",
    "read_target_already_available_status",
    "target_already_available_message",
]
