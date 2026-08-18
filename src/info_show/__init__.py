"""Read-only views over completed GLADE workflow artifacts."""

from __future__ import annotations

from typing import Any

from src.info_show.chassis_info import (
    format_chassis_info_zh,
    get_chassis_info,
    run_chassis_info,
)
from src.info_show.chassis_expansion_info import (
    get_chassis_expansion_info,
    run_chassis_expansion_info,
)
from src.info_show.gap_info import get_gap_info, run_gap_info
from src.pathway_analyze.list_solution_steps import run_solution_steps


def run_info(config: Any) -> dict[str, Any]:
    """Dispatch the unified ``info`` command."""

    if getattr(config, "chassis", False):
        raw_depth = getattr(config, "depth", None)
        depth = 0 if raw_depth is None else int(raw_depth)
        if depth < 0:
            raise ValueError("depth 必须大于等于 0")
        if depth == 0:
            return run_chassis_info(config)
        return run_chassis_expansion_info(config, depth)
    if getattr(config, "gap", False):
        return run_gap_info(config)
    if getattr(config, "solution", None) is not None:
        return run_solution_steps(config)
    raise ValueError("未指定信息查看类型，请使用 --chassis、--gap 或 --solution")


__all__ = [
    "format_chassis_info_zh",
    "get_chassis_expansion_info",
    "get_chassis_info",
    "get_gap_info",
    "run_chassis_expansion_info",
    "run_chassis_info",
    "run_gap_info",
    "run_info",
]
