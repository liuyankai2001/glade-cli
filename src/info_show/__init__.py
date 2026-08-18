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


def run_info(config: Any) -> dict[str, Any]:
    """Dispatch the unified ``info`` command."""

    if getattr(config, "show_chassis", False):
        raw_depth = getattr(config, "info_depth", None)
        depth = 0 if raw_depth is None else int(raw_depth)
        if depth < 0:
            raise ValueError("depth 必须大于等于 0")
        if depth == 0:
            return run_chassis_info(config)
        return run_chassis_expansion_info(config, depth)
    raise ValueError("未指定信息查看类型；--depth 必须与 --chassis 等类型参数配合使用")


__all__ = [
    "format_chassis_info_zh",
    "get_chassis_expansion_info",
    "get_chassis_info",
    "run_chassis_expansion_info",
    "run_chassis_info",
    "run_info",
]
