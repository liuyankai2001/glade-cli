"""Read-only views over completed GLADE workflow artifacts."""

from __future__ import annotations

from typing import Any

from src.info_show.chassis_info import (
    format_chassis_info_zh,
    get_chassis_info,
    run_chassis_info,
)


def run_info(config: Any) -> dict[str, Any]:
    """Dispatch the unified ``info`` command."""

    if getattr(config, "show_chassis", False):
        return run_chassis_info(config)
    raise ValueError("未指定信息查看类型，请使用 --chassis")


__all__ = [
    "format_chassis_info_zh",
    "get_chassis_info",
    "run_chassis_info",
    "run_info",
]
