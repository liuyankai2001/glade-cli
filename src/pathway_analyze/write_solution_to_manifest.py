"""Compatibility imports for internal callers of the retired ``select`` command.

The active implementation lives in :mod:`src.write_manifest.solution`.
"""

from __future__ import annotations

from copy import copy
from typing import Any

from src.write_manifest.solution import run_write_solution, write_solution


def _forward_config(config: Any) -> Any:
    forwarded_config = copy(config)
    if getattr(forwarded_config, "solution", None) is None:
        forwarded_config.solution = getattr(forwarded_config, "solution_id", None)
    return forwarded_config


def select_solution(config: Any) -> dict[str, Any]:
    """Forward a legacy internal call to the unified manifest writer."""

    return write_solution(_forward_config(config))


def run_select_solution(config: Any) -> dict[str, Any]:
    """Forward the retired handler to the unified manifest writer."""

    return run_write_solution(_forward_config(config))


__all__ = ["run_select_solution", "select_solution"]
