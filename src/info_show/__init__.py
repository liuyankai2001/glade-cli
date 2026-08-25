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
from src.info_show.main_enzyme_candidates_info import (
    get_main_enzyme_candidate_info,
    get_main_enzyme_candidates_info,
    run_main_enzyme_candidate_info,
    run_main_enzyme_candidates_info,
)
from src.info_show.main_enzyme_sets_info import (
    get_main_enzyme_set_info,
    get_main_enzyme_sets_info,
    run_main_enzyme_set_info,
    run_main_enzyme_sets_info,
)
from src.info_show.protein_info import (
    get_protein_info,
    get_proteins_info,
    run_protein_info,
    run_proteins_info,
)
from src.info_show.retropath_info import (
    get_retropath_candidate_info,
    get_retropath_info,
    run_retropath_candidate_info,
    run_retropath_info,
)
from src.info_show.solution_info import get_solution_info, run_solution_info


def run_info(config: Any) -> dict[str, Any]:
    """Dispatch the unified ``info`` command."""

    main_enzyme_candidate = getattr(config, "main_enzyme_candidate", None)
    main_enzyme_set = getattr(config, "main_enzyme_set", None)
    retropath_candidate = getattr(config, "retropath_candidate", None)
    if (
        main_enzyme_candidate is not None
        and getattr(config, "step", None) is None
    ):
        raise ValueError("--main-enzyme-candidate 必须与 --step 一起使用")
    if getattr(config, "step", None) is not None:
        supports_step = (
            getattr(config, "solution", None) is not None
            or getattr(config, "main_enzyme_candidates", False)
            or main_enzyme_candidate is not None
            or retropath_candidate is not None
        )
        if not supports_step:
            raise ValueError(
                "--step 必须与 --solution、--main-enzyme-candidates、"
                "--main-enzyme-candidate 或 --retropath-candidate 一起使用"
            )
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
    if getattr(config, "retropath", False):
        return run_retropath_info(config)
    if retropath_candidate is not None:
        return run_retropath_candidate_info(config)
    if getattr(config, "protein", None) is not None:
        return run_protein_info(config)
    if getattr(config, "proteins", False):
        return run_proteins_info(config)
    if main_enzyme_set is not None:
        return run_main_enzyme_set_info(config)
    if getattr(config, "main_enzyme_sets", False):
        return run_main_enzyme_sets_info(config)
    if main_enzyme_candidate is not None:
        return run_main_enzyme_candidate_info(config)
    if getattr(config, "main_enzyme_candidates", False):
        return run_main_enzyme_candidates_info(config)
    if getattr(config, "solution", None) is not None:
        return run_solution_info(config)
    raise ValueError(
        "未指定信息查看类型，请使用 --chassis、--gap、--solution、"
        "--retropath、--retropath-candidate、"
        "--main-enzyme-candidates、--main-enzyme-candidate、"
        "--main-enzyme-sets、--main-enzyme-set、--proteins 或 --protein"
    )


__all__ = [
    "format_chassis_info_zh",
    "get_chassis_expansion_info",
    "get_chassis_info",
    "get_gap_info",
    "get_main_enzyme_candidate_info",
    "get_main_enzyme_candidates_info",
    "get_main_enzyme_set_info",
    "get_main_enzyme_sets_info",
    "get_protein_info",
    "get_proteins_info",
    "get_retropath_candidate_info",
    "get_retropath_info",
    "get_solution_info",
    "run_chassis_expansion_info",
    "run_chassis_info",
    "run_gap_info",
    "run_main_enzyme_candidate_info",
    "run_main_enzyme_candidates_info",
    "run_main_enzyme_set_info",
    "run_main_enzyme_sets_info",
    "run_protein_info",
    "run_proteins_info",
    "run_retropath_candidate_info",
    "run_retropath_info",
    "run_solution_info",
    "run_info",
]
