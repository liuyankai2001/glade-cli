from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from src.main_protein_selection.common import (
    evidence_paths,
    file_summary,
    get_solution_steps,
    heterologous_requirements,
    read_manifest,
    unique_main_ecs,
)


def _solution_summary(manifest: dict[str, Any], solution_id: int, steps: list[dict[str, Any]]) -> dict[str, Any]:
    solution = manifest.get("solution") if isinstance(manifest.get("solution"), dict) else {}
    summary = solution.get("summary") if isinstance(solution.get("summary"), dict) else {}
    return {
        "source": solution.get("source"),
        "solution_id": solution_id,
        "target_compound_id": summary.get("target_compound_id"),
        "target_compound_name": summary.get("target_compound_name"),
        "step_count": len(steps),
        "gap_dir": solution.get("gap_dir"),
    }


def get_enzyme_system_context(
    *,
    manifest_path: str | Path,
    output_dir: str | Path,
    detail: Literal["compact", "full"] = "compact",
) -> dict[str, Any]:
    """
    读取当前路线与主酶候选选择上下文。

    限制：只读；不检索 UniProt、不分析风险、不写任何文件或 manifest。
    """

    manifest = read_manifest(manifest_path)
    solution_id, steps = get_solution_steps(manifest)
    requirements = heterologous_requirements(steps)
    protein_selection = manifest.get("protein_selection")
    protein_selection_summary = (
        protein_selection if isinstance(protein_selection, dict) else {}
    )
    paths = evidence_paths(output_dir)
    result = {
        "ok": True,
        "revision": manifest.get("revision"),
        "solution_summary": _solution_summary(manifest, solution_id, steps),
        "heterologous_step_count": len(requirements),
        "ec_numbers": unique_main_ecs(requirements),
        "protein_selection": {
            "available": isinstance(protein_selection, dict),
            "source": protein_selection_summary.get("source"),
            "selected_solution_id": protein_selection_summary.get("selected_solution_id"),
            "chassis_key": protein_selection_summary.get("chassis_key"),
            "recommended_design": protein_selection_summary.get(
                "recommended_design", {}
            ),
        },
        "evidence_files": {
            key: file_summary(path)
            for key, path in paths.items()
        },
        "next_actions": ["Run select_main_enzymes."],
    }
    if detail == "full":
        result["heterologous_requirements"] = requirements
        result["solution_steps"] = steps
    return result

