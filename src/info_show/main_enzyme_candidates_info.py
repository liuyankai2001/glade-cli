"""Compact read-only view of main-enzyme candidates grouped by route step."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from src.main_protein_selection.common import get_solution_steps, read_manifest
from src.main_protein_selection.models import MainEnzymeSelectionResult
from src.main_protein_selection.provenance import solution_fingerprint


FIT_STATUS_NAMES = {
    "verified": "已验证",
    "verified_with_risk": "已验证，但需要复核",
}

DIRECTION_NAMES = {
    "supported": "支持",
    "unknown": "待确认",
    "contradicted": "矛盾",
}

CONFIDENCE_NAMES = {
    "high": "高",
    "medium": "中",
    "low": "低",
}


def _selection_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "main_protein_selection"
        / "main_enzyme_selection.json"
    )


def _read_selection(path: Path) -> MainEnzymeSelectionResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到主酶候选，请先运行：main-enzyme -i <输入文件>"
        )
    try:
        return MainEnzymeSelectionResult.model_validate_json(
            path.read_text(encoding="utf-8-sig")
        )
    except ValueError as exc:
        raise ValueError(f"主酶候选结果格式无效：{path}") from exc


def get_main_enzyme_candidates_info(config: Any) -> dict[str, Any]:
    """读取当前已选路线的主酶候选，并按正向步骤编号分组。"""

    raw_step = getattr(config, "step", None)
    try:
        selected_step = int(raw_step) if raw_step is not None else None
    except (TypeError, ValueError) as exc:
        raise ValueError(f"step 必须是正整数，当前值：{raw_step!r}") from exc
    if selected_step is not None and selected_step < 1:
        raise ValueError("step 必须是正整数")

    selection = _read_selection(_selection_path(config))
    manifest = read_manifest(config.manifest_output_path)
    solution_id, steps = get_solution_steps(manifest)
    current_fingerprint = solution_fingerprint(solution_id, steps)

    if selection.selected_solution_id != solution_id:
        raise ValueError(
            "主酶候选对应的路线与 manifest 当前路线不一致，"
            "请重新运行 main-enzyme"
        )
    if selection.solution_fingerprint != current_fingerprint:
        raise ValueError(
            "manifest 中的路线已经发生变化，请重新运行 main-enzyme"
        )

    steps_by_index = {
        int(step.get("step_index") or 0): step
        for step in steps
        if int(step.get("step_index") or 0) > 0
    }
    if selected_step is not None and selected_step not in steps_by_index:
        raise ValueError(
            f"路线中不存在 Step {selected_step}；"
            f"可用步骤：{sorted(steps_by_index)}"
        )

    if selected_step is None:
        step_indexes = sorted(selection.candidates_by_step)
    else:
        step_indexes = [selected_step]

    grouped_candidates: list[dict[str, Any]] = []
    for step_index in step_indexes:
        step = steps_by_index.get(step_index)
        if step is None:
            raise ValueError(
                f"主酶候选引用了 manifest 中不存在的步骤：{step_index}"
            )
        candidates = sorted(
            selection.candidates_by_step.get(step_index, []),
            key=lambda candidate: candidate.candidate_rank,
        )
        candidate_note = ""
        if not candidates:
            if str(step.get("status") or "").strip() == "endogenous":
                candidate_note = "该步骤为内源反应，无需选择主酶"
            else:
                candidate_note = "该步骤没有可用的主酶候选"
        grouped_candidates.append({
            "步骤编号": step_index,
            "反应ID": step.get("reaction_id"),
            "候选数量": len(candidates),
            "说明": candidate_note,
            "候选酶": [
                {
                    "排名": candidate.candidate_rank,
                    "UniProt": candidate.accession,
                    "蛋白名称": candidate.protein_name,
                    "来源物种": candidate.organism_name,
                    "EC编号": candidate.ec_number,
                    "Reviewed": candidate.reviewed,
                    "蛋白评分": candidate.protein_score,
                    "反应匹配": FIT_STATUS_NAMES.get(
                        candidate.reaction_fit_status,
                        candidate.reaction_fit_status,
                    ),
                    "反应匹配分数": candidate.reaction_fit_score,
                    "方向判断": DIRECTION_NAMES.get(
                        candidate.direction_verdict,
                        candidate.direction_verdict,
                    ),
                    "方向置信度": CONFIDENCE_NAMES.get(
                        candidate.direction_confidence,
                        candidate.direction_confidence,
                    ),
                }
                for candidate in candidates
            ],
        })

    def filter_step_indexes(values: list[int]) -> list[int]:
        if selected_step is None:
            return values
        return [value for value in values if value == selected_step]

    result: dict[str, Any] = {
        "运行成功": True,
        "目标化合物": str(config.target_name),
        "路径编号": solution_id,
        "候选生成状态": selection.status,
    }
    if selected_step is not None:
        result["查看步骤"] = selected_step
    result.update({
        "未覆盖步骤": filter_step_indexes(selection.uncovered_step_indexes),
        "方向待确认步骤": filter_step_indexes(
            selection.direction_risk_step_indexes
        ),
        "步骤候选": grouped_candidates,
    })
    return result


def run_main_enzyme_candidates_info(config: Any) -> dict[str, Any]:
    """CLI-compatible entry point for ``info --main-enzyme-candidates``."""

    result = get_main_enzyme_candidates_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "get_main_enzyme_candidates_info",
    "run_main_enzyme_candidates_info",
]
