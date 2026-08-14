from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from langchain.tools import tool

from src.pathway_analyze.target_id import validate_target_compound_id
from src.runtime.monitor import monitor
from src.tools.common.session_paths import kegg_gap_dir, outputs_dir


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing all_solution_steps.csv: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _optional_number(value: str) -> int | float | str | bool | None:
    text = str(value or "").strip()
    if text == "":
        return None
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        number = float(text)
    except ValueError:
        return text
    return int(number) if number.is_integer() else number


def _normalize_row(row: dict[str, str]) -> dict[str, Any]:
    return {key: _optional_number(value) for key, value in row.items()}


def _resolve_gap_dir(
    gap_dir: str | None,
    target_compound_id: str | None,
    output_dir: str | None,
) -> Path:
    if gap_dir:
        path = Path(gap_dir)
        if path.is_absolute():
            return path.resolve()
        if path.exists():
            return path.resolve()
        return (outputs_dir(output_dir) / path).resolve()

    if not target_compound_id:
        raise ValueError("Provide gap_dir or target_compound_id.")

    target_compound_id = validate_target_compound_id(target_compound_id)
    return kegg_gap_dir(target_compound_id, output_dir).resolve()


def _to_int(value: Any, field_name: str) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer-compatible value: {value}") from exc


@tool
def list_solution_steps(
    solution_id: int,
    step_index: int | None = None,
    target_compound_id: str | None = None,
) -> dict[str, Any]:
    """
    读取某个候选通路方案的步骤明细。
    
    调用时机：kegg_gap_analyze 后，用户想查看 solution 的反应、化合物和酶信息。
    输入：target_compound_id、solution_id。
    返回：ok、步骤列表、EC/KO、化合物摘要和证据文件路径。
    限制：只读；不写 manifest，不做蛋白选择。
    """
    tool_name = "list_solution_steps"
    monitor.report_start(tool_name, {"solution_id": solution_id, "step_index": step_index, "target_compound_id": target_compound_id})
    try:
        selected_solution_id = int(solution_id)
        selected_step_index = int(step_index) if step_index is not None else None
        gap_dir_path = _resolve_gap_dir(None, target_compound_id, None)
        steps_path = gap_dir_path / "all_solution_steps.csv"

        rows = []
        for row in _read_csv_rows(steps_path):
            if _to_int(row.get("solution_id"), "solution_id") != selected_solution_id:
                continue
            if selected_step_index is not None and _to_int(row.get("step_index"), "step_index") != selected_step_index:
                continue
            rows.append(_normalize_row(row))

        rows.sort(key=lambda item: int(item.get("step_index") or 0))
        if not rows:
            target = f"solution_id={selected_solution_id}"
            if selected_step_index is not None:
                target += f", step_index={selected_step_index}"
            raise ValueError(f"No matching solution steps found in {steps_path}: {target}")

        result = {
            "ok": True,
            "gap_dir": str(gap_dir_path.resolve()),
            "all_solution_steps_csv": str(steps_path.resolve()),
            "solution_id": selected_solution_id,
            "step_index": selected_step_index,
            "step_count": len(rows),
            "steps": rows,
        }
        monitor.report_end(tool_name, {"step_count": len(rows)})
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise
