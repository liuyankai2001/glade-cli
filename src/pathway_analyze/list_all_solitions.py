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
        raise FileNotFoundError(f"Missing solutions.csv: {path}")
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


def _list_all_solutions_impl(
    gap_dir: str | None = None,
    target_compound_id: str | None = None,
    output_dir: str | None = None,
) -> dict[str, Any]:
    gap_dir_path = _resolve_gap_dir(gap_dir, target_compound_id, output_dir)
    solutions_path = gap_dir_path / "solutions.csv"
    rows = [_normalize_row(row) for row in _read_csv_rows(solutions_path)]

    return {
        "ok": True,
        "gap_dir": str(gap_dir_path.resolve()),
        "solutions_csv": str(solutions_path.resolve()),
        "solution_count": len(rows),
        "solutions": rows,
    }


@tool
def list_all_solutions(
    target_compound_id: str | None = None,
) -> dict[str, Any]:
    """
    列出某个目标化合物的候选通路方案摘要。
    
    调用时机：kegg_gap_analyze 完成后，用户要比较可选 solution。
    输入：target_compound_id。
    返回：ok、solution 列表、排名指标、输出文件路径和 next_actions。
    限制：只读；不写 manifest，不做步骤详情或蛋白选择。
    """

    tool_name = "list_all_solutions"
    monitor.report_start(tool_name, {"target_compound_id": target_compound_id})
    try:
        result = _list_all_solutions_impl(
            gap_dir=None,
            target_compound_id=target_compound_id,
            output_dir=None,
        )
        monitor.report_end(tool_name, {"solution_count": result.get("solution_count")})
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise
