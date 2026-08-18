from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir
from src.pathway_analyze.target_id import validate_target_compound_id


STEP_FIELD_NAMES = {
    "solution_id": "路径编号",
    "step_index": "步骤编号",
    "status": "反应类型",
    "produced_compound_id": "生成化合物ID",
    "produced_compound_name": "生成化合物名称",
    "reaction_id": "反应ID",
    "reaction_name": "反应名称",
    "equation": "反应方程",
    "direction": "反应方向",
    "oxygen_required": "是否需要氧气",
    "thermo_direction": "热力学倾向",
    "screening_rule_hits": "筛选规则",
    "precursor_compound_ids": "前体化合物ID",
    "precursor_compound_labels": "前体化合物",
    "ko_ids": "KO编号",
    "enzyme_ecs": "酶EC编号",
    "source_reaction_ids": "原始反应ID",
    "resolution_action": "反应解析操作",
    "resolution_evidence": "反应解析证据",
    "module_ids": "KEGG模块ID",
    "rhea_ids": "Rhea反应ID",
}

FIELD_VALUE_NAMES = {
    "status": {
        "heterologous": "异源反应",
        "endogenous": "内源反应",
    },
    "direction": {
        "left_to_right": "从左到右",
        "right_to_left": "从右到左",
    },
    "thermo_direction": {
        "favored": "有利",
        "neutral": "中性",
        "disfavored": "不利",
    },
    "resolution_action": {
        "none": "无",
        "explicit_multistep_component": "多步反应的已解析组分",
        "merged_complete_reaction": "已合并为完整反应",
    },
}


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


def _translate_row(row: dict[str, Any]) -> dict[str, Any]:
    """翻译用户可见字段和枚举值，保留 KEGG/Rhea 原始数据。"""

    translated: dict[str, Any] = {}
    for field, value in row.items():
        output_field = STEP_FIELD_NAMES.get(field, field)
        translated[output_field] = FIELD_VALUE_NAMES.get(field, {}).get(value, value)
    return translated


def _to_int(value: Any, field_name: str) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer-compatible value: {value}") from exc


def get_solution_info(config: Any) -> dict[str, Any]:
    """读取 ``config.solution`` 对应深度的候选路径步骤。"""

    target_compound = validate_target_compound_id(config.target_name)
    raw_solution = getattr(
        config,
        "solution",
        getattr(config, "solution_id", None),
    )
    if raw_solution is None:
        raise ValueError("未指定 solution，请使用 info --solution N")
    selected_solution_id = int(raw_solution)
    expansion_depth = int(getattr(config, "depth", 0))
    if expansion_depth < 0:
        raise ValueError("depth 必须大于等于 0")
    raw_step_index = getattr(config, "step_index", None)
    selected_step_index = int(raw_step_index) if raw_step_index is not None else None
    gap_dir = gap_depth_output_dir(
        Path(config.gap_output_path).expanduser().resolve(),
        expansion_depth,
    )
    steps_path = gap_dir / "all_solution_steps.csv"

    all_rows = _read_csv_rows(steps_path)
    if not all_rows:
        raise ValueError(f"Gap analysis produced no solution steps: {steps_path}")
    required_columns = {"solution_id", "step_index"}
    missing_columns = required_columns.difference(all_rows[0])
    if missing_columns:
        raise ValueError(
            f"Missing columns in {steps_path}: {sorted(missing_columns)}"
        )

    rows = []
    available_solution_ids: set[int] = set()
    for row in all_rows:
        row_solution_id = _to_int(row.get("solution_id"), "solution_id")
        available_solution_ids.add(row_solution_id)
        if row_solution_id != selected_solution_id:
            continue
        if (
            selected_step_index is not None
            and _to_int(row.get("step_index"), "step_index")
            != selected_step_index
        ):
            continue
        rows.append(_normalize_row(row))

    rows.sort(key=lambda item: int(item.get("step_index") or 0))
    if not rows:
        if selected_solution_id not in available_solution_ids:
            raise ValueError(
                f"solution {selected_solution_id} not found; available solutions: "
                f"{sorted(available_solution_ids)}"
            )
        raise ValueError(
            f"step {selected_step_index} not found in solution {selected_solution_id}"
        )

    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "Gap深度": expansion_depth,
        "路径编号": selected_solution_id,
        "指定步骤编号": selected_step_index,
        "步骤数量": len(rows),
        "反应步骤": [_translate_row(row) for row in rows],
    }


def run_solution_info(config: Any) -> dict[str, Any]:
    """CLI entry point for ``info --solution``."""

    result = get_solution_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["get_solution_info", "run_solution_info"]
