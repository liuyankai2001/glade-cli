from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.kegg_gap_analyze import (
    ELECTRON_INFERENCE_VERSION,
    gap_depth_output_dir,
)
from src.pathway_analyze.target_status import (
    TARGET_ALREADY_AVAILABLE_STATUS,
    read_target_already_available_status,
    target_already_available_message,
)
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
    "reaction_comment": "反应备注",
    "step_source": "步骤来源",
    "expansion_depth": "扩展深度",
    "expansion_anchor_compounds": "扩展锚点化合物",
    "retropath_step_id": "RetroPath步骤ID",
    "retropath_hypothesis_id": "计量假设ID",
    "retropath_rule_id": "RetroRules规则ID",
    "source_mnxr_id": "来源MNXR反应",
    "source_ec_numbers": "来源模板EC",
    "formal_mapping_exact": "是否精确映射已知反应",
    "full_reaction_smiles": "完整Reaction SMILES",
    "prediction_review_required": "预测风险待复核",
    "validation_status": "RetroPath验证状态",
    "stoichiometry_status": "计量补全状态",
    "gem_status": "GEM验证状态",
    "cofactor_mode": "辅因子验证模式",
    "cofactor_relaxed": "是否放宽通用载体",
    "opened_generic_compound_ids": "已开放通用载体",
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
        "retropath_prediction": "RetroPath预测反应",
    },
}


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"未找到路线结果文件，请先运行 gap：{path}")
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


def _split_values(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def _risk_name(value: Any) -> str:
    return {
        "none": "无",
        "low": "低",
        "medium": "中",
        "high": "高",
    }.get(str(value or "").strip().lower(), str(value or "").strip())


def _balance_name(value: Any) -> str:
    return {
        "not_applicable": "不涉及电子载体",
        "internally_balanced": "路线内已平衡",
        "unbalanced": "存在净消耗",
        "unresolved": "载体配对待确认",
    }.get(str(value or "").strip(), str(value or "").strip())


def _electron_status_name(value: Any) -> str:
    return {
        "not_required": "不需要额外电子系统",
        "internally_balanced": "电子载体已在路线内平衡",
        "internally_balanced_carrier_check_required": (
            "电子已内部平衡，但真实载体兼容性待确认"
        ),
        "carrier_check_required": "真实载体或电子伙伴待确认",
        "external_regeneration_required": "需要外部电子载体再生系统",
    }.get(str(value or "").strip(), str(value or "").strip())


def _select_solution_summary(
    rows: list[dict[str, str]],
    solution_id: int,
) -> dict[str, Any]:
    available_solution_ids: set[int] = set()
    for row in rows:
        row_solution_id = _to_int(row.get("solution_id"), "solution_id")
        available_solution_ids.add(row_solution_id)
        if row_solution_id == solution_id:
            return _normalize_row(row)
    raise ValueError(
        f"未找到路线 {solution_id}；可用路线：{sorted(available_solution_ids)}"
    )


def _step_overview(
    row: dict[str, Any],
    display_step_index: int,
) -> dict[str, Any]:
    status = row.get("status")
    result = {
        "步骤编号": display_step_index,
        "反应ID": row.get("reaction_id"),
        "反应类型": FIELD_VALUE_NAMES["status"].get(status, status),
    }
    if str(row.get("step_source") or "").strip() == "retropath":
        result.update({
            "RetroPath验证状态": row.get("validation_status"),
            "计量补全状态": row.get("stoichiometry_status"),
            "GEM验证状态": row.get("gem_status"),
        })
    return result


def _order_steps_forward(
    rows: list[dict[str, Any]],
    target_compound: str,
) -> list[dict[str, Any]]:
    """按底盘前体到目标产物的方向排列路线步骤。"""

    row_by_product = {
        str(row.get("produced_compound_id") or "").strip(): row
        for row in rows
        if str(row.get("produced_compound_id") or "").strip()
    }
    ordered_rows: list[dict[str, Any]] = []
    visited_products: set[str] = set()
    active_products: set[str] = set()

    def visit(product_id: str) -> None:
        if product_id in visited_products or product_id in active_products:
            return
        row = row_by_product.get(product_id)
        if row is None:
            return
        active_products.add(product_id)
        for precursor_id in _split_values(row.get("precursor_compound_ids")):
            visit(precursor_id)
        active_products.remove(product_id)
        visited_products.add(product_id)
        ordered_rows.append(row)

    visit(target_compound)
    for row in reversed(rows):
        product_id = str(row.get("produced_compound_id") or "").strip()
        if product_id:
            visit(product_id)
    return ordered_rows


def _build_path_chain(
    rows: list[dict[str, Any]],
    target_compound: str,
) -> str:
    """按“前体-反应->产物”格式生成从底盘前体到目标的路线链。"""

    row_by_product = {
        str(row.get("produced_compound_id") or "").strip(): row
        for row in rows
        if str(row.get("produced_compound_id") or "").strip()
    }
    active_compounds: set[str] = set()

    def render(compound_id: str) -> str:
        row = row_by_product.get(compound_id)
        if row is None or compound_id in active_compounds:
            return compound_id

        active_compounds.add(compound_id)
        precursor_ids = _split_values(row.get("precursor_compound_ids"))
        precursor_expressions = [render(item) for item in precursor_ids]
        active_compounds.remove(compound_id)

        if not precursor_expressions:
            precursor_text = "未知前体"
        else:
            precursor_text = "+".join(precursor_expressions)
            if len(precursor_expressions) > 1 and any(
                "->" in expression for expression in precursor_expressions
            ):
                precursor_text = f"({precursor_text})"

        reaction_id = str(row.get("reaction_id") or "未知反应").strip()
        return f"{precursor_text}-{reaction_id}->{compound_id}"

    return render(target_compound)


def get_solution_info(config: Any) -> dict[str, Any]:
    """读取路线概要，指定 ``config.step`` 时返回单步详情。"""

    target_compound = validate_target_compound_id(config.target_name)
    expansion_depth = int(getattr(config, "depth", 0))
    if expansion_depth < 0:
        raise ValueError("depth 必须大于等于 0")
    gap_dir = gap_depth_output_dir(
        Path(config.gap_output_path).expanduser().resolve(),
        expansion_depth,
    )
    target_status = read_target_already_available_status(
        gap_dir, target_compound
    )
    if target_status is not None:
        return {
            "运行成功": True,
            "运行状态": "目标化合物已在底盘细胞中",
            "原始状态代码": TARGET_ALREADY_AVAILABLE_STATUS,
            "提示": str(target_status.get("message") or "").strip()
            or target_already_available_message(target_compound),
            "目标化合物": target_compound,
            "Gap深度": expansion_depth,
            "需要新增合成路径": False,
        }
    raw_solution = getattr(config, "solution", None)
    if raw_solution is None:
        raise ValueError("未指定 solution，请使用 info --solution N")
    selected_solution_id = int(raw_solution)
    if selected_solution_id < 1:
        raise ValueError("solution 必须是正整数")
    raw_step_index = getattr(config, "step", None)
    selected_step_index = int(raw_step_index) if raw_step_index is not None else None
    if selected_step_index is not None and selected_step_index < 1:
        raise ValueError("step 必须是正整数")
    summaries_path = gap_dir / "solutions.csv"
    steps_path = gap_dir / "all_solution_steps.csv"

    summary_rows = _read_csv_rows(summaries_path)
    if not summary_rows:
        raise ValueError(f"gap 分析没有生成候选路线：{summaries_path}")
    if "solution_id" not in summary_rows[0]:
        raise ValueError(f"路线结果缺少 solution_id 字段：{summaries_path}")
    required_electron_columns = {
        "electron_balance_status",
        "electron_system_status",
        "requires_external_electron_regeneration",
        "requires_carrier_compatibility_check",
    }
    missing_electron_columns = required_electron_columns.difference(
        summary_rows[0]
    )
    summary = _select_solution_summary(summary_rows, selected_solution_id)
    solution_source = str(summary.get("solution_source") or "kegg").strip().lower()
    if solution_source == "retropath":
        from src.pathway_analyze.retropath_materialization import (
            MATERIALIZATION_MANIFEST_FILE_NAME,
            verify_retropath_solution_materialization,
        )
        if (gap_dir / "retropath" / MATERIALIZATION_MANIFEST_FILE_NAME).is_file():
            verify_retropath_solution_materialization(
                gap_dir=gap_dir,
                target_compound=target_compound,
                expansion_depth=expansion_depth,
                solution_id=selected_solution_id,
            )
        else:
            from src.pathway_analyze.retropath_promotion import (
                verify_retropath_solution_promotion,
            )

            verify_retropath_solution_promotion(
                gap_dir=gap_dir,
                target_compound=target_compound,
                expansion_depth=expansion_depth,
                solution_id=selected_solution_id,
            )
        if missing_electron_columns:
            raise ValueError(
                "RetroPath 路线缺少电子系统字段，请重新运行 gap --retropath"
            )
    else:
        run_config_path = gap_dir / "run_config.json"
        try:
            run_config = json.loads(run_config_path.read_text(encoding="utf-8-sig"))
        except (FileNotFoundError, json.JSONDecodeError) as exc:
            raise ValueError(f"无法读取 gap 运行配置：{run_config_path}") from exc
        if (
            not isinstance(run_config, dict)
            or run_config.get("electron_inference_version")
            != ELECTRON_INFERENCE_VERSION
            or missing_electron_columns
        ):
            raise ValueError(
                "gap 结果使用旧版电子推断，请重新运行："
                f"gap -i <输入文件> -d {expansion_depth}"
            )

    all_rows = _read_csv_rows(steps_path)
    if not all_rows:
        raise ValueError(f"gap 分析没有生成路线步骤：{steps_path}")
    required_columns = {"solution_id", "step_index"}
    missing_columns = required_columns.difference(all_rows[0])
    if missing_columns:
        raise ValueError(
            f"Missing columns in {steps_path}: {sorted(missing_columns)}"
        )

    rows: list[dict[str, Any]] = []
    for row in all_rows:
        row_solution_id = _to_int(row.get("solution_id"), "solution_id")
        if row_solution_id != selected_solution_id:
            continue
        rows.append(_normalize_row(row))

    rows.sort(key=lambda item: int(item.get("step_index") or 0))
    if not rows:
        raise ValueError(f"路线 {selected_solution_id} 没有反应步骤")
    forward_rows = _order_steps_forward(rows, target_compound)

    common = {
        "运行成功": True,
        "目标化合物": target_compound,
        "Gap深度": expansion_depth,
        "路径编号": selected_solution_id,
        "路线来源": "RetroPath预测" if solution_source == "retropath" else "KEGG",
    }
    if solution_source == "retropath":
        common.update({
            "RetroPath验证状态": summary.get("validation_status"),
            "计量补全状态": summary.get("stoichiometry_status"),
            "GEM验证状态": summary.get("gem_status"),
            "辅因子验证模式": summary.get("cofactor_mode"),
            "是否放宽通用载体": bool(summary.get("cofactor_relaxed")),
            "已开放通用载体": _split_values(
                summary.get("opened_generic_compound_ids")
            ),
            "验证风险": _split_values(summary.get("validation_issue")),
        })
    if selected_step_index is not None:
        if selected_step_index > len(forward_rows):
            raise ValueError(
                f"路线 {selected_solution_id} 中没有步骤 {selected_step_index}；"
                f"可用步骤：{list(range(1, len(forward_rows) + 1))}"
            )
        selected_row = forward_rows[selected_step_index - 1]
        step_detail = _translate_row(selected_row)
        step_detail["步骤编号"] = selected_step_index
        return {
            **common,
            "步骤编号": selected_step_index,
            "步骤详情": step_detail,
        }

    return {
        **common,
        "路径链": _build_path_chain(rows, target_compound),
        "步骤数量": _to_int(summary.get("total_steps") or len(rows), "total_steps"),
        "异源步骤数": _to_int(
            summary.get("heterologous_steps")
            or sum(row.get("status") == "heterologous" for row in rows),
            "heterologous_steps",
        ),
        "可达起始化合物": _split_values(
            summary.get("reachable_anchor_labels")
            or summary.get("reachable_anchor_compounds")
        ),
        "最大电子风险": _risk_name(summary.get("max_electron_risk_level")),
        "电子载体平衡": _balance_name(summary.get("electron_balance_status")),
        "电子系统结论": _electron_status_name(
            summary.get("electron_system_status")
        ),
        "需要额外电子再生系统": bool(
            summary.get("requires_external_electron_regeneration")
        ),
        "需要确认真实载体兼容性": bool(
            summary.get("requires_carrier_compatibility_check")
        ),
        "可以推荐": bool(summary.get("eligible_for_recommendation")),
        "预测风险待复核": bool(summary.get("prediction_review_required")),
        "RetroPath候选排名": summary.get("retropath_candidate_rank"),
        "RetroPath组合ID": summary.get("retropath_combination_id"),
        "反应步骤": [
            _step_overview(row, step_index)
            for step_index, row in enumerate(forward_rows, start=1)
        ],
    }


def run_solution_info(config: Any) -> dict[str, Any]:
    """CLI entry point for ``info --solution``."""

    result = get_solution_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["get_solution_info", "run_solution_info"]
