from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.kegg_gap_analyze import (
    ELECTRON_INFERENCE_VERSION,
    gap_depth_output_dir,
)
from src.pathway_analyze.retropath_mnxref import (
    MnxrefIndex,
    MnxrefIndexError,
    default_install_dir,
    default_rules_path,
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

VALIDATION_VALUE_NAMES = {
    "not_run": "未验证",
    "raw": "原始预测",
    "core_only": "仅核心反应，计量未补全",
    "complete": "已完成",
    "not_applicable": "不适用",
    "passed": "通过",
    "failed": "未通过",
    "strict": "严格模式",
    "relaxed": "宽松模式",
}

VALIDATION_RISK_NAMES = {
    "GEM validation has not been run": "尚未运行 GEM 验证",
    "this predicted route requires manual review before experimental use.": (
        "实验使用前必须人工复核"
    ),
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
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def _validation_name(value: Any) -> str:
    text = str(value or "").strip()
    return VALIDATION_VALUE_NAMES.get(text, text or "未知")


def _friendly_compound(
    value: Any,
    aliases: dict[str, str] | None = None,
) -> str:
    text = str(value or "").strip()
    if aliases and text in aliases:
        return aliases[text]
    if text.startswith("RP2CPD:"):
        identifier = text.split(":", 1)[1].split("-", 1)[0]
        return f"预测化合物 {identifier[:8]}"
    return text or "未知化合物"


def _friendly_compounds(
    values: Any,
    aliases: dict[str, str] | None = None,
) -> list[str]:
    if isinstance(values, list):
        items = values
    else:
        items = _split_values(values)
    return [_friendly_compound(item, aliases) for item in items]


def _friendly_named_compound(
    compound_id: Any,
    compound_name: Any,
    aliases: dict[str, str] | None = None,
) -> str:
    identifier = str(compound_id or "").strip()
    name = str(compound_name or "").strip()
    if aliases and identifier in aliases:
        return aliases[identifier]
    if identifier.startswith("RP2CPD:"):
        return _friendly_compound(identifier, aliases)
    if identifier and name and name != identifier:
        return f"{identifier} ({name})"
    return _friendly_compound(name or identifier, aliases)


def _retropath_anchor_aliases(
    gap_dir: Path,
    summary: dict[str, Any],
) -> dict[str, str]:
    candidate_rank = _to_int(
        summary.get("retropath_candidate_rank") or 0,
        "retropath_candidate_rank",
    )
    if candidate_rank < 1:
        return {}
    route_path = gap_dir / "retropath" / "candidate_routes.csv"
    if not route_path.is_file():
        return {}
    try:
        with route_path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, UnicodeError, csv.Error):
        return {}
    route = next(
        (
            row
            for row in rows
            if _to_int(row.get("candidate_rank") or 0, "candidate_rank")
            == candidate_rank
        ),
        None,
    )
    if route is None:
        return {}
    sink_ids = _split_values(route.get("sink_kegg_ids"))
    sink_keys = _split_values(route.get("sink_inchikeys"))
    anchor_ids = _split_values(summary.get("reachable_anchor_compounds"))
    anchor_labels = _split_values(summary.get("reachable_anchor_labels"))
    labels_by_id = dict(zip(anchor_ids, anchor_labels))
    return {
        f"RP2CPD:{inchikey}": labels_by_id.get(kegg_id, kegg_id)
        for kegg_id, inchikey in zip(sink_ids, sink_keys)
    }


def _retropath_compound_aliases(
    config: Any,
    gap_dir: Path,
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
) -> tuple[dict[str, str], bool]:
    aliases = _retropath_anchor_aliases(gap_dir, summary)
    predicted_ids = {
        compound_id
        for row in rows
        for field in ("produced_compound_id", "precursor_compound_ids")
        for compound_id in _split_values(row.get(field))
        if compound_id.startswith("RP2CPD:") and compound_id not in aliases
    }
    keys_by_id = {
        compound_id: compound_id.split(":", 1)[1].split("-", 1)[0]
        for compound_id in predicted_ids
        if len(compound_id.split(":", 1)[1].split("-", 1)[0]) == 14
    }
    if not keys_by_id:
        return aliases, False
    mnxref_dir = Path(
        getattr(
            config,
            "data_dir",
            default_install_dir().parents[2],
        )
    ) / "retropath" / "mnxref" / "3.0"
    rules_path = Path(
        getattr(config, "retropath_rules_path", default_rules_path())
    )
    try:
        with MnxrefIndex(mnxref_dir, rules_path) as index:
            matches = index.chemicals_by_connectivity_keys(
                keys_by_id.values()
            )
    except (OSError, MnxrefIndexError):
        return aliases, False

    enriched = False
    for compound_id, connectivity_key in sorted(keys_by_id.items()):
        candidates = matches.get(connectivity_key, tuple())
        names = {
            item.name.strip()
            for item in candidates
            if item.name.strip() and item.name.strip().lower() not in {"na", "none"}
        }
        if len(names) != 1:
            continue
        name = next(iter(names))
        exact = any(
            f"RP2CPD:{item.inchikey}" == compound_id for item in candidates
        )
        aliases[compound_id] = (
            name
            if exact
            else f"{name}*"
        )
        enriched = True
    return aliases, enriched


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
    compound_aliases: dict[str, str] | None = None,
) -> dict[str, Any]:
    status = row.get("status")
    result = {
        "步骤编号": display_step_index,
        "反应ID": row.get("reaction_id"),
        "反应类型": FIELD_VALUE_NAMES["status"].get(status, status),
        "步骤来源": (
            "RetroPath预测"
            if str(row.get("step_source") or "").strip() == "retropath"
            else "KEGG"
        ),
        "输入": _friendly_compounds(
            row.get("precursor_compound_labels")
            or row.get("precursor_compound_ids"),
            compound_aliases,
        ),
        "输出": _friendly_named_compound(
            row.get("produced_compound_id"),
            row.get("produced_compound_name"),
            compound_aliases,
        ),
        "反应名称": (
            "RetroPath预测反应"
            if str(row.get("reaction_name") or "").strip()
            == "RetroPath core prediction"
            else row.get("reaction_name")
        ),
        "酶EC编号": _split_values(
            row.get("enzyme_ecs") or row.get("source_ec_numbers")
        ),
    }
    if str(row.get("step_source") or "").strip() == "retropath":
        result.update({
            "RetroPath验证状态": _validation_name(row.get("validation_status")),
            "计量补全状态": _validation_name(row.get("stoichiometry_status")),
            "GEM验证状态": _validation_name(row.get("gem_status")),
        })
    return result


def _format_solution_info(result: dict[str, Any]) -> str:
    if result.get("需要新增合成路径") is False:
        return "\n".join(
            (
                "目标化合物已在底盘可达范围内",
                "",
                f"目标化合物：{result.get('目标化合物', '未知')}",
                str(result.get("提示") or "无需新增合成路线。"),
            )
        )

    solution_id = result.get("路径编号", "未知")
    target = result.get("目标化合物", "未知")
    depth = result.get("Gap深度", 0)
    source = result.get("路线来源", "未知")

    if "步骤详情" in result:
        detail = result["步骤详情"]
        step_index = result.get("步骤编号", "未知")
        inputs = _friendly_compounds(
            detail.get("前体化合物") or detail.get("前体化合物ID")
        )
        output = _friendly_compound(
            detail.get("生成化合物名称") or detail.get("生成化合物ID")
        )
        lines = [
            f"路线 {solution_id} · Step {step_index}",
            "",
            f"目标化合物：{target}",
            f"路线来源：{source}",
            f"输入：{' + '.join(inputs)}",
            f"输出：{output}",
        ]
        reaction_name = str(detail.get("反应名称") or "").strip()
        if reaction_name == "RetroPath core prediction":
            reaction_name = "RetroPath预测反应"
        if reaction_name:
            lines.append(f"反应：{reaction_name}")
        ecs = _split_values(
            detail.get("酶EC编号") or detail.get("来源模板EC")
        )
        if ecs:
            lines.append(f"酶 EC：{'、'.join(ecs)}")
        lines.extend(
            (
                "",
                "查看完整原始字段：",
                f"outputs/{target}/kegg_gap_{target}/depth{depth}/"
                "all_solution_steps.csv",
            )
        )
        return "\n".join(lines)

    total_steps = result.get("步骤数量", 0)
    heterologous_steps = result.get("异源步骤数", 0)
    starts = _friendly_compounds(result.get("可达起始化合物", []))
    lines = [
        f"路线 {solution_id}：{target}",
        "",
        f"路线来源：{source}",
        f"路线规模：共 {total_steps} 步，其中异源反应 {heterologous_steps} 步",
        f"起始代谢物：{'、'.join(starts)}",
        f"最大电子风险：{result.get('最大电子风险', '未知')}",
        f"电子系统：{result.get('电子系统结论', '未知')}",
    ]
    if source == "RetroPath预测":
        lines.append(
            "验证状态："
            f"{_validation_name(result.get('RetroPath验证状态'))}；"
            f"计量：{_validation_name(result.get('计量补全状态'))}；"
            f"GEM：{_validation_name(result.get('GEM验证状态'))}"
        )
        name_note = str(result.get("中间体名称说明") or "").strip()
        if name_note:
            lines.append(f"名称说明：{name_note}")
    risks = [
        VALIDATION_RISK_NAMES.get(str(item).strip(), str(item).strip())
        for item in result.get("验证风险", [])
        if str(item).strip()
    ]
    if risks:
        lines.append(f"验证提醒：{'；'.join(risks)}")

    lines.extend(("", "反应步骤："))
    for step in result.get("反应步骤", []):
        step_index = step.get("步骤编号", "未知")
        step_source = step.get("步骤来源") or step.get("反应类型") or "反应"
        reaction_name = str(step.get("反应名称") or "").strip()
        title = f"Step {step_index} · {reaction_name or step_source}"
        lines.append(title)
        inputs = _friendly_compounds(step.get("输入", []))
        lines.append(f"  输入：{' + '.join(inputs)}")
        lines.append(f"  输出：{_friendly_compound(step.get('输出'))}")
        ecs = _split_values(step.get("酶EC编号"))
        if ecs:
            lines.append(f"  酶 EC：{'、'.join(ecs)}")

    lines.extend(
        (
            "",
            "查看指定步骤：",
            "python main.py info -i <输入文件名> "
            f"--solution {solution_id} --step N -d {depth}",
        )
    )
    return "\n".join(lines)


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
    if solution_source == "retropath":
        compound_aliases, names_enriched = _retropath_compound_aliases(
            config,
            gap_dir,
            summary,
            forward_rows,
        )
    else:
        compound_aliases, names_enriched = {}, False

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
        if names_enriched:
            common["中间体名称说明"] = (
                "名称来自本地 MNXref 3.0；带 * 的名称仅为连接结构匹配，"
                "立体构型或质子化状态仍待确认"
            )
    if selected_step_index is not None:
        if selected_step_index > len(forward_rows):
            raise ValueError(
                f"路线 {selected_solution_id} 中没有步骤 {selected_step_index}；"
                f"可用步骤：{list(range(1, len(forward_rows) + 1))}"
            )
        selected_row = forward_rows[selected_step_index - 1]
        step_detail = _translate_row(selected_row)
        step_detail["步骤编号"] = selected_step_index
        step_detail["前体化合物"] = _friendly_compounds(
            selected_row.get("precursor_compound_labels")
            or selected_row.get("precursor_compound_ids"),
            compound_aliases,
        )
        step_detail["生成化合物名称"] = _friendly_named_compound(
            selected_row.get("produced_compound_id"),
            selected_row.get("produced_compound_name"),
            compound_aliases,
        )
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
            _step_overview(row, step_index, compound_aliases)
            for step_index, row in enumerate(forward_rows, start=1)
        ],
    }


def run_solution_info(config: Any) -> dict[str, Any]:
    """CLI entry point for ``info --solution``."""

    result = get_solution_info(config)
    print(_format_solution_info(result))
    return result


__all__ = ["get_solution_info", "run_solution_info"]
