from __future__ import annotations

import csv
import json
import math
import time
from collections import Counter
from pathlib import Path
from typing import Any

import cobra
from cobra.util.solver import fix_objective_as_constraint

from src.pathway_analyze.target_status import (
    TARGET_ALREADY_AVAILABLE_STATUS,
)


DEFAULT_GROWTH_FRACTION = 0.1
DEFAULT_FLUX_THRESHOLD = 1e-8
DEFAULT_TEST_COMPARTMENTS = ("c",)
DEFAULT_PROGRESS_INTERVAL = 250

TARGET_SUPPLY_CONFIRMED = "direct_supply_confirmed"
TARGET_SUPPLY_NOT_DETECTED = "no_direct_supply_detected"
TARGET_SUPPLY_INDETERMINATE = "indeterminate"


def _test_compartments(config: Any) -> set[str] | None:
    value = getattr(config, "test_compartments", DEFAULT_TEST_COMPARTMENTS)
    if value is None:
        return None
    if isinstance(value, str):
        value = value.replace(",", " ").split()
    compartments = {str(item).strip() for item in value if str(item).strip()}
    return compartments or None


def _load_medium(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Medium JSON root must be an object: {path}")

    try:
        return {str(reaction_id): float(bound) for reaction_id, bound in payload.items()}
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Medium bounds must be numeric: {path}") from exc


def _kegg_ids(metabolite: cobra.Metabolite) -> list[str]:
    annotation = getattr(metabolite, "annotation", {}) or {}
    values: list[str] = []
    for key in ("kegg.compound", "kegg.drug", "kegg.glycan"):
        annotation_values = annotation.get(key, [])
        if not isinstance(annotation_values, (list, tuple, set)):
            annotation_values = [annotation_values]
        values.extend(str(value) for value in annotation_values if value is not None)
    return list(dict.fromkeys(values))


def _has_kegg_id(metabolite: cobra.Metabolite, kegg_id: str) -> bool:
    """Return whether a metabolite has one specific KEGG annotation."""

    normalized = str(kegg_id or "").strip().upper()
    return bool(normalized) and any(
        value.strip().upper() == normalized for value in _kegg_ids(metabolite)
    )


def _analyze_producibility(
    model: cobra.Model,
    growth_fraction: float,
    flux_threshold: float,
    compartments: set[str] | None,
    progress_interval: int = DEFAULT_PROGRESS_INTERVAL,
    target_compound: str = "",
) -> tuple[list[dict[str, str]], float, dict[str, Any]]:
    """在维持最低生长率时，逐个检测胞内代谢物的最大 demand flux。"""

    baseline_growth = model.slim_optimize()
    if baseline_growth is None or not math.isfinite(float(baseline_growth)):
        raise ValueError("当前培养基下模型求解失败，请检查模型和培养基设置")
    baseline_growth = float(baseline_growth)
    if baseline_growth <= flux_threshold:
        raise ValueError("当前培养基下模型几乎不生长，请检查培养基设置")
    required_growth = baseline_growth * growth_fraction

    external_counts = Counter()
    for reaction in model.exchanges:
        exchange_metabolites = list(reaction.metabolites)
        if len(exchange_metabolites) == 1:
            external_counts[exchange_metabolites[0].compartment] += 1
    external = (
        {external_counts.most_common(1)[0][0]} if external_counts else {"e"}
    )
    metabolites = [
        metabolite
        for metabolite in model.metabolites
        if metabolite.compartment not in external
        and (compartments is None or metabolite.compartment in compartments)
    ]
    target_mapped_metabolites = [
        metabolite
        for metabolite in model.metabolites
        if _has_kegg_id(metabolite, target_compound)
    ]
    target_tested_metabolite_ids = {
        metabolite.id
        for metabolite in metabolites
        if _has_kegg_id(metabolite, target_compound)
    }
    kegg_rows: list[dict[str, str]] = []
    target_supply_evidence: list[dict[str, Any]] = []
    producible_count = 0
    producible_with_kegg_count = 0
    producible_without_kegg_count = 0
    optimization_failed_count = 0
    below_threshold_count = 0
    target_optimization_failed_count = 0
    target_below_threshold_count = 0
    loop_started = time.perf_counter()
    total_metabolites = len(metabolites)

    print(f"[INFO] detected external compartment(s): {','.join(sorted(external))}")
    print(
        "[INFO] chassis screening scope: "
        f"{total_metabolites} metabolite(s), "
        "compartments="
        f"{','.join(sorted(compartments)) if compartments else 'all non-external'}"
    )
    print(f"[INFO] baseline growth: {baseline_growth:.6g}")
    print(
        f"[INFO] required minimum growth: {required_growth:.6g} "
        f"({growth_fraction:.1%} of baseline)"
    )

    for index, metabolite in enumerate(metabolites, start=1):
        with model:
            fix_objective_as_constraint(model, fraction=growth_fraction)
            demand = model.add_boundary(metabolite, type="demand")
            model.objective = demand
            model.objective_direction = "max"
            flux = model.slim_optimize()

        flux_is_valid = flux is not None and math.isfinite(float(flux))
        flux_value = float(flux) if flux_is_valid else None
        is_target_metabolite = metabolite.id in target_tested_metabolite_ids
        if flux_is_valid and flux_value > flux_threshold:
            producible_count += 1
            kegg_ids = _kegg_ids(metabolite)
            if kegg_ids:
                producible_with_kegg_count += 1
            else:
                producible_without_kegg_count += 1
            for kegg_id in kegg_ids:
                kegg_rows.append(
                    {
                        "source": "producible",
                        "met_id": metabolite.id,
                        "met_name": metabolite.name,
                        "compartment": metabolite.compartment,
                        "kegg_id": kegg_id,
                    }
                )
            if is_target_metabolite:
                target_supply_evidence.append(
                    {
                        "met_id": metabolite.id,
                        "met_name": metabolite.name,
                        "compartment": metabolite.compartment,
                        "max_demand_flux": flux_value,
                    }
                )
        elif not flux_is_valid:
            optimization_failed_count += 1
            if is_target_metabolite:
                target_optimization_failed_count += 1
        else:
            below_threshold_count += 1
            if is_target_metabolite:
                target_below_threshold_count += 1

        if index % progress_interval == 0 or index == total_metabolites:
            elapsed = time.perf_counter() - loop_started
            rate = index / elapsed if elapsed > 0 else 0.0
            remaining = total_metabolites - index
            eta = remaining / rate if rate > 0 else 0.0
            percentage = 100.0 * index / total_metabolites if total_metabolites else 100.0
            print(
                f"[INFO] chassis progress: {index}/{total_metabolites} "
                f"({percentage:.1f}%), producible={producible_count}, "
                f"with_kegg={producible_with_kegg_count}, "
                f"elapsed={elapsed:.1f}s, rate={rate:.1f} metabolites/s, "
                f"ETA={eta:.1f}s"
            )

    kegg_rows.sort(key=lambda row: (row["kegg_id"], row["met_id"]))
    target_supply_evidence.sort(key=lambda row: (row["met_id"], row["compartment"]))
    unique_kegg_count = len({row["kegg_id"] for row in kegg_rows})
    target_tested_count = len(target_tested_metabolite_ids)
    target_out_of_scope_count = len(target_mapped_metabolites) - target_tested_count
    if target_supply_evidence:
        target_supply_status = TARGET_SUPPLY_CONFIRMED
    elif target_optimization_failed_count:
        target_supply_status = TARGET_SUPPLY_INDETERMINATE
    else:
        target_supply_status = TARGET_SUPPLY_NOT_DETECTED
    stats: dict[str, Any] = {
        "required_growth": required_growth,
        "external_compartments": sorted(external),
        "tested_metabolites": total_metabolites,
        "producible_metabolites": producible_count,
        "producible_with_kegg": producible_with_kegg_count,
        "producible_without_kegg": producible_without_kegg_count,
        "optimization_failed": optimization_failed_count,
        "below_flux_threshold": below_threshold_count,
        "kegg_mapping_rows": len(kegg_rows),
        "unique_kegg_compounds": unique_kegg_count,
        "screening_elapsed_seconds": time.perf_counter() - loop_started,
        "target_supply_status": target_supply_status,
        "target_supply_evidence": target_supply_evidence,
        "target_mapped_metabolites": len(target_mapped_metabolites),
        "target_tested_metabolites": target_tested_count,
        "target_out_of_scope_metabolites": target_out_of_scope_count,
        "target_below_threshold": target_below_threshold_count,
        "target_optimization_failed": target_optimization_failed_count,
    }
    return kegg_rows, baseline_growth, stats


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def analyze_chassis_metabolites(config: Any) -> dict[str, Any]:
    analysis_started = time.perf_counter()
    model_path = Path(config.model_path).expanduser().resolve()
    medium_path = Path(config.medium_path).expanduser().resolve()
    output_dir = Path(config.chassis_output_path).expanduser().resolve()
    producible_csv = Path(config.chassis_producible_csv).expanduser().resolve()
    summary_csv = Path(config.chassis_metabolites_summary_csv).expanduser().resolve()
    growth_fraction = float(
        getattr(config, "growth_fraction", DEFAULT_GROWTH_FRACTION)
    )
    flux_threshold = float(getattr(config, "flux_threshold", DEFAULT_FLUX_THRESHOLD))
    compartments = _test_compartments(config)
    progress_interval = int(
        getattr(config, "chassis_progress_interval", DEFAULT_PROGRESS_INTERVAL)
    )

    if not model_path.is_file():
        raise FileNotFoundError(f"GEM model file not found: {model_path}")
    if model_path.suffix.lower() != ".json":
        raise ValueError(f"Only JSON GEM models are supported: {model_path}")
    if not medium_path.is_file():
        raise FileNotFoundError(f"Medium file not found: {medium_path}")
    if not 0 < growth_fraction <= 1:
        raise ValueError("growth_fraction must be in the interval (0, 1]")
    if flux_threshold <= 0:
        raise ValueError("flux_threshold must be greater than 0")
    if progress_interval < 1:
        raise ValueError("chassis_progress_interval must be greater than or equal to 1")

    print("[INFO] starting chassis metabolite analysis")
    print(f"[INFO] GEM model: {model_path}")
    print(f"[INFO] medium file: {medium_path}")
    print(
        f"[INFO] parameters: growth_fraction={growth_fraction}, "
        f"flux_threshold={flux_threshold:g}, "
        "test_compartments="
        f"{','.join(sorted(compartments)) if compartments else 'all non-external'}, "
        f"progress_interval={progress_interval}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    model = cobra.io.load_json_model(str(model_path))
    medium = _load_medium(medium_path)
    model.medium = medium
    print(
        "[INFO] model loaded: "
        f"metabolites={len(model.metabolites)}, "
        f"reactions={len(model.reactions)}, genes={len(model.genes)}, "
        f"exchanges={len(model.exchanges)}"
    )
    print(f"[INFO] active medium reactions: {len(medium)}")
    target_compound = str(getattr(config, "target_name", "") or "").strip().upper()
    kegg_rows, baseline_growth, stats = _analyze_producibility(
        model,
        growth_fraction,
        flux_threshold,
        compartments,
        progress_interval,
        target_compound,
    )
    target_already_available = (
        stats["target_supply_status"] == TARGET_SUPPLY_CONFIRMED
    )
    status = (
        TARGET_ALREADY_AVAILABLE_STATUS if target_already_available else "complete"
    )
    if stats["target_supply_status"] == TARGET_SUPPLY_CONFIRMED:
        message = f"底盘分析完成；检测到目标化合物 {target_compound} 的直接供给通量。"
    elif stats["target_supply_status"] == TARGET_SUPPLY_NOT_DETECTED:
        message = f"底盘分析完成；未检测到目标化合物 {target_compound} 的直接供给通量。"
    else:
        message = f"底盘分析完成；目标化合物 {target_compound} 的供给结果不确定。"
    _write_csv(
        producible_csv,
        ["source", "met_id", "met_name", "compartment", "kegg_id"],
        kegg_rows,
    )

    summary_rows = [
        {"item": "baseline_growth", "value": baseline_growth},
        {"item": "required_growth", "value": stats["required_growth"]},
        {"item": "growth_fraction", "value": growth_fraction},
        {"item": "flux_threshold", "value": flux_threshold},
        {
            "item": "test_compartments",
            "value": ";".join(sorted(compartments or ())) or "non_external",
        },
        {
            "item": "detected_external_compartments",
            "value": ";".join(stats["external_compartments"]),
        },
        {"item": "progress_interval", "value": progress_interval},
        {"item": "tested_metabolites", "value": stats["tested_metabolites"]},
        {
            "item": "producible_metabolites",
            "value": stats["producible_metabolites"],
        },
        {
            "item": "producible_with_kegg",
            "value": stats["producible_with_kegg"],
        },
        {
            "item": "producible_without_kegg",
            "value": stats["producible_without_kegg"],
        },
        {
            "item": "optimization_failed",
            "value": stats["optimization_failed"],
        },
        {
            "item": "below_flux_threshold",
            "value": stats["below_flux_threshold"],
        },
        {"item": "kegg_mapping_rows", "value": stats["kegg_mapping_rows"]},
        {"item": "producible_kegg_compounds", "value": len(kegg_rows)},
        {
            "item": "unique_producible_kegg_compounds",
            "value": stats["unique_kegg_compounds"],
        },
        {"item": "reachable_kegg_compounds", "value": len(kegg_rows)},
        {
            "item": "screening_elapsed_seconds",
            "value": round(stats["screening_elapsed_seconds"], 3),
        },
        {"item": "target_compound", "value": target_compound},
        {"item": "target_supply_status", "value": stats["target_supply_status"]},
        {
            "item": "target_mapped_metabolites",
            "value": stats["target_mapped_metabolites"],
        },
        {
            "item": "target_tested_metabolites",
            "value": stats["target_tested_metabolites"],
        },
        {
            "item": "target_out_of_scope_metabolites",
            "value": stats["target_out_of_scope_metabolites"],
        },
        {
            "item": "target_below_threshold",
            "value": stats["target_below_threshold"],
        },
        {
            "item": "target_optimization_failed",
            "value": stats["target_optimization_failed"],
        },
        {
            "item": "target_supply_evidence_json",
            "value": json.dumps(stats["target_supply_evidence"], ensure_ascii=False),
        },
        {"item": "status", "value": status},
        {
            "item": "target_already_available_in_chassis",
            "value": target_already_available,
        },
    ]
    _write_csv(summary_csv, ["item", "value"], summary_rows)

    total_elapsed = time.perf_counter() - analysis_started
    return {
        "ok": True,
        "status": status,
        "message": message,
        "model_path": str(model_path),
        "medium_path": str(medium_path),
        "target_compound": target_compound,
        "target_already_available_in_chassis": target_already_available,
        "pathway_search_required": not target_already_available,
        "target_supply_status": stats["target_supply_status"],
        "target_supply_evidence": stats["target_supply_evidence"],
        "target_mapped_metabolites": stats["target_mapped_metabolites"],
        "target_tested_metabolites": stats["target_tested_metabolites"],
        "target_out_of_scope_metabolites": stats["target_out_of_scope_metabolites"],
        "target_below_threshold": stats["target_below_threshold"],
        "target_optimization_failed": stats["target_optimization_failed"],
        "baseline_growth": baseline_growth,
        "growth_fraction": growth_fraction,
        "flux_threshold": flux_threshold,
        "test_compartments": sorted(compartments or ()),
        "detected_external_compartments": stats["external_compartments"],
        "progress_interval": progress_interval,
        "required_growth": stats["required_growth"],
        "tested_metabolites": stats["tested_metabolites"],
        "producible_metabolites": stats["producible_metabolites"],
        "producible_with_kegg": stats["producible_with_kegg"],
        "producible_without_kegg": stats["producible_without_kegg"],
        "optimization_failed": stats["optimization_failed"],
        "below_flux_threshold": stats["below_flux_threshold"],
        "kegg_mapping_rows": stats["kegg_mapping_rows"],
        "producible_kegg_compounds": len(kegg_rows),
        "unique_producible_kegg_compounds": stats["unique_kegg_compounds"],
        "elapsed_seconds": total_elapsed,
        "reachable_kegg_file": str(producible_csv),
        "summary_file": str(summary_csv),
    }


def _format_number(value: Any) -> str:
    """Format numeric values compactly for the human-facing CLI summary."""

    if value is None:
        return "-"
    if isinstance(value, bool):
        return "是" if value else "否"
    if isinstance(value, int):
        return str(value)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(numeric):
        return str(value)
    return f"{numeric:.6g}"


def _target_supply_description(result: dict[str, Any]) -> str:
    """Explain the target conclusion without directing the user to another step."""

    status = result["target_supply_status"]
    if status == TARGET_SUPPLY_CONFIRMED:
        return "在本次约束下检测到超过通量阈值的目标供给通量。"
    if result["target_optimization_failed"]:
        return "目标关联代谢物的优化失败，无法据此判断是否可直接供给。"
    return "在本次 GEM、培养基和生长约束下，未检测到目标的直接供给通量。"


def format_chassis_result_zh(result: dict[str, Any]) -> str:
    """Render the completed chassis analysis as a concise Chinese CLI report."""

    status_labels = {
        TARGET_SUPPLY_CONFIRMED: "检测到直接供给",
        TARGET_SUPPLY_NOT_DETECTED: "未检测到直接供给",
        TARGET_SUPPLY_INDETERMINATE: "结果不确定",
    }
    compartments = result["test_compartments"]
    scope = ", ".join(compartments) if compartments else "全部非胞外区室"
    lines = [
        "",
        "底盘供给分析完成",
        "",
        "目标供给结论",
        f"  目标：{result['target_compound']}",
        f"  状态：{status_labels.get(result['target_supply_status'], result['target_supply_status'])}",
        f"  含义：{_target_supply_description(result)}",
    ]
    evidence = result["target_supply_evidence"]
    if evidence:
        lines.append("  供给证据：")
        for item in evidence:
            lines.extend(
                [
                    f"    - GEM 代谢物：{item['met_id']}",
                    f"      名称：{item['met_name'] or '-'}",
                    f"      区室：{item['compartment']}",
                    f"      最大 demand flux：{_format_number(item['max_demand_flux'])}",
                ]
            )

    lines.extend(
        [
            "",
            "分析条件",
            f"  底盘模型：{Path(result['model_path']).name}",
            f"  培养基：{Path(result['medium_path']).name}",
            f"  生长约束：维持基线生长的 {float(result['growth_fraction']):.1%}",
            f"  基线生长：{_format_number(result['baseline_growth'])}",
            f"  最低维持生长：{_format_number(result['required_growth'])}",
            f"  检测区室：{scope}",
            f"  通量阈值：{_format_number(result['flux_threshold'])}",
            "",
            "分析结果",
            f"  检测代谢物：{_format_number(result['tested_metabolites'])}",
            f"  可供给代谢物：{_format_number(result['producible_metabolites'])}",
            "  可供给 KEGG 化合物（去重）："
            f"{_format_number(result['unique_producible_kegg_compounds'])}",
            f"  优化失败：{_format_number(result['optimization_failed'])}",
            "  可供给但缺少 KEGG 注释："
            f"{_format_number(result['producible_without_kegg'])}",
            "",
            "结果文件",
            f"  供给化合物表：{result['reachable_kegg_file']}",
            f"  分析摘要：{result['summary_file']}",
        ]
    )
    return "\n".join(lines)


def run_chassis(config: Any) -> dict[str, Any]:
    """运行入口；JSON 读取和 RunConfig 构造由 ``main.py`` 负责。"""

    result = analyze_chassis_metabolites(config)
    print(format_chassis_result_zh(result))
    return result
