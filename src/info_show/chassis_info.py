"""Read-only presentation of completed chassis-analysis results."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.target_status import (
    TARGET_ALREADY_AVAILABLE_STATUS,
)


CHASSIS_INFO_SCHEMA_VERSION = "chassis_info.v2"

TARGET_SUPPLY_CONFIRMED = "direct_supply_confirmed"
TARGET_SUPPLY_NOT_DETECTED = "no_direct_supply_detected"
TARGET_SUPPLY_INDETERMINATE = "indeterminate"


def _convert_value(value: Any) -> Any:
    text = str(value or "").strip()
    if not text:
        return None
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        integer = int(text)
    except ValueError:
        try:
            return float(text)
        except ValueError:
            return text
    return integer


def _read_summary(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到底盘分析摘要，请先运行 chassis 命令："
            f"{path}"
        )

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"item", "value"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                f"底盘分析摘要必须包含 item 和 value 列：{path}"
            )
        summary: dict[str, Any] = {}
        for row in reader:
            item = str(row.get("item") or "").strip()
            if not item:
                continue
            if item in summary:
                raise ValueError(f"底盘分析摘要包含重复项目 {item!r}：{path}")
            summary[item] = _convert_value(row.get("value"))
    return summary


def _read_producible_compounds(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到底盘可生成化合物表，请先运行 chassis 命令："
            f"{path}"
        )

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "kegg_id" not in reader.fieldnames:
            raise ValueError(
                f"底盘可生成化合物表必须包含 kegg_id 列：{path}"
            )
        return [
            {
                str(key): str(value or "").strip()
                for key, value in row.items()
                if key is not None
            }
            for row in reader
        ]


def _int_or_none(value: Any) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _read_target_evidence(value: Any) -> list[dict[str, Any]]:
    """Read the target-specific FBA evidence persisted in the summary CSV."""

    if not isinstance(value, str) or not value.strip():
        return []
    try:
        payload = json.loads(value)
    except json.JSONDecodeError:
        return []
    if not isinstance(payload, list):
        return []
    return [item for item in payload if isinstance(item, dict)]


def _target_supply_message(
    target_compound_id: str,
    target_supply_status: str,
    optimization_failed: int | None,
) -> str:
    """Describe only the completed chassis-analysis conclusion."""

    if target_supply_status == TARGET_SUPPLY_CONFIRMED:
        return (
            f"在本次 GEM、培养基和生长约束下，检测到目标化合物 "
            f"{target_compound_id} 的直接供给通量。"
        )
    if target_supply_status == TARGET_SUPPLY_INDETERMINATE:
        if optimization_failed:
            return (
                f"目标化合物 {target_compound_id} 的关联代谢物优化失败，"
                "本次分析无法判断其是否可直接供给。"
            )
    return (
        f"在本次 GEM、培养基和生长约束下，未检测到目标化合物 "
        f"{target_compound_id} 的直接供给通量。"
    )


def get_chassis_info(config: Any) -> dict[str, Any]:
    """Return a concise, auditable view of existing chassis results."""

    summary_path = Path(config.chassis_metabolites_summary_csv).expanduser().resolve()
    compounds_path = Path(config.chassis_producible_csv).expanduser().resolve()
    summary = _read_summary(summary_path)
    compounds = _read_producible_compounds(compounds_path)

    target_compound_id = str(getattr(config, "target_name", "") or "").strip().upper()
    target_matches = [
        row
        for row in compounds
        if str(row.get("kegg_id") or "").strip().upper() == target_compound_id
    ]
    unique_kegg_ids = {
        str(row.get("kegg_id") or "").strip().upper()
        for row in compounds
        if str(row.get("kegg_id") or "").strip()
    }

    warnings: list[str] = []
    optimization_failed = _int_or_none(summary.get("optimization_failed")) or 0
    without_kegg = _int_or_none(summary.get("producible_without_kegg")) or 0
    if optimization_failed:
        warnings.append(
            f"有 {optimization_failed} 个候选代谢物优化失败"
        )
    if without_kegg:
        warnings.append(
            f"有 {without_kegg} 个可生成代谢物缺少 KEGG 注释"
        )

    reported_rows = _int_or_none(
        summary.get("kegg_mapping_rows")
        if summary.get("kegg_mapping_rows") is not None
        else summary.get("producible_kegg_compounds")
    )
    if reported_rows is not None and reported_rows != len(compounds):
        warnings.append(
            "摘要记录的 KEGG 映射行数与可生成化合物表不一致"
            f"({reported_rows} != {len(compounds)})"
        )

    target_already_available = bool(target_matches)
    target_supply_status = str(summary.get("target_supply_status") or "").strip()
    target_mapped_metabolites = _int_or_none(
        summary.get("target_mapped_metabolites")
    )
    target_tested_metabolites = _int_or_none(
        summary.get("target_tested_metabolites")
    )
    target_out_of_scope_metabolites = _int_or_none(
        summary.get("target_out_of_scope_metabolites")
    )
    target_optimization_failed = _int_or_none(
        summary.get("target_optimization_failed")
    )
    # Chassis output v2 initially used ``indeterminate`` when the target had
    # no GEM mapping or fell outside the configured screening compartments.
    # Both cases now have the same workflow meaning as any other negative
    # screen: no direct supply was detected and pathway search is required.
    if (
        target_supply_status == TARGET_SUPPLY_INDETERMINATE
        and not target_optimization_failed
    ):
        target_supply_status = TARGET_SUPPLY_NOT_DETECTED
    if target_supply_status not in {
        TARGET_SUPPLY_CONFIRMED,
        TARGET_SUPPLY_NOT_DETECTED,
        TARGET_SUPPLY_INDETERMINATE,
    }:
        # Results written before chassis output v2 cannot distinguish an
        # unmapped target from a screened-but-unsupplied one.
        target_supply_status = (
            TARGET_SUPPLY_CONFIRMED
            if target_already_available
            else TARGET_SUPPLY_NOT_DETECTED
        )
    target_supply_evidence = _read_target_evidence(
        summary.get("target_supply_evidence_json")
    )
    status = (
        TARGET_ALREADY_AVAILABLE_STATUS if target_already_available else "complete"
    )
    message = _target_supply_message(
        target_compound_id,
        target_supply_status,
        target_optimization_failed,
    )
    return {
        "schema_version": CHASSIS_INFO_SCHEMA_VERSION,
        "ok": True,
        "status": status,
        "message": message,
        "target_compound_id": target_compound_id,
        "target_producible_by_chassis": target_already_available,
        "pathway_search_required": not target_already_available,
        "target_matches": target_matches,
        "target_supply": {
            "status": target_supply_status,
            "mapped_metabolites": target_mapped_metabolites,
            "tested_metabolites": target_tested_metabolites,
            "out_of_scope_metabolites": target_out_of_scope_metabolites,
            "optimization_failed": target_optimization_failed,
            "evidence": target_supply_evidence,
        },
        "model_path": str(Path(config.model_path).expanduser().resolve()),
        "medium_path": str(Path(config.medium_path).expanduser().resolve()),
        "growth": {
            "baseline_growth": summary.get("baseline_growth"),
            "required_growth": summary.get("required_growth"),
            "growth_fraction": summary.get("growth_fraction"),
        },
        "screening": {
            "flux_threshold": summary.get("flux_threshold"),
            "test_compartments": summary.get("test_compartments"),
            "detected_external_compartments": summary.get(
                "detected_external_compartments"
            ),
            "tested_metabolites": summary.get("tested_metabolites"),
            "producible_metabolites": summary.get("producible_metabolites"),
            "producible_with_kegg": summary.get("producible_with_kegg"),
            "producible_without_kegg": summary.get("producible_without_kegg"),
            "below_flux_threshold": summary.get("below_flux_threshold"),
            "optimization_failed": summary.get("optimization_failed"),
            "screening_elapsed_seconds": summary.get(
                "screening_elapsed_seconds"
            ),
        },
        "kegg_compounds": {
            "mapping_rows": len(compounds),
            "unique_compounds": len(unique_kegg_ids),
            "reported_mapping_rows": reported_rows,
            "reported_unique_compounds": summary.get(
                "unique_producible_kegg_compounds"
            ),
        },
        "warnings": warnings,
        "files": {
            "summary_csv": str(summary_path),
            "producible_compounds_csv": str(compounds_path),
        },
    }


def format_chassis_info_zh(result: dict[str, Any]) -> dict[str, Any]:
    """Return a compact Chinese CLI view of the chassis result."""

    growth = result["growth"]
    screening = result["screening"]
    compounds = result["kegg_compounds"]
    target_supply = result["target_supply"]
    target_status_labels = {
        TARGET_SUPPLY_CONFIRMED: "检测到直接供给",
        TARGET_SUPPLY_NOT_DETECTED: "未检测到直接供给",
        TARGET_SUPPLY_INDETERMINATE: "结果不确定",
    }
    return {
        "运行成功": result["ok"],
        "运行状态": "底盘分析完成",
        "提示": result["message"],
        "目标化合物ID": result["target_compound_id"],
        "目标供给状态": target_status_labels.get(
            target_supply["status"], target_supply["status"]
        ),
        "目标映射GEM代谢物数": target_supply["mapped_metabolites"],
        "目标纳入检测的代谢物数": target_supply["tested_metabolites"],
        "目标供给证据": target_supply["evidence"],
        "底盘模型": Path(result["model_path"]).name,
        "培养基": Path(result["medium_path"]).name,
        "基线生长通量": growth["baseline_growth"],
        "检测代谢物数": screening["tested_metabolites"],
        "可生成代谢物数": screening["producible_metabolites"],
        "可生成KEGG化合物数": compounds["unique_compounds"],
        "警告": result["warnings"],
    }


def run_chassis_info(config: Any) -> dict[str, Any]:
    """CLI entry point for ``info --chassis``."""

    result = format_chassis_info_zh(get_chassis_info(config))
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "CHASSIS_INFO_SCHEMA_VERSION",
    "format_chassis_info_zh",
    "get_chassis_info",
    "run_chassis_info",
]
