"""Read-only views of the proteins currently selected in the manifest."""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.write_manifest.store import read_design_manifest


_PREVIEW_HEAD = 40
_PREVIEW_TAIL = 20


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少有效的 {field_name}")
    return value


def _text(value: Any) -> str:
    return str(value or "").strip()


def _string_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return list(dict.fromkeys(
        _text(item) for item in value if _text(item)
    ))


def _step_list(value: Any) -> list[int]:
    if not isinstance(value, list):
        return []
    steps: set[int] = set()
    for item in value:
        try:
            step = int(item)
        except (TypeError, ValueError):
            continue
        if step > 0:
            steps.add(step)
    return sorted(steps)


def _record(records: dict[str, dict[str, Any]], accession: str) -> dict[str, Any]:
    normalized = _text(accession).upper()
    if not normalized:
        raise ValueError("蛋白记录缺少 ID")
    return records.setdefault(
        normalized,
        {
            "accession": normalized,
            "roles": set(),
            "protein_name": "",
            "organism_name": "",
            "assigned_step_indexes": set(),
            "reviewed": None,
            "cofactors": [],
            "cds": None,
            "cds_failure": None,
        },
    )


def _merge_main_proteins(
    manifest: Mapping[str, Any],
    records: dict[str, dict[str, Any]],
) -> Mapping[str, Any]:
    selection = _mapping(
        manifest.get("main_enzyme_selection"),
        "main_enzyme_selection；请先写入主酶组合",
    )
    raw_proteins = selection.get("proteins")
    if not isinstance(raw_proteins, list) or not raw_proteins:
        raise ValueError("main_enzyme_selection.proteins 不能为空")
    for raw in raw_proteins:
        item = _mapping(raw, "main_enzyme_selection.proteins[]")
        current = _record(records, _text(item.get("accession")))
        current["roles"].add("main_enzyme")
        current["protein_name"] = _text(item.get("protein_name"))
        current["organism_name"] = _text(item.get("organism_name"))
        current["assigned_step_indexes"].update(
            _step_list(item.get("assigned_step_indexes"))
        )
        current["reviewed"] = item.get("reviewed")
        current["cofactors"] = _string_list(item.get("cofactors"))
    return selection


def _merge_cds_status(
    manifest: Mapping[str, Any],
    records: dict[str, dict[str, Any]],
) -> str:
    selection = manifest.get("cds_selection")
    if not isinstance(selection, Mapping):
        return "未运行"
    raw_proteins = selection.get("proteins")
    if isinstance(raw_proteins, list):
        for raw in raw_proteins:
            if not isinstance(raw, Mapping):
                continue
            accession = _text(raw.get("accession")).upper()
            if accession in records:
                records[accession]["cds"] = dict(raw)
    raw_failures = selection.get("failures")
    if isinstance(raw_failures, list):
        for raw in raw_failures:
            if not isinstance(raw, Mapping):
                continue
            accession = _text(raw.get("accession")).upper()
            if accession in records:
                records[accession]["cds_failure"] = dict(raw)
    return _text(selection.get("status")) or "未知"


def _current_records(
    manifest: Mapping[str, Any],
) -> tuple[Mapping[str, Any], str, list[dict[str, Any]]]:
    records: dict[str, dict[str, Any]] = {}
    main_selection = _merge_main_proteins(manifest, records)
    cds_status = _merge_cds_status(manifest, records)
    ordered = sorted(
        records.values(),
        key=lambda item: (
            0 if "main_enzyme" in item["roles"] else 1,
            min(item["assigned_step_indexes"] or {999999}),
            item["accession"],
        ),
    )
    return main_selection, cds_status, ordered


def _cds_summary(record: Mapping[str, Any]) -> dict[str, Any]:
    cds = record["cds"]
    failure = record["cds_failure"]
    if isinstance(cds, Mapping):
        optimized = cds.get("optimized_cds")
        optimized = optimized if isinstance(optimized, Mapping) else {}
        if optimized.get("optimization_skipped") is True:
            status = "用户CDS直接使用"
        elif optimized.get("processing_mode") == "codon_transformer_only":
            status = "CodonTransformer已生成（未修正）"
        else:
            status = "已优化"
        return {
            "状态": status,
            "长度_nt": optimized.get("length_nt"),
            "路径": _text(optimized.get("path")),
            "跳过优化": optimized.get("optimization_skipped") is True,
        }
    if isinstance(failure, Mapping):
        return {
            "状态": "处理失败",
            "长度_nt": None,
            "路径": "",
            "跳过优化": False,
        }
    return {
        "状态": "未处理",
        "长度_nt": None,
        "路径": "",
        "跳过优化": False,
    }


def _compact_record(config: Any, record: Mapping[str, Any]) -> dict[str, Any]:
    cds = _cds_summary(record)
    return {
        "蛋白ID": record["accession"],
        "类型": "主酶",
        "角色": ["主酶"],
        "蛋白名称": record["protein_name"],
        "来源生物": record["organism_name"],
        "负责步骤": sorted(record["assigned_step_indexes"]),
        "来源": "UniProt",
        "CDS状态": cds["状态"],
        "CDS长度_nt": cds["长度_nt"],
    }


def _safe_project_file(config: Any, path_value: str) -> tuple[Path | None, str]:
    if not path_value:
        return None, ""
    project_root = Path(config.project_output_path).expanduser().resolve()
    path = (project_root / path_value).resolve()
    try:
        path.relative_to(project_root)
    except ValueError:
        return None, "路径超出项目输出目录"
    return path, ""


def _preview(sequence: str) -> str:
    if len(sequence) <= _PREVIEW_HEAD + _PREVIEW_TAIL:
        return sequence
    return sequence[:_PREVIEW_HEAD] + "..." + sequence[-_PREVIEW_TAIL:]


def _sequence_file_info(config: Any, record: Mapping[str, Any]) -> dict[str, Any]:
    path_value = ""
    source_kind = ""
    if not path_value and isinstance(record["cds"], Mapping):
        protein_sequence = record["cds"].get("protein_sequence")
        if isinstance(protein_sequence, Mapping):
            path_value = _text(protein_sequence.get("path"))
            source_kind = "protein-to-CDS 蛋白文件"
    path, path_error = _safe_project_file(config, path_value)
    if path is None:
        return {
            "路径": path_value,
            "来源": source_kind,
            "存在": False,
            "长度": None,
            "预览": "",
            "读取错误": path_error,
        }
    if not path.is_file():
        return {
            "路径": path_value,
            "来源": source_kind,
            "存在": False,
            "长度": None,
            "预览": "",
            "读取错误": "文件不存在",
        }
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        sequence = "".join(
            re.sub(r"\s+", "", line)
            for line in lines
            if line.strip() and not line.lstrip().startswith(">")
        ).upper()
    except (OSError, UnicodeError) as exc:
        return {
            "路径": path_value,
            "来源": source_kind,
            "存在": True,
            "长度": None,
            "预览": "",
            "读取错误": str(exc),
        }
    return {
        "路径": path_value,
        "来源": source_kind,
        "存在": True,
        "长度": len(sequence),
        "预览": _preview(sequence),
        "读取错误": "",
    }


def _load_current(config: Any) -> tuple[dict[str, Any], Mapping[str, Any], str, list[dict[str, Any]]]:
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    main_selection, cds_status, records = _current_records(manifest)
    return manifest, main_selection, cds_status, records


def get_proteins_info(config: Any) -> dict[str, Any]:
    """Return a compact list of all currently selected proteins."""

    manifest, main_selection, cds_status, records = _load_current(config)
    return {
        "运行成功": True,
        "目标化合物": _text(manifest.get("target_compound_id")),
        "主酶组合编号": main_selection.get("selected_set_id"),
        "CDS整体状态": cds_status,
        "蛋白总数": len(records),
        "主酶数": sum("main_enzyme" in item["roles"] for item in records),
        "蛋白列表": [_compact_record(config, item) for item in records],
    }


def get_protein_info(config: Any) -> dict[str, Any]:
    """Return one current protein with sequence and CDS details."""

    manifest, main_selection, cds_status, records = _load_current(config)
    requested = _text(getattr(config, "protein", None)).upper()
    by_accession = {item["accession"]: item for item in records}
    if requested not in by_accession:
        raise ValueError(
            f"未找到当前蛋白 {requested or '空'}；"
            f"可查看 ID：{sorted(by_accession)}"
        )
    record = by_accession[requested]
    result = _compact_record(config, record)
    result.update({
        "运行成功": True,
        "目标化合物": _text(manifest.get("target_compound_id")),
        "主酶组合编号": main_selection.get("selected_set_id"),
        "CDS整体状态": cds_status,
        "来源生物Reviewed": record["reviewed"],
        "辅因子": list(record["cofactors"]),
        "序列文件": _sequence_file_info(config, record),
        "CDS详情": {
            **_cds_summary(record),
            "失败信息": (
                dict(record["cds_failure"])
                if isinstance(record["cds_failure"], Mapping)
                else None
            ),
        },
    })
    return result


def run_proteins_info(config: Any) -> dict[str, Any]:
    result = get_proteins_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def run_protein_info(config: Any) -> dict[str, Any]:
    result = get_protein_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "get_protein_info",
    "get_proteins_info",
    "run_protein_info",
    "run_proteins_info",
]
