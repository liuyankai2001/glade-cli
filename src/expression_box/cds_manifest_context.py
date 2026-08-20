from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from src.tools.common.session_paths import design_manifest_file, session_dir as resolve_session_dir
from src.tools.expression_cassette_assembly_tools.expression_host_context import (
    describe_manifest_host,
)


def json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def build_cds_manifest_context(
    session_dir: str | None = None,
    manifest_path: str | None = None,
) -> dict[str, Any]:
    session_root = _session_root(session_dir, manifest_path)
    manifest_file = (
        _safe_path(session_root, manifest_path)
        if manifest_path
        else design_manifest_file()
    )
    manifest = _read_json(manifest_file)

    has_protein_selection = isinstance(manifest.get("protein_selection"), dict)
    has_cds_selection = isinstance(manifest.get("cds_selection"), dict)
    recommended_accessions = _recommended_accessions(manifest)
    selected_cds = _selected_cds_items(manifest)
    cds_list = [
        _cds_summary(session_root, item, index)
        for index, item in enumerate(selected_cds, start=1)
    ]
    selected_cds_accessions = [
        str(item.get("protein", {}).get("accession") or "").strip()
        for item in cds_list
        if isinstance(item.get("protein"), dict)
    ]

    recommended_set = set(recommended_accessions)
    selected_cds_accession_set = {
        accession for accession in selected_cds_accessions if accession
    }
    missing_accessions = [
        accession
        for accession in recommended_accessions
        if accession not in selected_cds_accession_set
    ]
    extra_accessions = [
        accession
        for accession in selected_cds_accession_set
        if accession not in recommended_set
    ]

    recommended_count = len(recommended_accessions)
    selected_cds_count = len(cds_list)
    is_valid = (
        has_protein_selection
        and has_cds_selection
        and recommended_count > 0
        and selected_cds_count > 0
        and recommended_count == selected_cds_count
        and not missing_accessions
        and not extra_accessions
    )
    message = _validation_message(
        has_protein_selection=has_protein_selection,
        has_cds_selection=has_cds_selection,
        recommended_count=recommended_count,
        selected_cds_count=selected_cds_count,
        missing_accessions=missing_accessions,
        extra_accessions=extra_accessions,
        is_valid=is_valid,
    )

    return {
        "ok": True,
        "valid": is_valid,
        "message": message,
        "manifest_path": _relpath(session_root, manifest_file),
        "revision": manifest.get("revision"),
        "has_protein_selection": has_protein_selection,
        "has_cds_selection": has_cds_selection,
        "cds_selection_status": _cds_selection_status(manifest),
        "host_context": describe_manifest_host(manifest),
        "recommended_accession_count": recommended_count,
        "selected_cds_count": selected_cds_count,
        "cds_count": selected_cds_count,
        "recommended_accessions": recommended_accessions,
        "selected_cds_accessions": selected_cds_accessions,
        "missing_accessions": missing_accessions,
        "extra_selected_cds_accessions": extra_accessions,
        "cds_list": cds_list,
        "assembly_plan_hints": _assembly_plan_hints(selected_cds_count),
    }


def compact_validation_result(context: dict[str, Any]) -> dict[str, Any]:
    return {
        "ok": context.get("ok", True),
        "valid": context.get("valid", False),
        "message": context.get("message"),
        "manifest_path": context.get("manifest_path"),
        "revision": context.get("revision"),
        "has_protein_selection": context.get("has_protein_selection", False),
        "has_cds_selection": context.get("has_cds_selection", False),
        "host_context": context.get("host_context", {}),
        "recommended_accession_count": context.get("recommended_accession_count", 0),
        "selected_cds_count": context.get("selected_cds_count", 0),
        "recommended_accessions": context.get("recommended_accessions", []),
        "selected_cds_accessions": context.get("selected_cds_accessions", []),
        "missing_accessions": context.get("missing_accessions", []),
        "extra_selected_cds_accessions": context.get("extra_selected_cds_accessions", []),
    }


def _validation_message(
    *,
    has_protein_selection: bool,
    has_cds_selection: bool,
    recommended_count: int,
    selected_cds_count: int,
    missing_accessions: list[str],
    extra_accessions: list[str],
    is_valid: bool,
) -> str:
    if is_valid:
        return "cds信息已全部获取"
    if not has_protein_selection:
        return "缺少 protein_selection，请先完成候选蛋白推荐"
    if recommended_count <= 0:
        return "protein_selection 中没有推荐 accession，请先完成候选蛋白推荐"
    if not has_cds_selection:
        return "缺少 cds_selection，请先完成 protein_to_cds_agent 阶段"
    if selected_cds_count <= 0:
        return "cds_selection.selected_cds 为空，请先写入优化 CDS"
    if missing_accessions:
        return "cds_selection 缺少部分推荐 accession 的 CDS"
    if extra_accessions:
        return "cds_selection 中存在不属于推荐列表的 CDS"
    if recommended_count != selected_cds_count:
        return "推荐 accession 数量与已选 CDS 数量不一致"
    return "cds没有全部获取信息"


def _session_root(session_dir: str | None, manifest_path: str | None) -> Path:
    if session_dir:
        return Path(session_dir).resolve()
    if manifest_path:
        path = Path(manifest_path).resolve()
        if path.name == "design_manifest.json" and path.parent.name == "outputs":
            return path.parent.parent.resolve()
    return resolve_session_dir()


def _safe_path(session_root: Path, path_value: str | Path) -> Path:
    path = Path(path_value)
    if not path.is_absolute():
        path = session_root / path
    path = path.resolve()
    session_root = session_root.resolve()

    if path != session_root and session_root not in path.parents:
        raise ValueError(f"path escapes session_dir: {path}")
    return path


def _relpath(session_root: Path, path: Path) -> str:
    return path.resolve().relative_to(session_root.resolve()).as_posix()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return data


def _split_accessions(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        raw_values = [str(item).strip() for item in value]
    else:
        raw_values = [
            item.strip()
            for item in re.split(r"[;,|，、\s]+", str(value).strip())
        ]

    result = []
    seen = set()
    for accession in raw_values:
        if accession and accession not in seen:
            seen.add(accession)
            result.append(accession)
    return result


def _recommended_accessions(manifest: dict[str, Any]) -> list[str]:
    protein_selection = manifest.get("protein_selection")
    if not isinstance(protein_selection, dict):
        return []

    recommended_design = protein_selection.get("recommended_design")
    if not isinstance(recommended_design, dict):
        return []

    return _split_accessions(recommended_design.get("selected_accessions"))


def _selected_cds_items(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    cds_selection = manifest.get("cds_selection")
    if not isinstance(cds_selection, dict):
        return []

    selected_cds = cds_selection.get("selected_cds")
    if selected_cds is None:
        return []
    if not isinstance(selected_cds, list):
        raise ValueError('design_manifest.json field "cds_selection.selected_cds" must be a list')

    return [item for item in selected_cds if isinstance(item, dict)]


def _cds_selection_status(manifest: dict[str, Any]) -> str | None:
    cds_selection = manifest.get("cds_selection")
    if not isinstance(cds_selection, dict):
        return None
    status = cds_selection.get("status")
    return str(status) if status is not None else None


def _file_exists(session_root: Path, path_value: Any) -> bool | None:
    if not path_value:
        return None
    try:
        return _safe_path(session_root, str(path_value)).exists()
    except ValueError:
        return False


def _cds_summary(session_root: Path, item: dict[str, Any], index: int) -> dict[str, Any]:
    protein = item.get("protein") if isinstance(item.get("protein"), dict) else {}
    optimized_cds = item.get("optimized_cds") if isinstance(item.get("optimized_cds"), dict) else {}
    sequence_file = (
        optimized_cds.get("sequence_file")
        if isinstance(optimized_cds.get("sequence_file"), dict)
        else {}
    )
    sequence_path = sequence_file.get("path")

    return {
        "index": index,
        "cds_id": item.get("cds_id", ""),
        "step_index": item.get("step_index"),
        "reaction_id": item.get("reaction_id", ""),
        "ec_number": item.get("ec_number", ""),
        "protein": {
            "accession": protein.get("accession", ""),
            "protein_name": protein.get("protein_name", ""),
            "organism_name": protein.get("organism_name", ""),
        },
        "optimized_cds": {
            "sequence_file": {
                "path": sequence_path,
                "exists": _file_exists(session_root, sequence_path),
                "sha256": sequence_file.get("sha256", ""),
            },
            "length_nt": optimized_cds.get("length_nt"),
            "gc_percent": optimized_cds.get("gc_percent"),
            "start_codon": optimized_cds.get("start_codon", ""),
            "stop_codon": optimized_cds.get("stop_codon", ""),
            "length_multiple_of_3": optimized_cds.get("length_multiple_of_3"),
        },
    }


def _assembly_plan_hints(cds_count: int) -> list[dict[str, Any]]:
    if cds_count <= 0:
        return []

    return [
        {
            "cds_per_cassette": size,
            "cassette_count": (cds_count + size - 1) // size,
            "description": f"每个表达盒最多组装 {size} 个 CDS",
        }
        for size in range(1, cds_count + 1)
    ]
