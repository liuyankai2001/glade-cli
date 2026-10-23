"""Read-only validation of selected main enzymes and their generated CDS."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from src.expression_box.expression_host_context import describe_manifest_host
from src.write_manifest.store import read_design_manifest


def json_error(error: str, **details: Any) -> str:
    return json.dumps(
        {"ok": False, "error": error, **details}, ensure_ascii=False, default=str
    )


def build_cds_manifest_context(
    session_dir: str | None = None,
    manifest_path: str | None = None,
) -> dict[str, Any]:
    """Compare current CDS records with the explicitly selected main enzymes.

    ``session_dir`` identifies the project output directory when no manifest
    path is supplied. This adapter never resolves an external agent session.
    """
    from src.protein_to_cds.get_protein_selection_context import get_proteins_for_cds
    from src.protein_to_cds.write_protein_to_manifest import (
        CDS_SELECTION_SCHEMA_VERSION,
    )

    if manifest_path:
        path = Path(manifest_path).expanduser()
        if session_dir and not path.is_absolute():
            path = Path(session_dir).expanduser() / path
    elif session_dir:
        path = Path(session_dir).expanduser() / "design_manifest.json"
    else:
        raise ValueError("请提供 manifest_path，或用 session_dir 指定项目输出目录")
    path = path.resolve()
    context = get_proteins_for_cds(path)
    manifest = read_design_manifest(path)
    if manifest.get("revision", 0) != context.manifest_revision:
        raise ValueError("读取期间项目状态已变化，请重试")
    selection = manifest.get("cds_selection")
    selection = selection if isinstance(selection, dict) else {}
    if selection and selection.get("schema_version") != CDS_SELECTION_SCHEMA_VERSION:
        raise ValueError("CDS 结果已过期，请重新运行 protein-to-cds")
    rows = selection.get("proteins", [])
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise TypeError("cds_selection.proteins 必须是对象列表")
    expected = [protein.accession for protein in context.proteins]
    actual = [str(row.get("accession") or "").strip().upper() for row in rows]
    missing = sorted(set(expected) - set(actual))
    extra = sorted(set(actual) - set(expected))
    valid = bool(
        selection.get("status") == "complete"
        and selection.get("source_fingerprint") == context.source_fingerprint
        and len(actual) == len(set(actual)) == len(expected)
        and not missing
        and not extra
        and all(row.get("roles") == ["main_enzyme"] for row in rows)
    )
    cds_list = []
    for index, row in enumerate(rows, 1):
        optimized = row.get("optimized_cds", {})
        if not isinstance(optimized, dict):
            raise TypeError("optimized_cds 必须是对象")
        cds_path = str(optimized.get("path") or "")
        file_path = (path.parent / cds_path).resolve() if cds_path else None
        exists = bool(
            file_path and file_path.is_relative_to(path.parent) and file_path.is_file()
        )
        cds_list.append(
            {
                "index": index,
                "cds_id": row.get("accession", ""),
                "assigned_step_indexes": row.get("assigned_step_indexes", []),
                "protein": {
                    key: row.get(key, "")
                    for key in ("accession", "protein_name", "organism_name")
                },
                "optimized_cds": {
                    **optimized,
                    "sequence_file": {
                        "path": cds_path,
                        "exists": exists,
                        "sha256": optimized.get("file_sha256", ""),
                    },
                },
            }
        )
    return {
        "ok": True,
        "valid": valid,
        "message": "主酶 CDS 已全部生成"
        if valid
        else "CDS 与当前主酶选择不一致或尚未完成，请重新运行 protein-to-cds",
        "manifest_path": str(path),
        "revision": context.manifest_revision,
        "has_main_enzyme_selection": True,
        "has_cds_selection": bool(selection),
        "cds_selection_status": selection.get("status"),
        "host_context": describe_manifest_host(manifest),
        "recommended_accession_count": len(expected),
        "selected_cds_count": len(actual),
        "cds_count": len(actual),
        "recommended_accessions": expected,
        "selected_cds_accessions": actual,
        "missing_accessions": missing,
        "extra_selected_cds_accessions": extra,
        "cds_list": cds_list,
        "assembly_plan_hints": [
            {
                "cds_per_cassette": size,
                "cassette_count": (len(actual) + size - 1) // size,
                "description": f"每个表达盒最多组装 {size} 个 CDS",
            }
            for size in range(1, len(actual) + 1)
        ],
    }


def compact_validation_result(context: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in context.items()
        if key not in {"cds_list", "assembly_plan_hints"}
    }
