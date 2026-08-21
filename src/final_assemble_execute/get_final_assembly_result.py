"""Read-only inspection of a committed final-assembly output bundle."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.final_assemble_execute.common import inspect_file
from src.final_assemble_execute.config import FINAL_ASSEMBLY_SCHEMA_VERSION
from src.write_manifest.store import read_design_manifest


def get_final_assembly_result(config: Any) -> dict[str, Any]:
    project_root = Path(config.project_output_path).expanduser().resolve()
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    section = manifest.get("final_assembly")
    if not isinstance(section, Mapping):
        return {
            "ok": True,
            "available": False,
            "manifest_path": str(manifest_path),
            "manifest_revision": manifest.get("revision"),
        }
    if section.get("schema_version") != FINAL_ASSEMBLY_SCHEMA_VERSION:
        raise ValueError("manifest 中的 final_assembly 版本已过期")
    constructs = section.get("constructs")
    constructs = constructs if isinstance(constructs, list) else []
    inspected: list[dict[str, Any]] = []
    issue_count = 0
    for item in constructs:
        if not isinstance(item, Mapping):
            issue_count += 1
            continue
        files = item.get("files")
        files = files if isinstance(files, Mapping) else {}
        file_status = {
            key: inspect_file(project_root, value)
            for key, value in files.items()
            if isinstance(value, Mapping)
        }
        intact = bool(file_status) and all(
            value.get("exists") and value.get("file_sha256_matches")
            for value in file_status.values()
        )
        if not intact:
            issue_count += 1
        inspected.append(
            {
                "parts_design_id": item.get("parts_design_id"),
                "assembly_method": item.get("assembly_method"),
                "length_bp": item.get("length_bp"),
                "intact": intact,
                "files": file_status,
            }
        )
    top_level_files: dict[str, Any] = {}
    for key in ("run_summary_file", "report_file"):
        value = section.get(key)
        if isinstance(value, Mapping):
            top_level_files[key] = inspect_file(project_root, value)
            if not (
                top_level_files[key].get("exists")
                and top_level_files[key].get("file_sha256_matches")
            ):
                issue_count += 1
    return {
        "ok": issue_count == 0,
        "available": True,
        "status": section.get("status"),
        "target_compound_id": section.get("target_compound_id"),
        "counts": {
            "planned": section.get("planned_design_count"),
            "succeeded": section.get("succeeded_count"),
            "failed": section.get("failed_count"),
            "integrity_issues": issue_count,
        },
        "method_counts": section.get("method_counts"),
        "constructs": inspected,
        "files": top_level_files,
        "failures": section.get("failures"),
        "warnings": section.get("warnings"),
        "manifest_path": str(manifest_path),
        "manifest_revision": manifest.get("revision"),
    }


__all__ = ["get_final_assembly_result"]
