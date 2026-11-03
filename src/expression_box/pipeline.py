"""Pipeline for system-recommended expression-box protein groupings."""

from __future__ import annotations

import json
import os
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from src.expression_box.config import (
    BALANCED_MAX_MAIN_UNITS_PER_CASSETTE,
    COMPACT_MAX_CDS_LENGTH_NT,
    COMPACT_MAX_PROTEINS,
    EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION,
    GROUPING_ALGORITHM_VERSION,
)
from src.expression_box.manifest_adapter import load_expression_grouping_context
from src.expression_box.models import (
    ExpressionGroupingContext,
    ExpressionGroupingDesign,
    ExpressionProtein,
    GroupingGenerationResult,
)
from src.expression_box.protein_grouping import generate_grouping_designs

EXPRESSION_BOX_DESIGNS_FILENAME = "expression_box_designs.json"


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _protein_payload(protein: ExpressionProtein) -> dict[str, Any]:
    return {
        "accession": protein.accession,
        "roles": list(protein.roles),
        "assigned_step_indexes": list(protein.assigned_step_indexes),
        "optimized_cds_length_nt": protein.optimized_cds_length_nt,
    }


def _design_payload(
    design: ExpressionGroupingDesign,
    design_id: int,
) -> dict[str, Any]:
    cassettes = []
    for cassette_index, cassette in enumerate(design.cassettes, start=1):
        cassettes.append(
            {
                "cassette_index": cassette_index,
                "protein_accessions": [
                    protein.accession for protein in cassette.proteins
                ],
                "protein_count": len(cassette.proteins),
                "main_enzyme_count": cassette.main_enzyme_count,
                "total_cds_length_nt": cassette.total_cds_length_nt,
                "reason": cassette.reason,
                "proteins": [
                    _protein_payload(protein) for protein in cassette.proteins
                ],
            }
        )
    return {
        "design_id": design_id,
        "rank": design_id,
        "strategy": design.strategy,
        "name": design.name,
        "recommended": design.recommended,
        "cassette_count": len(cassettes),
        "protein_count": sum(item["protein_count"] for item in cassettes),
        "total_cds_length_nt": sum(item["total_cds_length_nt"] for item in cassettes),
        "cassettes": cassettes,
        "warnings": list(design.warnings),
    }


def _full_payload(
    context: ExpressionGroupingContext,
    generation: GroupingGenerationResult,
) -> dict[str, Any]:
    designs = [
        _design_payload(design, design_id)
        for design_id, design in enumerate(generation.designs, start=1)
    ]
    return {
        "schema_version": EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION,
        "status": "complete",
        "generated_at": _utc_now(),
        "algorithm_version": GROUPING_ALGORITHM_VERSION,
        "source": {
            "manifest_path": str(context.manifest_path),
            "manifest_revision": context.manifest_revision,
        },
        "constraints": {
            "balanced_max_main_units_per_cassette": (
                BALANCED_MAX_MAIN_UNITS_PER_CASSETTE
            ),
            "compact_max_proteins": COMPACT_MAX_PROTEINS,
            "compact_max_cds_length_nt": COMPACT_MAX_CDS_LENGTH_NT,
        },
        "target_compound_id": context.target_compound_id,
        "protein_count": len(context.proteins),
        "design_count": len(designs),
        "designs": designs,
        "skipped_strategies": [
            {
                "strategy": item.strategy,
                "reason": item.reason,
            }
            for item in generation.skipped_strategies
        ],
        "warnings": list(generation.warnings),
    }


def _summary(
    payload: dict[str, Any],
    output_path: Path,
    *,
    reused_existing: bool,
) -> dict[str, Any]:
    return {
        "ok": True,
        "status": payload["status"],
        "target_compound_id": payload["target_compound_id"],
        "protein_count": payload["protein_count"],
        "design_count": payload["design_count"],
        "designs": [
            {
                "design_id": item["design_id"],
                "name": item["name"],
                "strategy": item["strategy"],
                "recommended": item["recommended"],
                "cassette_count": item["cassette_count"],
                "cassettes": [
                    cassette["protein_accessions"] for cassette in item["cassettes"]
                ],
            }
            for item in payload["designs"]
        ],
        "skipped_strategies": payload["skipped_strategies"],
        "warnings": payload["warnings"],
        "output_path": str(output_path),
        "reused_existing": reused_existing,
        "manifest_modified": True,
    }



def _invalidate_calculation(manifest_path: Path, *, grouping: bool, preserve_draft: bool = False) -> None:
    from src.write_manifest.store import read_design_manifest, update_design_manifest
    from src.expression_box.parts_preparation import normalize_parts_draft, stale_parts_designs
    manifest = read_design_manifest(manifest_path)
    target = manifest.get("target_compound_id")
    if not target:
        return
    downstream = ("parts_selection", "assembled_expression_cassettes",
                  "assembled_expression_constructs", "plasmid_selection",
                  "final_assembly_plan", "final_assembly", "final_design_report")
    if grouping:
        downstream = ("expression_box_selection", "expression_cassette_assembly", *downstream)
    sections = {}
    raw_draft = manifest.get("expression_parts_draft")
    if isinstance(raw_draft, dict) and (preserve_draft or raw_draft.get("source_type") == "user_uploaded" or raw_draft.get("schema_version") == "expression_parts_draft.v1"):
        try:
            draft = normalize_parts_draft(raw_draft)
            stale_parts_designs(draft["designs"])
            sections["expression_parts_draft"] = draft
        except (KeyError, TypeError, ValueError):
            pass
    elif raw_draft is not None:
        downstream = (*downstream, "expression_parts_draft")
    update_design_manifest(manifest_path, target_compound_id=target, sections=sections,
                           discard_sections=downstream, expected_revision=int(manifest["revision"]))

def _run_expression_box_design(config: Any) -> dict[str, Any]:
    """Generate and persist system-recommended protein grouping candidates."""

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    _invalidate_calculation(manifest_path, grouping=True)
    output_path = project_root / "expression_box" / EXPRESSION_BOX_DESIGNS_FILENAME
    _write_json_atomic(output_path, {"schema_version": EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION, "status": "in_progress"})
    context = load_expression_grouping_context(manifest_path)
    current_manifest = manifest_path.read_bytes() if manifest_path.is_file() else None
    generation = generate_grouping_designs(context.proteins)
    payload = _full_payload(context, generation)
    if current_manifest is not None and manifest_path.read_bytes() != current_manifest:
        raise ValueError("manifest 在表达盒计算期间发生变化，请重试")
    _write_json_atomic(output_path, payload)
    return _summary(payload, output_path, reused_existing=False)


def run_expression_box_design(config: Any) -> dict[str, Any]:
    from src.protein_to_cds.artifacts import ArtifactTransaction
    root = Path(config.project_output_path).expanduser().resolve()
    with ArtifactTransaction([], lock_path=root / "protein_to_cds" / ".gc_optimization.lock"):
        return _run_expression_box_design(config)


__all__ = [
    "EXPRESSION_BOX_DESIGNS_FILENAME",
    "run_expression_box_design",
]
