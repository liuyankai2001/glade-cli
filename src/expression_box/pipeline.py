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
        "required_by_main_accessions": list(protein.required_by_main_accessions),
        "optimized_cds_length_nt": protein.optimized_cds_length_nt,
        "optimized_cds_sequence_sha256": (protein.optimized_cds_sequence_sha256),
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
            "cds_selection_source_fingerprint": (
                context.cds_selection_source_fingerprint
            ),
            "input_fingerprint": context.input_fingerprint,
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


def _load_cached_payload(
    path: Path,
    input_fingerprint: str,
) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    source = payload.get("source") if isinstance(payload, dict) else None
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION
        or payload.get("status") != "complete"
        or not isinstance(source, dict)
        or source.get("input_fingerprint") != input_fingerprint
        or not isinstance(payload.get("designs"), list)
        or not payload["designs"]
    ):
        return None
    return payload


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
        "manifest_modified": False,
    }


def run_expression_box_design(config: Any) -> dict[str, Any]:
    """Generate and persist system-recommended protein grouping candidates."""

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    context = load_expression_grouping_context(manifest_path)
    output_path = project_root / "expression_box" / EXPRESSION_BOX_DESIGNS_FILENAME
    cached = _load_cached_payload(output_path, context.input_fingerprint)
    if cached is not None:
        try:
            return _summary(cached, output_path, reused_existing=True)
        except (KeyError, TypeError, ValueError):
            pass

    generation = generate_grouping_designs(context.proteins)
    payload = _full_payload(context, generation)
    _write_json_atomic(output_path, payload)
    return _summary(payload, output_path, reused_existing=False)


__all__ = [
    "EXPRESSION_BOX_DESIGNS_FILENAME",
    "run_expression_box_design",
]
