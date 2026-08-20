"""Read the current CDS-selection schema for expression-box grouping."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.expression_box.config import (
    BALANCED_MAX_MAIN_UNITS_PER_CASSETTE,
    COMPACT_MAX_CDS_LENGTH_NT,
    COMPACT_MAX_PROTEINS,
    GROUPING_ALGORITHM_VERSION,
)
from src.expression_box.models import ExpressionGroupingContext, ExpressionProtein
from src.write_manifest.store import read_design_manifest

_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_SUPPORTED_ROLES = {"main_enzyme", "auxiliary_protein"}


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest missing object: {field_name}")
    return value


def _nonempty_text(value: Any, field_name: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"manifest field must not be empty: {field_name}")
    return text


def _positive_integer(value: Any, field_name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"manifest field must be a positive integer: {field_name}")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"manifest field must be a positive integer: {field_name}"
        ) from exc
    if parsed < 1 or isinstance(value, float) and not value.is_integer():
        raise ValueError(f"manifest field must be a positive integer: {field_name}")
    return parsed


def _string_tuple(value: Any, field_name: str) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise ValueError(f"manifest field must be a list: {field_name}")
    return tuple(
        sorted({str(item or "").strip() for item in value if str(item or "").strip()})
    )


def _step_indexes(value: Any, field_name: str) -> tuple[int, ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"manifest field must be a non-empty list: {field_name}")
    return tuple(sorted({_positive_integer(item, field_name) for item in value}))


def _sequence_sha256(value: Any, field_name: str) -> str:
    normalized = str(value or "").strip().lower()
    if _SHA256_RE.fullmatch(normalized) is None:
        raise ValueError(f"manifest field is not a SHA-256 digest: {field_name}")
    return normalized


def _stable_fingerprint(payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _protein_from_manifest(raw: Any, index: int) -> ExpressionProtein:
    field = f"cds_selection.proteins[{index}]"
    item = _mapping(raw, field)
    accession = _nonempty_text(item.get("accession"), f"{field}.accession").upper()
    roles = _string_tuple(item.get("roles"), f"{field}.roles")
    if not roles or not set(roles).issubset(_SUPPORTED_ROLES):
        raise ValueError(
            f"manifest field contains unsupported expression roles: {field}.roles"
        )
    optimized_cds = _mapping(item.get("optimized_cds"), f"{field}.optimized_cds")
    return ExpressionProtein(
        accession=accession,
        roles=roles,
        assigned_step_indexes=_step_indexes(
            item.get("assigned_step_indexes"),
            f"{field}.assigned_step_indexes",
        ),
        required_by_main_accessions=tuple(
            value.upper()
            for value in _string_tuple(
                item.get("required_by_main_accessions", []),
                f"{field}.required_by_main_accessions",
            )
        ),
        optimized_cds_length_nt=_positive_integer(
            optimized_cds.get("length_nt"),
            f"{field}.optimized_cds.length_nt",
        ),
        optimized_cds_sequence_sha256=_sequence_sha256(
            optimized_cds.get("sequence_sha256"),
            f"{field}.optimized_cds.sequence_sha256",
        ),
    )


def load_expression_grouping_context(
    manifest_path: str | Path,
) -> ExpressionGroupingContext:
    """Return the minimal current-manifest input needed for grouping."""

    path = Path(manifest_path).expanduser().resolve()
    manifest = read_design_manifest(path)
    target_compound_id = _nonempty_text(
        manifest.get("target_compound_id"),
        "target_compound_id",
    )
    cds_selection = _mapping(manifest.get("cds_selection"), "cds_selection")
    if cds_selection.get("status") != "complete":
        raise ValueError("cds_selection must be complete before expression-box design")
    raw_proteins = cds_selection.get("proteins")
    if not isinstance(raw_proteins, list) or not raw_proteins:
        raise ValueError("cds_selection.proteins must be a non-empty list")

    proteins = tuple(
        sorted(
            (
                _protein_from_manifest(raw, index)
                for index, raw in enumerate(raw_proteins)
            ),
            key=lambda item: item.accession,
        )
    )
    accessions = [item.accession for item in proteins]
    if len(accessions) != len(set(accessions)):
        raise ValueError("cds_selection.proteins contains duplicate accessions")
    main_accessions = {item.accession for item in proteins if item.is_main_enzyme}
    if not main_accessions:
        raise ValueError("cds_selection.proteins contains no main enzymes")
    for protein in proteins:
        missing_dependencies = sorted(
            set(protein.required_by_main_accessions) - main_accessions
        )
        if missing_dependencies:
            raise ValueError(
                f"auxiliary protein {protein.accession} references unknown main enzymes: "
                + ", ".join(missing_dependencies)
            )

    try:
        revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision must be an integer") from exc
    cds_source_fingerprint = str(cds_selection.get("source_fingerprint") or "").strip()
    fingerprint_payload = {
        "algorithm_version": GROUPING_ALGORITHM_VERSION,
        "constraints": {
            "balanced_max_main_units_per_cassette": (
                BALANCED_MAX_MAIN_UNITS_PER_CASSETTE
            ),
            "compact_max_proteins": COMPACT_MAX_PROTEINS,
            "compact_max_cds_length_nt": COMPACT_MAX_CDS_LENGTH_NT,
        },
        "target_compound_id": target_compound_id,
        "cds_selection_source_fingerprint": cds_source_fingerprint,
        "proteins": [
            {
                "accession": item.accession,
                "roles": item.roles,
                "assigned_step_indexes": item.assigned_step_indexes,
                "required_by_main_accessions": (item.required_by_main_accessions),
                "optimized_cds_length_nt": item.optimized_cds_length_nt,
                "optimized_cds_sequence_sha256": (item.optimized_cds_sequence_sha256),
            }
            for item in proteins
        ],
    }
    return ExpressionGroupingContext(
        manifest_path=path,
        manifest_revision=revision,
        target_compound_id=target_compound_id,
        cds_selection_source_fingerprint=cds_source_fingerprint,
        input_fingerprint=_stable_fingerprint(fingerprint_payload),
        proteins=proteins,
    )


__all__ = ["load_expression_grouping_context"]
