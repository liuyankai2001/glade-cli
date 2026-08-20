"""Typed records used by expression-part recommendation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class ExpressionCds:
    accession: str
    sequence: str
    sequence_sha256: str
    path: Path


@dataclass(frozen=True, slots=True)
class ExpressionPartsCassette:
    cassette_index: int
    cds: tuple[ExpressionCds, ...]


@dataclass(frozen=True, slots=True)
class ExpressionPartsContext:
    manifest_path: Path
    manifest_revision: int
    project_root: Path
    target_compound_id: str
    expression_box_selection_fingerprint: str
    cds_selection_source_fingerprint: str
    host_name: str
    host_key: str
    host_labels: tuple[str, ...]
    input_fingerprint: str
    cassettes: tuple[ExpressionPartsCassette, ...]


@dataclass(frozen=True, slots=True)
class PartCandidate:
    part_id: str
    role: str
    sequence: str
    sequence_sha256: str
    sequence_type: str
    host_match_kind: str
    strength: str
    regulation: str
    direction: str
    regulator_required: str
    evidence_grade: str
    role_confidence: str
    sequence_type_confidence: str
    host_confidence: str
    strength_confidence: str
    regulation_confidence: str
    activity_value: float | None
    activity_percentile: float | None
    activity_dataset: str
    activity_unit: str
    activity_context: str
    activity_source: str
    source: str
    evidence: str
    warnings: tuple[str, ...]
    updated_at: str

    def fingerprint_payload(self) -> dict[str, Any]:
        return {
            "part_id": self.part_id,
            "role": self.role,
            "sequence_sha256": self.sequence_sha256,
            "sequence_type": self.sequence_type,
            "host_match_kind": self.host_match_kind,
            "strength": self.strength,
            "regulation": self.regulation,
            "direction": self.direction,
            "regulator_required": self.regulator_required,
            "evidence_grade": self.evidence_grade,
            "role_confidence": self.role_confidence,
            "sequence_type_confidence": self.sequence_type_confidence,
            "host_confidence": self.host_confidence,
            "strength_confidence": self.strength_confidence,
            "regulation_confidence": self.regulation_confidence,
            "activity_value": self.activity_value,
            "activity_percentile": self.activity_percentile,
            "activity_dataset": self.activity_dataset,
            "activity_unit": self.activity_unit,
            "activity_context": self.activity_context,
            "activity_source": self.activity_source,
            "source": self.source,
            "evidence": self.evidence,
            "warnings": list(self.warnings),
            "updated_at": self.updated_at,
        }


@dataclass(frozen=True, slots=True)
class PartsSnapshot:
    collection_name: str
    collection_row_count: int
    schema_fingerprint: str
    candidate_fingerprint: str
    candidates: tuple[PartCandidate, ...]
    raw_counts: dict[str, int]
    accepted_counts: dict[str, int]
    rejected_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class RbsPrediction:
    part_id: str
    accession: str
    expression: float
    d_g_total: float
    intended_start_position: int
    unintended_start_count: int
    context_sha256: str


__all__ = [
    "ExpressionCds",
    "ExpressionPartsCassette",
    "ExpressionPartsContext",
    "PartCandidate",
    "PartsSnapshot",
    "RbsPrediction",
]
