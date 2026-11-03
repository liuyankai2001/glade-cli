"""Typed records used by expression-part recommendation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ExpressionCds:
    accession: str
    sequence: str
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
    host_name: str
    host_key: str
    host_labels: tuple[str, ...]
    cassettes: tuple[ExpressionPartsCassette, ...]
    restriction_enzymes: tuple[str, ...] = ()
    homopolymer_max: int | None = None


@dataclass(frozen=True, slots=True)
class PartCandidate:
    part_id: str
    role: str
    sequence: str
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

@dataclass(frozen=True, slots=True)
class PartsSnapshot:
    collection_name: str
    collection_row_count: int
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


__all__ = [
    "ExpressionCds",
    "ExpressionPartsCassette",
    "ExpressionPartsContext",
    "PartCandidate",
    "PartsSnapshot",
    "RbsPrediction",
]
