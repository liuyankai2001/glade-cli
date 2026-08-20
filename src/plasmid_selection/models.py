"""Small immutable models shared by the plasmid-selection pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class ExpressionBurden:
    schema_version: str
    model_version: str
    score: float
    level: str
    confidence: str
    raw_load_units: float
    reference_load_units: float
    gene_count: int
    cassette_count: int
    total_cds_length_nt: int
    minimum_ostir_reference_count: int
    fingerprint: str
    warnings: tuple[str, ...]

    def fingerprint_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "model_version": self.model_version,
            "score": self.score,
            "level": self.level,
            "confidence": self.confidence,
            "raw_load_units": self.raw_load_units,
            "reference_load_units": self.reference_load_units,
            "gene_count": self.gene_count,
            "cassette_count": self.cassette_count,
            "total_cds_length_nt": self.total_cds_length_nt,
            "minimum_ostir_reference_count": (
                self.minimum_ostir_reference_count
            ),
            "fingerprint": self.fingerprint,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True)
class ExpressionConstruct:
    design_id: int
    rank: int
    length_bp: int
    cassette_count: int
    component_count: int
    sequence_sha256: str
    file_sha256: str
    path: Path
    burden: ExpressionBurden

    def fingerprint_payload(self) -> dict[str, Any]:
        return {
            "design_id": self.design_id,
            "rank": self.rank,
            "length_bp": self.length_bp,
            "cassette_count": self.cassette_count,
            "component_count": self.component_count,
            "sequence_sha256": self.sequence_sha256,
            "file_sha256": self.file_sha256,
            "expression_burden": self.burden.fingerprint_payload(),
        }


@dataclass(frozen=True, slots=True)
class PlasmidContext:
    target_compound_id: str
    manifest_path: Path
    project_output_path: Path
    manifest_revision: int
    host_key: str
    host_name: str
    parts_selection_fingerprint: str
    assembled_constructs_fingerprint: str
    input_fingerprint: str
    selected_design_ids: tuple[int, ...]
    constructs: tuple[ExpressionConstruct, ...]

    @property
    def minimum_insert_length_bp(self) -> int:
        return min(item.length_bp for item in self.constructs)

    @property
    def maximum_insert_length_bp(self) -> int:
        return max(item.length_bp for item in self.constructs)

    @property
    def maximum_cassette_count(self) -> int:
        return max(item.cassette_count for item in self.constructs)

    @property
    def minimum_burden_score(self) -> float:
        return min(item.burden.score for item in self.constructs)

    @property
    def maximum_burden_score(self) -> float:
        return max(item.burden.score for item in self.constructs)


@dataclass(frozen=True, slots=True)
class PlasmidSnapshot:
    collection_name: str
    collection_row_count: int
    schema_fingerprint: str
    candidate_fingerprint: str
    candidates: tuple[dict[str, Any], ...]


@dataclass(frozen=True, slots=True)
class DownloadedSequence:
    content: bytes
    length_bp: int
    sequence_content_sha256: str
    canonical_sequence_sha256: str
    file_sha256: str
    source_download_url: str


__all__ = [
    "DownloadedSequence",
    "ExpressionBurden",
    "ExpressionConstruct",
    "PlasmidContext",
    "PlasmidSnapshot",
]
