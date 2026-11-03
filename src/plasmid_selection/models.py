"""Small immutable models shared by the plasmid-selection pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
from pathlib import Path
from typing import Any


def without_audit_hashes(value: Any) -> Any:
    """Remove historical digest metadata while retaining actual source data."""
    if isinstance(value, Mapping):
        return {key: without_audit_hashes(item) for key, item in value.items()
                if "fingerprint" not in str(key).lower()
                and not str(key).lower().endswith("sha256")}
    if isinstance(value, (list, tuple)):
        return [without_audit_hashes(item) for item in value]
    return value


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
    warnings: tuple[str, ...]

    def to_payload(self) -> dict[str, Any]:
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
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True, slots=True)
class ExpressionConstruct:
    design_id: int
    rank: int
    length_bp: int
    cassette_count: int
    component_count: int
    path: Path
    burden: ExpressionBurden


@dataclass(frozen=True, slots=True)
class PlasmidContext:
    target_compound_id: str
    manifest_path: Path
    project_output_path: Path
    manifest_revision: int
    host_key: str
    host_name: str
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
    candidates: tuple[dict[str, Any], ...]


@dataclass(frozen=True, slots=True)
class DownloadedSequence:
    content: bytes
    length_bp: int
    source_download_url: str


__all__ = [
    "DownloadedSequence",
    "ExpressionBurden",
    "ExpressionConstruct",
    "PlasmidContext",
    "PlasmidSnapshot",
]
