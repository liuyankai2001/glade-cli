"""Typed models for expression-box protein grouping."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ExpressionProtein:
    accession: str
    roles: tuple[str, ...]
    assigned_step_indexes: tuple[int, ...]
    required_by_main_accessions: tuple[str, ...]
    optimized_cds_length_nt: int
    optimized_cds_sequence_sha256: str

    @property
    def is_main_enzyme(self) -> bool:
        return "main_enzyme" in self.roles

    @property
    def is_auxiliary_protein(self) -> bool:
        return "auxiliary_protein" in self.roles

    @property
    def first_step_index(self) -> int:
        return min(self.assigned_step_indexes)


@dataclass(frozen=True, slots=True)
class ExpressionGroupingContext:
    manifest_path: Path
    manifest_revision: int
    target_compound_id: str
    cds_selection_source_fingerprint: str
    input_fingerprint: str
    proteins: tuple[ExpressionProtein, ...]


@dataclass(frozen=True, slots=True)
class ProteinUnit:
    main_accession: str
    proteins: tuple[ExpressionProtein, ...]
    first_step_index: int

    @property
    def total_cds_length_nt(self) -> int:
        return sum(item.optimized_cds_length_nt for item in self.proteins)


@dataclass(frozen=True, slots=True)
class ExpressionCassette:
    proteins: tuple[ExpressionProtein, ...]
    reason: str

    @property
    def total_cds_length_nt(self) -> int:
        return sum(item.optimized_cds_length_nt for item in self.proteins)

    @property
    def main_enzyme_count(self) -> int:
        return sum(item.is_main_enzyme for item in self.proteins)

    @property
    def signature(self) -> tuple[str, ...]:
        return tuple(item.accession for item in self.proteins)


@dataclass(frozen=True, slots=True)
class ExpressionGroupingDesign:
    strategy: str
    name: str
    recommended: bool
    cassettes: tuple[ExpressionCassette, ...]
    warnings: tuple[str, ...] = ()

    @property
    def signature(self) -> tuple[tuple[str, ...], ...]:
        return tuple(cassette.signature for cassette in self.cassettes)


@dataclass(frozen=True, slots=True)
class SkippedGroupingStrategy:
    strategy: str
    reason: str


@dataclass(frozen=True, slots=True)
class GroupingGenerationResult:
    designs: tuple[ExpressionGroupingDesign, ...]
    skipped_strategies: tuple[SkippedGroupingStrategy, ...]
    warnings: tuple[str, ...]


__all__ = [
    "ExpressionCassette",
    "ExpressionGroupingContext",
    "ExpressionGroupingDesign",
    "ExpressionProtein",
    "GroupingGenerationResult",
    "ProteinUnit",
    "SkippedGroupingStrategy",
]
