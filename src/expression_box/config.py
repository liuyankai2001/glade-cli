"""Configuration for deterministic expression-box grouping."""

from __future__ import annotations

EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION = "expression_box_designs.v1"
GROUPING_ALGORITHM_VERSION = "expression_box_grouping.v1"

BALANCED_MAX_MAIN_UNITS_PER_CASSETTE = 2
COMPACT_MAX_PROTEINS = 6
COMPACT_MAX_CDS_LENGTH_NT = 12_000


__all__ = [
    "BALANCED_MAX_MAIN_UNITS_PER_CASSETTE",
    "COMPACT_MAX_CDS_LENGTH_NT",
    "COMPACT_MAX_PROTEINS",
    "EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION",
    "GROUPING_ALGORITHM_VERSION",
]
