"""Unified writers for the project-level design manifest."""

from src.write_manifest.store import (
    SCHEMA_VERSION,
    read_design_manifest,
    update_design_manifest,
)
from src.write_manifest.solution import run_write_solution, write_solution
from src.write_manifest.expression_parts import (
    run_write_expression_parts_selection,
    write_expression_parts_selection,
)

__all__ = [
    "SCHEMA_VERSION",
    "read_design_manifest",
    "run_write_expression_parts_selection",
    "run_write_solution",
    "update_design_manifest",
    "write_expression_parts_selection",
    "write_solution",
]
"""Native writers for sections of ``design_manifest.json``."""
