"""Small immutable models shared by final-assembly execution modules."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.SeqRecord import SeqRecord

from src.final_assemble_plan.models import FinalAssemblyContext


@dataclass(frozen=True, slots=True)
class FinalAssemblyExecutionContext:
    manifest: dict[str, Any]
    manifest_path: Path
    project_output_path: Path
    manifest_revision: int
    source_context: FinalAssemblyContext
    plan_selection_fingerprint: str
    plans: tuple[dict[str, Any], ...]
    backbone_record: SeqRecord
    insert_records: dict[int, SeqRecord]


@dataclass(frozen=True, slots=True)
class SequenceAssemblyResult:
    sequence: str
    inserted_start_bp: int
    inserted_end_bp: int
    vector_replace_start_index: int
    vector_replace_end_index: int
    replacement_payload_length_bp: int
    target_audit: dict[str, Any]


__all__ = [
    "FinalAssemblyExecutionContext",
    "SequenceAssemblyResult",
]
