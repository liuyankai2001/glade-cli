"""Immutable models for final-assembly planning."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class AssemblyConstruct:
    design_id: int
    rank: int
    path: Path
    length_bp: int
    sequence: str


@dataclass(frozen=True, slots=True)
class AssemblyBackbone:
    plasmid_id: str
    name: str
    path: Path
    length_bp: int
    sequence: str
    topology: str
    assembly_policy: str
    insertion_regions: tuple[dict[str, Any], ...]
    protected_features: tuple[dict[str, Any], ...]
    genbank_features: tuple[dict[str, Any], ...]



@dataclass(frozen=True, slots=True)
class FinalAssemblyContext:
    target_compound_id: str
    manifest_path: Path
    project_output_path: Path
    constructs: tuple[AssemblyConstruct, ...]
    backbone: AssemblyBackbone


__all__ = [
    "AssemblyBackbone",
    "AssemblyConstruct",
    "FinalAssemblyContext",
]
