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
    sequence_sha256: str
    file_sha256: str

    def fingerprint_payload(self) -> dict[str, Any]:
        return {
            "parts_design_id": self.design_id,
            "rank": self.rank,
            "length_bp": self.length_bp,
            "sequence_sha256": self.sequence_sha256,
            "file_sha256": self.file_sha256,
        }


@dataclass(frozen=True, slots=True)
class AssemblyBackbone:
    plasmid_id: str
    name: str
    path: Path
    length_bp: int
    sequence: str
    file_sha256: str
    sequence_content_sha256: str
    topology: str
    assembly_policy: str
    insertion_regions: tuple[dict[str, Any], ...]
    protected_features: tuple[dict[str, Any], ...]
    genbank_features: tuple[dict[str, Any], ...]

    def fingerprint_payload(self) -> dict[str, Any]:
        return {
            "plasmid_id": self.plasmid_id,
            "name": self.name,
            "length_bp": self.length_bp,
            "file_sha256": self.file_sha256,
            "sequence_content_sha256": self.sequence_content_sha256,
            "topology": self.topology,
            "assembly_policy": self.assembly_policy,
            "insertion_regions": list(self.insertion_regions),
            "protected_features": list(self.protected_features),
        }


@dataclass(frozen=True, slots=True)
class FinalAssemblyContext:
    target_compound_id: str
    manifest_path: Path
    project_output_path: Path
    manifest_revision: int
    parts_selection_fingerprint: str
    assembled_constructs_fingerprint: str
    plasmid_selection_fingerprint: str
    input_fingerprint: str
    constructs: tuple[AssemblyConstruct, ...]
    backbone: AssemblyBackbone


__all__ = [
    "AssemblyBackbone",
    "AssemblyConstruct",
    "FinalAssemblyContext",
]
