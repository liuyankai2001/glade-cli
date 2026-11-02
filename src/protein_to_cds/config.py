"""Configuration owned by the manifest-driven protein-to-CDS workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from src.main_protein_selection.host_registry import get_host_preset
from src.protein_to_cds.organisms import resolve_organism

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CODON_TRANSFORMER_MODEL_DIR = PROJECT_ROOT / "model" / "CodonTransformer"


@dataclass(frozen=True, slots=True)
class HostProfile:
    chassis_key: str
    host_name: str
    codon_transformer_organism_id: int
    taxon_id: int | None = None
    codon_transformer_organism_name: str = ''

    def to_dict(self) -> dict:
        return {
            'chassis_key': self.chassis_key,
            'name': self.host_name,
            'taxon_id': self.taxon_id,
            'codon_transformer_organism_id': self.codon_transformer_organism_id,
            'codon_transformer_organism_name': self.codon_transformer_organism_name,
        }


HOST_PROFILES = {
    "ecoli_mg1655": HostProfile(
        chassis_key="ecoli_mg1655",
        host_name="Escherichia coli str. K-12 substr. MG1655",
        codon_transformer_organism_id=52,
        taxon_id=511145,
        codon_transformer_organism_name="Escherichia coli str. K-12 substr. MG1655",
    ),
}


@dataclass(frozen=True, slots=True)
class CdsConstraintConfig:
    policy_version: str = "cds_constraints.mg1655.v1"
    local_gc_window_nt: int = 50
    local_gc_min: float = 0.20
    local_gc_max: float = 0.80
    homopolymer_limit: int = 7
    rare_codon_weight_lt: float = 0.10
    rare_cluster_window_codons: int = 10
    rare_cluster_min_codons: int = 3


CDS_CONSTRAINT_CONFIG = CdsConstraintConfig()


def host_profile_for_chassis(
    chassis_key: str, *, codon_transformer_organism: int | None = None,
) -> HostProfile:
    normalized = str(chassis_key or "").strip()
    preset = get_host_preset(normalized)
    if codon_transformer_organism is None:
        if preset['taxon_id'] != 511145:
            raise ValueError(
                'codon_transformer_organism is required for this chassis; '
                f'unsupported chassis_key without an explicit generation category: {normalized!r}'
            )
        codon_transformer_organism = 52
    identifier, name = resolve_organism(codon_transformer_organism)
    return HostProfile(
        chassis_key=normalized, host_name=preset['name'], taxon_id=preset['taxon_id'],
        codon_transformer_organism_id=identifier, codon_transformer_organism_name=name,
    )


__all__ = [
    "CDS_CONSTRAINT_CONFIG",
    "CODON_TRANSFORMER_MODEL_DIR",
    "HOST_PROFILES",
    "CdsConstraintConfig",
    "HostProfile",
    "host_profile_for_chassis",
]
