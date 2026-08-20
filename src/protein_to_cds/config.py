"""Configuration owned by the manifest-driven protein-to-CDS workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CODON_TRANSFORMER_MODEL_DIR = PROJECT_ROOT / "model" / "CodonTransformer"


@dataclass(frozen=True, slots=True)
class HostProfile:
    chassis_key: str
    host_name: str
    codon_transformer_organism_id: int


HOST_PROFILES = {
    "ecoli_mg1655": HostProfile(
        chassis_key="ecoli_mg1655",
        host_name="Escherichia coli str. K-12 substr. MG1655",
        codon_transformer_organism_id=52,
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


def host_profile_for_chassis(chassis_key: str) -> HostProfile:
    normalized = str(chassis_key or "").strip()
    try:
        return HOST_PROFILES[normalized]
    except KeyError as exc:
        supported = ", ".join(sorted(HOST_PROFILES))
        raise ValueError(
            f"unsupported chassis_key for CDS optimization: {normalized!r}; "
            f"supported: {supported}"
        ) from exc


__all__ = [
    "CDS_CONSTRAINT_CONFIG",
    "CODON_TRANSFORMER_MODEL_DIR",
    "HOST_PROFILES",
    "CdsConstraintConfig",
    "HostProfile",
    "host_profile_for_chassis",
]
