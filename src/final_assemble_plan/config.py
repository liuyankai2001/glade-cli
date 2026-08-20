"""Configuration for final-assembly plan recommendation."""

from __future__ import annotations


ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION = (
    "assembly_plan_recommendations.v1"
)
FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION = "final_assembly_plan.v2"
ASSEMBLY_PLAN_ALGORITHM_VERSION = "final_assembly_planning.v1.0.0"
ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME = "assembly_plan_recommendations.json"
SUPPORTED_METHODS = ("restriction", "gibson")

DEFAULT_HOMOLOGY_ARM_LENGTH = 30
MIN_HOMOLOGY_ARM_GC_PERCENT = 30.0
MAX_HOMOLOGY_ARM_GC_PERCENT = 70.0
MAX_HOMOLOGY_ARM_HOMOPOLYMER = 7

CURATED_RESTRICTION_ENZYMES = (
    "EcoRI",
    "BamHI",
    "HindIII",
    "XhoI",
    "XbaI",
    "SpeI",
    "NheI",
    "PstI",
    "KpnI",
    "SacI",
    "SalI",
    "NotI",
    "SphI",
    "BglII",
    "AgeI",
    "AvrII",
    "BstBI",
    "ClaI",
    "MluI",
    "NcoI",
    "NdeI",
    "PacI",
    "AscI",
    "SbfI",
    "SmaI",
    "ApaI",
    "EagI",
)

ENZYME_PRIORITIES = {
    "EcoRI": 95,
    "BamHI": 95,
    "HindIII": 95,
    "XhoI": 92,
    "XbaI": 90,
    "PstI": 90,
    "SpeI": 88,
    "NheI": 88,
    "KpnI": 88,
    "SacI": 88,
    "SalI": 86,
    "NotI": 86,
    "SphI": 84,
    "BglII": 84,
    "AgeI": 82,
    "AvrII": 82,
    "BstBI": 80,
    "ClaI": 80,
    "MluI": 80,
    "NcoI": 78,
    "NdeI": 78,
    "PacI": 76,
    "AscI": 76,
    "SbfI": 76,
    "SmaI": 72,
    "ApaI": 72,
    "EagI": 72,
}

PLAN_SCORE_WEIGHTS = {
    "insertion_region_safety": 40.0,
    "junction_quality": 25.0,
    "experimental_simplicity": 20.0,
    "source_integrity": 10.0,
    "estimated_final_size": 5.0,
}


__all__ = [
    "ASSEMBLY_PLAN_ALGORITHM_VERSION",
    "ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME",
    "ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION",
    "CURATED_RESTRICTION_ENZYMES",
    "DEFAULT_HOMOLOGY_ARM_LENGTH",
    "ENZYME_PRIORITIES",
    "FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION",
    "MAX_HOMOLOGY_ARM_GC_PERCENT",
    "MAX_HOMOLOGY_ARM_HOMOPOLYMER",
    "MIN_HOMOLOGY_ARM_GC_PERCENT",
    "PLAN_SCORE_WEIGHTS",
    "SUPPORTED_METHODS",
]
