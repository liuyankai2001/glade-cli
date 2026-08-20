"""Configuration for manifest-driven plasmid backbone recommendation."""

from __future__ import annotations


PLASMID_COLLECTION = "plasmid_templates_v2"
PLASMID_CANDIDATES_SCHEMA_VERSION = "plasmid_candidates.v1"
PLASMID_SELECTION_SCHEMA_VERSION = "plasmid_selection.v1"
PLASMID_RECOMMENDATION_ALGORITHM_VERSION = "plasmid_recommendation.v1.0.0"
PLASMID_CANDIDATES_FILENAME = "plasmid_candidates.json"
DEFAULT_CANDIDATE_COUNT = 5
MIN_CANDIDATE_COUNT = 1
MAX_CANDIDATE_COUNT = 20
MILVUS_QUERY_LIMIT = 1_000
MILVUS_TIMEOUT_SECONDS = 30
DOWNLOAD_TIMEOUT_SECONDS = 30

SUPPORTED_PRIORITIES = ("stability", "balanced", "expression")
SUPPORTED_ASSEMBLY_POLICIES = {
    "insert_into_mcs",
    "replace_seva_cargo_paci_spei",
}

# The score is deliberately transparent and is not an experimental-success
# probability. Each row sums to the 35-point copy/load component.
COPY_CLASS_SCORES = {
    "stability": {
        "low": 35.0,
        "medium": 28.0,
        "context_dependent": 14.0,
        "high": 8.0,
    },
    "balanced": {
        "medium": 35.0,
        "low": 31.0,
        "context_dependent": 18.0,
        "high": 14.0,
    },
    "expression": {
        "high": 35.0,
        "medium": 30.0,
        "low": 18.0,
        "context_dependent": 14.0,
    },
}

MARKER_SCORES = {
    "kanamycin": 10.0,
    "chloramphenicol": 8.0,
    "gentamicin": 7.0,
    "streptomycin/spectinomycin": 6.0,
    "tetracycline": 5.0,
    "ampicillin": 3.0,
}

MARKER_ALIASES = {
    "kan": "kanamycin",
    "km": "kanamycin",
    "kanamycin": "kanamycin",
    "neo": "kanamycin",
    "chloramphenicol": "chloramphenicol",
    "cm": "chloramphenicol",
    "cat": "chloramphenicol",
    "gentamicin": "gentamicin",
    "gm": "gentamicin",
    "aacc1": "gentamicin",
    "streptomycin": "streptomycin/spectinomycin",
    "spectinomycin": "streptomycin/spectinomycin",
    "streptomycin/spectinomycin": "streptomycin/spectinomycin",
    "sm/sp": "streptomycin/spectinomycin",
    "sp/sm": "streptomycin/spectinomycin",
    "aada": "streptomycin/spectinomycin",
    "tetracycline": "tetracycline",
    "tet": "tetracycline",
    "tc": "tetracycline",
    "teta": "tetracycline",
    "ampicillin": "ampicillin",
    "amp": "ampicillin",
    "ap": "ampicillin",
    "bla": "ampicillin",
}


__all__ = [
    "COPY_CLASS_SCORES",
    "DEFAULT_CANDIDATE_COUNT",
    "DOWNLOAD_TIMEOUT_SECONDS",
    "MARKER_ALIASES",
    "MARKER_SCORES",
    "MAX_CANDIDATE_COUNT",
    "MILVUS_QUERY_LIMIT",
    "MILVUS_TIMEOUT_SECONDS",
    "MIN_CANDIDATE_COUNT",
    "PLASMID_CANDIDATES_FILENAME",
    "PLASMID_CANDIDATES_SCHEMA_VERSION",
    "PLASMID_COLLECTION",
    "PLASMID_RECOMMENDATION_ALGORITHM_VERSION",
    "PLASMID_SELECTION_SCHEMA_VERSION",
    "SUPPORTED_ASSEMBLY_POLICIES",
    "SUPPORTED_PRIORITIES",
]
