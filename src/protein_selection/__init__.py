"""Protein Supply application package with lazy public exports."""

from __future__ import annotations

from importlib import import_module
from typing import Any


_PUBLIC_EXPORTS = {
    "AUXILIARY_PROTEIN_PIPELINE_VERSION": (
        "src.protein_selection.models",
        "AUXILIARY_PROTEIN_PIPELINE_VERSION",
    ),
    "AUXILIARY_PROTEIN_RESEARCH_SCHEMA_VERSION": (
        "src.protein_selection.models",
        "AUXILIARY_PROTEIN_RESEARCH_SCHEMA_VERSION",
    ),
    "AuxiliaryProteinCombinationResult": (
        "src.protein_selection.models",
        "AuxiliaryProteinCombinationResult",
    ),
    "MainEnzymeAuxiliaryResearchResult": (
        "src.protein_selection.models",
        "MainEnzymeAuxiliaryResearchResult",
    ),
    "auxiliary_protein_result_path": (
        "src.protein_selection.pipeline",
        "auxiliary_protein_result_path",
    ),
    "execute_research_units": (
        "src.protein_selection.pipeline",
        "execute_research_units",
    ),
    "run_auxiliary_protein_pipeline": (
        "src.protein_selection.pipeline",
        "run_auxiliary_protein_pipeline",
    ),
    "write_auxiliary_protein_result": (
        "src.protein_selection.pipeline",
        "write_auxiliary_protein_result",
    ),
}

__all__ = sorted(_PUBLIC_EXPORTS)


def __getattr__(name: str) -> Any:
    target = _PUBLIC_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attribute_name = target
    value = getattr(import_module(module_name), attribute_name)
    globals()[name] = value
    return value
