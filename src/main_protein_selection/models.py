"""Validated public output models for main-enzyme selection."""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


MAIN_ENZYME_SELECTION_SCHEMA_VERSION = "main_enzyme_selection.v1"
AcceptedReactionFitStatus = Literal["verified", "verified_with_risk"]
SelectionStatus = Literal["complete", "source_unavailable"]


def _string(value: Any) -> str:
    return str(value or "").strip()


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _optional_int(value: Any) -> int | None:
    text = _string(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return _string(value).lower() in {"1", "true", "yes"}


def _values(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        raw_values = [str(item).strip() for item in value]
    else:
        raw_values = [
            item.strip()
            for item in re.split(r"\s*[;|]\s*", _string(value))
        ]
    return list(dict.fromkeys(item for item in raw_values if item))


class MainEnzymeCandidate(BaseModel):
    """One reaction-verified candidate for one selected-route step."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    ec_number: str | None = None
    accession: str = Field(min_length=1)
    protein_name: str = ""
    organism_name: str = ""
    reviewed: bool
    length: int | None = Field(default=None, ge=1)
    candidate_rank: int = Field(ge=1)
    protein_score: float = Field(ge=0.0)
    reaction_fit_status: AcceptedReactionFitStatus
    reaction_fit_score: float = Field(ge=0.0)
    direction_verdict: str = ""
    direction_confidence: str = ""
    retrieval_strategies: list[str] = Field(default_factory=list)
    matched_rhea_ids: list[str] = Field(default_factory=list)
    matched_ko_ids: list[str] = Field(default_factory=list)
    reaction_confidence: str = ""
    sequence_version: int | None = Field(default=None, ge=1)
    sequence_sha256: str = ""
    sequence: str = ""
    warnings: list[str] = Field(default_factory=list)
    reasons: list[str] = Field(default_factory=list)

    @classmethod
    def from_candidate_row(
        cls,
        row: Mapping[str, Any],
    ) -> "MainEnzymeCandidate":
        """Validate and project one in-memory candidate table row."""

        return cls(
            step_index=int(float(_string(row.get("step_index")) or 0)),
            reaction_id=_string(row.get("reaction_id")),
            ec_number=_string(row.get("ec_number")) or None,
            accession=_string(row.get("accession")).upper(),
            protein_name=_string(row.get("protein_name")),
            organism_name=_string(row.get("organism_name")),
            reviewed=_bool(row.get("reviewed")),
            length=_optional_int(row.get("length")),
            candidate_rank=int(float(_string(row.get("candidate_rank")) or 0)),
            protein_score=_float(row.get("score")),
            reaction_fit_status=_string(row.get("reaction_fit_status")),
            reaction_fit_score=_float(row.get("reaction_fit_score")),
            direction_verdict=_string(row.get("direction_verdict")),
            direction_confidence=_string(row.get("direction_confidence")),
            retrieval_strategies=_values(row.get("retrieval_strategy")),
            matched_rhea_ids=_values(row.get("matched_rhea_ids")),
            matched_ko_ids=_values(row.get("matched_ko_ids")),
            reaction_confidence=_string(row.get("reaction_confidence")),
            sequence_version=_optional_int(row.get("sequence_version")),
            sequence_sha256=_string(row.get("sequence_sha256")),
            sequence=_string(row.get("sequence")),
            warnings=_values(row.get("warnings")),
            reasons=_values(row.get("reasons")),
        )


class MainEnzymeSelectionParameters(BaseModel):
    """Selection controls that materially affect the candidate set."""

    model_config = ConfigDict(extra="forbid")

    top_n: int = Field(ge=1)
    max_results: int = Field(ge=1)
    allow_transmembrane: bool
    fetch_proteins: bool


class MainEnzymeSelectionResult(BaseModel):
    """Canonical machine-readable result for one selected route."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    schema_version: Literal["main_enzyme_selection.v1"] = (
        MAIN_ENZYME_SELECTION_SCHEMA_VERSION
    )
    ok: bool
    status: SelectionStatus
    selected_solution_id: int = Field(ge=1)
    expansion_depth: int = Field(ge=0)
    solution_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    chassis_key: str = Field(min_length=1)
    parameters: MainEnzymeSelectionParameters
    candidates_by_step: dict[int, list[MainEnzymeCandidate]]
    uncovered_step_indexes: list[int] = Field(default_factory=list)
    direction_rejected_step_indexes: list[int] = Field(default_factory=list)
    direction_risk_step_indexes: list[int] = Field(default_factory=list)
    evidence_files: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_consistency(self) -> "MainEnzymeSelectionResult":
        if self.ok != (self.status == "complete"):
            raise ValueError("ok must agree with selection status")
        for step_index, candidates in self.candidates_by_step.items():
            if step_index < 1:
                raise ValueError("candidate step indexes must be positive")
            if any(item.step_index != step_index for item in candidates):
                raise ValueError("candidate does not match its step group")
            ranks = [item.candidate_rank for item in candidates]
            if ranks != sorted(ranks) or len(ranks) != len(set(ranks)):
                raise ValueError("candidate ranks must be unique and sorted per step")
        return self


__all__ = [
    "MAIN_ENZYME_SELECTION_SCHEMA_VERSION",
    "MainEnzymeCandidate",
    "MainEnzymeSelectionParameters",
    "MainEnzymeSelectionResult",
]
