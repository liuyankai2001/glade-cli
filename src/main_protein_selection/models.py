"""Validated public output models for main-enzyme selection."""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from math import isclose
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


MAIN_ENZYME_SELECTION_SCHEMA_VERSION = "main_enzyme_selection.v3"
MAIN_ENZYME_SETS_SCHEMA_VERSION = "main_enzyme_sets.v3"
MAIN_ENZYME_SETS_ALGORITHM_VERSION = (
    "evidence_constrained_assignment.v3_taxonomy_ranked"
)
AcceptedReactionFitStatus = Literal[
    "verified",
    "verified_with_risk",
    "manual_review",
]
SelectionStatus = Literal["complete", "source_unavailable"]
MainEnzymeSetStatus = Literal["complete", "review_required"]
MainEnzymeSetsStatus = Literal[
    "complete",
    "review_required",
    "infeasible",
    "truncated",
    "stale_input",
    "source_unavailable",
]
AuxiliarySelectionStatus = Literal[
    "pending_user_selection",
    "integrated_in_main_enzyme",
]
AuxiliaryRequirementStatus = Literal[
    "not_required",
    "pending_user_selection",
    "integrated",
    "mixed",
]


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


class AuxiliaryRoleRequirement(BaseModel):
    """One still-unselected or main-enzyme-integrated helper role."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    role: str = Field(min_length=1)
    necessity: Literal["required", "possibly_required"]
    confidence: Literal["high", "medium", "low"]
    selection_status: AuxiliarySelectionStatus
    carrier_ids: list[str] = Field(default_factory=list)
    evidence: list[str] = Field(default_factory=list)
    step_indexes: list[int] = Field(default_factory=list)
    main_enzyme_accessions: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_stable_lists(self) -> "AuxiliaryRoleRequirement":
        if self.carrier_ids != sorted(set(self.carrier_ids)):
            raise ValueError("carrier_ids must be sorted and unique")
        if self.evidence != list(dict.fromkeys(self.evidence)):
            raise ValueError("evidence must be unique in stable order")
        if self.step_indexes != sorted(set(self.step_indexes)):
            raise ValueError("step_indexes must be sorted and unique")
        accessions = [item.upper() for item in self.main_enzyme_accessions]
        if accessions != sorted(set(accessions)):
            raise ValueError(
                "main_enzyme_accessions must be uppercase, sorted and unique"
            )
        return self


def _auxiliary_requirements(
    value: Any,
    *,
    step_index: int | None = None,
    accession: str | None = None,
) -> list[AuxiliaryRoleRequirement]:
    if isinstance(value, list):
        raw = value
    else:
        text = _string(value)
        if not text:
            raw = []
        else:
            try:
                raw = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError("invalid auxiliary requirements JSON") from exc
    if not isinstance(raw, list):
        raise ValueError("auxiliary requirements must be a list")
    result: list[AuxiliaryRoleRequirement] = []
    for item in raw:
        payload = dict(item)
        if step_index is not None:
            payload["step_indexes"] = sorted(set([
                *list(payload.get("step_indexes") or []),
                step_index,
            ]))
        if accession:
            payload["main_enzyme_accessions"] = sorted(set([
                *[str(value).upper() for value in payload.get(
                    "main_enzyme_accessions", []
                )],
                accession.upper(),
            ]))
        payload["carrier_ids"] = sorted(set(payload.get("carrier_ids") or []))
        payload["evidence"] = list(dict.fromkeys(payload.get("evidence") or []))
        result.append(AuxiliaryRoleRequirement.model_validate(payload))
    return sorted(result, key=lambda item: item.role)


class MainEnzymeCandidate(BaseModel):
    """One reaction-verified candidate for one selected-route step."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    ec_number: str | None = None
    accession: str = Field(min_length=1)
    protein_name: str = ""
    organism_name: str = ""
    organism_id: int | None = Field(default=None, ge=1)
    taxonomic_shared_taxon_id: int | None = Field(default=None, ge=1)
    taxonomic_shared_name: str = ""
    taxonomic_shared_rank: str = ""
    taxonomic_fit_status: str = Field(default="unknown", min_length=1)
    taxonomic_fit_score: float = Field(default=50.0, ge=0.0, le=100.0)
    taxonomy_evidence_source: str = ""
    reviewed: bool
    length: int | None = Field(default=None, ge=1)
    candidate_rank: int = Field(ge=1)
    protein_score: float = Field(ge=0.0)
    reaction_fit_status: AcceptedReactionFitStatus
    reaction_fit_score: float = Field(ge=0.0)
    direction_verdict: str = ""
    direction_confidence: str = ""
    retrieval_strategies: list[str] = Field(default_factory=list)
    retrieval_query_ids: list[str] = Field(default_factory=list)
    publication_ids: list[str] = Field(default_factory=list)
    matched_rhea_ids: list[str] = Field(default_factory=list)
    matched_ko_ids: list[str] = Field(default_factory=list)
    reaction_confidence: str = ""
    enzyme_system_type: str = ""
    auxiliary_requirement_status: AuxiliaryRequirementStatus = "not_required"
    auxiliary_requirements: list[AuxiliaryRoleRequirement] = Field(
        default_factory=list
    )
    sequence_version: int | None = Field(default=None, ge=1)
    sequence_sha256: str = ""
    sequence: str = ""
    warnings: list[str] = Field(default_factory=list)
    reasons: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_auxiliary_status(self) -> "MainEnzymeCandidate":
        statuses = {
            item.selection_status for item in self.auxiliary_requirements
        }
        expected: AuxiliaryRequirementStatus
        if not statuses:
            expected = "not_required"
        elif statuses == {"integrated_in_main_enzyme"}:
            expected = "integrated"
        elif statuses == {"pending_user_selection"}:
            expected = "pending_user_selection"
        else:
            expected = "mixed"
        if self.auxiliary_requirement_status != expected:
            raise ValueError(
                "auxiliary_requirement_status does not match requirements"
            )
        return self

    @classmethod
    def from_candidate_row(
        cls,
        row: Mapping[str, Any],
    ) -> "MainEnzymeCandidate":
        """Validate and project one in-memory candidate table row."""

        step_index = int(float(_string(row.get("step_index")) or 0))
        accession = _string(row.get("accession")).upper()
        return cls(
            step_index=step_index,
            reaction_id=_string(row.get("reaction_id")),
            ec_number=_string(row.get("ec_number")) or None,
            accession=accession,
            protein_name=_string(row.get("protein_name")),
            organism_name=_string(row.get("organism_name")),
            organism_id=_optional_int(row.get("organism_id")),
            taxonomic_shared_taxon_id=_optional_int(
                row.get("taxonomic_shared_taxon_id")
            ),
            taxonomic_shared_name=_string(row.get("taxonomic_shared_name")),
            taxonomic_shared_rank=_string(row.get("taxonomic_shared_rank")),
            taxonomic_fit_status=(
                _string(row.get("taxonomic_fit_status")) or "unknown"
            ),
            taxonomic_fit_score=_float(
                row.get("taxonomic_fit_score")
                if _string(row.get("taxonomic_fit_score"))
                else 50.0
            ),
            taxonomy_evidence_source=_string(
                row.get("taxonomy_evidence_source")
            ),
            reviewed=_bool(row.get("reviewed")),
            length=_optional_int(row.get("length")),
            candidate_rank=int(float(_string(row.get("candidate_rank")) or 0)),
            protein_score=_float(row.get("score")),
            reaction_fit_status=_string(row.get("reaction_fit_status")),
            reaction_fit_score=_float(row.get("reaction_fit_score")),
            direction_verdict=_string(row.get("direction_verdict")),
            direction_confidence=_string(row.get("direction_confidence")),
            retrieval_strategies=_values(row.get("retrieval_strategy")),
            retrieval_query_ids=_values(row.get("retrieval_query_id")),
            publication_ids=_values(row.get("publication_ids")),
            matched_rhea_ids=_values(row.get("matched_rhea_ids")),
            matched_ko_ids=_values(row.get("matched_ko_ids")),
            reaction_confidence=_string(row.get("reaction_confidence")),
            enzyme_system_type=_string(row.get("enzyme_system_type")),
            auxiliary_requirement_status=(
                _string(row.get("auxiliary_requirement_status"))
                or "not_required"
            ),
            auxiliary_requirements=_auxiliary_requirements(
                row.get("auxiliary_requirements_json"),
                step_index=step_index,
                accession=accession,
            ),
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
    literature_search: bool = False


class MainEnzymeSelectionResult(BaseModel):
    """Canonical machine-readable result for one selected route."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    schema_version: Literal["main_enzyme_selection.v3"] = (
        MAIN_ENZYME_SELECTION_SCHEMA_VERSION
    )
    ok: bool
    status: SelectionStatus
    selected_solution_id: int = Field(ge=1)
    expansion_depth: int = Field(ge=0)
    solution_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    chassis_key: str = Field(min_length=1)
    host_taxon_id: int = Field(default=511145, ge=1)
    taxonomy_status: Literal["resolved", "unknown"] = "resolved"
    taxonomy_source: str = Field(default="builtin_snapshot", min_length=1)
    taxonomy_scoring_policy_version: str = Field(
        default="taxonomy_lca_rank.v1",
        min_length=1,
    )
    taxonomy_fingerprint: str = Field(
        default="0" * 64,
        pattern=r"^[0-9a-f]{64}$",
    )
    scoring_weights: dict[str, float] = Field(default_factory=lambda: {
        "function": 0.40,
        "evidence": 0.25,
        "expression": 0.20,
        "host": 0.15,
    })
    parameters: MainEnzymeSelectionParameters
    shortlist_decision_fingerprint: str | None = Field(
        default=None,
        pattern=r"^[0-9a-f]{64}$",
    )
    candidates_by_step: dict[int, list[MainEnzymeCandidate]]
    uncovered_step_indexes: list[int] = Field(default_factory=list)
    direction_rejected_step_indexes: list[int] = Field(default_factory=list)
    direction_risk_step_indexes: list[int] = Field(default_factory=list)
    evidence_files: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_consistency(self) -> "MainEnzymeSelectionResult":
        expected_weights = {"function", "evidence", "expression", "host"}
        if set(self.scoring_weights) != expected_weights:
            raise ValueError("scoring_weights must define function/evidence/expression/host")
        if not isclose(sum(self.scoring_weights.values()), 1.0, abs_tol=1e-9):
            raise ValueError("scoring_weights must sum to 1")
        if self.ok != (self.status == "complete"):
            raise ValueError("ok must agree with selection status")
        for step_index, candidates in self.candidates_by_step.items():
            if step_index < 1:
                raise ValueError("candidate step indexes must be positive")
            if any(item.step_index != step_index for item in candidates):
                raise ValueError("candidate does not match its step group")
            ranks = [item.candidate_rank for item in candidates]
            if ranks != list(range(1, len(candidates) + 1)):
                raise ValueError(
                    "candidate ranks must be contiguous from 1 within each step"
                )
            if len(candidates) > self.parameters.top_n:
                raise ValueError("candidate count exceeds parameters.top_n")
        return self


def _require_sorted_unique_positive(
    values: list[int],
    field_name: str,
    *,
    allow_empty: bool = True,
) -> None:
    """Validate a deterministic list of positive step indexes."""

    if not allow_empty and not values:
        raise ValueError(f"{field_name} must not be empty")
    if any(value < 1 for value in values):
        raise ValueError(f"{field_name} must contain positive integers")
    if values != sorted(set(values)):
        raise ValueError(f"{field_name} must be sorted and unique")


class MainEnzymeSetParameters(BaseModel):
    """Controls that materially affect enzyme-set enumeration."""

    model_config = ConfigDict(extra="forbid")

    max_sets: int = Field(ge=1)
    max_search_nodes: int = Field(ge=1)
    candidate_scope: Literal["top_n_shortlist"] = "top_n_shortlist"


class MainEnzymeStepAssignment(BaseModel):
    """The one protein assigned as the primary catalyst for one step."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    accession: str = Field(min_length=1)
    candidate_rank: int = Field(ge=1)
    protein_score: float = Field(ge=0.0)
    host_fit_score: float = Field(ge=0.0)
    reaction_fit_status: AcceptedReactionFitStatus
    reaction_fit_score: float = Field(ge=0.0)
    direction_verdict: str = ""
    direction_confidence: str = ""
    specificity_status: str = Field(min_length=1)


class MainEnzymeSetProtein(BaseModel):
    """One distinct protein used by an enzyme set."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    accession: str = Field(min_length=1)
    protein_name: str = ""
    organism_name: str = ""
    reviewed: bool
    sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    cofactors: list[str] = Field(default_factory=list)
    capable_step_indexes: list[int]
    assigned_step_indexes: list[int]
    enzyme_system_types: list[str] = Field(default_factory=list)
    auxiliary_requirements: list[AuxiliaryRoleRequirement] = Field(
        default_factory=list
    )

    @model_validator(mode="after")
    def validate_consistency(self) -> "MainEnzymeSetProtein":
        _require_sorted_unique_positive(
            self.capable_step_indexes,
            "capable_step_indexes",
            allow_empty=False,
        )
        _require_sorted_unique_positive(
            self.assigned_step_indexes,
            "assigned_step_indexes",
            allow_empty=False,
        )
        if not set(self.assigned_step_indexes).issubset(
            self.capable_step_indexes
        ):
            raise ValueError(
                "assigned_step_indexes must be a subset of "
                "capable_step_indexes"
            )
        if any(not cofactor for cofactor in self.cofactors):
            raise ValueError("cofactors must not contain empty values")
        if self.cofactors != sorted(set(self.cofactors)):
            raise ValueError("cofactors must be sorted and unique")
        if self.enzyme_system_types != sorted(set(self.enzyme_system_types)):
            raise ValueError("enzyme_system_types must be sorted and unique")
        for requirement in self.auxiliary_requirements:
            if requirement.main_enzyme_accessions != [self.accession]:
                raise ValueError(
                    "protein auxiliary requirement must reference its accession"
                )
            if not set(requirement.step_indexes).issubset(
                self.assigned_step_indexes
            ):
                raise ValueError(
                    "protein auxiliary requirement references unassigned steps"
                )
        return self


class MainEnzymeSetMetrics(BaseModel):
    """Deterministic ranking and review metrics for one enzyme set."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    protein_count: int = Field(ge=1)
    organism_count: int = Field(ge=0)
    min_reaction_fit_score: float = Field(ge=0.0)
    mean_reaction_fit_score: float = Field(ge=0.0)
    min_protein_score: float = Field(ge=0.0)
    mean_protein_score: float = Field(ge=0.0)
    min_host_fit_score: float = Field(ge=0.0)
    mean_host_fit_score: float = Field(ge=0.0)
    reviewed_fraction: float = Field(ge=0.0, le=1.0)
    reaction_fit_risk_count: int = Field(ge=0)
    direction_risk_count: int = Field(ge=0)
    low_direction_confidence_count: int = Field(ge=0)
    specificity_risk_count: int = Field(ge=0)
    exact_specificity_count: int = Field(ge=0)
    warning_count: int = Field(ge=0)
    carrier_compatibility_status: str = Field(min_length=1)
    electron_reassessment_status: str = Field(min_length=1)
    auxiliary_requirement_status: AuxiliaryRequirementStatus
    pending_auxiliary_role_count: int = Field(ge=0)

    @model_validator(mode="after")
    def validate_score_ranges(self) -> "MainEnzymeSetMetrics":
        if self.min_reaction_fit_score > self.mean_reaction_fit_score:
            raise ValueError(
                "min_reaction_fit_score must not exceed "
                "mean_reaction_fit_score"
            )
        if self.min_protein_score > self.mean_protein_score:
            raise ValueError(
                "min_protein_score must not exceed mean_protein_score"
            )
        if self.min_host_fit_score > self.mean_host_fit_score:
            raise ValueError(
                "min_host_fit_score must not exceed mean_host_fit_score"
            )
        return self


class MainEnzymeSet(BaseModel):
    """One complete or review-required cover of the selected-route steps."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    set_id: int = Field(ge=1)
    set_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    status: MainEnzymeSetStatus
    protein_count: int = Field(ge=1)
    coverage_complete: bool
    covered_step_indexes: list[int]
    uncovered_step_indexes: list[int] = Field(default_factory=list)
    proteins: list[MainEnzymeSetProtein]
    step_assignments: list[MainEnzymeStepAssignment]
    review_required_step_indexes: list[int] = Field(default_factory=list)
    electron_assessment: str = Field(min_length=1)
    auxiliary_requirement_status: AuxiliaryRequirementStatus
    auxiliary_requirements: list[AuxiliaryRoleRequirement] = Field(
        default_factory=list
    )
    metrics: MainEnzymeSetMetrics
    reasons: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_consistency(self) -> "MainEnzymeSet":
        _require_sorted_unique_positive(
            self.covered_step_indexes,
            "covered_step_indexes",
            allow_empty=False,
        )
        _require_sorted_unique_positive(
            self.uncovered_step_indexes,
            "uncovered_step_indexes",
        )
        _require_sorted_unique_positive(
            self.review_required_step_indexes,
            "review_required_step_indexes",
        )
        covered = set(self.covered_step_indexes)
        uncovered = set(self.uncovered_step_indexes)
        if covered & uncovered:
            raise ValueError(
                "covered_step_indexes and uncovered_step_indexes "
                "must be disjoint"
            )
        if self.coverage_complete != (not self.uncovered_step_indexes):
            raise ValueError(
                "coverage_complete must agree with uncovered_step_indexes"
            )
        if not self.coverage_complete:
            raise ValueError(
                "main-enzyme sets must cover every required step"
            )
        if not set(self.review_required_step_indexes).issubset(covered):
            raise ValueError(
                "review_required_step_indexes must be covered steps"
            )

        protein_accessions = [protein.accession for protein in self.proteins]
        if not protein_accessions:
            raise ValueError("proteins must not be empty")
        if protein_accessions != sorted(set(protein_accessions)):
            raise ValueError("proteins must be sorted by unique accession")
        if self.protein_count != len(self.proteins):
            raise ValueError("protein_count must equal len(proteins)")
        if self.metrics.protein_count != self.protein_count:
            raise ValueError(
                "metrics.protein_count must equal protein_count"
            )

        assignment_steps = [
            assignment.step_index for assignment in self.step_assignments
        ]
        if assignment_steps != sorted(set(assignment_steps)):
            raise ValueError(
                "step_assignments must be sorted with one assignment per step"
            )
        if assignment_steps != self.covered_step_indexes:
            raise ValueError(
                "step_assignments must cover exactly covered_step_indexes"
            )
        known_accessions = set(protein_accessions)
        if any(
            assignment.accession not in known_accessions
            for assignment in self.step_assignments
        ):
            raise ValueError(
                "every step assignment must reference a set protein"
            )
        for requirement in self.auxiliary_requirements:
            if not set(requirement.main_enzyme_accessions).issubset(
                known_accessions
            ):
                raise ValueError(
                    "set auxiliary requirement references an unknown protein"
                )
            if not set(requirement.step_indexes).issubset(covered):
                raise ValueError(
                    "set auxiliary requirement references an uncovered step"
                )

        assignments_by_accession: dict[str, list[int]] = {
            accession: [] for accession in protein_accessions
        }
        for assignment in self.step_assignments:
            assignments_by_accession[assignment.accession].append(
                assignment.step_index
            )
        for protein in self.proteins:
            if (
                assignments_by_accession[protein.accession]
                != protein.assigned_step_indexes
            ):
                raise ValueError(
                    "protein assigned_step_indexes must agree with "
                    "step_assignments"
                )

        organisms = {
            protein.organism_name
            for protein in self.proteins
            if protein.organism_name
        }
        if self.metrics.organism_count != len(organisms):
            raise ValueError(
                "metrics.organism_count must equal the distinct "
                "non-empty organism count"
            )
        reviewed_fraction = sum(
            protein.reviewed for protein in self.proteins
        ) / len(self.proteins)
        if not isclose(
            self.metrics.reviewed_fraction,
            reviewed_fraction,
            rel_tol=1e-9,
            abs_tol=1e-9,
        ):
            raise ValueError(
                "metrics.reviewed_fraction must agree with proteins"
            )
        if self.metrics.warning_count != len(self.warnings):
            raise ValueError("metrics.warning_count must equal len(warnings)")
        if (
            self.metrics.auxiliary_requirement_status
            != self.auxiliary_requirement_status
        ):
            raise ValueError(
                "metrics auxiliary status must agree with the set"
            )
        pending_auxiliary_count = sum(
            item.selection_status == "pending_user_selection"
            for item in self.auxiliary_requirements
        )
        if self.metrics.pending_auxiliary_role_count != pending_auxiliary_count:
            raise ValueError(
                "pending_auxiliary_role_count must agree with requirements"
            )

        reaction_fit_scores = [
            assignment.reaction_fit_score
            for assignment in self.step_assignments
        ]
        protein_scores = [
            assignment.protein_score for assignment in self.step_assignments
        ]
        host_fit_scores = [
            assignment.host_fit_score
            for assignment in self.step_assignments
        ]
        score_checks = (
            (
                self.metrics.min_reaction_fit_score,
                min(reaction_fit_scores),
                "metrics.min_reaction_fit_score",
            ),
            (
                self.metrics.mean_reaction_fit_score,
                sum(reaction_fit_scores) / len(reaction_fit_scores),
                "metrics.mean_reaction_fit_score",
            ),
            (
                self.metrics.min_protein_score,
                min(protein_scores),
                "metrics.min_protein_score",
            ),
            (
                self.metrics.mean_protein_score,
                sum(protein_scores) / len(protein_scores),
                "metrics.mean_protein_score",
            ),
            (
                self.metrics.min_host_fit_score,
                min(host_fit_scores),
                "metrics.min_host_fit_score",
            ),
            (
                self.metrics.mean_host_fit_score,
                sum(host_fit_scores) / len(host_fit_scores),
                "metrics.mean_host_fit_score",
            ),
        )
        for actual, expected, field_name in score_checks:
            if not isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-6):
                raise ValueError(f"{field_name} must agree with assignments")

        supported_directions = {
            "supported",
            "verified",
            "compatible",
            "forward",
            "reversible",
        }
        direction_risk_steps = {
            assignment.step_index
            for assignment in self.step_assignments
            if assignment.direction_verdict.lower()
            not in supported_directions
        }
        specificity_risk_steps = {
            assignment.step_index
            for assignment in self.step_assignments
            if assignment.specificity_status.lower()
            not in {"exact", "supported"}
        }
        fit_risk_steps = {
            assignment.step_index
            for assignment in self.step_assignments
            if assignment.reaction_fit_status != "verified"
        }
        low_direction_confidence_steps = {
            assignment.step_index
            for assignment in self.step_assignments
            if assignment.direction_confidence.lower() not in {"high", "medium"}
        }
        exact_specificity_steps = {
            assignment.step_index
            for assignment in self.step_assignments
            if assignment.specificity_status.lower() == "exact"
        }
        if self.metrics.reaction_fit_risk_count != len(fit_risk_steps):
            raise ValueError(
                "metrics.reaction_fit_risk_count must agree with assignments"
            )
        if self.metrics.direction_risk_count != len(direction_risk_steps):
            raise ValueError(
                "metrics.direction_risk_count must agree with assignments"
            )
        if self.metrics.specificity_risk_count != len(
            specificity_risk_steps
        ):
            raise ValueError(
                "metrics.specificity_risk_count must agree with assignments"
            )
        if self.metrics.low_direction_confidence_count != len(
            low_direction_confidence_steps
        ):
            raise ValueError(
                "metrics.low_direction_confidence_count must agree with "
                "assignments"
            )
        if self.metrics.exact_specificity_count != len(
            exact_specificity_steps
        ):
            raise ValueError(
                "metrics.exact_specificity_count must agree with assignments"
            )
        assignment_review_steps = (
            direction_risk_steps
            | low_direction_confidence_steps
            | specificity_risk_steps
            | fit_risk_steps
        )
        if not assignment_review_steps.issubset(
            self.review_required_step_indexes
        ):
            raise ValueError(
                "review_required_step_indexes must include assignment risks"
            )

        requires_review = bool(
            self.review_required_step_indexes
            or self.warnings
            or self.metrics.direction_risk_count
            or self.metrics.low_direction_confidence_count
            or self.metrics.specificity_risk_count
            or self.metrics.electron_reassessment_status
            == "review_required"
            or any(
                assignment.reaction_fit_status != "verified"
                for assignment in self.step_assignments
            )
        )
        if self.status == "complete" and requires_review:
            raise ValueError(
                "complete set cannot contain unresolved review indicators"
            )
        if (
            self.status == "complete"
            and self.metrics.electron_reassessment_status
            not in {"not_required", "auxiliary_role_identified"}
        ):
            raise ValueError(
                "complete set has an unresolved electron reassessment"
            )
        if self.status == "review_required" and not requires_review:
            raise ValueError(
                "review_required set must contain a review indicator"
            )
        return self


class MainEnzymeSetsResult(BaseModel):
    """Canonical machine-readable enzyme-set enumeration result."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    schema_version: Literal["main_enzyme_sets.v3"] = (
        MAIN_ENZYME_SETS_SCHEMA_VERSION
    )
    algorithm_version: Literal[
        "evidence_constrained_assignment.v3_taxonomy_ranked"
    ] = (
        MAIN_ENZYME_SETS_ALGORITHM_VERSION
    )
    ok: bool
    status: MainEnzymeSetsStatus
    selected_solution_id: int = Field(ge=1)
    expansion_depth: int = Field(ge=0)
    solution_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    chassis_key: str = Field(min_length=1)
    parameters: MainEnzymeSetParameters
    candidate_pool_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    input_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    required_step_indexes: list[int]
    minimum_protein_count: int | None = Field(default=None, ge=1)
    search_complete: bool
    search_nodes: int = Field(ge=0)
    sets: list[MainEnzymeSet] = Field(default_factory=list)
    uncovered_step_indexes: list[int] = Field(default_factory=list)
    blocking_reasons: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    source_artifacts: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_consistency(self) -> "MainEnzymeSetsResult":
        _require_sorted_unique_positive(
            self.required_step_indexes,
            "required_step_indexes",
        )
        _require_sorted_unique_positive(
            self.uncovered_step_indexes,
            "uncovered_step_indexes",
        )
        required = set(self.required_step_indexes)
        if not set(self.uncovered_step_indexes).issubset(required):
            raise ValueError(
                "uncovered_step_indexes must be required steps"
            )

        set_ids = [enzyme_set.set_id for enzyme_set in self.sets]
        if set_ids != list(range(1, len(self.sets) + 1)):
            raise ValueError("set_id values must be contiguous from 1")
        set_fingerprints = [
            enzyme_set.set_fingerprint for enzyme_set in self.sets
        ]
        if len(set_fingerprints) != len(set(set_fingerprints)):
            raise ValueError("set_fingerprint values must be unique")

        for enzyme_set in self.sets:
            set_steps = set(enzyme_set.covered_step_indexes) | set(
                enzyme_set.uncovered_step_indexes
            )
            if set_steps != required:
                raise ValueError(
                    "each set's covered and uncovered steps must partition "
                    "required_step_indexes"
                )
            if any(
                not set(protein.capable_step_indexes).issubset(required)
                for protein in enzyme_set.proteins
            ):
                raise ValueError(
                    "protein capable steps must be required route steps"
                )

        covered_by_any = {
            step_index
            for enzyme_set in self.sets
            for step_index in enzyme_set.covered_step_indexes
        }
        expected_uncovered = sorted(required - covered_by_any)
        if self.uncovered_step_indexes != expected_uncovered:
            raise ValueError(
                "uncovered_step_indexes must be steps uncovered by all sets"
            )

        if self.sets:
            smallest_reported = min(
                enzyme_set.protein_count for enzyme_set in self.sets
            )
            if (
                self.minimum_protein_count is None
                or self.minimum_protein_count > smallest_reported
            ):
                raise ValueError(
                    "minimum_protein_count must not exceed the smallest "
                    "reported set"
                )
        elif self.minimum_protein_count is not None:
            raise ValueError(
                "minimum_protein_count must be null when no sets exist"
            )

        successful = self.status in {"complete", "review_required"}
        failed = self.status in {
            "infeasible",
            "stale_input",
            "source_unavailable",
        }
        if successful and not self.ok:
            raise ValueError("ok must be true for successful statuses")
        if successful and self.required_step_indexes and not self.sets:
            raise ValueError(
                "successful result with required steps must contain a set"
            )
        if successful and self.uncovered_step_indexes:
            raise ValueError(
                "successful result cannot contain uncovered required steps"
            )
        if failed and self.ok:
            raise ValueError("ok must be false for failure statuses")
        if self.status in {"stale_input", "source_unavailable"} and self.sets:
            raise ValueError(f"{self.status} result must not contain sets")
        if self.status == "truncated" and self.ok != bool(self.sets):
            raise ValueError(
                "truncated results are ok only when partial sets exist"
            )
        if self.status == "truncated" and self.search_complete:
            raise ValueError("truncated status requires search_complete=false")
        if self.status in {"complete", "review_required", "infeasible"}:
            if not self.search_complete:
                raise ValueError(
                    f"{self.status} status requires search_complete=true"
                )
        if (
            self.status == "complete"
            and self.sets
            and self.sets[0].status != "complete"
        ):
            raise ValueError("complete result requires a complete top set")
        if (
            self.status == "review_required"
            and self.sets
            and self.sets[0].status != "review_required"
        ):
            raise ValueError(
                "review_required result requires a review-required top set"
            )
        if self.status == "infeasible":
            if self.sets:
                raise ValueError("infeasible result must not contain sets")
            if self.uncovered_step_indexes != self.required_step_indexes:
                raise ValueError(
                    "infeasible result must leave every required step "
                    "uncovered"
                )
        return self


__all__ = [
    "AuxiliaryRoleRequirement",
    "MAIN_ENZYME_SELECTION_SCHEMA_VERSION",
    "MAIN_ENZYME_SETS_ALGORITHM_VERSION",
    "MAIN_ENZYME_SETS_SCHEMA_VERSION",
    "MainEnzymeCandidate",
    "MainEnzymeSet",
    "MainEnzymeSetMetrics",
    "MainEnzymeSetParameters",
    "MainEnzymeSetProtein",
    "MainEnzymeSetsResult",
    "MainEnzymeStepAssignment",
    "MainEnzymeSelectionParameters",
    "MainEnzymeSelectionResult",
]
