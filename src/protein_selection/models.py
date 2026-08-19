"""Public result models for main-enzyme-combination protein research."""

from __future__ import annotations

import re
from typing import Literal, Self

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
    model_validator,
)

from src.main_protein_selection.provenance import stable_json_hash
from src.protein_selection.agents.main_research_agent import MainResearchResult
from src.protein_selection.reaction_scope import (
    normalize_reaction_ids,
    require_exact_reaction_scope,
)
from src.protein_selection.state import ReactionScope, ResearchMode


AUXILIARY_PROTEIN_RESEARCH_SCHEMA_VERSION = "auxiliary_protein_research.v3"
AUXILIARY_PROTEIN_PIPELINE_VERSION = "auxiliary_protein_pipeline.v3"

AuxiliaryRequirementClassification = Literal[
    "confirmed_required",
    "possibly_required",
    "possibly_not_required",
    "confirmed_not_required",
    "unknown",
]

MainEnzymeResearchStatus = Literal[
    "complete",
    "review_required",
    "blocked",
    "service_error",
]
CombinationResearchStatus = MainEnzymeResearchStatus

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class AuxiliaryRequirementAssessment(BaseModel):
    """Deterministic five-level auxiliary-protein conclusion."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    classification: AuxiliaryRequirementClassification
    reason: str = Field(min_length=1)
    evidence_ids: list[str] = Field(default_factory=list)
    candidate_auxiliary_accessions: list[str] = Field(default_factory=list)
    evidence_conflict: bool = False

    @field_validator("evidence_ids")
    @classmethod
    def normalize_unique_evidence_ids(cls, value: list[str]) -> list[str]:
        normalized = [item.strip() for item in value if item.strip()]
        if len(normalized) != len(value):
            raise ValueError("assessment identifiers cannot be empty")
        if normalized != list(dict.fromkeys(normalized)):
            raise ValueError("assessment identifiers must be unique")
        return normalized

    @field_validator("candidate_auxiliary_accessions")
    @classmethod
    def normalize_accessions(cls, value: list[str]) -> list[str]:
        normalized = [item.strip().upper() for item in value if item.strip()]
        if len(normalized) != len(value):
            raise ValueError("candidate auxiliary accessions cannot be empty")
        if normalized != list(dict.fromkeys(normalized)):
            raise ValueError("candidate auxiliary accessions must be unique")
        return normalized

    @model_validator(mode="after")
    def validate_classification_evidence(self) -> Self:
        if self.evidence_conflict and self.classification != "unknown":
            raise ValueError("conflicting evidence must classify as unknown")
        if self.classification != "unknown" and not self.evidence_ids:
            raise ValueError("non-unknown classifications require evidence")
        if self.classification in {
            "confirmed_required",
            "possibly_required",
        } and not self.candidate_auxiliary_accessions:
            raise ValueError(
                "required classifications need candidate auxiliary accessions"
            )
        return self


def assess_auxiliary_requirement(
    result: MainResearchResult | None,
    *,
    service_error: str | None = None,
) -> AuxiliaryRequirementAssessment:
    """Derive the five-level conclusion from validated evidence only."""

    if result is None:
        return AuxiliaryRequirementAssessment(
            classification="unknown",
            reason=(
                f"研究服务未完成：{service_error}"
                if service_error
                else "研究服务未返回可用结果。"
            ),
        )

    all_evidence_ids = list(
        dict.fromkeys(item.citation_id for item in result.evidence)
    )
    if result.reaction_match != "matched":
        return AuxiliaryRequirementAssessment(
            classification="unknown",
            reason="主酶与目标反应尚未确认匹配，不能判断辅助蛋白需求。",
            evidence_ids=all_evidence_ids,
            evidence_conflict=bool(result.conflicting_evidence),
        )

    if result.conflicting_evidence:
        return AuxiliaryRequirementAssessment(
            classification="unknown",
            reason="辅助蛋白依赖证据存在明确冲突，不能给出方向性结论。",
            evidence_ids=all_evidence_ids,
            candidate_auxiliary_accessions=list(
                dict.fromkeys(
                    [
                        *(item.uniprot_id for item in result.auxiliary_proteins),
                        *(
                            item.uniprot_id
                            for item in result.candidate_auxiliary_proteins
                        ),
                    ]
                )
            ),
            evidence_conflict=True,
        )

    confirmed_required = [
        item
        for item in result.auxiliary_proteins
        if item.necessity == "required"
    ]
    if result.outcome in {"host_supported", "supplement_required"}:
        return AuxiliaryRequirementAssessment(
            classification="confirmed_required",
            reason=result.auxiliary_requirement_reason,
            evidence_ids=list(
                dict.fromkeys(
                    evidence_id
                    for item in confirmed_required
                    for evidence_id in item.evidence_citation_ids
                )
            ),
            candidate_auxiliary_accessions=[
                item.uniprot_id for item in confirmed_required
            ],
        )

    if result.outcome == "independent":
        return AuxiliaryRequirementAssessment(
            classification="confirmed_not_required",
            reason=result.auxiliary_requirement_reason,
            evidence_ids=all_evidence_ids,
        )

    citations = {item.citation_id: item for item in result.evidence}
    possible_required_candidates = []
    possible_required_evidence_ids: list[str] = []
    for candidate in result.candidate_auxiliary_proteins:
        if candidate.proposed_necessity != "required":
            continue
        supported_ids = [
            citation_id
            for citation_id in candidate.evidence_citation_ids
            if (citation := citations.get(citation_id)) is not None
            and citation.direction == "supports"
            and citation.strength
            in {"direct_experimental", "curated", "indirect"}
        ]
        if not supported_ids:
            continue
        possible_required_candidates.append(candidate.uniprot_id)
        possible_required_evidence_ids.extend(supported_ids)

    independence = result.independence_assessment
    possible_not_required_evidence_ids = list(
        dict.fromkeys(
            evidence_id
            for item in independence.evidence
            for evidence_id in item.evidence_ids
        )
    )
    has_possible_not_required = (
        independence.confidence == "medium"
        and independence.conclusion == "likely_independent"
        and bool(possible_not_required_evidence_ids)
    )
    if possible_required_candidates and has_possible_not_required:
        return AuxiliaryRequirementAssessment(
            classification="unknown",
            reason=(
                "同时存在可能依赖和可能独立的证据，当前不能确定辅助蛋白需求。"
            ),
            evidence_ids=list(
                dict.fromkeys(
                    [
                        *possible_required_evidence_ids,
                        *possible_not_required_evidence_ids,
                    ]
                )
            ),
            candidate_auxiliary_accessions=list(
                dict.fromkeys(possible_required_candidates)
            ),
            evidence_conflict=True,
        )
    if possible_required_candidates:
        return AuxiliaryRequirementAssessment(
            classification="possibly_required",
            reason=(
                "存在支持辅助蛋白依赖的方向性证据，但尚未达到确认门槛。"
            ),
            evidence_ids=list(dict.fromkeys(possible_required_evidence_ids)),
            candidate_auxiliary_accessions=list(
                dict.fromkeys(possible_required_candidates)
            ),
        )
    if has_possible_not_required:
        return AuxiliaryRequirementAssessment(
            classification="possibly_not_required",
            reason=(
                independence.reasons[0]
                if independence.reasons
                else "存在支持主酶可能独立工作的实验依据。"
            ),
            evidence_ids=possible_not_required_evidence_ids,
        )
    return AuxiliaryRequirementAssessment(
        classification="unknown",
        reason=result.auxiliary_requirement_reason,
        evidence_ids=all_evidence_ids,
    )


def status_for_auxiliary_assessment(
    result: MainResearchResult,
    assessment: AuxiliaryRequirementAssessment,
) -> MainEnzymeResearchStatus:
    """Map one five-level assessment to workflow status."""

    if result.outcome == "reaction_mismatch":
        return "blocked"
    if assessment.classification in {
        "confirmed_required",
        "possibly_not_required",
        "confirmed_not_required",
    }:
        return "complete"
    return "review_required"


def _stable_unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _canonical_unit_result(
    item: "MainEnzymeAuxiliaryResearchResult",
) -> dict:
    payload = item.model_dump(mode="json")
    research_result = payload.get("research_result")
    if isinstance(research_result, dict):
        research_result.pop("stage_timings_seconds", None)
        research_result.pop("retrieval_stats", None)
    return payload


def auxiliary_protein_result_fingerprint(
    *,
    input_fingerprint: str,
    status: CombinationResearchStatus,
    main_enzyme_results: list["MainEnzymeAuxiliaryResearchResult"],
    main_enzyme_accessions: list[str],
    required_auxiliary_protein_accessions: list[str],
    recommended_auxiliary_protein_accessions: list[str],
    host_available_auxiliary_protein_accessions: list[str],
    auxiliary_proteins_to_introduce: list[str],
    candidate_auxiliary_protein_accessions: list[str],
    complete_protein_list: list[str],
    blocking_reasons: list[str],
    warnings: list[str],
) -> str:
    """Hash canonical biological results while excluding runtime counters."""

    return stable_json_hash(
        {
            "schema_version": AUXILIARY_PROTEIN_RESEARCH_SCHEMA_VERSION,
            "pipeline_version": AUXILIARY_PROTEIN_PIPELINE_VERSION,
            "input_fingerprint": input_fingerprint,
            "status": status,
            "main_enzyme_results": [
                _canonical_unit_result(item) for item in main_enzyme_results
            ],
            "main_enzyme_accessions": main_enzyme_accessions,
            "required_auxiliary_protein_accessions": (
                required_auxiliary_protein_accessions
            ),
            "recommended_auxiliary_protein_accessions": (
                recommended_auxiliary_protein_accessions
            ),
            "host_available_auxiliary_protein_accessions": (
                host_available_auxiliary_protein_accessions
            ),
            "auxiliary_proteins_to_introduce": auxiliary_proteins_to_introduce,
            "candidate_auxiliary_protein_accessions": (
                candidate_auxiliary_protein_accessions
            ),
            "complete_protein_list": complete_protein_list,
            "blocking_reasons": blocking_reasons,
            "warnings": warnings,
        }
    )


class MainEnzymeAuxiliaryResearchResult(BaseModel):
    """One selected main enzyme and its complete auxiliary-protein decision."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    accession: str = Field(min_length=1)
    sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    reaction_scope: ReactionScope
    assigned_step_indexes: list[int] = Field(min_length=1)
    reaction_ids: list[str] = Field(min_length=1)
    status: MainEnzymeResearchStatus
    auxiliary_requirement_assessment: AuxiliaryRequirementAssessment
    research_result: MainResearchResult | None = None
    error: str | None = None

    @field_validator("accession")
    @classmethod
    def normalize_accession(cls, value: str) -> str:
        return value.upper()

    @field_validator("sequence_sha256")
    @classmethod
    def normalize_sequence_sha256(cls, value: str) -> str:
        return value.lower()

    @field_validator("assigned_step_indexes")
    @classmethod
    def validate_step_indexes(cls, value: list[int]) -> list[int]:
        if any(isinstance(item, bool) or item < 1 for item in value):
            raise ValueError("assigned_step_indexes must contain positive integers")
        if value != sorted(set(value)):
            raise ValueError("assigned_step_indexes must be sorted and unique")
        return value

    @field_validator("reaction_ids")
    @classmethod
    def validate_reaction_ids(cls, value: list[str]) -> list[str]:
        return normalize_reaction_ids(value)

    @field_serializer("research_result")
    def serialize_research_result(
        self,
        value: MainResearchResult | None,
    ) -> dict | None:
        if value is None:
            return None
        return value.model_dump(
            mode="json",
            exclude=set(MainResearchResult.model_computed_fields),
        )

    @model_validator(mode="after")
    def validate_result_scope_and_status(self) -> Self:
        expected_scope = (
            "single_step"
            if len(self.assigned_step_indexes) == 1
            else "multi_step"
        )
        if self.reaction_scope != expected_scope:
            raise ValueError(
                "reaction_scope does not match assigned_step_indexes"
            )
        if self.status == "service_error":
            if self.research_result is not None:
                raise ValueError("service_error units cannot contain a result")
            if not self.error:
                raise ValueError("service_error units require an error message")
            expected_assessment = assess_auxiliary_requirement(
                None,
                service_error=self.error,
            )
            if self.auxiliary_requirement_assessment != expected_assessment:
                raise ValueError(
                    "service_error auxiliary assessment is inconsistent"
                )
            return self

        if self.error is not None:
            raise ValueError("non-error units cannot contain an error message")
        if self.research_result is None:
            raise ValueError("non-error units require a research result")
        if self.research_result.input_uniprot_id.upper() != self.accession:
            raise ValueError("research result changed the main-enzyme accession")
        require_exact_reaction_scope(
            self.research_result.reaction_ids,
            self.reaction_ids,
            label=f"main enzyme {self.accession}",
        )
        expected_assessment = assess_auxiliary_requirement(
            self.research_result
        )
        if self.auxiliary_requirement_assessment != expected_assessment:
            raise ValueError(
                "auxiliary requirement assessment does not match evidence"
            )
        expected_status = status_for_auxiliary_assessment(
            self.research_result,
            expected_assessment,
        )
        if self.status != expected_status:
            raise ValueError(
                "unit status does not match the main research outcome"
            )
        return self


class AuxiliaryProteinCombinationResult(BaseModel):
    """Canonical research result for one user-selected main-enzyme set."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    schema_version: Literal["auxiliary_protein_research.v3"] = (
        AUXILIARY_PROTEIN_RESEARCH_SCHEMA_VERSION
    )
    pipeline_version: Literal["auxiliary_protein_pipeline.v3"] = (
        AUXILIARY_PROTEIN_PIPELINE_VERSION
    )
    generated_at: str = Field(min_length=1)
    source_manifest: str = Field(min_length=1)
    source_manifest_revision: int = Field(ge=0)
    target_compound_id: str = Field(min_length=1)
    selected_solution_id: int = Field(ge=1)
    expansion_depth: int = Field(ge=0)
    selected_set_id: int = Field(ge=1)
    selected_set_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    solution_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    chassis_key: str = Field(min_length=1)
    research_mode: ResearchMode
    input_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    result_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    status: CombinationResearchStatus
    can_advance: bool
    main_enzyme_results: list[MainEnzymeAuxiliaryResearchResult] = Field(
        min_length=1
    )
    main_enzyme_accessions: list[str] = Field(min_length=1)
    required_auxiliary_protein_accessions: list[str] = Field(
        default_factory=list
    )
    recommended_auxiliary_protein_accessions: list[str] = Field(
        default_factory=list
    )
    host_available_auxiliary_protein_accessions: list[str] = Field(
        default_factory=list
    )
    auxiliary_proteins_to_introduce: list[str] = Field(default_factory=list)
    candidate_auxiliary_protein_accessions: list[str] = Field(
        default_factory=list
    )
    complete_protein_list: list[str] = Field(min_length=1)
    blocking_reasons: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    @field_validator("selected_set_fingerprint", "solution_fingerprint")
    @classmethod
    def normalize_source_fingerprint(cls, value: str) -> str:
        normalized = value.lower()
        if not _SHA256_PATTERN.fullmatch(normalized):
            raise ValueError("source fingerprint must be a SHA-256 digest")
        return normalized

    @field_validator(
        "main_enzyme_accessions",
        "required_auxiliary_protein_accessions",
        "recommended_auxiliary_protein_accessions",
        "host_available_auxiliary_protein_accessions",
        "auxiliary_proteins_to_introduce",
        "candidate_auxiliary_protein_accessions",
        "complete_protein_list",
    )
    @classmethod
    def validate_protein_list(cls, value: list[str]) -> list[str]:
        normalized = [item.upper() for item in value]
        if any(not item for item in normalized):
            raise ValueError("protein accession lists cannot contain empty values")
        if normalized != _stable_unique(normalized):
            raise ValueError("protein accession lists must be unique")
        return normalized

    @model_validator(mode="after")
    def validate_combination_consistency(self) -> Self:
        accessions = [item.accession for item in self.main_enzyme_results]
        if len(accessions) != len(set(accessions)):
            raise ValueError("main_enzyme_results contains duplicate accessions")
        expected_order = sorted(
            self.main_enzyme_results,
            key=lambda item: (
                min(item.assigned_step_indexes),
                item.accession,
            ),
        )
        if self.main_enzyme_results != expected_order:
            raise ValueError("main_enzyme_results is not in route order")
        if self.main_enzyme_accessions != accessions:
            raise ValueError(
                "main_enzyme_accessions does not match main_enzyme_results"
            )

        service_errors = [
            item for item in self.main_enzyme_results
            if item.status == "service_error"
        ]
        blocked = [
            item for item in self.main_enzyme_results
            if item.status == "blocked"
        ]
        reviews = [
            item for item in self.main_enzyme_results
            if item.status == "review_required"
        ]
        expected_status: CombinationResearchStatus
        if service_errors:
            expected_status = "service_error"
        elif blocked:
            expected_status = "blocked"
        elif reviews:
            expected_status = "review_required"
        else:
            expected_status = "complete"
        if self.status != expected_status:
            raise ValueError("combination status does not match unit results")
        if self.can_advance != (self.status == "complete"):
            raise ValueError("can_advance must be true only for complete results")
        if self.status == "complete" and self.blocking_reasons:
            raise ValueError("complete results cannot contain blocking reasons")
        if self.status != "complete" and not self.blocking_reasons:
            raise ValueError("non-complete results require blocking reasons")
        for item in self.main_enzyme_results:
            if (
                item.auxiliary_requirement_assessment.classification
                != "possibly_not_required"
            ):
                continue
            if not any(
                item.accession in warning and "可能不需要辅助蛋白" in warning
                for warning in self.warnings
            ):
                raise ValueError(
                    "possibly_not_required units require an explicit warning"
                )

        required: list[str] = []
        recommended: list[str] = []
        host_available: list[str] = []
        to_introduce: list[str] = []
        candidates: list[str] = []
        for item in self.main_enzyme_results:
            result = item.research_result
            if result is None:
                continue
            required.extend(
                protein.uniprot_id
                for protein in result.auxiliary_proteins
                if protein.necessity == "required"
            )
            recommended.extend(
                protein.uniprot_id
                for protein in result.auxiliary_proteins
                if protein.necessity == "enhancing"
            )
            host_available.extend(result.host_available_auxiliary_proteins)
            to_introduce.extend(result.auxiliary_proteins_to_introduce)
            candidates.extend(result.candidate_proteins)

        expected_lists = {
            "required_auxiliary_protein_accessions": _stable_unique(required),
            "recommended_auxiliary_protein_accessions": _stable_unique(
                recommended
            ),
            "host_available_auxiliary_protein_accessions": _stable_unique(
                host_available
            ),
            "auxiliary_proteins_to_introduce": _stable_unique(to_introduce),
            "candidate_auxiliary_protein_accessions": _stable_unique(
                candidates
            ),
        }
        for field_name, expected in expected_lists.items():
            if getattr(self, field_name) != expected:
                raise ValueError(f"{field_name} does not match unit results")

        expected_complete = _stable_unique(
            [*accessions, *expected_lists["required_auxiliary_protein_accessions"]]
        )
        if self.complete_protein_list != expected_complete:
            raise ValueError(
                "complete_protein_list does not match main and required proteins"
            )
        expected_fingerprint = auxiliary_protein_result_fingerprint(
            input_fingerprint=self.input_fingerprint,
            status=self.status,
            main_enzyme_results=self.main_enzyme_results,
            main_enzyme_accessions=self.main_enzyme_accessions,
            required_auxiliary_protein_accessions=(
                self.required_auxiliary_protein_accessions
            ),
            recommended_auxiliary_protein_accessions=(
                self.recommended_auxiliary_protein_accessions
            ),
            host_available_auxiliary_protein_accessions=(
                self.host_available_auxiliary_protein_accessions
            ),
            auxiliary_proteins_to_introduce=(
                self.auxiliary_proteins_to_introduce
            ),
            candidate_auxiliary_protein_accessions=(
                self.candidate_auxiliary_protein_accessions
            ),
            complete_protein_list=self.complete_protein_list,
            blocking_reasons=self.blocking_reasons,
            warnings=self.warnings,
        )
        if self.result_fingerprint != expected_fingerprint:
            raise ValueError("result_fingerprint does not match result content")
        return self


__all__ = [
    "AUXILIARY_PROTEIN_PIPELINE_VERSION",
    "AUXILIARY_PROTEIN_RESEARCH_SCHEMA_VERSION",
    "AuxiliaryRequirementAssessment",
    "AuxiliaryRequirementClassification",
    "AuxiliaryProteinCombinationResult",
    "CombinationResearchStatus",
    "MainEnzymeAuxiliaryResearchResult",
    "MainEnzymeResearchStatus",
    "assess_auxiliary_requirement",
    "auxiliary_protein_result_fingerprint",
    "status_for_auxiliary_assessment",
]
