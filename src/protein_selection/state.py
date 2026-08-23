"""Serializable state types shared by workflow nodes."""

from typing import Any, Literal, NotRequired, TypedDict


ValidationStatus = Literal["valid", "invalid", "service_error"]
ResearchMode = Literal["balanced", "deep"]
RequirementStatus = Literal[
    "required",
    "enhancing",
    "not_required",
    "undetermined",
]
ResearchOutcome = Literal[
    "independent",
    "host_supported",
    "supplement_required",
    "reaction_mismatch",
    "unresolved",
]
CompletionStatus = Literal["success", "unresolved", "service_error"]
ValidationField = Literal["uniprot_id", "reaction_id"]
ReactionScope = Literal["single_step", "multi_step"]
WholeReactionEvidenceStatus = Literal[
    "supported",
    "uncertain",
    "unavailable",
]
ValidationErrorCode = Literal[
    "missing_input",
    "invalid_format",
    "not_found",
    "inactive_record",
    "service_error",
    "invalid_response",
]


class AssignedReactionStep(TypedDict):
    """One selected-route step assigned to a main enzyme."""

    step_index: int
    reaction_id: str
    reaction_name: str
    equation: str
    direction: str
    precursor_compound_ids: list[str]
    produced_compound_id: str
    produced_compound_name: str
    rhea_ids: list[str]


class WholeReactionContext(TypedDict):
    """Whole-reaction evidence for a single- or multi-step main enzyme."""

    equation: str | None
    rhea_ids: list[str]
    start_compound_ids: list[str]
    end_compound_ids: list[str]
    intermediate_compound_ids: list[str]
    evidence_status: WholeReactionEvidenceStatus
    evidence: list[str]


class MainEnzymeResearchUnit(TypedDict):
    """Serializable research input for one unique main-enzyme accession."""

    accession: str
    reaction_scope: ReactionScope
    assigned_step_indexes: list[int]
    reaction_steps: list[AssignedReactionStep]
    protein_name: NotRequired[str]
    organism_name: NotRequired[str]
    whole_reaction: NotRequired[WholeReactionContext | None]


class ValidationIssue(TypedDict):
    """One actionable input-validation problem."""

    code: ValidationErrorCode
    field: ValidationField
    message: str


class UniProtRecord(TypedDict):
    """Small UniProt record retained after existence validation."""

    primary_accession: str
    entry_name: str
    entry_type: str
    organism_name: str
    taxon_id: int
    source_url: str


class ReactionRecord(TypedDict):
    """Small KEGG reaction record retained after existence validation."""

    reaction_id: str
    name: str | None
    source_url: str
    names: NotRequired[list[str]]
    definition: NotRequired[str | None]
    equation: NotRequired[str | None]
    enzyme_ids: NotRequired[list[str]]
    orthology: NotRequired[list[dict[str, str]]]
    rhea_ids: NotRequired[list[str]]


class ValidationState(TypedDict, total=False):
    """LangGraph-compatible state slice used by the validation node."""

    uniprot_id: str
    reaction_id: str
    research_mode: ResearchMode
    validation_status: ValidationStatus
    validation_errors: list[ValidationIssue]
    uniprot_record: UniProtRecord | None
    reaction_record: ReactionRecord | None


class AuxiliaryRequirementAssessment(TypedDict):
    """Structured result produced by the LLM requirement node."""

    status: RequirementStatus
    reason: str
    supporting_annotations: list[str]


class ProteinSupplyState(ValidationState, total=False):
    """Application state extended with auxiliary-requirement fields."""

    research_unit: MainEnzymeResearchUnit
    reaction_records: dict[str, ReactionRecord]
    uniprot_annotation: dict[str, Any]
    research_context: dict[str, Any] | None
    main_enzyme_research_context: dict[str, Any] | None
    preliminary_reaction_match: Literal["matched", "mismatched", "uncertain"]
    preliminary_reaction_match_reason: str
    requirement_status: RequirementStatus
    requirement_reason: str
    requirement_evidence: list[str]
    requirement_assessment: AuxiliaryRequirementAssessment
    dependency_status: RequirementStatus
    can_catalyze_independently: bool | None
    completion_status: CompletionStatus
    main_research_result: dict[str, Any] | None
    research_outcome: ResearchOutcome | None
    final_protein_list: list[str]
    recommended_proteins: list[str]
    host_available_auxiliary_proteins: list[str]
    auxiliary_proteins_to_introduce: list[str]
    recommended_proteins_to_introduce: list[str]
    candidate_proteins: list[str]
    research_timings_seconds: dict[str, float]
    retrieval_stats: dict[str, int]
    research_error: str | None
