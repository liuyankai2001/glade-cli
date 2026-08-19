"""Validated models for literature-backed non-standard enzyme activity."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Protocol

from pydantic import BaseModel, ConfigDict, Field, model_validator

from src.main_protein_selection.uniprot_protein_candidates import ProteinCandidate


LITERATURE_ACTIVITY_SCHEMA_VERSION = "literature_activity_evidence.v1"
LITERATURE_ACTIVITY_ALGORITHM_VERSION = "bounded_literature_activity.v1"

LiteratureActivityStatus = Literal[
    "disabled",
    "not_needed",
    "complete",
    "not_found",
    "partial",
    "source_unavailable",
]
LiteratureSource = Literal["PubMed", "Europe PMC", "Crossref"]
QueryStatus = Literal["success", "not_found", "error", "timeout", "cache_hit"]
RetractionStatus = Literal[
    "not_retracted",
    "retracted",
    "expression_of_concern",
    "not_checked",
    "uncertain",
]
AccessLevel = Literal[
    "full_text_snippets",
    "abstract_only",
    "metadata_only",
]
ActivityDirection = Literal[
    "left_to_right",
    "right_to_left",
    "bidirectional",
    "unknown",
]
ActivityRelationship = Literal["supports", "contradicts", "context_only"]
ProteinIdentifierKind = Literal[
    "accession",
    "external_protein_id",
    "external_nucleotide_id",
    "gene_or_protein_name",
]
AssayType = Literal[
    "purified_enzyme",
    "biochemical_reconstitution",
    "whole_cell_overexpression",
    "engineered_whole_cell",
    "genetic_knockout",
    "genetic_complementation",
    "cell_free_extract",
    "review_statement",
    "homology_inference",
    "computational_prediction",
    "unknown",
]
EvidenceLevel = Literal["A", "B", "C", "Reject"]
EvidenceFitStatus = Literal["verified_with_risk", "audit_only", "rejected"]
EvidenceReviewStatus = Literal["pending", "rejected"]
IdentityStatus = Literal[
    "not_attempted",
    "unique",
    "ambiguous",
    "no_hit",
    "source_unavailable",
]


class ReactionCompound(BaseModel):
    """One KEGG compound expected on one side of a selected reaction."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    compound_id: str = Field(pattern=r"^C\d{5}$")
    name: str = ""
    core: bool = False


class LiteratureActivityRequirement(BaseModel):
    """Normalized route-step context supplied to literature retrieval."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    reaction_name: str = ""
    equation: str = ""
    expected_direction: ActivityDirection
    ec_numbers: list[str] = Field(default_factory=list)
    substrates: list[ReactionCompound] = Field(default_factory=list)
    products: list[ReactionCompound] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_reaction_scope(self) -> "LiteratureActivityRequirement":
        if not any(item.core for item in self.substrates):
            raise ValueError("at least one core substrate is required")
        if not any(item.core for item in self.products):
            raise ValueError("at least one core product is required")
        return self

    @property
    def core_substrate_ids(self) -> list[str]:
        return [item.compound_id for item in self.substrates if item.core]

    @property
    def core_product_ids(self) -> list[str]:
        return [item.compound_id for item in self.products if item.core]


class LiteratureSearchQuery(BaseModel):
    """One deterministic, source-specific search query."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    query_id: str = Field(min_length=1)
    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    source: Literal["PubMed", "Europe PMC"]
    query: str = Field(min_length=1)
    rationale: str = Field(min_length=1)


class LiteratureQueryAudit(BaseModel):
    """Auditable outcome of one bounded external-source call."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    query_id: str = Field(min_length=1)
    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    source: LiteratureSource
    operation: str = Field(min_length=1)
    query: str = ""
    status: QueryStatus
    result_count: int = Field(default=0, ge=0)
    cache_hit: bool = False
    response_sha256: str = ""
    source_record_ids: list[str] = Field(default_factory=list)
    error: str = ""


class RetrievedLiteraturePaper(BaseModel):
    """Internal bounded paper record retained only for extraction and cache."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    paper_id: str = Field(min_length=1)
    title: str = ""
    abstract: str = ""
    full_text_snippets: str = ""
    authors: list[str] = Field(default_factory=list)
    journal: str = ""
    year: int | None = Field(default=None, ge=1800, le=2200)
    doi: str = ""
    pmid: str = ""
    pmcid: str = ""
    source_databases: list[LiteratureSource] = Field(default_factory=list)
    source_urls: list[str] = Field(default_factory=list)
    access_level: AccessLevel = "metadata_only"
    retraction_status: RetractionStatus = "not_checked"
    raw_evidence_ids: list[str] = Field(default_factory=list)
    source_text_sha256: str = ""

    @property
    def source_text(self) -> str:
        values = [self.title, self.abstract, self.full_text_snippets]
        return "\n\n".join(value for value in values if value)


class ExtractedActivityClaim(BaseModel):
    """Untrusted structured claim emitted by the extraction model."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    protein_identifier: str = ""
    protein_identifier_kind: ProteinIdentifierKind = "gene_or_protein_name"
    gene_name: str = ""
    protein_name: str = ""
    organism_name: str = ""
    taxon_id: int | None = Field(default=None, ge=1)
    external_identifiers: list[str] = Field(default_factory=list)
    tested_substrates: list[str] = Field(default_factory=list)
    tested_products: list[str] = Field(default_factory=list)
    matched_substrate_ids: list[str] = Field(default_factory=list)
    matched_product_ids: list[str] = Field(default_factory=list)
    assay_type: AssayType = "unknown"
    direct_activity_measured: bool = False
    direction: ActivityDirection = "unknown"
    relationship: ActivityRelationship = "context_only"
    evidence_summary: str = ""
    evidence_excerpt: str = Field(default="", max_length=400)
    source_locator: str = ""
    limitations: list[str] = Field(default_factory=list)


class PaperActivityExtraction(BaseModel):
    """Bounded result for one paper and one route step."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    claims: list[ExtractedActivityClaim] = Field(default_factory=list, max_length=5)
    extraction_notes: list[str] = Field(default_factory=list)


class LiteratureActivityEvidence(BaseModel):
    """One deterministic, traceable enzyme-activity evidence record."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    evidence_id: str = Field(min_length=1)
    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    equation: str = ""
    expected_direction: ActivityDirection
    substrate_ids: list[str] = Field(default_factory=list)
    substrate_names: list[str] = Field(default_factory=list)
    product_ids: list[str] = Field(default_factory=list)
    product_names: list[str] = Field(default_factory=list)
    protein_identifier: str = ""
    protein_identifier_kind: ProteinIdentifierKind = "gene_or_protein_name"
    gene_name: str = ""
    protein_name: str = ""
    organism_name: str = ""
    taxon_id: int | None = Field(default=None, ge=1)
    resolved_accession: str = ""
    identity_status: IdentityStatus = "not_attempted"
    identity_match_basis: list[str] = Field(default_factory=list)
    identity_query_ids: list[str] = Field(default_factory=list)
    tested_substrates: list[str] = Field(default_factory=list)
    tested_products: list[str] = Field(default_factory=list)
    matched_substrate_ids: list[str] = Field(default_factory=list)
    matched_product_ids: list[str] = Field(default_factory=list)
    evidence_type: str = "nonstandard_enzyme_activity"
    evidence_level: EvidenceLevel
    assay_type: AssayType
    direct_activity_measured: bool
    direction: ActivityDirection
    relationship: ActivityRelationship
    source_databases: list[LiteratureSource] = Field(default_factory=list)
    title: str = ""
    year: int | None = Field(default=None, ge=1800, le=2200)
    doi: str = ""
    pmid: str = ""
    pmcid: str = ""
    source_urls: list[str] = Field(default_factory=list)
    source_locator: str = ""
    access_level: AccessLevel
    retraction_status: RetractionStatus
    evidence_excerpt: str = Field(default="", max_length=400)
    evidence_summary: str = ""
    fit_status: EvidenceFitStatus
    review_status: EvidenceReviewStatus
    limitations: list[str] = Field(default_factory=list)
    validation_checks: dict[str, Any] = Field(default_factory=dict)
    rejection_reasons: list[str] = Field(default_factory=list)
    raw_evidence_ids: list[str] = Field(default_factory=list)
    source_text_sha256: str = ""


class LiteratureActivityFailure(BaseModel):
    """One non-fatal retrieval, extraction, or identity-resolution failure."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    failure_id: str = Field(min_length=1)
    step_index: int = Field(ge=1)
    reaction_id: str = Field(min_length=1)
    stage: Literal["retrieval", "extraction", "identity", "validation"]
    source: str = ""
    query_id: str = ""
    message: str = Field(min_length=1)
    retryable: bool = False


class LiteratureActivitySummary(BaseModel):
    """Small aggregate used by CLI and info consumers."""

    model_config = ConfigDict(extra="forbid")

    requirement_count: int = Field(ge=0)
    query_count: int = Field(ge=0)
    paper_count: int = Field(ge=0)
    evidence_count: int = Field(ge=0)
    grade_a_count: int = Field(ge=0)
    grade_b_count: int = Field(ge=0)
    grade_c_count: int = Field(ge=0)
    rejected_count: int = Field(ge=0)
    candidate_count: int = Field(ge=0)
    eligible_step_indexes: list[int] = Field(default_factory=list)
    unresolved_step_indexes: list[int] = Field(default_factory=list)
    failure_count: int = Field(ge=0)


class LiteratureActivityArtifact(BaseModel):
    """Canonical JSON artifact for one bounded literature search."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    schema_version: Literal["literature_activity_evidence.v1"] = (
        LITERATURE_ACTIVITY_SCHEMA_VERSION
    )
    algorithm_version: Literal["bounded_literature_activity.v1"] = (
        LITERATURE_ACTIVITY_ALGORITHM_VERSION
    )
    generated_at: str = Field(min_length=1)
    status: LiteratureActivityStatus
    enabled: bool
    cache_hit: bool = False
    request_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    artifact_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    chassis_key: str = Field(min_length=1)
    component_versions: dict[str, str] = Field(default_factory=dict)
    requirements: list[LiteratureActivityRequirement] = Field(default_factory=list)
    queries: list[LiteratureQueryAudit] = Field(default_factory=list)
    evidence: list[LiteratureActivityEvidence] = Field(default_factory=list)
    failures: list[LiteratureActivityFailure] = Field(default_factory=list)
    summary: LiteratureActivitySummary

    @model_validator(mode="after")
    def validate_consistency(self) -> "LiteratureActivityArtifact":
        if self.enabled != (self.status != "disabled"):
            raise ValueError("enabled must agree with status")
        if self.status == "disabled" and (self.queries or self.evidence or self.failures):
            raise ValueError("disabled artifact must not contain online activity")
        evidence_ids = [item.evidence_id for item in self.evidence]
        if len(evidence_ids) != len(set(evidence_ids)):
            raise ValueError("evidence_id values must be unique")
        requirement_steps = {item.step_index for item in self.requirements}
        if any(item.step_index not in requirement_steps for item in self.evidence):
            raise ValueError("evidence references an unknown requirement step")
        return self


@dataclass(slots=True)
class LiteratureActivitySearchResult:
    """Runtime result returned to the main-enzyme selector."""

    status: LiteratureActivityStatus
    candidates_by_step: dict[int, list[ProteinCandidate]]
    artifact: LiteratureActivityArtifact
    json_path: Path
    csv_path: Path
    query_errors: dict[str, str] = field(default_factory=dict)


class LiteratureRetriever(Protocol):
    """Injectable, async paper retrieval boundary used by offline tests."""

    async def retrieve(
        self,
        requirement: LiteratureActivityRequirement,
        queries: list[LiteratureSearchQuery],
        *,
        max_papers: int,
        max_full_texts: int,
    ) -> "LiteratureRetrievalBatch": ...


class ActivityExtractor(Protocol):
    """Injectable, async structured extraction boundary used by offline tests."""

    async def extract(
        self,
        requirement: LiteratureActivityRequirement,
        paper: RetrievedLiteraturePaper,
    ) -> PaperActivityExtraction: ...


class LiteratureRetrievalBatch(BaseModel):
    """Internal retriever response; failures are data, never exceptions."""

    model_config = ConfigDict(extra="forbid")

    papers: list[RetrievedLiteraturePaper] = Field(default_factory=list)
    queries: list[LiteratureQueryAudit] = Field(default_factory=list)
    failures: list[LiteratureActivityFailure] = Field(default_factory=list)


__all__ = [
    "ActivityExtractor",
    "ExtractedActivityClaim",
    "LITERATURE_ACTIVITY_ALGORITHM_VERSION",
    "LITERATURE_ACTIVITY_SCHEMA_VERSION",
    "LiteratureActivityArtifact",
    "LiteratureActivityEvidence",
    "LiteratureActivityFailure",
    "LiteratureActivityRequirement",
    "LiteratureActivitySearchResult",
    "LiteratureActivityStatus",
    "LiteratureActivitySummary",
    "LiteratureQueryAudit",
    "LiteratureRetrievalBatch",
    "LiteratureRetriever",
    "LiteratureSearchQuery",
    "PaperActivityExtraction",
    "ReactionCompound",
    "RetrievedLiteraturePaper",
]
