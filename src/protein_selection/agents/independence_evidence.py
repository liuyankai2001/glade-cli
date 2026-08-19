"""Structured evidence that a main enzyme can catalyze without another protein."""

from __future__ import annotations

from typing import Literal, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from src.protein_selection.reaction_scope import (
    ReactionScopedModel,
    normalize_reaction_ids,
)


IndependenceAssayType = Literal[
    "purified_single_protein_activity",
    "defined_reconstitution_without_partner",
    "single_gene_heterologous_activity",
    "active_homomeric_unit",
]

DIRECT_INDEPENDENCE_ASSAYS = {
    "purified_single_protein_activity",
    "defined_reconstitution_without_partner",
}
SUPPORTIVE_INDEPENDENCE_ASSAYS = {
    "single_gene_heterologous_activity",
    "active_homomeric_unit",
}


class IndependentCatalysisEvidence(ReactionScopedModel):
    """One source-grounded assay relevant to protein independence."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    independence_id: str = Field(min_length=1)
    input_uniprot_id: str = Field(min_length=1)
    assay_type: IndependenceAssayType
    evidence_ids: list[str] = Field(min_length=1)
    activity_observed: bool
    input_protein_tested: bool
    defined_protein_components: bool
    different_protein_present: bool
    protein_components: list[str] = Field(default_factory=list)
    experimental_context: str = Field(min_length=1)
    source_urls: list[str] = Field(default_factory=list)
    supporting_excerpt: str | None = None
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_assay_claim(self) -> Self:
        self.reaction_ids = normalize_reaction_ids(self.reaction_ids)
        if len(self.evidence_ids) != len(set(self.evidence_ids)):
            raise ValueError("independence evidence IDs must be unique")
        if len(self.protein_components) != len(set(self.protein_components)):
            raise ValueError("protein components must be unique")
        if len(self.source_urls) != len(set(self.source_urls)):
            raise ValueError("independence source URLs must be unique")
        if self.assay_type in DIRECT_INDEPENDENCE_ASSAYS:
            if not self.defined_protein_components:
                raise ValueError(
                    "direct independence assays require defined protein components"
                )
            if self.different_protein_present:
                raise ValueError(
                    "direct independence assays cannot contain another protein"
                )
        return self


class IndependenceAssessment(BaseModel):
    """Deterministically derived independence conclusion for one main enzyme."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    confidence: Literal["none", "medium", "high"] = "none"
    conclusion: Literal[
        "unresolved",
        "likely_independent",
        "independent",
    ] = "unresolved"
    evidence: list[IndependentCatalysisEvidence] = Field(default_factory=list)
    reasons: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_consistency(self) -> Self:
        if self.confidence == "none":
            if self.conclusion != "unresolved" or self.evidence:
                raise ValueError("no-confidence assessment cannot carry a conclusion")
        elif not self.evidence:
            raise ValueError("independence confidence requires evidence")
        if self.confidence == "high" and self.conclusion != "independent":
            raise ValueError("high-confidence evidence must conclude independent")
        if self.confidence == "medium" and self.conclusion != "likely_independent":
            raise ValueError(
                "medium-confidence evidence must conclude likely_independent"
            )
        return self


__all__ = [
    "DIRECT_INDEPENDENCE_ASSAYS",
    "SUPPORTIVE_INDEPENDENCE_ASSAYS",
    "IndependenceAssessment",
    "IndependentCatalysisEvidence",
    "IndependenceAssayType",
]
