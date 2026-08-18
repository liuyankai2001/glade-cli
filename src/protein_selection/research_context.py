"""Deterministic extraction of protein and reaction research context.

The research agents should receive a compact, auditable identity packet instead
of independently interpreting a complete UniProtKB JSON document.  This module
contains no model calls and deliberately treats missing annotations as
``uncertain`` rather than negative biological evidence.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import re
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from src.protein_selection.state import (
    AssignedReactionStep,
    MainEnzymeResearchUnit,
    ReactionScope,
    WholeReactionEvidenceStatus,
)


ReactionMatch = Literal["matched", "mismatched", "uncertain"]

_RELEVANT_CROSS_REFERENCE_DATABASES = frozenset(
    {
        "BioCyc",
        "ComplexPortal",
        "GO",
        "IntAct",
        "KEGG",
        "OMA",
        "Rhea",
        "STRING",
    }
)
_RHEA_ID_PATTERN = re.compile(r"(?:RHEA:)?(\d+)", re.IGNORECASE)
_PMID_PATTERN = re.compile(r"(?:(?:PMID|PUBMED):)?(\d+)", re.IGNORECASE)
_PARENTHETICAL_ALIAS_PATTERN = re.compile(r"\(([^()]+)\)")
_BINOMIAL_NAME_PATTERN = re.compile(r"\b([A-Z][a-z]+\s+[a-z][a-z-]+)\b")
_REFERENCE_PMID_LIMIT = 12


class ProteinCatalyticActivity(BaseModel):
    """One curated catalytic-activity statement from UniProtKB."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    name: str = Field(min_length=1)
    ec_number: str | None = None
    rhea_ids: list[str] = Field(default_factory=list)
    evidence_codes: list[str] = Field(default_factory=list)


class ProteinCrossReference(BaseModel):
    """One research-relevant cross-reference from UniProtKB."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    database: str = Field(min_length=1)
    record_id: str = Field(min_length=1)
    properties: dict[str, str] = Field(default_factory=dict)
    evidence_ids: list[str] = Field(default_factory=list)


class ProteinResearchIdentity(BaseModel):
    """Canonical protein identifiers and annotations used for retrieval."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    primary_accession: str = Field(min_length=1)
    entry_name: str | None = None
    protein_names: list[str] = Field(default_factory=list)
    gene_names: list[str] = Field(default_factory=list)
    locus_tags: list[str] = Field(default_factory=list)
    organism_name: str | None = None
    organism_aliases: list[str] = Field(default_factory=list)
    taxon_id: int | None = None
    ec_numbers: list[str] = Field(default_factory=list)
    catalytic_activities: list[ProteinCatalyticActivity] = Field(
        default_factory=list
    )
    subunit_annotations: list[str] = Field(default_factory=list)
    cofactor_annotations: list[str] = Field(default_factory=list)
    cross_references: list[ProteinCrossReference] = Field(default_factory=list)
    reference_pmids: list[str] = Field(default_factory=list)


class ReactionOrthology(BaseModel):
    """One KEGG orthology entry associated with a reaction."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    orthology_id: str = Field(min_length=1)
    description: str = Field(min_length=1)


class ReactionResearchIdentity(BaseModel):
    """Normalized official KEGG reaction fields used for matching."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    reaction_id: str = Field(min_length=1)
    names: list[str] = Field(default_factory=list)
    definition: str | None = None
    equation: str | None = None
    ec_numbers: list[str] = Field(default_factory=list)
    rhea_ids: list[str] = Field(default_factory=list)
    orthology: list[ReactionOrthology] = Field(default_factory=list)
    source_url: str | None = None


class ResearchContext(BaseModel):
    """Serializable identity packet shared by every research stage."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    protein: ProteinResearchIdentity
    reaction: ReactionResearchIdentity
    preliminary_reaction_match: ReactionMatch
    preliminary_reaction_match_reason: str = Field(min_length=1)


class ReactionStepResearchContext(BaseModel):
    """One route-assigned Step plus its official reaction identity."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    step_index: int = Field(ge=1)
    route_step: AssignedReactionStep
    reaction: ReactionResearchIdentity
    preliminary_reaction_match: ReactionMatch
    preliminary_reaction_match_reason: str = Field(min_length=1)


class WholeReactionResearchContext(BaseModel):
    """Structured evidence describing one main enzyme's whole reaction."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    equation: str | None = None
    rhea_ids: list[str] = Field(default_factory=list)
    start_compound_ids: list[str] = Field(default_factory=list)
    end_compound_ids: list[str] = Field(default_factory=list)
    intermediate_compound_ids: list[str] = Field(default_factory=list)
    evidence_status: WholeReactionEvidenceStatus
    evidence: list[str] = Field(default_factory=list)


class MainEnzymeResearchContext(BaseModel):
    """Auditable protein context spanning one or more selected-route Steps."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    protein: ProteinResearchIdentity
    reaction_scope: ReactionScope
    assigned_step_indexes: list[int]
    reaction_steps: list[ReactionStepResearchContext]
    whole_reaction: WholeReactionResearchContext | None = None
    preliminary_reaction_match: ReactionMatch
    preliminary_reaction_match_reason: str = Field(min_length=1)


def build_research_context(
    uniprot_annotation: Mapping[str, Any],
    reaction_record: Mapping[str, Any],
) -> ResearchContext:
    """Build a deterministic research packet from official input records."""

    protein = _extract_protein_identity(uniprot_annotation)
    reaction = _extract_reaction_identity(reaction_record)
    status, reason = assess_preliminary_reaction_match(protein, reaction)
    return ResearchContext(
        protein=protein,
        reaction=reaction,
        preliminary_reaction_match=status,
        preliminary_reaction_match_reason=reason,
    )


def build_main_enzyme_research_context(
    uniprot_annotation: Mapping[str, Any],
    research_unit: MainEnzymeResearchUnit,
    reaction_records: Mapping[str, Mapping[str, Any]],
) -> MainEnzymeResearchContext:
    """Build one context for a main enzyme assigned to one or more Steps."""

    protein = _extract_protein_identity(uniprot_annotation)
    accession = _required_string(
        research_unit.get("accession"),
        "research_unit.accession",
    ).upper()
    if protein.primary_accession.upper() != accession:
        raise ValueError(
            "UniProt annotation accession does not match the main-enzyme "
            f"research unit: {protein.primary_accession} != {accession}"
        )

    reaction_scope = research_unit.get("reaction_scope")
    if reaction_scope not in {"single_step", "multi_step"}:
        raise ValueError("research_unit.reaction_scope is invalid")
    assigned_step_indexes = _validated_step_indexes(
        research_unit.get("assigned_step_indexes"),
        "research_unit.assigned_step_indexes",
    )
    raw_steps = research_unit.get("reaction_steps")
    if not isinstance(raw_steps, list) or any(
        not isinstance(step, Mapping) for step in raw_steps
    ):
        raise ValueError("research_unit.reaction_steps must be a list")
    if reaction_scope == "single_step" and len(assigned_step_indexes) != 1:
        raise ValueError("single_step research units must contain one Step")
    if reaction_scope == "multi_step" and len(assigned_step_indexes) < 2:
        raise ValueError("multi_step research units must contain multiple Steps")

    normalized_records: dict[str, Mapping[str, Any]] = {}
    for key, record in reaction_records.items():
        normalized_key = str(key or "").strip().upper()
        if not normalized_key:
            raise ValueError("reaction_records contains an empty reaction ID")
        if normalized_key in normalized_records:
            raise ValueError(
                f"reaction_records contains duplicate ID {normalized_key}"
            )
        if not isinstance(record, Mapping):
            raise ValueError(
                f"reaction_records[{normalized_key}] must be an object"
            )
        normalized_records[normalized_key] = record

    steps_by_index: dict[int, Mapping[str, Any]] = {}
    for raw_step in raw_steps:
        step_index = _positive_integer(
            raw_step.get("step_index"),
            "research_unit.reaction_steps[].step_index",
        )
        if step_index in steps_by_index:
            raise ValueError(
                f"research_unit contains duplicate Step {step_index}"
            )
        steps_by_index[step_index] = raw_step
    if sorted(steps_by_index) != assigned_step_indexes:
        raise ValueError(
            "research_unit assigned_step_indexes and reaction_steps do not "
            "describe the same Steps"
        )

    reaction_steps: list[ReactionStepResearchContext] = []
    for step_index in assigned_step_indexes:
        raw_step = steps_by_index[step_index]
        reaction_id = _required_string(
            raw_step.get("reaction_id"),
            f"research_unit Step {step_index} reaction_id",
        ).upper()
        raw_record = normalized_records.get(reaction_id)
        if raw_record is None:
            raise ValueError(
                f"reaction_records is missing {reaction_id} for Step "
                f"{step_index}"
            )
        reaction = _extract_reaction_identity(raw_record)
        if reaction.reaction_id.upper() != reaction_id:
            raise ValueError(
                f"reaction_records[{reaction_id}] contains "
                f"{reaction.reaction_id}"
            )
        status, reason = assess_preliminary_reaction_match(
            protein,
            reaction,
        )
        route_step = AssignedReactionStep(**dict(raw_step))
        reaction_steps.append(
            ReactionStepResearchContext(
                step_index=step_index,
                route_step=route_step,
                reaction=reaction,
                preliminary_reaction_match=status,
                preliminary_reaction_match_reason=reason,
            )
        )

    raw_whole_reaction = research_unit.get("whole_reaction")
    if raw_whole_reaction is not None and not isinstance(
        raw_whole_reaction,
        Mapping,
    ):
        raise ValueError("research_unit.whole_reaction must be an object or null")
    whole_reaction = (
        WholeReactionResearchContext.model_validate(raw_whole_reaction)
        if isinstance(raw_whole_reaction, Mapping)
        else None
    )
    overall_status, overall_reason = assess_main_enzyme_reaction_match(
        protein,
        reaction_scope,
        reaction_steps,
        whole_reaction,
    )
    return MainEnzymeResearchContext(
        protein=protein,
        reaction_scope=reaction_scope,
        assigned_step_indexes=assigned_step_indexes,
        reaction_steps=reaction_steps,
        whole_reaction=whole_reaction,
        preliminary_reaction_match=overall_status,
        preliminary_reaction_match_reason=overall_reason,
    )


def assess_main_enzyme_reaction_match(
    protein: ProteinResearchIdentity,
    reaction_scope: ReactionScope,
    reaction_steps: Sequence[ReactionStepResearchContext],
    whole_reaction: WholeReactionResearchContext | None,
) -> tuple[ReactionMatch, str]:
    """Conservatively assess one protein against its complete route scope."""

    if reaction_scope == "single_step":
        if len(reaction_steps) != 1:
            raise ValueError("single_step context must contain exactly one Step")
        step = reaction_steps[0]
        return (
            step.preliminary_reaction_match,
            "Single-step main-enzyme context: "
            + step.preliminary_reaction_match_reason,
        )
    if reaction_scope != "multi_step":
        raise ValueError(f"unsupported reaction_scope: {reaction_scope}")
    if len(reaction_steps) < 2:
        raise ValueError("multi_step context must contain at least two Steps")

    if whole_reaction is not None and whole_reaction.evidence_status == "supported":
        protein_rhea_ids = {
            numeric_id
            for activity in protein.catalytic_activities
            for rhea_id in activity.rhea_ids
            if (numeric_id := _parse_rhea_numeric_id(rhea_id)) is not None
        }
        whole_rhea_ids = {
            numeric_id
            for rhea_id in whole_reaction.rhea_ids
            if (numeric_id := _parse_rhea_numeric_id(rhea_id)) is not None
        }
        common_rhea_ids = sorted(protein_rhea_ids & whole_rhea_ids)
        if common_rhea_ids:
            return (
                "matched",
                "The supported whole-reaction context and UniProt catalytic "
                "activity share exact Rhea ID(s): "
                + ", ".join(f"RHEA:{value}" for value in common_rhea_ids),
            )
        if whole_reaction.evidence:
            return (
                "matched",
                "The structured whole-reaction context is marked supported "
                "and contains explicit source evidence for the selected "
                "multi-step enzyme.",
            )

    matched_count = sum(
        step.preliminary_reaction_match == "matched"
        for step in reaction_steps
    )
    uncertain_count = sum(
        step.preliminary_reaction_match == "uncertain"
        for step in reaction_steps
    )
    mismatched_count = sum(
        step.preliminary_reaction_match == "mismatched"
        for step in reaction_steps
    )
    return (
        "uncertain",
        "The main enzyme spans multiple decomposed route Steps, but no "
        "supported whole-reaction evidence has been established. Step-level "
        f"results are matched={matched_count}, uncertain={uncertain_count}, "
        f"mismatched={mismatched_count}; decomposed Step agreement or "
        "disagreement alone cannot resolve the whole enzyme reaction.",
    )


def assess_preliminary_reaction_match(
    protein: ProteinResearchIdentity,
    reaction: ReactionResearchIdentity,
) -> tuple[ReactionMatch, str]:
    """Return a conservative positive match from exact official mappings.

    Rhea assigns a master plus three direction-specific identifiers to one
    reaction family.  KEGG links the bidirectional member, so a UniProt
    catalytic activity is matched against that exact four-ID family rather than
    by raw string equality or global modulo arithmetic.  EC agreement alone is
    intentionally insufficient because one EC can cover multiple reactions.
    """

    protein_rhea_ids = {
        numeric_id
        for activity in protein.catalytic_activities
        for rhea_id in activity.rhea_ids
        if (numeric_id := _parse_rhea_numeric_id(rhea_id)) is not None
    }
    reaction_rhea_families = {
        family
        for rhea_id in reaction.rhea_ids
        if (family := normalize_rhea_family_id(rhea_id)) is not None
    }
    common_families = sorted(
        family
        for family in reaction_rhea_families
        if any(
            _rhea_id_belongs_to_family(numeric_id, family)
            for numeric_id in protein_rhea_ids
        )
    )
    if common_families:
        return (
            "matched",
            "UniProt and KEGG map to the same normalized Rhea reaction "
            f"family ({', '.join(common_families)}).",
        )

    protein_ec = set(protein.ec_numbers)
    reaction_ec = set(reaction.ec_numbers)
    common_ec = sorted(protein_ec & reaction_ec)
    if protein_rhea_ids and reaction_rhea_families:
        if protein_ec and reaction_ec and not common_ec:
            return (
                "mismatched",
                "UniProt and KEGG have disjoint curated Rhea mappings and "
                "disjoint EC annotations; the input protein is mapped to a "
                "different biochemical reaction than the KEGG target.",
            )
        if common_ec:
            return (
                "uncertain",
                "UniProt and KEGG share EC "
                f"{', '.join(common_ec)}, but their curated Rhea mappings "
                "do not overlap, so EC agreement alone cannot resolve the "
                "exact biochemical reaction.",
            )
        return (
            "uncertain",
            "The curated Rhea mappings do not overlap, but a definitive "
            "mismatch requires non-overlapping EC annotations on both input "
            "records; targeted verification is required.",
        )

    if common_ec:
        return (
            "uncertain",
            "UniProt and KEGG share EC "
            f"{', '.join(common_ec)}, but EC agreement alone does not prove "
            "that the exact biochemical reaction matches.",
        )

    return (
        "uncertain",
        "The official input records do not provide an exact shared Rhea "
        "mapping; targeted database or literature verification is required.",
    )


def normalize_rhea_family_id(value: str) -> str | None:
    """Return the master family for a KEGG-linked Rhea identifier.

    Rhea allocates one master reaction followed by its LR, RL and BI variants,
    but family starts are not globally aligned to a fixed modulo.  KEGG reaction
    cross-references are attached to the bidirectional (fourth) member, so its
    master is exactly three identifiers earlier.
    """

    numeric_id = _parse_rhea_numeric_id(value)
    if numeric_id is None or numeric_id <= 3:
        return None
    return f"RHEA-FAMILY:{numeric_id - 3}"


def _parse_rhea_numeric_id(value: str) -> int | None:
    match = _RHEA_ID_PATTERN.fullmatch(value.strip())
    if match is None:
        return None
    numeric_id = int(match.group(1))
    return numeric_id if numeric_id > 0 else None


def _rhea_id_belongs_to_family(
    numeric_id: int,
    family: str,
) -> bool:
    prefix = "RHEA-FAMILY:"
    if not family.startswith(prefix):
        return False
    family_start = int(family[len(prefix) :])
    return family_start <= numeric_id <= family_start + 3


def annotation_contains_text(annotation: Mapping[str, Any], text: str) -> bool:
    """Verify that a claimed annotation fragment occurs in the source JSON."""

    needle = _normalize_text(text)
    if not needle:
        return False
    return any(
        needle in _normalize_text(candidate)
        for candidate in _iter_strings(annotation)
    )


def _extract_protein_identity(
    annotation: Mapping[str, Any],
) -> ProteinResearchIdentity:
    accession = _as_nonempty_string(annotation.get("primaryAccession"))
    if accession is None:
        raise ValueError("UniProt annotation omitted primaryAccession")

    protein_names: list[str] = []
    description = annotation.get("proteinDescription")
    if isinstance(description, Mapping):
        for name_field in ("recommendedName", "submissionNames"):
            value = description.get(name_field)
            values = value if isinstance(value, list) else [value]
            for item in values:
                if isinstance(item, Mapping):
                    _append_name_value(protein_names, item.get("fullName"))
        alternative_names = description.get("alternativeNames")
        if isinstance(alternative_names, list):
            for item in alternative_names:
                if isinstance(item, Mapping):
                    _append_name_value(protein_names, item.get("fullName"))

    gene_names: list[str] = []
    locus_tags: list[str] = []
    genes = annotation.get("genes")
    if isinstance(genes, list):
        for gene in genes:
            if not isinstance(gene, Mapping):
                continue
            _append_name_value(gene_names, gene.get("geneName"))
            for field in ("synonyms",):
                values = gene.get(field)
                if isinstance(values, list):
                    for value in values:
                        _append_name_value(gene_names, value)
            for field in ("orderedLocusNames", "orfNames"):
                values = gene.get(field)
                if isinstance(values, list):
                    for value in values:
                        _append_name_value(locus_tags, value)

    activities: list[ProteinCatalyticActivity] = []
    subunit_annotations: list[str] = []
    cofactor_annotations: list[str] = []
    comments = annotation.get("comments")
    if isinstance(comments, list):
        for comment in comments:
            if not isinstance(comment, Mapping):
                continue
            comment_type = comment.get("commentType")
            if comment_type == "CATALYTIC ACTIVITY":
                reaction = comment.get("reaction")
                if not isinstance(reaction, Mapping):
                    continue
                reaction_name = _as_nonempty_string(reaction.get("name"))
                if reaction_name is None:
                    continue
                rhea_ids: list[str] = []
                crossrefs = reaction.get("reactionCrossReferences")
                if isinstance(crossrefs, list):
                    for crossref in crossrefs:
                        if (
                            isinstance(crossref, Mapping)
                            and crossref.get("database") == "Rhea"
                        ):
                            _append_string(rhea_ids, crossref.get("id"))
                evidence_codes = [
                    evidence_code
                    for evidence in reaction.get("evidences", [])
                    if isinstance(evidence, Mapping)
                    and (
                        evidence_code := _as_nonempty_string(
                            evidence.get("evidenceCode")
                        )
                    )
                    is not None
                ]
                activities.append(
                    ProteinCatalyticActivity(
                        name=reaction_name,
                        ec_number=_as_nonempty_string(
                            reaction.get("ecNumber")
                        ),
                        rhea_ids=_deduplicate(rhea_ids),
                        evidence_codes=_deduplicate(evidence_codes),
                    )
                )
            elif comment_type == "SUBUNIT":
                subunit_annotations.extend(_extract_comment_texts(comment))
            elif comment_type == "COFACTOR":
                cofactor_annotations.extend(_extract_cofactor_texts(comment))

    ec_numbers: list[str] = []
    if isinstance(description, Mapping):
        recommended = description.get("recommendedName")
        if isinstance(recommended, Mapping):
            for ec_number in recommended.get("ecNumbers", []):
                if isinstance(ec_number, Mapping):
                    _append_string(ec_numbers, ec_number.get("value"))
    for activity in activities:
        _append_string(ec_numbers, activity.ec_number)

    cross_references: list[ProteinCrossReference] = []
    raw_crossrefs = annotation.get("uniProtKBCrossReferences")
    if isinstance(raw_crossrefs, list):
        for item in raw_crossrefs:
            if not isinstance(item, Mapping):
                continue
            database = _as_nonempty_string(item.get("database"))
            record_id = _as_nonempty_string(item.get("id"))
            if (
                database not in _RELEVANT_CROSS_REFERENCE_DATABASES
                or record_id is None
            ):
                continue
            properties = {
                key: value
                for raw_property in item.get("properties", [])
                if isinstance(raw_property, Mapping)
                and (key := _as_nonempty_string(raw_property.get("key")))
                is not None
                and (
                    value := _as_nonempty_string(raw_property.get("value"))
                )
                is not None
            }
            evidence_ids = [
                f"{source}:{source_id}"
                for evidence in item.get("evidences", [])
                if isinstance(evidence, Mapping)
                and (source := _as_nonempty_string(evidence.get("source")))
                is not None
                and (
                    source_id := _as_nonempty_string(evidence.get("id"))
                )
                is not None
            ]
            cross_references.append(
                ProteinCrossReference(
                    database=database,
                    record_id=record_id,
                    properties=properties,
                    evidence_ids=_deduplicate(evidence_ids),
                )
            )

    organism = annotation.get("organism")
    organism_name = None
    organism_aliases: list[str] = []
    taxon_id = None
    if isinstance(organism, Mapping):
        organism_name = _as_nonempty_string(organism.get("scientificName"))
        _append_string(organism_aliases, organism.get("commonName"))
        raw_synonyms = organism.get("synonyms")
        if isinstance(raw_synonyms, list):
            for synonym in raw_synonyms:
                if isinstance(synonym, Mapping):
                    _append_string(organism_aliases, synonym.get("value"))
                else:
                    _append_string(organism_aliases, synonym)
        if organism_name is not None:
            organism_aliases.extend(
                match.group(1).strip()
                for match in _PARENTHETICAL_ALIAS_PATTERN.finditer(
                    organism_name
                )
                if match.group(1).strip()
            )
            organism_aliases.extend(
                _extract_reference_organism_aliases(
                    annotation,
                    organism_name,
                    gene_names,
                )
            )
        raw_taxon_id = organism.get("taxonId")
        if isinstance(raw_taxon_id, int):
            taxon_id = raw_taxon_id

    return ProteinResearchIdentity(
        primary_accession=accession,
        entry_name=_as_nonempty_string(annotation.get("uniProtkbId")),
        protein_names=_deduplicate(protein_names),
        gene_names=_deduplicate(gene_names),
        locus_tags=_deduplicate(locus_tags),
        organism_name=organism_name,
        organism_aliases=_deduplicate_casefold(organism_aliases),
        taxon_id=taxon_id,
        ec_numbers=_deduplicate(ec_numbers),
        catalytic_activities=activities,
        subunit_annotations=_deduplicate(subunit_annotations),
        cofactor_annotations=_deduplicate(cofactor_annotations),
        cross_references=cross_references,
        reference_pmids=_extract_reference_pmids(annotation),
    )


def _extract_reaction_identity(
    reaction_record: Mapping[str, Any],
) -> ReactionResearchIdentity:
    reaction_id = _as_nonempty_string(reaction_record.get("reaction_id"))
    if reaction_id is None:
        raise ValueError("KEGG reaction record omitted reaction_id")
    raw_orthology = reaction_record.get("orthology", [])
    orthology = [
        ReactionOrthology.model_validate(item)
        for item in raw_orthology
        if isinstance(item, Mapping)
    ]
    names = reaction_record.get("names", [])
    if not isinstance(names, list):
        names = []
    legacy_name = _as_nonempty_string(reaction_record.get("name"))
    normalized_names = [
        value
        for value in names
        if isinstance(value, str) and value.strip()
    ]
    if legacy_name is not None:
        normalized_names.insert(0, legacy_name)
    return ReactionResearchIdentity(
        reaction_id=reaction_id,
        names=_deduplicate(normalized_names),
        definition=_as_nonempty_string(reaction_record.get("definition")),
        equation=_as_nonempty_string(reaction_record.get("equation")),
        ec_numbers=_string_list(reaction_record.get("enzyme_ids")),
        rhea_ids=_string_list(reaction_record.get("rhea_ids")),
        orthology=orthology,
        source_url=_as_nonempty_string(reaction_record.get("source_url")),
    )


def _extract_comment_texts(comment: Mapping[str, Any]) -> list[str]:
    texts: list[str] = []
    for item in comment.get("texts", []):
        if isinstance(item, Mapping):
            _append_string(texts, item.get("value"))
    note = comment.get("note")
    if isinstance(note, Mapping):
        for item in note.get("texts", []):
            if isinstance(item, Mapping):
                _append_string(texts, item.get("value"))
    return texts


def _extract_cofactor_texts(comment: Mapping[str, Any]) -> list[str]:
    texts = _extract_comment_texts(comment)
    for cofactor in comment.get("cofactors", []):
        if isinstance(cofactor, Mapping):
            _append_string(texts, cofactor.get("name"))
    return texts


def _extract_reference_pmids(annotation: Mapping[str, Any]) -> list[str]:
    """Return up to twelve unique PubMed IDs in UniProt reference order."""

    pmids: list[str] = []
    seen: set[str] = set()
    references = annotation.get("references")
    if not isinstance(references, list):
        return pmids

    for reference in references:
        if not isinstance(reference, Mapping):
            continue
        citation = reference.get("citation")
        if not isinstance(citation, Mapping):
            continue
        crossrefs = citation.get("citationCrossReferences")
        if not isinstance(crossrefs, list):
            continue
        for crossref in crossrefs:
            if not isinstance(crossref, Mapping):
                continue
            database = _as_nonempty_string(crossref.get("database"))
            raw_id = _as_nonempty_string(crossref.get("id"))
            if database is None or database.casefold() != "pubmed":
                continue
            if raw_id is None:
                continue
            match = _PMID_PATTERN.fullmatch(raw_id)
            if match is None or int(match.group(1)) <= 0:
                continue
            pmid = match.group(1)
            if pmid in seen:
                continue
            seen.add(pmid)
            pmids.append(pmid)
            if len(pmids) == _REFERENCE_PMID_LIMIT:
                return pmids
    return pmids


def _extract_reference_organism_aliases(
    annotation: Mapping[str, Any],
    organism_name: str,
    gene_names: Iterable[str],
) -> list[str]:
    """Recover source-backed historical names from gene-specific references."""

    genus = organism_name.split(maxsplit=1)[0].casefold()
    gene_roots: set[str] = set()
    for gene in gene_names:
        normalized_gene = str(gene).strip()
        match = re.match(r"[A-Za-z]{3,}", normalized_gene)
        if match is None:
            continue
        gene_roots.add(match.group(0).casefold())
        operon_match = re.fullmatch(
            r"([a-z]{3,})[A-Z](?:\d+)?",
            normalized_gene,
        )
        if operon_match is not None:
            gene_roots.add(operon_match.group(1).casefold())
    aliases: list[str] = []
    references = annotation.get("references")
    if not isinstance(references, list):
        return aliases
    for reference in references:
        if not isinstance(reference, Mapping):
            continue
        citation = reference.get("citation")
        if not isinstance(citation, Mapping):
            continue
        title = _as_nonempty_string(citation.get("title"))
        if title is None:
            continue
        normalized_title = title.casefold()
        if gene_roots and not any(
            _reference_title_mentions_gene_root(normalized_title, root)
            for root in gene_roots
        ):
            continue
        for match in _BINOMIAL_NAME_PATTERN.finditer(title):
            alias = match.group(1).strip()
            if alias.split(maxsplit=1)[0].casefold() == genus:
                aliases.append(alias)
    return _deduplicate_casefold(aliases)


def _reference_title_mentions_gene_root(title: str, root: str) -> bool:
    """Match a gene/operon token without accepting an arbitrary substring."""

    return bool(
        re.search(
            rf"(?<![a-z0-9]){re.escape(root)}[a-z0-9]{{0,4}}(?![a-z0-9])",
            title,
            re.IGNORECASE,
        )
    )


def _required_string(value: Any, field_name: str) -> str:
    normalized = _as_nonempty_string(value)
    if normalized is None:
        raise ValueError(f"{field_name} must be a non-empty string")
    return normalized


def _positive_integer(value: Any, field_name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be a positive integer")
    try:
        normalized = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{field_name} must be a positive integer"
        ) from exc
    if normalized < 1:
        raise ValueError(f"{field_name} must be a positive integer")
    return normalized


def _validated_step_indexes(value: Any, field_name: str) -> list[int]:
    if not isinstance(value, list):
        raise ValueError(f"{field_name} must be a list")
    result = [_positive_integer(item, field_name) for item in value]
    if result != sorted(set(result)):
        raise ValueError(f"{field_name} must be sorted and unique")
    if not result:
        raise ValueError(f"{field_name} must not be empty")
    return result


def _append_name_value(target: list[str], raw: Any) -> None:
    if isinstance(raw, Mapping):
        _append_string(target, raw.get("value"))


def _append_string(target: list[str], raw: Any) -> None:
    value = _as_nonempty_string(raw)
    if value is not None:
        target.append(value)


def _as_nonempty_string(raw: Any) -> str | None:
    if not isinstance(raw, str):
        return None
    value = raw.strip()
    return value or None


def _string_list(raw: Any) -> list[str]:
    if not isinstance(raw, list):
        return []
    return _deduplicate(
        value.strip()
        for value in raw
        if isinstance(value, str) and value.strip()
    )


def _deduplicate(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _deduplicate_casefold(values: Iterable[str]) -> list[str]:
    """Deduplicate strings case-insensitively while retaining first spelling."""

    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        key = value.casefold()
        if key in seen:
            continue
        seen.add(key)
        result.append(value)
    return result


def _iter_strings(value: Any) -> Iterable[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _iter_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _iter_strings(item)


def _normalize_text(value: str) -> str:
    return " ".join(value.split()).casefold()
