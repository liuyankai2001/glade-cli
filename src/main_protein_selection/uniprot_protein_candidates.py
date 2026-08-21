"""UniProt 蛋白候选的检索、解析、过滤与评分。"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass, field, replace
from typing import Any, Literal

import requests

from src.main_protein_selection.settings import (
    UNIPROT_HTTP_CONFIG,
    UNIPROT_PAGE_SIZE as CONFIGURED_UNIPROT_PAGE_SIZE,
    UNIPROT_SEARCH_URL as CONFIGURED_UNIPROT_SEARCH_URL,
)
# UniProt 服务地址和请求参数统一记录在 service_config.py。
UNIPROT_SEARCH_URL = CONFIGURED_UNIPROT_SEARCH_URL
UNIPROT_PAGE_SIZE = CONFIGURED_UNIPROT_PAGE_SIZE
REQUEST_TIMEOUT = UNIPROT_HTTP_CONFIG.timeout_seconds
UNIPROT_IDENTITY_RESOLVER_VERSION = "literature_name_uniprot_resolver.v1"
# The long-form alias makes the artifact version easy to discover while the
# shorter name remains convenient for callers in this module.
LITERATURE_NAME_UNIPROT_RESOLVER_VERSION = UNIPROT_IDENTITY_RESOLVER_VERSION
UNIPROT_FIELDS = ",".join([
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "taxonomic_lineage",
    "reviewed",
    "protein_existence",
    "annotation_score",
    "ec",
    "length",
    "lineage",
    "sequence",
    "cc_catalytic_activity",
    "cc_cofactor",
    "cc_function",
    "cc_ptm",
    "cc_subunit",
    "cc_subcellular_location",
    "rhea",
    "ft_transmem",
    "ft_signal",
    "ft_transit",
    "ft_mod_res",
    "ft_binding",
    "ft_act_site",
    "ft_domain",
    "ft_crosslnk",
    "xref_interpro",
    "xref_pfam",
    "lit_pubmed_id",
    "xref_embl",
    "xref_refseq",
    "keyword",
    "fragment",
    "cc_sequence_caution",
])


CHASSIS_TAXON_PRESETS = {
    "ecoli": {
        "name": "Escherichia coli",
        "species_taxon_id": 562,
        "strain_taxon_id": None,
        "domain": "Bacteria",
        "preferred_lineage_terms": [
            "Escherichia",
            "Enterobacteriaceae",
            "Enterobacterales",
            "Gammaproteobacteria",
            "Pseudomonadota",
            "Bacteria",
        ],
    },
    "ecoli_mg1655": {
        "name": "Escherichia coli K-12 MG1655",
        "species_taxon_id": 562,
        "strain_taxon_id": 511145,
        "domain": "Bacteria",
        "preferred_lineage_terms": [
            "Escherichia coli",
            "Escherichia",
            "Enterobacteriaceae",
            "Enterobacterales",
            "Gammaproteobacteria",
            "Pseudomonadota",
            "Bacteria",
        ],
    },
}


STEP_CANDIDATE_COLUMNS = [
    "solution_id",
    "step_index",
    "reaction_id",
    "reaction_name",
    "produced_compound_id",
    "produced_compound_name",
    "role",
    "enzyme_system_type",
    "required_auxiliary_roles",
    "auxiliary_requirement_status",
    "auxiliary_requirements_json",
    "ec_number",
    "accession",
    "entry_name",
    "protein_name",
    "organism_name",
    "organism_id",
    "taxonomic_lineage",
    "reviewed",
    "length",
    "score",
    "score_breakdown",
    "evaluation_rank",
    "candidate_rank",
    "selection_status",
    "reaction_fit_status",
    "reaction_fit_score",
    "reaction_fit_rule_ids",
    "reaction_fit_evidence",
    "candidate_role",
    "component_type",
    "system_anchor_accession",
    "complex_portal_ac",
    "complex_portal_name",
    "complex_match_basis",
    "complex_evidence_type",
    "component_stoichiometry",
    "retrieval_strategy",
    "retrieval_query_id",
    "matched_rhea_ids",
    "matched_ko_ids",
    "kegg_gene_ids",
    "direction_support",
    "direction_verdict",
    "direction_confidence",
    "direction_evidence_level",
    "direction_evidence_source_ids",
    "direction_evidence",
    "required_rhea_direction_ids",
    "reaction_confidence",
    "selenzyme_rank",
    "selenzyme_score",
    "selenzyme_sim_rf",
    "selenzyme_sim_2018",
    "selenzyme_reaction_similarity",
    "selenzyme_risk_status",
    "selenzyme_matched_reaction_id",
    "selenzyme_taxonomic_distance",
    "selenzyme_direction_used",
    "selenzyme_direction_preferred",
    "gene_names",
    "ec_numbers",
    "catalytic_activities",
    "catalytic_activity_records_json",
    "cofactors",
    "subunit",
    "function_comments",
    "ptm_comments",
    "feature_annotations",
    "domain_ids",
    "keywords",
    "protein_existence",
    "sequence_version",
    "sequence_sha256",
    "aliases",
    "publication_ids",
    "cross_references",
    "subcellular_locations",
    "rhea_ids",
    "reasons",
    "warnings",
    "sequence",
]

PROTEIN_CANDIDATE_COLUMNS = [
    "accession",
    "entry_name",
    "protein_name",
    "organism_name",
    "organism_id",
    "reviewed",
    "length",
    "covered_step_indexes",
    "covered_reaction_ids",
    "covered_ec_numbers",
    "roles",
    "enzyme_system_types",
    "required_auxiliary_roles",
    "auxiliary_requirement_statuses",
    "auxiliary_requirements_json",
    "best_score",
    "mean_score",
    "min_score",
    "match_count",
    "candidate_roles",
    "component_types",
    "system_anchor_accessions",
    "complex_portal_acs",
    "complex_portal_names",
    "complex_match_bases",
    "complex_evidence_types",
    "component_stoichiometries",
    "retrieval_strategies",
    "retrieval_query_ids",
    "matched_rhea_ids",
    "matched_ko_ids",
    "kegg_gene_ids",
    "direction_support",
    "direction_verdict",
    "direction_confidence",
    "direction_evidence_level",
    "direction_evidence_source_ids",
    "required_rhea_direction_ids",
    "reaction_confidence",
    "best_selenzyme_rank",
    "selenzyme_ranks",
    "best_selenzyme_score",
    "best_selenzyme_reaction_similarity",
    "selenzyme_match_types",
    "selenzyme_selection_evidence_json",
    "selenzyme_risk_statuses",
    "selenzyme_risk_evidence_json",
    "gene_names",
    "catalytic_activities",
    "catalytic_activity_records_json",
    "cofactors",
    "subunit",
    "function_comments",
    "ptm_comments",
    "feature_annotations",
    "domain_ids",
    "keywords",
    "protein_existence",
    "sequence_version",
    "sequence_sha256",
    "subcellular_locations",
    "rhea_ids",
    "warnings",
    "sequence",
]

@dataclass
class ProteinCandidate:
    accession: str
    entry_name: str
    protein_name: str
    organism_name: str
    organism_id: int | None
    reviewed: bool
    length: int | None
    ec_numbers: list[str]
    score: float
    reasons: list[str]
    sequence: str = ""
    taxonomic_lineage: list[str] = field(default_factory=list)
    gene_names: list[str] = field(default_factory=list)
    catalytic_activities: list[str] = field(default_factory=list)
    catalytic_activity_records: list[dict[str, Any]] = field(default_factory=list)
    cofactors: list[str] = field(default_factory=list)
    subunit: list[str] = field(default_factory=list)
    function_comments: list[str] = field(default_factory=list)
    ptm_comments: list[str] = field(default_factory=list)
    feature_annotations: list[str] = field(default_factory=list)
    domain_ids: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    protein_existence: str = ""
    sequence_version: int | None = None
    sequence_sha256: str = ""
    aliases: list[str] = field(default_factory=list)
    publication_ids: list[str] = field(default_factory=list)
    cross_references: list[str] = field(default_factory=list)
    subcellular_locations: list[str] = field(default_factory=list)
    rhea_ids: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    score_breakdown: dict[str, float] = field(default_factory=dict)
    candidate_role: str = "catalytic_main"
    component_type: str = ""
    system_anchor_accession: str = ""
    complex_portal_ac: str = ""
    complex_portal_name: str = ""
    complex_match_basis: str = ""
    complex_evidence_type: str = ""
    component_stoichiometry: dict[str, float | None] = field(default_factory=dict)
    retrieval_strategy: str = "ec_exact"
    retrieval_query_id: str = ""
    matched_rhea_ids: list[str] = field(default_factory=list)
    matched_ko_ids: list[str] = field(default_factory=list)
    kegg_gene_ids: list[str] = field(default_factory=list)
    direction_support: str = ""
    direction_verdict: str = ""
    direction_confidence: str = ""
    direction_evidence_level: str = ""
    direction_evidence_source_ids: list[str] = field(default_factory=list)
    direction_evidence: list[str] = field(default_factory=list)
    required_rhea_direction_ids: list[str] = field(default_factory=list)
    reaction_confidence: str = ""
    selenzyme_rank: int | None = None
    selenzyme_score: float | None = None
    selenzyme_sim_rf: float | None = None
    selenzyme_sim_2018: float | None = None
    selenzyme_reaction_similarity: float | None = None
    selenzyme_risk_status: str = ""
    selenzyme_matched_reaction_id: str = ""
    selenzyme_taxonomic_distance: float | None = None
    selenzyme_direction_used: str = ""
    selenzyme_direction_preferred: str = ""


class _DataclassMapping(Mapping[str, Any]):
    """Small compatibility layer for typed results consumed as JSON mappings."""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)  # type: ignore[arg-type]

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __iter__(self) -> Iterator[str]:
        return (item.name for item in self.__dataclass_fields__.values())  # type: ignore[attr-defined]

    def __len__(self) -> int:
        return len(self.__dataclass_fields__)  # type: ignore[attr-defined]


@dataclass(frozen=True)
class UniProtIdentityHit(_DataclassMapping):
    """Auditable projection of one sequence-valid UniProt identity hit."""

    accession: str
    entry_name: str
    protein_name: str
    gene_names: list[str]
    aliases: list[str]
    organism_name: str
    taxon_id: int | None
    lineage: list[str]
    reviewed: bool
    sequence: str
    sequence_length: int | None
    sequence_version: int | None
    sequence_sha256: str
    publication_ids: list[str]
    cross_references: list[str]
    external_protein_ids: list[str]
    external_nucleotide_ids: list[str]
    match_basis: list[str]
    matched_role_terms: list[str]
    matched_taxon_ids: list[int]
    matched_publication_ids: list[str]
    matched_external_identifiers: list[str]
    matched_organism_names: list[str]
    matched_strain_names: list[str]
    role_match: bool
    taxon_match: bool
    identity_score: float
    query_ids: list[str]
    qualifying: bool = True
    rejection_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class UniProtIdentityResolution(_DataclassMapping):
    """Conservative resolution result; ambiguous hits are never auto-selected."""

    resolver_version: str
    status: Literal["unique", "ambiguous", "no_hit", "source_unavailable"]
    identifier: str
    identifier_kind: str
    hits: list[UniProtIdentityHit]
    queries: list[str]
    query_ids: list[str]
    match_basis: list[str]
    truncated: bool
    selected_hit: UniProtIdentityHit | None = None
    selected_accession: str = ""
    rejected_hits: list[UniProtIdentityHit] = field(default_factory=list)
    query_metadata: list[dict[str, Any]] = field(default_factory=list)
    source_errors: dict[str, str] = field(default_factory=dict)


def _unique(values: list[str] | tuple[str, ...]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def _join(values: list[Any] | tuple[Any, ...]) -> str:
    return ";".join(str(value) for value in values if str(value or "").strip())


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _split_list_field(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return _unique([str(item) for item in value])
    text = str(value).strip()
    if not text:
        return []
    return _unique([
        part.strip()
        for chunk in text.split("|")
        for part in chunk.split(";")
        if part.strip()
    ])


def _value_from_name(name: Any) -> str:
    if isinstance(name, dict):
        return str(name.get("value", "")).strip()
    return str(name or "").strip()


def _entry_is_reviewed(entry: dict[str, Any]) -> bool:
    return entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)" or entry.get("reviewed") is True


def _comment_type(comment: dict[str, Any]) -> str:
    return str(comment.get("commentType", "")).strip().upper()


def _comment_texts(comment: dict[str, Any]) -> list[str]:
    return _unique([
        str(text.get("value", "")).strip()
        for text in comment.get("texts", [])
        if isinstance(text, dict)
    ])


def _feature_type(feature: dict[str, Any]) -> str:
    return str(feature.get("type", "")).strip().lower()


def _has_feature(entry: dict[str, Any], names: set[str]) -> bool:
    return any(_feature_type(feature) in names for feature in entry.get("features", []))


def _clamp(value: float, low: float = 0.0, high: float = 100.0) -> float:
    return max(low, min(high, value))


def _format_evidence_values(values: list[str], max_items: int = 5) -> str:
    shown = values[:max_items]
    suffix = f" and {len(values) - max_items} more" if len(values) > max_items else ""
    return "; ".join(shown) + suffix


def normalize_ec_number(ec_number: str) -> str:
    return str(ec_number or "").strip().removeprefix("EC:").removeprefix("ec:")


def search_uniprot_by_ec(
    ec_number: str,
    reviewed_only: bool = False,
    size: int = UNIPROT_PAGE_SIZE,
    max_results: int = 1000,
    session: requests.Session | None = None,
) -> list[dict[str, Any]]:
    ec_number = normalize_ec_number(ec_number)
    if max_results <= 0:
        return []

    query_parts = [f"ec:{ec_number}"]
    if reviewed_only:
        query_parts.append("reviewed:true")

    params = {
        "query": " AND ".join(f"({part})" for part in query_parts),
        "format": "json",
        "fields": UNIPROT_FIELDS,
        "size": min(max(1, size), UNIPROT_PAGE_SIZE, max_results),
    }

    http = session or requests.Session()
    results: list[dict[str, Any]] = []
    next_url: str | None = UNIPROT_SEARCH_URL
    next_params: dict[str, Any] | None = params

    while next_url and len(results) < max_results:
        response = http.get(next_url, params=next_params, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        data = response.json()
        remaining = max_results - len(results)
        results.extend(data.get("results", [])[:remaining])
        response_links = getattr(response, "links", {}) or {}
        next_link = response_links.get("next", {}).get("url")
        next_url = next_link if next_link and len(results) < max_results else None
        next_params = None

    return results


def search_uniprot_by_query(
    query: str,
    reviewed_only: bool = False,
    size: int = UNIPROT_PAGE_SIZE,
    max_results: int = 1000,
    session: requests.Session | None = None,
) -> list[dict[str, Any]]:
    """Run one auditable UniProtKB query using the standard protein fields."""
    normalized_query = str(query or "").strip()
    if not normalized_query or max_results <= 0:
        return []
    query_parts = [normalized_query]
    if reviewed_only:
        query_parts.append("reviewed:true")
    params = {
        "query": " AND ".join(f"({part})" for part in query_parts),
        "format": "json",
        "fields": UNIPROT_FIELDS,
        "size": min(max(1, size), UNIPROT_PAGE_SIZE, max_results),
    }
    http = session or requests.Session()
    results: list[dict[str, Any]] = []
    next_url: str | None = UNIPROT_SEARCH_URL
    next_params: dict[str, Any] | None = params
    while next_url and len(results) < max_results:
        response = http.get(next_url, params=next_params, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        data = response.json()
        remaining = max_results - len(results)
        results.extend(data.get("results", [])[:remaining])
        response_links = getattr(response, "links", {}) or {}
        next_link = response_links.get("next", {}).get("url")
        next_url = next_link if next_link and len(results) < max_results else None
        next_params = None
    return results


def extract_ec_numbers(entry: dict[str, Any]) -> list[str]:
    ecs = set()
    protein_desc = entry.get("proteinDescription", {})
    rec = protein_desc.get("recommendedName", {})
    for ec in rec.get("ecNumbers", []):
        if ec.get("value"):
            ecs.add(str(ec["value"]))
    for alt in protein_desc.get("alternativeNames", []):
        for ec in alt.get("ecNumbers", []):
            if ec.get("value"):
                ecs.add(str(ec["value"]))
    # UniProtKB/TrEMBL entries commonly store the queried EC only under
    # submissionNames. Ignoring this branch makes exact-EC search results look
    # EC-less and causes the hard filter to discard every unreviewed candidate.
    for submission in protein_desc.get("submissionNames", []):
        for ec in submission.get("ecNumbers", []):
            if ec.get("value"):
                ecs.add(str(ec["value"]))
    return sorted(ecs)


def get_protein_name(entry: dict[str, Any]) -> str:
    protein_desc = entry.get("proteinDescription", {})
    rec = protein_desc.get("recommendedName", {})
    full_name = _value_from_name(rec.get("fullName", {}))
    if full_name:
        return full_name
    for alt in protein_desc.get("alternativeNames", []):
        full_name = _value_from_name(alt.get("fullName", {}))
        if full_name:
            return full_name
    for submission in protein_desc.get("submissionNames", []):
        full_name = _value_from_name(submission.get("fullName", {}))
        if full_name:
            return full_name
    return ""


def extract_gene_names(entry: dict[str, Any]) -> list[str]:
    names = []
    for gene in entry.get("genes", []):
        if not isinstance(gene, dict):
            continue
        names.append(_value_from_name(gene.get("geneName", {})))
        for field_name in ("synonyms", "orderedLocusNames", "orfNames"):
            for value in gene.get(field_name, []):
                names.append(_value_from_name(value))
    return _unique(names)


def extract_aliases(entry: dict[str, Any]) -> list[str]:
    """Return stable protein/gene aliases useful for exact literature lookup."""

    values = [*extract_gene_names(entry)]
    protein_desc = entry.get("proteinDescription", {})
    for key in ("recommendedName",):
        name = protein_desc.get(key, {})
        if isinstance(name, dict):
            values.append(_value_from_name(name.get("fullName", {})))
            for short in name.get("shortNames", []):
                values.append(_value_from_name(short))
    for key in ("alternativeNames", "submissionNames"):
        for name in protein_desc.get(key, []):
            if not isinstance(name, dict):
                continue
            values.append(_value_from_name(name.get("fullName", {})))
            for short in name.get("shortNames", []):
                values.append(_value_from_name(short))
    return _unique(values)


def extract_cross_references(entry: dict[str, Any]) -> list[str]:
    """Project exact identifier mappings needed by the evidence-graph audit."""

    values: list[str] = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if not isinstance(xref, dict):
            continue
        database = str(xref.get("database") or "").strip()
        identifier = str(xref.get("id") or "").strip()
        if database and identifier:
            values.append(f"{database}:{identifier}")
        for prop in xref.get("properties", []):
            if not isinstance(prop, dict):
                continue
            key = str(prop.get("key") or "").strip()
            value = str(prop.get("value") or "").strip()
            if database and key and value:
                values.append(f"{database}:{key}:{value}")
    return _unique(values)


def extract_publication_ids(entry: dict[str, Any]) -> list[str]:
    """Return PubMed/DOI anchors without relying on free-text citation parsing."""

    values: list[str] = []
    for reference in entry.get("references", []):
        citation = reference.get("citation", {}) if isinstance(reference, dict) else {}
        if not isinstance(citation, dict):
            continue
        for xref in citation.get("citationCrossReferences", []):
            if not isinstance(xref, dict):
                continue
            database = str(xref.get("database") or "").strip()
            identifier = str(xref.get("id") or "").strip()
            if database and identifier:
                values.append(f"{database}:{identifier}")
    for value in extract_cross_references(entry):
        if value.lower().startswith(("pubmed:", "doi:")):
            values.append(value)
    return _unique(values)


def get_sequence_value(entry: dict[str, Any]) -> str:
    sequence = entry.get("sequence", {})
    if isinstance(sequence, dict):
        return str(sequence.get("value", "")).strip()
    return str(sequence or "").strip()


def get_sequence_length(entry: dict[str, Any]) -> int | None:
    sequence = entry.get("sequence", {})
    length = sequence.get("length") if isinstance(sequence, dict) else None
    if length is not None:
        return int(length)
    sequence_value = get_sequence_value(entry)
    return len(sequence_value) if sequence_value else None


def has_transmembrane(entry: dict[str, Any]) -> bool:
    return _has_feature(entry, {"transmembrane"})


def has_signal_peptide(entry: dict[str, Any]) -> bool:
    return _has_feature(entry, {"signal", "signal peptide"})


def has_transit_peptide(entry: dict[str, Any]) -> bool:
    return _has_feature(entry, {"transit peptide"})


def has_sequence_caution(entry: dict[str, Any]) -> bool:
    return any(comment.get("commentType") == "SEQUENCE CAUTION" for comment in entry.get("comments", []))


def is_fragment(entry: dict[str, Any]) -> bool:
    flag = entry.get("proteinDescription", {}).get("flag")
    if flag is None:
        return False
    if isinstance(flag, str):
        return flag.lower() == "fragment"
    if isinstance(flag, dict):
        return str(flag.get("type", "")).lower() == "fragment"
    return False


def get_lineage(entry: dict[str, Any]) -> list[str]:
    return entry.get("organism", {}).get("lineage", [])


def extract_catalytic_activities(entry: dict[str, Any]) -> list[str]:
    activities = []
    for comment in entry.get("comments", []):
        if _comment_type(comment) != "CATALYTIC ACTIVITY":
            continue
        reaction = comment.get("reaction", {})
        parts = []
        if reaction.get("name"):
            parts.append(str(reaction["name"]))
        if reaction.get("ecNumber"):
            parts.append(f"EC {reaction['ecNumber']}")
        xrefs = []
        for xref in reaction.get("reactionCrossReferences", []):
            database = xref.get("database")
            identifier = xref.get("id")
            if database and identifier:
                xrefs.append(f"{database}:{identifier}")
        if xrefs:
            parts.append(", ".join(xrefs))
        if parts:
            activities.append("; ".join(parts))
    return _unique(activities)


def _reaction_cross_references(value: Any) -> list[dict[str, str]]:
    if not isinstance(value, list):
        return []
    records: list[dict[str, str]] = []
    for item in value:
        if not isinstance(item, dict):
            continue
        database = str(item.get("database") or "").strip()
        identifier = str(item.get("id") or "").strip()
        if database and identifier:
            records.append({"database": database, "id": identifier})
    return records


def _evidence_records(value: Any) -> list[dict[str, str]]:
    if not isinstance(value, list):
        return []
    records: list[dict[str, str]] = []
    for item in value:
        if not isinstance(item, dict):
            continue
        record = {
            key: str(item.get(key) or "").strip()
            for key in ("evidenceCode", "source", "id")
            if str(item.get(key) or "").strip()
        }
        if record:
            records.append(record)
    return records


def extract_catalytic_activity_records(
    entry: dict[str, Any],
) -> list[dict[str, Any]]:
    """Preserve catalytic reactions and physiological directions from UniProt."""

    records: list[dict[str, Any]] = []
    for comment in entry.get("comments", []):
        if (
            not isinstance(comment, dict)
            or _comment_type(comment) != "CATALYTIC ACTIVITY"
        ):
            continue
        reaction = comment.get("reaction")
        reaction = reaction if isinstance(reaction, dict) else {}
        physiological_records: list[dict[str, Any]] = []
        physiological_values = (
            comment.get("physiologicalReactions")
            or reaction.get("physiologicalReactions")
            or []
        )
        for physiological in physiological_values:
            if not isinstance(physiological, dict):
                continue
            physiological_records.append({
                "direction": str(
                    physiological.get("direction") or ""
                ).strip(),
                "reaction_cross_references": _reaction_cross_references(
                    physiological.get("reactionCrossReferences")
                ),
                "evidences": _evidence_records(
                    physiological.get("evidences")
                ),
            })
        records.append({
            "equation": str(reaction.get("name") or "").strip(),
            "ec_number": str(reaction.get("ecNumber") or "").strip(),
            "reaction_cross_references": _reaction_cross_references(
                reaction.get("reactionCrossReferences")
            ),
            "evidences": _evidence_records(reaction.get("evidences")),
            "physiological_reactions": physiological_records,
        })
    return records


def extract_rhea_ids(entry: dict[str, Any]) -> list[str]:
    rhea_ids = []
    for comment in entry.get("comments", []):
        if _comment_type(comment) != "CATALYTIC ACTIVITY":
            continue
        reaction = comment.get("reaction", {})
        for xref in reaction.get("reactionCrossReferences", []):
            if str(xref.get("database", "")).lower() == "rhea" and xref.get("id"):
                rhea_ids.append(str(xref["id"]))
        physiological_values = (
            comment.get("physiologicalReactions")
            or reaction.get("physiologicalReactions")
            or []
        )
        for physiological in physiological_values:
            if not isinstance(physiological, dict):
                continue
            for xref in physiological.get("reactionCrossReferences", []):
                if (
                    str(xref.get("database", "")).lower() == "rhea"
                    and xref.get("id")
                ):
                    rhea_ids.append(str(xref["id"]))
    for xref in entry.get("uniProtKBCrossReferences", []):
        if str(xref.get("database", "")).lower() == "rhea" and xref.get("id"):
            rhea_ids.append(str(xref["id"]))
    return _unique(rhea_ids)


def extract_cofactors(entry: dict[str, Any]) -> list[str]:
    cofactors = []
    for comment in entry.get("comments", []):
        if _comment_type(comment) != "COFACTOR":
            continue
        for cofactor in comment.get("cofactors", []):
            name = cofactor.get("name")
            cofactors.append(_value_from_name(name) if isinstance(name, dict) else str(name or "").strip())
        cofactors.extend(_comment_texts(comment))
    return _unique(cofactors)


def extract_subunit_comments(entry: dict[str, Any]) -> list[str]:
    subunits = []
    for comment in entry.get("comments", []):
        if _comment_type(comment) == "SUBUNIT":
            subunits.extend(_comment_texts(comment))
    return _unique(subunits)


def extract_comment_texts(entry: dict[str, Any], comment_type: str) -> list[str]:
    """Return exact UniProt comment text for one allow-listed comment type."""

    expected = str(comment_type or "").strip().upper()
    return _unique([
        text
        for comment in entry.get("comments", [])
        if isinstance(comment, dict) and _comment_type(comment) == expected
        for text in _comment_texts(comment)
    ])


def extract_feature_annotations(entry: dict[str, Any]) -> list[str]:
    """Project selection-relevant UniProt features into stable short strings."""

    allowed = {
        "active site",
        "binding site",
        "cross-link",
        "domain",
        "modified residue",
    }
    annotations: list[str] = []
    for feature in entry.get("features", []):
        if not isinstance(feature, dict):
            continue
        feature_type = _feature_type(feature)
        if feature_type not in allowed:
            continue
        location = feature.get("location", {})
        start = location.get("start", {}) if isinstance(location, dict) else {}
        end = location.get("end", {}) if isinstance(location, dict) else {}
        start_value = start.get("value") if isinstance(start, dict) else None
        end_value = end.get("value") if isinstance(end, dict) else None
        span = ""
        if start_value is not None:
            span = str(start_value)
            if end_value is not None and end_value != start_value:
                span += f"-{end_value}"
        description = str(feature.get("description") or "").strip()
        evidence_ids = _unique([
            str(item.get("evidenceCode") or item.get("source") or "").strip()
            for item in feature.get("evidences", [])
            if isinstance(item, dict)
        ])
        annotations.append(" | ".join(filter(None, [
            feature_type,
            span,
            description,
            ",".join(evidence_ids),
        ])))
    return _unique(annotations)


def extract_domain_ids(entry: dict[str, Any]) -> list[str]:
    values: list[str] = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if not isinstance(xref, dict):
            continue
        database = str(xref.get("database") or "").strip()
        identifier = str(xref.get("id") or "").strip()
        if database.lower() in {"interpro", "pfam"} and identifier:
            values.append(f"{database}:{identifier}")
    return _unique(values)


def extract_keywords(entry: dict[str, Any]) -> list[str]:
    return _unique([
        str(item.get("name") or item.get("id") or "").strip()
        for item in entry.get("keywords", [])
        if isinstance(item, dict)
    ])


def get_sequence_version(entry: dict[str, Any]) -> int | None:
    sequence = entry.get("sequence", {})
    value = sequence.get("version") if isinstance(sequence, dict) else None
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def get_sequence_sha256(entry: dict[str, Any]) -> str:
    value = get_sequence_value(entry)
    return hashlib.sha256(value.encode("ascii")).hexdigest() if value else ""


def extract_subcellular_locations(entry: dict[str, Any]) -> list[str]:
    locations = []
    for comment in entry.get("comments", []):
        if _comment_type(comment) != "SUBCELLULAR LOCATION":
            continue
        for location in comment.get("subcellularLocations", []):
            for key in ("location", "topology", "orientation"):
                value = _value_from_name(location.get(key, {}))
                if value:
                    locations.append(value)
        locations.extend(_comment_texts(comment))
    return _unique(locations)


def get_protein_existence(entry: dict[str, Any]) -> str:
    return str(entry.get("proteinExistence", "")).strip()


def get_annotation_score(entry: dict[str, Any]) -> float | None:
    value = entry.get("annotationScore")
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def hard_filter_candidate(
    entry: dict[str, Any],
    target_ec: str,
    allow_transmembrane: bool = False,
) -> tuple[bool, list[str], list[str]]:
    target_ec = normalize_ec_number(target_ec)
    reasons = []
    warnings = []
    ec_numbers = extract_ec_numbers(entry)
    if target_ec not in ec_numbers:
        return False, [f"filtered: EC is not an exact match; candidate_ecs={ec_numbers or ['none']}"], warnings
    reasons.append("EC exact match")
    if not get_sequence_value(entry):
        return False, ["filtered: missing amino acid sequence"], warnings
    if is_fragment(entry):
        return False, ["filtered: sequence is a fragment"], warnings
    if has_sequence_caution(entry):
        warnings.append("sequence caution present; verify protein sequence before CDS optimization")
    if has_transmembrane(entry):
        if not allow_transmembrane:
            return False, ["filtered: transmembrane region present and allow_transmembrane=False"], warnings
        warnings.append("transmembrane region present")
    if has_signal_peptide(entry):
        warnings.append("signal peptide present")
    if has_transit_peptide(entry):
        warnings.append("transit peptide present")
    return True, reasons, warnings


def hard_filter_candidate_without_ec(
    entry: dict[str, Any],
    allow_transmembrane: bool = False,
) -> tuple[bool, list[str], list[str]]:
    """Apply sequence/expression safety filters when reaction evidence replaces EC."""
    warnings: list[str] = []
    if not get_sequence_value(entry):
        return False, ["filtered: missing amino acid sequence"], warnings
    if is_fragment(entry):
        return False, ["filtered: sequence is a fragment"], warnings
    if has_sequence_caution(entry):
        warnings.append("sequence caution present; verify protein sequence before CDS optimization")
    if has_transmembrane(entry):
        if not allow_transmembrane:
            return False, ["filtered: transmembrane region present and allow_transmembrane=False"], warnings
        warnings.append("transmembrane region present")
    if has_signal_peptide(entry):
        warnings.append("signal peptide present")
    if has_transit_peptide(entry):
        warnings.append("transit peptide present")
    return True, ["reaction-backed candidate passed sequence safety filters"], warnings


def _score_function(entry: dict[str, Any], target_ec: str) -> tuple[float, list[str], list[str]]:
    reasons = ["function: EC exact match"]
    warnings = []
    score = 70.0
    catalytic_activities = extract_catalytic_activities(entry)
    rhea_ids = extract_rhea_ids(entry)
    if catalytic_activities:
        score += 10
        reasons.append("function: catalytic activity annotation present")
        if any(target_ec in activity for activity in catalytic_activities):
            score += 10
            reasons.append("function: catalytic activity contains target EC")
    else:
        warnings.append("missing catalytic activity text annotation")
    if rhea_ids:
        score += 10
        reasons.append("function: Rhea reaction linked")
    return _clamp(score), reasons, warnings


def _score_evidence(entry: dict[str, Any]) -> tuple[float, list[str], list[str]]:
    reasons = []
    warnings = []
    reviewed = _entry_is_reviewed(entry)
    score = 40.0 if reviewed else 15.0
    if reviewed:
        reasons.append("evidence: UniProt reviewed / Swiss-Prot")
    else:
        reasons.append("evidence: UniProt unreviewed / TrEMBL")
        warnings.append("unreviewed UniProt entry")
    annotation_score = get_annotation_score(entry)
    if annotation_score is not None:
        score += min(annotation_score, 5.0) * 7
        reasons.append(f"evidence: annotation score={annotation_score:g}/5")
    existence = get_protein_existence(entry).lower()
    if "protein level" in existence:
        score += 25
        reasons.append("evidence: protein-level existence")
    elif "transcript level" in existence:
        score += 15
        reasons.append("evidence: transcript-level existence")
    elif "homology" in existence:
        score += 10
        reasons.append("evidence: homology-based existence")
    elif "predicted" in existence:
        score += 5
        reasons.append("evidence: predicted existence")
    elif "uncertain" in existence:
        warnings.append("protein existence uncertain")
    if has_sequence_caution(entry):
        score -= 30
        warnings.append("sequence caution present; verify protein sequence before CDS optimization")
    return _clamp(score), reasons, warnings


def _score_expression(entry: dict[str, Any]) -> tuple[float, list[str], list[str]]:
    reasons = []
    warnings = []
    score = 60.0
    length = get_sequence_length(entry)
    if length is None:
        warnings.append("missing protein length")
    elif 80 <= length <= 700:
        score += 20
        reasons.append("expression: length is in common soluble-enzyme range")
    elif 700 < length <= 1200:
        score += 12
        reasons.append("expression: length acceptable but long")
    elif length > 2000:
        score -= 25
        warnings.append("protein very long")
    elif length < 50:
        score -= 25
        warnings.append("protein very short")
    else:
        score -= 8
        warnings.append("protein length outside common range")
    if has_transmembrane(entry):
        score -= 35
        warnings.append("transmembrane region increases expression risk")
    else:
        score += 10
        reasons.append("expression: no transmembrane region detected")
    if has_signal_peptide(entry):
        score -= 8
    if has_transit_peptide(entry):
        score -= 8
    cofactors = extract_cofactors(entry)
    if cofactors:
        reasons.append(f"expression: cofactor annotated ({_format_evidence_values(cofactors)})")
    return _clamp(score), reasons, warnings


def _score_host(entry: dict[str, Any], chassis_profile: dict[str, Any]) -> tuple[float, list[str]]:
    reasons = []
    organism = entry.get("organism", {})
    organism_id = organism.get("taxonId")
    lineage = get_lineage(entry)
    strain_taxon_id = chassis_profile.get("strain_taxon_id")
    species_taxon_id = chassis_profile.get("species_taxon_id")
    if strain_taxon_id and organism_id == strain_taxon_id:
        return 100.0, ["host: source strain matches chassis"]
    if organism_id == species_taxon_id:
        return 90.0, ["host: source species matches chassis"]
    preferred_terms = chassis_profile.get("preferred_lineage_terms", [])
    matched_terms = [term for term in preferred_terms if term in lineage]
    if matched_terms:
        best_index = min(preferred_terms.index(term) for term in matched_terms)
        score = max(35.0, 80.0 - best_index * 8)
        reasons.append(f"host: phylogenetically close to chassis via {matched_terms[0]}")
        return score, reasons
    return 20.0, ["host: source organism is distant from chassis"]


def score_candidate(
    entry: dict[str, Any],
    target_ec: str,
    chassis_profile: dict[str, Any],
) -> tuple[float, list[str], list[str], dict[str, float]]:
    target_ec = normalize_ec_number(target_ec)
    function_score, function_reasons, function_warnings = _score_function(entry, target_ec)
    evidence_score, evidence_reasons, evidence_warnings = _score_evidence(entry)
    expression_score, expression_reasons, expression_warnings = _score_expression(entry)
    host_score, host_reasons = _score_host(entry, chassis_profile)
    total = (
        function_score * 0.45
        + evidence_score * 0.25
        + expression_score * 0.20
        + host_score * 0.10
    )
    return (
        round(total, 2),
        _unique(function_reasons + evidence_reasons + expression_reasons + host_reasons),
        _unique(function_warnings + evidence_warnings + expression_warnings),
        {
            "function": round(function_score, 2),
            "evidence": round(evidence_score, 2),
            "expression": round(expression_score, 2),
            "host": round(host_score, 2),
            "total": round(total, 2),
        },
    )


_UNIPROT_IDENTITY_KINDS = {
    "accession",
    "external_protein_id",
    "external_nucleotide_id",
    "gene_or_protein_name",
}
_IDENTITY_QUERY_PRIORITIES = {
    "accession_exact": 0,
    "external_protein_xref": 1,
    "gene_or_protein_name_exact_role": 2,
    "external_nucleotide_xref": 3,
    "external_xref": 4,
}


@dataclass(frozen=True)
class _IdentityQuerySpec:
    strategy: str
    query: str
    identifier: str
    identifier_kind: str


def _identity_text(value: Any) -> str:
    """Normalize biological prose for exact-name and contained-role matching."""

    return " ".join(
        re.sub(r"[^\w]+", " ", str(value or "").casefold(), flags=re.UNICODE).split()
    )


def _identity_id(value: Any) -> str:
    return str(value or "").strip().casefold()


def _coerce_identity_values(values: Any) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        return _unique([values])
    if isinstance(values, Mapping):
        flattened: list[str] = []
        for nested in values.values():
            flattened.extend(_coerce_identity_values(nested))
        return _unique(flattened)
    if isinstance(values, Sequence):
        return _unique([str(value) for value in values])
    return _unique([str(values)])


def _strip_external_namespace(value: str) -> str:
    text = str(value or "").strip()
    namespace, separator, identifier = text.partition(":")
    if not separator:
        return text
    if namespace.strip().casefold() in {
        "embl",
        "ena",
        "genbank",
        "ddbj",
        "refseq",
        "proteinid",
        "protein_id",
        "nucleotidesequenceid",
        "nucleotide_id",
    }:
        return identifier.strip()
    return text


def _infer_external_identifier_kind(value: str) -> str:
    identifier = _strip_external_namespace(value).upper()
    if re.fullmatch(r"(?:NP|XP|YP|WP|AP|ZP)_\d+(?:\.\d+)?", identifier):
        return "protein"
    if re.fullmatch(r"[A-Z]{3}\d{5,}(?:\.\d+)?", identifier):
        return "protein"
    if re.fullmatch(
        r"(?:NM|NR|XM|XR|NC|NG|NT|NW|NZ)_\d+(?:\.\d+)?",
        identifier,
    ):
        return "nucleotide"
    if re.fullmatch(r"[A-Z]{1,2}\d{5,}(?:\.\d+)?", identifier):
        return "nucleotide"
    return "unknown"


def _external_identifier_items(values: Any) -> list[tuple[str, str]]:
    """Return identifier/kind pairs while accepting lists or namespaced maps."""

    items: list[tuple[str, str]] = []
    if isinstance(values, Mapping):
        for namespace, nested in values.items():
            key = str(namespace or "").casefold()
            hint = (
                "nucleotide"
                if any(term in key for term in ("nucleotide", "dna", "rna"))
                else "protein"
                if any(term in key for term in ("protein", "peptide"))
                else "unknown"
            )
            for value in _coerce_identity_values(nested):
                normalized = _strip_external_namespace(value)
                items.append((
                    normalized,
                    hint if hint != "unknown" else _infer_external_identifier_kind(normalized),
                ))
    else:
        for value in _coerce_identity_values(values):
            normalized = _strip_external_namespace(value)
            items.append((normalized, _infer_external_identifier_kind(normalized)))
    seen: set[tuple[str, str]] = set()
    answer: list[tuple[str, str]] = []
    for value, kind in items:
        key = (_identity_id(value), kind)
        if value and key not in seen:
            seen.add(key)
            answer.append((value, kind))
    return answer


def _quoted_uniprot_value(value: str) -> str:
    escaped = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _identity_query_specs(
    identifier: str,
    identifier_kind: str,
    expected_role_terms: Sequence[str],
    external_identifiers: Any,
) -> list[_IdentityQuerySpec]:
    specs: list[_IdentityQuerySpec] = []

    def add(strategy: str, value: str, kind: str, query: str) -> None:
        specs.append(_IdentityQuerySpec(
            strategy=strategy,
            query=query,
            identifier=value,
            identifier_kind=kind,
        ))

    if identifier_kind == "accession":
        add(
            "accession_exact",
            identifier,
            identifier_kind,
            f"accession:{_quoted_uniprot_value(identifier)}",
        )
    elif identifier_kind == "external_protein_id":
        add(
            "external_protein_xref",
            identifier,
            identifier_kind,
            f"xref:{_quoted_uniprot_value(identifier)}",
        )
    elif identifier_kind == "gene_or_protein_name":
        quoted_name = _quoted_uniprot_value(identifier)
        query = f"((gene_exact:{quoted_name}) OR (protein_name:{quoted_name}))"
        roles = _unique([str(term) for term in expected_role_terms])
        if roles:
            role_query = " OR ".join(_quoted_uniprot_value(term) for term in roles)
            query += f" AND ({role_query})"
        add("gene_or_protein_name_exact_role", identifier, identifier_kind, query)
    else:
        add(
            "external_nucleotide_xref",
            identifier,
            identifier_kind,
            f"xref:{_quoted_uniprot_value(identifier)}",
        )

    for external_identifier, kind in _external_identifier_items(external_identifiers):
        strategy = {
            "protein": "external_protein_xref",
            "nucleotide": "external_nucleotide_xref",
        }.get(kind, "external_xref")
        add(
            strategy,
            external_identifier,
            f"external_{kind}_id" if kind != "unknown" else "external_identifier",
            f"xref:{_quoted_uniprot_value(external_identifier)}",
        )

    specs.sort(key=lambda item: _IDENTITY_QUERY_PRIORITIES[item.strategy])
    seen_queries: set[str] = set()
    deduplicated: list[_IdentityQuerySpec] = []
    for spec in specs:
        if spec.query not in seen_queries:
            seen_queries.add(spec.query)
            deduplicated.append(spec)
    return deduplicated


def _entry_external_identifiers(
    entry: dict[str, Any],
) -> tuple[list[str], list[str], list[str]]:
    protein_ids: list[str] = []
    nucleotide_ids: list[str] = []
    all_ids: list[str] = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if not isinstance(xref, dict):
            continue
        database = str(xref.get("database") or "").strip()
        database_key = database.casefold()
        identifier = str(xref.get("id") or "").strip()
        if identifier:
            all_ids.append(identifier)
            if database_key in {"embl", "ena", "genbank", "ddbj"}:
                nucleotide_ids.append(identifier)
            elif database_key in {
                "refseq",
                "pdb",
                "ensembl",
                "ensemblplants",
                "ensemblbacteria",
            }:
                protein_ids.append(identifier)
        for prop in xref.get("properties", []):
            if not isinstance(prop, dict):
                continue
            key = str(prop.get("key") or "").casefold()
            value = str(prop.get("value") or "").strip()
            if not value:
                continue
            all_ids.append(value)
            if "protein" in key or "translation" in key:
                protein_ids.append(value)
            elif any(term in key for term in ("nucleotide", "sequenceid", "genomic")):
                nucleotide_ids.append(value)
    return _unique(protein_ids), _unique(nucleotide_ids), _unique(all_ids)


def _entry_identity_names(entry: dict[str, Any]) -> list[str]:
    return _unique([
        get_protein_name(entry),
        *extract_aliases(entry),
    ])


def _entry_organism_names(entry: dict[str, Any]) -> list[str]:
    organism = entry.get("organism", {})
    if not isinstance(organism, dict):
        return []
    values = [
        str(organism.get("scientificName") or ""),
        str(organism.get("commonName") or ""),
    ]
    for key in ("synonyms", "names"):
        nested = organism.get(key, [])
        if not isinstance(nested, list):
            continue
        for value in nested:
            values.append(_value_from_name(value))
    return _unique(values)


def _normalize_publication_id(value: str) -> tuple[str, str]:
    text = str(value or "").strip()
    namespace, separator, identifier = text.partition(":")
    namespace_key = namespace.strip().casefold()
    if separator and namespace_key in {"pubmed", "pmid"}:
        return "pubmed", identifier.strip().casefold()
    if separator and namespace_key in {"doi"}:
        return "doi", identifier.strip().casefold()
    if separator and namespace_key in {"pmc", "pmcid"}:
        return "pmcid", identifier.strip().casefold()
    if text.isdigit():
        return "pubmed", text
    if text.casefold().startswith("pmc"):
        return "pmcid", text.casefold()
    if text.casefold().startswith("10."):
        return "doi", text.casefold()
    return namespace_key if separator else "other", identifier.strip().casefold() if separator else text.casefold()


def _entry_role_text(entry: dict[str, Any]) -> str:
    values = [
        *_entry_identity_names(entry),
        *extract_catalytic_activities(entry),
        *extract_subunit_comments(entry),
        *extract_comment_texts(entry, "FUNCTION"),
        *extract_cofactors(entry),
        *extract_keywords(entry),
    ]
    return _identity_text(" ".join(values))


def _identity_quality_score(entry: dict[str, Any], match_basis: Sequence[str]) -> tuple[float, list[str]]:
    evidence_score, _, evidence_warnings = _score_evidence(entry)
    expression_score, _, expression_warnings = _score_expression(entry)
    exact_evidence_count = len([
        basis
        for basis in match_basis
        if basis.endswith("_exact") or basis.endswith("_match")
    ])
    score = _clamp(evidence_score * 0.55 + expression_score * 0.35 + min(10, exact_evidence_count * 2))
    return round(score, 2), _unique(evidence_warnings + expression_warnings)


def _strong_identity_basis(hit: UniProtIdentityHit) -> list[str]:
    direct = [
        basis
        for basis in hit.match_basis
        if basis in {
            "accession_exact",
            "external_protein_id_exact",
            "source_publication_id_exact",
        }
    ]
    if direct:
        return direct
    comprehensive_categories = {
        basis
        for basis in hit.match_basis
        if basis in {
            "gene_or_protein_name_exact",
            "external_nucleotide_id_exact",
            "external_identifier_exact",
            "expected_role_match",
            "taxon_id_exact",
            "organism_name_exact",
            "strain_name_match",
        }
    }
    return ["unique_comprehensive_qualifier"] if len(comprehensive_categories) >= 2 else []


def resolve_uniprot_identity(
    identifier: str,
    identifier_kind: Literal[
        "accession",
        "external_protein_id",
        "external_nucleotide_id",
        "gene_or_protein_name",
    ],
    expected_role_terms: Sequence[str] = (),
    taxon_ids: Sequence[int] = (),
    max_results: int = 25,
    session: requests.Session | None = None,
    source_publication_ids: Sequence[str] = (),
    external_identifiers: Any = (),
    organism_names: Sequence[str] = (),
    strain_names: Sequence[str] = (),
) -> UniProtIdentityResolution:
    """Resolve a paper/database identifier to UniProt without guessing.

    The resolver deliberately separates retrieval/ranking from identity proof.
    A reviewed entry or a higher annotation score may order hits, but neither
    can turn multiple hits into a unique result.  A unique result requires one
    sequence-valid qualifying hit plus a strong exact or combined anchor.
    """

    normalized_identifier = str(identifier or "").strip()
    normalized_kind = str(identifier_kind or "").strip()
    if not normalized_identifier:
        raise ValueError("identifier must not be empty")
    if normalized_kind not in _UNIPROT_IDENTITY_KINDS:
        allowed = ", ".join(sorted(_UNIPROT_IDENTITY_KINDS))
        raise ValueError(f"identifier_kind must be one of: {allowed}")
    if isinstance(max_results, bool) or int(max_results) <= 0:
        raise ValueError("max_results must be a positive integer")
    result_limit = int(max_results)

    roles = _coerce_identity_values(expected_role_terms)
    publication_ids = _coerce_identity_values(source_publication_ids)
    requested_organisms = _coerce_identity_values(organism_names)
    requested_strains = _coerce_identity_values(strain_names)
    requested_external = _external_identifier_items(external_identifiers)
    normalized_taxa: list[int] = []
    for value in _coerce_identity_values(taxon_ids):
        try:
            taxon_id = int(value)
        except (TypeError, ValueError):
            raise ValueError(f"invalid taxon ID: {value}") from None
        if taxon_id > 0 and taxon_id not in normalized_taxa:
            normalized_taxa.append(taxon_id)

    specs = _identity_query_specs(
        normalized_identifier,
        normalized_kind,
        roles,
        external_identifiers,
    )
    query_ids = [
        f"uniprot_identity_q{index}_{spec.strategy}"
        for index, spec in enumerate(specs, start=1)
    ]
    source_errors: dict[str, str] = {}
    query_metadata: list[dict[str, Any]] = []
    entry_by_accession: dict[str, dict[str, Any]] = {}
    entry_query_ids: dict[str, list[str]] = {}
    entry_priority: dict[str, int] = {}
    successful_queries = 0
    truncated = False

    for query_id, spec in zip(query_ids, specs):
        try:
            entries = search_uniprot_by_query(
                spec.query,
                max_results=result_limit + 1,
                session=session,
            )
            successful_queries += 1
            query_was_truncated = len(entries) > result_limit
            truncated = truncated or query_was_truncated
            query_metadata.append({
                "query_id": query_id,
                "strategy": spec.strategy,
                "query": spec.query,
                "status": "ok",
                "returned_count": len(entries),
                "truncated": query_was_truncated,
            })
        except Exception as exc:  # requests and injectable client failures
            source_errors[query_id] = f"{type(exc).__name__}: {exc}"
            query_metadata.append({
                "query_id": query_id,
                "strategy": spec.strategy,
                "query": spec.query,
                "status": "source_unavailable",
                "returned_count": 0,
                "truncated": False,
            })
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            accession = str(entry.get("primaryAccession") or "").strip().upper()
            if not accession:
                continue
            entry_by_accession.setdefault(accession, entry)
            entry_query_ids.setdefault(accession, [])
            if query_id not in entry_query_ids[accession]:
                entry_query_ids[accession].append(query_id)
            entry_priority[accession] = min(
                entry_priority.get(accession, 99),
                _IDENTITY_QUERY_PRIORITIES[spec.strategy],
            )

    if successful_queries == 0:
        return UniProtIdentityResolution(
            resolver_version=UNIPROT_IDENTITY_RESOLVER_VERSION,
            status="source_unavailable",
            identifier=normalized_identifier,
            identifier_kind=normalized_kind,
            hits=[],
            queries=[spec.query for spec in specs],
            query_ids=query_ids,
            match_basis=[],
            truncated=False,
            query_metadata=query_metadata,
            source_errors=source_errors,
        )

    ordered_accessions = sorted(
        entry_by_accession,
        key=lambda accession: (entry_priority.get(accession, 99), accession),
    )
    if len(ordered_accessions) > result_limit:
        truncated = True
        ordered_accessions = ordered_accessions[:result_limit]

    wanted_publications = {
        _normalize_publication_id(value): value for value in publication_ids
    }
    wanted_external = {
        _identity_id(value): (value, kind) for value, kind in requested_external
    }
    wanted_organisms = {
        _identity_text(value): value for value in requested_organisms
    }
    wanted_strains = {
        _identity_text(value): value for value in requested_strains
    }
    qualifying_hits: list[UniProtIdentityHit] = []
    rejected_hits: list[UniProtIdentityHit] = []

    for accession in ordered_accessions:
        entry = entry_by_accession[accession]
        protein_ids, nucleotide_ids, all_external_ids = _entry_external_identifiers(entry)
        protein_keys = {_identity_id(value) for value in protein_ids}
        nucleotide_keys = {_identity_id(value) for value in nucleotide_ids}
        external_keys = {_identity_id(value) for value in all_external_ids}
        names = _entry_identity_names(entry)
        name_keys = {_identity_text(value) for value in names}
        organism_values = _entry_organism_names(entry)
        organism_keys = {_identity_text(value) for value in organism_values}
        organism_text = _identity_text(" ".join([
            *organism_values,
            *[str(value) for value in get_lineage(entry)],
        ]))
        organism = entry.get("organism", {})
        try:
            entry_taxon_id = int(organism.get("taxonId"))
        except (TypeError, ValueError, AttributeError):
            entry_taxon_id = None

        secondary_accessions = {
            str(value or "").strip().upper()
            for value in entry.get("secondaryAccessions", [])
            if str(value or "").strip()
        }
        identifier_key = _identity_id(_strip_external_namespace(normalized_identifier))
        match_basis: list[str] = []
        main_match = False
        if normalized_kind == "accession":
            main_match = normalized_identifier.upper() in {accession, *secondary_accessions}
            if main_match:
                match_basis.append("accession_exact")
        elif normalized_kind == "external_protein_id":
            main_match = identifier_key in protein_keys
            if main_match:
                match_basis.append("external_protein_id_exact")
        elif normalized_kind == "external_nucleotide_id":
            main_match = identifier_key in nucleotide_keys
            if main_match:
                match_basis.append("external_nucleotide_id_exact")
        else:
            main_match = _identity_text(normalized_identifier) in name_keys
            if main_match:
                match_basis.append("gene_or_protein_name_exact")

        role_text = _entry_role_text(entry)
        matched_roles = [
            role for role in roles
            if _identity_text(role) and _identity_text(role) in role_text
        ]
        role_match = not roles or bool(matched_roles)
        if roles and role_match:
            match_basis.append("expected_role_match")

        matched_taxa = (
            [entry_taxon_id]
            if entry_taxon_id is not None and entry_taxon_id in normalized_taxa
            else []
        )
        taxon_match = not normalized_taxa or bool(matched_taxa)
        if normalized_taxa and taxon_match:
            match_basis.append("taxon_id_exact")

        entry_publications = extract_publication_ids(entry)
        entry_publication_map = {
            _normalize_publication_id(value): value for value in entry_publications
        }
        matched_publications = [
            entry_publication_map[key]
            for key in wanted_publications
            if key in entry_publication_map
        ]
        publication_match = not wanted_publications or bool(matched_publications)
        if wanted_publications and publication_match:
            match_basis.append("source_publication_id_exact")

        matched_external: list[str] = []
        matched_external_kinds: set[str] = set()
        for key, (value, kind) in wanted_external.items():
            kind_keys = (
                protein_keys
                if kind == "protein"
                else nucleotide_keys
                if kind == "nucleotide"
                else external_keys
            )
            if key in kind_keys:
                matched_external.append(value)
                matched_external_kinds.add(kind)
        external_match = not wanted_external or bool(matched_external)
        if wanted_external and external_match:
            if "protein" in matched_external_kinds:
                match_basis.append("external_protein_id_exact")
            if "nucleotide" in matched_external_kinds:
                match_basis.append("external_nucleotide_id_exact")
            if matched_external_kinds == {"unknown"}:
                match_basis.append("external_identifier_exact")

        matched_organisms = [
            original
            for key, original in wanted_organisms.items()
            if key in organism_keys
        ]
        organism_match = not wanted_organisms or bool(matched_organisms)
        if wanted_organisms and organism_match:
            match_basis.append("organism_name_exact")

        matched_strains = [
            original
            for key, original in wanted_strains.items()
            if len(key) >= 4
            and re.search(
                rf"(?<!\w){re.escape(key)}(?!\w)",
                organism_text,
            )
        ]
        strain_match = not wanted_strains or bool(matched_strains)
        if wanted_strains and strain_match:
            match_basis.append("strain_name_match")

        sequence_passed, filter_reasons, filter_warnings = hard_filter_candidate_without_ec(
            entry,
            # Identity resolution must not reject a real integral-membrane
            # identity; the expression-design gate can apply that policy later.
            allow_transmembrane=True,
        )
        rejection_reasons: list[str] = []
        if not main_match:
            rejection_reasons.append("identifier did not match the requested identifier kind exactly")
        if not role_match:
            rejection_reasons.append("expected role terms did not match UniProt annotation")
        if not taxon_match:
            rejection_reasons.append("taxon constraint did not match")
        # Literature organism strings are disambiguation hints, not hard
        # identity constraints.  Historical taxonomic renames are common
        # (for example, the same deposited strain can move to another species
        # name), so a missing literal name/strain match must not discard an
        # otherwise exact accession/publication mapping.  When at least one
        # hit does match, the positive hint is applied below to narrow the
        # candidate set.
        if not sequence_passed:
            rejection_reasons.extend(filter_reasons)

        match_basis = _unique(match_basis)
        identity_score, score_warnings = _identity_quality_score(entry, match_basis)
        hit = UniProtIdentityHit(
            accession=accession,
            entry_name=str(entry.get("uniProtkbId") or ""),
            protein_name=get_protein_name(entry),
            gene_names=extract_gene_names(entry),
            aliases=extract_aliases(entry),
            organism_name=str(organism.get("scientificName") or "") if isinstance(organism, dict) else "",
            taxon_id=entry_taxon_id,
            lineage=[str(value) for value in get_lineage(entry)],
            reviewed=_entry_is_reviewed(entry),
            sequence=get_sequence_value(entry),
            sequence_length=get_sequence_length(entry),
            sequence_version=get_sequence_version(entry),
            sequence_sha256=get_sequence_sha256(entry),
            publication_ids=entry_publications,
            cross_references=extract_cross_references(entry),
            external_protein_ids=protein_ids,
            external_nucleotide_ids=nucleotide_ids,
            match_basis=match_basis,
            matched_role_terms=matched_roles,
            matched_taxon_ids=matched_taxa,
            matched_publication_ids=matched_publications,
            matched_external_identifiers=matched_external,
            matched_organism_names=matched_organisms,
            matched_strain_names=matched_strains,
            role_match=role_match,
            taxon_match=taxon_match,
            identity_score=identity_score,
            query_ids=entry_query_ids.get(accession, []),
            qualifying=not rejection_reasons,
            rejection_reasons=_unique(rejection_reasons),
            warnings=_unique(filter_warnings + score_warnings),
        )
        if hit.qualifying:
            qualifying_hits.append(hit)
        else:
            rejected_hits.append(hit)

    # Publication and secondary external identifiers are positive
    # disambiguators, not negative evidence.  If none of the otherwise valid
    # candidates shares one, retain the ambiguous name-derived candidates.
    # If one or more do share it, the exact shared anchor may safely narrow the
    # pool.  Conflicting positive anchors never authorize a selection.
    base_hits = list(qualifying_hits)
    active_disambiguators: list[tuple[str, Any]] = []
    if wanted_publications and any(hit.matched_publication_ids for hit in base_hits):
        active_disambiguators.append((
            "source publication identifier was not shared",
            lambda hit: bool(hit.matched_publication_ids),
        ))
    if wanted_external and any(hit.matched_external_identifiers for hit in base_hits):
        active_disambiguators.append((
            "external identifier was not shared",
            lambda hit: bool(hit.matched_external_identifiers),
        ))
    if wanted_organisms and any(hit.matched_organism_names for hit in base_hits):
        active_disambiguators.append((
            "organism name hint was not shared",
            lambda hit: bool(hit.matched_organism_names),
        ))
    if wanted_strains and any(hit.matched_strain_names for hit in base_hits):
        active_disambiguators.append((
            "strain name hint was not shared",
            lambda hit: bool(hit.matched_strain_names),
        ))
    if active_disambiguators:
        narrowed_hits = [
            hit
            for hit in base_hits
            if all(predicate(hit) for _, predicate in active_disambiguators)
        ]
        if narrowed_hits:
            qualifying_hits = narrowed_hits
            narrowed_accessions = {hit.accession for hit in narrowed_hits}
            for hit in base_hits:
                if hit.accession in narrowed_accessions:
                    continue
                reasons = [
                    reason
                    for reason, predicate in active_disambiguators
                    if not predicate(hit)
                ]
                rejected_hits.append(replace(
                    hit,
                    qualifying=False,
                    rejection_reasons=_unique(hit.rejection_reasons + reasons),
                ))

    qualifying_hits.sort(key=lambda hit: (
        min((
            int(query_id.split("_q", 1)[1].split("_", 1)[0])
            for query_id in hit.query_ids
        ), default=999),
        -hit.identity_score,
        hit.accession,
    ))
    rejected_hits.sort(key=lambda hit: hit.accession)

    selected_hit: UniProtIdentityHit | None = None
    selected_accession = ""
    strong_basis: list[str] = []
    if not qualifying_hits:
        status: Literal["unique", "ambiguous", "no_hit", "source_unavailable"] = (
            "ambiguous" if truncated else "no_hit"
        )
    elif len(qualifying_hits) == 1 and not truncated:
        strong_basis = _strong_identity_basis(qualifying_hits[0])
        if strong_basis:
            status = "unique"
            selected_hit = qualifying_hits[0]
            selected_accession = selected_hit.accession
        else:
            status = "ambiguous"
    else:
        status = "ambiguous"

    result_match_basis = _unique([
        basis for hit in qualifying_hits for basis in hit.match_basis
    ] + strong_basis)
    return UniProtIdentityResolution(
        resolver_version=UNIPROT_IDENTITY_RESOLVER_VERSION,
        status=status,
        identifier=normalized_identifier,
        identifier_kind=normalized_kind,
        hits=qualifying_hits,
        queries=[spec.query for spec in specs],
        query_ids=query_ids,
        match_basis=result_match_basis,
        truncated=truncated,
        selected_hit=selected_hit,
        selected_accession=selected_accession,
        rejected_hits=rejected_hits,
        query_metadata=query_metadata,
        source_errors=source_errors,
    )


def _merge_entries(entry_groups: list[list[dict[str, Any]]]) -> list[dict[str, Any]]:
    merged = {}
    for entries in entry_groups:
        for entry in entries:
            accession = entry.get("primaryAccession")
            if accession and accession not in merged:
                merged[accession] = entry
    return list(merged.values())


def _candidate_from_entry(
    entry: dict[str, Any],
    score: float,
    reasons: list[str],
    warnings: list[str],
    score_breakdown: dict[str, float],
) -> ProteinCandidate:
    organism = entry.get("organism", {})
    return ProteinCandidate(
        accession=entry.get("primaryAccession", ""),
        entry_name=entry.get("uniProtkbId", ""),
        protein_name=get_protein_name(entry),
        organism_name=organism.get("scientificName", ""),
        organism_id=organism.get("taxonId"),
        taxonomic_lineage=[str(value) for value in get_lineage(entry)],
        reviewed=_entry_is_reviewed(entry),
        length=get_sequence_length(entry),
        ec_numbers=extract_ec_numbers(entry),
        score=score,
        reasons=reasons,
        sequence=get_sequence_value(entry),
        gene_names=extract_gene_names(entry),
        catalytic_activities=extract_catalytic_activities(entry),
        catalytic_activity_records=extract_catalytic_activity_records(entry),
        cofactors=extract_cofactors(entry),
        subunit=extract_subunit_comments(entry),
        function_comments=extract_comment_texts(entry, "FUNCTION"),
        ptm_comments=extract_comment_texts(entry, "PTM"),
        feature_annotations=extract_feature_annotations(entry),
        domain_ids=extract_domain_ids(entry),
        keywords=extract_keywords(entry),
        protein_existence=get_protein_existence(entry),
        sequence_version=get_sequence_version(entry),
        sequence_sha256=get_sequence_sha256(entry),
        aliases=extract_aliases(entry),
        publication_ids=extract_publication_ids(entry),
        cross_references=extract_cross_references(entry),
        subcellular_locations=extract_subcellular_locations(entry),
        rhea_ids=extract_rhea_ids(entry),
        warnings=warnings,
        score_breakdown=score_breakdown,
    )


def candidate_from_reaction_entry(
    entry: dict[str, Any],
    chassis_key: str,
    *,
    retrieval_strategy: str,
    retrieval_query_id: str,
    matched_rhea_ids: list[str] | None = None,
    candidate_role: str = "catalytic_main",
    component_type: str = "",
    system_anchor_accession: str = "",
    allow_transmembrane: bool = False,
    function_evidence_reason: str = "",
    complex_portal_ac: str = "",
    complex_portal_name: str = "",
    complex_match_basis: str = "",
    complex_evidence_type: str = "",
    component_stoichiometry: dict[str, float | None] | None = None,
) -> ProteinCandidate | None:
    if chassis_key not in CHASSIS_TAXON_PRESETS:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    passed, filter_reasons, filter_warnings = hard_filter_candidate_without_ec(
        entry,
        allow_transmembrane=allow_transmembrane,
    )
    if not passed:
        return None
    entry_ecs = extract_ec_numbers(entry)
    scoring_ec = entry_ecs[0] if entry_ecs else ""
    if scoring_ec:
        score, reasons, warnings, score_breakdown = score_candidate(
            entry=entry,
            target_ec=scoring_ec,
            chassis_profile=CHASSIS_TAXON_PRESETS[chassis_key],
        )
    else:
        evidence_score, evidence_reasons, evidence_warnings = _score_evidence(entry)
        expression_score, expression_reasons, expression_warnings = _score_expression(entry)
        host_score, host_reasons = _score_host(entry, CHASSIS_TAXON_PRESETS[chassis_key])
        function_score = 100.0
        score_breakdown = {
            "function": function_score,
            "evidence": evidence_score,
            "expression": expression_score,
            "host": host_score,
            "total": round(
                function_score * 0.45
                + evidence_score * 0.25
                + expression_score * 0.20
                + host_score * 0.10,
                2,
            ),
        }
        score = score_breakdown["total"]
        reasons = [
            function_evidence_reason or "function: exact Rhea reaction match"
        ] + evidence_reasons + expression_reasons + host_reasons
        warnings = evidence_warnings + expression_warnings
    if function_evidence_reason and function_evidence_reason not in reasons:
        reasons = [function_evidence_reason] + reasons
    candidate = _candidate_from_entry(
        entry=entry,
        score=score,
        reasons=_unique(filter_reasons + reasons),
        warnings=_unique(filter_warnings + warnings),
        score_breakdown=score_breakdown,
    )
    candidate.candidate_role = candidate_role
    candidate.component_type = component_type
    candidate.system_anchor_accession = system_anchor_accession
    candidate.complex_portal_ac = complex_portal_ac
    candidate.complex_portal_name = complex_portal_name
    candidate.complex_match_basis = complex_match_basis
    candidate.complex_evidence_type = complex_evidence_type
    candidate.component_stoichiometry = component_stoichiometry or {}
    candidate.retrieval_strategy = retrieval_strategy
    candidate.retrieval_query_id = retrieval_query_id
    candidate.matched_rhea_ids = _unique(matched_rhea_ids or [])
    candidate.direction_support = "unassessed"
    candidate.reaction_confidence = "medium"
    return candidate


def recommend_uniprot_proteins(
    ec_number: str,
    chassis_key: str = "ecoli_mg1655",
    top_n: int = 10,
    max_results: int = 1000,
    allow_transmembrane: bool = False,
    session: requests.Session | None = None,
) -> list[ProteinCandidate]:
    ec_number = normalize_ec_number(ec_number)
    if chassis_key not in CHASSIS_TAXON_PRESETS:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    chassis_profile = CHASSIS_TAXON_PRESETS[chassis_key]
    http = session or requests.Session()
    reviewed_entries = search_uniprot_by_ec(
        ec_number=ec_number,
        reviewed_only=True,
        max_results=min(max_results, 500),
        session=http,
    )
    all_entries = search_uniprot_by_ec(
        ec_number=ec_number,
        reviewed_only=False,
        max_results=max_results,
        session=http,
    )
    entries = _merge_entries([reviewed_entries, all_entries])
    candidates = []
    for entry in entries:
        passed, filter_reasons, filter_warnings = hard_filter_candidate(
            entry=entry,
            target_ec=ec_number,
            allow_transmembrane=allow_transmembrane,
        )
        if not passed:
            continue
        score, reasons, warnings, score_breakdown = score_candidate(
            entry=entry,
            target_ec=ec_number,
            chassis_profile=chassis_profile,
        )
        candidates.append(_candidate_from_entry(
            entry=entry,
            score=score,
            reasons=_unique(filter_reasons + reasons),
            warnings=_unique(filter_warnings + warnings),
            score_breakdown=score_breakdown,
        ))
    candidates.sort(key=lambda candidate: candidate.score, reverse=True)
    return candidates[:top_n]


def _step_candidate_row(
    requirement: dict[str, Any],
    ec_number: str,
    candidate: ProteinCandidate,
) -> dict[str, Any]:
    return {
        "solution_id": requirement["solution_id"],
        "step_index": requirement["step_index"],
        "reaction_id": requirement["reaction_id"],
        "reaction_name": requirement["reaction_name"],
        "produced_compound_id": requirement["produced_compound_id"],
        "produced_compound_name": requirement["produced_compound_name"],
        "role": "main",
        "ec_number": ec_number,
        "accession": candidate.accession,
        "entry_name": candidate.entry_name,
        "protein_name": candidate.protein_name,
        "organism_name": candidate.organism_name,
        "organism_id": candidate.organism_id or "",
        "taxonomic_lineage": _join(candidate.taxonomic_lineage),
        "reviewed": candidate.reviewed,
        "length": candidate.length or "",
        "score": candidate.score,
        "score_breakdown": _json(candidate.score_breakdown),
        "evaluation_rank": "",
        "candidate_rank": "",
        "selection_status": "",
        "candidate_role": candidate.candidate_role,
        "component_type": candidate.component_type,
        "system_anchor_accession": candidate.system_anchor_accession,
        "complex_portal_ac": candidate.complex_portal_ac,
        "complex_portal_name": candidate.complex_portal_name,
        "complex_match_basis": candidate.complex_match_basis,
        "complex_evidence_type": candidate.complex_evidence_type,
        "component_stoichiometry": _json(candidate.component_stoichiometry),
        "retrieval_strategy": candidate.retrieval_strategy,
        "retrieval_query_id": candidate.retrieval_query_id,
        "matched_rhea_ids": _join(candidate.matched_rhea_ids),
        "matched_ko_ids": _join(candidate.matched_ko_ids),
        "kegg_gene_ids": _join(candidate.kegg_gene_ids),
        "direction_support": candidate.direction_support,
        "direction_verdict": candidate.direction_verdict,
        "direction_confidence": candidate.direction_confidence,
        "direction_evidence_level": candidate.direction_evidence_level,
        "direction_evidence_source_ids": _join(candidate.direction_evidence_source_ids),
        "direction_evidence": " | ".join(candidate.direction_evidence),
        "required_rhea_direction_ids": _join(candidate.required_rhea_direction_ids),
        "reaction_confidence": candidate.reaction_confidence,
        "selenzyme_rank": candidate.selenzyme_rank or "",
        "selenzyme_score": (
            candidate.selenzyme_score if candidate.selenzyme_score is not None else ""
        ),
        "selenzyme_sim_rf": (
            candidate.selenzyme_sim_rf if candidate.selenzyme_sim_rf is not None else ""
        ),
        "selenzyme_sim_2018": (
            candidate.selenzyme_sim_2018
            if candidate.selenzyme_sim_2018 is not None
            else ""
        ),
        "selenzyme_reaction_similarity": (
            candidate.selenzyme_reaction_similarity
            if candidate.selenzyme_reaction_similarity is not None
            else ""
        ),
        "selenzyme_risk_status": candidate.selenzyme_risk_status,
        "selenzyme_matched_reaction_id": candidate.selenzyme_matched_reaction_id,
        "selenzyme_taxonomic_distance": (
            candidate.selenzyme_taxonomic_distance
            if candidate.selenzyme_taxonomic_distance is not None
            else ""
        ),
        "selenzyme_direction_used": candidate.selenzyme_direction_used,
        "selenzyme_direction_preferred": candidate.selenzyme_direction_preferred,
        "gene_names": _join(candidate.gene_names),
        "ec_numbers": _join(candidate.ec_numbers),
        "catalytic_activities": " | ".join(candidate.catalytic_activities),
        "catalytic_activity_records_json": _json(
            candidate.catalytic_activity_records
        ),
        "cofactors": " | ".join(candidate.cofactors),
        "subunit": " | ".join(candidate.subunit),
        "function_comments": " | ".join(candidate.function_comments),
        "ptm_comments": " | ".join(candidate.ptm_comments),
        "feature_annotations": " | ".join(candidate.feature_annotations),
        "domain_ids": _join(candidate.domain_ids),
        "keywords": _join(candidate.keywords),
        "protein_existence": candidate.protein_existence,
        "sequence_version": candidate.sequence_version or "",
        "sequence_sha256": candidate.sequence_sha256,
        "aliases": _join(candidate.aliases),
        "publication_ids": _join(candidate.publication_ids),
        "cross_references": _join(candidate.cross_references),
        "subcellular_locations": " | ".join(candidate.subcellular_locations),
        "rhea_ids": _join(candidate.rhea_ids),
        "reasons": " | ".join(candidate.reasons),
        "warnings": " | ".join(candidate.warnings),
        "sequence": candidate.sequence,
    }
