"""UniProt 蛋白候选的检索、解析、过滤与评分。"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import requests

from src.main_protein_selection.sequence_quality import (
    unsupported_amino_acid_reason,
)
from src.main_protein_selection.settings import (
    UNIPROT_HTTP_CONFIG,
    UNIPROT_PAGE_SIZE as CONFIGURED_UNIPROT_PAGE_SIZE,
    UNIPROT_SEARCH_URL as CONFIGURED_UNIPROT_SEARCH_URL,
)
from src.main_protein_selection.taxonomy_compatibility import (
    CHASSIS_TAXON_PRESETS,
    ChassisTaxonomyProfile,
    TaxonomyFit,
    normalize_scoring_weights,
    resolve_chassis_taxonomy,
    score_taxonomic_fit,
)
# UniProt 服务地址和请求参数统一记录在 service_config.py。
UNIPROT_SEARCH_URL = CONFIGURED_UNIPROT_SEARCH_URL
UNIPROT_PAGE_SIZE = CONFIGURED_UNIPROT_PAGE_SIZE
REQUEST_TIMEOUT = UNIPROT_HTTP_CONFIG.timeout_seconds
UNIPROT_ACCESSION_BATCH_SIZE = 50
UNIPROT_FIELDS = ",".join([
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "reviewed",
    "protein_existence",
    "annotation_score",
    "ec",
    "length",
    "lineage",
    "lineage_ids",
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


@dataclass(frozen=True, slots=True)
class UniProtAccessionBatchRecord:
    """One auditable UniProt accession-batch request."""

    query_id: str
    requested_accessions: tuple[str, ...]
    returned_accessions: tuple[str, ...] = ()
    missing_accessions: tuple[str, ...] = ()
    error: str = ""
    recovered_by_split: bool = False


@dataclass(frozen=True, slots=True)
class UniProtAccessionBatchResolution:
    """Result of populating a shared accession-to-entry cache."""

    records: tuple[UniProtAccessionBatchRecord, ...]
    accession_errors: dict[str, str]

    @property
    def query_ids(self) -> list[str]:
        return [record.query_id for record in self.records]

    @property
    def query_errors(self) -> dict[str, str]:
        return {
            record.query_id: record.error
            for record in self.records
            if record.error and not record.recovered_by_split
        }


STEP_CANDIDATE_COLUMNS = [
    "solution_id",
    "step_index",
    "reaction_id",
    "reaction_name",
    "produced_compound_id",
    "produced_compound_name",
    "role",
    "ec_number",
    "accession",
    "entry_name",
    "protein_name",
    "organism_name",
    "organism_id",
    "taxonomic_lineage",
    "taxonomic_lineage_ids",
    "taxonomic_shared_taxon_id",
    "taxonomic_shared_name",
    "taxonomic_shared_rank",
    "taxonomic_fit_status",
    "taxonomic_fit_score",
    "taxonomy_evidence_source",
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
    "taxonomic_shared_taxon_id",
    "taxonomic_shared_name",
    "taxonomic_shared_rank",
    "taxonomic_fit_status",
    "taxonomic_fit_score",
    "taxonomy_evidence_source",
    "reviewed",
    "length",
    "covered_step_indexes",
    "covered_reaction_ids",
    "covered_ec_numbers",
    "roles",
    "best_score",
    "mean_score",
    "min_score",
    "match_count",
    "candidate_roles",
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
    taxonomic_lineage_ids: list[int] = field(default_factory=list)
    taxonomic_shared_taxon_id: int | None = None
    taxonomic_shared_name: str = ""
    taxonomic_shared_rank: str = ""
    taxonomic_fit_status: str = "unknown"
    taxonomic_fit_score: float = 50.0
    taxonomy_evidence_source: str = ""
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


def _normalized_accessions(accessions: Sequence[str]) -> list[str]:
    return list(dict.fromkeys(
        str(accession or "").strip().upper()
        for accession in accessions
        if str(accession or "").strip()
    ))


def _accession_batch_query_id(prefix: str, accessions: Sequence[str]) -> str:
    query_hash = hashlib.sha256(
        ";".join(accessions).encode("utf-8")
    ).hexdigest()[:12]
    return f"{str(prefix or 'uniprot_accession_batch').strip()}_{query_hash}"


def _http_error_status(exc: Exception) -> int | None:
    response = getattr(exc, "response", None)
    status_code = getattr(response, "status_code", None)
    try:
        return int(status_code) if status_code is not None else None
    except (TypeError, ValueError):
        return None


def resolve_uniprot_accession_batches(
    accessions: Sequence[str],
    *,
    session: requests.Session,
    entry_cache: MutableMapping[str, dict[str, Any] | None],
    query_id_prefix: str = "uniprot_accession_batch",
    batch_size: int = UNIPROT_ACCESSION_BATCH_SIZE,
) -> UniProtAccessionBatchResolution:
    """Populate ``entry_cache`` through deterministic accession OR queries.

    Successful no-hits are negative-cached as ``None``. Request failures do
    not mutate the affected cache entries, so a later retrieval stage may
    retry them. Bad-request and URI-too-long responses are bisected to keep
    one malformed accession from hiding the remaining valid entries.
    """

    normalized = _normalized_accessions(accessions)
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    missing = [accession for accession in normalized if accession not in entry_cache]
    records: list[UniProtAccessionBatchRecord] = []
    accession_errors: dict[str, str] = {}

    def resolve_batch(batch: list[str]) -> None:
        if not batch:
            return
        query_id = _accession_batch_query_id(query_id_prefix, batch)
        query = " OR ".join(f"accession:{accession}" for accession in batch)
        try:
            entries = search_uniprot_by_query(
                query,
                reviewed_only=False,
                max_results=len(batch),
                session=session,
            )
        except Exception as exc:
            error = str(exc)
            recovered_by_split = (
                _http_error_status(exc) in {400, 414} and len(batch) > 1
            )
            records.append(UniProtAccessionBatchRecord(
                query_id=query_id,
                requested_accessions=tuple(batch),
                error=error,
                recovered_by_split=recovered_by_split,
            ))
            if recovered_by_split:
                midpoint = len(batch) // 2
                resolve_batch(batch[:midpoint])
                resolve_batch(batch[midpoint:])
                return
            for accession in batch:
                accession_errors[accession] = error
            return

        requested = set(batch)
        returned: list[str] = []
        for entry in entries:
            accession = str(entry.get("primaryAccession") or "").strip().upper()
            if accession not in requested or accession in returned:
                continue
            entry_cache[accession] = entry
            returned.append(accession)
        returned_set = set(returned)
        absent = [accession for accession in batch if accession not in returned_set]
        for accession in absent:
            entry_cache[accession] = None
        records.append(UniProtAccessionBatchRecord(
            query_id=query_id,
            requested_accessions=tuple(batch),
            returned_accessions=tuple(returned),
            missing_accessions=tuple(absent),
        ))

    for index in range(0, len(missing), batch_size):
        resolve_batch(missing[index : index + batch_size])
    return UniProtAccessionBatchResolution(
        records=tuple(records),
        accession_errors=accession_errors,
    )


def extract_ec_numbers(entry: dict[str, Any]) -> list[str]:
    ecs = set()

    def collect(description: dict[str, Any]) -> None:
        names = [description.get("recommendedName", {})]
        names.extend(description.get("alternativeNames", []))
        names.extend(description.get("submissionNames", []))
        for name in names:
            for ec in name.get("ecNumbers", []):
                if ec.get("value"):
                    ecs.add(str(ec["value"]))
        # Multi-domain and processed proteins can place all catalytic ECs
        # on their included domains/chains rather than on the whole protein.
        # UniProt EC searches include these annotations, so local exact-EC
        # filtering must read them as well.
        for branch in ("includes", "contains"):
            for component in description.get(branch, []):
                collect(component)

    collect(entry.get("proteinDescription", {}))
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
    """Return stable protein/gene aliases from database annotations."""

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


def get_lineage_ids(entry: dict[str, Any]) -> list[int]:
    values: list[int] = []
    for item in entry.get("lineages", []):
        if not isinstance(item, Mapping):
            continue
        try:
            taxon_id = int(item.get("taxonId") or 0)
        except (TypeError, ValueError):
            continue
        if taxon_id > 0:
            values.append(taxon_id)
    return list(dict.fromkeys(values))


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
    sequence = get_sequence_value(entry)
    if not sequence:
        return False, ["filtered: missing amino acid sequence"], warnings
    sequence_rejection = unsupported_amino_acid_reason(sequence)
    if sequence_rejection:
        return False, [f"filtered: {sequence_rejection}"], warnings
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
    sequence = get_sequence_value(entry)
    if not sequence:
        return False, ["filtered: missing amino acid sequence"], warnings
    sequence_rejection = unsupported_amino_acid_reason(sequence)
    if sequence_rejection:
        return False, [f"filtered: {sequence_rejection}"], warnings
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


def _score_host(
    entry: dict[str, Any],
    chassis_profile: ChassisTaxonomyProfile,
) -> tuple[float, list[str], TaxonomyFit]:
    fit = score_taxonomic_fit(entry, chassis_profile)
    if fit.status == "unknown":
        reason = "host taxonomy: unavailable; neutral taxonomy score applied"
    elif fit.status in {"cross_domain", "root_only"}:
        reason = "host taxonomy: source shares no ranked lineage with chassis"
    elif fit.status == "exact_host_taxon":
        reason = "host taxonomy: source taxon exactly matches chassis"
    else:
        reason = (
            "host taxonomy: lowest common ancestor "
            f"{fit.shared_name or 'unknown'} ({fit.shared_rank or 'unknown'})"
        )
    return fit.score, [reason], fit


def score_candidate(
    entry: dict[str, Any],
    target_ec: str,
    chassis_profile: ChassisTaxonomyProfile,
    *,
    scoring_weights: Mapping[str, Any] | None = None,
) -> tuple[float, list[str], list[str], dict[str, float], TaxonomyFit]:
    active_weights = normalize_scoring_weights(scoring_weights)
    target_ec = normalize_ec_number(target_ec)
    function_score, function_reasons, function_warnings = _score_function(entry, target_ec)
    evidence_score, evidence_reasons, evidence_warnings = _score_evidence(entry)
    expression_score, expression_reasons, expression_warnings = _score_expression(entry)
    host_score, host_reasons, taxonomy_fit = _score_host(entry, chassis_profile)
    total = (
        function_score * active_weights["function"]
        + evidence_score * active_weights["evidence"]
        + expression_score * active_weights["expression"]
        + host_score * active_weights["host"]
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
        taxonomy_fit,
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
    taxonomy_fit: TaxonomyFit | None = None,
) -> ProteinCandidate:
    organism = entry.get("organism", {})
    taxonomy_fit = taxonomy_fit or TaxonomyFit(
        status="unknown",
        score=50.0,
        candidate_taxon_id=organism.get("taxonId"),
        evidence_source="not_scored",
    )
    return ProteinCandidate(
        accession=entry.get("primaryAccession", ""),
        entry_name=entry.get("uniProtkbId", ""),
        protein_name=get_protein_name(entry),
        organism_name=organism.get("scientificName", ""),
        organism_id=organism.get("taxonId"),
        taxonomic_lineage=[str(value) for value in get_lineage(entry)],
        taxonomic_lineage_ids=get_lineage_ids(entry),
        taxonomic_shared_taxon_id=taxonomy_fit.shared_taxon_id,
        taxonomic_shared_name=taxonomy_fit.shared_name,
        taxonomic_shared_rank=taxonomy_fit.shared_rank,
        taxonomic_fit_status=taxonomy_fit.status,
        taxonomic_fit_score=taxonomy_fit.score,
        taxonomy_evidence_source=taxonomy_fit.evidence_source,
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
    allow_transmembrane: bool = False,
    function_evidence_reason: str = "",
    taxonomy_profile: ChassisTaxonomyProfile | None = None,
    scoring_weights: Mapping[str, Any] | None = None,
) -> ProteinCandidate | None:
    active_weights = normalize_scoring_weights(scoring_weights)
    if chassis_key not in CHASSIS_TAXON_PRESETS:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    profile = taxonomy_profile or resolve_chassis_taxonomy(
        chassis_key,
        allow_network=False,
    )
    passed, filter_reasons, filter_warnings = hard_filter_candidate_without_ec(
        entry,
        allow_transmembrane=allow_transmembrane,
    )
    if not passed:
        return None
    entry_ecs = extract_ec_numbers(entry)
    scoring_ec = entry_ecs[0] if entry_ecs else ""
    if scoring_ec:
        score, reasons, warnings, score_breakdown, taxonomy_fit = score_candidate(
            entry=entry,
            target_ec=scoring_ec,
            chassis_profile=profile,
            scoring_weights=active_weights,
        )
    else:
        evidence_score, evidence_reasons, evidence_warnings = _score_evidence(entry)
        expression_score, expression_reasons, expression_warnings = _score_expression(entry)
        host_score, host_reasons, taxonomy_fit = _score_host(entry, profile)
        function_score = 100.0
        score_breakdown = {
            "function": function_score,
            "evidence": evidence_score,
            "expression": expression_score,
            "host": host_score,
            "total": round(
                function_score * active_weights["function"]
                + evidence_score * active_weights["evidence"]
                + expression_score * active_weights["expression"]
                + host_score * active_weights["host"],
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
        taxonomy_fit=taxonomy_fit,
    )
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
    taxonomy_profile: ChassisTaxonomyProfile | None = None,
    scoring_weights: Mapping[str, Any] | None = None,
) -> list[ProteinCandidate]:
    active_weights = normalize_scoring_weights(scoring_weights)
    ec_number = normalize_ec_number(ec_number)
    if chassis_key not in CHASSIS_TAXON_PRESETS:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    http = session or requests.Session()
    chassis_profile = taxonomy_profile or resolve_chassis_taxonomy(
        chassis_key,
        session=http,
    )
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
        score, reasons, warnings, score_breakdown, taxonomy_fit = score_candidate(
            entry=entry,
            target_ec=ec_number,
            chassis_profile=chassis_profile,
            scoring_weights=active_weights,
        )
        candidates.append(_candidate_from_entry(
            entry=entry,
            score=score,
            reasons=_unique(filter_reasons + reasons),
            warnings=_unique(filter_warnings + warnings),
            score_breakdown=score_breakdown,
            taxonomy_fit=taxonomy_fit,
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
        "taxonomic_lineage_ids": _join(
            [str(value) for value in candidate.taxonomic_lineage_ids]
        ),
        "taxonomic_shared_taxon_id": candidate.taxonomic_shared_taxon_id or "",
        "taxonomic_shared_name": candidate.taxonomic_shared_name,
        "taxonomic_shared_rank": candidate.taxonomic_shared_rank,
        "taxonomic_fit_status": candidate.taxonomic_fit_status,
        "taxonomic_fit_score": candidate.taxonomic_fit_score,
        "taxonomy_evidence_source": candidate.taxonomy_evidence_source,
        "reviewed": candidate.reviewed,
        "length": candidate.length or "",
        "score": candidate.score,
        "score_breakdown": _json(candidate.score_breakdown),
        "evaluation_rank": "",
        "candidate_rank": "",
        "selection_status": "",
        "candidate_role": candidate.candidate_role,
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
