"""Conservative UniProt identity resolution and candidate conversion."""

from __future__ import annotations

import asyncio
import hashlib
import inspect
import re
from collections.abc import Mapping, Sequence
from typing import Any, Callable

from src.main_protein_selection.literature_activity.models import (
    LiteratureActivityEvidence,
    LiteratureActivityFailure,
)
from src.main_protein_selection.uniprot_protein_candidates import (
    ProteinCandidate,
    resolve_uniprot_identity,
)
from src.main_protein_selection.taxonomy_compatibility import (
    SCORING_WEIGHTS,
    ChassisTaxonomyProfile,
    score_taxonomic_fit,
)


IDENTITY_ADAPTER_VERSION = "literature_activity_identity.v2_taxonomy_ranked"
_UNIPROT_ACCESSION = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]{5}|[A-Z0-9]{10})$",
    re.IGNORECASE,
)


def _unique(values: Sequence[str]) -> list[str]:
    answer: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            answer.append(text)
    return answer


def _field(value: Any, name: str, default: Any = None) -> Any:
    if isinstance(value, Mapping):
        return value.get(name, default)
    return getattr(value, name, default)


def _failure(
    evidence: LiteratureActivityEvidence,
    message: str,
    *,
    retryable: bool,
) -> LiteratureActivityFailure:
    digest = hashlib.sha256(
        f"{evidence.evidence_id}|{message}".encode("utf-8")
    ).hexdigest()[:16]
    return LiteratureActivityFailure(
        failure_id=f"failure_identity_{digest}",
        step_index=evidence.step_index,
        reaction_id=evidence.reaction_id,
        stage="identity",
        source="UniProt",
        query_id=evidence.evidence_id,
        message=message,
        retryable=retryable,
    )


async def _call_resolver(
    resolver: Callable[..., Any],
    **kwargs: Any,
) -> Any:
    if inspect.iscoroutinefunction(resolver):
        return await resolver(**kwargs)
    result = await asyncio.to_thread(resolver, **kwargs)
    return await result if inspect.isawaitable(result) else result


def _publication_ids(evidence: LiteratureActivityEvidence) -> list[str]:
    return _unique([
        f"DOI:{evidence.doi}" if evidence.doi else "",
        f"PubMed:{evidence.pmid}" if evidence.pmid else "",
        f"PMC:{evidence.pmcid}" if evidence.pmcid else "",
    ])


def _candidate_from_hit(
    evidence: LiteratureActivityEvidence,
    hit: Any,
    taxonomy_profile: ChassisTaxonomyProfile,
) -> ProteinCandidate:
    identity_score = float(_field(hit, "identity_score", 0.0) or 0.0)
    level_base = 82.0 if evidence.evidence_level == "A" else 74.0
    function_score = 90.0 if evidence.evidence_level == "A" else 82.0
    evidence_score = 90.0 if evidence.evidence_level == "A" else 75.0
    pseudo_entry = {
        "organism": {
            "scientificName": str(_field(hit, "organism_name", "")),
            "taxonId": _field(hit, "taxon_id"),
            "lineage": [str(item) for item in _field(hit, "lineage", [])],
        },
        "lineages": [
            dict(item)
            for item in _field(hit, "ranked_lineage", [])
            if isinstance(item, Mapping)
        ],
    }
    taxonomy_fit = score_taxonomic_fit(pseudo_entry, taxonomy_profile)
    score = round(
        function_score * SCORING_WEIGHTS["function"]
        + evidence_score * SCORING_WEIGHTS["evidence"]
        + identity_score * SCORING_WEIGHTS["expression"]
        + taxonomy_fit.score * SCORING_WEIGHTS["host"],
        2,
    )
    hit_publications = [str(item) for item in _field(hit, "publication_ids", [])]
    limitations = list(evidence.limitations)
    evidence_reason = (
        f"literature activity evidence {evidence.evidence_id}: "
        f"level {evidence.evidence_level}, assay {evidence.assay_type}"
    )
    warnings = _unique([
        *[str(item) for item in _field(hit, "warnings", [])],
        "literature-derived non-standard activity requires human review",
        *(f"literature limitation: {item}" for item in limitations),
    ])
    confidence = (
        "literature_grade_a"
        if evidence.evidence_level == "A"
        else "literature_grade_b"
    )
    return ProteinCandidate(
        accession=str(_field(hit, "accession", "")).upper(),
        entry_name=str(_field(hit, "entry_name", "")),
        protein_name=str(_field(hit, "protein_name", "")) or evidence.protein_name,
        organism_name=str(_field(hit, "organism_name", "")) or evidence.organism_name,
        organism_id=_field(hit, "taxon_id"),
        reviewed=bool(_field(hit, "reviewed", False)),
        length=_field(hit, "sequence_length"),
        ec_numbers=[],
        score=score,
        reasons=[
            evidence_reason,
            "function: experimentally supported non-standard activity",
            f"literature review status: {evidence.review_status}",
        ],
        sequence=str(_field(hit, "sequence", "")),
        taxonomic_lineage=[str(item) for item in _field(hit, "lineage", [])],
        taxonomic_lineage_ids=[
            int(item.get("taxonId"))
            for item in _field(hit, "ranked_lineage", [])
            if isinstance(item, Mapping) and str(item.get("taxonId") or "").isdigit()
        ],
        taxonomic_shared_taxon_id=taxonomy_fit.shared_taxon_id,
        taxonomic_shared_name=taxonomy_fit.shared_name,
        taxonomic_shared_rank=taxonomy_fit.shared_rank,
        taxonomic_fit_status=taxonomy_fit.status,
        taxonomic_fit_score=taxonomy_fit.score,
        taxonomy_evidence_source=taxonomy_fit.evidence_source,
        gene_names=[str(item) for item in _field(hit, "gene_names", [])],
        protein_existence="",
        sequence_version=_field(hit, "sequence_version"),
        sequence_sha256=str(_field(hit, "sequence_sha256", "")),
        aliases=[str(item) for item in _field(hit, "aliases", [])],
        publication_ids=_unique([*_publication_ids(evidence), *hit_publications]),
        cross_references=[str(item) for item in _field(hit, "cross_references", [])],
        warnings=warnings,
        score_breakdown={
            "function": function_score,
            "evidence": evidence_score,
            "expression": identity_score,
            "host": taxonomy_fit.score,
            "total": score,
        },
        retrieval_strategy="literature_experimental_activity",
        retrieval_query_id=evidence.evidence_id,
        direction_support="literature_experimental",
        direction_verdict="supported",
        direction_confidence="high" if evidence.evidence_level == "A" else "medium",
        direction_evidence_level="experimental_literature",
        direction_evidence_source_ids=[evidence.evidence_id],
        direction_evidence=[
            f"{evidence.source_locator}: {evidence.evidence_summary}"
        ],
        reaction_confidence=confidence,
    )


def _merge_candidates(candidates: Sequence[ProteinCandidate]) -> list[ProteinCandidate]:
    by_accession: dict[str, ProteinCandidate] = {}
    for candidate in candidates:
        accession = candidate.accession.upper()
        existing = by_accession.get(accession)
        if existing is None:
            by_accession[accession] = candidate
            continue
        if candidate.reaction_confidence == "literature_grade_a":
            existing.reaction_confidence = "literature_grade_a"
            existing.direction_confidence = "high"
        retrieval_query_ids = sorted({
            value
            for raw in (
                existing.retrieval_query_id,
                candidate.retrieval_query_id,
            )
            for value in str(raw or "").split(";")
            if value
        })
        existing.retrieval_query_id = ";".join(retrieval_query_ids)
        existing.score = max(existing.score, candidate.score)
        existing.reasons = _unique([*existing.reasons, *candidate.reasons])
        existing.warnings = _unique([*existing.warnings, *candidate.warnings])
        existing.publication_ids = _unique([
            *existing.publication_ids,
            *candidate.publication_ids,
        ])
        existing.direction_evidence_source_ids = _unique([
            *existing.direction_evidence_source_ids,
            *candidate.direction_evidence_source_ids,
        ])
        existing.direction_evidence = _unique([
            *existing.direction_evidence,
            *candidate.direction_evidence,
        ])
    return sorted(
        by_accession.values(),
        key=lambda item: (
            0 if item.reaction_confidence == "literature_grade_a" else 1,
            -item.score,
            0 if item.reviewed else 1,
            item.accession,
        ),
    )


async def resolve_evidence_identities(
    evidence_records: Sequence[LiteratureActivityEvidence],
    *,
    max_results: int,
    allow_transmembrane: bool,
    session: Any = None,
    resolver: Callable[..., Any] = resolve_uniprot_identity,
    taxonomy_profile: ChassisTaxonomyProfile,
) -> tuple[
    list[LiteratureActivityEvidence],
    dict[int, list[ProteinCandidate]],
    list[LiteratureActivityFailure],
]:
    """Resolve only A/B evidence; ambiguity can never create a candidate."""

    updated_records: list[LiteratureActivityEvidence] = []
    candidates_by_step: dict[int, list[ProteinCandidate]] = {}
    failures: list[LiteratureActivityFailure] = []
    for evidence in evidence_records:
        if evidence.evidence_level not in {"A", "B"}:
            updated_records.append(evidence)
            continue
        identifier = evidence.protein_identifier
        identifier_kind = evidence.protein_identifier_kind
        if not identifier:
            updated_records.append(evidence.model_copy(update={
                "identity_status": "no_hit",
                "fit_status": "audit_only",
                "limitations": _unique([
                    *evidence.limitations,
                    "protein identity could not be resolved because no identifier was extracted",
                ]),
            }))
            continue
        if _UNIPROT_ACCESSION.fullmatch(identifier):
            identifier_kind = "accession"
        try:
            resolution = await _call_resolver(
                resolver,
                identifier=identifier,
                identifier_kind=identifier_kind,
                expected_role_terms=(),
                taxon_ids=[evidence.taxon_id] if evidence.taxon_id else (),
                max_results=max_results,
                session=session,
                source_publication_ids=_publication_ids(evidence),
                external_identifiers=(),
                organism_names=[evidence.organism_name] if evidence.organism_name else (),
                strain_names=(),
            )
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            failures.append(_failure(evidence, message, retryable=True))
            updated_records.append(evidence.model_copy(update={
                "identity_status": "source_unavailable",
                "fit_status": "audit_only",
                "limitations": _unique([
                    *evidence.limitations,
                    "UniProt identity service failed; activity evidence was retained for review",
                ]),
            }))
            continue

        status = str(_field(resolution, "status", "source_unavailable"))
        selected = _field(resolution, "selected_hit")
        match_basis = [str(item) for item in _field(resolution, "match_basis", [])]
        query_ids = [str(item) for item in _field(resolution, "query_ids", [])]
        if status != "unique" or selected is None:
            message = f"UniProt identity resolution returned {status} for {identifier}"
            failures.append(_failure(
                evidence,
                message,
                retryable=status == "source_unavailable",
            ))
            normalized_status = (
                status if status in {"ambiguous", "no_hit", "source_unavailable"}
                else "source_unavailable"
            )
            updated_records.append(evidence.model_copy(update={
                "identity_status": normalized_status,
                "identity_match_basis": match_basis,
                "identity_query_ids": query_ids,
                "fit_status": "audit_only",
                "limitations": _unique([
                    *evidence.limitations,
                    f"UniProt identity was not unique ({normalized_status})",
                ]),
            }))
            continue

        hit_warnings = [str(item) for item in _field(selected, "warnings", [])]
        if (
            not allow_transmembrane
            and any("transmembrane region present" in item for item in hit_warnings)
        ):
            updated_records.append(evidence.model_copy(update={
                "resolved_accession": str(_field(selected, "accession", "")).upper(),
                "identity_status": "unique",
                "identity_match_basis": [str(item) for item in _field(selected, "match_basis", match_basis)],
                "identity_query_ids": [str(item) for item in _field(selected, "query_ids", query_ids)],
                "fit_status": "audit_only",
                "limitations": _unique([
                    *evidence.limitations,
                    "resolved protein contains a transmembrane region and allow_transmembrane=False",
                ]),
            }))
            continue

        resolved = evidence.model_copy(update={
            "resolved_accession": str(_field(selected, "accession", "")).upper(),
            "identity_status": "unique",
            "identity_match_basis": [str(item) for item in _field(selected, "match_basis", match_basis)],
            "identity_query_ids": [str(item) for item in _field(selected, "query_ids", query_ids)],
        })
        updated_records.append(resolved)
        candidates_by_step.setdefault(evidence.step_index, []).append(
            _candidate_from_hit(resolved, selected, taxonomy_profile)
        )

    return (
        updated_records,
        {
            step: _merge_candidates(candidates)
            for step, candidates in candidates_by_step.items()
        },
        failures,
    )


__all__ = [
    "IDENTITY_ADAPTER_VERSION",
    "resolve_evidence_identities",
]
