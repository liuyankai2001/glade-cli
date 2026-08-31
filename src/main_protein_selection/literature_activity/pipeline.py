"""Public synchronous pipeline for literature-backed enzyme candidates."""

from __future__ import annotations

import asyncio
import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable

from pydantic import ValidationError

from src.main_protein_selection.literature_activity.extractor import (
    build_default_extractor,
    literature_model_identity,
)
from src.main_protein_selection.literature_activity.identity import (
    resolve_evidence_identities,
)
from src.main_protein_selection.literature_activity.models import (
    ActivityExtractor,
    LiteratureActivityArtifact,
    LiteratureActivityEvidence,
    LiteratureActivityFailure,
    LiteratureActivityRequirement,
    LiteratureActivitySearchResult,
    LiteratureActivitySummary,
    LiteratureRetriever,
)
from src.main_protein_selection.literature_activity.query_builder import (
    build_literature_queries,
    literature_component_versions,
    normalize_requirements,
    request_fingerprint,
)
from src.main_protein_selection.literature_activity.storage import (
    finalize_artifact,
    load_cached_result,
    store_cached_result,
    utc_now_text,
    write_artifact,
)
from src.main_protein_selection.literature_activity.validator import (
    validate_activity_claim,
)
from src.main_protein_selection.uniprot_protein_candidates import (
    resolve_uniprot_identity,
)
from src.main_protein_selection.taxonomy_compatibility import (
    ChassisTaxonomyProfile,
    resolve_chassis_taxonomy,
)


MAX_EXTRACTION_PAPERS_PER_STEP = 6
MAX_FULL_TEXT_PAPERS_PER_STEP = 3
EXTRACTION_CONCURRENCY = 3


def _retryable_structured_output_error(exc: Exception) -> bool:
    if isinstance(exc, (ValidationError, json.JSONDecodeError)):
        return True
    class_name = type(exc).__name__.casefold()
    if any(term in class_name for term in (
        "validationerror", "jsondecode", "outputparser",
        "apiconnectionerror", "apitimeouterror",
    )):
        return True
    message = str(exc).casefold()
    return bool(
        re.search(r"\b(?:invalid|malformed)\s+json\b", message)
        or re.search(r"\beof\b.*\b(?:json|pars\w*|value)\b", message)
        or re.search(r"\b(?:json|structured output)\b.*\b(?:pars\w*|decod\w*|valid\w*)\b", message)
        or re.search(r"\b(?:api |network )?connection\b.*\b(?:error|failed|reset|closed)\b", message)
    )


def _paper_prefilter_text(value: str) -> str:
    translated = str(value or "").translate(str.maketrans({
        "−": "-", "–": "-", "—": "-", "＋": "+",
    }))
    return " ".join(
        re.sub(r"[^\w]+", " ", translated.casefold()).split()
    )


def _core_compound_terms(compounds: Sequence[Any]) -> list[str]:
    terms: list[str] = []
    for item in compounds:
        if not item.core:
            continue
        name = re.sub(
            r"^\s*\([+\-RSDL]\)\s*-\s*",
            "",
            str(item.name or ""),
            flags=re.IGNORECASE,
        )
        normalized = _paper_prefilter_text(name or item.compound_id)
        if normalized:
            terms.append(normalized)
    return list(dict.fromkeys(terms))


def _papers_for_extraction(
    requirement: LiteratureActivityRequirement,
    papers: Sequence[Any],
) -> list[Any]:
    """Preserve retriever rank while cheaply rejecting irrelevant papers."""

    substrate_terms = _core_compound_terms(requirement.substrates)
    product_terms = _core_compound_terms(requirement.products)
    eligible: list[Any] = []
    for paper in papers:
        raw_text = str(paper.source_text or "")
        text = _paper_prefilter_text(raw_text)
        product_match = any(term in text for term in product_terms)
        substrate_match = any(term in text for term in substrate_terms)
        reaction_context = substrate_match or bool(re.search(
            r"\b(?:reaction|activity|cataly\w*|convert\w*|dephosphorylat\w*|"
            r"hydroly\w*|biosynth\w*|bioconversion|production|formation)\b",
            text,
        ))
        protein_context = bool(re.search(
            r"\b(?:protein|enzyme|gene|phosphatase|hydrolase|reductase|synthase|"
            r"overexpress\w*|express\w*|knockout|complement\w*|recombinant|purified)\b",
            text,
        )) or bool(re.search(r"\b[A-Za-z]{2,}\d*[A-Z]\b", raw_text))
        if product_match and reaction_context and protein_context:
            eligible.append(paper)
    selected = eligible if eligible else list(papers)
    return selected[:MAX_EXTRACTION_PAPERS_PER_STEP]


def _failure(
    requirement: LiteratureActivityRequirement,
    *,
    stage: str,
    message: str,
    source: str,
    query_id: str = "",
    retryable: bool,
) -> LiteratureActivityFailure:
    digest = hashlib.sha256(
        f"{requirement.step_index}|{stage}|{source}|{query_id}|{message}".encode("utf-8")
    ).hexdigest()[:16]
    return LiteratureActivityFailure(
        failure_id=f"failure_{stage}_{digest}",
        step_index=requirement.step_index,
        reaction_id=requirement.reaction_id,
        stage=stage,
        source=source,
        query_id=query_id,
        message=message,
        retryable=retryable,
    )


def _summary(
    requirements: Sequence[LiteratureActivityRequirement],
    *,
    query_count: int,
    paper_count: int,
    evidence: Sequence[LiteratureActivityEvidence],
    candidates_by_step: Mapping[int, Sequence[Any]],
    failures: Sequence[LiteratureActivityFailure],
) -> LiteratureActivitySummary:
    eligible_steps = sorted(
        step for step, candidates in candidates_by_step.items() if candidates
    )
    required_steps = [item.step_index for item in requirements]
    return LiteratureActivitySummary(
        requirement_count=len(requirements),
        query_count=query_count,
        paper_count=paper_count,
        evidence_count=len(evidence),
        grade_a_count=sum(item.evidence_level == "A" for item in evidence),
        grade_b_count=sum(item.evidence_level == "B" for item in evidence),
        grade_c_count=sum(item.evidence_level == "C" for item in evidence),
        rejected_count=sum(item.evidence_level == "Reject" for item in evidence),
        candidate_count=sum(len(items) for items in candidates_by_step.values()),
        eligible_step_indexes=eligible_steps,
        unresolved_step_indexes=[
            step for step in required_steps if step not in eligible_steps
        ],
        failure_count=len(failures),
    )


def _status(
    *,
    candidate_count: int,
    evidence_count: int,
    successful_query_count: int,
    failure_count: int,
) -> str:
    if successful_query_count == 0 and failure_count:
        return "source_unavailable"
    if candidate_count:
        return "partial" if failure_count else "complete"
    if failure_count:
        return "partial"
    return "not_found"


def _disabled_request_fingerprint(
    requirements: Sequence[Mapping[str, Any] | LiteratureActivityRequirement],
    *,
    chassis_key: str,
    top_n: int,
    max_results: int,
    allow_transmembrane: bool,
) -> str:
    """Fingerprint disabled input without applying literature-only validation."""

    raw_requirements = [
        item.model_dump(mode="json")
        if isinstance(item, LiteratureActivityRequirement)
        else dict(item)
        for item in requirements
    ]
    canonical = json.dumps(
        {
            "policy": "literature_activity_disabled.v1",
            "component_versions": literature_component_versions(),
            "requirements": raw_requirements,
            "chassis_key": chassis_key,
            "top_n": top_n,
            "max_results": max_results,
            "allow_transmembrane": allow_transmembrane,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


async def _run_enabled(
    requirements: list[LiteratureActivityRequirement],
    *,
    chassis_key: str,
    top_n: int,
    max_results: int,
    allow_transmembrane: bool,
    session: Any,
    retriever: LiteratureRetriever | None,
    extractor: ActivityExtractor | None,
    identity_resolver: Callable[..., Any],
    model: Any,
    env_path: str | Path | None,
    taxonomy_profile: ChassisTaxonomyProfile,
) -> tuple[
    list[Any],
    list[LiteratureActivityEvidence],
    dict[int, list[Any]],
    list[LiteratureActivityFailure],
    int,
]:
    if not requirements:
        return [], [], {}, [], 0
    if retriever is None:
        from src.main_protein_selection.literature_activity.retriever import (
            ToolUniverseLiteratureRetriever,
        )

        retriever = ToolUniverseLiteratureRetriever()

    all_queries: list[Any] = []
    failures: list[LiteratureActivityFailure] = []
    papers_by_step: dict[int, list[Any]] = {}
    paper_count = 0
    for requirement in requirements:
        search_queries = build_literature_queries(
            requirement,
            chassis_key=chassis_key,
        )
        try:
            batch = await retriever.retrieve(
                requirement,
                search_queries,
                max_papers=min(20, max_results),
                max_full_texts=MAX_FULL_TEXT_PAPERS_PER_STEP,
            )
        except Exception as exc:
            failures.append(_failure(
                requirement,
                stage="retrieval",
                message=f"{type(exc).__name__}: {exc}",
                source=type(retriever).__name__,
                retryable=True,
            ))
            papers_by_step[requirement.step_index] = []
            continue
        all_queries.extend(batch.queries)
        failures.extend(batch.failures)
        papers_by_step[requirement.step_index] = _papers_for_extraction(
            requirement,
            batch.papers,
        )
        paper_count += len(batch.papers)

    if extractor is None and any(papers_by_step.values()):
        try:
            extractor = build_default_extractor(model=model, env_path=env_path)
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            for requirement in requirements:
                if papers_by_step.get(requirement.step_index):
                    failures.append(_failure(
                        requirement,
                        stage="extraction",
                        message=message,
                        source="configured language model",
                        retryable=False,
                    ))

    evidence_by_id: dict[str, LiteratureActivityEvidence] = {}
    if extractor is not None:
        semaphore = asyncio.Semaphore(EXTRACTION_CONCURRENCY)

        async def extract_one(requirement, paper):
            async with semaphore:
                for attempt in range(2):
                    try:
                        return await extractor.extract(requirement, paper), None
                    except Exception as exc:
                        if attempt == 0 and _retryable_structured_output_error(exc):
                            continue
                        return None, _failure(
                            requirement,
                            stage="extraction",
                            message=f"{type(exc).__name__}: {exc}",
                            source=paper.paper_id,
                            query_id=paper.paper_id,
                            retryable=_retryable_structured_output_error(exc),
                        )
                raise AssertionError("unreachable extraction retry loop")

        for requirement in requirements:
            papers = papers_by_step.get(requirement.step_index, [])
            extracted = await asyncio.gather(*(
                extract_one(requirement, paper) for paper in papers
            )) if papers else []
            for paper, (result, extraction_failure) in zip(papers, extracted):
                if extraction_failure is not None:
                    failures.append(extraction_failure)
                    continue
                for claim in result.claims:
                    validated = validate_activity_claim(requirement, paper, claim)
                    evidence_by_id.setdefault(validated.evidence_id, validated)

    evidence, candidates_by_step, identity_failures = await resolve_evidence_identities(
        list(evidence_by_id.values()),
        max_results=min(25, max_results),
        allow_transmembrane=allow_transmembrane,
        session=session,
        resolver=identity_resolver,
        taxonomy_profile=taxonomy_profile,
    )
    failures.extend(identity_failures)
    candidates_by_step = {
        step: candidates[:top_n]
        for step, candidates in candidates_by_step.items()
    }
    return all_queries, evidence, candidates_by_step, failures, paper_count


def _run_coroutine(coroutine):
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coroutine)
    # The public API is intentionally synchronous for the CLI. If an embedding
    # application already owns an event loop, isolate our loop in one worker.
    with ThreadPoolExecutor(max_workers=1) as executor:
        return executor.submit(asyncio.run, coroutine).result()


def run_literature_activity_search(
    requirements: Sequence[Mapping[str, Any] | LiteratureActivityRequirement],
    *,
    enabled: bool,
    output_dir: str | Path,
    cache_dir: str | Path,
    chassis_key: str,
    top_n: int = 5,
    max_results: int = 25,
    allow_transmembrane: bool = False,
    session: Any = None,
    retriever: LiteratureRetriever | None = None,
    extractor: ActivityExtractor | None = None,
    identity_resolver: Callable[..., Any] = resolve_uniprot_identity,
    model: Any = None,
    env_path: str | Path | None = None,
    taxonomy_profile: ChassisTaxonomyProfile | None = None,
) -> LiteratureActivitySearchResult:
    """Search unresolved steps and return literature-backed candidates.

    Disabled mode is deliberately handled before construction of the default
    retriever/extractor. It writes an empty audit artifact but cannot load a
    model, start ToolUniverse, call UniProt, or perform any other network I/O.
    """

    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    if max_results < 1:
        raise ValueError("max_results must be at least 1")
    if not enabled:
        fingerprint = _disabled_request_fingerprint(
            requirements,
            chassis_key=chassis_key,
            top_n=top_n,
            max_results=max_results,
            allow_transmembrane=allow_transmembrane,
        )
        summary = _summary(
            [],
            query_count=0,
            paper_count=0,
            evidence=[],
            candidates_by_step={},
            failures=[],
        )
        artifact = finalize_artifact({
            "generated_at": utc_now_text(),
            "status": "disabled",
            "enabled": False,
            "cache_hit": False,
            "request_fingerprint": fingerprint,
            "chassis_key": chassis_key,
            "component_versions": literature_component_versions(),
            "requirements": [],
            "queries": [],
            "evidence": [],
            "failures": [],
            "summary": summary,
        })
        json_path, csv_path = write_artifact(artifact, output_dir)
        return LiteratureActivitySearchResult(
            status="disabled",
            candidates_by_step={},
            artifact=artifact,
            json_path=json_path,
            csv_path=csv_path,
            query_errors={},
        )

    normalized = normalize_requirements(requirements)
    taxonomy_profile = taxonomy_profile or resolve_chassis_taxonomy(
        chassis_key,
        session=session,
        cache_root=cache_dir,
    )
    model_cache_identity = literature_model_identity(
        model=model,
        env_path=env_path,
    )
    component_versions = literature_component_versions(
        model_identity=model_cache_identity
    )
    fingerprint = request_fingerprint(
        normalized,
        chassis_key=chassis_key,
        top_n=top_n,
        max_results=max_results,
        allow_transmembrane=allow_transmembrane,
        model_identity=model_cache_identity,
    )
    if not normalized:
        summary = _summary(
            [],
            query_count=0,
            paper_count=0,
            evidence=[],
            candidates_by_step={},
            failures=[],
        )
        artifact = finalize_artifact({
            "generated_at": utc_now_text(),
            "status": "not_needed",
            "enabled": True,
            "cache_hit": False,
            "request_fingerprint": fingerprint,
            "chassis_key": chassis_key,
            "component_versions": component_versions,
            "requirements": [],
            "queries": [],
            "evidence": [],
            "failures": [],
            "summary": summary,
        })
        json_path, csv_path = write_artifact(artifact, output_dir)
        return LiteratureActivitySearchResult(
            status="not_needed",
            candidates_by_step={},
            artifact=artifact,
            json_path=json_path,
            csv_path=csv_path,
            query_errors={},
        )

    cached = load_cached_result(
        cache_dir=cache_dir,
        output_dir=output_dir,
        request_fingerprint=fingerprint,
    )
    if cached is not None:
        return cached

    queries, evidence, candidates_by_step, failures, paper_count = _run_coroutine(
        _run_enabled(
            normalized,
            chassis_key=chassis_key,
            top_n=top_n,
            max_results=max_results,
            allow_transmembrane=allow_transmembrane,
            session=session,
            retriever=retriever,
            extractor=extractor,
            identity_resolver=identity_resolver,
            model=model,
            env_path=env_path,
            taxonomy_profile=taxonomy_profile,
        )
    )
    successful_queries = sum(
        item.status in {"success", "not_found", "cache_hit"} for item in queries
    )
    candidate_count = sum(len(items) for items in candidates_by_step.values())
    status = _status(
        candidate_count=candidate_count,
        evidence_count=len(evidence),
        successful_query_count=successful_queries,
        failure_count=len(failures),
    )
    summary = _summary(
        normalized,
        query_count=len(queries),
        paper_count=paper_count,
        evidence=evidence,
        candidates_by_step=candidates_by_step,
        failures=failures,
    )
    artifact = finalize_artifact({
        "generated_at": utc_now_text(),
        "status": status,
        "enabled": True,
        "cache_hit": False,
        "request_fingerprint": fingerprint,
        "chassis_key": chassis_key,
        "component_versions": component_versions,
        "requirements": normalized,
        "queries": queries,
        "evidence": evidence,
        "failures": failures,
        "summary": summary,
    })
    fingerprint_reason = f"literature_artifact_sha256:{artifact.artifact_fingerprint}"
    for candidates in candidates_by_step.values():
        for candidate in candidates:
            if fingerprint_reason not in candidate.reasons:
                candidate.reasons.append(fingerprint_reason)
    json_path, csv_path = write_artifact(artifact, output_dir)
    query_errors = {
        failure.query_id or failure.failure_id: failure.message
        for failure in failures
    }
    result = LiteratureActivitySearchResult(
        status=status,
        candidates_by_step=candidates_by_step,
        artifact=artifact,
        json_path=json_path,
        csv_path=csv_path,
        query_errors=query_errors,
    )
    store_cached_result(result, cache_dir=cache_dir)
    return result


def write_source_unavailable_artifact(
    requirements: Sequence[Mapping[str, Any] | LiteratureActivityRequirement],
    *,
    output_dir: str | Path,
    chassis_key: str,
    message: str,
    top_n: int = 5,
    max_results: int = 25,
    allow_transmembrane: bool = False,
) -> LiteratureActivitySearchResult:
    """Overwrite stale literature outputs after an unexpected pipeline failure.

    This emergency audit path deliberately avoids requirement normalization,
    cache access, model construction, ToolUniverse and every network client.
    """

    error_message = str(message or "").strip()
    if not error_message:
        raise ValueError("message must not be empty")
    base_fingerprint = _disabled_request_fingerprint(
        requirements,
        chassis_key=chassis_key,
        top_n=top_n,
        max_results=max_results,
        allow_transmembrane=allow_transmembrane,
    )
    fingerprint = hashlib.sha256(
        f"literature_activity_source_unavailable.v1|{base_fingerprint}".encode("utf-8")
    ).hexdigest()
    digest = hashlib.sha256(error_message.encode("utf-8")).hexdigest()[:16]
    failure = LiteratureActivityFailure(
        failure_id=f"failure_pipeline_{digest}",
        step_index=1,
        reaction_id="UNKNOWN",
        stage="retrieval",
        source="literature_activity_pipeline",
        query_id="literature_activity_pipeline",
        message=error_message,
        retryable=True,
    )
    summary = _summary(
        [],
        query_count=0,
        paper_count=0,
        evidence=[],
        candidates_by_step={},
        failures=[failure],
    )
    artifact = finalize_artifact({
        "generated_at": utc_now_text(),
        "status": "source_unavailable",
        "enabled": True,
        "cache_hit": False,
        "request_fingerprint": fingerprint,
        "chassis_key": chassis_key,
        "component_versions": literature_component_versions(),
        "requirements": [],
        "queries": [],
        "evidence": [],
        "failures": [failure],
        "summary": summary,
    })
    json_path, csv_path = write_artifact(artifact, output_dir)
    return LiteratureActivitySearchResult(
        status="source_unavailable",
        candidates_by_step={},
        artifact=artifact,
        json_path=json_path,
        csv_path=csv_path,
        query_errors={"literature_activity_pipeline": error_message},
    )


__all__ = [
    "run_literature_activity_search",
    "write_source_unavailable_artifact",
]
