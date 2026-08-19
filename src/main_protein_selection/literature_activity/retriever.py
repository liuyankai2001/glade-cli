"""Bounded ToolUniverse retrieval for experimental enzyme-activity papers."""

from __future__ import annotations

import asyncio
import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from src.main_protein_selection.literature_activity.models import (
    LiteratureActivityFailure,
    LiteratureActivityRequirement,
    LiteratureQueryAudit,
    LiteratureRetrievalBatch,
    LiteratureSearchQuery,
    RetrievedLiteraturePaper,
)


RETRIEVER_VERSION = "tooluniverse_literature_activity.v3"
_PMID_PATTERN = re.compile(r"(?<!\d)\d{5,9}(?!\d)")
_PMCID_PATTERN = re.compile(r"\bPMC\d+\b", re.IGNORECASE)
_DOI_PATTERN = re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)
_RETRYABLE_PATTERN = re.compile(
    r"timeout|temporar|rate.?limit|429|502|503|504|connection|reset",
    re.IGNORECASE,
)


@dataclass(slots=True)
class _ToolCallResult:
    audit: LiteratureQueryAudit
    content: str
    raw_evidence_id: str
    failure: LiteratureActivityFailure | None = None


def _unique(values: Sequence[str]) -> list[str]:
    answer: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        key = text.casefold()
        if text and key not in seen:
            seen.add(key)
            answer.append(text)
    return answer


def _normalise_doi(value: Any) -> str:
    text = str(value or "").strip()
    text = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", text, flags=re.I)
    match = _DOI_PATTERN.search(text)
    return match.group(0).rstrip(".,;)").lower() if match else ""


def _normalise_pmid(value: Any) -> str:
    text = str(value or "").strip()
    match = _PMID_PATTERN.fullmatch(text.removeprefix("PMID:").strip())
    return match.group(0) if match else ""


def _normalise_pmcid(value: Any) -> str:
    match = _PMCID_PATTERN.search(str(value or ""))
    return match.group(0).upper() if match else ""


def _serialize_tool_result(raw: Any) -> str:
    if raw is None:
        return ""
    if isinstance(raw, str):
        return raw
    if hasattr(raw, "model_dump"):
        try:
            raw = raw.model_dump(mode="json")
        except TypeError:
            raw = raw.model_dump()
    if isinstance(raw, Mapping):
        content = raw.get("content")
        if isinstance(content, list):
            texts = [
                str(item.get("text") or "")
                for item in content
                if isinstance(item, Mapping) and item.get("type") == "text"
            ]
            if any(texts):
                return "\n".join(texts)
        return json.dumps(raw, ensure_ascii=False, default=str)
    if isinstance(raw, Sequence) and not isinstance(raw, (bytes, bytearray)):
        blocks = list(raw)
        if blocks and all(
            isinstance(item, Mapping)
            and item.get("type") == "text"
            and isinstance(item.get("text"), str)
            for item in blocks
        ):
            return "\n".join(str(item["text"]) for item in blocks)
        return json.dumps(blocks, ensure_ascii=False, default=str)
    return str(raw)


def _mapping_from_raw(raw: Any) -> Mapping[str, Any] | None:
    if isinstance(raw, Mapping):
        return raw
    if hasattr(raw, "model_dump"):
        try:
            dumped = raw.model_dump(mode="json")
        except TypeError:
            dumped = raw.model_dump()
        return dumped if isinstance(dumped, Mapping) else None
    return None


def _error_flag(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().casefold() in {"1", "true", "yes"}


def _envelope_error(mapping: Mapping[str, Any]) -> str:
    status = str(mapping.get("status") or "").strip().casefold()
    flagged = _error_flag(mapping.get("isError")) or _error_flag(
        mapping.get("is_error")
    )
    error_value = mapping.get("error")
    has_error_value = (
        error_value is not None
        and error_value is not False
        and error_value != ""
        and error_value != 0
    )
    if status not in {
        "error",
        "failed",
        "failure",
        "timeout",
        "timed_out",
        "timed out",
    } and not flagged and not has_error_value:
        return ""
    message = (
        mapping.get("message")
        or mapping.get("detail")
        or error_value
        or status
        or "external tool returned an error envelope"
    )
    if isinstance(message, Mapping):
        message = message.get("message") or message.get("detail") or str(message)
    return _text_value(message) or "external tool returned an error envelope"


def _tool_result_error(raw: Any, content: str) -> str:
    """Recognize MCP/tool error envelopes without misclassifying empty hits."""

    raw_mapping = _mapping_from_raw(raw)
    if raw_mapping is not None:
        message = _envelope_error(raw_mapping)
        if message:
            return message
        nested = raw_mapping.get("result")
        if isinstance(nested, Mapping):
            message = _envelope_error(nested)
            if message:
                return message

    decoded = _decode_json(content)
    if isinstance(decoded, Mapping):
        message = _envelope_error(decoded)
        if message:
            return message
        nested = decoded.get("result")
        if isinstance(nested, Mapping):
            message = _envelope_error(nested)
            if message:
                return message

    prefix = re.match(
        r"^\s*(?:error|failed|failure|timeout|timed\s+out)\b\s*[:\-]?\s*(.*)$",
        str(content or ""),
        re.IGNORECASE | re.DOTALL,
    )
    if prefix is not None:
        return prefix.group(0).strip()[:1000]
    return ""


def _decode_json(content: str) -> Any:
    text = str(content or "").strip()
    if not text:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        match = re.search(r"(?:\{[\s\S]*\}|\[[\s\S]*\])", text)
        if match is None:
            return None
        try:
            return json.loads(match.group(0))
        except json.JSONDecodeError:
            return None


def _walk_mappings(value: Any):
    if isinstance(value, Mapping):
        yield value
        for nested in value.values():
            yield from _walk_mappings(nested)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for nested in value:
            yield from _walk_mappings(nested)


def _mapping_value(mapping: Mapping[str, Any], *names: str) -> Any:
    lowered = {str(key).casefold(): value for key, value in mapping.items()}
    for name in names:
        if name.casefold() in lowered:
            return lowered[name.casefold()]
    return None


def _text_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return " ".join(value.split())
    if isinstance(value, Mapping):
        for key in ("value", "text", "name", "title", "abstract"):
            candidate = value.get(key)
            if isinstance(candidate, str):
                return " ".join(candidate.split())
        return ""
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        return " ".join(filter(None, (_text_value(item) for item in value)))
    return str(value).strip()


def _authors(value: Any) -> list[str]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        text = _text_value(value)
        return [text] if text else []
    answer: list[str] = []
    for item in value:
        if isinstance(item, Mapping):
            name = _text_value(
                item.get("name")
                or item.get("full_name")
                or " ".join(
                    filter(None, [
                        _text_value(item.get("given")),
                        _text_value(item.get("family")),
                    ])
                )
            )
        else:
            name = _text_value(item)
        if name:
            answer.append(name)
    return _unique(answer)


def _year(value: Any) -> int | None:
    match = re.search(r"\b(?:18|19|20|21)\d{2}\b", _text_value(value))
    return int(match.group(0)) if match else None


def _paper_key(*, doi: str, pmid: str, pmcid: str, title: str) -> str:
    if doi:
        return f"doi:{doi}"
    if pmid:
        return f"pmid:{pmid}"
    if pmcid:
        return f"pmcid:{pmcid}"
    normalized_title = " ".join(re.sub(r"[^\w]+", " ", title.casefold()).split())
    if normalized_title:
        return "title:" + hashlib.sha256(normalized_title.encode("utf-8")).hexdigest()[:20]
    return ""


def _paper_from_mapping(
    mapping: Mapping[str, Any],
    *,
    source: str,
    raw_evidence_id: str,
) -> RetrievedLiteraturePaper | None:
    title = _text_value(_mapping_value(mapping, "title", "article_title", "name"))
    abstract = _text_value(
        _mapping_value(mapping, "abstract", "abstractText", "abstract_text", "summary")
    )
    pmid = _normalise_pmid(
        _mapping_value(mapping, "pmid", "pubmed_id", "pubmedid", "uid")
    )
    pmcid = _normalise_pmcid(
        _mapping_value(mapping, "pmcid", "pmc_id", "pmc", "fullTextId")
    )
    doi = _normalise_doi(
        _mapping_value(mapping, "doi", "DOI", "digital_object_identifier")
    )
    if not pmid:
        direct = re.search(r"\bPMID\s*:?\s*(\d{5,9})\b", json.dumps(mapping, default=str), re.I)
        pmid = direct.group(1) if direct else ""
    if not doi:
        doi_match = _DOI_PATTERN.search(json.dumps(mapping, default=str))
        doi = _normalise_doi(doi_match.group(0)) if doi_match else ""
    if not pmcid:
        pmcid = _normalise_pmcid(json.dumps(mapping, default=str))
    if not (title and (pmid or pmcid or doi or abstract)):
        return None
    key = _paper_key(doi=doi, pmid=pmid, pmcid=pmcid, title=title)
    if not key:
        return None
    url_values = _mapping_value(mapping, "url", "source_url", "links", "urls")
    urls = re.findall(r"https?://[^\s\"'<>]+", _text_value(url_values))
    journal = _text_value(
        _mapping_value(mapping, "journal", "journal_title", "container-title", "source")
    )
    year = _year(
        _mapping_value(mapping, "year", "publication_year", "pubdate", "published", "date")
    )
    source_name = "Europe PMC" if source.startswith("Europe") else "PubMed"
    return RetrievedLiteraturePaper(
        paper_id=key,
        title=title,
        abstract=abstract,
        authors=_authors(_mapping_value(mapping, "authors", "authorList", "author")),
        journal=journal,
        year=year,
        doi=doi,
        pmid=pmid,
        pmcid=pmcid,
        source_databases=[source_name],
        source_urls=_unique(urls),
        access_level="abstract_only" if abstract else "metadata_only",
        raw_evidence_ids=[raw_evidence_id],
        source_text_sha256=hashlib.sha256(
            "\n".join([title, abstract]).encode("utf-8")
        ).hexdigest(),
    )


def _papers_from_content(
    content: str,
    *,
    source: str,
    raw_evidence_id: str,
) -> list[RetrievedLiteraturePaper]:
    decoded = _decode_json(content)
    papers: list[RetrievedLiteraturePaper] = []
    seen: set[str] = set()
    if decoded is not None:
        for mapping in _walk_mappings(decoded):
            paper = _paper_from_mapping(
                mapping,
                source=source,
                raw_evidence_id=raw_evidence_id,
            )
            if paper is not None and paper.paper_id not in seen:
                seen.add(paper.paper_id)
                papers.append(paper)
    if papers:
        return papers

    # Conservative fallback: IDs without a title are sufficient only to fetch
    # PubMed detail; they are never passed directly to the extraction model.
    for pmid in _unique(re.findall(r"\bPMID\s*:?\s*(\d{5,9})\b", content, re.I)):
        papers.append(RetrievedLiteraturePaper(
            paper_id=f"pmid:{pmid}",
            pmid=pmid,
            source_databases=["PubMed"],
            raw_evidence_ids=[raw_evidence_id],
        ))
    return papers


def _merge_papers(
    first: RetrievedLiteraturePaper,
    second: RetrievedLiteraturePaper,
) -> RetrievedLiteraturePaper:
    def longer(left: str, right: str) -> str:
        return right if len(right) > len(left) else left

    title = longer(first.title, second.title)
    abstract = longer(first.abstract, second.abstract)
    snippets = longer(first.full_text_snippets, second.full_text_snippets)
    access_order = {"metadata_only": 0, "abstract_only": 1, "full_text_snippets": 2}
    access_level = max(
        (first.access_level, second.access_level),
        key=lambda item: access_order[item],
    )
    retraction_order = {
        "not_checked": 0,
        "uncertain": 1,
        "not_retracted": 2,
        "expression_of_concern": 3,
        "retracted": 4,
    }
    retraction_status = max(
        (first.retraction_status, second.retraction_status),
        key=lambda item: retraction_order[item],
    )
    source_text = "\n".join([title, abstract, snippets])
    return first.model_copy(update={
        "title": title,
        "abstract": abstract,
        "full_text_snippets": snippets,
        "authors": _unique([*first.authors, *second.authors]),
        "journal": longer(first.journal, second.journal),
        "year": first.year or second.year,
        "doi": first.doi or second.doi,
        "pmid": first.pmid or second.pmid,
        "pmcid": first.pmcid or second.pmcid,
        "source_databases": _unique([
            *first.source_databases,
            *second.source_databases,
        ]),
        "source_urls": _unique([*first.source_urls, *second.source_urls]),
        "access_level": access_level,
        "retraction_status": retraction_status,
        "raw_evidence_ids": _unique([
            *first.raw_evidence_ids,
            *second.raw_evidence_ids,
        ]),
        "source_text_sha256": hashlib.sha256(source_text.encode("utf-8")).hexdigest(),
    })


def _discovery_text(value: str) -> str:
    """Normalize discovery text while deliberately ignoring stereochemistry."""

    text = str(value or "").translate(str.maketrans({
        "−": "-", "–": "-", "—": "-", "＋": "+",
    }))
    text = re.sub(
        r"\(\s*[+\-RSDL]\s*\)\s*-\s*",
        "",
        text,
        flags=re.IGNORECASE,
    )
    return " ".join(re.sub(r"[^\w]+", " ", text.casefold()).split())


def _relevance_score(
    paper: RetrievedLiteraturePaper,
    requirement: LiteratureActivityRequirement,
) -> tuple[int, str]:
    terms = _unique([
        *(item.name for item in requirement.substrates if item.name),
        *(item.name for item in requirement.products if item.name),
        requirement.reaction_name,
        requirement.reaction_id,
    ])
    normalized_terms = _unique([_discovery_text(term) for term in terms])
    title = _discovery_text(paper.title)
    body = _discovery_text(paper.source_text)
    score = sum(5 for term in normalized_terms if term and term in title)
    score += sum(2 for term in normalized_terms if term and term in body)
    substrate_terms = _unique([
        _discovery_text(item.name or item.compound_id)
        for item in requirement.substrates if item.core
    ])
    product_terms = _unique([
        _discovery_text(item.name or item.compound_id)
        for item in requirement.products if item.core
    ])
    substrate_match = any(term and term in body for term in substrate_terms)
    product_match = any(term and term in body for term in product_terms)
    if substrate_match and product_match:
        score += 8
    if re.search(r"\boverexpress\w*\b", body):
        score += 6
    if re.search(r"\b(?:dephosphorylat\w*|bioconversion|enzyme activity)\b", body):
        score += 5
    if "escherichia coli" in body or re.search(r"\be coli\b", body):
        score += 3
    if re.search(r"\b[A-Z][a-z]{1,8}[A-Z0-9]\b", paper.source_text):
        score += 3
    score += min(8, 2 * len(re.findall(
        r"\b(?:enzyme|activity|cataly\w*|knockout|complement\w*|purif\w*)\b",
        body,
        re.IGNORECASE,
    )))
    return score, paper.paper_id


def _source_record_ids(content: str) -> list[str]:
    values = [
        *(f"PMID:{value}" for value in _unique(re.findall(r"\bPMID\s*:?\s*(\d{5,9})\b", content, re.I))),
        *(f"PMC:{value.upper()}" for value in _unique(_PMCID_PATTERN.findall(content))),
        *(f"DOI:{_normalise_doi(value)}" for value in _unique(_DOI_PATTERN.findall(content))),
    ]
    return _unique([value for value in values if not value.endswith(":")])


class ToolUniverseLiteratureRetriever:
    """Retrieve a small evidence set from a strict ToolUniverse allowlist."""

    def __init__(
        self,
        *,
        tools: Sequence[Any] | None = None,
        timeout_seconds: float = 45.0,
        max_attempts: int = 2,
    ) -> None:
        self._provided_tools = list(tools) if tools is not None else None
        self._tools: dict[str, Any] | None = None
        self._tools_lock = asyncio.Lock()
        self._timeout_seconds = float(timeout_seconds)
        self._max_attempts = int(max_attempts)
        if self._timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive")
        if self._max_attempts < 1:
            raise ValueError("max_attempts must be positive")

    async def _ensure_tools(self) -> dict[str, Any]:
        if self._tools is not None:
            return self._tools
        async with self._tools_lock:
            if self._tools is not None:
                return self._tools
            if self._provided_tools is None:
                # Deliberately lazy: importing this class or executing disabled
                # mode cannot start ToolUniverse.
                from src.protein_selection.agents.literature_researcher import (
                    load_literature_researcher_tools,
                )

                available = await load_literature_researcher_tools()
            else:
                available = self._provided_tools
            self._tools = {tool.name: tool for tool in available}
            return self._tools

    async def _call(
        self,
        requirement: LiteratureActivityRequirement,
        *,
        query_id: str,
        source: str,
        operation: str,
        arguments: dict[str, Any],
        query: str = "",
    ) -> _ToolCallResult:
        tools = await self._ensure_tools()
        tool = tools.get(operation)
        if tool is None:
            message = f"ToolUniverse did not expose required tool {operation}"
            return self._failed_call(
                requirement,
                query_id=query_id,
                source=source,
                operation=operation,
                query=query,
                message=message,
                retryable=False,
            )
        for attempt in range(1, self._max_attempts + 1):
            try:
                async with asyncio.timeout(self._timeout_seconds):
                    raw = await tool.ainvoke(arguments)
                content = _serialize_tool_result(raw)
                tool_error = _tool_result_error(raw, content)
                if tool_error:
                    message = tool_error
                    retryable = bool(_RETRYABLE_PATTERN.search(message))
                    if retryable and attempt < self._max_attempts:
                        await asyncio.sleep(0.5 * attempt)
                        continue
                    return self._failed_call(
                        requirement,
                        query_id=query_id,
                        source=source,
                        operation=operation,
                        query=query,
                        message=message,
                        retryable=retryable,
                        timeout=bool(re.match(r"^\s*(?:timeout|timed\s+out)\b", message, re.I)),
                    )
                digest = hashlib.sha256(content.encode("utf-8")).hexdigest()
                raw_id = f"raw_{operation}_{digest[:16]}"
                papers = _papers_from_content(
                    content,
                    source=source,
                    raw_evidence_id=raw_id,
                )
                audit = LiteratureQueryAudit(
                    query_id=query_id,
                    step_index=requirement.step_index,
                    reaction_id=requirement.reaction_id,
                    source=source,
                    operation=operation,
                    query=query,
                    status="success" if (content or papers) else "not_found",
                    result_count=len(papers),
                    response_sha256=digest,
                    source_record_ids=_source_record_ids(content),
                )
                return _ToolCallResult(audit=audit, content=content, raw_evidence_id=raw_id)
            except TimeoutError:
                message = f"timed out after {self._timeout_seconds:g} seconds"
                retryable = True
            except Exception as exc:  # external tools must fail closed into data
                message = f"{type(exc).__name__}: {exc}"
                retryable = bool(_RETRYABLE_PATTERN.search(message))
            if retryable and attempt < self._max_attempts:
                await asyncio.sleep(0.5 * attempt)
                continue
            return self._failed_call(
                requirement,
                query_id=query_id,
                source=source,
                operation=operation,
                query=query,
                message=message,
                retryable=retryable,
                timeout=message.startswith("timed out"),
            )
        raise AssertionError("unreachable retry loop")

    @staticmethod
    def _failed_call(
        requirement: LiteratureActivityRequirement,
        *,
        query_id: str,
        source: str,
        operation: str,
        query: str,
        message: str,
        retryable: bool,
        timeout: bool = False,
    ) -> _ToolCallResult:
        digest = hashlib.sha256(
            f"{query_id}|{operation}|{message}".encode("utf-8")
        ).hexdigest()
        failure = LiteratureActivityFailure(
            failure_id=f"failure_{digest[:16]}",
            step_index=requirement.step_index,
            reaction_id=requirement.reaction_id,
            stage="retrieval",
            source=source,
            query_id=query_id,
            message=message,
            retryable=retryable,
        )
        audit = LiteratureQueryAudit(
            query_id=query_id,
            step_index=requirement.step_index,
            reaction_id=requirement.reaction_id,
            source=source,
            operation=operation,
            query=query,
            status="timeout" if timeout else "error",
            error=message,
        )
        return _ToolCallResult(
            audit=audit,
            content="",
            raw_evidence_id=f"raw_error_{digest[:16]}",
            failure=failure,
        )

    async def retrieve(
        self,
        requirement: LiteratureActivityRequirement,
        queries: list[LiteratureSearchQuery],
        *,
        max_papers: int,
        max_full_texts: int,
    ) -> LiteratureRetrievalBatch:
        if max_papers < 1 or max_full_texts < 0:
            raise ValueError("invalid literature retrieval limits")
        # Initialize once before concurrent search calls; otherwise every first
        # query could race to start its own stdio server process.
        await self._ensure_tools()
        audits: list[LiteratureQueryAudit] = []
        failures: list[LiteratureActivityFailure] = []
        papers_by_key: dict[str, RetrievedLiteraturePaper] = {}

        async def run_search(item: LiteratureSearchQuery) -> _ToolCallResult:
            if item.source == "PubMed":
                operation = "PubMed_search_articles"
                arguments = {
                    "query": item.query,
                    "limit": min(5, max_papers),
                    "include_abstract": True,
                }
            else:
                operation = "EuropePMC_search_articles"
                arguments = {
                    "query": item.query,
                    "limit": min(5, max_papers),
                    "require_has_ft": False,
                    "enrich_missing_abstract": True,
                }
            return await self._call(
                requirement,
                query_id=item.query_id,
                source=item.source,
                operation=operation,
                arguments=arguments,
                query=item.query,
            )

        search_results = await asyncio.gather(*(run_search(item) for item in queries))
        for result in search_results:
            audits.append(result.audit)
            if result.failure is not None:
                failures.append(result.failure)
                continue
            for paper in _papers_from_content(
                result.content,
                source=result.audit.source,
                raw_evidence_id=result.raw_evidence_id,
            ):
                existing = papers_by_key.get(paper.paper_id)
                papers_by_key[paper.paper_id] = (
                    _merge_papers(existing, paper) if existing else paper
                )

        ranked = sorted(
            papers_by_key.values(),
            key=lambda paper: (-_relevance_score(paper, requirement)[0], paper.paper_id),
        )[:max_papers]

        detail_targets = _unique([paper.pmid for paper in ranked if paper.pmid])
        detail_results = await asyncio.gather(*(
            self._call(
                requirement,
                query_id=f"lit_s{requirement.step_index}_pmid_{pmid}",
                source="PubMed",
                operation="PubMed_get_article",
                arguments={"pmid": pmid},
                query=f"PMID:{pmid}",
            )
            for pmid in detail_targets[:max_papers]
        )) if detail_targets else []
        for result in detail_results:
            audits.append(result.audit)
            if result.failure is not None:
                failures.append(result.failure)
                continue
            for paper in _papers_from_content(
                result.content,
                source="PubMed",
                raw_evidence_id=result.raw_evidence_id,
            ):
                existing = papers_by_key.get(paper.paper_id)
                if existing is None:
                    # Cross-identifier merge when search and detail use DOI vs PMID.
                    existing_key = next((
                        key for key, value in papers_by_key.items()
                        if (paper.pmid and value.pmid == paper.pmid)
                        or (paper.doi and value.doi == paper.doi)
                    ), None)
                    if existing_key is not None:
                        existing = papers_by_key.pop(existing_key)
                papers_by_key[paper.paper_id] = (
                    _merge_papers(existing, paper) if existing else paper
                )

        ranked = sorted(
            papers_by_key.values(),
            key=lambda paper: (-_relevance_score(paper, requirement)[0], paper.paper_id),
        )[:max_papers]
        fulltext_targets = [paper for paper in ranked if paper.pmcid][:max_full_texts]
        terms = _unique([
            *(item.name for item in requirement.substrates if item.name),
            *(item.name for item in requirement.products if item.name),
            "enzyme activity",
            "overexpression",
            "purified",
        ])[:12]
        fulltext_results = await asyncio.gather(*(
            self._call(
                requirement,
                query_id=f"lit_s{requirement.step_index}_fulltext_{paper.pmcid}",
                source="Europe PMC",
                operation="EuropePMC_get_fulltext_snippets",
                arguments={
                    "pmcid": paper.pmcid,
                    "terms": terms,
                    "window_chars": 1200,
                    "max_snippets_per_term": 3,
                    "max_total_chars": 20_000,
                },
                query=f"PMCID:{paper.pmcid}",
            )
            for paper in fulltext_targets
        )) if fulltext_targets else []
        for paper, result in zip(fulltext_targets, fulltext_results):
            audits.append(result.audit)
            if result.failure is not None:
                failures.append(result.failure)
                continue
            snippets = result.content[:20_000]
            if snippets:
                updated = paper.model_copy(update={
                    "full_text_snippets": snippets,
                    "access_level": "full_text_snippets",
                    "source_databases": _unique([*paper.source_databases, "Europe PMC"]),
                    "raw_evidence_ids": _unique([
                        *paper.raw_evidence_ids,
                        result.raw_evidence_id,
                    ]),
                    "source_text_sha256": hashlib.sha256(
                        (paper.title + "\n" + paper.abstract + "\n" + snippets).encode("utf-8")
                    ).hexdigest(),
                })
                papers_by_key[paper.paper_id] = updated

        ranked = sorted(
            papers_by_key.values(),
            key=lambda paper: (-_relevance_score(paper, requirement)[0], paper.paper_id),
        )[:max_papers]
        # Crossref is metadata validation only, bounded to the ten most relevant
        # papers. It never creates a biological claim.
        for paper in ranked[:10]:
            if not paper.doi:
                continue
            metadata_result, retract_result = await asyncio.gather(
                self._call(
                    requirement,
                    query_id=f"lit_s{requirement.step_index}_crossref_{hashlib.sha256(paper.doi.encode()).hexdigest()[:10]}",
                    source="Crossref",
                    operation="Crossref_get_work",
                    arguments={"doi": paper.doi},
                    query=f"DOI:{paper.doi}",
                ),
                self._call(
                    requirement,
                    query_id=f"lit_s{requirement.step_index}_retraction_{hashlib.sha256(paper.doi.encode()).hexdigest()[:10]}",
                    source="Crossref",
                    operation="Crossref_check_retraction",
                    arguments={"doi": paper.doi},
                    query=f"DOI:{paper.doi}",
                ),
            )
            for result in (metadata_result, retract_result):
                audits.append(result.audit)
                if result.failure is not None:
                    failures.append(result.failure)
            status = _retraction_status(retract_result.content)
            key = next((
                candidate_key for candidate_key, value in papers_by_key.items()
                if value.doi == paper.doi
            ), paper.paper_id)
            current = papers_by_key.get(key, paper)
            crossref_has_content = any(
                result.failure is None and bool(result.content.strip())
                for result in (metadata_result, retract_result)
            )
            papers_by_key[key] = current.model_copy(update={
                "source_databases": _unique([
                    *current.source_databases,
                    *(["Crossref"] if crossref_has_content else []),
                ]),
                "retraction_status": status,
                "raw_evidence_ids": _unique([
                    *current.raw_evidence_ids,
                    *(
                        [metadata_result.raw_evidence_id]
                        if metadata_result.content else []
                    ),
                    *(
                        [retract_result.raw_evidence_id]
                        if retract_result.content else []
                    ),
                ]),
            })

        ranked = sorted(
            papers_by_key.values(),
            key=lambda paper: (-_relevance_score(paper, requirement)[0], paper.paper_id),
        )[:max_papers]
        return LiteratureRetrievalBatch(
            papers=[paper for paper in ranked if paper.source_text.strip()],
            queries=audits,
            failures=failures,
        )


def _boolean_value(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().casefold()
    if text in {"1", "true", "yes"}:
        return True
    if text in {"0", "false", "no"}:
        return False
    return None


def _structured_retraction_flags(value: Any) -> tuple[bool, bool, bool]:
    """Return concern/retracted/explicit-not-retracted flags recursively."""

    concern = False
    retracted = False
    not_retracted = False
    for mapping in _walk_mappings(value):
        for raw_key, raw_value in mapping.items():
            key = re.sub(r"[^a-z0-9]+", "_", str(raw_key).casefold()).strip("_")
            boolean = _boolean_value(raw_value)
            scalar = (
                str(raw_value or "").strip().casefold()
                if not isinstance(raw_value, (Mapping, list, tuple))
                else ""
            )
            normalized_scalar = re.sub(r"[_\-]+", " ", scalar)
            if (
                "expression" in key and "concern" in key and boolean is True
            ) or (
                key in {"status", "retraction_status", "publication_status", "relation_type"}
                and "expression of concern" in normalized_scalar
            ):
                concern = True
            if (
                ("retracted" in key or key in {"is_retraction", "has_retraction"})
                and boolean is True
            ) or (
                key in {"status", "retraction_status", "publication_status", "relation_type"}
                and normalized_scalar.strip() == "retracted"
            ):
                retracted = True
            if (
                ("retracted" in key or key in {"is_retraction", "has_retraction"})
                and boolean is False
            ) or (
                key in {"status", "retraction_status", "publication_status"}
                and normalized_scalar.strip() in {"not retracted", "notretracted"}
            ):
                not_retracted = True
    return concern, retracted, not_retracted


def _retraction_status(content: str) -> str:
    text = str(content or "").strip()
    if not text:
        return "not_checked"
    decoded = _decode_json(text)
    if decoded is not None:
        concern, retracted, not_retracted = _structured_retraction_flags(decoded)
        # An expression of concern is material even when the same envelope also
        # reports is_retracted=false (or, defensively, contradictory true).
        if concern:
            return "expression_of_concern"
        if retracted:
            return "retracted"
        if not_retracted:
            return "not_retracted"

    normalized = re.sub(r"[_\-]+", " ", text.casefold())
    if "expression of concern" in normalized:
        return "expression_of_concern"
    if re.match(r"^\s*retracted\b", normalized) or re.search(
        r'"(?:is )?retracted"\s*:\s*true|"status"\s*:\s*"retracted"',
        normalized,
    ):
        return "retracted"
    if "not retracted" in normalized or re.search(
        r'"(?:is )?retracted"\s*:\s*false',
        normalized,
    ):
        return "not_retracted"
    return "uncertain"


__all__ = ["RETRIEVER_VERSION", "ToolUniverseLiteratureRetriever"]
