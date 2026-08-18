"""Deterministic retrieval plans and an auditable raw-evidence ledger."""

from __future__ import annotations

import asyncio
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import hashlib
import json
import re
from typing import Any, Literal

from langchain_core.tools import BaseTool
from pydantic import BaseModel, ConfigDict, Field

from src.protein_selection.agents.bio_database_researcher import (
    BioDatabaseResearchResult,
)
from src.protein_selection.agents.research_policy import ResearchMode, ResearchPolicy
from src.protein_selection.research_context import ResearchContext
from src.protein_selection.services.cache import (
    PersistentTTLCache,
    RetrievalStats,
    ttl_for_tool,
)


EvidenceStage = Literal[
    "validated_input",
    "database",
    "literature",
    "web",
    "host_compatibility",
]
RetrievalStatus = Literal["success", "error", "timeout", "skipped"]

_PMID_PATTERNS = (
    re.compile(r'"pmid"\s*:\s*"?(\d{5,9})"?', re.IGNORECASE),
    re.compile(r'"uid"\s*:\s*"?(\d{5,9})"?', re.IGNORECASE),
    re.compile(
        r'"database"\s*:\s*"PubMed"\s*,\s*"id"\s*:\s*"?(\d{5,9})"?',
        re.IGNORECASE,
    ),
    re.compile(
        r'"source"\s*:\s*"PubMed"\s*,\s*"id"\s*:\s*"?(\d{5,9})"?',
        re.IGNORECASE,
    ),
    re.compile(r"pubmed(?:\.ncbi\.nlm\.nih\.gov)?/(\d{5,9})", re.IGNORECASE),
)
_URL_PATTERN = re.compile(r"https?://[^\s\"'<>]+")
_UNIPROT_ACCESSION_PATTERN = re.compile(
    r"(?<![A-Z0-9])(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?![A-Z0-9])",
    re.IGNORECASE,
)
_KO_PATTERN = re.compile(r"(?<![A-Z0-9])(?:ko:)?(K\d{5})(?![A-Z0-9])", re.I)
_KEGG_GENE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9_])([a-z][a-z0-9_]{1,9}:[A-Za-z0-9_.-]+)"
)
_REFERENCE_PMID_PATTERN = re.compile(
    r"(?:PMID\s*[:#-]?\s*|pubmed(?:\.ncbi\.nlm\.nih\.gov)?/)"
    r"(\d{5,9})",
    re.IGNORECASE,
)
_GENE_TOKEN_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:[a-z]{2,10}[A-Z]\d?|"
    r"[A-Z][a-z]{1,9}[A-Z]\d?)(?![A-Za-z0-9])"
)

_CURATED_CANDIDATE_TOOLS = frozenset(
    {
        "validated_input",
        "UniProt_validated_input",
        "KEGG_validated_input",
        "ComplexPortal_get_complex",
        "ComplexPortal_search_complexes",
        "UniProt_get_entry_by_accession",
        "UniProt_get_function_by_accession",
        "KEGG_get_reaction",
        "KEGG_get_enzyme",
        "Rhea_get_reaction",
        "Rhea_get_reaction_participants",
    }
)
_CANDIDATE_SEARCH_FIELDS = [
    "accession",
    "id",
    "gene_names",
    "protein_name",
    "organism_name",
    "organism_id",
    "reviewed",
]


class RawResearchEvidence(BaseModel):
    """One raw tool response stored before any model interpretation."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    evidence_id: str = Field(min_length=1)
    stage: EvidenceStage
    tool_name: str = Field(min_length=1)
    query_arguments: dict[str, Any]
    status: RetrievalStatus
    content: str = ""
    content_sha256: str | None = None
    source_record_ids: list[str] = Field(default_factory=list)
    source_urls: list[str] = Field(default_factory=list)
    error: str | None = None
    truncated: bool = False
    original_char_count: int = 0
    cache_hit: bool = False
    attempt_count: int = Field(default=1, ge=0)


@dataclass(frozen=True, slots=True)
class PlannedToolCall:
    """One code-generated tool invocation."""

    stage: EvidenceStage
    tool_name: str
    arguments: dict[str, Any]
    rationale: str


@dataclass(frozen=True, slots=True)
class PlannedLiteratureQuery:
    """One bounded source-specific literature query."""

    source: Literal["PubMed", "Europe PMC", "Semantic Scholar"]
    query: str


CandidateSeedKind = Literal[
    "accession",
    "gene",
    "locus",
    "protein_name",
    "ko",
]
CandidateResolutionStatus = Literal[
    "verified",
    "ambiguous",
    "not_found",
    "verified_empty",
    "source_error",
]


@dataclass(frozen=True, slots=True)
class CandidateSeed:
    """One curated, provenance-preserving clue to a partner identity.

    A seed is deliberately weaker than a resolved UniProt identity.  In
    particular, gene symbols and protein names are only searched within the
    exact source taxon, while a KEGG orthology seed must retain the source
    organism code needed for the KO -> gene -> UniProt conversion ladder.
    """

    kind: CandidateSeedKind
    value: str
    source_evidence_id: str
    source_span: str
    role_hint: str | None = None
    taxon_id: int | None = None
    kegg_organism_code: str | None = None


@dataclass(frozen=True, slots=True)
class VerifiedCandidateHint:
    """A partner identity verified by an exact UniProt lookup.

    The hint is used only to construct a literature query and to resolve an
    identifier after a paper has independently established the biological
    role.  It is not itself evidence that the partner is required.
    """

    uniprot_id: str
    protein_name: str | None = None
    gene_names: tuple[str, ...] = ()
    organism_name: str | None = None
    taxon_id: int | None = None
    aliases: tuple[str, ...] = ()
    ko_ids: tuple[str, ...] = ()
    seed_evidence_ids: tuple[str, ...] = ()
    reference_pmids: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class CandidateResolution:
    """Auditable outcome for one candidate seed.

    ``verified_empty`` means the authoritative source answered successfully
    and returned no matching records.  It must never be synthesized from a
    timeout or transport/server failure; those are ``source_error``.
    """

    status: CandidateResolutionStatus
    seed: CandidateSeed
    candidates: tuple[VerifiedCandidateHint, ...]
    reason: str


class ResearchEvidenceLedger:
    """In-memory append-only ledger for one workflow invocation."""

    def __init__(self) -> None:
        self._records: list[RawResearchEvidence] = []
        self._by_cache_key: dict[str, RawResearchEvidence] = {}
        self._inflight: dict[str, asyncio.Future[RawResearchEvidence]] = {}
        self._lock = asyncio.Lock()

    @property
    def records(self) -> tuple[RawResearchEvidence, ...]:
        return tuple(self._records)

    def records_for(self, stage: EvidenceStage) -> list[RawResearchEvidence]:
        return [record for record in self._records if record.stage == stage]

    async def add_static(self, record: RawResearchEvidence) -> None:
        """Append a deterministic non-tool source such as validated input."""

        cache_key = f"static:{record.evidence_id}"
        async with self._lock:
            if cache_key in self._by_cache_key:
                return
            self._records.append(record)
            self._by_cache_key[cache_key] = record

    async def cached(self, cache_key: str) -> RawResearchEvidence | None:
        async with self._lock:
            return self._by_cache_key.get(cache_key)

    async def claim(
        self,
        cache_key: str,
    ) -> tuple[asyncio.Future[RawResearchEvidence], bool]:
        async with self._lock:
            completed = self._by_cache_key.get(cache_key)
            if completed is not None:
                future = asyncio.get_running_loop().create_future()
                future.set_result(completed)
                return future, False
            existing = self._inflight.get(cache_key)
            if existing is not None:
                return existing, False
            future = asyncio.get_running_loop().create_future()
            self._inflight[cache_key] = future
            return future, True

    async def complete(
        self,
        cache_key: str,
        record: RawResearchEvidence,
    ) -> None:
        async with self._lock:
            self._records.append(record)
            self._by_cache_key[cache_key] = record
            future = self._inflight.pop(cache_key)
            if not future.done():
                future.set_result(record)


class ResearchToolRunner:
    """Invoke only named tools from a deterministic plan and retain results."""

    def __init__(
        self,
        *,
        tools: Sequence[BaseTool],
        ledger: ResearchEvidenceLedger,
        policy: ResearchPolicy,
        stage_to_role: Mapping[EvidenceStage, str],
        max_content_chars: int = 24_000,
        response_cache: PersistentTTLCache | None = None,
        stats: RetrievalStats | None = None,
        max_attempts: int = 3,
    ) -> None:
        self._tools = {tool.name: tool for tool in tools}
        self._ledger = ledger
        self._policy = policy
        self._stage_to_role = dict(stage_to_role)
        self._max_content_chars = max_content_chars
        self._response_cache = response_cache
        self._stats = stats
        self._max_attempts = max_attempts
        self._global_semaphore = asyncio.Semaphore(4)
        self._source_semaphores = {
            "pubmed": asyncio.Semaphore(2),
            "europe_pmc": asyncio.Semaphore(2),
            "database": asyncio.Semaphore(2),
            "semantic_scholar": asyncio.Semaphore(1),
            "web": asyncio.Semaphore(1),
        }

    async def run(self, call: PlannedToolCall) -> RawResearchEvidence:
        cache_key = _cache_key(call.tool_name, call.arguments)
        future, owns_call = await self._ledger.claim(cache_key)
        if not owns_call:
            return await future

        tool = self._tools.get(call.tool_name)
        if tool is None:
            record = self._failure_record(
                call,
                "error",
                f"planned tool {call.tool_name} is unavailable",
            )
            await self._ledger.complete(cache_key, record)
            return record

        cached = self._cached_record(call, cache_key)
        if cached is not None:
            await self._ledger.complete(cache_key, cached)
            if self._stats is not None:
                self._stats.cache_hits += 1
            return cached

        timeout_seconds = self._timeout_for(call)
        record: RawResearchEvidence | None = None
        for attempt in range(1, self._max_attempts + 1):
            if self._stats is not None and attempt == 1:
                self._stats.live_calls += 1
            try:
                raw = await self._invoke_limited(
                    tool,
                    call,
                    timeout_seconds,
                )
                content = _serialize_tool_result(raw)
                tool_error = _tool_result_error(
                    raw,
                    content,
                    tool_name=call.tool_name,
                )
                if tool_error is not None:
                    if (
                        attempt < self._max_attempts
                        and _is_transient_tool_error(tool_error)
                    ):
                        await self._retry_delay(attempt)
                        continue
                    record = self._failure_record(
                        call,
                        "error",
                        tool_error,
                        attempt_count=attempt,
                    )
                    break
                content = _compact_known_tool_content(
                    call.tool_name,
                    content,
                )
                original_char_count = len(content)
                truncated = original_char_count > self._max_content_chars
                bounded_content = content[: self._max_content_chars]
                digest = hashlib.sha256(content.encode("utf-8")).hexdigest()
                record = RawResearchEvidence(
                    evidence_id=_evidence_id(call, digest),
                    stage=call.stage,
                    tool_name=call.tool_name,
                    query_arguments=call.arguments,
                    status="success",
                    content=bounded_content,
                    content_sha256=digest,
                    source_record_ids=_extract_source_record_ids(
                        call.tool_name,
                        call.arguments,
                        content,
                    ),
                    source_urls=_extract_urls(content),
                    truncated=truncated,
                    original_char_count=original_char_count,
                    attempt_count=attempt,
                )
                self._cache_record(call, cache_key, record)
                break
            except TimeoutError:
                if attempt < self._max_attempts:
                    await self._retry_delay(attempt)
                    continue
                record = self._failure_record(
                    call,
                    "timeout",
                    f"timed out after {timeout_seconds:g} seconds",
                    attempt_count=attempt,
                )
            except Exception as exc:
                if (
                    attempt < self._max_attempts
                    and _is_transient_tool_error(str(exc))
                ):
                    await self._retry_delay(attempt)
                    continue
                record = self._failure_record(
                    call,
                    "error",
                    str(exc),
                    attempt_count=attempt,
                )
                break

        if record is None:
            record = self._failure_record(
                call,
                "error",
                "tool retry loop ended unexpectedly",
                attempt_count=self._max_attempts,
            )
        if record.status in {"error", "timeout"} and self._stats is not None:
            self._stats.failed_calls += 1

        await self._ledger.complete(cache_key, record)
        return record

    async def _invoke_limited(
        self,
        tool: BaseTool,
        call: PlannedToolCall,
        timeout_seconds: float,
    ) -> Any:
        source = _tool_source(call.tool_name)
        semaphore = self._source_semaphores[source]
        async with self._global_semaphore, semaphore:
            async with asyncio.timeout(timeout_seconds):
                return await tool.ainvoke(call.arguments)

    async def _retry_delay(self, attempt: int) -> None:
        if self._stats is not None:
            self._stats.retries += 1
        await asyncio.sleep(0.5 * (2 ** (attempt - 1)))

    def _cached_record(
        self,
        call: PlannedToolCall,
        cache_key: str,
    ) -> RawResearchEvidence | None:
        if self._response_cache is None:
            return None
        cached = self._response_cache.get("research.tool", cache_key)
        if not isinstance(cached, Mapping):
            return None
        content = cached.get("content")
        digest = cached.get("content_sha256")
        if not isinstance(content, str) or not isinstance(digest, str):
            return None
        return RawResearchEvidence(
            evidence_id=_evidence_id(call, digest),
            stage=call.stage,
            tool_name=call.tool_name,
            query_arguments=call.arguments,
            status="success",
            content=content,
            content_sha256=digest,
            source_record_ids=[
                value
                for value in cached.get("source_record_ids", [])
                if isinstance(value, str)
            ],
            source_urls=[
                value
                for value in cached.get("source_urls", [])
                if isinstance(value, str)
            ],
            truncated=bool(cached.get("truncated", False)),
            original_char_count=int(
                cached.get("original_char_count", len(content))
            ),
            cache_hit=True,
            attempt_count=0,
        )

    def _cache_record(
        self,
        call: PlannedToolCall,
        cache_key: str,
        record: RawResearchEvidence,
    ) -> None:
        if self._response_cache is None or record.status != "success":
            return
        self._response_cache.set(
            "research.tool",
            cache_key,
            {
                "content": record.content,
                "content_sha256": record.content_sha256,
                "source_record_ids": record.source_record_ids,
                "source_urls": record.source_urls,
                "truncated": record.truncated,
                "original_char_count": record.original_char_count,
            },
            ttl_seconds=ttl_for_tool(call.tool_name),
        )

    async def run_many(
        self,
        calls: Sequence[PlannedToolCall],
    ) -> list[RawResearchEvidence]:
        if not calls:
            return []
        return list(await asyncio.gather(*(self.run(call) for call in calls)))

    def _timeout_for(self, call: PlannedToolCall) -> float:
        role = self._stage_to_role.get(call.stage)
        if role is None:
            return self._policy.model_timeout_seconds
        budget = self._policy.budget_for(role)  # type: ignore[arg-type]
        if budget is None:
            return 45.0
        return budget.per_tool_timeouts.get(
            call.tool_name,
            budget.default_tool_timeout_seconds,
        )

    @staticmethod
    def _failure_record(
        call: PlannedToolCall,
        status: Literal["error", "timeout"],
        error: str,
        *,
        attempt_count: int = 1,
    ) -> RawResearchEvidence:
        digest = hashlib.sha256(error.encode("utf-8")).hexdigest()
        return RawResearchEvidence(
            evidence_id=_evidence_id(call, digest),
            stage=call.stage,
            tool_name=call.tool_name,
            query_arguments=call.arguments,
            status=status,
            error=error,
            attempt_count=attempt_count,
        )


def build_database_call_plan(
    context: ResearchContext,
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Build exact database calls from official identifiers and annotations."""

    calls: list[PlannedToolCall] = []
    protein = context.protein
    reaction = context.reaction
    complex_limit = 2 if mode == "balanced" else 5

    if context.preliminary_reaction_match != "matched":
        for rhea_id in reaction.rhea_ids[:1 if mode == "balanced" else 2]:
            calls.append(
                PlannedToolCall(
                    stage="database",
                    tool_name="Rhea_get_reaction",
                    arguments={"rhea_id": rhea_id},
                    rationale="Resolve the exact KEGG reaction in Rhea.",
                )
            )
        for ec_number in reaction.ec_numbers[:1]:
            calls.append(
                PlannedToolCall(
                    stage="database",
                    tool_name="Rhea_search_by_ec",
                    arguments={
                        "ec_number": ec_number,
                        "limit": 10 if mode == "balanced" else 20,
                    },
                    rationale="Enumerate reaction ambiguity for the shared EC.",
                )
            )

    complex_ids = [
        crossref.record_id
        for crossref in protein.cross_references
        if crossref.database == "ComplexPortal"
    ]
    for complex_id in complex_ids[:complex_limit]:
        calls.append(
            PlannedToolCall(
                stage="database",
                tool_name="ComplexPortal_get_complex",
                arguments={"complex_id": complex_id},
                rationale="Read a complex explicitly linked by UniProt.",
            )
        )

    if not complex_ids and protein.subunit_annotations:
        query = _first_search_term(protein)
        if query is not None and protein.taxon_id is not None:
            calls.append(
                PlannedToolCall(
                    stage="database",
                    tool_name="ComplexPortal_search_complexes",
                    arguments={
                        "query": query,
                        "species": str(protein.taxon_id),
                        "number": 5 if mode == "balanced" else 10,
                    },
                    rationale="Resolve an annotated heteromer using exact taxonomy.",
                )
            )

    has_interaction_clue = any(
        crossref.database == "IntAct" for crossref in protein.cross_references
    )
    if has_interaction_clue or (mode == "deep" and protein.subunit_annotations):
        calls.append(
            PlannedToolCall(
                stage="database",
                tool_name="intact_get_interactions",
                arguments={
                    "identifier": protein.primary_accession,
                    "format": "json",
                },
                rationale="Retrieve interactions for the exact accession once.",
            )
        )

    return _deduplicate_calls(calls)


def extract_candidate_seeds(
    context: ResearchContext,
    records: Sequence[RawResearchEvidence],
) -> list[CandidateSeed]:
    """Extract only curated partner clues, retaining their exact provenance.

    Free-form web and literature text is intentionally excluded.  Exact
    accessions, structured participant fields, curated subunit annotations,
    and KEGG reaction orthology are useful discovery clues; none of them is
    promoted to a verified dependency identity at this stage.
    """

    protein = context.protein
    reaction = context.reaction
    taxon_id = getattr(protein, "taxon_id", None)
    kegg_organism_code = _context_kegg_organism_code(context)
    input_accession = str(protein.primary_accession).upper()
    excluded = {
        "accession": {input_accession.casefold()},
        "gene": {
            str(value).casefold()
            for value in getattr(protein, "gene_names", ())
            if isinstance(value, str)
        },
        "locus": {
            str(value).casefold()
            for value in getattr(protein, "locus_tags", ())
            if isinstance(value, str)
        },
        "protein_name": {
            _normalize_candidate_text(value)
            for value in getattr(protein, "protein_names", ())
            if isinstance(value, str)
        },
        "ko": set(),
    }
    seeds: list[CandidateSeed] = []
    seen: set[tuple[str, str, str, int | None]] = set()

    def add_seed(
        kind: CandidateSeedKind,
        raw_value: Any,
        evidence_id: str,
        source_span: str,
        role_hint: str | None = None,
        *,
        seed_taxon_id: int | None = taxon_id,
        organism_code: str | None = kegg_organism_code,
    ) -> None:
        value = _normalize_seed_value(kind, raw_value)
        if value is None:
            return
        comparison = (
            _normalize_candidate_text(value)
            if kind == "protein_name"
            else value.casefold()
        )
        if comparison in excluded[kind]:
            return
        key = (kind, comparison, evidence_id, seed_taxon_id)
        if key in seen:
            return
        seen.add(key)
        seeds.append(
            CandidateSeed(
                kind=kind,
                value=value,
                source_evidence_id=evidence_id,
                source_span=_bounded_source_span(source_span),
                role_hint=_bounded_source_span(role_hint)
                if role_hint
                else None,
                taxon_id=seed_taxon_id,
                kegg_organism_code=(
                    organism_code.casefold() if organism_code else None
                ),
            )
        )

    # Future context builders may expose explicit curated candidate fields.
    # getattr keeps this planner compatible with the current strict model.
    explicit_fields: tuple[tuple[CandidateSeedKind, tuple[str, ...]], ...] = (
        (
            "accession",
            (
                "candidate_accessions",
                "curated_candidate_accessions",
                "partner_accessions",
            ),
        ),
        (
            "gene",
            ("candidate_genes", "candidate_gene_names", "partner_genes"),
        ),
        (
            "locus",
            ("candidate_loci", "candidate_locus_tags", "partner_loci"),
        ),
        (
            "protein_name",
            ("candidate_protein_names", "partner_protein_names"),
        ),
        (
            "ko",
            ("candidate_ko_ids", "candidate_kos", "partner_ko_ids"),
        ),
    )
    for owner_name, owner in (
        ("context", context),
        ("context.protein", protein),
        ("context.reaction", reaction),
    ):
        for kind, field_names in explicit_fields:
            for field_name in field_names:
                for value in _iter_candidate_values(
                    getattr(owner, field_name, ())
                ):
                    add_seed(
                        kind,
                        value,
                        "CONTEXT-CURATED",
                        f"{owner_name}.{field_name}: {value}",
                    )

    for item in getattr(reaction, "orthology", ()):
        orthology_id = getattr(item, "orthology_id", None)
        if orthology_id is None and isinstance(item, Mapping):
            orthology_id = item.get("orthology_id") or item.get("id")
        description = getattr(item, "description", None)
        if description is None and isinstance(item, Mapping):
            description = item.get("description")
        add_seed(
            "ko",
            orthology_id,
            "CONTEXT-REACTION",
            f"reaction.orthology: {orthology_id} {description or ''}".strip(),
            str(description) if description else None,
        )

    for annotation in getattr(protein, "subunit_annotations", ()):
        if not isinstance(annotation, str):
            continue
        for match in _UNIPROT_ACCESSION_PATTERN.finditer(annotation):
            add_seed(
                "accession",
                match.group(0),
                "CONTEXT-PROTEIN",
                annotation,
                annotation,
            )
        for match in _KO_PATTERN.finditer(annotation):
            add_seed(
                "ko",
                match.group(1),
                "CONTEXT-PROTEIN",
                annotation,
                annotation,
            )
        for match in _GENE_TOKEN_PATTERN.finditer(annotation):
            gene_token = match.group(0)
            if gene_token[0].isupper():
                gene_token = gene_token[0].lower() + gene_token[1:]
            add_seed(
                "gene",
                gene_token,
                "CONTEXT-PROTEIN",
                annotation,
                annotation,
            )

    for record in records:
        if (
            record.status != "success"
            or record.stage not in {"validated_input", "database"}
            or record.tool_name not in _CURATED_CANDIDATE_TOOLS
        ):
            continue
        payload = _json_value(record.content)
        if payload is not None:
            for kind, value, span, role_hint in _structured_candidate_values(
                payload
            ):
                add_seed(
                    kind,
                    value,
                    record.evidence_id,
                    span,
                    role_hint,
                )
        if record.tool_name.startswith("ComplexPortal_"):
            for match in _UNIPROT_ACCESSION_PATTERN.finditer(record.content):
                add_seed(
                    "accession",
                    match.group(0),
                    record.evidence_id,
                    match.group(0),
                    "curated complex participant",
                )
        for match in _KO_PATTERN.finditer(record.content):
            add_seed(
                "ko",
                match.group(1),
                record.evidence_id,
                match.group(0),
                "curated reaction orthology",
            )

    return seeds


def candidate_resolution_calls(
    seeds: Sequence[CandidateSeed],
    input_uniprot_id: str,
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Plan the first bounded identity-resolution rung for curated seeds."""

    calls: list[PlannedToolCall] = []
    limit = 5 if mode == "balanced" else 10
    for seed in _select_candidate_seeds(seeds, mode):
        value = seed.value.strip()
        if seed.kind == "accession":
            accession = value.upper()
            if (
                accession == input_uniprot_id.upper()
                or not _UNIPROT_ACCESSION_PATTERN.fullmatch(accession)
            ):
                continue
            calls.append(
                PlannedToolCall(
                    stage="database",
                    tool_name="UniProt_get_entry_by_accession",
                    arguments={"accession": accession, "compact": True},
                    rationale=(
                        "Verify an exact accession from curated candidate "
                        f"evidence {seed.source_evidence_id}."
                    ),
                )
            )
            continue

        if seed.kind == "ko":
            if seed.kegg_organism_code is None:
                continue
            calls.append(
                PlannedToolCall(
                    stage="database",
                    tool_name="KEGG_link_entries",
                    arguments={
                        "source": f"ko:{value.upper()}",
                        "target": seed.kegg_organism_code,
                    },
                    rationale=(
                        "Resolve a curated KO to genes in the exact KEGG "
                        "source organism before ID conversion."
                    ),
                )
            )
            continue

        if seed.taxon_id is None:
            # Gene/name searches without an exact source taxon are too prone
            # to silently returning a human or unrelated bacterial protein.
            continue
        escaped = _uniprot_query_value(value)
        field = (
            "protein_name"
            if seed.kind == "protein_name"
            else "gene_exact"
        )
        query = (
            f"{field}:{escaped} AND taxonomy_id:{seed.taxon_id} "
            "AND reviewed:true"
        )
        calls.extend(
            [
                PlannedToolCall(
                    stage="database",
                    tool_name="UniProt_search",
                    arguments={
                        "query": query,
                        "limit": limit,
                        "fields": list(_CANDIDATE_SEARCH_FIELDS),
                    },
                    rationale=(
                        "Search reviewed UniProt entries for the curated "
                        f"{seed.kind} in exact taxon {seed.taxon_id}."
                    ),
                ),
                PlannedToolCall(
                    stage="database",
                    tool_name="proteins_api_search",
                    arguments={
                        "query": value,
                        "taxid": str(seed.taxon_id),
                        "reviewed": True,
                        "size": limit,
                        "format": "json",
                    },
                    rationale=(
                        "Cross-check the same curated identity through EBI "
                        "Proteins with an exact taxon and reviewed filter."
                    ),
                ),
            ]
        )
    return _deduplicate_calls(calls)


def candidate_kegg_conversion_calls(
    link_records: Sequence[RawResearchEvidence],
    seeds: Sequence[CandidateSeed],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Plan KO-linked organism gene -> UniProt conversions.

    This is deliberately a second phase: KEGG conversion accepts organism
    gene identifiers, not KO identifiers.  Failed links yield no conversion
    calls and remain distinguishable from successful empty link responses in
    :func:`resolve_candidate_identities`.
    """

    calls: list[PlannedToolCall] = []
    per_seed_limit = 3 if mode == "balanced" else 6
    ko_seeds = [seed for seed in seeds if seed.kind == "ko"]
    for record in link_records:
        if record.tool_name != "KEGG_link_entries" or record.status != "success":
            continue
        source = str(record.query_arguments.get("source") or "")
        target = str(record.query_arguments.get("target") or "").casefold()
        seed = next(
            (
                item
                for item in ko_seeds
                if source.casefold() == f"ko:{item.value}".casefold()
                and item.kegg_organism_code is not None
                and target == item.kegg_organism_code.casefold()
            ),
            None,
        )
        if seed is None:
            continue
        prefix = f"{seed.kegg_organism_code}:".casefold()
        gene_ids = [
            value
            for value in _extract_kegg_gene_ids(record.content)
            if value.casefold().startswith(prefix)
        ]
        for gene_id in gene_ids[:per_seed_limit]:
            calls.append(
                PlannedToolCall(
                    stage="database",
                    tool_name="KEGG_convert_ids",
                    arguments={
                        "kegg_id": gene_id,
                        "target_db": "uniprot",
                    },
                    rationale=(
                        "Convert one exact-organism KEGG gene linked from "
                        f"{seed.value} to UniProt."
                    ),
                )
            )
    return _deduplicate_calls(calls)


def candidate_resolution_verification_calls(
    resolution_records: Sequence[RawResearchEvidence],
    seeds: Sequence[CandidateSeed],
    input_uniprot_id: str,
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Verify accessions returned by search or KEGG conversion exactly."""

    accessions: list[str] = []
    max_candidates = 6 if mode == "balanced" else 12
    for record in resolution_records:
        if record.status != "success":
            continue
        if record.tool_name in {"UniProt_search", "proteins_api_search"}:
            if not any(
                seed.kind in {"gene", "locus", "protein_name"}
                and _record_matches_seed_search(record, seed)
                for seed in seeds
            ):
                continue
            values = [
                str(item.get("accession") or "")
                for item in _extract_search_candidate_mappings(record.content)
            ]
        elif record.tool_name == "KEGG_convert_ids":
            if not any(seed.kind == "ko" for seed in seeds):
                continue
            values = _extract_converted_accessions(record.content)
        else:
            continue
        for raw_accession in values:
            accession = raw_accession.upper()
            if (
                accession == input_uniprot_id.upper()
                or accession in accessions
                or not _UNIPROT_ACCESSION_PATTERN.fullmatch(accession)
            ):
                continue
            accessions.append(accession)
            if len(accessions) >= max_candidates:
                break
        if len(accessions) >= max_candidates:
            break
    return [
        PlannedToolCall(
            stage="database",
            tool_name="UniProt_get_entry_by_accession",
            arguments={"accession": accession, "compact": True},
            rationale="Verify a search- or KEGG-resolved candidate exactly.",
        )
        for accession in accessions
    ]


def resolve_candidate_identities(
    seeds: Sequence[CandidateSeed],
    records: Sequence[RawResearchEvidence],
    input_uniprot_id: str,
) -> list[CandidateResolution]:
    """Resolve seeds without conflating ambiguity, emptiness, and failure."""

    provisional: list[CandidateResolution] = []
    for seed in seeds:
        related, lineage_accessions = _candidate_seed_lineage(seed, records)
        successful = [record for record in related if record.status == "success"]
        failures = [record for record in related if record.status != "success"]
        not_found_failures = [
            record
            for record in failures
            if _is_explicit_not_found(record.error or record.content)
        ]
        source_failures = [
            record for record in failures if record not in not_found_failures
        ]

        hints: list[VerifiedCandidateHint] = []
        seen_accessions: set[str] = set()
        allowed = {
            value.upper()
            for value in lineage_accessions
            if value.upper() != input_uniprot_id.upper()
        }
        if seed.kind == "accession":
            allowed.add(seed.value.upper())
        for record in records:
            if (
                record.status != "success"
                or record.tool_name != "UniProt_get_entry_by_accession"
            ):
                continue
            hint = _verified_hint_from_exact_record(record, seed, allowed)
            if (
                hint is None
                or hint.uniprot_id.upper() == input_uniprot_id.upper()
                or hint.uniprot_id.upper() in seen_accessions
            ):
                continue
            seen_accessions.add(hint.uniprot_id.upper())
            hints.append(hint)

        unresolved_lineage = allowed - seen_accessions
        if len(allowed) > 1 or len(hints) > 1:
            status: CandidateResolutionStatus = "ambiguous"
            reason = (
                "The curated seed maps to multiple candidate accessions; "
                "ambiguity is retained pending disambiguating evidence."
            )
        elif len(hints) == 1 and not unresolved_lineage:
            status = "verified"
            reason = "One exact UniProt identity matched the seed and source taxon."
        elif source_failures:
            status = "source_error"
            reason = (
                "At least one required identity source failed; no empty-result "
                "claim was inferred from that failure."
            )
        elif not_found_failures and not successful:
            status = "not_found"
            reason = "The exact identifier was explicitly reported as not found."
        elif _has_verified_empty_response(seed, successful):
            status = "verified_empty"
            reason = (
                "The authoritative lookup completed successfully and returned "
                "no matching entries."
            )
        elif successful:
            status = "not_found"
            reason = (
                "Lookup responses completed, but no candidate passed exact "
                "identifier and taxon verification."
            )
        elif failures:
            status = "not_found" if not_found_failures else "source_error"
            reason = (
                "The exact identifier was explicitly reported as not found."
                if status == "not_found"
                else "All relevant identity sources failed."
            )
        else:
            status = "not_found"
            reason = "No identity-resolution response was available for this seed."
        provisional.append(
            CandidateResolution(
                status=status,
                seed=seed,
                candidates=tuple(hints),
                reason=reason,
            )
        )

    merged_hints = _merge_candidate_hints(
        hint
        for resolution in provisional
        for hint in resolution.candidates
    )
    return [
        CandidateResolution(
            status=resolution.status,
            seed=resolution.seed,
            candidates=tuple(
                merged_hints[hint.uniprot_id.upper()]
                for hint in resolution.candidates
            ),
            reason=resolution.reason,
        )
        for resolution in provisional
    ]


def build_literature_query_plan(
    context: ResearchContext,
    database_report: BioDatabaseResearchResult | None,
    mode: ResearchMode,
    candidate_hints: Sequence[VerifiedCandidateHint] = (),
) -> list[PlannedLiteratureQuery]:
    """Build a fair, bounded candidate discovery-query ladder.

    Candidate ladders are interleaved by depth instead of exhausting the
    first candidate.  Balanced mode covers at most three candidates/eight
    plans; deep mode covers at most six candidates/eighteen plans.  Curated
    PMIDs are planned separately by :func:`linked_article_detail_calls` so
    they cannot crowd discovery searches out of this budget.
    """

    protein = context.protein
    primary_gene = next(
        (
            value
            for value in getattr(protein, "gene_names", ())
            if isinstance(value, str) and value.strip()
        ),
        None,
    )
    primary_name = next(
        (
            value
            for value in getattr(protein, "protein_names", ())
            if isinstance(value, str) and value.strip()
        ),
        None,
    )
    identity_terms = [
        f'"{protein.primary_accession}"',
        primary_gene,
        _quote_term(primary_name),
    ]
    identity = " OR ".join(term for term in identity_terms if term)

    report_candidates = (
        getattr(database_report, "candidate_protein_dependencies", ())
        if database_report is not None
        else ()
    )
    candidates: list[dict[str, Any]] = []
    by_accession: dict[str, dict[str, Any]] = {}
    for candidate in report_candidates:
        accession = getattr(candidate, "uniprot_id", None)
        item = {
            "uniprot_id": accession,
            "protein_name": getattr(candidate, "protein_name", None),
            "gene_names": (),
            "aliases": (),
            "ko_ids": (),
            "organism_name": getattr(candidate, "organism_name", None),
            "role": getattr(candidate, "role", None),
            "reference_pmids": (),
        }
        candidates.append(item)
        if isinstance(accession, str):
            by_accession[accession.upper()] = item

    for hint in candidate_hints:
        item = by_accession.get(hint.uniprot_id.upper())
        if item is None:
            item = {
                "uniprot_id": hint.uniprot_id,
                "protein_name": hint.protein_name,
                "gene_names": hint.gene_names,
                "aliases": hint.aliases,
                "ko_ids": hint.ko_ids,
                "organism_name": hint.organism_name,
                "role": hint.protein_name or "candidate protein subunit",
                "reference_pmids": hint.reference_pmids,
            }
            candidates.append(item)
            by_accession[hint.uniprot_id.upper()] = item
        else:
            item["protein_name"] = item.get("protein_name") or hint.protein_name
            item["gene_names"] = _unique_strings(
                (*item.get("gene_names", ()), *hint.gene_names)
            )
            item["aliases"] = _unique_strings(
                (*item.get("aliases", ()), *hint.aliases)
            )
            item["ko_ids"] = _unique_strings(
                (*item.get("ko_ids", ()), *hint.ko_ids)
            )
            item["reference_pmids"] = _unique_strings(
                (*item.get("reference_pmids", ()), *hint.reference_pmids)
            )

    budget = 8 if mode == "balanced" else 18
    max_candidates = 3 if mode == "balanced" else 6
    planned: list[PlannedLiteratureQuery] = []

    reaction_identity = _reaction_literature_identity(context)
    ladders: list[list[PlannedLiteratureQuery]] = []
    for candidate in candidates[:max_candidates]:
        accession = candidate.get("uniprot_id")
        protein_name = candidate.get("protein_name")
        genes = tuple(candidate.get("gene_names", ()))
        aliases = tuple(candidate.get("aliases", ()))
        ko_ids = tuple(candidate.get("ko_ids", ()))
        candidate_terms = _unique_strings(
            (
                *genes,
                *aliases,
                *ko_ids,
                str(accession) if accession else "",
                str(protein_name) if protein_name else "",
            )
        )
        if not candidate_terms:
            continue
        candidate_identity = " OR ".join(
            _quote_term(term) or term for term in candidate_terms
        )
        input_anchor = primary_gene or f'"{protein.primary_accession}"'

        # Gene pairs are the highest-recall first rung for older articles;
        # long modern protein names are often absent from their metadata.
        unscoped_candidate = (
            genes[0]
            if genes
            else (
                _quote_term(str(protein_name))
                if protein_name
                else _quote_term(candidate_terms[0])
            )
        )
        unscoped = f"({input_anchor} AND {unscoped_candidate})"
        paired = (
            f"({input_anchor} AND {genes[0]})"
            if genes
            else f"(({identity}) AND ({candidate_identity}))"
        )
        organism_query = _organism_alias_query(
            candidate.get("organism_name")
            or getattr(protein, "organism_name", None),
            getattr(protein, "organism_aliases", ()),
        )
        scoped = paired
        if organism_query:
            scoped += f" AND ({organism_query})"
        mechanism = _mechanism_terms(str(candidate.get("role") or ""))

        ladder = [PlannedLiteratureQuery("PubMed", unscoped)]
        if aliases:
            alias_identity = " OR ".join(
                _quote_term(term) or term
                for term in _unique_strings(
                    (*genes, *aliases, str(protein_name or ""))
                )
            )
            alias_query = f"({alias_identity})"
            if organism_query:
                alias_query += f" AND ({organism_query})"
            ladder.append(PlannedLiteratureQuery("PubMed", alias_query))
        role_text = str(candidate.get("role") or "")
        if mode == "deep" and re.search(
            r"(?:electron|redox|reductase|ferredoxin|flavodoxin)",
            role_text,
            re.IGNORECASE,
        ):
            component_anchor = (
                genes[0]
                if genes
                else _quote_term(str(protein_name or candidate_terms[0]))
            )
            ladder.append(
                PlannedLiteratureQuery(
                    "PubMed",
                    f"({component_anchor}) AND ({mechanism})",
                )
            )
        if reaction_identity:
            ladder.append(
                PlannedLiteratureQuery(
                    "Europe PMC",
                    f"{scoped} AND ({reaction_identity}) AND "
                    "(activity OR reconstitution OR complementation)",
                )
            )
        ladder.append(
            PlannedLiteratureQuery(
                "PubMed",
                f"{scoped} AND ({mechanism})",
            )
        )
        if mode == "deep":
            ladder.extend(
                [
                    PlannedLiteratureQuery(
                        "Semantic Scholar",
                        f"{input_anchor} {candidate_terms[0]} "
                        f"{context.reaction.reaction_id} biochemical complex",
                    ),
                    PlannedLiteratureQuery(
                        "Europe PMC",
                        f"{scoped} AND (subunit OR complex OR assembly OR "
                        "protein-protein interaction)",
                    ),
                ]
            )
        ladders.append(ladder)

    if ladders:
        for depth in range(max(len(ladder) for ladder in ladders)):
            for ladder in ladders:
                if depth < len(ladder):
                    planned.append(ladder[depth])
    else:
        organism_query = _organism_alias_query(
            getattr(protein, "organism_name", None),
            getattr(protein, "organism_aliases", ()),
        )
        broad = f"({identity})"
        if organism_query:
            broad += f" AND ({organism_query})"
        planned.append(
            PlannedLiteratureQuery(
                "PubMed",
                broad + " AND (subunit OR complex OR reconstitution OR activity)",
            )
        )
        if reaction_identity:
            planned.append(
                PlannedLiteratureQuery(
                    "Europe PMC",
                    f"({identity}) AND ({reaction_identity}) AND "
                    "(subunit OR complex OR activity)",
                )
            )
        if mode == "deep":
            planned.extend(
                [
                    PlannedLiteratureQuery(
                        "PubMed",
                        f"({identity}) AND (knockout OR complementation OR "
                        "purified OR residual activity)",
                    ),
                    PlannedLiteratureQuery(
                        "Semantic Scholar",
                        f"{protein.primary_accession} protein subunit "
                        "biochemical activity",
                    ),
                ]
            )

    unique: list[PlannedLiteratureQuery] = []
    seen: set[tuple[str, str]] = set()
    for query in planned:
        key = (query.source, query.query.casefold())
        if key in seen:
            continue
        seen.add(key)
        unique.append(query)
        if len(unique) >= budget:
            break
    return unique


def linked_article_detail_calls(
    context: ResearchContext,
    hints: Sequence[VerifiedCandidateHint],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Fetch curated PMIDs directly without consuming search-query budget.

    Context-level references are considered first.  Candidate references are
    then interleaved so a long bibliography on the first candidate cannot
    crowd every other candidate out of the bounded direct-detail allowance.
    """

    limit = 6 if mode == "balanced" else 12
    pmids = list(_context_reference_pmids(context))
    hint_pmids = [list(_normalize_pmids(hint.reference_pmids)) for hint in hints]
    max_depth = max((len(values) for values in hint_pmids), default=0)
    for depth in range(max_depth):
        for values in hint_pmids:
            if depth < len(values):
                pmids.append(values[depth])
    return [
        PlannedToolCall(
            stage="literature",
            tool_name="PubMed_get_article",
            arguments={"pmid": pmid},
            rationale=(
                "Retrieve a curated linked PMID directly before literature "
                "discovery searches."
            ),
        )
        for pmid in _normalize_pmids(pmids)[:limit]
    ]


def verified_candidate_hints(
    records: Sequence[RawResearchEvidence],
    input_uniprot_id: str,
) -> list[VerifiedCandidateHint]:
    """Extract exact candidate identities from successful UniProt records."""

    hints: list[VerifiedCandidateHint] = []
    seen: set[str] = set()
    for record in records:
        if (
            record.status != "success"
            or record.tool_name != "UniProt_get_entry_by_accession"
        ):
            continue
        requested = record.query_arguments.get("accession")
        if not isinstance(requested, str):
            continue
        requested = requested.upper()
        if requested == input_uniprot_id.upper() or requested in seen:
            continue
        payload = _json_object(record.content)
        if payload is None:
            continue
        data = payload.get("data") if payload.get("status") == "success" else payload
        if not isinstance(data, Mapping):
            continue
        returned = data.get("primaryAccession")
        if not isinstance(returned, str) or returned.upper() != requested:
            continue
        protein_name = data.get("protein_name")
        if not isinstance(protein_name, str):
            protein_name = None
        raw_genes = data.get("gene_names")
        genes = tuple(
            value
            for value in raw_genes
            if isinstance(value, str) and value.strip()
        ) if isinstance(raw_genes, list) else ()
        organism = data.get("organism")
        organism_name = (
            organism.get("scientificName")
            if isinstance(organism, Mapping)
            and isinstance(organism.get("scientificName"), str)
            else None
        )
        taxon_id = (
            organism.get("taxonId")
            if isinstance(organism, Mapping)
            and isinstance(organism.get("taxonId"), int)
            else None
        )
        seen.add(requested)
        hints.append(
            VerifiedCandidateHint(
                uniprot_id=requested,
                protein_name=protein_name,
                gene_names=genes,
                organism_name=organism_name,
                taxon_id=taxon_id,
                aliases=_unique_strings(
                    (*genes, protein_name or "")
                ),
                seed_evidence_ids=(record.evidence_id,),
                reference_pmids=tuple(extract_pmids(record.content)),
            )
        )
    return hints


def candidate_identity_calls(
    records: Sequence[RawResearchEvidence],
    input_uniprot_id: str,
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Verify a few exact accessions discovered in curated complex records."""

    max_candidates = 3 if mode == "balanced" else 6
    accessions: list[str] = []
    for record in records:
        if (
            record.status != "success"
            or "ComplexPortal" not in record.tool_name
        ):
            continue
        for match in _UNIPROT_ACCESSION_PATTERN.finditer(record.content):
            accession = match.group(0).upper()
            if accession == input_uniprot_id.upper() or accession in accessions:
                continue
            accessions.append(accession)
            if len(accessions) >= max_candidates:
                break
        if len(accessions) >= max_candidates:
            break
    return [
        PlannedToolCall(
            stage="database",
            tool_name="UniProt_get_entry_by_accession",
            arguments={"accession": accession, "compact": True},
            rationale="Verify a candidate accession found in a curated complex.",
        )
        for accession in accessions
    ]


def report_candidate_identity_calls(
    report: Any | None,
    known_hints: Sequence[VerifiedCandidateHint],
    input_uniprot_id: str,
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Verify exact accessions proposed by any grounded dependency report."""

    if report is None:
        return []
    excluded = {
        input_uniprot_id.upper(),
        *(hint.uniprot_id.upper() for hint in known_hints),
    }
    limit = 3 if mode == "balanced" else 6
    accessions: list[str] = []
    for candidate in report.candidate_protein_dependencies:
        accession = candidate.uniprot_id
        if (
            not isinstance(accession, str)
            or accession.upper() in excluded
            or not _UNIPROT_ACCESSION_PATTERN.fullmatch(accession)
        ):
            continue
        accessions.append(accession.upper())
        excluded.add(accession.upper())
        if len(accessions) >= limit:
            break
    return [
        PlannedToolCall(
            stage="database",
            tool_name="UniProt_get_entry_by_accession",
            arguments={"accession": accession, "compact": True},
            rationale="Verify an exact candidate proposed by database analysis.",
        )
        for accession in accessions
    ]


def literature_tool_calls(
    queries: Sequence[PlannedLiteratureQuery],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Translate planned literature queries to exact MCP tool arguments."""

    result_limit = 5 if mode == "balanced" else 10
    calls: list[PlannedToolCall] = []
    for item in queries:
        if item.source == "PubMed":
            direct_pmid = _direct_pmid_query(item.query)
            if direct_pmid is not None:
                calls.append(
                    PlannedToolCall(
                        stage="literature",
                        tool_name="PubMed_get_article",
                        arguments={"pmid": direct_pmid},
                        rationale=(
                            "Retrieve a curated reference PMID directly before "
                            "discovery searches."
                        ),
                    )
                )
                continue
            calls.append(
                PlannedToolCall(
                    stage="literature",
                    tool_name="PubMed_search_articles",
                    arguments={
                        "query": item.query,
                        "limit": result_limit,
                        "include_abstract": True,
                    },
                    rationale="Search a bounded candidate-specific PubMed query.",
                )
            )
        elif item.source == "Europe PMC":
            calls.append(
                PlannedToolCall(
                    stage="literature",
                    tool_name="EuropePMC_search_articles",
                    arguments={
                        "query": item.query,
                        "limit": result_limit,
                        "require_has_ft": False,
                        "enrich_missing_abstract": True,
                    },
                    rationale="Search Europe PMC for the same dependency hypothesis.",
                )
            )
        else:
            calls.append(
                PlannedToolCall(
                    stage="literature",
                    tool_name="SemanticScholar_search_papers",
                    arguments={"query": item.query, "limit": result_limit},
                    rationale="Deep-mode fallback discovery only.",
                )
            )
    return calls


def article_detail_calls(
    search_records: Sequence[RawResearchEvidence],
    mode: ResearchMode,
    context: ResearchContext | None = None,
    candidate_hints: Sequence[VerifiedCandidateHint] = (),
) -> list[PlannedToolCall]:
    """Rank individual papers and reserve coverage for verified candidates."""

    max_articles = 8 if mode == "balanced" else 18
    ranked: dict[str, tuple[int, int, set[str]]] = {}
    sequence = 0
    for record in search_records:
        if record.status != "success":
            continue
        if record.tool_name == "PubMed_get_article":
            # A curated PMID rung already fetched the full detail record.
            continue
        articles = _article_mappings(record.content)
        if not articles:
            articles = [{"pmid": pmid} for pmid in extract_pmids(record.content)]
        for article in articles:
            pmid = str(article.get("pmid") or article.get("uid") or "")
            if not re.fullmatch(r"\d{5,9}", pmid):
                continue
            score, matched_candidates = _article_priority(
                article,
                record,
                context,
                candidate_hints,
            )
            existing = ranked.get(pmid)
            if existing is None or score > existing[0]:
                ranked[pmid] = (score, sequence, matched_candidates)
            sequence += 1

    ordered = sorted(
        ranked,
        key=lambda pmid: (-ranked[pmid][0], ranked[pmid][1], pmid),
    )
    selected: list[str] = []
    for hint in candidate_hints:
        best = next(
            (
                pmid
                for pmid in ordered
                if hint.uniprot_id in ranked[pmid][2]
                and pmid not in selected
            ),
            None,
        )
        if best is not None:
            selected.append(best)
        if len(selected) >= max_articles:
            break
    selected.extend(
        pmid
        for pmid in ordered
        if pmid not in selected
    )
    pmids = selected[:max_articles]
    return [
        PlannedToolCall(
            stage="literature",
            tool_name="PubMed_get_article",
            arguments={"pmid": pmid},
            rationale="Read one top-ranked unique PubMed article.",
        )
        for pmid in pmids
    ]


def fulltext_metadata_calls(
    literature_records: Sequence[RawResearchEvidence],
    context: ResearchContext,
    candidate_hints: Sequence[VerifiedCandidateHint],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Resolve PMCIDs for paired incomplete PubMed abstracts by exact PMID."""

    if not candidate_hints:
        return []
    pmids_with_pmcid = {
        str(article.get("pmid") or article.get("uid"))
        for record in literature_records
        if record.status == "success"
        for article in _article_mappings(record.content)
        if _normalize_pmcid(article.get("pmcid")) is not None
    }
    attempted_pmids = {
        match.group(1)
        for record in literature_records
        if record.tool_name == "EuropePMC_search_articles"
        and (
            match := re.fullmatch(
                r"\s*EXT_ID\s*:\s*(\d{5,9})\s*",
                str(record.query_arguments.get("query") or ""),
                re.I,
            )
        )
        is not None
    }
    ranked: dict[str, tuple[int, int]] = {}
    sequence = 0
    for record in literature_records:
        if record.status != "success" or "PubMed" not in record.tool_name:
            continue
        for article in _article_mappings(record.content):
            pmid = str(article.get("pmid") or article.get("uid") or "")
            if (
                not re.fullmatch(r"\d{5,9}", pmid)
                or pmid in pmids_with_pmcid
                or pmid in attempted_pmids
            ):
                continue
            text = " ".join(
                str(article.get(field) or "")
                for field in ("title", "abstract")
            )
            if not _abstract_needs_fulltext_escalation(text, context):
                continue
            scores = [
                _article_priority(article, record, context, (hint,))[0]
                for hint in candidate_hints
                if _fulltext_article_matches_pair(text, context, hint)
            ]
            if not scores:
                continue
            score = max(scores)
            existing = ranked.get(pmid)
            if existing is None or score > existing[0]:
                ranked[pmid] = (score, sequence)
            sequence += 1

    max_articles = 2 if mode == "balanced" else 6
    selected = sorted(
        ranked,
        key=lambda pmid: (-ranked[pmid][0], ranked[pmid][1], pmid),
    )[:max_articles]
    return [
        PlannedToolCall(
            stage="literature",
            tool_name="EuropePMC_search_articles",
            arguments={
                "query": f"EXT_ID:{pmid}",
                "limit": 1,
                "require_has_ft": False,
                "enrich_missing_abstract": False,
            },
            rationale=(
                "Resolve a PMCID by exact PMID before attempting bounded "
                "full-text snippets for an incomplete paired abstract."
            ),
        )
        for pmid in selected
    ]


def fulltext_snippet_calls(
    literature_records: Sequence[RawResearchEvidence],
    context: ResearchContext,
    candidate_hints: Sequence[VerifiedCandidateHint],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Escalate relevant but incomplete abstracts to bounded full-text snippets.

    Europe PMC full text is queried only when its metadata supplies a PMCID,
    the abstract names both verified protein identities, and the abstract is
    still missing an exact organism or a complete experimental dependency
    relation.  A failed snippet call is considered attempted for this workflow
    and is never silently retried or replaced by a weaker source.
    """

    if not candidate_hints:
        return []
    attempted_pmcids = {
        pmcid
        for record in literature_records
        if record.tool_name == "EuropePMC_get_fulltext_snippets"
        and (
            pmcid := _normalize_pmcid(
                record.query_arguments.get("pmcid")
                or record.query_arguments.get("article_id")
            )
        )
        is not None
    }
    ranked: dict[str, tuple[int, int, tuple[str, ...]]] = {}
    sequence = 0
    for record in literature_records:
        if record.status != "success":
            continue
        for article in _article_mappings(record.content):
            pmcid = _normalize_pmcid(article.get("pmcid"))
            if pmcid is None or pmcid in attempted_pmcids:
                continue
            text = " ".join(
                str(article.get(field) or "")
                for field in ("title", "abstract")
            )
            for hint in candidate_hints:
                if not _fulltext_article_matches_pair(text, context, hint):
                    continue
                if not _abstract_needs_fulltext_escalation(
                    text,
                    context,
                ):
                    continue
                score, _ = _article_priority(
                    article,
                    record,
                    context,
                    (hint,),
                )
                terms = _fulltext_query_terms(context, hint, text)
                if not terms:
                    continue
                existing = ranked.get(pmcid)
                if existing is None:
                    ranked[pmcid] = (score, sequence, terms)
                else:
                    ranked[pmcid] = (
                        max(score, existing[0]),
                        existing[1],
                        _unique_strings((*existing[2], *terms))[:12],
                    )
                sequence += 1

    max_articles = 2 if mode == "balanced" else 6
    selected = sorted(
        ranked,
        key=lambda pmcid: (-ranked[pmcid][0], ranked[pmcid][1], pmcid),
    )[:max_articles]
    return [
        PlannedToolCall(
            stage="literature",
            tool_name="EuropePMC_get_fulltext_snippets",
            arguments={
                "pmcid": pmcid,
                "terms": list(ranked[pmcid][2]),
                "window_chars": 1200,
                "max_snippets_per_term": 3,
                "max_total_chars": 20_000,
            },
            rationale=(
                "Read bounded verbatim full-text windows because the paired "
                "abstract does not establish the exact organism and complete "
                "protein-level dependency experiment."
            ),
        )
        for pmcid in selected
    ]


def _normalize_pmcid(value: Any) -> str | None:
    if not isinstance(value, (str, int)):
        return None
    match = re.fullmatch(r"\s*(?:PMC)?(\d{3,12})\s*", str(value), re.I)
    return f"PMC{match.group(1)}" if match is not None else None


def _fulltext_article_matches_pair(
    text: str,
    context: ResearchContext,
    hint: VerifiedCandidateHint,
) -> bool:
    normalized = text.casefold()
    input_terms = (
        context.protein.primary_accession,
        *context.protein.gene_names,
        *context.protein.protein_names,
    )
    candidate_terms = (
        hint.uniprot_id,
        *hint.gene_names,
        *hint.aliases,
        *((hint.protein_name,) if hint.protein_name else ()),
    )
    return _any_search_term_occurs(
        normalized,
        input_terms,
    ) and _any_search_term_occurs(normalized, candidate_terms)


def _abstract_needs_fulltext_escalation(
    text: str,
    context: ResearchContext,
) -> bool:
    normalized = text.casefold()
    exact_organisms = _unique_strings(
        (
            _species_search_name(context.protein.organism_name) or "",
            *(
                alias
                for alias in context.protein.organism_aliases
                if len(alias.split()) >= 2
                and not alias.casefold().startswith("strain ")
            ),
        )
    )
    has_exact_organism = any(
        _any_search_term_occurs(normalized, (organism,))
        for organism in exact_organisms
    )
    has_complete_reconstitution = bool(
        re.search(r"(?:reconstitut|complement|add(?:ed|ition))", normalized)
        and re.search(r"(?:activity|cataly|holoenzyme)", normalized)
        and re.search(r"(?:alone|without|absence|delet|knockout)", normalized)
    )
    has_complete_genetic_loss = bool(
        re.search(
            r"(?:null\s+mutation|gene\s+(?:delet|knockout|disrupt)|"
            r"blocking\s+(?:the\s+)?expression)",
            normalized,
        )
        and re.search(
            r"(?:no\s+(?:detectable\s+)?activity|abolish|inactivat|"
            r"unable\s+to\s+grow|growth\s+(?:was\s+)?(?:completely\s+)?blocked)",
            normalized,
        )
    )
    return not (
        has_exact_organism
        and (has_complete_reconstitution or has_complete_genetic_loss)
    )


def _fulltext_query_terms(
    context: ResearchContext,
    hint: VerifiedCandidateHint,
    article_text: str,
) -> tuple[str, ...]:
    """Build compact identity/lineage/mechanism terms without case hardcoding."""

    terms: list[str] = []
    organism_name = context.protein.organism_name or ""
    for match in re.finditer(r"\(([^()]*)\)", organism_name):
        parenthetical = re.sub(
            r"\b(?:strain|substrain|isolate)\b",
            " ",
            match.group(1),
            flags=re.I,
        )
        terms.extend(
            token
            for token in re.findall(r"[A-Za-z0-9][A-Za-z0-9_.-]*", parenthetical)
            if len(token) >= 2
        )

    canonical = _species_search_name(organism_name)
    canonical_key = canonical.casefold() if canonical else ""
    for alias in context.protein.organism_aliases:
        alias_species = _species_search_name(alias)
        if (
            alias_species is None
            or len(alias_species.split()) < 2
            or alias_species.casefold() == canonical_key
        ):
            continue
        terms.append(alias_species.split()[1])

    terms.extend(hint.gene_names[:2] or (hint.uniprot_id,))
    terms.extend(
        context.protein.gene_names[:2]
        or (context.protein.primary_accession,)
    )
    normalized = article_text.casefold()
    mechanism_terms = (
        "expression",
        "heat-sensitive",
        "temperature-sensitive",
        "knockout",
        "gene deletion",
        "loss-of-function",
        "complementation",
        "reconstitution",
        "no activity",
    )
    matched_mechanisms = [
        term
        for term in mechanism_terms
        if term in normalized
        or term.replace("-", " ") in normalized
    ][:2]
    terms.extend(matched_mechanisms or ["activity"])
    return _unique_strings(terms)[:12]


def _article_mappings(content: str) -> list[Mapping[str, Any]]:
    try:
        payload = json.loads(content)
    except (TypeError, ValueError):
        return []
    found: list[Mapping[str, Any]] = []

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            if value.get("pmid") is not None or value.get("uid") is not None:
                found.append(value)
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(payload)
    return found


def _article_priority(
    article: Mapping[str, Any],
    record: RawResearchEvidence,
    context: ResearchContext | None,
    candidate_hints: Sequence[VerifiedCandidateHint],
) -> tuple[int, set[str]]:
    title = str(article.get("title") or "")
    abstract = str(article.get("abstract") or "")
    text = f"{title} {abstract}".casefold()
    score = _literature_detail_priority(record)
    input_match = False
    if context is not None:
        input_terms = [
            context.protein.primary_accession,
            *context.protein.gene_names,
            *(context.protein.protein_names[:1]),
            *(
                alias
                for annotation in context.protein.subunit_annotations
                for alias in _gene_cluster_aliases(annotation)
            ),
        ]
        input_match = _any_search_term_occurs(text, input_terms)
        if input_match:
            score += 5
        species = _species_search_name(context.protein.organism_name)
        if species and species.casefold() in text:
            score += 3
    matched: set[str] = set()
    for hint in candidate_hints:
        terms = [
            hint.uniprot_id,
            *hint.gene_names,
            *hint.aliases,
            *hint.ko_ids,
        ]
        if hint.protein_name:
            terms.append(hint.protein_name)
        if _any_search_term_occurs(text, terms):
            matched.add(hint.uniprot_id)
            score += 5
            if input_match:
                score += 8
    if re.search(
        r"(?:reconstitut|complement|addback|add(?:ed|ition)|knockout|"
        r"purif|isolat|activity|cataly)",
        text,
    ):
        score += 3
    if "reconstitution" in title.casefold():
        score += 3
    return score, matched


def _any_search_term_occurs(text: str, terms: Sequence[str]) -> bool:
    return any(
        term
        and re.search(
            rf"(?<![a-z0-9]){re.escape(term.casefold())}(?![a-z0-9])",
            text,
        )
        for term in terms
    )


def _literature_detail_priority(record: RawResearchEvidence) -> int:
    """Prefer candidate-pair and explicit dependency search results."""

    if record.status != "success" or "PubMed" not in record.tool_name:
        return -1
    query = str(record.query_arguments.get("query") or "").casefold()
    content = record.content.casefold()
    score = 0
    if "residual activity" in query or "reconstitution" in query:
        score += 4
    if query.count(" and ") >= 2:
        score += 2
    if "residual activity" in content or "approximately 5%" in content:
        score += 3
    if "small subunit" in content and "large subunit" in content:
        score += 2
    return score


def build_web_call_plan(
    context: ResearchContext,
    database_report: BioDatabaseResearchResult | None,
    mode: ResearchMode,
    literature_report: Any | None = None,
) -> list[PlannedToolCall]:
    """Build a last-resort exact web search; pages are fetched separately."""

    protein = context.protein
    candidate_terms: list[str] = []
    if database_report is not None:
        for candidate in database_report.candidate_protein_dependencies:
            candidate_terms.extend(
                value
                for value in (candidate.uniprot_id, candidate.protein_name)
                if value
            )
    if literature_report is not None:
        for candidate in literature_report.candidate_protein_dependencies:
            candidate_terms.extend(
                value
                for value in (candidate.uniprot_id, candidate.protein_name)
                if value
            )
    identity = " ".join(
        [
            protein.primary_accession,
            *(protein.gene_names[:1]),
            *(candidate_terms[:2]),
            "protein complex subunit UniProt",
        ]
    )
    limit = 5 if mode == "balanced" else 8
    return [
        PlannedToolCall(
            stage="web",
            tool_name="search",
            arguments={
                "query": identity,
                "limit": limit,
                "engines": ["bing"],
                "searchMode": "request",
            },
            rationale="Resolve only an identifier or official-source conflict.",
        )
    ]


def build_host_call_plan(
    context: ResearchContext,
    candidate_roles: Sequence[Mapping[str, Any]],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Build bounded host-mapping calls for already supported native roles."""

    calls: list[PlannedToolCall] = [
        PlannedToolCall(
            stage="host_compatibility",
            tool_name="find_iml1515_reactions",
            arguments={
                "kegg_reaction_id": context.reaction.reaction_id,
                "max_results": 10,
            },
            rationale="Read exact iML1515 reaction and raw GPR context.",
        )
    ]
    max_candidates = 3 if mode == "balanced" else 6
    for role in candidate_roles[:max_candidates]:
        uniprot_id = role.get("uniprot_id")
        if not isinstance(uniprot_id, str) or not uniprot_id:
            continue
        calls.extend(
            [
                PlannedToolCall(
                    stage="host_compatibility",
                    tool_name="UniProt_get_entry_by_accession",
                    arguments={"accession": uniprot_id, "compact": True},
                    rationale="Verify the candidate identity and organism.",
                ),
                PlannedToolCall(
                    stage="host_compatibility",
                    tool_name="get_iml1515_gene",
                    arguments={"identifier": uniprot_id, "max_results": 10},
                    rationale="Check exact UniProt-to-MG1655 model mapping.",
                ),
            ]
        )
        for ko_id in role.get("ko_ids", [])[:3]:
            if not isinstance(ko_id, str) or not _KO_PATTERN.fullmatch(ko_id):
                continue
            calls.append(
                PlannedToolCall(
                    stage="host_compatibility",
                    tool_name="KEGG_link_entries",
                    arguments={
                        "source": f"ko:{ko_id.upper()}",
                        "target": "eco",
                    },
                    rationale=(
                        "Check whether the exact native component KO has an "
                        "E. coli gene mapping; transport errors are not absence."
                    ),
                )
            )
        taxon_id = role.get("taxon_id")
        if taxon_id not in {511145, 83333, 562}:
            calls.extend(
                [
                    PlannedToolCall(
                        stage="host_compatibility",
                        tool_name="OMA_resolve_xref",
                        arguments={"search": uniprot_id, "limit": 10},
                        rationale="Resolve the native protein in OMA.",
                    ),
                    PlannedToolCall(
                        stage="host_compatibility",
                        tool_name="OMA_get_orthologs",
                        arguments={
                            "protein_id": uniprot_id,
                            "per_page": 30 if mode == "balanced" else 60,
                        },
                        rationale="Find ortholog candidates before any BLAST fallback.",
                    ),
                ]
            )
    return _deduplicate_calls(calls)


def deep_host_sequence_calls(
    candidate_roles: Sequence[Mapping[str, Any]],
    host_records: Sequence[RawResearchEvidence],
) -> list[PlannedToolCall]:
    """Request native sequences only when OMA found no E. coli mapping."""

    oma_content = "\n".join(
        record.content
        for record in host_records
        if record.status == "success" and record.tool_name.startswith("OMA_")
    ).casefold()
    has_ecoli_oma_hit = any(
        marker in oma_content
        for marker in (
            "escherichia coli",
            '"taxon_id": 562',
            '"taxon_id": 511145',
            '"species": "ecoli"',
        )
    )
    if has_ecoli_oma_hit:
        return []
    calls: list[PlannedToolCall] = []
    for role in candidate_roles[:3]:
        uniprot_id = role.get("uniprot_id")
        taxon_id = role.get("taxon_id")
        if (
            isinstance(uniprot_id, str)
            and uniprot_id
            and taxon_id not in {511145, 83333, 562}
        ):
            calls.append(
                PlannedToolCall(
                    stage="host_compatibility",
                    tool_name="UniProt_get_sequence_by_accession",
                    arguments={"accession": uniprot_id},
                    rationale="Obtain native sequence for a deep-mode BLAST fallback.",
                )
            )
    return _deduplicate_calls(calls)


def blast_calls_from_sequences(
    sequence_records: Sequence[RawResearchEvidence],
) -> list[PlannedToolCall]:
    """Build BLAST calls from successfully retrieved native sequences."""

    calls: list[PlannedToolCall] = []
    for record in sequence_records:
        if (
            record.status != "success"
            or record.tool_name != "UniProt_get_sequence_by_accession"
        ):
            continue
        sequence = _extract_protein_sequence(record.content)
        if sequence is None:
            continue
        calls.append(
            PlannedToolCall(
                stage="host_compatibility",
                tool_name="BLAST_protein_search",
                arguments={
                    "sequence": sequence,
                    "database": "swissprot",
                    "expect": 1e-5,
                    "hitlist_size": 20,
                },
                rationale="Deep-mode fallback after OMA coverage was insufficient.",
            )
        )
    return calls


def host_candidate_verification_calls(
    host_records: Sequence[RawResearchEvidence],
    known_native_ids: Sequence[str],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Verify exact accessions discovered by OMA or BLAST before analysis."""

    max_candidates = 3 if mode == "balanced" else 6
    excluded = {value.upper() for value in known_native_ids}
    accessions: list[str] = []
    for record in host_records:
        if (
            record.status != "success"
            or record.tool_name
            not in {"OMA_get_orthologs", "BLAST_protein_search"}
        ):
            continue
        for match in _UNIPROT_ACCESSION_PATTERN.finditer(record.content):
            accession = match.group(0).upper()
            if accession in excluded or accession in accessions:
                continue
            accessions.append(accession)
            if len(accessions) >= max_candidates:
                break
        if len(accessions) >= max_candidates:
            break
    calls: list[PlannedToolCall] = []
    for accession in accessions:
        calls.extend(
            [
                PlannedToolCall(
                    stage="host_compatibility",
                    tool_name="UniProt_get_entry_by_accession",
                    arguments={"accession": accession, "compact": True},
                    rationale="Verify an OMA/BLAST host candidate in UniProt.",
                ),
                PlannedToolCall(
                    stage="host_compatibility",
                    tool_name="get_iml1515_gene",
                    arguments={"identifier": accession, "max_results": 10},
                    rationale="Map the verified host candidate to iML1515.",
                ),
            ]
        )
    return _deduplicate_calls(calls)


def web_fetch_calls(
    search_records: Sequence[RawResearchEvidence],
    mode: ResearchMode,
) -> list[PlannedToolCall]:
    """Fetch a small number of official or scholarly search results."""

    max_pages = 2 if mode == "balanced" else 4
    urls: list[str] = []
    preferred_hosts = (
        "uniprot.org",
        "kegg.jp",
        "rhea-db.org",
        "ebi.ac.uk",
        "ncbi.nlm.nih.gov",
        "europepmc.org",
    )
    discovered = [
        url
        for record in search_records
        if record.status == "success"
        for url in _extract_urls(record.content)
    ]
    for url in sorted(
        dict.fromkeys(discovered),
        key=lambda value: (
            not any(host in value.casefold() for host in preferred_hosts),
            value,
        ),
    ):
        urls.append(url)
        if len(urls) >= max_pages:
            break
    return [
        PlannedToolCall(
            stage="web",
            tool_name="fetchWebContent",
            arguments={
                "url": url.rstrip(".,);]"),
                "maxChars": 20_000,
                "readability": True,
                "includeLinks": True,
            },
            rationale="Read the body of a selected official or scholarly page.",
        )
        for url in urls
    ]


def extract_pmids(content: str) -> list[str]:
    values: list[str] = []
    for pattern in _PMID_PATTERNS:
        for value in pattern.findall(content):
            if value not in values:
                values.append(value)
    return values


def evidence_bundle(
    records: Sequence[RawResearchEvidence],
) -> list[dict[str, Any]]:
    return [record.model_dump(mode="json") for record in records]


def _context_kegg_organism_code(context: Any) -> str | None:
    for owner in (context, getattr(context, "protein", None)):
        if owner is None:
            continue
        explicit = getattr(owner, "kegg_organism_code", None)
        if isinstance(explicit, str) and re.fullmatch(
            r"[a-z][a-z0-9_]{1,9}", explicit, re.IGNORECASE
        ):
            return explicit.casefold()
    protein = getattr(context, "protein", None)
    for crossref in getattr(protein, "cross_references", ()):
        database = getattr(crossref, "database", None)
        record_id = getattr(crossref, "record_id", None)
        if isinstance(crossref, Mapping):
            database = database or crossref.get("database")
            record_id = record_id or crossref.get("record_id") or crossref.get("id")
        if (
            not isinstance(database, str)
            or database.casefold() != "kegg"
            or not isinstance(record_id, str)
            or ":" not in record_id
        ):
            continue
        prefix = record_id.split(":", 1)[0].strip()
        if re.fullmatch(r"[a-z][a-z0-9_]{1,9}", prefix, re.IGNORECASE):
            return prefix.casefold()
    return None


def _normalize_candidate_text(value: str) -> str:
    return " ".join(value.replace('"', " ").split()).casefold()


def _gene_cluster_aliases(value: str) -> tuple[str, ...]:
    """Derive compact historical operon tokens such as ``benABC``.

    Older abstracts often mention only a concatenated gene-cluster name while
    modern curated records list the individual subunits.  The derivation is
    restricted to two or more explicit mixed-case gene symbols sharing the
    same stem in one curated span; arbitrary prose is never combined.
    """

    grouped: dict[str, tuple[str, list[str]]] = {}
    for match in _GENE_TOKEN_PATTERN.finditer(value):
        gene = match.group(0)
        if gene[0].isupper():
            gene = gene[0].lower() + gene[1:]
        parts = re.fullmatch(r"(.+?)([A-Z]\d?)", gene)
        if parts is None:
            continue
        stem, suffix = parts.groups()
        key = stem.casefold()
        if key not in grouped:
            grouped[key] = (stem, [])
        suffixes = grouped[key][1]
        if suffix not in suffixes:
            suffixes.append(suffix)
    return tuple(
        stem + "".join(suffixes)
        for stem, suffixes in grouped.values()
        if len(suffixes) >= 2
    )


def _normalize_seed_value(
    kind: CandidateSeedKind,
    raw_value: Any,
) -> str | None:
    if not isinstance(raw_value, str):
        return None
    value = " ".join(raw_value.strip().split())
    if not value:
        return None
    if kind == "accession":
        match = _UNIPROT_ACCESSION_PATTERN.fullmatch(value)
        return match.group(0).upper() if match is not None else None
    if kind == "ko":
        match = re.fullmatch(r"(?:ko:)?(K\d{5})", value, re.IGNORECASE)
        return match.group(1).upper() if match is not None else None
    if kind in {"gene", "locus"}:
        value = value.strip(".,;:()[]{}")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{0,39}", value):
            return None
        return value
    return value if len(value) <= 300 else None


def _bounded_source_span(value: str | None, limit: int = 400) -> str:
    normalized = " ".join(str(value or "").split())
    return normalized[:limit]


def _iter_candidate_values(value: Any):
    if isinstance(value, str):
        if value.strip():
            yield value.strip()
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _iter_candidate_values(item)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for item in value:
            yield from _iter_candidate_values(item)


def _structured_candidate_values(
    payload: Any,
) -> list[tuple[CandidateSeedKind, str, str, str | None]]:
    """Read candidate fields from structured curated records only."""

    results: list[tuple[CandidateSeedKind, str, str, str | None]] = []
    accession_keys = {
        "accession",
        "primaryaccession",
        "uniprotaccession",
        "uniprotid",
    }
    gene_keys = {
        "gene",
        "genename",
        "genenames",
        "genesymbol",
        "genesymbols",
        "symbol",
    }
    locus_keys = {"locus", "locustag", "locustags", "orderedlocusname"}
    ko_keys = {
        "ko",
        "koid",
        "koids",
        "orthology",
        "orthologyid",
        "keggorthology",
    }
    protein_keys = {
        "proteinname",
        "proteinfullname",
        "recommendedname",
        "fullname",
    }

    def normalized_key(value: Any) -> str:
        return re.sub(r"[^a-z0-9]", "", str(value).casefold())

    def role_for(mapping: Mapping[str, Any]) -> str | None:
        for raw_key, raw_value in mapping.items():
            if normalized_key(raw_key) not in {
                "role",
                "function",
                "description",
                "proteinname",
                "name",
                "comment",
            }:
                continue
            value = next(_iter_candidate_values(raw_value), None)
            if isinstance(value, str):
                return _bounded_source_span(value)
        return None

    def add_values(
        kind: CandidateSeedKind,
        raw: Any,
        path: str,
        role_hint: str | None,
    ) -> None:
        for value in _iter_candidate_values(raw):
            if kind == "accession":
                matches = [match.group(0) for match in _UNIPROT_ACCESSION_PATTERN.finditer(value)]
            elif kind == "ko":
                matches = [match.group(1) for match in _KO_PATTERN.finditer(value)]
            elif kind in {"gene", "locus"}:
                matches = [
                    token.strip(".,;:()[]{}")
                    for token in re.split(r"[\s,;/|]+", value)
                    if token.strip(".,;:()[]{}")
                ]
            else:
                matches = [value]
            for match in matches:
                if _normalize_seed_value(kind, match) is not None:
                    results.append(
                        (
                            kind,
                            match,
                            _bounded_source_span(f"{path}: {value}"),
                            role_hint,
                        )
                    )

    def visit(value: Any, path: str) -> None:
        if isinstance(value, Mapping):
            role_hint = role_for(value)
            mapping_keys = {normalized_key(key) for key in value}
            participant_like = bool(
                mapping_keys
                & (accession_keys | gene_keys | locus_keys | {"role", "function"})
            ) or any(
                marker in path.casefold()
                for marker in ("participant", "component", "subunit", "interactor")
            )
            for raw_key, child in value.items():
                key = normalized_key(raw_key)
                child_path = f"{path}.{raw_key}"
                if key in accession_keys:
                    add_values("accession", child, child_path, role_hint)
                elif key in gene_keys:
                    add_values("gene", child, child_path, role_hint)
                elif key in locus_keys:
                    add_values("locus", child, child_path, role_hint)
                elif key in ko_keys:
                    add_values("ko", child, child_path, role_hint)
                elif key in protein_keys and participant_like:
                    add_values("protein_name", child, child_path, role_hint)
                visit(child, child_path)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            for index, child in enumerate(value):
                visit(child, f"{path}[{index}]")

    visit(payload, "$")
    return results


def _select_candidate_seeds(
    seeds: Sequence[CandidateSeed],
    mode: ResearchMode,
) -> list[CandidateSeed]:
    limit = 3 if mode == "balanced" else 6
    priorities = {
        "accession": 0,
        "gene": 1,
        "locus": 2,
        "protein_name": 3,
        "ko": 4,
    }
    ordered = sorted(
        enumerate(seeds),
        key=lambda item: (priorities[item[1].kind], item[0]),
    )
    selected: list[CandidateSeed] = []
    seen: set[tuple[str, str, int | None, str | None]] = set()
    for _, seed in ordered:
        key = (
            seed.kind,
            _normalize_candidate_text(seed.value),
            seed.taxon_id,
            seed.kegg_organism_code,
        )
        if key in seen:
            continue
        seen.add(key)
        selected.append(seed)
        if len(selected) >= limit:
            break
    return selected


def _uniprot_query_value(value: str) -> str:
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{0,39}", value):
        return value
    return '"' + value.replace('"', " ").strip() + '"'


def _extract_kegg_gene_ids(content: str) -> list[str]:
    values: list[str] = []
    payload = _json_value(content)

    def visit(value: Any) -> None:
        if isinstance(value, str):
            for match in _KEGG_GENE_PATTERN.finditer(value):
                gene_id = match.group(1)
                if gene_id not in values:
                    values.append(gene_id)
        elif isinstance(value, Mapping):
            for child in value.values():
                visit(child)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            for child in value:
                visit(child)

    visit(payload if payload is not None else content)
    return values


def _mapping_gene_names(mapping: Mapping[str, Any]) -> tuple[str, ...]:
    values: list[str] = []
    for key in ("gene_names", "geneNames", "genes", "gene"):
        raw = mapping.get(key)
        if raw is None:
            continue
        if isinstance(raw, str):
            values.extend(raw.split())
            continue
        for item in _iter_candidate_values(raw):
            # Raw Proteins API gene objects also contain type/status strings;
            # keep only compact identifier-shaped values.
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{0,39}", item):
                values.append(item)
    return _unique_strings(values)


def _mapping_protein_name(mapping: Mapping[str, Any]) -> str | None:
    direct = mapping.get("protein_name")
    if isinstance(direct, str) and direct.strip():
        return direct.strip()
    for key in ("proteinDescription", "protein", "recommendedName"):
        raw = mapping.get(key)
        if raw is None:
            continue
        if isinstance(raw, str) and raw.strip():
            return raw.strip()
        strings = [
            value
            for value in _iter_candidate_values(raw)
            if isinstance(value, str)
            and len(value.strip()) > 3
            and value.casefold() not in {"recommended", "alternative"}
        ]
        if strings:
            return strings[0].strip()
    return None


def _mapping_organism(
    mapping: Mapping[str, Any],
) -> tuple[str | None, int | None]:
    name: str | None = None
    taxon_id: int | None = None
    raw = mapping.get("organism")
    if isinstance(raw, str):
        name = raw
    elif isinstance(raw, Mapping):
        for key in ("scientificName", "name"):
            value = raw.get(key)
            if isinstance(value, str) and value.strip():
                name = value.strip()
                break
        for key in ("taxonId", "taxonomy", "taxonomyId", "id"):
            value = raw.get(key)
            if isinstance(value, int):
                taxon_id = value
                break
            if isinstance(value, str) and value.isdigit():
                taxon_id = int(value)
                break
    if name is None:
        value = mapping.get("organism_name")
        if isinstance(value, str) and value.strip():
            name = value.strip()
    if taxon_id is None:
        for key in ("organism_id", "taxon_id", "taxonId"):
            value = mapping.get(key)
            if isinstance(value, int):
                taxon_id = value
                break
            if isinstance(value, str) and value.isdigit():
                taxon_id = int(value)
                break
    return name, taxon_id


def _extract_search_candidate_mappings(
    content: str,
) -> list[dict[str, Any]]:
    payload = _json_value(content)
    if payload is None:
        return []
    found: list[dict[str, Any]] = []
    seen: set[str] = set()

    def visit(value: Any) -> None:
        if isinstance(value, Mapping):
            accession = value.get("primaryAccession") or value.get("accession")
            if isinstance(accession, str) and _UNIPROT_ACCESSION_PATTERN.fullmatch(
                accession
            ):
                upper = accession.upper()
                if upper not in seen:
                    organism_name, taxon_id = _mapping_organism(value)
                    raw_reviewed = value.get("reviewed")
                    entry_type = value.get("entryType")
                    reviewed: bool | None = None
                    if isinstance(raw_reviewed, bool):
                        reviewed = raw_reviewed
                    elif isinstance(raw_reviewed, str):
                        reviewed = raw_reviewed.casefold() in {
                            "true",
                            "reviewed",
                            "yes",
                        }
                    elif isinstance(entry_type, str):
                        reviewed = "unreviewed" not in entry_type.casefold()
                    found.append(
                        {
                            "accession": upper,
                            "gene_names": _mapping_gene_names(value),
                            "protein_name": _mapping_protein_name(value),
                            "organism_name": organism_name,
                            "taxon_id": taxon_id,
                            "reviewed": reviewed,
                        }
                    )
                    seen.add(upper)
            for child in value.values():
                visit(child)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            for child in value:
                visit(child)

    visit(payload)
    return found


def _extract_converted_accessions(content: str) -> list[str]:
    payload = _json_value(content)
    values: list[str] = []

    def visit(value: Any) -> None:
        if isinstance(value, str):
            for match in _UNIPROT_ACCESSION_PATTERN.finditer(value):
                accession = match.group(0).upper()
                if accession not in values:
                    values.append(accession)
        elif isinstance(value, Mapping):
            for key, child in value.items():
                if str(key).casefold() in {
                    "external_id",
                    "externalid",
                    "uniprot",
                    "accession",
                }:
                    visit(child)
                elif isinstance(child, (Mapping, list, tuple)):
                    visit(child)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            for child in value:
                visit(child)

    visit(payload if payload is not None else content)
    return values


def _record_matches_seed_search(
    record: RawResearchEvidence,
    seed: CandidateSeed,
) -> bool:
    if record.tool_name == "proteins_api_search":
        query = record.query_arguments.get("query")
        taxid = record.query_arguments.get("taxid")
        return (
            isinstance(query, str)
            and _normalize_candidate_text(query)
            == _normalize_candidate_text(seed.value)
            and (
                seed.taxon_id is None
                or str(taxid) == str(seed.taxon_id)
            )
        )
    if record.tool_name == "UniProt_search":
        query = str(record.query_arguments.get("query") or "")
        return (
            _normalize_candidate_text(seed.value)
            in _normalize_candidate_text(query)
            and (
                seed.taxon_id is None
                or f"taxonomy_id:{seed.taxon_id}" in query.casefold()
            )
        )
    return False


def _search_mapping_matches_seed(
    mapping: Mapping[str, Any],
    seed: CandidateSeed,
) -> bool:
    taxon_id = mapping.get("taxon_id")
    if seed.taxon_id is not None and taxon_id is not None and taxon_id != seed.taxon_id:
        return False
    if mapping.get("reviewed") is False:
        return False
    target = _normalize_candidate_text(seed.value)
    if seed.kind in {"gene", "locus"}:
        genes = mapping.get("gene_names") or ()
        return any(_normalize_candidate_text(str(gene)) == target for gene in genes)
    if seed.kind == "protein_name":
        protein_name = mapping.get("protein_name")
        if not isinstance(protein_name, str):
            return False
        normalized = _normalize_candidate_text(protein_name)
        return target == normalized or target in normalized
    return True


def _candidate_seed_lineage(
    seed: CandidateSeed,
    records: Sequence[RawResearchEvidence],
) -> tuple[list[RawResearchEvidence], set[str]]:
    related: list[RawResearchEvidence] = []
    lineage: set[str] = set()
    linked_genes: set[str] = set()

    if seed.kind == "accession":
        lineage.add(seed.value.upper())
    for record in records:
        if seed.kind == "accession":
            requested = record.query_arguments.get("accession")
            if (
                record.tool_name == "UniProt_get_entry_by_accession"
                and isinstance(requested, str)
                and requested.upper() == seed.value.upper()
            ):
                related.append(record)
        elif seed.kind in {"gene", "locus", "protein_name"}:
            if not _record_matches_seed_search(record, seed):
                continue
            related.append(record)
            if record.status == "success":
                for mapping in _extract_search_candidate_mappings(record.content):
                    if _search_mapping_matches_seed(mapping, seed):
                        lineage.add(str(mapping["accession"]).upper())
        elif seed.kind == "ko" and record.tool_name == "KEGG_link_entries":
            source = str(record.query_arguments.get("source") or "")
            target = str(record.query_arguments.get("target") or "")
            if source.casefold() != f"ko:{seed.value}".casefold():
                continue
            if (
                seed.kegg_organism_code is not None
                and target.casefold() != seed.kegg_organism_code.casefold()
            ):
                continue
            related.append(record)
            if record.status == "success":
                linked_genes.update(_extract_kegg_gene_ids(record.content))

    if seed.kind == "ko":
        for record in records:
            if record.tool_name != "KEGG_convert_ids":
                continue
            gene_id = str(record.query_arguments.get("kegg_id") or "")
            if gene_id not in linked_genes:
                continue
            related.append(record)
            if record.status == "success":
                lineage.update(_extract_converted_accessions(record.content))

    for record in records:
        if record.tool_name != "UniProt_get_entry_by_accession":
            continue
        requested = record.query_arguments.get("accession")
        if isinstance(requested, str) and requested.upper() in lineage:
            if record not in related:
                related.append(record)
    return related, lineage


def _is_explicit_not_found(message: str) -> bool:
    return bool(
        re.search(
            r"(?:\bnot found\b|\bno results?(?: found)?\b|"
            r"\bno such (?:entry|record|identifier)\b|"
            r"\bunknown (?:accession|identifier)\b|\b404\b)",
            message,
            re.IGNORECASE,
        )
    )


def _verified_hint_from_exact_record(
    record: RawResearchEvidence,
    seed: CandidateSeed,
    allowed_accessions: set[str],
) -> VerifiedCandidateHint | None:
    requested = record.query_arguments.get("accession")
    if not isinstance(requested, str) or requested.upper() not in allowed_accessions:
        return None
    mappings = _extract_search_candidate_mappings(record.content)
    mapping = next(
        (
            item
            for item in mappings
            if str(item.get("accession") or "").upper() == requested.upper()
        ),
        None,
    )
    if mapping is None:
        return None
    taxon_id = mapping.get("taxon_id")
    if seed.taxon_id is not None and taxon_id != seed.taxon_id:
        return None
    genes = tuple(mapping.get("gene_names") or ())
    protein_name = mapping.get("protein_name")
    aliases = list(genes)
    if isinstance(protein_name, str):
        aliases.append(protein_name)
    if seed.kind in {"gene", "locus", "protein_name"}:
        aliases.append(seed.value)
    aliases.extend(_gene_cluster_aliases(seed.source_span))
    if seed.role_hint:
        aliases.extend(_gene_cluster_aliases(seed.role_hint))
    reference_pmids = _normalize_pmids(
        [
            *_REFERENCE_PMID_PATTERN.findall(seed.source_span),
            *(
                _REFERENCE_PMID_PATTERN.findall(seed.role_hint)
                if seed.role_hint
                else []
            ),
            *extract_pmids(record.content),
        ]
    )
    return VerifiedCandidateHint(
        uniprot_id=requested.upper(),
        protein_name=protein_name if isinstance(protein_name, str) else None,
        gene_names=genes,
        organism_name=(
            mapping.get("organism_name")
            if isinstance(mapping.get("organism_name"), str)
            else None
        ),
        taxon_id=taxon_id if isinstance(taxon_id, int) else None,
        aliases=_unique_strings(aliases),
        ko_ids=(seed.value.upper(),) if seed.kind == "ko" else (),
        seed_evidence_ids=(seed.source_evidence_id,),
        reference_pmids=reference_pmids,
    )


def _content_is_explicit_empty(tool_name: str, content: str) -> bool:
    payload = _json_value(content)
    if payload is None:
        return False
    if isinstance(payload, list):
        return not payload
    if not isinstance(payload, Mapping):
        return False
    if payload.get("error"):
        return False
    for key in ("total_results", "totalResults", "returned", "count"):
        value = payload.get(key)
        if isinstance(value, (int, float)) and value == 0:
            return True
    for key in ("results", "entries", "items"):
        value = payload.get(key)
        if isinstance(value, list):
            return not value
    if "data" in payload:
        data = payload["data"]
        if isinstance(data, list):
            return not data
        if isinstance(data, Mapping):
            serialized = json.dumps(data, ensure_ascii=False, default=str)
            return _content_is_explicit_empty(tool_name, serialized)
        return data in (None, "")
    return False


def _has_verified_empty_response(
    seed: CandidateSeed,
    successful_records: Sequence[RawResearchEvidence],
) -> bool:
    relevant_tools = {
        "accession": {"UniProt_get_entry_by_accession"},
        "gene": {"UniProt_search", "proteins_api_search"},
        "locus": {"UniProt_search", "proteins_api_search"},
        "protein_name": {"UniProt_search", "proteins_api_search"},
        "ko": {"KEGG_link_entries", "KEGG_convert_ids"},
    }[seed.kind]
    relevant = [
        record for record in successful_records if record.tool_name in relevant_tools
    ]
    return bool(relevant) and all(
        _content_is_explicit_empty(record.tool_name, record.content)
        for record in relevant
    )


def _merge_candidate_hints(hints: Any) -> dict[str, VerifiedCandidateHint]:
    merged: dict[str, VerifiedCandidateHint] = {}
    for hint in hints:
        key = hint.uniprot_id.upper()
        existing = merged.get(key)
        if existing is None:
            merged[key] = hint
            continue
        merged[key] = VerifiedCandidateHint(
            uniprot_id=key,
            protein_name=existing.protein_name or hint.protein_name,
            gene_names=_unique_strings((*existing.gene_names, *hint.gene_names)),
            organism_name=existing.organism_name or hint.organism_name,
            taxon_id=(
                existing.taxon_id
                if existing.taxon_id is not None
                else hint.taxon_id
            ),
            aliases=_unique_strings((*existing.aliases, *hint.aliases)),
            ko_ids=_unique_strings((*existing.ko_ids, *hint.ko_ids)),
            seed_evidence_ids=_unique_strings(
                (*existing.seed_evidence_ids, *hint.seed_evidence_ids)
            ),
            reference_pmids=_unique_strings(
                (*existing.reference_pmids, *hint.reference_pmids)
            ),
        )
    return merged


def _unique_strings(values: Sequence[str] | Any) -> tuple[str, ...]:
    unique: list[str] = []
    seen: set[str] = set()
    for value in values:
        if not isinstance(value, str):
            continue
        normalized = " ".join(value.split())
        key = normalized.casefold()
        if not normalized or key in seen:
            continue
        seen.add(key)
        unique.append(normalized)
    return tuple(unique)


def _normalize_pmids(values: Sequence[Any]) -> tuple[str, ...]:
    pmids: list[str] = []
    for value in values:
        text = str(value)
        match = re.fullmatch(r"(?:PMID\s*[:#-]?\s*)?(\d{5,9})", text, re.I)
        if match is None:
            match = _REFERENCE_PMID_PATTERN.search(text)
        if match is not None and match.group(1) not in pmids:
            pmids.append(match.group(1))
    return tuple(pmids)


def _context_reference_pmids(context: Any) -> tuple[str, ...]:
    values: list[Any] = []
    for owner in (
        context,
        getattr(context, "protein", None),
        getattr(context, "reaction", None),
    ):
        if owner is None:
            continue
        for field_name in (
            "reference_pmids",
            "pmids",
            "publication_pmids",
            "curated_pmids",
        ):
            raw = getattr(owner, field_name, ())
            if isinstance(raw, str):
                values.append(raw)
            elif isinstance(raw, Sequence):
                values.extend(raw)
    protein = getattr(context, "protein", None)
    for crossref in getattr(protein, "cross_references", ()):
        database = getattr(crossref, "database", None)
        record_id = getattr(crossref, "record_id", None)
        if isinstance(crossref, Mapping):
            database = database or crossref.get("database")
            record_id = record_id or crossref.get("record_id") or crossref.get("id")
        if isinstance(database, str) and database.casefold() == "pubmed":
            values.append(record_id)
    return _normalize_pmids(values)


def _reaction_literature_identity(context: Any) -> str:
    reaction = context.reaction
    terms: list[str] = [str(reaction.reaction_id)]
    names = getattr(reaction, "names", ())
    if names:
        terms.append(str(names[0]))
    ec_numbers = getattr(reaction, "ec_numbers", ())
    if ec_numbers:
        terms.append(str(ec_numbers[0]))
    for item in getattr(reaction, "orthology", ())[:2]:
        value = getattr(item, "orthology_id", None)
        if value is None and isinstance(item, Mapping):
            value = item.get("orthology_id") or item.get("id")
        if isinstance(value, str):
            terms.append(value)
    return " OR ".join(_quote_term(value) or value for value in _unique_strings(terms))


def _organism_alias_query(
    value: str | None,
    extra_aliases: Sequence[str] = (),
) -> str | None:
    species = _species_search_name(value)
    if species is None:
        return None
    aliases = [species]
    aliases.extend(
        alias
        for raw_alias in extra_aliases
        if (alias := _species_search_name(raw_alias)) is not None
        and len(alias.split()) >= 2
        and not alias.casefold().startswith("strain ")
    )
    for alias in list(aliases):
        parts = alias.split()
        if len(parts) >= 2 and len(parts[0]) > 1:
            aliases.append(f"{parts[0][0]}. {' '.join(parts[1:])}")
    return " OR ".join(_quote_term(alias) or alias for alias in _unique_strings(aliases))


def _direct_pmid_query(query: str) -> str | None:
    match = re.fullmatch(
        r"\s*(?:PMID\s*[:#-]?\s*)?(\d{5,9})(?:\s*\[PMID\])?\s*",
        query,
        re.IGNORECASE,
    )
    return match.group(1) if match is not None else None


def _first_search_term(context: Any) -> str | None:
    for values in (context.gene_names, context.protein_names):
        if values:
            return values[0]
    return None


def _mechanism_terms(role: str) -> str:
    normalized = role.casefold()
    if any(
        marker in normalized
        for marker in (
            "electron",
            "redox",
            "reductase",
            "ferredoxin",
            "flavodoxin",
        )
    ):
        return (
            "electron transfer OR NADH OR reductase OR ferredoxin OR "
            "flavodoxin OR reconstitution"
        )
    if "matur" in normalized or "assembl" in normalized:
        return "maturation OR assembly OR activation OR complementation"
    return "reconstitution OR residual activity OR purified OR knockout"


def _quote_term(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = " ".join(value.replace('"', " ").split())
    return f'"{normalized}"' if normalized else None


def _species_search_name(value: str | None) -> str | None:
    """Return a literature-friendly species name without strain qualifiers."""

    if value is None:
        return None
    species = re.sub(
        r"\s*\([^)]*(?:strain|substrain)[^)]*\)\s*",
        " ",
        value,
        flags=re.IGNORECASE,
    )
    normalized = " ".join(species.split())
    return normalized or None


def _deduplicate_calls(calls: Sequence[PlannedToolCall]) -> list[PlannedToolCall]:
    unique: list[PlannedToolCall] = []
    seen: set[str] = set()
    for call in calls:
        key = _cache_key(call.tool_name, call.arguments)
        if key in seen:
            continue
        seen.add(key)
        unique.append(call)
    return unique


def _cache_key(tool_name: str, arguments: Mapping[str, Any]) -> str:
    serialized = json.dumps(
        arguments,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return f"{tool_name}:{serialized}"


def _tool_source(tool_name: str) -> str:
    if tool_name.startswith("PubMed_"):
        return "pubmed"
    if tool_name.startswith("EuropePMC_"):
        return "europe_pmc"
    if tool_name.startswith("SemanticScholar_"):
        return "semantic_scholar"
    if tool_name in {"search", "fetchWebContent"}:
        return "web"
    return "database"


_TRANSIENT_TOOL_ERROR_PATTERN = re.compile(
    r"(?:ssl|unexpected_eof|eof occurred|connection|temporar|rate.?limit|"
    r"too many requests|\b429\b|\b500\b|\b502\b|\b503\b|\b504\b|"
    r"timeout|timed out|server error|service unavailable)",
    re.IGNORECASE,
)


def _is_transient_tool_error(message: str) -> bool:
    normalized = message.strip()
    if not normalized:
        return False
    if re.search(
        r"(?:not found|invalid (?:id|identifier|request)|no results?)",
        normalized,
        re.IGNORECASE,
    ):
        return False
    return _TRANSIENT_TOOL_ERROR_PATTERN.search(normalized) is not None


def _evidence_id(call: PlannedToolCall, digest: str) -> str:
    seed = f"{call.stage}:{_cache_key(call.tool_name, call.arguments)}:{digest}"
    return "RAW-" + hashlib.sha256(seed.encode("utf-8")).hexdigest()[:16]


def _serialize_tool_result(raw: Any) -> str:
    text_payload = _mcp_text_payload(raw)
    if text_payload is not None:
        parsed = _json_value(text_payload)
        if parsed is not None:
            return json.dumps(
                parsed,
                ensure_ascii=False,
                sort_keys=True,
                default=str,
            )
        return text_payload
    if isinstance(raw, str):
        return raw
    if hasattr(raw, "model_dump"):
        raw = raw.model_dump(mode="json")
    try:
        return json.dumps(raw, ensure_ascii=False, sort_keys=True, default=str)
    except (TypeError, ValueError):
        return str(raw)


def _compact_known_tool_content(tool_name: str, content: str) -> str:
    """Keep only decision-relevant fields from verbose official responses."""

    if tool_name != "UniProt_get_entry_by_accession":
        return content
    payload = _json_object(content)
    if payload is None:
        return content
    data = payload.get("data")
    if not isinstance(data, Mapping):
        return content
    compact_data = {
        key: data[key]
        for key in (
            "entryType",
            "primaryAccession",
            "uniProtkbId",
            "protein_name",
            "gene_names",
            "organism",
            "proteinExistence",
            "comments",
        )
        if key in data
    }
    compact_payload = {
        key: payload[key]
        for key in ("status", "metadata")
        if key in payload
    }
    compact_payload["data"] = compact_data
    return json.dumps(
        compact_payload,
        ensure_ascii=False,
        sort_keys=True,
        default=str,
    )


def _mcp_text_payload(raw: Any) -> str | None:
    """Unwrap MCP text content blocks before indexing and grounding.

    ``langchain-mcp-adapters`` commonly returns ``[{type: text, text: ...}]``.
    Serializing that wrapper would escape every quote in the actual JSON,
    preventing PMID/URL extraction and causing exact excerpts to fail their
    grounding check.
    """

    blocks: Any = None
    if isinstance(raw, list):
        blocks = raw
    elif isinstance(raw, Mapping) and isinstance(raw.get("content"), list):
        blocks = raw["content"]
    else:
        candidate = getattr(raw, "content", None)
        if isinstance(candidate, list):
            blocks = candidate
    if not isinstance(blocks, list) or not blocks:
        return None

    texts: list[str] = []
    for block in blocks:
        if isinstance(block, Mapping):
            block_type = block.get("type")
            block_text = block.get("text")
        else:
            block_type = getattr(block, "type", None)
            block_text = getattr(block, "text", None)
        if block_type != "text" or not isinstance(block_text, str):
            return None
        texts.append(block_text)
    return "\n".join(texts)


def _json_value(content: str) -> Any | None:
    try:
        return json.loads(content)
    except (TypeError, ValueError):
        return None


def _json_object(content: str) -> Mapping[str, Any] | None:
    value = _json_value(content)
    return value if isinstance(value, Mapping) else None


def _tool_result_error(
    raw: Any,
    content: str,
    *,
    tool_name: str | None = None,
) -> str | None:
    """Recognize tool-level failures returned as ordinary values.

    MCP adapters do not always raise Python exceptions.  Some instead return a
    ToolMessage/CallToolResult with an error status, or a short error string.
    Such values must never enter the ledger as successful biological evidence.
    """

    if (
        tool_name == "KEGG_link_entries"
        and re.search(
            r"\bno\s+[a-z0-9_]+\s+entries\s+linked\s+to\s+ko:K\d{5}\b",
            content,
            re.IGNORECASE,
        )
    ):
        # ToolUniverse 1.4.1 reports a genuine empty KEGG mapping with an
        # error-shaped envelope.  It is a successful, auditable empty result,
        # unlike connection/SSL/server failures.
        return None

    status = getattr(raw, "status", None)
    if isinstance(status, str) and status.casefold() in {
        "error",
        "failed",
        "failure",
    }:
        return content or f"tool returned status {status}"

    is_error = getattr(raw, "isError", None)
    if is_error is None:
        is_error = getattr(raw, "is_error", None)
    if is_error is True:
        return content or "tool returned an MCP error result"

    raw_mapping: Mapping[str, Any] | None = None
    if isinstance(raw, Mapping):
        raw_mapping = raw
    elif hasattr(raw, "model_dump"):
        dumped = raw.model_dump(mode="json")
        if isinstance(dumped, Mapping):
            raw_mapping = dumped
    if raw_mapping is not None:
        mapping_status = raw_mapping.get("status")
        if (
            isinstance(mapping_status, str)
            and mapping_status.casefold() in {"error", "failed", "failure"}
        ) or raw_mapping.get("isError") is True or raw_mapping.get(
            "is_error"
        ) is True:
            return content or "tool returned an error result"

    content_mapping = _json_object(content)
    if content_mapping is not None:
        content_status = content_mapping.get("status")
        if isinstance(content_status, str) and content_status.casefold() in {
            "error",
            "failed",
            "failure",
        }:
            message = content_mapping.get("message") or content_mapping.get(
                "error"
            )
            return str(message) if message else content

    if re.match(
        r"^\s*(?:error\b|failed\b|failure\b|timeout\b|timed out\b|"
        r"tool\s+.+?\s+timed out\b)",
        content,
        re.IGNORECASE,
    ):
        return content
    return None


def _extract_urls(content: str) -> list[str]:
    return list(dict.fromkeys(_URL_PATTERN.findall(content)))


def _extract_source_record_ids(
    tool_name: str,
    arguments: Mapping[str, Any],
    content: str,
) -> list[str]:
    ids: list[str] = []
    for key in (
        "accession",
        "complex_id",
        "identifier",
        "pmid",
        "rhea_id",
    ):
        value = arguments.get(key)
        if isinstance(value, (str, int)):
            ids.append(str(value))
    if "PubMed" in tool_name:
        ids.extend(extract_pmids(content))
    return list(dict.fromkeys(ids))


def _extract_protein_sequence(content: str) -> str | None:
    patterns = (
        re.compile(r'"sequence"\s*:\s*"([A-Z]{10,})"', re.IGNORECASE),
        re.compile(r"(?<![A-Z])[ACDEFGHIKLMNPQRSTVWY]{30,}(?![A-Z])"),
    )
    for pattern in patterns:
        match = pattern.search(content)
        if match is not None:
            sequence = match.group(1) if match.lastindex else match.group(0)
            return sequence.upper()
    return None
