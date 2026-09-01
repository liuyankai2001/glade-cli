"""Exact enzyme retrieval through KEGG Orthology before SelenzymeRF fallback."""

from __future__ import annotations

import hashlib
import json
import tempfile
import time
from pathlib import Path
from typing import Any

import requests

from src.main_protein_selection.taxonomy_compatibility import ChassisTaxonomyProfile

from src.main_protein_selection.settings import KEGG_HTTP_CONFIG, KEGG_REST_BASE_URL
from src.main_protein_selection.uniprot_protein_candidates import (
    ProteinCandidate,
    UNIPROT_ACCESSION_BATCH_SIZE,
    candidate_from_reaction_entry,
    resolve_uniprot_accession_batches,
)


KO_TO_UNIPROT_BATCH_SIZE = 10


class KeggKOSourceUnavailable(RuntimeError):
    """Raised when a KO query cannot use either KEGG or a local cache."""


def _safe_token(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in value
    )


def _json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
        suffix=".tmp",
    ) as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def _chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


class KeggKOClient:
    def __init__(
        self,
        session: requests.Session | None = None,
        cache_root: Path | None = None,
        base_url: str | None = None,
    ) -> None:
        self.session = session or requests.Session()
        if cache_root is None:
            raise ValueError("cache_root is required")
        self.cache_root = Path(cache_root)
        self.base_url = (base_url or KEGG_REST_BASE_URL).rstrip("/")

    def _get_text(self, url: str) -> str:
        last_error: Exception | None = None
        for attempt in range(KEGG_HTTP_CONFIG.retries):
            try:
                response = self.session.get(
                    url,
                    timeout=KEGG_HTTP_CONFIG.timeout_seconds,
                )
                response.raise_for_status()
                return response.text
            except requests.RequestException as exc:
                last_error = exc
                if attempt + 1 < KEGG_HTTP_CONFIG.retries:
                    time.sleep(KEGG_HTTP_CONFIG.sleep_seconds * (2**attempt))
        raise KeggKOSourceUnavailable(
            f"KEGG KO request failed: {url}: {last_error or 'unknown error'}"
        )

    def proteins_for_ko(self, ko_id: str) -> dict[str, Any]:
        normalized = str(ko_id or "").strip().upper().removeprefix("KO:")
        if not normalized.startswith("K") or not normalized[1:].isdigit():
            raise ValueError(f"Invalid KEGG KO identifier: {ko_id!r}")

        query_url = f"{self.base_url}/link/genes/ko:{normalized}"
        cache_key = hashlib.sha256(
            json.dumps(
                {
                    "cache_schema": "kegg_ko_client.v1",
                    "query_url": query_url,
                },
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest()[:16]
        query_id = f"kegg_ko_{normalized}_{cache_key}"
        cache_path = self.cache_root / f"{_safe_token(query_id)}.json"
        if cache_path.exists():
            try:
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                cached = None
            if (
                isinstance(cached, dict)
                and cached.get("cache_schema") == "kegg_ko_client.v1"
                and cached.get("ko_id") == normalized
                and cached.get("query_type") == "uniprot_by_kegg_ko"
            ):
                return {**cached, "cache_hit": True}

        gene_text = self._get_text(query_url)
        gene_ids: list[str] = []
        for line in gene_text.splitlines():
            if "\t" not in line:
                continue
            raw_ko, raw_gene = line.split("\t", 1)
            if raw_ko.strip().upper().removeprefix("KO:") != normalized:
                continue
            gene_id = raw_gene.strip()
            if gene_id and gene_id not in gene_ids:
                gene_ids.append(gene_id)

        mappings: list[dict[str, str]] = []
        allowed_gene_ids = set(gene_ids)
        conversion_urls: list[str] = []
        conversion_texts: list[str] = []
        for batch in _chunks(gene_ids, KO_TO_UNIPROT_BATCH_SIZE):
            conversion_url = f"{self.base_url}/conv/uniprot/{'+'.join(batch)}"
            conversion_urls.append(conversion_url)
            conversion_text = self._get_text(conversion_url)
            conversion_texts.append(conversion_text)
            for line in conversion_text.splitlines():
                if "\t" not in line:
                    continue
                raw_gene, raw_accession = line.split("\t", 1)
                gene_id = raw_gene.strip()
                accession = raw_accession.strip().removeprefix("up:")
                if gene_id in allowed_gene_ids and accession:
                    mappings.append({
                        "gene_id": gene_id,
                        "accession": accession,
                    })

        unique_mappings: list[dict[str, str]] = []
        seen_mappings: set[tuple[str, str]] = set()
        for mapping in mappings:
            key = (mapping["gene_id"], mapping["accession"])
            if key not in seen_mappings:
                seen_mappings.add(key)
                unique_mappings.append(mapping)

        status = "ok" if unique_mappings else (
            "no_uniprot_mapping" if gene_ids else "no_hit"
        )
        raw_evidence = "\n".join([gene_text, *conversion_texts])
        result = {
            "cache_schema": "kegg_ko_client.v1",
            "query_id": query_id,
            "query_type": "uniprot_by_kegg_ko",
            "ko_id": normalized,
            "status": status,
            "gene_ids": gene_ids,
            "mappings": unique_mappings,
            "request": {
                "gene_link_url": query_url,
                "conversion_urls": conversion_urls,
            },
            "response_sha256": hashlib.sha256(
                raw_evidence.encode("utf-8")
            ).hexdigest(),
            "cache_hit": False,
        }
        _json_atomic(cache_path, result)
        return result


def retrieve_ko_candidates(
    requirement: dict[str, Any],
    query_result: dict[str, Any],
    chassis_key: str,
    *,
    top_n: int,
    max_results: int,
    allow_transmembrane: bool,
    session: requests.Session,
    entry_cache: dict[str, dict[str, Any] | None] | None = None,
    taxonomy_profile: ChassisTaxonomyProfile | None = None,
) -> tuple[list[ProteinCandidate], list[dict[str, Any]], list[str], dict[str, str]]:
    """Resolve an exact KEGG KO mapping to filtered UniProt candidates."""

    entry_cache = entry_cache if entry_cache is not None else {}
    ko_id = str(query_result.get("ko_id") or "").strip().upper()
    requirement_ko_ids = {
        str(value or "").strip().upper().removeprefix("KO:")
        for value in requirement.get("ko_ids", [])
        if str(value or "").strip()
    }
    if ko_id not in requirement_ko_ids:
        raise ValueError(
            f"KO query {ko_id or '<missing>'} does not match requirement KO IDs"
        )
    mappings = [
        mapping
        for mapping in query_result.get("mappings", [])
        if isinstance(mapping, dict)
        and str(mapping.get("accession") or "").strip()
    ]
    genes_by_accession: dict[str, list[str]] = {}
    for mapping in mappings:
        accession = str(mapping.get("accession") or "").strip().upper()
        gene_id = str(mapping.get("gene_id") or "").strip()
        genes_by_accession.setdefault(accession, [])
        if gene_id and gene_id not in genes_by_accession[accession]:
            genes_by_accession[accession].append(gene_id)

    accessions = list(genes_by_accession)[: max(0, int(max_results))]
    query_ids: list[str] = []
    errors: dict[str, str] = {}
    resolution = resolve_uniprot_accession_batches(
        accessions,
        session=session,
        entry_cache=entry_cache,
        query_id_prefix=f"uniprot_kegg_ko_{ko_id}",
        batch_size=UNIPROT_ACCESSION_BATCH_SIZE,
    )
    query_ids.extend(resolution.query_ids)
    errors.update(resolution.query_errors)

    candidates: list[ProteinCandidate] = []
    audit_rows: list[dict[str, Any]] = []
    for accession in accessions:
        entry = entry_cache.get(accession)
        if not isinstance(entry, dict):
            audit_rows.append({
                "ko_id": ko_id,
                "accession": accession,
                "kegg_gene_ids": genes_by_accession.get(accession, []),
                "status": "uniprot_unresolved",
            })
            continue
        candidate = candidate_from_reaction_entry(
            entry,
            chassis_key,
            retrieval_strategy="kegg_ko_exact",
            retrieval_query_id=str(query_result.get("query_id") or ""),
            allow_transmembrane=allow_transmembrane,
            function_evidence_reason=f"function: exact KEGG KO match {ko_id}",
            taxonomy_profile=taxonomy_profile,
        )
        if candidate is None:
            audit_rows.append({
                "ko_id": ko_id,
                "accession": accession,
                "kegg_gene_ids": genes_by_accession.get(accession, []),
                "status": "sequence_filter_rejected",
            })
            continue
        candidate.matched_ko_ids = [ko_id]
        candidate.kegg_gene_ids = genes_by_accession.get(accession, [])
        candidate.reaction_confidence = "ko_exact"
        candidates.append(candidate)
        audit_rows.append({
            "ko_id": ko_id,
            "accession": accession,
            "kegg_gene_ids": candidate.kegg_gene_ids,
            "status": "candidate",
            "score": candidate.score,
        })

    candidates.sort(key=lambda candidate: candidate.score, reverse=True)
    return candidates[:top_n], audit_rows, query_ids, errors


__all__ = [
    "KeggKOClient",
    "KeggKOSourceUnavailable",
    "retrieve_ko_candidates",
]
