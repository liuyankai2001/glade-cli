from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from pathlib import Path
from typing import Any

import requests

from src.main_protein_selection.settings import (
    GLADE_CONTACT_EMAIL_ENV,
    RHEA_HTTP_CONFIG,
    RHEA_REST_URL as CONFIGURED_RHEA_REST_URL,
    env_value,
)
from src.main_protein_selection.uniprot_protein_candidates import (
    ProteinCandidate,
    _merge_entries,
    candidate_from_reaction_entry,
    extract_rhea_ids,
    search_uniprot_by_query,
)


# 兼容原有模块字段，便于测试按模块临时替换。
RHEA_REST_URL = CONFIGURED_RHEA_REST_URL
RHEA_TIMEOUT_SECONDS = RHEA_HTTP_CONFIG.timeout_seconds
RHEA_RETRIES = RHEA_HTTP_CONFIG.retries


def _safe_token(value: str) -> str:
    return "".join(char if char.isalnum() or char in "._-" else "_" for char in value)


def _json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _user_agent() -> str:
    contact = env_value(GLADE_CONTACT_EMAIL_ENV)
    return f"GLADE/2.2 ({contact})" if contact else "GLADE/2.2"


class RheaClient:
    def __init__(
        self,
        session: requests.Session | None = None,
        cache_root: Path | None = None,
    ) -> None:
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", _user_agent())
        if cache_root is None:
            raise ValueError("cache_root is required")
        self.cache_root = Path(cache_root)

    def reactions_for_kegg(self, reaction_id: str) -> dict[str, Any]:
        return self.reactions_for_kegg_many([reaction_id]).get(
            str(reaction_id or "").strip().upper(),
            {},
        )

    def reactions_for_kegg_many(
        self,
        reaction_ids: list[str],
    ) -> dict[str, dict[str, Any]]:
        normalized_ids = list(dict.fromkeys(
            str(value or "").strip().upper()
            for value in reaction_ids
            if str(value or "").strip()
        ))
        results: dict[str, dict[str, Any]] = {}
        missing: list[str] = []
        for normalized in normalized_ids:
            query_id = f"rhea_kegg_{normalized}"
            cache_path = self.cache_root / f"{_safe_token(query_id)}.json"
            if cache_path.exists():
                cached = json.loads(cache_path.read_text(encoding="utf-8"))
                if isinstance(cached, dict):
                    cached["cache_hit"] = True
                    results[normalized] = cached
                    continue
            missing.append(normalized)
        if not missing:
            return results

        batch_query_id = "rhea_kegg_batch_" + hashlib.sha256(
            ";".join(missing).encode("utf-8")
        ).hexdigest()[:16]
        params = {
            "query": " OR ".join(missing),
            "columns": "rhea-id,equation,ec,reaction-xref(KEGG),uniprot",
            "format": "tsv",
            "limit": 100,
        }
        last_error: Exception | None = None
        response: requests.Response | None = None
        for attempt in range(RHEA_RETRIES):
            try:
                response = self.session.get(
                    RHEA_REST_URL,
                    params=params,
                    timeout=RHEA_TIMEOUT_SECONDS,
                )
                response.raise_for_status()
                break
            except requests.RequestException as exc:
                last_error = exc
                if attempt + 1 < RHEA_RETRIES:
                    time.sleep(RHEA_HTTP_CONFIG.sleep_seconds * (2**attempt))
        if response is None:
            for normalized in missing:
                results[normalized] = {
                    "query_id": f"rhea_kegg_{normalized}",
                    "batch_query_id": batch_query_id,
                    "query_type": "rhea_by_kegg_reaction",
                    "reaction_id": normalized,
                    "status": "source_unavailable",
                    "records": [],
                    "error": str(last_error or "Rhea request failed"),
                    "cache_hit": False,
                }
            return results

        raw_text = response.text
        records: list[dict[str, Any]] = []
        reader = csv.DictReader(io.StringIO(raw_text), delimiter="\t")
        for row in reader:
            raw_rhea_id = str(row.get("Reaction identifier") or "").strip()
            if not raw_rhea_id:
                continue
            records.append({
                "rhea_id": raw_rhea_id,
                "equation": str(row.get("Equation") or "").strip(),
                "ec_numbers": [
                    item.strip().removeprefix("EC:")
                    for item in str(row.get("EC number") or "").split(";")
                    if item.strip()
                ],
                "kegg_reactions": [
                    item.strip().removeprefix("KEGG:")
                    for item in str(row.get("Cross-reference (KEGG)") or "").split(";")
                    if item.strip()
                ],
                "uniprot_count": int(str(row.get("Enzymes") or "0").strip() or 0),
            })
        for normalized in missing:
            matched_records = [
                record
                for record in records
                if normalized in set(record.get("kegg_reactions", []))
            ]
            payload = {
                "query_id": f"rhea_kegg_{normalized}",
                "batch_query_id": batch_query_id,
                "query_type": "rhea_by_kegg_reaction",
                "reaction_id": normalized,
                "status": "ok",
                "records": matched_records,
                "request": {"url": RHEA_REST_URL, "params": params},
                "response_sha256": hashlib.sha256(raw_text.encode("utf-8")).hexdigest(),
                "raw_tsv": raw_text,
                "cache_hit": False,
            }
            cache_path = self.cache_root / f"{_safe_token(payload['query_id'])}.json"
            _json_atomic(cache_path, payload)
            results[normalized] = payload
        return results


def reaction_evidence_for_requirements(
    requirements: list[dict[str, Any]],
    client: RheaClient | None = None,
) -> list[dict[str, Any]]:
    rhea = client or RheaClient()
    evidence: list[dict[str, Any]] = []
    reaction_ids = list(dict.fromkeys(
        str(requirement.get("reaction_id") or "").strip().upper()
        for requirement in requirements
        if str(requirement.get("reaction_id") or "").strip()
    ))
    by_reaction = rhea.reactions_for_kegg_many(reaction_ids)
    for requirement in requirements:
        reaction_id = str(requirement.get("reaction_id") or "").strip().upper()
        source = by_reaction.get(reaction_id, {
            "query_id": f"rhea_kegg_{reaction_id}",
            "status": "source_unavailable",
            "records": [],
            "error": "Reaction identifier was not queried",
        })
        rhea_master_ids = [
            str(record.get("rhea_id") or "")
            for record in source.get("records", [])
            if str(record.get("rhea_id") or "")
        ]
        # Keep the directional/bidirectional IDs supplied by the locked KEGG
        # route.  The Rhea REST result is normally the undirected master and
        # must not silently replace that more specific route evidence.
        kegg_rhea_ids = [
            str(value) for value in requirement.get("rhea_ids", []) if str(value)
        ]
        requirement["kegg_rhea_ids"] = kegg_rhea_ids
        requirement["rhea_master_ids"] = rhea_master_ids
        requirement["rhea_master_equations"] = [
            str(record.get("equation") or "")
            for record in source.get("records", [])
            if str(record.get("equation") or "").strip()
        ]
        requirement["reaction_evidence_query_id"] = source.get("query_id", "")
        evidence.append({
            "step_index": int(requirement.get("step_index") or 0),
            "reaction_id": reaction_id,
            "direction": str(requirement.get("direction") or ""),
            "ec_status": str(requirement.get("ec_status") or ""),
            "kegg_rhea_ids": kegg_rhea_ids,
            "rhea_master_ids": rhea_master_ids,
            "source_status": source.get("status", ""),
            "query_id": source.get("query_id", ""),
            "records": source.get("records", []),
            "error": source.get("error", ""),
        })
    return evidence


def _rhea_candidate_entries(
    rhea_ids: list[str],
    session: requests.Session,
    max_results: int,
) -> tuple[list[dict[str, Any]], list[str], dict[str, str]]:
    groups: list[list[dict[str, Any]]] = []
    query_ids: list[str] = []
    errors: dict[str, str] = {}
    for raw_id in rhea_ids:
        normalized = str(raw_id).upper().removeprefix("RHEA:")
        query_id = f"uniprot_rhea_{normalized}"
        query = f'cc_catalytic_activity:"rhea:{normalized}"'
        query_ids.append(query_id)
        try:
            reviewed = search_uniprot_by_query(
                query,
                reviewed_only=True,
                max_results=min(max_results, 500),
                session=session,
            )
            all_entries = search_uniprot_by_query(
                query,
                reviewed_only=False,
                max_results=max_results,
                session=session,
            )
            groups.extend([reviewed, all_entries])
        except Exception as exc:
            errors[query_id] = str(exc)
    return _merge_entries(groups), query_ids, errors


def retrieve_rhea_candidates_for_requirement(
    requirement: dict[str, Any],
    chassis_key: str,
    *,
    max_results: int,
    top_n: int | None = None,
    allow_transmembrane: bool,
    session: requests.Session,
) -> tuple[list[ProteinCandidate], list[str], dict[str, str]]:
    rhea_ids = [
        str(value)
        for value in (
            requirement.get("rhea_retrieval_ids")
            or requirement.get("rhea_ids", [])
        )
        if str(value)
    ]
    if not rhea_ids:
        return [], [], {}
    entries, query_ids, errors = _rhea_candidate_entries(rhea_ids, session, max_results)
    candidates: list[ProteinCandidate] = []
    for entry in entries:
        candidate_rhea_ids = {
            str(value).upper()
            for value in extract_rhea_ids(entry)
        }
        matched = [rhea_id for rhea_id in rhea_ids if rhea_id.upper() in candidate_rhea_ids]
        if not matched:
            continue
        candidate = candidate_from_reaction_entry(
            entry,
            chassis_key,
            retrieval_strategy="rhea_reaction",
            retrieval_query_id=";".join(query_ids),
            matched_rhea_ids=matched,
            allow_transmembrane=allow_transmembrane,
        )
        if candidate is not None:
            candidates.append(candidate)
    candidates.sort(key=lambda row: row.score, reverse=True)
    return candidates[:top_n] if top_n is not None else candidates, query_ids, errors


__all__ = [
    "RheaClient",
    "reaction_evidence_for_requirements",
    "retrieve_rhea_candidates_for_requirement",
]
