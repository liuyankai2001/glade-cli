from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import requests
from .rhea_http import new_rhea_session, request_response, reject_html

from .transmembrane import validate_options

from src.main_protein_selection.reaction_direction_verifier import normalize_rhea_ids
from src.main_protein_selection.taxonomy_compatibility import (
    ChassisTaxonomyProfile,
    normalize_scoring_weights,
)

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


def _parse_rhea_tsv(text: str) -> list[dict[str, Any]]:
    reject_html(text)
    reader = csv.DictReader(io.StringIO(text.lstrip("\ufeff")), delimiter="\t")
    required = {"Reaction identifier", "Equation", "EC number", "Cross-reference (KEGG)", "Enzymes"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames) or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("Rhea TSV is missing required columns")
    records = []
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("Rhea TSV row has an invalid number of fields")
        identifier = row["Reaction identifier"].strip()
        if not re.fullmatch(r"(?:RHEA:)?\d+", identifier) or not row["Equation"].strip():
            raise ValueError("Invalid Rhea reaction row")
        enzyme_count = int(row["Enzymes"].strip() or "0")
        if enzyme_count < 0:
            raise ValueError("Invalid Rhea enzyme count")
        records.append({"rhea_id": identifier, "equation": row["Equation"].strip(),
            "ec_numbers": [v.strip().removeprefix("EC:") for v in row["EC number"].split(";") if v.strip()],
            "kegg_reactions": [v.strip().removeprefix("KEGG:") for v in row["Cross-reference (KEGG)"].split(";") if v.strip()],
            "uniprot_count": enzyme_count})
    return records


class RheaClient:
    def __init__(
        self,
        session: requests.Session | None = None,
        cache_root: Path | None = None,
    ) -> None:
        if cache_root is None:
            raise ValueError("cache_root is required")
        self._owns_session = session is None
        self.session = session if session is not None else new_rhea_session()
        self.cache_root = Path(cache_root)

    def close(self):
        if self._owns_session:
            self.session.close()
            self._owns_session = False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

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
                try:
                    cached = json.loads(cache_path.read_text(encoding="utf-8"))
                    if not isinstance(cached, dict) or cached.get("status") != "ok" or cached.get("reaction_id") != normalized:
                        raise ValueError("Rhea cache identity/status mismatch")
                    if cached.get("query_id") != query_id or cached.get("query_type") != "rhea_by_kegg_reaction":
                        raise ValueError("Rhea cache query mismatch")
                    request = cached.get("request", {})
                    if request.get("url", "").rstrip("/") != RHEA_REST_URL.rstrip("/") or not re.search(rf"\b{normalized}\b", request.get("params", {}).get("query", "")):
                        raise ValueError("Rhea cache request mismatch")
                    records = _parse_rhea_tsv(cached["raw_tsv"])
                    cached["records"] = [r for r in records if normalized in r["kegg_reactions"]]
                    cached.pop("response_sha256", None)
                    cached["cache_hit"] = True
                    results[normalized] = cached
                    continue
                except (OSError, ValueError, TypeError, KeyError, AttributeError):
                    # Keep the old file until a valid refresh can replace it.
                    pass
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
        try:
            response = request_response(self.session, RHEA_REST_URL, params=params,
                                        timeout=RHEA_TIMEOUT_SECONDS, retries=RHEA_RETRIES)
            raw_text = response.text
            records = _parse_rhea_tsv(raw_text)
        except (requests.RequestException, ValueError, TypeError) as exc:
            for normalized in missing:
                results[normalized] = {
                    "query_id": f"rhea_kegg_{normalized}",
                    "batch_query_id": batch_query_id,
                    "query_type": "rhea_by_kegg_reaction",
                    "reaction_id": normalized,
                    "status": "source_unavailable",
                    "records": [],
                    "error": str(exc),
                    "cache_hit": False,
                }
            return results

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
    for rhea_id in rhea_ids:
        query_id = f"uniprot_rhea_{rhea_id}"
        query = f'cc_catalytic_activity:"rhea:{rhea_id}"'
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
    taxonomy_profile: ChassisTaxonomyProfile | None = None,
    scoring_weights: Mapping[str, Any] | None = None,
    transmembrane_penalty: int = 35,
    no_transmembrane_bonus: int = 10,
    filter_audit: list[dict[str, Any]] | None = None,
) -> tuple[list[ProteinCandidate], list[str], dict[str, str]]:
    validate_options(allow_transmembrane, transmembrane_penalty, no_transmembrane_bonus)
    active_weights = normalize_scoring_weights(scoring_weights)
    rhea_ids = normalize_rhea_ids(
        requirement.get("rhea_retrieval_ids") or requirement.get("rhea_ids")
    )
    if not rhea_ids:
        return [], [], {}
    entries, query_ids, errors = _rhea_candidate_entries(rhea_ids, session, max_results)
    candidates: list[ProteinCandidate] = []
    for entry in entries:
        candidate_rhea_ids = set(normalize_rhea_ids(extract_rhea_ids(entry)))
        matched = [rhea_id for rhea_id in rhea_ids if rhea_id in candidate_rhea_ids]
        if not matched:
            continue
        candidate = candidate_from_reaction_entry(
            entry,
            chassis_key,
            retrieval_strategy="rhea_reaction",
            retrieval_query_id=";".join(query_ids),
            matched_rhea_ids=matched,
            allow_transmembrane=allow_transmembrane,
            transmembrane_penalty=transmembrane_penalty,
            no_transmembrane_bonus=no_transmembrane_bonus,
            filter_audit=filter_audit,
            taxonomy_profile=taxonomy_profile,
            scoring_weights=active_weights,
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
