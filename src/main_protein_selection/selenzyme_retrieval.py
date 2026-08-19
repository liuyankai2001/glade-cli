"""Reaction-first enzyme retrieval through the SelenzymeRF REST service."""

from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any

import requests

from src.main_protein_selection.settings import (
    SELENZYME_HTTP_CONFIG,
    get_selenzyme_rest_url,
)
from src.main_protein_selection.uniprot_protein_candidates import (
    CHASSIS_TAXON_PRESETS,
    ProteinCandidate,
    candidate_from_reaction_entry,
    hard_filter_candidate_without_ec,
    search_uniprot_by_query,
)


EXACT_SIMILARITY_TOLERANCE = 1e-6
COMPLETE_EC_PATTERN = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
EC_RELATION_EXACT = "exact_current_uniprot_ec"
EC_RELATION_SHARED_REACTION = "shared_reaction_ec_overlap"
EC_RELATION_UNANNOTATED = "unannotated_current_uniprot_ec"
EC_RELATION_CONTRADICTED = "contradicted_ec"


class SelenzymeSourceUnavailable(RuntimeError):
    """Raised when a required Selenzyme query cannot use network or cache."""


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


def _decode_selenzyme_rows(value: Any) -> list[dict[str, Any]]:
    """Decode Selenzyme's JSON-string, record, or pandas column orientation."""

    data = value
    if isinstance(data, str):
        stripped = data.strip()
        if not stripped:
            return []
        data = json.loads(stripped)
    if isinstance(data, list):
        return [dict(item) for item in data if isinstance(item, dict)]
    if not isinstance(data, dict) or not data:
        return []

    if all(isinstance(item, dict) for item in data.values()):
        keys = list(data)
        if all(str(key).isdigit() for key in keys):
            return [
                dict(data[key])
                for key in sorted(keys, key=lambda item: int(str(item)))
                if isinstance(data[key], dict)
            ]

        indexes: set[str] = set()
        for values in data.values():
            indexes.update(str(index) for index in values)
        rows: list[dict[str, Any]] = []
        for index in sorted(
            indexes,
            key=lambda item: (not item.isdigit(), int(item) if item.isdigit() else item),
        ):
            row = {
                str(column): values.get(index, values.get(int(index)) if index.isdigit() else None)
                for column, values in data.items()
                if isinstance(values, dict)
            }
            rows.append(row)
        return rows
    return [dict(data)]


def _first(row: dict[str, Any], *names: str) -> Any:
    lowered = {str(key).strip().lower(): value for key, value in row.items()}
    for name in names:
        if name in row:
            return row[name]
        match = lowered.get(name.strip().lower())
        if match is not None:
            return match
    return None


def _float_or_none(value: Any) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _ec_numbers(value: Any) -> set[str]:
    if isinstance(value, (list, tuple, set)):
        values = value
    else:
        values = str(value or "").replace("|", ";").split(";")
    return {
        str(item).strip()
        for item in values
        if COMPLETE_EC_PATTERN.fullmatch(str(item).strip())
    }


def classify_selenzyme_ec_relation(
    *,
    required_ecs: Any,
    candidate_ecs: Any,
    reported_ecs: Any,
) -> str:
    """Classify protein EC evidence without confusing reaction and protein annotations."""

    required = _ec_numbers(required_ecs)
    candidate = _ec_numbers(candidate_ecs)
    reported = _ec_numbers(reported_ecs)
    if candidate & required:
        return EC_RELATION_EXACT
    if not candidate:
        return EC_RELATION_UNANNOTATED
    if candidate & reported and required & reported:
        return EC_RELATION_SHARED_REACTION
    return EC_RELATION_CONTRADICTED


def _normalize_row(row: dict[str, Any], rank: int) -> dict[str, Any]:
    accession = str(
        _first(row, "Seq. ID", "Seq ID", "Sequence ID", "UniProt", "Entry") or ""
    ).strip()
    return {
        "rank": rank,
        "accession": accession,
        "score": _float_or_none(_first(row, "Score")),
        "description": str(_first(row, "Description") or "").strip(),
        "organism_source": str(_first(row, "Organism Source") or "").strip(),
        "matched_reaction_id": str(_first(row, "Rxn. ID", "Rxn ID") or "").strip(),
        "ec_number": str(_first(row, "EC Number", "EC") or "").strip(),
        "sim_rf": _float_or_none(
            _first(row, "sim_RF", "sim RF", "Rxn Sim RF.", "Rxn Sim RF")
        ),
        "sim_2018": _float_or_none(
            _first(row, "sim_2018", "sim 2018", "Rxn Sim.", "Rxn Sim")
        ),
        "reaction_similarity": _float_or_none(
            _first(
                row,
                "Reaction similarity",
                "Combined score",
                "Rxn Sim.",
                "Rxn Sim",
            )
        ),
        "taxonomic_distance": _float_or_none(
            _first(row, "Tax. distance", "Tax distance")
        ),
        "protein_evidence": str(
            _first(row, "Uniprot protein evidence", "UniProt protein evidence") or ""
        ).strip(),
        "direction_used": str(_first(row, "Direction Used") or "").strip(),
        "direction_preferred": str(_first(row, "Direction Preferred") or "").strip(),
    }


def selenzyme_match_type(row: dict[str, Any]) -> str:
    """Classify SelenzymeRF by its combined reaction-similarity score.

    ``sim_RF`` and ``sim_2018`` remain useful audit evidence, but the v2.6
    fallback policy deliberately uses only the combined score as its gate.
    """

    value = row.get("reaction_similarity")
    if not isinstance(value, (int, float)):
        return "invalid"
    similarity = float(value)
    if similarity < -EXACT_SIMILARITY_TOLERANCE or similarity > 1.0 + EXACT_SIMILARITY_TOLERANCE:
        return "invalid"
    if abs(similarity - 1.0) <= EXACT_SIMILARITY_TOLERANCE:
        return "exact"
    return "risk"


def selenzyme_target_count(top_n: int) -> int:
    return min(max(50, max(1, int(top_n)) * 10), 200)


def chassis_host_taxon_id(chassis_key: str) -> int:
    if chassis_key not in CHASSIS_TAXON_PRESETS:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    profile = CHASSIS_TAXON_PRESETS[chassis_key]
    return int(profile.get("strain_taxon_id") or profile["species_taxon_id"])


class SelenzymeClient:
    def __init__(
        self,
        session: requests.Session | None = None,
        cache_root: Path | None = None,
        rest_url: str | None = None,
    ) -> None:
        self.session = session or requests.Session()
        if cache_root is None:
            raise ValueError("cache_root is required")
        self.cache_root = Path(cache_root)
        self.rest_url = (rest_url or get_selenzyme_rest_url()).rstrip("/")
        self.query_url = (
            self.rest_url
            if self.rest_url.lower().endswith("/query")
            else f"{self.rest_url}/Query"
        )

    def _query_identifier(
        self,
        identifier: str,
        *,
        database: str,
        query_type: str,
        result_field: str,
        host_taxon_id: int,
        targets: int,
    ) -> dict[str, Any]:
        normalized = str(identifier or "").strip()
        if not normalized:
            raise ValueError(f"Selenzyme {database} query requires an identifier")
        payload = {
            "db": database,
            "rxnid": normalized,
            "targets": int(targets),
            "host": str(int(host_taxon_id)),
            "noMSA": True,
            "direction": 0,
        }
        cache_key = hashlib.sha256(
            json.dumps(
                {
                    "cache_schema": "selenzyme_client.v3",
                    "url": self.query_url,
                    "payload": payload,
                },
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest()[:20]
        query_id = f"selenzyme_{database}_{normalized}_{cache_key}"
        cache_path = self.cache_root / f"{_safe_token(query_id)}.json"
        if cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if isinstance(cached, dict):
                cached["cache_hit"] = True
                return cached

        last_error: Exception | None = None
        for attempt in range(SELENZYME_HTTP_CONFIG.retries):
            try:
                response = self.session.request(
                    "POST",
                    self.query_url,
                    json=payload,
                    timeout=SELENZYME_HTTP_CONFIG.timeout_seconds,
                )
                response.raise_for_status()
                envelope = response.json()
                if not isinstance(envelope, dict):
                    raise ValueError("Selenzyme response root is not an object")
                rows = [
                    _normalize_row(row, index)
                    for index, row in enumerate(
                        _decode_selenzyme_rows(envelope.get("data")),
                        start=1,
                    )
                ]
                raw_text = getattr(response, "text", "") or json.dumps(
                    envelope,
                    ensure_ascii=False,
                    sort_keys=True,
                )
                result = {
                    "cache_schema": "selenzyme_client.v3",
                    "query_id": query_id,
                    "query_type": query_type,
                    "query_database": database,
                    "query_value": normalized,
                    result_field: normalized,
                    "status": "ok" if rows else "no_hit",
                    "request": {"url": self.query_url, "payload": payload},
                    "app": str(envelope.get("app") or ""),
                    "version": str(envelope.get("version") or ""),
                    "rows": rows,
                    "response_sha256": hashlib.sha256(
                        raw_text.encode("utf-8")
                    ).hexdigest(),
                    "cache_hit": False,
                }
                _json_atomic(cache_path, result)
                return result
            except (requests.RequestException, ValueError, json.JSONDecodeError) as exc:
                last_error = exc
                if attempt + 1 < SELENZYME_HTTP_CONFIG.retries:
                    time.sleep(SELENZYME_HTTP_CONFIG.sleep_seconds * (2**attempt))
        raise SelenzymeSourceUnavailable(
            f"Selenzyme query failed for {normalized}: {last_error or 'unknown error'}"
        )

    def query_kegg_reaction(
        self,
        reaction_id: str,
        *,
        host_taxon_id: int,
        targets: int,
    ) -> dict[str, Any]:
        normalized = str(reaction_id or "").strip().upper()
        if not normalized:
            raise ValueError("Selenzyme query requires a KEGG reaction ID")
        return self._query_identifier(
            normalized,
            database="kegg",
            query_type="selenzyme_by_kegg_reaction",
            result_field="reaction_id",
            host_taxon_id=host_taxon_id,
            targets=targets,
        )

    def query_ec_number(
        self,
        ec_number: str,
        *,
        host_taxon_id: int,
        targets: int,
    ) -> dict[str, Any]:
        normalized = str(ec_number or "").strip()
        if not COMPLETE_EC_PATTERN.fullmatch(normalized):
            raise ValueError(
                "Selenzyme EC query requires a complete four-level EC number"
            )
        return self._query_identifier(
            normalized,
            database="ec",
            query_type="selenzyme_by_ec_number",
            result_field="ec_number",
            host_taxon_id=host_taxon_id,
            targets=targets,
        )


def retrieve_selenzyme_candidates(
    requirement: dict[str, Any],
    query_result: dict[str, Any],
    chassis_key: str,
    *,
    top_n: int,
    allow_transmembrane: bool,
    session: requests.Session,
    entry_cache: dict[str, dict[str, Any] | None] | None = None,
) -> tuple[list[ProteinCandidate], list[dict[str, Any]], list[str], dict[str, str]]:
    """Resolve ranked Selenzyme accessions to filtered UniProt candidates."""

    entry_cache = entry_cache if entry_cache is not None else {}
    candidates: list[ProteinCandidate] = []
    audit_rows: list[dict[str, Any]] = []
    query_ids: list[str] = []
    errors: dict[str, str] = {}
    seen_accessions: set[str] = set()
    query_type = str(query_result.get("query_type") or "")
    is_ec_query = query_type == "selenzyme_by_ec_number"
    queried_ec = str(query_result.get("ec_number") or "").strip()
    required_ecs = _ec_numbers(
        requirement.get("ec_numbers")
        or requirement.get("locked_ec_numbers")
        or requirement.get("locked_enzyme_ecs")
    )

    source_rows = [
        row
        for row in query_result.get("rows", [])
        if isinstance(row, dict)
    ]
    source_rows.sort(
        key=lambda row: (
            -float(row.get("reaction_similarity"))
            if isinstance(row.get("reaction_similarity"), (int, float))
            else float("inf"),
            int(row.get("rank") or 10**9),
        )
    )
    for source_row in source_rows:
        row = dict(source_row)
        accession = str(row.get("accession") or "").strip().upper()
        row["match_type"] = selenzyme_match_type(row)
        row["gate_status"] = "rejected"
        row["rejection_reasons"] = []
        if row["match_type"] == "invalid":
            row["rejection_reasons"] = ["missing_or_invalid_combined_reaction_similarity"]
            audit_rows.append(row)
            continue
        if is_ec_query and queried_ec not in _ec_numbers(row.get("ec_number")):
            row["rejection_reasons"] = [
                "selenzyme_ec_query_row_missing_requested_ec"
            ]
            audit_rows.append(row)
            continue
        if not accession:
            row["rejection_reasons"] = ["missing_uniprot_accession"]
            audit_rows.append(row)
            continue
        if accession in seen_accessions:
            row["rejection_reasons"] = ["duplicate_uniprot_accession"]
            audit_rows.append(row)
            continue
        seen_accessions.add(accession)
        uniprot_query_id = f"uniprot_accession_{accession}"
        query_ids.append(uniprot_query_id)
        if accession not in entry_cache:
            try:
                entries = search_uniprot_by_query(
                    f"accession:{accession}",
                    max_results=1,
                    session=session,
                )
                entry_cache[accession] = next(
                    (
                        entry
                        for entry in entries
                        if str(entry.get("primaryAccession") or "").strip().upper()
                        == accession
                    ),
                    None,
                )
            except Exception as exc:
                entry_cache[accession] = None
                errors[uniprot_query_id] = str(exc)
        entry = entry_cache.get(accession)
        if not entry:
            row["rejection_reasons"] = ["uniprot_entry_not_found"]
            audit_rows.append(row)
            continue

        passed, filter_reasons, filter_warnings = hard_filter_candidate_without_ec(
            entry,
            allow_transmembrane=allow_transmembrane,
        )
        if not passed:
            row["rejection_reasons"] = filter_reasons
            row["warnings"] = filter_warnings
            audit_rows.append(row)
            continue

        match_type = str(row["match_type"])
        retrieval_strategy = (
            "selenzyme_ec_risk"
            if is_ec_query
            else f"selenzyme_kegg_{match_type}"
        )
        candidate = candidate_from_reaction_entry(
            entry,
            chassis_key,
            retrieval_strategy=retrieval_strategy,
            retrieval_query_id=str(query_result.get("query_id") or ""),
            matched_rhea_ids=[],
            allow_transmembrane=allow_transmembrane,
            function_evidence_reason=(
                "function: SelenzymeRF EC association; locked substrate specificity unverified"
                if is_ec_query
                else "function: exact SelenzymeRF reaction match"
                if match_type == "exact"
                else "function: SelenzymeRF similar-reaction candidate"
            ),
        )
        if candidate is None:
            row["rejection_reasons"] = ["uniprot_sequence_safety_filter_failed"]
            audit_rows.append(row)
            continue
        candidate.retrieval_strategy = retrieval_strategy
        candidate.retrieval_query_id = str(query_result.get("query_id") or "")
        candidate_ecs = set(candidate.ec_numbers)
        reported_ecs = _ec_numbers(row.get("ec_number"))
        ec_relation = (
            classify_selenzyme_ec_relation(
                required_ecs=required_ecs,
                candidate_ecs=candidate_ecs,
                reported_ecs=reported_ecs,
            )
            if is_ec_query
            else ""
        )
        if is_ec_query:
            row["ec_relation"] = ec_relation
            row["candidate_ec_numbers"] = sorted(candidate_ecs)
            row["required_ec_numbers"] = sorted(required_ecs)
            row["reported_ec_numbers"] = sorted(reported_ecs)
        if is_ec_query and ec_relation == EC_RELATION_CONTRADICTED:
            row["rejection_reasons"] = [
                "selenzyme_candidate_ec_contradicts_requirement"
            ]
            audit_rows.append(row)
            continue
        candidate.reaction_confidence = (
            "selenzyme_ec_risk"
            if is_ec_query
            else "selenzyme_exact"
            if match_type == "exact"
            else "selenzyme_risk"
        )
        candidate.direction_support = "selenzyme_direction_neutral_no_known_conflict"
        candidate.selenzyme_rank = int(row.get("rank") or 0) or None
        candidate.selenzyme_score = row.get("score")
        candidate.selenzyme_sim_rf = row.get("sim_rf")
        candidate.selenzyme_sim_2018 = row.get("sim_2018")
        candidate.selenzyme_reaction_similarity = row.get("reaction_similarity")
        candidate.selenzyme_matched_reaction_id = str(
            row.get("matched_reaction_id") or ""
        )
        candidate.selenzyme_taxonomic_distance = row.get("taxonomic_distance")
        candidate.selenzyme_direction_used = str(row.get("direction_used") or "")
        candidate.selenzyme_direction_preferred = str(
            row.get("direction_preferred") or ""
        )
        if is_ec_query:
            candidate.selenzyme_risk_status = ec_relation
            relation_warning = {
                EC_RELATION_EXACT: (
                    "Current UniProt annotation contains the queried EC, but the locked substrate/product specificity remains unverified"
                ),
                EC_RELATION_SHARED_REACTION: (
                    "Current UniProt EC and the queried EC overlap only through the same Selenzyme reaction record; the locked substrate/product specificity remains unverified"
                ),
                EC_RELATION_UNANNOTATED: (
                    "SelenzymeRF EC association is not confirmed by the current UniProt EC annotation"
                ),
            }[ec_relation]
            current_activity = (
                f"Current UniProt catalytic activity: {candidate.catalytic_activities[0]}"
                if candidate.catalytic_activities
                else ""
            )
            candidate.reasons = list(dict.fromkeys(
                candidate.reasons
                + [
                    f"selenzyme_ec_relation:{ec_relation}",
                    "selenzyme_reported_ec_numbers:"
                    + ";".join(sorted(reported_ecs)),
                ]
            ))
            candidate.warnings = list(dict.fromkeys(
                candidate.warnings
                + [
                    "SelenzymeRF EC association does not establish the locked substrate/product reaction",
                    relation_warning,
                    *([current_activity] if current_activity else []),
                ]
            ))
        else:
            candidate.selenzyme_risk_status = (
                "" if match_type == "exact" else "combined_reaction_similarity_below_1"
            )
        if not is_ec_query and match_type != "exact":
            similarity = float(row["reaction_similarity"])
            candidate.warnings = list(dict.fromkeys(
                candidate.warnings
                + [
                    "SelenzymeRF risk fallback accepted: combined reaction similarity="
                    f"{similarity:.10g} (<1)"
                ]
            ))
        row["gate_status"] = "passed"
        row["rejection_reasons"] = []
        row["warnings"] = candidate.warnings
        audit_rows.append(row)
        candidates.append(candidate)
        if len(candidates) >= max(1, int(top_n)):
            break
    return candidates, audit_rows, query_ids, errors


__all__ = [
    "COMPLETE_EC_PATTERN",
    "EC_RELATION_CONTRADICTED",
    "EC_RELATION_EXACT",
    "EC_RELATION_SHARED_REACTION",
    "EC_RELATION_UNANNOTATED",
    "EXACT_SIMILARITY_TOLERANCE",
    "SelenzymeClient",
    "SelenzymeSourceUnavailable",
    "_decode_selenzyme_rows",
    "chassis_host_taxon_id",
    "classify_selenzyme_ec_relation",
    "retrieve_selenzyme_candidates",
    "selenzyme_match_type",
    "selenzyme_target_count",
]
