"""Retrieve and validate expression-part facts from the remote Milvus store."""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dotenv import load_dotenv

from src.config.config import ROOT
from src.expression_box.config import (
    EXPRESSION_PARTS_COLLECTION,
    MILVUS_QUERY_LIMIT,
    MILVUS_TIMEOUT_SECONDS,
)
from src.expression_box.expression_host_context import candidate_host_match_kind
from src.expression_box.parts_models import (
    ExpressionPartsContext,
    PartCandidate,
    PartsSnapshot,
)


_DNA_RE = re.compile(r"[ACGT]+")
_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_ROLES = ("promoter", "rbs", "terminator")
_STRENGTHS = {"low", "medium", "high"}
_FORWARD_DIRECTIONS = {"forward", "both", "bidirectional"}
_NO_REGULATOR = {"", "none", "not_required", "not_applicable"}
_REQUIRED_FIELDS = {
    "part_id",
    "role",
    "role_confidence",
    "sequence_type",
    "sequence_type_confidence",
    "sequence",
    "sequence_sha256",
    "sequence_available",
    "length_bp",
    "host",
    "host_confidence",
    "strength",
    "strength_confidence",
    "regulation",
    "regulation_confidence",
    "warnings",
    "registry_metadata",
    "activity_value",
    "activity_percentile",
    "activity_dataset",
    "activity_unit",
    "activity_context",
    "activity_source",
    "evidence_grade",
    "direction",
    "regulator_required",
    "source",
    "evidence",
    "updated_at",
}
_OUTPUT_FIELDS = sorted(_REQUIRED_FIELDS | {"doc_id"})


class MilvusPartsError(RuntimeError):
    """Raised when the authoritative remote parts store cannot be used."""


def _stable_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _clean(value: Any) -> str:
    return str(value or "").strip()


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _float_or_none(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 0 else None


def _list_of_text(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [_clean(item) for item in value if _clean(item)]


def _client_from_environment():
    load_dotenv(Path(ROOT) / ".env")
    host = _clean(os.getenv("MILVUS_HOST"))
    port = _clean(os.getenv("MILVUS_PORT")) or "19530"
    if not host:
        raise MilvusPartsError("MILVUS_HOST is not configured in .env")
    uri = host.rstrip("/") if host.startswith(("http://", "https://")) else f"http://{host}:{port}"
    kwargs: dict[str, Any] = {"uri": uri}
    token = _clean(os.getenv("MILVUS_TOKEN"))
    if token:
        kwargs["token"] = token
    db_name = _clean(os.getenv("MILVUS_DB_NAME"))
    if db_name:
        kwargs["db_name"] = db_name
    try:
        from pymilvus import MilvusClient
    except ImportError as exc:
        raise MilvusPartsError(
            "pymilvus is unavailable; install pymilvus==2.6.17"
        ) from exc
    try:
        return MilvusClient(**kwargs)
    except Exception as exc:
        raise MilvusPartsError(
            f"could not connect to the configured Milvus service: {type(exc).__name__}"
        ) from exc


def _schema_description(client: Any, collection_name: str) -> tuple[dict[str, Any], str]:
    try:
        if not client.has_collection(collection_name):
            raise MilvusPartsError(
                f"Milvus collection not found: {collection_name}"
            )
        client.load_collection(collection_name)
        description = client.describe_collection(collection_name)
    except MilvusPartsError:
        raise
    except Exception as exc:
        raise MilvusPartsError(
            f"could not inspect Milvus collection {collection_name}: {type(exc).__name__}"
        ) from exc
    fields = {
        str(field.get("name") or ""): {
            "type": str(field.get("type") or field.get("data_type") or ""),
            "params": field.get("params") or {},
        }
        for field in description.get("fields", [])
        if field.get("name")
    }
    missing = sorted(_REQUIRED_FIELDS - set(fields))
    if missing:
        raise MilvusPartsError(
            f"Milvus collection {collection_name} is missing fields: {missing}"
        )
    return description, _stable_hash(fields)


def _query_role(
    client: Any,
    *,
    collection_name: str,
    role: str,
    host_labels: tuple[str, ...],
) -> list[dict[str, Any]]:
    encoded_labels = json.dumps(list(host_labels), ensure_ascii=False)
    expression = (
        f'role == "{role}" and sequence_available == true '
        f"and ARRAY_CONTAINS_ANY(host, {encoded_labels})"
    )
    try:
        rows = client.query(
            collection_name=collection_name,
            filter=expression,
            output_fields=_OUTPUT_FIELDS,
            limit=MILVUS_QUERY_LIMIT,
            timeout=MILVUS_TIMEOUT_SECONDS,
        )
    except Exception as exc:
        raise MilvusPartsError(
            f"Milvus query failed for role {role}: {type(exc).__name__}"
        ) from exc
    return [dict(row) for row in rows if isinstance(row, Mapping)]


def _rejection_reason(
    row: Mapping[str, Any],
    *,
    expected_role: str,
    context: ExpressionPartsContext,
) -> str | None:
    part_id = _clean(row.get("part_id"))
    if not part_id:
        return "missing_part_id"
    role = _lower(row.get("role"))
    if role != expected_role:
        return "role_mismatch"
    if _lower(row.get("sequence_type")) != role:
        return "sequence_type_mismatch"
    sequence = re.sub(r"\s+", "", _clean(row.get("sequence"))).upper()
    if _DNA_RE.fullmatch(sequence) is None:
        return "invalid_sequence"
    try:
        length_bp = int(row.get("length_bp") or 0)
    except (TypeError, ValueError):
        return "invalid_length"
    if length_bp != len(sequence):
        return "length_mismatch"
    recorded_sha = _lower(row.get("sequence_sha256"))
    computed_sha = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    if _SHA256_RE.fullmatch(recorded_sha) is None or recorded_sha != computed_sha:
        return "sequence_hash_mismatch"
    if candidate_host_match_kind(_list_of_text(row.get("host")), _host_resolution(context)) == "none":
        return "host_mismatch"
    if _lower(row.get("strength")) not in _STRENGTHS:
        return "unsupported_strength"
    if _lower(row.get("direction")) not in _FORWARD_DIRECTIONS:
        return "unsupported_direction"
    registry = row.get("registry_metadata")
    registry = registry if isinstance(registry, Mapping) else {}
    failure_text = " ".join(
        [
            _lower(registry.get("part_results")),
            *(_lower(value) for value in _list_of_text(row.get("warnings"))),
        ]
    )
    if "fails" in failure_text or "should not work" in failure_text:
        return "failed_evidence"
    if role == "promoter":
        if _lower(row.get("regulation")) != "constitutive":
            return "non_constitutive_promoter"
        if _lower(row.get("regulator_required")) not in _NO_REGULATOR:
            return "regulator_required"
    return None


def _host_resolution(context: ExpressionPartsContext):
    """Return the small protocol object expected by candidate_host_match_kind."""

    from src.expression_box.expression_host_context import ExpressionHostResolution

    return ExpressionHostResolution(
        codon_transformer_name=context.host_name,
        codon_transformer_organism_id=52,
        requested_name="",
        expression_host_key=context.host_key,
        milvus_host_labels=context.host_labels,
        query_name=context.host_name,
    )


def _candidate(row: Mapping[str, Any], context: ExpressionPartsContext) -> PartCandidate:
    sequence = re.sub(r"\s+", "", _clean(row.get("sequence"))).upper()
    host_match = candidate_host_match_kind(
        _list_of_text(row.get("host")),
        _host_resolution(context),
    )
    return PartCandidate(
        part_id=_clean(row.get("part_id")),
        role=_lower(row.get("role")),
        sequence=sequence,
        sequence_sha256=hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
        sequence_type=_lower(row.get("sequence_type")),
        host_match_kind=host_match,
        strength=_lower(row.get("strength")),
        regulation=_lower(row.get("regulation")),
        direction=_lower(row.get("direction")),
        regulator_required=_lower(row.get("regulator_required")),
        evidence_grade=_clean(row.get("evidence_grade")).upper() or "C",
        role_confidence=_lower(row.get("role_confidence")),
        sequence_type_confidence=_lower(row.get("sequence_type_confidence")),
        host_confidence=_lower(row.get("host_confidence")),
        strength_confidence=_lower(row.get("strength_confidence")),
        regulation_confidence=_lower(row.get("regulation_confidence")),
        activity_value=_float_or_none(row.get("activity_value")),
        activity_percentile=_float_or_none(row.get("activity_percentile")),
        activity_dataset=_clean(row.get("activity_dataset")),
        activity_unit=_clean(row.get("activity_unit")),
        activity_context=_clean(row.get("activity_context")),
        activity_source=_clean(row.get("activity_source")),
        source=_clean(row.get("source")),
        evidence=_clean(row.get("evidence")),
        warnings=tuple(_list_of_text(row.get("warnings"))),
        updated_at=_clean(row.get("updated_at")),
    )


def fetch_expression_part_candidates(
    context: ExpressionPartsContext,
    *,
    client: Any | None = None,
    collection_name: str = EXPRESSION_PARTS_COLLECTION,
) -> PartsSnapshot:
    """Return a normalized snapshot of authoritative remote candidates."""

    owned_client = client is None
    milvus = client if client is not None else _client_from_environment()
    try:
        _, schema_fingerprint = _schema_description(milvus, collection_name)
        raw_by_role = {
            role: _query_role(
                milvus,
                collection_name=collection_name,
                role=role,
                host_labels=context.host_labels,
            )
            for role in _ROLES
        }
        try:
            stats = milvus.get_collection_stats(collection_name)
            row_count = int(stats.get("row_count") or 0)
        except Exception:
            row_count = sum(len(rows) for rows in raw_by_role.values())
    finally:
        if owned_client:
            close = getattr(milvus, "close", None)
            if callable(close):
                close()

    accepted: list[PartCandidate] = []
    rejected: Counter[str] = Counter()
    raw_counts = {role: len(rows) for role, rows in raw_by_role.items()}
    for role, rows in raw_by_role.items():
        for row in rows:
            reason = _rejection_reason(row, expected_role=role, context=context)
            if reason is not None:
                rejected[reason] += 1
                continue
            accepted.append(_candidate(row, context))
    accepted.sort(key=lambda item: (item.role, item.part_id, item.sequence_sha256))
    identities = [(item.role, item.part_id) for item in accepted]
    if len(identities) != len(set(identities)):
        raise MilvusPartsError("Milvus returned duplicate role/part_id candidates")
    accepted_counts = Counter(item.role for item in accepted)
    missing_roles = [role for role in _ROLES if not accepted_counts[role]]
    if missing_roles:
        raise MilvusPartsError(
            "Milvus has no usable MG1655 candidates for roles: "
            + ", ".join(missing_roles)
        )
    fingerprint = _stable_hash(
        [item.fingerprint_payload() for item in accepted]
    )
    return PartsSnapshot(
        collection_name=collection_name,
        collection_row_count=row_count,
        schema_fingerprint=schema_fingerprint,
        candidate_fingerprint=fingerprint,
        candidates=tuple(accepted),
        raw_counts=raw_counts,
        accepted_counts={role: accepted_counts[role] for role in _ROLES},
        rejected_counts=dict(sorted(rejected.items())),
    )


__all__ = [
    "MilvusPartsError",
    "fetch_expression_part_candidates",
]
