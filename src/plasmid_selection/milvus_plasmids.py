"""Scalar retrieval of audited plasmid-template cards from remote Milvus."""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dotenv import load_dotenv

from src.config.config import ROOT
from src.plasmid_selection.config import (
    MILVUS_QUERY_LIMIT,
    MILVUS_TIMEOUT_SECONDS,
    PLASMID_COLLECTION,
)
from src.plasmid_selection.models import PlasmidSnapshot


REQUIRED_FIELDS = {
    "plasmid_id",
    "name",
    "schema_version",
    "source",
    "source_url",
    "sequence_url",
    "source_record_id",
    "source_record_version",
    "sequence_file",
    "source_file_sha256",
    "sequence_content_sha256",
    "canonical_sequence_sha256",
    "length_bp",
    "topology",
    "replicon_family",
    "copy_number_class",
    "cargo_type",
    "assembly_policy",
    "audit_status",
    "audit_passed",
    "mg1655_compatible",
    "has_expression_cargo",
    "origins",
    "resistance_markers",
    "selection_markers",
    "host_compatibility",
    "replication_dependencies",
    "insertion_regions",
    "protected_features",
    "source_provenance",
    "evidence_refs",
}

OPTIONAL_FIELDS = {
    "addgene_id",
    "description",
    "depositor",
    "article",
    "pubmed_id",
    "vector_type",
    "bacterial_resistance",
    "growth_strain",
    "copy_number",
    "cloning_method",
    "backbone",
    "requires_cargo_replacement",
    "sequence_sha256",
    "module_boundaries",
    "regulatory_contexts",
    "mcs_features",
    "features",
    "audit_version",
    "audit_checks",
    "normalization_applied",
    "normalization_events",
}


class MilvusPlasmidError(RuntimeError):
    """The authoritative remote plasmid collection could not be used."""


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


def _client_from_environment():
    load_dotenv(Path(ROOT) / ".env")
    host = _clean(os.getenv("MILVUS_HOST"))
    port = _clean(os.getenv("MILVUS_PORT")) or "19530"
    if not host:
        raise MilvusPlasmidError("MILVUS_HOST is not configured in .env")
    uri = (
        host.rstrip("/")
        if host.startswith(("http://", "https://"))
        else f"http://{host}:{port}"
    )
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
        raise MilvusPlasmidError(
            "pymilvus is unavailable; install pymilvus==2.6.17"
        ) from exc
    try:
        return MilvusClient(**kwargs)
    except Exception as exc:
        raise MilvusPlasmidError(
            "could not connect to the configured Milvus service: "
            f"{type(exc).__name__}"
        ) from exc


def _collection_schema(
    client: Any,
    collection_name: str,
) -> tuple[set[str], str]:
    try:
        if not client.has_collection(collection_name):
            raise MilvusPlasmidError(
                f"Milvus collection not found: {collection_name}"
            )
        client.load_collection(collection_name)
        description = client.describe_collection(collection_name)
    except MilvusPlasmidError:
        raise
    except Exception as exc:
        raise MilvusPlasmidError(
            f"could not inspect Milvus collection {collection_name}: "
            f"{type(exc).__name__}"
        ) from exc
    field_descriptions = {
        str(field.get("name") or ""): {
            "type": str(field.get("type") or field.get("data_type") or ""),
            "params": field.get("params") or {},
        }
        for field in description.get("fields", [])
        if isinstance(field, Mapping) and field.get("name")
    }
    missing = sorted(REQUIRED_FIELDS - set(field_descriptions))
    if missing:
        raise MilvusPlasmidError(
            f"Milvus collection {collection_name} is missing fields: {missing}"
        )
    return set(field_descriptions), _stable_hash(field_descriptions)


def _json_like(value: Any, expected: type) -> Any:
    if isinstance(value, expected):
        return value
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError:
            return expected()
        return parsed if isinstance(parsed, expected) else expected()
    return expected()


def _normalize_row(row: Mapping[str, Any], fields: set[str]) -> dict[str, Any]:
    normalized = {
        key: row.get(key)
        for key in sorted((REQUIRED_FIELDS | OPTIONAL_FIELDS) & fields)
    }
    for key in (
        "origins",
        "resistance_markers",
        "selection_markers",
        "host_compatibility",
        "replication_dependencies",
        "insertion_regions",
        "protected_features",
        "evidence_refs",
        "regulatory_contexts",
        "mcs_features",
        "normalization_events",
    ):
        if key in normalized:
            normalized[key] = _json_like(normalized.get(key), list)
    for key in (
        "source_provenance",
        "module_boundaries",
        "features",
        "audit_checks",
    ):
        if key in normalized:
            normalized[key] = _json_like(normalized.get(key), dict)
    for key in (
        "plasmid_id",
        "name",
        "schema_version",
        "source",
        "source_url",
        "sequence_url",
        "source_record_id",
        "source_record_version",
        "sequence_file",
        "source_file_sha256",
        "sequence_sha256",
        "sequence_content_sha256",
        "canonical_sequence_sha256",
        "topology",
        "replicon_family",
        "copy_number_class",
        "cargo_type",
        "assembly_policy",
        "audit_status",
        "bacterial_resistance",
    ):
        if key in normalized:
            normalized[key] = _clean(normalized.get(key))
    try:
        normalized["length_bp"] = int(normalized.get("length_bp") or 0)
    except (TypeError, ValueError):
        normalized["length_bp"] = 0
    return normalized


def fetch_plasmid_snapshot(
    *,
    client: Any | None = None,
    collection_name: str = PLASMID_COLLECTION,
) -> PlasmidSnapshot:
    """Query all v2 cards using scalar fields; no embedding model is used."""

    owned_client = client is None
    milvus = client if client is not None else _client_from_environment()
    try:
        fields, schema_fingerprint = _collection_schema(milvus, collection_name)
        output_fields = sorted((REQUIRED_FIELDS | OPTIONAL_FIELDS) & fields)
        try:
            raw_rows = milvus.query(
                collection_name=collection_name,
                filter='schema_version == "plasmid_template.v2"',
                output_fields=output_fields,
                limit=MILVUS_QUERY_LIMIT,
                timeout=MILVUS_TIMEOUT_SECONDS,
            )
        except Exception as exc:
            raise MilvusPlasmidError(
                f"Milvus query failed: {type(exc).__name__}"
            ) from exc
        try:
            stats = milvus.get_collection_stats(collection_name)
            row_count = int(stats.get("row_count") or 0)
        except Exception:
            row_count = len(raw_rows)
    finally:
        if owned_client:
            close = getattr(milvus, "close", None)
            if callable(close):
                close()

    candidates = [
        _normalize_row(row, fields)
        for row in raw_rows
        if isinstance(row, Mapping)
    ]
    candidates.sort(key=lambda item: (item.get("plasmid_id", ""), item.get("name", "")))
    identities = [item.get("plasmid_id") for item in candidates]
    if any(not item for item in identities):
        raise MilvusPlasmidError("Milvus returned a plasmid card without plasmid_id")
    if len(identities) != len(set(identities)):
        raise MilvusPlasmidError("Milvus returned duplicate plasmid_id values")
    if not candidates:
        raise MilvusPlasmidError(
            f"Milvus collection {collection_name} has no plasmid_template.v2 cards"
        )
    return PlasmidSnapshot(
        collection_name=collection_name,
        collection_row_count=row_count,
        schema_fingerprint=schema_fingerprint,
        candidate_fingerprint=_stable_hash(candidates),
        candidates=tuple(candidates),
    )


__all__ = [
    "MilvusPlasmidError",
    "REQUIRED_FIELDS",
    "fetch_plasmid_snapshot",
]
