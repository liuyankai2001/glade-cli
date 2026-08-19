"""Atomic artifact and TTL-cache storage for literature activity research."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from dataclasses import asdict
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

from src.main_protein_selection.literature_activity.models import (
    LiteratureActivityArtifact,
    LiteratureActivitySearchResult,
)
from src.main_protein_selection.uniprot_protein_candidates import ProteinCandidate


JSON_FILENAME = "literature_activity_evidence.json"
CSV_FILENAME = "literature_activity_evidence.csv"
_CACHE_SCHEMA_VERSION = "literature_activity_cache.v1"
_CSV_COLUMNS = (
    "evidence_id",
    "step_index",
    "reaction_id",
    "resolved_accession",
    "gene_name",
    "protein_name",
    "organism_name",
    "evidence_level",
    "assay_type",
    "direct_activity_measured",
    "direction",
    "identity_status",
    "fit_status",
    "review_status",
    "doi",
    "pmid",
    "pmcid",
    "title",
    "source_locator",
    "evidence_summary",
    "limitations",
    "rejection_reasons",
)


def utc_now_text() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def artifact_fingerprint(payload: dict[str, Any] | LiteratureActivityArtifact) -> str:
    """Hash stable scientific content, excluding runtime/cache bookkeeping."""

    data = (
        payload.model_dump(mode="json")
        if isinstance(payload, LiteratureActivityArtifact)
        else dict(payload)
    )
    data.pop("artifact_fingerprint", None)
    data.pop("generated_at", None)
    data.pop("cache_hit", None)
    canonical = json.dumps(
        data,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=lambda value: value.model_dump(mode="json"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def finalize_artifact(payload: dict[str, Any]) -> LiteratureActivityArtifact:
    """Validate an artifact after deriving its stable content fingerprint."""

    data = dict(payload)
    # Validate once with a placeholder so schema defaults participate in the
    # canonical hash exactly as they do when an artifact is read from disk.
    data["artifact_fingerprint"] = "0" * 64
    provisional = LiteratureActivityArtifact.model_validate(data)
    normalized = provisional.model_dump(mode="json")
    normalized["artifact_fingerprint"] = artifact_fingerprint(provisional)
    return LiteratureActivityArtifact.model_validate(normalized)


def _write_text_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=str(path.parent),
        text=True,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def write_artifact(
    artifact: LiteratureActivityArtifact,
    output_dir: str | Path,
) -> tuple[Path, Path]:
    """Atomically write canonical JSON and a compact human-auditable CSV."""

    directory = Path(output_dir).expanduser().resolve(strict=False)
    json_path = directory / JSON_FILENAME
    csv_path = directory / CSV_FILENAME
    _write_text_atomic(
        json_path,
        json.dumps(
            artifact.model_dump(mode="json"),
            ensure_ascii=False,
            indent=2,
        ) + "\n",
    )

    with tempfile.TemporaryFile(mode="w+", encoding="utf-8", newline="") as buffer:
        writer = csv.DictWriter(buffer, fieldnames=_CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for evidence in artifact.evidence:
            row = evidence.model_dump(mode="json")
            for field_name in ("limitations", "rejection_reasons"):
                row[field_name] = " | ".join(row.get(field_name, []))
            writer.writerow(row)
        buffer.seek(0)
        _write_text_atomic(csv_path, buffer.read())
    return json_path, csv_path


def _cache_path(cache_dir: str | Path, request_fingerprint: str) -> Path:
    return (
        Path(cache_dir).expanduser().resolve(strict=False)
        / "literature_activity"
        / f"{request_fingerprint}.json"
    )


def load_cached_result(
    *,
    cache_dir: str | Path,
    output_dir: str | Path,
    request_fingerprint: str,
) -> LiteratureActivitySearchResult | None:
    """Load one unexpired result without invoking a model or external source."""

    path = _cache_path(cache_dir, request_fingerprint)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != _CACHE_SCHEMA_VERSION
        or payload.get("request_fingerprint") != request_fingerprint
    ):
        return None
    try:
        expires_at = datetime.fromisoformat(
            str(payload.get("expires_at") or "").replace("Z", "+00:00")
        )
    except ValueError:
        return None
    if expires_at <= datetime.now(UTC):
        return None

    try:
        artifact_data = dict(payload["artifact"])
        artifact_data["generated_at"] = utc_now_text()
        artifact_data["cache_hit"] = True
        artifact_data["artifact_fingerprint"] = artifact_fingerprint(artifact_data)
        artifact = LiteratureActivityArtifact.model_validate(artifact_data)
        candidates_by_step = {
            int(step): [ProteinCandidate(**item) for item in candidates]
            for step, candidates in dict(payload.get("candidates_by_step") or {}).items()
        }
    except (KeyError, TypeError, ValueError):
        return None

    json_path, csv_path = write_artifact(artifact, output_dir)
    query_errors = {
        failure.query_id or failure.failure_id: failure.message
        for failure in artifact.failures
    }
    return LiteratureActivitySearchResult(
        status=artifact.status,
        candidates_by_step=candidates_by_step,
        artifact=artifact,
        json_path=json_path,
        csv_path=csv_path,
        query_errors=query_errors,
    )


def store_cached_result(
    result: LiteratureActivitySearchResult,
    *,
    cache_dir: str | Path,
) -> None:
    """Cache successful science longer than transient source failures."""

    if result.status == "disabled":
        return
    if result.status == "source_unavailable":
        ttl = timedelta(hours=1)
    elif result.status == "not_found":
        ttl = timedelta(days=7)
    else:
        ttl = timedelta(days=30)
    now = datetime.now(UTC)
    payload = {
        "schema_version": _CACHE_SCHEMA_VERSION,
        "request_fingerprint": result.artifact.request_fingerprint,
        "created_at": now.isoformat().replace("+00:00", "Z"),
        "expires_at": (now + ttl).isoformat().replace("+00:00", "Z"),
        "artifact": result.artifact.model_dump(mode="json"),
        "candidates_by_step": {
            str(step): [asdict(candidate) for candidate in candidates]
            for step, candidates in result.candidates_by_step.items()
        },
    }
    _write_text_atomic(
        _cache_path(cache_dir, result.artifact.request_fingerprint),
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
    )


__all__ = [
    "CSV_FILENAME",
    "JSON_FILENAME",
    "artifact_fingerprint",
    "finalize_artifact",
    "load_cached_result",
    "store_cached_result",
    "utc_now_text",
    "write_artifact",
]
