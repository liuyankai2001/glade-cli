from __future__ import annotations

import hashlib
import re
from pathlib import Path


BIOBRICK_PART_ID_RE = re.compile(r"^BBa_[A-Za-z0-9_]+$", re.IGNORECASE)
SOURCE_SCOPED_PART_ID_RE = re.compile(
    r"^[A-Za-z][A-Za-z0-9_.-]{1,31}:[A-Za-z0-9][A-Za-z0-9_.-]{0,95}$"
)


def normalize_part_id(part_id: str) -> str:
    """Validate a BioBrick or source-scoped identifier without losing provenance."""

    normalized = str(part_id or "").strip()
    if normalized.lower().startswith("igem:"):
        normalized = normalized.split(":", 1)[1]
    if BIOBRICK_PART_ID_RE.fullmatch(normalized):
        return normalized
    if SOURCE_SCOPED_PART_ID_RE.fullmatch(normalized):
        namespace, local_id = normalized.split(":", 1)
        return f"{namespace.upper()}:{local_id}"
    raise ValueError(
        "part_id must be a BioBrick BBa_* id or a source-scoped id such as "
        f"KOSURI2013:P001: {part_id}"
    )


def part_storage_key(part_id: str) -> str:
    """Return a deterministic Windows-safe key while preserving legacy BBa paths."""

    normalized = normalize_part_id(part_id)
    if BIOBRICK_PART_ID_RE.fullmatch(normalized):
        return normalized
    readable = re.sub(r"[^A-Za-z0-9_.-]+", "__", normalized)
    digest = hashlib.sha256(normalized.encode("utf-8")).hexdigest()[:10]
    return f"{readable}__{digest}"


def part_sequence_path(parts_dir: Path, part_id: str) -> Path:
    return parts_dir / f"{part_storage_key(part_id)}.txt"
