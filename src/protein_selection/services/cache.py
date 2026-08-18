"""Small versioned TTL cache for successful public-database responses."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import time
from typing import Any, Callable


CACHE_SCHEMA_VERSION = 1
DAY_SECONDS = 24 * 60 * 60
SEARCH_TTL_SECONDS = DAY_SECONDS
CURATED_RECORD_TTL_SECONDS = 7 * DAY_SECONDS
ARTICLE_TTL_SECONDS = 30 * DAY_SECONDS


@dataclass(slots=True)
class RetrievalStats:
    """Mutable query-scoped counters exposed only through audit output."""

    cache_hits: int = 0
    live_calls: int = 0
    retries: int = 0
    failed_calls: int = 0

    def as_dict(self) -> dict[str, int]:
        return {
            "cache_hits": self.cache_hits,
            "live_calls": self.live_calls,
            "retries": self.retries,
            "failed_calls": self.failed_calls,
        }


class PersistentTTLCache:
    """Store JSON-serializable values in atomic, content-keyed files."""

    def __init__(
        self,
        root: str | Path,
        *,
        enabled: bool = True,
        refresh: bool = False,
        clock: Callable[[], float] = time.time,
    ) -> None:
        self.root = Path(root)
        self.enabled = enabled
        self.refresh = refresh
        self._clock = clock

    def get(self, namespace: str, key: str) -> Any | None:
        if not self.enabled or self.refresh:
            return None
        path = self._path(namespace, key)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, TypeError, ValueError):
            return None
        if not isinstance(payload, dict):
            return None
        if payload.get("schema_version") != CACHE_SCHEMA_VERSION:
            return None
        expires_at = payload.get("expires_at")
        if not isinstance(expires_at, (int, float)):
            return None
        if expires_at <= self._clock():
            return None
        return payload.get("value")

    def set(
        self,
        namespace: str,
        key: str,
        value: Any,
        *,
        ttl_seconds: float,
    ) -> None:
        if not self.enabled or ttl_seconds <= 0:
            return
        path = self._path(namespace, key)
        now = self._clock()
        payload = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "created_at": now,
            "expires_at": now + ttl_seconds,
            "value": value,
        }
        temporary: Path | None = None
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(
                f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
            )
            temporary.write_text(
                json.dumps(
                    payload,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                    default=str,
                ),
                encoding="utf-8",
            )
            os.replace(temporary, path)
        except OSError:
            # A cache is an optimization; filesystem failures must not make a
            # biological query fail.
            if temporary is not None:
                try:
                    temporary.unlink(missing_ok=True)
                except OSError:
                    pass
            return

    def _path(self, namespace: str, key: str) -> Path:
        safe_namespace = re.sub(r"[^a-zA-Z0-9_.-]+", "_", namespace)
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.root / safe_namespace / digest[:2] / f"{digest}.json"


def ttl_for_tool(tool_name: str) -> float:
    """Return a conservative freshness window for one public source call."""

    normalized = tool_name.casefold()
    if "search" in normalized or tool_name in {
        "STRING_map_identifiers",
        "STRING_get_protein_interactions",
        "intact_search_interactions",
        "intact_get_interactions",
        "fetchWebContent",
    }:
        return SEARCH_TTL_SECONDS
    if tool_name in {
        "PubMed_get_article",
        "SemanticScholar_get_paper",
        "SemanticScholar_get_pdf_snippets",
    }:
        return ARTICLE_TTL_SECONDS
    return CURATED_RECORD_TTL_SECONDS
