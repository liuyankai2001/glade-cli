"""Minimal KEGG reaction client used during input validation."""

from collections.abc import Callable
import json
import re
import time

import httpx

from src.services.errors import DatabaseLookupError, LookupErrorKind
from src.services.cache import (
    CURATED_RECORD_TTL_SECONDS,
    PersistentTTLCache,
    RetrievalStats,
)
from src.services.http import get_with_retry
from src.state import ReactionRecord


class KeggClient:
    """Validate reaction identifiers against the official KEGG REST API."""

    SERVICE_NAME = "KEGG"
    BASE_URL = "https://rest.kegg.jp"
    USER_AGENT = "protein-supply/0.1.0"
    ENTRY_PATTERN = re.compile(r"^ENTRY\s+(R\d{5})\s+Reaction\s*$", re.MULTILINE)
    FIELD_PATTERN = re.compile(r"^([A-Z][A-Z0-9_]*)\s+(.+)$")
    ORTHOLOGY_PATTERN = re.compile(r"^(K\d{5})\s+(.+)$")

    def __init__(
        self,
        client: httpx.Client,
        *,
        sleep: Callable[[float], None] = time.sleep,
        cache: PersistentTTLCache | None = None,
        stats: RetrievalStats | None = None,
        attempts: int = 2,
    ) -> None:
        self._client = client
        self._sleep = sleep
        self._cache = cache
        self._stats = stats
        self._attempts = attempts

    def get_reaction(self, reaction_id: str) -> ReactionRecord:
        """Return a compact record for one KEGG reaction."""

        url = f"{self.BASE_URL}/get/rn:{reaction_id}"
        cache_key = json.dumps(
            {"reaction_id": reaction_id},
            sort_keys=True,
            separators=(",", ":"),
        )
        if self._cache is not None:
            cached = self._cache.get("input.kegg", cache_key)
            if isinstance(cached, dict):
                if self._stats is not None:
                    self._stats.cache_hits += 1
                return cached  # type: ignore[return-value]
        if self._stats is not None:
            self._stats.live_calls += 1
        try:
            response = get_with_retry(
                self._client,
                url,
                headers={"Accept": "text/plain", "User-Agent": self.USER_AGENT},
                sleep=self._sleep,
                on_retry=self._record_retry,
                attempts=self._attempts,
            )
        except httpx.RequestError as exc:
            if self._stats is not None:
                self._stats.failed_calls += 1
            raise DatabaseLookupError(
                kind=LookupErrorKind.SERVICE_ERROR,
                service=self.SERVICE_NAME,
                identifier=reaction_id,
                message=f"KEGG request failed: {exc}",
            ) from exc

        if response.status_code == 404:
            self._record_failure()
            raise DatabaseLookupError(
                kind=LookupErrorKind.NOT_FOUND,
                service=self.SERVICE_NAME,
                identifier=reaction_id,
                message=f"KEGG reaction {reaction_id} was not found",
            )
        if response.status_code >= 400:
            self._record_failure()
            raise DatabaseLookupError(
                kind=LookupErrorKind.SERVICE_ERROR,
                service=self.SERVICE_NAME,
                identifier=reaction_id,
                message=f"KEGG returned HTTP {response.status_code}",
            )

        entry_match = self.ENTRY_PATTERN.search(response.text)
        if entry_match is None or entry_match.group(1) != reaction_id:
            self._record_failure()
            raise DatabaseLookupError(
                kind=LookupErrorKind.INVALID_RESPONSE,
                service=self.SERVICE_NAME,
                identifier=reaction_id,
                message="KEGG response did not contain the requested reaction entry",
            )

        fields = self._parse_flat_file(response.text)
        names = self._split_semicolon_values(fields.get("NAME", []))
        enzyme_ids = [
            value
            for line in fields.get("ENZYME", [])
            for value in line.split()
            if value
        ]
        orthology = []
        for line in fields.get("ORTHOLOGY", []):
            match = self.ORTHOLOGY_PATTERN.fullmatch(line)
            if match is not None:
                orthology.append(
                    {
                        "orthology_id": match.group(1),
                        "description": match.group(2).strip(),
                    }
                )
        rhea_ids: list[str] = []
        for line in fields.get("DBLINKS", []):
            if not line.upper().startswith("RHEA:"):
                continue
            rhea_ids.extend(
                f"RHEA:{token}"
                for token in line.split(":", maxsplit=1)[1].split()
                if token.isdigit()
            )
        result: ReactionRecord = {
            "reaction_id": entry_match.group(1),
            "name": names[0] if names else None,
            "names": names,
            "definition": self._join_field(fields.get("DEFINITION", [])),
            "equation": self._join_field(fields.get("EQUATION", [])),
            "enzyme_ids": list(dict.fromkeys(enzyme_ids)),
            "orthology": orthology,
            "rhea_ids": list(dict.fromkeys(rhea_ids)),
            "source_url": str(response.url),
        }
        if self._cache is not None:
            self._cache.set(
                "input.kegg",
                cache_key,
                result,
                ttl_seconds=CURATED_RECORD_TTL_SECONDS,
            )
        return result

    def _record_retry(self) -> None:
        if self._stats is not None:
            self._stats.retries += 1

    def _record_failure(self) -> None:
        if self._stats is not None:
            self._stats.failed_calls += 1

    @classmethod
    def _parse_flat_file(cls, text: str) -> dict[str, list[str]]:
        """Parse the continuation-line format used by KEGG flat files."""

        fields: dict[str, list[str]] = {}
        current_field: str | None = None
        for raw_line in text.splitlines():
            if raw_line == "///":
                break
            match = cls.FIELD_PATTERN.match(raw_line)
            if match is not None:
                current_field = match.group(1)
                fields.setdefault(current_field, []).append(
                    match.group(2).strip()
                )
                continue
            if current_field is None or not raw_line[:1].isspace():
                continue
            continuation = raw_line.strip()
            if continuation:
                fields[current_field].append(continuation)
        return fields

    @staticmethod
    def _split_semicolon_values(lines: list[str]) -> list[str]:
        combined = " ".join(lines)
        return [
            value.strip()
            for value in combined.split(";")
            if value.strip()
        ]

    @staticmethod
    def _join_field(lines: list[str]) -> str | None:
        combined = " ".join(lines).strip()
        return combined or None
