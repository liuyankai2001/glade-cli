"""Minimal UniProtKB client used during input validation."""

from collections.abc import Callable
import json
import time
from typing import Any

import httpx

from src.services.errors import DatabaseLookupError, LookupErrorKind
from src.services.cache import (
    CURATED_RECORD_TTL_SECONDS,
    PersistentTTLCache,
    RetrievalStats,
)
from src.services.http import get_with_retry
from src.state import UniProtRecord


class UniProtClient:
    """Validate accessions against the official UniProtKB REST API."""

    SERVICE_NAME = "UniProt"
    BASE_URL = "https://rest.uniprot.org"
    USER_AGENT = "protein-supply/0.1.0"

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
        self._full_record_cache: dict[str, tuple[dict[str, Any], str]] = {}

    def get_record(self, accession: str) -> UniProtRecord:
        """Return a compact active record, resolving to its primary accession."""

        payload, source_url = self._get_payload(
            accession,
            fields="accession,id,organism_name",
        )
        return self._compact_record(payload, source_url, accession)

    def get_record_and_full(
        self,
        accession: str,
    ) -> tuple[UniProtRecord, dict[str, Any]]:
        """Fetch once and return both compact and complete UniProt records."""

        payload, source_url = self._get_payload(accession)
        return self._compact_record(payload, source_url, accession), payload

    def _compact_record(
        self,
        payload: dict[str, Any],
        source_url: str,
        accession: str,
    ) -> UniProtRecord:
        """Validate and project the fields used by workflow state."""

        entry_type = payload.get("entryType")
        primary_accession = payload.get("primaryAccession")
        entry_name = payload.get("uniProtkbId")
        organism = payload.get("organism")
        if (
            not isinstance(primary_accession, str)
            or not isinstance(entry_name, str)
            or not isinstance(entry_type, str)
            or not isinstance(organism, dict)
            or not isinstance(organism.get("scientificName"), str)
            or not isinstance(organism.get("taxonId"), int)
        ):
            raise self._invalid_response(accession, "required fields were missing")

        return {
            "primary_accession": primary_accession,
            "entry_name": entry_name,
            "entry_type": entry_type,
            "organism_name": organism["scientificName"],
            "taxon_id": organism["taxonId"],
            "source_url": source_url,
        }

    def get_full_record(self, accession: str) -> dict[str, Any]:
        """Return the complete active UniProtKB JSON record for LLM analysis."""

        cached = self._full_record_cache.get(accession.upper())
        if cached is not None:
            return cached[0]
        payload, _ = self._get_payload(accession)
        return payload

    def _get_payload(
        self,
        accession: str,
        *,
        fields: str | None = None,
    ) -> tuple[dict[str, Any], str]:
        url = f"{self.BASE_URL}/uniprotkb/{accession}.json"
        params = {"fields": fields} if fields is not None else None
        cache_key = json.dumps(
            {"accession": accession.upper(), "fields": fields},
            sort_keys=True,
            separators=(",", ":"),
        )
        if self._cache is not None:
            cached = self._cache.get("input.uniprot", cache_key)
            if (
                isinstance(cached, dict)
                and isinstance(cached.get("payload"), dict)
                and isinstance(cached.get("source_url"), str)
            ):
                if self._stats is not None:
                    self._stats.cache_hits += 1
                payload = cached["payload"]
                source_url = cached["source_url"]
                if fields is None:
                    self._full_record_cache[accession.upper()] = (
                        payload,
                        source_url,
                    )
                return payload, source_url
        if self._stats is not None:
            self._stats.live_calls += 1
        try:
            response = get_with_retry(
                self._client,
                url,
                params=params,
                headers={"Accept": "application/json", "User-Agent": self.USER_AGENT},
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
                identifier=accession,
                message=f"UniProt request failed: {exc}",
            ) from exc

        if response.status_code == 404:
            self._record_failure()
            raise DatabaseLookupError(
                kind=LookupErrorKind.NOT_FOUND,
                service=self.SERVICE_NAME,
                identifier=accession,
                message=f"UniProt accession {accession} was not found",
            )
        if response.status_code >= 400:
            self._record_failure()
            raise DatabaseLookupError(
                kind=LookupErrorKind.SERVICE_ERROR,
                service=self.SERVICE_NAME,
                identifier=accession,
                message=f"UniProt returned HTTP {response.status_code}",
            )

        try:
            payload: Any = response.json()
        except ValueError as exc:
            self._record_failure()
            raise self._invalid_response(accession, "response was not valid JSON") from exc

        if not isinstance(payload, dict):
            self._record_failure()
            raise self._invalid_response(accession, "response was not a JSON object")

        entry_type = payload.get("entryType")
        if entry_type == "Inactive":
            self._record_failure()
            reason = payload.get("inactiveReason")
            detail = ""
            if isinstance(reason, dict) and isinstance(reason.get("inactiveReasonType"), str):
                detail = f" ({reason['inactiveReasonType']})"
            raise DatabaseLookupError(
                kind=LookupErrorKind.INACTIVE_RECORD,
                service=self.SERVICE_NAME,
                identifier=accession,
                message=f"UniProt accession {accession} is inactive{detail}",
            )

        source_url = str(response.url)
        if self._cache is not None:
            self._cache.set(
                "input.uniprot",
                cache_key,
                {"payload": payload, "source_url": source_url},
                ttl_seconds=CURATED_RECORD_TTL_SECONDS,
            )
        if fields is None:
            self._full_record_cache[accession.upper()] = (payload, source_url)
        return payload, source_url

    def _record_retry(self) -> None:
        if self._stats is not None:
            self._stats.retries += 1

    def _record_failure(self) -> None:
        if self._stats is not None:
            self._stats.failed_calls += 1

    def _invalid_response(self, accession: str, detail: str) -> DatabaseLookupError:
        return DatabaseLookupError(
            kind=LookupErrorKind.INVALID_RESPONSE,
            service=self.SERVICE_NAME,
            identifier=accession,
            message=f"UniProt {detail}",
        )
