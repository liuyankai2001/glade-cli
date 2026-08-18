"""Shared HTTP behavior for external database clients."""

from collections.abc import Callable, Mapping
from email.utils import parsedate_to_datetime
import time

import httpx


DEFAULT_TIMEOUT_SECONDS = 10.0
DEFAULT_ATTEMPTS = 3
DEFAULT_RETRY_DELAY_SECONDS = 0.5
RETRYABLE_STATUS_CODES = frozenset({408, 429, 500, 502, 503, 504})


def get_with_retry(
    client: httpx.Client,
    url: str,
    *,
    params: Mapping[str, str] | None = None,
    headers: Mapping[str, str] | None = None,
    timeout: float = DEFAULT_TIMEOUT_SECONDS,
    attempts: int = DEFAULT_ATTEMPTS,
    retry_delay: float = DEFAULT_RETRY_DELAY_SECONDS,
    sleep: Callable[[float], None] = time.sleep,
    on_retry: Callable[[], None] | None = None,
) -> httpx.Response:
    """Perform a GET with bounded retries for transient transport failures."""

    if attempts < 1:
        raise ValueError("attempts must be at least 1")

    for attempt in range(attempts):
        response: httpx.Response | None = None
        try:
            response = client.get(
                url,
                params=params,
                headers=headers,
                timeout=timeout,
            )
        except httpx.RequestError:
            if attempt == attempts - 1:
                raise
        else:
            retryable_status = response.status_code in RETRYABLE_STATUS_CODES
            if not retryable_status or attempt == attempts - 1:
                return response

        if on_retry is not None:
            on_retry()
        delay = retry_delay * (2**attempt)
        if response is not None:
            delay = max(delay, min(30.0, _retry_after_seconds(response)))
        sleep(delay)

    raise RuntimeError("HTTP retry loop ended unexpectedly")


def _retry_after_seconds(response: httpx.Response) -> float:
    value = response.headers.get("Retry-After")
    if not value:
        return 0.0
    try:
        return max(0.0, float(value))
    except ValueError:
        try:
            retry_at = parsedate_to_datetime(value).timestamp()
        except (TypeError, ValueError, OverflowError):
            return 0.0
        return max(0.0, retry_at - time.time())
