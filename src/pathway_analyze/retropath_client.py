"""Synchronous client for the loopback RetroPath Docker HTTP service."""

from __future__ import annotations

import hashlib
import ipaddress
import json
import math
import re
import tempfile
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime
from pathlib import Path, PurePosixPath
from typing import Any, Optional
from urllib.parse import quote, urlparse

import httpx

from src.pathway_analyze.retropath_input import RetroPathInputBundle
from src.pathway_analyze.retropath_models import (
    RETROPATH_RUN_STATUSES,
    RetroPathRunResult,
    RetroPathRuntimeProvenance,
)


DEFAULT_RETROPATH_SERVICE_URL = "http://127.0.0.1:8765"
DEFAULT_REQUEST_TIMEOUT_SECONDS = 30.0
DEFAULT_GET_ATTEMPTS = 3
DEFAULT_RETRY_BACKOFF_SECONDS = 0.5
DEFAULT_POLL_INTERVAL_SECONDS = 1.0
DEFAULT_WAIT_TIMEOUT_SECONDS = 3900.0

RETROPATH_CLIENT_CONTRACT_VERSION = "retropath_http.v1"
RETROPATH_CLIENT_MANIFEST_SCHEMA = "retropath_client_run.v1"
RETROPATH_CLIENT_STATE_SCHEMA = "retropath_client_state.v1"
CLIENT_MANIFEST_FILE_NAME = "run_manifest.json"
CLIENT_STATE_FILE_NAME = "client_state.json"
RAW_DIRECTORY_NAME = "raw"
SERVICE_RESULTS_FILE_NAME = "service_results.json"
SERVICE_MANIFEST_FILE_NAME = "service_run_manifest.json"

TERMINAL_STATUSES = frozenset(
    {"succeeded", "no_solution", "source_in_sink", "failed", "timed_out"}
)
CACHEABLE_STATUSES = frozenset({"succeeded", "no_solution", "source_in_sink"})
RETRYABLE_STATUS_CODES = frozenset({408, 429, 500, 502, 503, 504})
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
JOB_ID_PATTERN = re.compile(r"^rp2-[0-9a-f]{32}$")


class RetroPathClientError(RuntimeError):
    """Stable error returned for transport, protocol, or local client failures."""

    def __init__(
        self,
        code: str,
        detail: str,
        *,
        status_code: Optional[int] = None,
        job_id: Optional[str] = None,
    ) -> None:
        self.code = str(code).strip()
        self.detail = str(detail).strip()
        self.status_code = status_code
        self.job_id = job_id
        suffix = f"; job_id={job_id}" if job_id else ""
        super().__init__(f"{self.code}: {self.detail}{suffix}")


def _validated_integer(
    value: Any,
    field_name: str,
    *,
    minimum: int,
    maximum: int,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field_name} must be an integer")
    if not minimum <= value <= maximum:
        raise ValueError(f"{field_name} must be between {minimum} and {maximum}")
    return value


@dataclass(frozen=True)
class RetroPathJobParameters:
    """Parameters accepted by the P0 RetroPath service."""

    max_steps: int = 3
    topx: int = 100
    dmin: int = 2
    dmax: int = 16
    mwmax_source: int = 1000
    msc_timeout: int = 10

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "max_steps",
            _validated_integer(self.max_steps, "max_steps", minimum=1, maximum=10),
        )
        object.__setattr__(
            self,
            "topx",
            _validated_integer(self.topx, "topx", minimum=1, maximum=1000),
        )
        object.__setattr__(
            self,
            "dmin",
            _validated_integer(self.dmin, "dmin", minimum=0, maximum=16),
        )
        object.__setattr__(
            self,
            "dmax",
            _validated_integer(self.dmax, "dmax", minimum=2, maximum=16),
        )
        object.__setattr__(
            self,
            "mwmax_source",
            _validated_integer(
                self.mwmax_source,
                "mwmax_source",
                minimum=1,
                maximum=5000,
            ),
        )
        object.__setattr__(
            self,
            "msc_timeout",
            _validated_integer(
                self.msc_timeout,
                "msc_timeout",
                minimum=1,
                maximum=60,
            ),
        )
        if self.dmin > self.dmax:
            raise ValueError("dmin must be less than or equal to dmax")

    def to_dict(self) -> dict[str, int]:
        return {
            "max_steps": self.max_steps,
            "topx": self.topx,
            "dmin": self.dmin,
            "dmax": self.dmax,
            "mwmax_source": self.mwmax_source,
            "msc_timeout": self.msc_timeout,
        }

    def to_form(self) -> dict[str, str]:
        return {key: str(value) for key, value in self.to_dict().items()}

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> "RetroPathJobParameters":
        if not isinstance(payload, Mapping):
            raise ValueError("RetroPath job parameters must be an object")
        required = (
            "max_steps",
            "topx",
            "dmin",
            "dmax",
            "mwmax_source",
            "msc_timeout",
        )
        missing = [key for key in required if key not in payload]
        if missing:
            raise ValueError(f"RetroPath job parameters are missing: {missing}")
        return cls(**{key: payload[key] for key in required})


def _required_text(payload: Mapping[str, Any], key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{key} must be a non-empty string")
    return value.strip()


def _optional_text(payload: Mapping[str, Any], key: str) -> Optional[str]:
    value = payload.get(key)
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"{key} must be a string or null")
    return value.strip() or None


@dataclass(frozen=True)
class RetroPathServiceHealth:
    """Validated service health and runtime provenance."""

    service_version: str
    provenance: RetroPathRuntimeProvenance
    queue_active: int
    worker_concurrency: int

    def fingerprint_dict(self) -> dict[str, Any]:
        return {
            "service_version": self.service_version,
            **self.provenance.to_dict(),
            "worker_concurrency": self.worker_concurrency,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "ready": True,
            "service_version": self.service_version,
            "queue_active": self.queue_active,
            "worker_concurrency": self.worker_concurrency,
            **self.provenance.to_dict(),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "RetroPathServiceHealth":
        if payload.get("ready") is not True:
            errors = payload.get("errors", [])
            detail = json.dumps(errors, ensure_ascii=False) if errors else "runtime is not ready"
            raise RetroPathClientError("runtime_not_ready", detail)
        queue_active = payload.get("queue_active")
        worker_concurrency = payload.get("worker_concurrency")
        if (
            isinstance(queue_active, bool)
            or not isinstance(queue_active, int)
            or queue_active < 0
        ):
            raise ValueError("queue_active must be a non-negative integer")
        if (
            isinstance(worker_concurrency, bool)
            or not isinstance(worker_concurrency, int)
            or worker_concurrency < 1
        ):
            raise ValueError("worker_concurrency must be a positive integer")
        provenance = RetroPathRuntimeProvenance(
            wrapper_version=_required_text(payload, "wrapper_version"),
            wrapper_reported_version=_optional_text(
                payload,
                "wrapper_reported_version",
            ),
            workflow_version=_required_text(payload, "workflow_version"),
            knime_version=_required_text(payload, "knime_version"),
            rdkit_plugin_version=_required_text(payload, "rdkit_plugin_version"),
            rules_version=_required_text(payload, "rules_version"),
            rules_sha256=_required_text(payload, "rules_sha256"),
        )
        return cls(
            service_version=_required_text(payload, "service_version"),
            provenance=provenance,
            queue_active=queue_active,
            worker_concurrency=worker_concurrency,
        )


@dataclass(frozen=True)
class RetroPathJobState:
    """Validated public state of one service job."""

    job_id: str
    status: str
    parameters: RetroPathJobParameters
    created_at: str
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    return_code: Optional[int] = None
    failure_code: Optional[str] = None
    error: Optional[str] = None

    @property
    def is_terminal(self) -> bool:
        return self.status in TERMINAL_STATUSES

    def to_dict(self) -> dict[str, Any]:
        return {
            "job_id": self.job_id,
            "status": self.status,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "return_code": self.return_code,
            "failure_code": self.failure_code,
            "error": self.error,
            "parameters": self.parameters.to_dict(),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "RetroPathJobState":
        job_id = _required_text(payload, "job_id")
        if not JOB_ID_PATTERN.fullmatch(job_id):
            raise ValueError("job_id does not use the expected rp2 UUID format")
        status = _required_text(payload, "status").lower()
        if status not in RETROPATH_RUN_STATUSES:
            raise ValueError(f"unknown RetroPath job status: {status}")
        return_code = payload.get("return_code")
        if return_code is not None and (
            isinstance(return_code, bool) or not isinstance(return_code, int)
        ):
            raise ValueError("return_code must be an integer or null")
        return cls(
            job_id=job_id,
            status=status,
            parameters=RetroPathJobParameters.from_mapping(
                payload.get("parameters", {})
            ),
            created_at=_required_text(payload, "created_at"),
            started_at=_optional_text(payload, "started_at"),
            finished_at=_optional_text(payload, "finished_at"),
            return_code=return_code,
            failure_code=_optional_text(payload, "failure_code"),
            error=_optional_text(payload, "error"),
        )


@dataclass(frozen=True)
class RetroPathClientRun:
    """Final client result and the local audit files that produced it."""

    result: RetroPathRunResult
    request_fingerprint: str
    output_dir: Path
    raw_dir: Path
    run_manifest_path: Path
    client_state_path: Path
    cache_hit: bool


def _normalize_base_url(value: str) -> str:
    raw = str(value).strip().rstrip("/")
    parsed = urlparse(raw)
    if parsed.scheme != "http" or not parsed.hostname:
        raise ValueError("RetroPath service URL must be an HTTP loopback URL")
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ValueError("RetroPath service URL must not contain credentials or query data")
    if parsed.path not in {"", "/"}:
        raise ValueError("RetroPath service URL must not contain a path")
    hostname = parsed.hostname.lower()
    is_loopback = hostname == "localhost"
    if not is_loopback:
        try:
            is_loopback = ipaddress.ip_address(hostname).is_loopback
        except ValueError:
            is_loopback = False
    if not is_loopback:
        raise ValueError("RetroPath service URL must resolve explicitly to loopback")
    return raw


def _positive_finite(value: Any, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field_name} must be a positive finite number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized <= 0:
        raise ValueError(f"{field_name} must be a positive finite number")
    return normalized


def _nonnegative_finite(value: Any, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field_name} must be a non-negative finite number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0:
        raise ValueError(f"{field_name} must be a non-negative finite number")
    return normalized


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(text)
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    _atomic_write_text(
        path,
        json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        + "\n",
    )


def _read_json_object(path: Path) -> Optional[dict[str, Any]]:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _error_detail(response: httpx.Response) -> str:
    try:
        payload = response.json()
    except ValueError:
        return response.text.strip() or f"HTTP {response.status_code}"
    if isinstance(payload, Mapping) and "detail" in payload:
        detail = payload["detail"]
        if isinstance(detail, str):
            return detail
        return json.dumps(detail, ensure_ascii=False, sort_keys=True)
    return json.dumps(payload, ensure_ascii=False, sort_keys=True)


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


def _safe_relative_path(value: str, field_name: str) -> PurePosixPath:
    if (
        not isinstance(value, str)
        or not value.strip()
        or "\\" in value
        or "\x00" in value
    ):
        raise RetroPathClientError(
            "artifact_path_unsafe",
            f"{field_name} is empty or contains a backslash: {value!r}",
        )
    path = PurePosixPath(value.strip())
    if path.is_absolute() or any(
        part in {"", ".", ".."} or ":" in part for part in path.parts
    ):
        raise RetroPathClientError(
            "artifact_path_unsafe",
            f"{field_name} escapes the artifact root: {value!r}",
        )
    return path


def _local_artifact_relative(remote_path: str) -> PurePosixPath:
    remote = _safe_relative_path(remote_path, "service artifact path")
    if remote.as_posix() == "run_manifest.json":
        return PurePosixPath(RAW_DIRECTORY_NAME, SERVICE_MANIFEST_FILE_NAME)
    if remote.as_posix() in {"stdout.log", "stderr.log"}:
        return PurePosixPath(RAW_DIRECTORY_NAME, remote.name)
    if remote.parts[0] == "raw":
        if len(remote.parts) == 1:
            raise RetroPathClientError(
                "artifact_path_unsafe",
                "service artifact path points to the raw directory itself",
            )
        return PurePosixPath(RAW_DIRECTORY_NAME, *remote.parts[1:])
    return PurePosixPath(RAW_DIRECTORY_NAME, "service", *remote.parts)


def _resolved_local_path(root: Path, relative: str | PurePosixPath) -> Path:
    normalized = _safe_relative_path(str(relative), "local artifact path")
    root = root.resolve()
    candidate = (root / Path(*normalized.parts)).resolve()
    if root not in candidate.parents:
        raise RetroPathClientError(
            "artifact_path_unsafe",
            f"local artifact path escapes output directory: {relative!s}",
        )
    return candidate


class RetroPathHttpClient:
    """Submit, monitor, download, and cache jobs from the local service."""

    def __init__(
        self,
        base_url: str = DEFAULT_RETROPATH_SERVICE_URL,
        *,
        request_timeout_seconds: float = DEFAULT_REQUEST_TIMEOUT_SECONDS,
        get_attempts: int = DEFAULT_GET_ATTEMPTS,
        retry_backoff_seconds: float = DEFAULT_RETRY_BACKOFF_SECONDS,
        poll_interval_seconds: float = DEFAULT_POLL_INTERVAL_SECONDS,
        wait_timeout_seconds: float = DEFAULT_WAIT_TIMEOUT_SECONDS,
        client: Optional[httpx.Client] = None,
        sleep: Callable[[float], None] = time.sleep,
        monotonic: Callable[[], float] = time.monotonic,
        now: Callable[[], datetime] = lambda: datetime.now(UTC),
    ) -> None:
        self.base_url = _normalize_base_url(base_url)
        self.request_timeout_seconds = _positive_finite(
            request_timeout_seconds,
            "request_timeout_seconds",
        )
        if isinstance(get_attempts, bool) or not isinstance(get_attempts, int):
            raise ValueError("get_attempts must be an integer")
        if get_attempts < 1:
            raise ValueError("get_attempts must be greater than or equal to 1")
        self.get_attempts = get_attempts
        self.retry_backoff_seconds = _nonnegative_finite(
            retry_backoff_seconds,
            "retry_backoff_seconds",
        )
        self.poll_interval_seconds = _positive_finite(
            poll_interval_seconds,
            "poll_interval_seconds",
        )
        self.wait_timeout_seconds = _positive_finite(
            wait_timeout_seconds,
            "wait_timeout_seconds",
        )
        self._sleep = sleep
        self._monotonic = monotonic
        self._now = now
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=self.request_timeout_seconds,
            follow_redirects=False,
            trust_env=False,
        )

    def close(self) -> None:
        if self._owns_client:
            self._client.close()

    def __enter__(self) -> "RetroPathHttpClient":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def _url(self, path: str) -> str:
        return f"{self.base_url}{path}"

    def _sleep_before_retry(
        self,
        attempt_index: int,
        response: Optional[httpx.Response],
    ) -> None:
        delay = self.retry_backoff_seconds * (2**attempt_index)
        if response is not None:
            delay = max(delay, min(30.0, _retry_after_seconds(response)))
        if delay > 0:
            self._sleep(delay)

    def _http_error(
        self,
        response: httpx.Response,
        *,
        operation: str,
        job_id: Optional[str] = None,
    ) -> RetroPathClientError:
        detail = _error_detail(response)
        status_code = response.status_code
        if status_code == 400:
            code = "input_invalid"
        elif status_code == 404:
            code = "job_not_found"
        elif status_code == 409:
            code = "protocol_error"
        elif status_code == 503 and operation == "submit" and "queue" in detail.lower():
            code = "queue_full"
        elif status_code in RETRYABLE_STATUS_CODES or status_code >= 500:
            code = "service_unavailable"
        else:
            code = "protocol_error"
        return RetroPathClientError(
            code,
            f"{operation} failed with HTTP {status_code}: {detail}",
            status_code=status_code,
            job_id=job_id,
        )

    def _request_json(
        self,
        method: str,
        path: str,
        *,
        operation: str,
        job_id: Optional[str] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        normalized_method = method.upper()
        attempts = self.get_attempts if normalized_method == "GET" else 1
        for attempt in range(attempts):
            response: Optional[httpx.Response] = None
            try:
                response = self._client.request(
                    normalized_method,
                    self._url(path),
                    timeout=self.request_timeout_seconds,
                    **kwargs,
                )
            except httpx.RequestError as exc:
                if normalized_method == "POST":
                    raise RetroPathClientError(
                        "submission_uncertain",
                        f"submission transport failed after the request may have been sent: {exc}",
                    ) from exc
                if attempt == attempts - 1:
                    raise RetroPathClientError(
                        "service_unavailable",
                        f"{operation} transport failed: {exc}",
                        job_id=job_id,
                    ) from exc
                self._sleep_before_retry(attempt, None)
                continue

            if (
                normalized_method == "GET"
                and response.status_code in RETRYABLE_STATUS_CODES
                and attempt < attempts - 1
            ):
                self._sleep_before_retry(attempt, response)
                continue
            if response.status_code >= 400:
                raise self._http_error(
                    response,
                    operation=operation,
                    job_id=job_id,
                )
            try:
                payload = response.json()
            except ValueError as exc:
                raise RetroPathClientError(
                    "protocol_error",
                    f"{operation} returned invalid JSON",
                    status_code=response.status_code,
                    job_id=job_id,
                ) from exc
            if not isinstance(payload, dict):
                raise RetroPathClientError(
                    "protocol_error",
                    f"{operation} JSON root must be an object",
                    status_code=response.status_code,
                    job_id=job_id,
                )
            return payload
        raise RuntimeError("HTTP request retry loop ended unexpectedly")

    def health(self) -> RetroPathServiceHealth:
        payload = self._request_json("GET", "/health", operation="health")
        try:
            return RetroPathServiceHealth.from_payload(payload)
        except RetroPathClientError:
            raise
        except ValueError as exc:
            raise RetroPathClientError(
                "protocol_error",
                f"health payload is invalid: {exc}",
            ) from exc

    def submit(
        self,
        input_bundle: RetroPathInputBundle,
        parameters: RetroPathJobParameters,
    ) -> RetroPathJobState:
        if not isinstance(input_bundle, RetroPathInputBundle):
            raise ValueError("input_bundle must be a RetroPathInputBundle")
        if not isinstance(parameters, RetroPathJobParameters):
            raise ValueError("parameters must be RetroPathJobParameters")
        try:
            with input_bundle.target_source_path.open("rb") as source_handle:
                with input_bundle.chassis_sink_path.open("rb") as sink_handle:
                    payload = self._request_json(
                        "POST",
                        "/v1/jobs",
                        operation="submit",
                        files={
                            "source_file": (
                                "source.csv",
                                source_handle,
                                "text/csv",
                            ),
                            "sink_file": ("sink.csv", sink_handle, "text/csv"),
                        },
                        data=parameters.to_form(),
                    )
        except OSError as exc:
            raise RetroPathClientError(
                "input_modified",
                f"cannot read P2 input files: {exc}",
            ) from exc
        state = self._parse_job(payload, operation="submit")
        if state.parameters != parameters:
            raise RetroPathClientError(
                "protocol_error",
                "submitted job parameters do not match the request",
                job_id=state.job_id,
            )
        return state

    def _parse_job(
        self,
        payload: Mapping[str, Any],
        *,
        operation: str,
        expected_job_id: Optional[str] = None,
    ) -> RetroPathJobState:
        try:
            state = RetroPathJobState.from_payload(payload)
        except ValueError as exc:
            raise RetroPathClientError(
                "protocol_error",
                f"{operation} job payload is invalid: {exc}",
                job_id=expected_job_id,
            ) from exc
        if expected_job_id is not None and state.job_id != expected_job_id:
            raise RetroPathClientError(
                "protocol_error",
                f"{operation} returned job_id {state.job_id}, expected {expected_job_id}",
                job_id=expected_job_id,
            )
        return state

    def get_job(self, job_id: str) -> RetroPathJobState:
        normalized_job_id = str(job_id).strip()
        if not JOB_ID_PATTERN.fullmatch(normalized_job_id):
            raise ValueError("job_id does not use the expected rp2 UUID format")
        payload = self._request_json(
            "GET",
            f"/v1/jobs/{quote(normalized_job_id, safe='')}",
            operation="get_job",
            job_id=normalized_job_id,
        )
        return self._parse_job(
            payload,
            operation="get_job",
            expected_job_id=normalized_job_id,
        )

    def wait_for_job(
        self,
        job_id: str,
        *,
        timeout_seconds: Optional[float] = None,
        poll_interval_seconds: Optional[float] = None,
        on_status: Optional[Callable[[RetroPathJobState], None]] = None,
    ) -> RetroPathJobState:
        timeout = (
            self.wait_timeout_seconds
            if timeout_seconds is None
            else _positive_finite(timeout_seconds, "timeout_seconds")
        )
        interval = (
            self.poll_interval_seconds
            if poll_interval_seconds is None
            else _positive_finite(poll_interval_seconds, "poll_interval_seconds")
        )
        deadline = self._monotonic() + timeout
        last_status: Optional[str] = None
        while True:
            state = self.get_job(job_id)
            if on_status is not None and state.status != last_status:
                on_status(state)
            last_status = state.status
            if state.is_terminal:
                return state
            remaining = deadline - self._monotonic()
            if remaining <= 0:
                raise RetroPathClientError(
                    "client_poll_timeout",
                    f"job did not reach a terminal state within {timeout:g} seconds",
                    job_id=job_id,
                )
            self._sleep(min(interval, remaining))

    def get_results(self, job_id: str) -> dict[str, Any]:
        payload = self._request_json(
            "GET",
            f"/v1/jobs/{quote(job_id, safe='')}/results",
            operation="get_results",
            job_id=job_id,
        )
        state = self._parse_job(
            payload,
            operation="get_results",
            expected_job_id=job_id,
        )
        if not state.is_terminal:
            raise RetroPathClientError(
                "protocol_error",
                "results endpoint returned a non-terminal job",
                job_id=job_id,
            )
        artifacts = payload.get("artifacts")
        if (
            isinstance(artifacts, (str, bytes))
            or not isinstance(artifacts, Sequence)
            or not all(isinstance(item, str) and item.strip() for item in artifacts)
        ):
            raise RetroPathClientError(
                "protocol_error",
                "results artifacts must be a list of non-empty strings",
                job_id=job_id,
            )
        if len(set(artifacts)) != len(artifacts):
            raise RetroPathClientError(
                "protocol_error",
                "results artifacts contain duplicate paths",
                job_id=job_id,
            )
        manifest = payload.get("manifest")
        if manifest is not None and not isinstance(manifest, Mapping):
            raise RetroPathClientError(
                "protocol_error",
                "results manifest must be an object or null",
                job_id=job_id,
            )
        return payload

    def _download_one(
        self,
        job_id: str,
        remote_path: str,
        local_path: Path,
    ) -> str:
        encoded_path = quote(remote_path, safe="/")
        url = self._url(
            f"/v1/jobs/{quote(job_id, safe='')}/artifacts/{encoded_path}"
        )
        for attempt in range(self.get_attempts):
            response: Optional[httpx.Response] = None
            temporary_path: Optional[Path] = None
            try:
                with self._client.stream(
                    "GET",
                    url,
                    timeout=self.request_timeout_seconds,
                ) as response:
                    if (
                        response.status_code in RETRYABLE_STATUS_CODES
                        and attempt < self.get_attempts - 1
                    ):
                        self._sleep_before_retry(attempt, response)
                        continue
                    if response.status_code >= 400:
                        error = self._http_error(
                            response,
                            operation="download_artifact",
                            job_id=job_id,
                        )
                        raise RetroPathClientError(
                            "artifact_download_failed",
                            error.detail,
                            status_code=error.status_code,
                            job_id=job_id,
                        )
                    local_path.parent.mkdir(parents=True, exist_ok=True)
                    digest = hashlib.sha256()
                    with tempfile.NamedTemporaryFile(
                        mode="wb",
                        dir=local_path.parent,
                        prefix=f".{local_path.name}.",
                        suffix=".tmp",
                        delete=False,
                    ) as handle:
                        temporary_path = Path(handle.name)
                        for chunk in response.iter_bytes(1024 * 1024):
                            digest.update(chunk)
                            handle.write(chunk)
                    temporary_path.replace(local_path)
                    temporary_path = None
                    expected_hash = digest.hexdigest()
                    if _sha256_file(local_path) != expected_hash:
                        raise RetroPathClientError(
                            "artifact_checksum_mismatch",
                            f"downloaded artifact failed local SHA-256 verification: {remote_path}",
                            job_id=job_id,
                        )
                    return expected_hash
            except RetroPathClientError:
                raise
            except (httpx.HTTPError, OSError) as exc:
                if attempt == self.get_attempts - 1:
                    raise RetroPathClientError(
                        "artifact_download_failed",
                        f"failed to download {remote_path}: {exc}",
                        job_id=job_id,
                    ) from exc
                self._sleep_before_retry(attempt, response)
            finally:
                if temporary_path is not None:
                    temporary_path.unlink(missing_ok=True)
        raise RuntimeError("artifact retry loop ended unexpectedly")

    def download_artifacts(
        self,
        job_id: str,
        artifact_paths: Sequence[str],
        output_dir: str | Path,
    ) -> tuple[dict[str, str], dict[str, str]]:
        root = Path(output_dir).expanduser().resolve()
        mapping: dict[str, str] = {}
        hashes: dict[str, str] = {}
        reserved = PurePosixPath(RAW_DIRECTORY_NAME, SERVICE_RESULTS_FILE_NAME)
        local_paths: set[str] = {reserved.as_posix()}
        for remote_path in sorted(set(artifact_paths)):
            local_relative = _local_artifact_relative(remote_path)
            local_key = local_relative.as_posix()
            if local_key in local_paths:
                raise RetroPathClientError(
                    "artifact_path_unsafe",
                    f"multiple service artifacts map to {local_key}",
                    job_id=job_id,
                )
            local_paths.add(local_key)
            local_path = _resolved_local_path(root, local_relative)
            hashes[local_key] = self._download_one(
                job_id,
                remote_path,
                local_path,
            )
            mapping[remote_path] = local_key
        return mapping, hashes

    def _validate_input_hashes(
        self,
        input_bundle: RetroPathInputBundle,
    ) -> tuple[str, str]:
        for path, field_name in (
            (input_bundle.target_source_path, "target_source.csv"),
            (input_bundle.chassis_sink_path, "chassis_sink.csv"),
        ):
            if not path.is_file():
                raise RetroPathClientError(
                    "input_modified",
                    f"P2 input file is missing: {field_name}: {path}",
                )
        source_sha256 = _sha256_file(input_bundle.target_source_path)
        sink_sha256 = _sha256_file(input_bundle.chassis_sink_path)
        if source_sha256 != input_bundle.target_source_sha256:
            raise RetroPathClientError(
                "input_modified",
                "target_source.csv no longer matches its P2 SHA-256",
            )
        if sink_sha256 != input_bundle.chassis_sink_sha256:
            raise RetroPathClientError(
                "input_modified",
                "chassis_sink.csv no longer matches its P2 SHA-256",
            )
        return source_sha256, sink_sha256

    def _request_fingerprint(
        self,
        input_bundle: RetroPathInputBundle,
        parameters: RetroPathJobParameters,
        health: RetroPathServiceHealth,
        source_sha256: str,
        sink_sha256: str,
    ) -> str:
        return _canonical_sha256(
            {
                "client_contract_version": RETROPATH_CLIENT_CONTRACT_VERSION,
                "service_url": self.base_url,
                "health": health.fingerprint_dict(),
                "target_compound_id": input_bundle.target_compound.compound_id,
                "expansion_depth": input_bundle.expansion_depth,
                "source_sha256": source_sha256,
                "sink_sha256": sink_sha256,
                "parameters": parameters.to_dict(),
            }
        )

    def _state_payload(
        self,
        fingerprint: str,
        phase: str,
        *,
        parameters: RetroPathJobParameters,
        job: Optional[RetroPathJobState] = None,
        error: Optional[RetroPathClientError] = None,
    ) -> dict[str, Any]:
        return {
            "schema_version": RETROPATH_CLIENT_STATE_SCHEMA,
            "updated_at": self._now().isoformat(),
            "request_fingerprint": fingerprint,
            "phase": phase,
            "parameters": parameters.to_dict(),
            "job": job.to_dict() if job is not None else None,
            "error": (
                None
                if error is None
                else {
                    "code": error.code,
                    "detail": error.detail,
                    "status_code": error.status_code,
                    "job_id": error.job_id,
                }
            ),
        }

    def _write_state(
        self,
        path: Path,
        fingerprint: str,
        phase: str,
        *,
        parameters: RetroPathJobParameters,
        job: Optional[RetroPathJobState] = None,
        error: Optional[RetroPathClientError] = None,
    ) -> None:
        _atomic_write_json(
            path,
            self._state_payload(
                fingerprint,
                phase,
                parameters=parameters,
                job=job,
                error=error,
            ),
        )

    def _try_cached_run(
        self,
        output_dir: Path,
        fingerprint: str,
    ) -> Optional[RetroPathClientRun]:
        manifest_path = output_dir / CLIENT_MANIFEST_FILE_NAME
        state_path = output_dir / CLIENT_STATE_FILE_NAME
        payload = _read_json_object(manifest_path)
        if (
            payload is None
            or payload.get("schema_version") != RETROPATH_CLIENT_MANIFEST_SCHEMA
            or payload.get("request_fingerprint") != fingerprint
        ):
            return None
        try:
            result = RetroPathRunResult.from_dict(payload.get("run_result", {}))
        except ValueError:
            return None
        if result.status not in CACHEABLE_STATUSES:
            return None
        artifact_hashes = payload.get("artifact_sha256")
        if not isinstance(artifact_hashes, Mapping):
            return None
        if set(result.artifacts) != set(artifact_hashes):
            return None
        for relative, expected_hash in artifact_hashes.items():
            if (
                not isinstance(relative, str)
                or not isinstance(expected_hash, str)
                or not SHA256_PATTERN.fullmatch(expected_hash)
            ):
                return None
            try:
                path = _resolved_local_path(output_dir, relative)
            except RetroPathClientError:
                return None
            if not path.is_file() or _sha256_file(path) != expected_hash:
                return None
        return RetroPathClientRun(
            result=result,
            request_fingerprint=fingerprint,
            output_dir=output_dir,
            raw_dir=output_dir / RAW_DIRECTORY_NAME,
            run_manifest_path=manifest_path,
            client_state_path=state_path,
            cache_hit=True,
        )

    def _resumable_job(
        self,
        state_path: Path,
        fingerprint: str,
        parameters: RetroPathJobParameters,
    ) -> Optional[RetroPathJobState]:
        payload = _read_json_object(state_path)
        if (
            payload is None
            or payload.get("schema_version") != RETROPATH_CLIENT_STATE_SCHEMA
            or payload.get("request_fingerprint") != fingerprint
        ):
            return None
        phase = payload.get("phase")
        if phase == "submission_uncertain":
            raise RetroPathClientError(
                "submission_uncertain",
                "a previous submission may have reached the service; use "
                "force=True only after checking the service",
            )
        job_payload = payload.get("job")
        if not isinstance(job_payload, Mapping):
            return None
        try:
            saved = RetroPathJobState.from_payload(job_payload)
        except ValueError:
            return None
        if saved.parameters != parameters:
            return None
        if phase == "completed" and saved.status not in CACHEABLE_STATUSES:
            return None
        try:
            return self.get_job(saved.job_id)
        except RetroPathClientError as exc:
            if exc.code == "job_not_found":
                return None
            raise

    def _validate_service_manifest(
        self,
        manifest: Optional[Mapping[str, Any]],
        *,
        job: RetroPathJobState,
        health: RetroPathServiceHealth,
        parameters: RetroPathJobParameters,
        source_sha256: str,
        sink_sha256: str,
        result_artifacts: Sequence[str],
    ) -> None:
        if manifest is None:
            if job.status in CACHEABLE_STATUSES:
                raise RetroPathClientError(
                    "protocol_error",
                    f"terminal status {job.status} is missing the service run manifest",
                    job_id=job.job_id,
                )
            return
        manifest_schema = manifest.get("schema_version")
        if manifest_schema not in {1, 2}:
            raise RetroPathClientError(
                "protocol_error",
                "service manifest schema_version is not supported",
                job_id=job.job_id,
            )
        expected_scalars = {
            "job_id": job.job_id,
            "status": job.status,
            "return_code": job.return_code,
        }
        if manifest_schema == 2:
            expected_scalars["failure_code"] = job.failure_code
        for key, expected in expected_scalars.items():
            if manifest.get(key) != expected:
                raise RetroPathClientError(
                    "protocol_error",
                    f"service manifest {key} does not match terminal job state",
                    job_id=job.job_id,
                )
        if manifest.get("parameters") != parameters.to_dict():
            raise RetroPathClientError(
                "protocol_error",
                "service manifest parameters do not match submission",
                job_id=job.job_id,
            )
        input_hashes = manifest.get("input_sha256")
        if not isinstance(input_hashes, Mapping) or input_hashes != {
            "source.csv": source_sha256,
            "sink.csv": sink_sha256,
        }:
            raise RetroPathClientError(
                "protocol_error",
                "service manifest input SHA-256 values do not match P2 inputs",
                job_id=job.job_id,
            )
        versions = manifest.get("versions")
        expected_versions = {
            "retropath2_wrapper": health.provenance.wrapper_version,
            "retropath2_wrapper_reported": (
                health.provenance.wrapper_reported_version
            ),
            "workflow": health.provenance.workflow_version,
            "knime": health.provenance.knime_version,
            "knime_rdkit_nodes": health.provenance.rdkit_plugin_version,
            "rules": health.provenance.rules_version,
        }
        if versions != expected_versions:
            raise RetroPathClientError(
                "protocol_error",
                "service manifest versions do not match the health response",
                job_id=job.job_id,
            )
        if manifest.get("rules_sha256") != health.provenance.rules_sha256:
            raise RetroPathClientError(
                "protocol_error",
                "service manifest rules SHA-256 does not match health",
                job_id=job.job_id,
            )
        if manifest_schema == 2 and not isinstance(
            manifest.get("resource_telemetry"),
            Mapping,
        ):
            raise RetroPathClientError(
                "protocol_error",
                "service manifest resource_telemetry is missing",
                job_id=job.job_id,
            )
        manifest_artifacts = manifest.get("artifacts")
        if (
            isinstance(manifest_artifacts, (str, bytes))
            or not isinstance(manifest_artifacts, Sequence)
            or not all(
                isinstance(item, str) and item.strip()
                for item in manifest_artifacts
            )
            or list(result_artifacts)
            != ["run_manifest.json", *list(manifest_artifacts)]
        ):
            raise RetroPathClientError(
                "protocol_error",
                "results artifact list does not match the service manifest",
                job_id=job.job_id,
            )

    def _build_run_result(
        self,
        input_bundle: RetroPathInputBundle,
        job: RetroPathJobState,
        health: RetroPathServiceHealth,
        parameters: RetroPathJobParameters,
        fingerprint: str,
        source_sha256: str,
        sink_sha256: str,
        local_artifacts: Sequence[str],
    ) -> RetroPathRunResult:
        errors: tuple[str, ...] = tuple()
        if job.error:
            errors = (job.error,)
        elif job.status == "failed":
            errors = ("RetroPath service reported a failed job without error text",)
        audit_parameters: dict[str, Any] = {
            **parameters.to_dict(),
            "service_url": self.base_url,
            "request_fingerprint": fingerprint,
            "source_sha256": source_sha256,
            "sink_sha256": sink_sha256,
            "target_compound_id": input_bundle.target_compound.compound_id,
            "expansion_depth": input_bundle.expansion_depth,
        }
        return RetroPathRunResult(
            job_id=job.job_id,
            status=job.status,
            return_code=job.return_code,
            failure_code=job.failure_code,
            provenance=health.provenance,
            parameters=tuple(audit_parameters.items()),
            artifacts=tuple(local_artifacts),
            errors=errors,
        )

    def run(
        self,
        input_bundle: RetroPathInputBundle,
        output_dir: str | Path,
        *,
        parameters: Optional[RetroPathJobParameters] = None,
        force: bool = False,
        timeout_seconds: Optional[float] = None,
        on_status: Optional[Callable[[RetroPathJobState], None]] = None,
    ) -> RetroPathClientRun:
        if not isinstance(input_bundle, RetroPathInputBundle):
            raise ValueError("input_bundle must be a RetroPathInputBundle")
        job_parameters = parameters or RetroPathJobParameters()
        if not isinstance(job_parameters, RetroPathJobParameters):
            raise ValueError("parameters must be RetroPathJobParameters or null")
        if not isinstance(force, bool):
            raise ValueError("force must be a boolean")
        resolved_output_dir = Path(output_dir).expanduser().resolve()
        resolved_output_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = resolved_output_dir / CLIENT_MANIFEST_FILE_NAME
        state_path = resolved_output_dir / CLIENT_STATE_FILE_NAME
        raw_dir = resolved_output_dir / RAW_DIRECTORY_NAME

        source_sha256, sink_sha256 = self._validate_input_hashes(input_bundle)
        health = self.health()
        fingerprint = self._request_fingerprint(
            input_bundle,
            job_parameters,
            health,
            source_sha256,
            sink_sha256,
        )
        if not force:
            cached = self._try_cached_run(resolved_output_dir, fingerprint)
            if cached is not None:
                return cached

        job: Optional[RetroPathJobState] = None
        if not force:
            job = self._resumable_job(state_path, fingerprint, job_parameters)
        if job is None:
            self._write_state(
                state_path,
                fingerprint,
                "submitting",
                parameters=job_parameters,
            )
            try:
                job = self.submit(input_bundle, job_parameters)
            except RetroPathClientError as exc:
                phase = (
                    "submission_uncertain"
                    if exc.code == "submission_uncertain"
                    else "submission_failed"
                )
                self._write_state(
                    state_path,
                    fingerprint,
                    phase,
                    parameters=job_parameters,
                    error=exc,
                )
                raise
            self._write_state(
                state_path,
                fingerprint,
                "submitted",
                parameters=job_parameters,
                job=job,
            )

        if job.parameters != job_parameters:
            raise RetroPathClientError(
                "protocol_error",
                "resumed job parameters do not match the current request",
                job_id=job.job_id,
            )
        if not job.is_terminal:
            self._write_state(
                state_path,
                fingerprint,
                "polling",
                parameters=job_parameters,
                job=job,
            )
            latest_job = job

            def record_status(state: RetroPathJobState) -> None:
                nonlocal latest_job
                latest_job = state
                self._write_state(
                    state_path,
                    fingerprint,
                    "polling",
                    parameters=job_parameters,
                    job=state,
                )
                if on_status is not None:
                    on_status(state)

            try:
                job = self.wait_for_job(
                    job.job_id,
                    timeout_seconds=timeout_seconds,
                    on_status=record_status,
                )
            except RetroPathClientError as exc:
                job = latest_job
                phase = (
                    "client_poll_timeout"
                    if exc.code == "client_poll_timeout"
                    else "poll_interrupted"
                )
                self._write_state(
                    state_path,
                    fingerprint,
                    phase,
                    parameters=job_parameters,
                    job=job,
                    error=exc,
                )
                raise
        self._write_state(
            state_path,
            fingerprint,
            "terminal",
            parameters=job_parameters,
            job=job,
        )

        try:
            service_results = self.get_results(job.job_id)
            results_job = self._parse_job(
                service_results,
                operation="get_results",
                expected_job_id=job.job_id,
            )
            if results_job != job:
                raise RetroPathClientError(
                    "protocol_error",
                    "results job state does not match the final polled state",
                    job_id=job.job_id,
                )
            service_manifest = service_results.get("manifest")
            self._validate_service_manifest(
                service_manifest,
                job=job,
                health=health,
                parameters=job_parameters,
                source_sha256=source_sha256,
                sink_sha256=sink_sha256,
                result_artifacts=service_results["artifacts"],
            )
            raw_dir.mkdir(parents=True, exist_ok=True)
            service_results_path = raw_dir / SERVICE_RESULTS_FILE_NAME
            service_results_text = json.dumps(
                service_results,
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
                allow_nan=False,
            ) + "\n"
            _atomic_write_text(service_results_path, service_results_text)
            service_results_relative = PurePosixPath(
                RAW_DIRECTORY_NAME,
                SERVICE_RESULTS_FILE_NAME,
            ).as_posix()
            artifact_hashes = {
                service_results_relative: hashlib.sha256(
                    service_results_text.encode("utf-8")
                ).hexdigest()
            }
            artifact_mapping, downloaded_hashes = self.download_artifacts(
                job.job_id,
                tuple(service_results["artifacts"]),
                resolved_output_dir,
            )
            artifact_hashes.update(downloaded_hashes)
        except RetroPathClientError as exc:
            self._write_state(
                state_path,
                fingerprint,
                "result_failed",
                parameters=job_parameters,
                job=job,
                error=exc,
            )
            raise

        run_result = self._build_run_result(
            input_bundle,
            job,
            health,
            job_parameters,
            fingerprint,
            source_sha256,
            sink_sha256,
            tuple(sorted(artifact_hashes)),
        )
        manifest_payload = {
            "schema_version": RETROPATH_CLIENT_MANIFEST_SCHEMA,
            "client_contract_version": RETROPATH_CLIENT_CONTRACT_VERSION,
            "recorded_at": self._now().isoformat(),
            "request_fingerprint": fingerprint,
            "cacheable": job.status in CACHEABLE_STATUSES,
            "service_url": self.base_url,
            "target_compound_id": input_bundle.target_compound.compound_id,
            "expansion_depth": input_bundle.expansion_depth,
            "input_sha256": {
                "target_source.csv": source_sha256,
                "chassis_sink.csv": sink_sha256,
            },
            "parameters": job_parameters.to_dict(),
            "health": health.to_dict(),
            "job": job.to_dict(),
            "artifact_mapping": dict(sorted(artifact_mapping.items())),
            "artifact_sha256": dict(sorted(artifact_hashes.items())),
            "run_result": run_result.to_dict(),
        }
        _atomic_write_json(manifest_path, manifest_payload)
        self._write_state(
            state_path,
            fingerprint,
            "completed",
            parameters=job_parameters,
            job=job,
        )
        return RetroPathClientRun(
            result=run_result,
            request_fingerprint=fingerprint,
            output_dir=resolved_output_dir,
            raw_dir=raw_dir,
            run_manifest_path=manifest_path,
            client_state_path=state_path,
            cache_hit=False,
        )


__all__ = [
    "CACHEABLE_STATUSES",
    "CLIENT_MANIFEST_FILE_NAME",
    "CLIENT_STATE_FILE_NAME",
    "DEFAULT_GET_ATTEMPTS",
    "DEFAULT_POLL_INTERVAL_SECONDS",
    "DEFAULT_REQUEST_TIMEOUT_SECONDS",
    "DEFAULT_RETROPATH_SERVICE_URL",
    "DEFAULT_RETRY_BACKOFF_SECONDS",
    "DEFAULT_WAIT_TIMEOUT_SECONDS",
    "RETROPATH_CLIENT_CONTRACT_VERSION",
    "RETROPATH_CLIENT_MANIFEST_SCHEMA",
    "RETROPATH_CLIENT_STATE_SCHEMA",
    "RetroPathClientError",
    "RetroPathClientRun",
    "RetroPathHttpClient",
    "RetroPathJobParameters",
    "RetroPathJobState",
    "RetroPathServiceHealth",
    "TERMINAL_STATUSES",
]
