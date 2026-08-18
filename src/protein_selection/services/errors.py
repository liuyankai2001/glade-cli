"""Typed failures raised by biological database clients."""

from enum import StrEnum


class LookupErrorKind(StrEnum):
    """Stable error categories consumed by workflow nodes."""

    NOT_FOUND = "not_found"
    INACTIVE_RECORD = "inactive_record"
    SERVICE_ERROR = "service_error"
    INVALID_RESPONSE = "invalid_response"


class DatabaseLookupError(Exception):
    """A database lookup failed in an expected, classifiable way."""

    def __init__(
        self,
        *,
        kind: LookupErrorKind,
        service: str,
        identifier: str,
        message: str,
    ) -> None:
        super().__init__(message)
        self.kind = kind
        self.service = service
        self.identifier = identifier
