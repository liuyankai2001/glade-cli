"""User-facing design failures shared by the local web API and DNA service."""

from __future__ import annotations


class DesignError(ValueError):
    def __init__(
        self,
        message: str,
        *,
        code: str = "invalid_design",
        issues: list[dict] | None = None,
    ):
        super().__init__(message)
        self.code = code
        self.issues = issues or []
