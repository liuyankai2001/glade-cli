"""Canonical reaction-scope helpers shared by protein researchers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import re
from typing import Any, Self

from pydantic import BaseModel, Field, field_validator, model_validator


_KEGG_REACTION_ID_PATTERN = re.compile(r"^R\d{5}$")


def normalize_reaction_ids(values: Sequence[Any]) -> list[str]:
    """Return a deterministic, non-empty KEGG reaction-ID set."""

    normalized: list[str] = []
    for value in values:
        if not isinstance(value, str):
            raise ValueError("reaction_ids must contain strings")
        reaction_id = value.strip().upper()
        if not _KEGG_REACTION_ID_PATTERN.fullmatch(reaction_id):
            raise ValueError(f"invalid KEGG reaction ID: {value}")
        normalized.append(reaction_id)
    if not normalized:
        raise ValueError("reaction_ids must not be empty")
    if len(normalized) != len(set(normalized)):
        raise ValueError("reaction_ids must be unique")
    return sorted(normalized)


class ReactionScopedModel(BaseModel):
    """Base model that emits plural reaction scope only.

    Old single-step cached objects containing ``reaction_id`` remain readable,
    but serialization always uses the canonical ``reaction_ids`` field.
    """

    reaction_ids: list[str] = Field(min_length=1)

    @model_validator(mode="before")
    @classmethod
    def migrate_legacy_reaction_id(cls, value: Any) -> Any:
        if not isinstance(value, Mapping):
            return value
        data = dict(value)
        has_legacy = "reaction_id" in data
        has_plural = "reaction_ids" in data
        if has_legacy and has_plural:
            raise ValueError(
                "provide reaction_ids only; reaction_id is a legacy input field"
            )
        if has_legacy:
            data["reaction_ids"] = [data.pop("reaction_id")]
        return data

    @field_validator("reaction_ids")
    @classmethod
    def validate_reaction_ids(cls, value: list[str]) -> list[str]:
        return normalize_reaction_ids(value)


def context_reactions(context: Any) -> tuple[Any, ...]:
    """Return every official reaction identity represented by a context."""

    reaction_steps = getattr(context, "reaction_steps", None)
    if isinstance(reaction_steps, Sequence) and not isinstance(
        reaction_steps,
        (str, bytes),
    ):
        reactions = tuple(
            reaction
            for step in reaction_steps
            if (reaction := getattr(step, "reaction", None)) is not None
        )
        if reactions:
            return reactions
    reaction = getattr(context, "reaction", None)
    if reaction is not None:
        return (reaction,)
    raise ValueError("research context does not contain reaction identities")


def context_reaction_ids(context: Any) -> list[str]:
    """Return the canonical reaction scope represented by a context."""

    values = [
        getattr(reaction, "reaction_id", None)
        for reaction in context_reactions(context)
    ]
    return normalize_reaction_ids(list(dict.fromkeys(values)))


def require_exact_reaction_scope(
    actual: Sequence[str],
    expected: Sequence[str],
    *,
    label: str,
) -> None:
    """Reject a report that drops or invents requested reactions."""

    if normalize_reaction_ids(actual) != normalize_reaction_ids(expected):
        raise ValueError(f"{label} changed the requested reaction scope")


def require_reaction_subset(
    actual: Sequence[str],
    allowed: Sequence[str],
    *,
    label: str,
) -> None:
    """Reject evidence whose reaction scope escapes its parent report."""

    actual_ids = set(normalize_reaction_ids(actual))
    allowed_ids = set(normalize_reaction_ids(allowed))
    unknown = sorted(actual_ids - allowed_ids)
    if unknown:
        raise ValueError(
            f"{label} references reactions outside the report scope: "
            + ", ".join(unknown)
        )
