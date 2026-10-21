"""Shared direction policy for KEGG expansion and pathway discovery."""

from __future__ import annotations

from typing import Any


DEFAULT_UNPRE_RE_REACTION = 3
REACTION_DIRECTION_POLICY_VERSION = "kegg_direction_policy.v2"


def validate_reaction_direction_mode(value: Any) -> int:
    """Reject booleans and coercions that could silently select another policy."""

    if isinstance(value, bool) or not isinstance(value, int) or value not in range(4):
        raise ValueError(
            "unpre_re_reaction must be an integer in {0, 1, 2, 3}; "
            f"got {value!r}"
        )
    return value


def reaction_direction_allowed(
    direction: str,
    thermo_direction: str,
    unpre_re_reaction: int,
) -> bool:
    """Apply progressively stricter direction policies, without fallback."""

    mode = validate_reaction_direction_mode(unpre_re_reaction)
    if mode == 0:
        return True
    if mode == 1:
        return thermo_direction != "disfavored"
    if mode == 2:
        return thermo_direction == "favored"
    return direction == "left_to_right" and thermo_direction == "favored"
