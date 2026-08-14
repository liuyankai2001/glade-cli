from __future__ import annotations

import re


KEGG_COMPOUND_ID_PATTERN = re.compile(r"^C\d{5}$")


def validate_target_compound_id(value: str) -> str:
    """Validate and return a target KEGG Compound ID in ``Cxxxxx`` format."""

    target_compound_id = str(value).strip()
    if not KEGG_COMPOUND_ID_PATTERN.fullmatch(target_compound_id):
        raise ValueError(
            "target must be a KEGG Compound ID in Cxxxxx format, "
            f"got: {value!r}"
        )
    return target_compound_id
