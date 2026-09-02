"""Shared amino-acid sequence quality checks for main-enzyme candidates."""

from __future__ import annotations

import re
from dataclasses import dataclass


STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
UNSUPPORTED_AMINO_ACID_REASON_CODE = "unsupported_amino_acid_characters"


@dataclass(frozen=True)
class ProteinSequenceQuality:
    """Normalized sequence plus unsupported-character locations."""

    normalized_sequence: str
    unsupported_positions: dict[str, tuple[int, ...]]

    @property
    def is_constructable(self) -> bool:
        return bool(self.normalized_sequence) and not self.unsupported_positions

    @property
    def unsupported_characters(self) -> tuple[str, ...]:
        return tuple(sorted(self.unsupported_positions))

    def rejection_reason(self) -> str:
        details = " ".join(
            f"{character}[count={len(positions)},positions="
            f"{','.join(map(str, positions))}]"
            for character, positions in sorted(self.unsupported_positions.items())
        )
        return f"{UNSUPPORTED_AMINO_ACID_REASON_CODE}: {details}"


def analyze_protein_sequence(sequence: object) -> ProteinSequenceQuality:
    """Return constructability information using 1-based normalized positions."""

    normalized = re.sub(r"\s+", "", str(sequence or "")).upper()
    unsupported: dict[str, list[int]] = {}
    for position, character in enumerate(normalized, start=1):
        if character not in STANDARD_AMINO_ACIDS:
            unsupported.setdefault(character, []).append(position)
    return ProteinSequenceQuality(
        normalized_sequence=normalized,
        unsupported_positions={
            character: tuple(positions)
            for character, positions in unsupported.items()
        },
    )


def unsupported_amino_acid_reason(sequence: object) -> str:
    """Return a stable rejection reason, or an empty string when acceptable."""

    quality = analyze_protein_sequence(sequence)
    if not quality.unsupported_positions:
        return ""
    return quality.rejection_reason()


__all__ = [
    "ProteinSequenceQuality",
    "STANDARD_AMINO_ACIDS",
    "UNSUPPORTED_AMINO_ACID_REASON_CODE",
    "analyze_protein_sequence",
    "unsupported_amino_acid_reason",
]
