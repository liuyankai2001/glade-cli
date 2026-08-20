"""Pure sequence and coordinate helpers used by assembly planners."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio import Restriction, SeqIO
from Bio.Seq import Seq

from src.final_assemble_plan.config import (
    CURATED_RESTRICTION_ENZYMES,
    ENZYME_PRIORITIES,
    MAX_HOMOLOGY_ARM_GC_PERCENT,
    MAX_HOMOLOGY_ARM_HOMOPOLYMER,
    MIN_HOMOLOGY_ARM_GC_PERCENT,
)


DNA_ALPHABET = frozenset("ACGT")


def stable_json_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def resolve_project_file(project_output: Path, value: Any) -> Path:
    text = str(value or "").strip()
    if not text:
        raise ValueError("sequence file path is empty")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = project_output / path
    path = path.resolve()
    try:
        path.relative_to(project_output.resolve())
    except ValueError as exc:
        raise ValueError(f"sequence file is outside the project output: {path}") from exc
    return path


def read_genbank(path: Path) -> tuple[Any, str, bytes]:
    if not path.is_file():
        raise FileNotFoundError(path)
    content = path.read_bytes()
    try:
        record = SeqIO.read(path, "genbank")
    except Exception as exc:
        raise ValueError(f"could not parse GenBank file: {path}") from exc
    sequence = str(record.seq).upper()
    if not sequence or set(sequence) - DNA_ALPHABET:
        raise ValueError(f"GenBank sequence must contain only A/C/G/T: {path}")
    return record, sequence, content


def feature_payloads(record: Any) -> tuple[dict[str, Any], ...]:
    result: list[dict[str, Any]] = []
    for feature in record.features:
        labels: list[str] = []
        for key in ("label", "gene", "product", "note"):
            labels.extend(str(item) for item in feature.qualifiers.get(key, []))
        result.append(
            {
                "type": str(feature.type),
                "start_bp": int(feature.location.start) + 1,
                "end_bp": int(feature.location.end),
                "strand": int(feature.location.strand or 0),
                "labels": labels,
            }
        )
    return tuple(result)


def relative_project_path(project_output: Path, path: Path) -> str:
    return path.resolve().relative_to(project_output.resolve()).as_posix()


def max_homopolymer(sequence: str) -> int:
    if not sequence:
        return 0
    return max(len(match.group(0)) for match in re.finditer(r"([ACGT])\1*", sequence))


def gc_percent(sequence: str) -> float:
    return 100.0 * sum(base in "GC" for base in sequence) / len(sequence)


def circular_segment(
    sequence: str,
    start_1based: int,
    length: int,
) -> str:
    if not sequence or length < 1:
        return ""
    size = len(sequence)
    start = (start_1based - 1) % size
    repeats = math.ceil((start + length) / size) + 1
    expanded = sequence * repeats
    return expanded[start : start + length]


def circular_occurrences(sequence: str, motif: str) -> list[int]:
    motif = motif.upper()
    if not motif or len(motif) > len(sequence):
        return []
    search = sequence + sequence[: len(motif) - 1]
    return [
        match.start() + 1
        for match in re.finditer(f"(?={re.escape(motif)})", search)
        if match.start() < len(sequence)
    ]


def arm_quality(sequence: str, backbone_sequence: str) -> dict[str, Any]:
    gc = gc_percent(sequence)
    homopolymer = max_homopolymer(sequence)
    occurrences = len(circular_occurrences(backbone_sequence, sequence))
    passes = (
        MIN_HOMOLOGY_ARM_GC_PERCENT <= gc <= MAX_HOMOLOGY_ARM_GC_PERCENT
        and homopolymer <= MAX_HOMOLOGY_ARM_HOMOPOLYMER
        and occurrences == 1
    )
    closeness = max(0.0, 1.0 - abs(gc - 50.0) / 20.0)
    return {
        "sequence": sequence,
        "length_bp": len(sequence),
        "gc_percent": round(gc, 4),
        "max_homopolymer": homopolymer,
        "backbone_occurrence_count": occurrences,
        "passes": passes,
        "quality": round(closeness, 4),
    }


def enzyme_catalog() -> tuple[dict[str, Any], ...]:
    enzymes: list[dict[str, Any]] = []
    for name in CURATED_RESTRICTION_ENZYMES:
        enzyme = getattr(Restriction, name, None)
        if enzyme is None:
            continue
        enzymes.append(
            {
                "name": name,
                "enzyme": enzyme,
                "site": str(enzyme.site),
                "site_length": len(str(enzyme.site)),
                "overhang": (
                    "five_prime"
                    if enzyme.is_5overhang()
                    else "three_prime"
                    if enzyme.is_3overhang()
                    else "blunt"
                ),
                "priority": int(ENZYME_PRIORITIES.get(name, 0)),
            }
        )
    return tuple(enzymes)


def enzyme_cut_positions(enzyme: Any, sequence: str, *, circular: bool) -> list[int]:
    return [
        int(position) - 1
        for position in enzyme.search(Seq(sequence), linear=not circular)
    ]


def enzyme_site_start(enzyme: Any, cut_after_bp: int, sequence_length: int) -> int:
    start = cut_after_bp + 1 - int(enzyme.fst5)
    return ((start - 1) % sequence_length) + 1


def region_bounds(region: Mapping[str, Any]) -> tuple[int, int] | None:
    try:
        start = int(region.get("start_bp"))
        end = int(region.get("end_bp"))
    except (TypeError, ValueError):
        return None
    if start < 1 or end < start:
        return None
    return start, end


def span_inside_region(start: int, end: int, region: Mapping[str, Any]) -> bool:
    bounds = region_bounds(region)
    return bool(bounds and bounds[0] <= start <= end <= bounds[1])


def protected_coordinate_spans(
    protected_features: tuple[dict[str, Any], ...],
    genbank_features: tuple[dict[str, Any], ...],
) -> tuple[tuple[int, int], ...]:
    protected_names = {
        str(item.get("name") or "").strip().lower()
        for item in protected_features
        if str(item.get("name") or "").strip()
    }
    protected_roles = {
        str(item.get("role") or "").strip().lower()
        for item in protected_features
        if str(item.get("role") or "").strip()
    }
    spans: list[tuple[int, int]] = []
    for feature in genbank_features:
        haystack = " ".join(
            [str(feature.get("type") or ""), *(feature.get("labels") or [])]
        ).lower()
        feature_type = str(feature.get("type") or "").lower()
        name_match = any(name in haystack for name in protected_names)
        role_match = (
            "replication" in protected_roles
            and feature_type in {"rep_origin", "origin_of_replication"}
        ) or (
            "selection" in protected_roles
            and feature_type in {"cds", "gene"}
            and any(token in haystack for token in ("resistance", "bla", "neo", "cat", "aada", "teta", "aacc1"))
        )
        if name_match or role_match:
            spans.append(
                (int(feature["start_bp"]), int(feature["end_bp"]))
            )
    return tuple(sorted(set(spans)))


def overlaps_spans(start: int, end: int, spans: tuple[tuple[int, int], ...]) -> bool:
    return any(start <= protected_end and end >= protected_start for protected_start, protected_end in spans)


__all__ = [
    "arm_quality",
    "circular_occurrences",
    "circular_segment",
    "enzyme_catalog",
    "enzyme_cut_positions",
    "enzyme_site_start",
    "feature_payloads",
    "gc_percent",
    "max_homopolymer",
    "overlaps_spans",
    "protected_coordinate_spans",
    "read_genbank",
    "region_bounds",
    "relative_project_path",
    "resolve_project_file",
    "sha256_bytes",
    "sha256_sequence",
    "span_inside_region",
    "stable_json_hash",
]
