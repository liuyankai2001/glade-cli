"""Build deterministic main-enzyme combinations for one selected route.

The candidate-selection stage answers "which proteins can catalyse each
step?".  This module answers the next, deliberately separate question:
"which smallest, evidence-supported protein sets cover every required
heterologous step?".

Only rows already selected into the visible Top-N shortlist are consumed.
The full candidate audit is intentionally not searched, so ``--top-n`` keeps
its user-facing meaning and every protein in a combination can be inspected
before it is selected.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from math import isclose
from pathlib import Path
from typing import Any, Mapping, Sequence

from .common import (
    get_solution_steps,
    heterologous_requirements,
    read_csv,
    read_manifest,
    rel_or_abs,
    write_csv,
    write_json_atomic,
)
from .models import (
    MAIN_ENZYME_SETS_ALGORITHM_VERSION,
    MAIN_ENZYME_SETS_SCHEMA_VERSION,
    MainEnzymeSelectionResult,
    MainEnzymeSetsResult,
)
from .provenance import file_sha256, solution_fingerprint, stable_json_hash


MAIN_ENZYME_SETS_FILENAME = "main_enzyme_sets.json"
MAIN_ENZYME_SETS_SUMMARY_FILENAME = "main_enzyme_sets.csv"
MAIN_ENZYME_SET_MEMBERS_FILENAME = "main_enzyme_set_members.csv"

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_FORMING_RE = re.compile(r"\(([^()]*(?:forming|producing)[^()]*)\)", re.I)
_ACCEPTED_FIT = {"verified", "verified_with_risk"}
_DIRECTION_CONTRADICTED = {"contradicted", "unsupported", "rejected"}
_DIRECTION_SUPPORTED = {
    "supported",
    "verified",
    "compatible",
    "forward",
    "reversible",
}
_DIRECT_SPECIFICITY = {"exact", "supported"}

_SUMMARY_COLUMNS = [
    "set_id",
    "status",
    "protein_count",
    "accessions",
    "covered_step_indexes",
    "review_required_step_indexes",
    "organism_count",
    "min_reaction_fit_score",
    "mean_reaction_fit_score",
    "min_protein_score",
    "mean_protein_score",
    "min_host_fit_score",
    "mean_host_fit_score",
    "reviewed_fraction",
    "reaction_fit_risk_count",
    "direction_risk_count",
    "low_direction_confidence_count",
    "specificity_risk_count",
    "exact_specificity_count",
    "carrier_compatibility_status",
    "electron_reassessment_status",
    "set_fingerprint",
    "warnings",
]

_MEMBER_COLUMNS = [
    "set_id",
    "set_status",
    "accession",
    "protein_name",
    "organism_name",
    "reviewed",
    "sequence_sha256",
    "cofactors",
    "capable_step_indexes",
    "assigned_step_indexes",
    "assigned_reaction_ids",
]


def _text(value: Any) -> str:
    return str(value or "").strip()


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(float(_text(value)))
    except (TypeError, ValueError):
        return default


def _float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return _text(value).lower() in {"1", "true", "yes", "y"}


def _split_values(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        values = [_text(item) for item in value]
    else:
        values = [
            item.strip()
            for item in re.split(r"\s*[;|]\s*", _text(value))
        ]
    return sorted({item for item in values if item})


def _score_breakdown(value: Any) -> dict[str, float]:
    if isinstance(value, Mapping):
        raw = value
    else:
        try:
            parsed = json.loads(_text(value))
        except (TypeError, ValueError, json.JSONDecodeError):
            parsed = {}
        raw = parsed if isinstance(parsed, Mapping) else {}
    return {
        key: _float(raw.get(key))
        for key in ("evidence", "expression", "function", "host", "total")
    }


def _normalise_words(value: Any) -> str:
    text = _text(value).lower().replace("α", "alpha").replace("β", "beta")
    return " ".join(re.findall(r"[a-z0-9]+", text))


def _sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def main_enzyme_set_paths(output_dir: str | Path) -> dict[str, Path]:
    """Return output paths without extending candidate evidence paths.

    Keeping these paths separate prevents the candidate-selection artifact
    from claiming not-yet-created combination files as its own evidence.
    """

    root = Path(output_dir)
    return {
        "main_enzyme_sets_json": root / MAIN_ENZYME_SETS_FILENAME,
        "main_enzyme_sets_csv": root / MAIN_ENZYME_SETS_SUMMARY_FILENAME,
        "main_enzyme_set_members_csv": root / MAIN_ENZYME_SET_MEMBERS_FILENAME,
    }


def _canonical_required_steps(
    required_steps: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    canonical: list[dict[str, Any]] = []
    seen: set[int] = set()
    for raw in required_steps:
        status = _text(raw.get("status")).lower()
        if status and status != "heterologous":
            continue
        spontaneous = _bool(raw.get("spontaneous"))
        annotation = " ".join(
            _text(raw.get(key)).lower()
            for key in ("reaction_name", "reaction_comment", "resolution_action")
        )
        if spontaneous or any(
            marker in annotation
            for marker in ("spontaneous", "non-enzymatic", "nonenzymatic")
        ):
            continue
        step_index = _int(raw.get("step_index"))
        if step_index < 1:
            raise ValueError("required step_index must be a positive integer")
        if step_index in seen:
            raise ValueError(f"duplicate required Step {step_index}")
        seen.add(step_index)
        canonical.append(
            {
                "step_index": step_index,
                "reaction_id": _text(raw.get("reaction_id")).upper(),
                "reaction_name": _text(raw.get("reaction_name")),
                "produced_compound_id": _text(
                    raw.get("produced_compound_id")
                ).upper(),
                "produced_compound_name": _text(
                    raw.get("produced_compound_name")
                ),
                "precursor_compound_ids": _split_values(
                    raw.get("precursor_compound_ids")
                ),
                "precursor_compound_labels": _split_values(
                    raw.get("precursor_compound_labels")
                ),
            }
        )
    return sorted(canonical, key=lambda item: item["step_index"])


def _specificity_status(
    row: Mapping[str, Any],
    requirement: Mapping[str, Any],
    *,
    terminal_step: bool,
) -> tuple[str, str]:
    """Classify product specificity conservatively from local annotations."""

    explicit = _text(row.get("specificity_status")).lower()
    if explicit in {"exact", "supported", "unknown", "conflict"}:
        if explicit == "conflict":
            return explicit, "candidate specificity_status is conflict"
        # Still apply known, explicit product conflicts below.

    target_name = _normalise_words(requirement.get("produced_compound_name"))
    target_id = _text(requirement.get("produced_compound_id")).lower()
    evidence = " ".join(
        _text(row.get(key))
        for key in (
            "protein_name",
            "function_comments",
            "catalytic_activities",
            "reaction_fit_evidence",
        )
    )
    normalised = _normalise_words(evidence)

    # A synthase explicitly annotated to continue beyond GGPP is unsuitable as
    # an automatic primary catalyst for a GGPP endpoint, even if its broad EC
    # annotation contains the intermediate reaction.
    if (
        "geranylgeranyl diphosphate" in target_name
        and "nonaprenyl" in normalised
    ):
        return "conflict", "annotation continues GGPP to nonaprenyl diphosphate"

    if terminal_step:
        forming_products = [
            _normalise_words(match)
            for match in _FORMING_RE.findall(evidence)
        ]
        if forming_products and target_name:
            ignored = {"all", "forming", "producing"}
            stereo_tokens = {"cis", "trans"}
            target_all_tokens = set(target_name.split()) - ignored
            target_stereo = target_all_tokens & stereo_tokens
            target_positions = {
                token for token in target_all_tokens if token.isdigit()
            }
            target_core = {
                token
                for token in target_all_tokens - stereo_tokens
                if not token.isdigit()
            }
            product_match_without_stereo = False
            conflicting_product_seen = False
            for product in forming_products:
                product_all_tokens = set(product.split()) - ignored
                product_stereo = product_all_tokens & stereo_tokens
                product_positions = {
                    token for token in product_all_tokens if token.isdigit()
                }
                product_core = {
                    token
                    for token in product_all_tokens - stereo_tokens
                    if not token.isdigit()
                }
                same_core = bool(target_core) and (
                    target_core <= product_core or product_core <= target_core
                )
                positions_match = (
                    not target_positions
                    or target_positions <= product_positions
                )
                if not same_core or not positions_match:
                    conflicting_product_seen = True
                    continue
                if target_stereo and product_stereo:
                    if target_stereo != product_stereo:
                        conflicting_product_seen = True
                        continue
                    return (
                        "exact",
                        "annotated forming product and stereochemistry match route endpoint",
                    )
                if target_stereo and not product_stereo:
                    product_match_without_stereo = True
                    continue
                return "exact", "annotated forming product matches route endpoint"
            if product_match_without_stereo:
                return (
                    "supported",
                    "annotated product matches but stereochemistry is unspecified",
                )
            if conflicting_product_seen:
                return (
                    "conflict",
                    "annotated forming product differs from route endpoint",
                )
            return "conflict", "annotated forming product differs from route endpoint"

    if target_name and target_name in normalised:
        return "exact", "route product is named in candidate annotation"
    if target_id and target_id in evidence.lower():
        return "exact", "route product identifier is named in candidate annotation"
    if explicit in _DIRECT_SPECIFICITY:
        return explicit, "candidate carries explicit specificity evidence"
    if _split_values(row.get("matched_rhea_ids")) or _split_values(
        row.get("matched_ko_ids")
    ):
        return "supported", "candidate has reaction-specific mapping evidence"
    if _text(row.get("reaction_fit_status")) == "verified":
        return "supported", "reaction fit is verified"
    return "unknown", "product specificity is not established"


@dataclass(frozen=True)
class _Edge:
    step_index: int
    reaction_id: str
    accession: str
    candidate_rank: int
    protein_score: float
    host_fit_score: float
    reaction_fit_status: str
    reaction_fit_score: float
    direction_verdict: str
    direction_confidence: str
    specificity_status: str
    specificity_reason: str
    warnings: tuple[str, ...]

    def assignment_order(self) -> tuple[Any, ...]:
        fit_level = 0 if self.reaction_fit_status == "verified" else 1
        direction_level = (
            0 if self.direction_verdict.lower() in _DIRECTION_SUPPORTED else 1
        )
        specificity_level = {
            "exact": 0,
            "supported": 1,
            "unknown": 2,
        }.get(self.specificity_status, 3)
        return (
            specificity_level,
            fit_level,
            direction_level,
            -self.reaction_fit_score,
            -self.host_fit_score,
            -self.protein_score,
            self.candidate_rank,
            self.accession,
        )


@dataclass
class _Protein:
    accession: str
    protein_name: str
    organism_name: str
    reviewed: bool
    sequence_sha256: str
    sequence: str
    cofactors: set[str] = field(default_factory=set)
    edges: dict[int, _Edge] = field(default_factory=dict)
    annotation_texts: set[str] = field(default_factory=set)
    multistep_supported: bool = False
    multistep_evidence_reason: str = ""

    @property
    def capable_steps(self) -> tuple[int, ...]:
        return tuple(sorted(self.edges))


def _names_from_labels(labels: Sequence[str]) -> list[str]:
    names: list[str] = []
    for label in labels:
        match = re.search(r"\((.+)\)\s*$", label)
        value = match.group(1) if match else label
        normalised = _normalise_words(value)
        if normalised and not re.fullmatch(r"c\d{5}", normalised):
            names.append(normalised)
    return names


def _annotation_mentions(annotation: str, compound_name: str) -> bool:
    target = _normalise_words(compound_name)
    return bool(target and target in _normalise_words(annotation))


def _multistep_annotation_evidence(
    protein: _Protein,
    requirements: Mapping[int, Mapping[str, Any]],
) -> tuple[bool, str]:
    steps = list(protein.capable_steps)
    if len(steps) < 2:
        return True, "single-step candidate"

    first = requirements[steps[0]]
    last = requirements[steps[-1]]
    start_names = _names_from_labels(first.get("precursor_compound_labels", []))
    end_name = _normalise_words(last.get("produced_compound_name"))
    intermediate_names = [
        _normalise_words(requirements[step].get("produced_compound_name"))
        for step in steps[:-1]
        if _normalise_words(
            requirements[step].get("produced_compound_name")
        )
    ]

    # Strongest evidence: one recorded catalytic activity or function statement
    # names both ends of the collapsed interval.  Exact normalised phrases keep
    # positional/stereochemical labels (for example 15-cis vs 9-cis).
    for annotation in protein.annotation_texts:
        if (
            end_name
            and _annotation_mentions(annotation, end_name)
            and any(
                _annotation_mentions(annotation, start_name)
                for start_name in start_names
            )
        ):
            return True, "one annotation names the whole-reaction start and end"

    corpus = " ".join(sorted(protein.annotation_texts))
    chain_named = (
        bool(start_names)
        and any(
            _annotation_mentions(corpus, start_name)
            for start_name in start_names
        )
        and bool(end_name)
        and _annotation_mentions(corpus, end_name)
        and all(
            _annotation_mentions(corpus, name)
            for name in intermediate_names
        )
    )
    if chain_named:
        return True, "annotations name the route start, end, and intermediates"

    # This branch is primarily useful for curated/offline fixtures and future
    # sources that cite KEGG reactions directly.  It is stronger than merely
    # sharing an EC because every Step and every Step product is named.
    reactions_named = all(
        _text(requirements[step].get("reaction_id")).lower()
        in corpus.lower()
        for step in steps
    )
    products_named = all(
        _annotation_mentions(
            corpus, _text(requirements[step].get("produced_compound_name"))
        )
        for step in steps
    )
    if reactions_named and products_named:
        return True, "annotations explicitly name every compressed Step"
    return False, "no auditable whole-reaction or complete chain annotation"


def _canonical_candidate_projection(row: Mapping[str, Any]) -> dict[str, Any]:
    """Project only fields that affect filtering, scoring, or audit output."""

    return {
        "step_index": _int(row.get("step_index")),
        "reaction_id": _text(row.get("reaction_id")).upper(),
        "candidate_rank": _int(row.get("candidate_rank")),
        "selection_status": _text(row.get("selection_status")).lower(),
        "candidate_role": _text(
            row.get("candidate_role") or row.get("role")
        ).lower(),
        "accession": _text(row.get("accession")).upper(),
        "protein_name": _text(row.get("protein_name")),
        "organism_name": _text(row.get("organism_name")),
        "reviewed": _bool(row.get("reviewed")),
        "score": _float(row.get("score") or row.get("protein_score")),
        "score_breakdown": _score_breakdown(row.get("score_breakdown")),
        "reaction_fit_status": _text(row.get("reaction_fit_status")).lower(),
        "reaction_fit_score": _float(row.get("reaction_fit_score")),
        "reaction_fit_evidence": _text(row.get("reaction_fit_evidence")),
        "direction_verdict": _text(row.get("direction_verdict")).lower(),
        "direction_confidence": _text(row.get("direction_confidence")).lower(),
        "specificity_status": _text(row.get("specificity_status")).lower(),
        "matched_rhea_ids": _split_values(row.get("matched_rhea_ids")),
        "matched_ko_ids": _split_values(row.get("matched_ko_ids")),
        "catalytic_activities": _text(row.get("catalytic_activities")),
        "function_comments": _text(row.get("function_comments")),
        "cofactors": _split_values(row.get("cofactors")),
        "warnings": _split_values(row.get("warnings")),
        "sequence_sha256": _text(row.get("sequence_sha256")).lower(),
        "sequence": _text(row.get("sequence")),
    }


def shortlist_decision_fingerprint(
    candidate_rows: Sequence[Mapping[str, Any]],
) -> str:
    """Fingerprint all shortlisted fields that can affect combinations."""

    projected = sorted(
        (
            _canonical_candidate_projection(row)
            for row in candidate_rows
            if _text(row.get("selection_status")).lower() == "selected"
        ),
        key=lambda item: (
            item["step_index"],
            item["candidate_rank"],
            item["accession"],
        ),
    )
    return stable_json_hash(projected)


def _build_protein_pool(
    required_steps: list[dict[str, Any]],
    candidate_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[_Protein], list[str]]:
    requirements = {item["step_index"]: item for item in required_steps}
    terminal_step = max(requirements, default=0)
    proteins: dict[str, _Protein] = {}
    rejected_reasons: list[str] = []
    seen_keys: set[tuple[int, int, str]] = set()
    seen_ranks: set[tuple[int, int]] = set()
    seen_step_accessions: set[tuple[int, str]] = set()

    projected_rows = sorted(
        (_canonical_candidate_projection(row) for row in candidate_rows),
        key=lambda row: (
            row["step_index"],
            row["candidate_rank"],
            row["accession"],
        ),
    )
    for row in projected_rows:
        # The visible Top-N shortlist is the complete search boundary.
        if row["selection_status"] != "selected":
            continue
        role = row["candidate_role"]
        if role and role not in {"main", "catalytic_main"}:
            continue
        step_index = row["step_index"]
        if step_index not in requirements:
            continue
        accession = row["accession"]
        rank = row["candidate_rank"]
        key = (step_index, rank, accession)
        if key in seen_keys:
            raise ValueError(
                f"duplicate shortlist candidate Step {step_index} rank {rank} "
                f"accession {accession}"
            )
        seen_keys.add(key)
        rank_key = (step_index, rank)
        if rank_key in seen_ranks:
            raise ValueError(
                f"duplicate shortlist rank {rank} within Step {step_index}"
            )
        seen_ranks.add(rank_key)
        step_accession_key = (step_index, accession)
        if step_accession_key in seen_step_accessions:
            raise ValueError(
                f"duplicate accession {accession} within Step {step_index}"
            )
        seen_step_accessions.add(step_accession_key)
        fit_status = row["reaction_fit_status"]
        direction = row["direction_verdict"]
        sequence = row["sequence"]
        sequence_sha256 = row["sequence_sha256"]

        invalid: list[str] = []
        if not accession:
            invalid.append("missing accession")
        if rank < 1:
            invalid.append("invalid candidate rank")
        if fit_status not in _ACCEPTED_FIT:
            invalid.append(f"reaction fit {fit_status or 'missing'}")
        if direction in _DIRECTION_CONTRADICTED:
            invalid.append(f"direction {direction}")
        if not sequence:
            invalid.append("missing sequence")
        if not _SHA256_RE.fullmatch(sequence_sha256):
            invalid.append("invalid sequence SHA256")
        elif sequence:
            try:
                actual_digest = _sequence_digest(sequence)
            except UnicodeEncodeError:
                invalid.append("sequence is not ASCII")
            else:
                if actual_digest != sequence_sha256:
                    invalid.append("sequence SHA256 mismatch")
        if row["reaction_id"] != requirements[step_index]["reaction_id"]:
            invalid.append("reaction ID mismatch")
        if invalid:
            rejected_reasons.append(
                f"Step {step_index} {accession or '<missing>'}: "
                + ", ".join(invalid)
            )
            continue

        specificity, specificity_reason = _specificity_status(
            row,
            requirements[step_index],
            terminal_step=step_index == terminal_step,
        )
        if specificity == "conflict":
            rejected_reasons.append(
                f"Step {step_index} {accession}: {specificity_reason}"
            )
            continue

        edge = _Edge(
            step_index=step_index,
            reaction_id=row["reaction_id"],
            accession=accession,
            candidate_rank=rank,
            protein_score=row["score"],
            host_fit_score=row["score_breakdown"]["host"],
            reaction_fit_status=fit_status,
            reaction_fit_score=row["reaction_fit_score"],
            direction_verdict=direction,
            direction_confidence=row["direction_confidence"],
            specificity_status=specificity,
            specificity_reason=specificity_reason,
            warnings=tuple(row["warnings"]),
        )

        if accession not in proteins:
            proteins[accession] = _Protein(
                accession=accession,
                protein_name=row["protein_name"],
                organism_name=row["organism_name"],
                reviewed=row["reviewed"],
                sequence_sha256=sequence_sha256,
                sequence=sequence,
            )
        protein = proteins[accession]
        if protein.sequence_sha256 != sequence_sha256 or protein.sequence != sequence:
            raise ValueError(
                f"accession {accession} has inconsistent sequences across steps"
            )
        if protein.protein_name != row["protein_name"]:
            raise ValueError(
                f"accession {accession} has inconsistent protein names"
            )
        if protein.organism_name != row["organism_name"]:
            raise ValueError(
                f"accession {accession} has inconsistent organism names"
            )
        protein.reviewed = protein.reviewed and row["reviewed"]
        protein.cofactors.update(row["cofactors"])
        protein.annotation_texts.update(
            value
            for value in (
                row["protein_name"],
                row["catalytic_activities"],
                row["function_comments"],
            )
            if value
        )
        existing = protein.edges.get(step_index)
        if existing is None or edge.assignment_order() < existing.assignment_order():
            protein.edges[step_index] = edge

    for protein in proteins.values():
        (
            protein.multistep_supported,
            protein.multistep_evidence_reason,
        ) = _multistep_annotation_evidence(protein, requirements)
    return sorted(proteins.values(), key=lambda item: item.accession), rejected_reasons


def _edge_quality(edge: _Edge) -> tuple[Any, ...]:
    return edge.assignment_order()


def _protein_branch_order(protein: _Protein) -> tuple[Any, ...]:
    edges = list(protein.edges.values())
    return (
        -len(edges),
        sum(1 for edge in edges if edge.specificity_status != "exact"),
        sum(1 for edge in edges if edge.reaction_fit_status != "verified"),
        sum(
            1
            for edge in edges
            if edge.direction_verdict not in _DIRECTION_SUPPORTED
        ),
        -min((edge.reaction_fit_score for edge in edges), default=0.0),
        -min((edge.protein_score for edge in edges), default=0.0),
        protein.accession,
    )


@dataclass
class _SearchResult:
    assignments: list[tuple[tuple[int, str], ...]]
    minimum_count: int | None
    nodes: int
    complete: bool


def _assignment_quality(
    assignment: tuple[tuple[int, str], ...],
    protein_by_accession: Mapping[str, _Protein],
) -> tuple[Any, ...]:
    edges = [
        protein_by_accession[accession].edges[step_index]
        for step_index, accession in assignment
    ]
    return (
        sum(edge.specificity_status not in _DIRECT_SPECIFICITY for edge in edges),
        -sum(edge.specificity_status == "exact" for edge in edges),
        sum(edge.reaction_fit_status != "verified" for edge in edges),
        sum(edge.direction_verdict not in _DIRECTION_SUPPORTED for edge in edges),
        sum(
            edge.direction_confidence not in {"high", "medium"}
            for edge in edges
        ),
        -min((edge.reaction_fit_score for edge in edges), default=0.0),
        -sum(edge.reaction_fit_score for edge in edges),
        -min((edge.host_fit_score for edge in edges), default=0.0),
        -sum(edge.host_fit_score for edge in edges),
        -min((edge.protein_score for edge in edges), default=0.0),
        -sum(edge.protein_score for edge in edges),
        assignment,
    )


def _enumerate_valid_assignments(
    required_step_indexes: list[int],
    proteins: list[_Protein],
    *,
    max_search_nodes: int,
) -> _SearchResult:
    """Enumerate primary Step assignments under multi-step evidence gates."""

    if not required_step_indexes:
        return _SearchResult([tuple()], 0, 0, True)

    protein_by_accession = {item.accession: item for item in proteins}
    by_step = {
        step_index: sorted(
            [item.accession for item in proteins if step_index in item.edges],
            key=lambda accession: (
                _edge_quality(
                    protein_by_accession[accession].edges[step_index]
                ),
                _protein_branch_order(protein_by_accession[accession]),
            ),
        )
        for step_index in required_step_indexes
    }
    if any(not by_step[step] for step in required_step_indexes):
        return _SearchResult([], None, 0, True)

    step_order = sorted(
        required_step_indexes,
        key=lambda step_index: (len(by_step[step_index]), step_index),
    )
    # One best primary assignment is retained for each distinct protein set.
    # This includes higher-cardinality sets when they reduce biochemical risk.
    found: dict[tuple[str, ...], tuple[tuple[int, str], ...]] = {}
    nodes = 0
    truncated = False

    def visit(
        position: int,
        assignment_by_step: dict[int, str],
        assigned_by_accession: dict[str, list[int]],
    ) -> None:
        nonlocal nodes, truncated
        if truncated:
            return
        if nodes >= max_search_nodes:
            truncated = True
            return
        nodes += 1

        if position == len(step_order):
            assignment = tuple(sorted(assignment_by_step.items()))
            accession_set = tuple(sorted(assigned_by_accession))
            current = found.get(accession_set)
            if current is None or _assignment_quality(
                assignment, protein_by_accession
            ) < _assignment_quality(current, protein_by_accession):
                found[accession_set] = assignment
            return

        step_index = step_order[position]
        candidates = list(by_step[step_index])
        candidates.sort(
            key=lambda accession: (
                0 if accession in assigned_by_accession else 1,
                _edge_quality(
                    protein_by_accession[accession].edges[step_index]
                ),
                _protein_branch_order(protein_by_accession[accession]),
            )
        )
        for accession in candidates:
            protein = protein_by_accession[accession]
            previous_steps = assigned_by_accession.get(accession, [])
            if previous_steps and not protein.multistep_supported:
                continue
            assignment_by_step[step_index] = accession
            if accession not in assigned_by_accession:
                assigned_by_accession[accession] = []
            assigned_by_accession[accession].append(step_index)
            visit(position + 1, assignment_by_step, assigned_by_accession)
            assigned_by_accession[accession].pop()
            if not assigned_by_accession[accession]:
                del assigned_by_accession[accession]
            del assignment_by_step[step_index]
            if truncated:
                return

    visit(0, {}, {})
    assignments = sorted(
        found.values(),
        key=lambda item: (
            tuple(sorted(accession for _, accession in item)),
            item,
        ),
    )
    minimum_count = min(
        (len({accession for _, accession in item}) for item in assignments),
        default=None,
    )
    return _SearchResult(
        assignments=assignments,
        minimum_count=minimum_count,
        nodes=nodes,
        complete=not truncated,
    )


def _electron_risk_steps(electron_inference: Mapping[str, Any]) -> list[int]:
    direct = electron_inference.get("electron_risk_step_indexes")
    indexes = {_int(value) for value in _split_values(direct)}
    risk_steps = electron_inference.get("risk_steps")
    if isinstance(risk_steps, list):
        indexes.update(
            _int(item.get("step_index"))
            for item in risk_steps
            if isinstance(item, Mapping)
        )
    step_requirements = electron_inference.get("step_requirements")
    if isinstance(step_requirements, list):
        indexes.update(
            _int(item.get("step_index"))
            for item in step_requirements
            if isinstance(item, Mapping)
        )
    return sorted(index for index in indexes if index > 0)


def _electron_flags(electron_inference: Mapping[str, Any]) -> tuple[bool, bool]:
    summary = electron_inference.get("summary")
    merged = dict(summary) if isinstance(summary, Mapping) else {}
    merged.update(electron_inference)
    return (
        _bool(merged.get("requires_external_electron_regeneration")),
        _bool(merged.get("requires_carrier_compatibility_check")),
    )


def _build_set_payload(
    primary_assignment: tuple[tuple[int, str], ...],
    protein_by_accession: Mapping[str, _Protein],
    requirements: Mapping[int, Mapping[str, Any]],
    electron_inference: Mapping[str, Any],
    candidate_pool_fingerprint: str,
) -> dict[str, Any]:
    accessions = tuple(sorted({item[1] for item in primary_assignment}))
    selected = [protein_by_accession[accession] for accession in accessions]
    assignments: list[dict[str, Any]] = []
    assigned_by_accession: dict[str, list[int]] = defaultdict(list)
    review_steps: set[int] = set()
    set_warnings: set[str] = set()

    assignment_lookup = dict(primary_assignment)
    for step_index in sorted(requirements):
        accession = assignment_lookup[step_index]
        edge = protein_by_accession[accession].edges[step_index]
        assigned_by_accession[edge.accession].append(step_index)
        assignments.append(
            {
                "step_index": step_index,
                "reaction_id": edge.reaction_id,
                "accession": edge.accession,
                "candidate_rank": edge.candidate_rank,
                "protein_score": edge.protein_score,
                "host_fit_score": edge.host_fit_score,
                "reaction_fit_status": edge.reaction_fit_status,
                "reaction_fit_score": edge.reaction_fit_score,
                "direction_verdict": edge.direction_verdict,
                "direction_confidence": edge.direction_confidence,
                "specificity_status": edge.specificity_status,
            }
        )
        if edge.reaction_fit_status != "verified":
            review_steps.add(step_index)
        if edge.direction_verdict not in _DIRECTION_SUPPORTED:
            review_steps.add(step_index)
        if edge.direction_confidence not in {"high", "medium"}:
            review_steps.add(step_index)
        if edge.specificity_status not in _DIRECT_SPECIFICITY:
            review_steps.add(step_index)
        for warning in edge.warnings:
            set_warnings.add(f"Step {step_index} {edge.accession}: {warning}")
            review_steps.add(step_index)

    risk_steps = [
        step for step in _electron_risk_steps(electron_inference) if step in requirements
    ]
    external_required, carrier_check = _electron_flags(electron_inference)
    if not risk_steps and not external_required and not carrier_check:
        carrier_status = "not_required"
        electron_status = "not_required"
        electron_assessment = "No route-level electron-system reassessment is required."
    else:
        risk_accessions = {
            assignment["accession"]
            for assignment in assignments
            if assignment["step_index"] in risk_steps
        }
        if external_required:
            carrier_status = "external_regeneration_required"
            electron_assessment = (
                "The selected set still requires an external electron "
                "regeneration mechanism."
            )
        elif len(risk_steps) > 1 and len(risk_accessions) == 1:
            carrier_status = "single_multistep_enzyme"
            electron_assessment = (
                "The selected multi-step enzyme replaces decomposed carrier "
                "steps; reassess its whole-reaction electron acceptor and "
                "cofactor reoxidation before downstream design."
            )
        else:
            carrier_status = "compatibility_review_required"
            electron_assessment = (
                "Electron-carrier compatibility remains unresolved across "
                "the selected proteins."
            )
        electron_status = "review_required"
        review_steps.update(risk_steps)
        set_warnings.add(electron_assessment)

    proteins_payload: list[dict[str, Any]] = []
    for protein in sorted(selected, key=lambda item: item.accession):
        proteins_payload.append(
            {
                "accession": protein.accession,
                "protein_name": protein.protein_name,
                "organism_name": protein.organism_name,
                "reviewed": protein.reviewed,
                "sequence_sha256": protein.sequence_sha256,
                "cofactors": sorted(protein.cofactors),
                "capable_step_indexes": list(protein.capable_steps),
                "assigned_step_indexes": sorted(
                    assigned_by_accession[protein.accession]
                ),
            }
        )

    fit_scores = [item["reaction_fit_score"] for item in assignments]
    protein_scores = [item["protein_score"] for item in assignments]
    host_fit_scores = [item["host_fit_score"] for item in assignments]
    warnings = sorted(set_warnings)
    review_indexes = sorted(review_steps)
    status = "review_required" if review_indexes or warnings else "complete"
    metrics = {
        "protein_count": len(proteins_payload),
        "organism_count": len(
            {
                item["organism_name"]
                for item in proteins_payload
                if item["organism_name"]
            }
        ),
        "min_reaction_fit_score": min(fit_scores),
        "mean_reaction_fit_score": sum(fit_scores) / len(fit_scores),
        "min_protein_score": min(protein_scores),
        "mean_protein_score": sum(protein_scores) / len(protein_scores),
        "min_host_fit_score": min(host_fit_scores),
        "mean_host_fit_score": sum(host_fit_scores) / len(host_fit_scores),
        "reviewed_fraction": sum(
            bool(item["reviewed"]) for item in proteins_payload
        )
        / len(proteins_payload),
        "reaction_fit_risk_count": sum(
            item["reaction_fit_status"] != "verified"
            for item in assignments
        ),
        "direction_risk_count": sum(
            item["direction_verdict"] not in _DIRECTION_SUPPORTED
            for item in assignments
        ),
        "low_direction_confidence_count": sum(
            item["direction_confidence"] not in {"high", "medium"}
            for item in assignments
        ),
        "specificity_risk_count": sum(
            item["specificity_status"] not in _DIRECT_SPECIFICITY
            for item in assignments
        ),
        "exact_specificity_count": sum(
            item["specificity_status"] == "exact"
            for item in assignments
        ),
        "warning_count": len(warnings),
        "carrier_compatibility_status": carrier_status,
        "electron_reassessment_status": electron_status,
    }
    set_fingerprint = stable_json_hash(
        {
            "candidate_pool_fingerprint": candidate_pool_fingerprint,
            "accessions": sorted(accessions),
            "assignments": [
                {
                    "step_index": item["step_index"],
                    "accession": item["accession"],
                    "candidate_rank": item["candidate_rank"],
                }
                for item in assignments
            ],
        }
    )
    reasons = [
        "covers every required heterologous step",
        f"uses {len(proteins_payload)} distinct protein(s)",
    ]
    reasons.extend(
        f"{protein.accession}: {protein.multistep_evidence_reason}"
        for protein in selected
        if len(assigned_by_accession[protein.accession]) > 1
    )
    return {
        "set_id": 1,  # Reassigned after deterministic ranking.
        "set_fingerprint": set_fingerprint,
        "status": status,
        "protein_count": len(proteins_payload),
        "coverage_complete": True,
        "covered_step_indexes": sorted(requirements),
        "uncovered_step_indexes": [],
        "proteins": proteins_payload,
        "step_assignments": assignments,
        "review_required_step_indexes": review_indexes,
        "electron_assessment": electron_assessment,
        "metrics": metrics,
        "reasons": reasons,
        "warnings": warnings,
    }


def _set_ranking_key(payload: Mapping[str, Any]) -> tuple[Any, ...]:
    metrics = payload["metrics"]
    accessions = tuple(item["accession"] for item in payload["proteins"])
    electron_severity = {
        "not_required": 0,
        "single_multistep_enzyme": 1,
        "compatibility_review_required": 2,
        "external_regeneration_required": 3,
    }.get(metrics["carrier_compatibility_status"], 4)
    return (
        0 if payload["status"] == "complete" else 1,
        metrics["specificity_risk_count"],
        -metrics["exact_specificity_count"],
        metrics["reaction_fit_risk_count"],
        metrics["direction_risk_count"],
        metrics["low_direction_confidence_count"],
        electron_severity,
        payload["protein_count"],
        -metrics["min_reaction_fit_score"],
        -metrics["mean_reaction_fit_score"],
        len(payload["warnings"]),
        -metrics["reviewed_fraction"],
        metrics["organism_count"],
        -metrics["min_host_fit_score"],
        -metrics["mean_host_fit_score"],
        -metrics["min_protein_score"],
        -metrics["mean_protein_score"],
        accessions,
    )


def _candidate_pool_fingerprint(
    *,
    solution_fingerprint_value: str,
    chassis_key: str,
    required_steps: list[dict[str, Any]],
    candidate_rows: Sequence[Mapping[str, Any]],
    electron_inference: Mapping[str, Any],
) -> str:
    candidates = sorted(
        (
            _canonical_candidate_projection(row)
            for row in candidate_rows
            if _text(row.get("selection_status")).lower() == "selected"
        ),
        key=lambda item: (
            item["step_index"],
            item["candidate_rank"],
            item["accession"],
        ),
    )
    return stable_json_hash(
        {
            "algorithm_version": MAIN_ENZYME_SETS_ALGORITHM_VERSION,
            "solution_fingerprint": solution_fingerprint_value,
            "chassis_key": chassis_key,
            "required_steps": required_steps,
            "candidate_scope": "top_n_shortlist",
            "candidates": candidates,
            "electron_inference": electron_inference,
        }
    )


def candidate_pool_fingerprint_from_rows(
    *,
    solution_fingerprint_value: str,
    chassis_key: str,
    required_steps: Sequence[Mapping[str, Any]],
    candidate_rows: Sequence[Mapping[str, Any]],
    electron_inference: Mapping[str, Any] | None,
) -> str:
    """Return the path-independent semantic fingerprint used by set results."""

    return _candidate_pool_fingerprint(
        solution_fingerprint_value=solution_fingerprint_value,
        chassis_key=chassis_key,
        required_steps=_canonical_required_steps(required_steps),
        candidate_rows=candidate_rows,
        electron_inference=dict(electron_inference or {}),
    )


def _input_fingerprint(
    candidate_pool_fingerprint: str,
    *,
    max_sets: int,
    max_search_nodes: int,
) -> str:
    return stable_json_hash(
        {
            "candidate_pool_fingerprint": candidate_pool_fingerprint,
            "parameters": {
                "max_sets": max_sets,
                "max_search_nodes": max_search_nodes,
            },
        }
    )


def build_main_enzyme_sets_from_rows(
    *,
    solution_id: int,
    expansion_depth: int,
    solution_fingerprint_value: str,
    chassis_key: str,
    required_steps: Sequence[Mapping[str, Any]],
    candidate_rows: Sequence[Mapping[str, Any]],
    electron_inference: Mapping[str, Any] | None,
    max_sets: int = 20,
    max_search_nodes: int = 1_000_000,
    source_artifacts: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Pure, offline enzyme-set search over validated shortlist-like rows."""

    if solution_id < 1:
        raise ValueError("solution_id must be positive")
    if expansion_depth < 0:
        raise ValueError("expansion_depth must not be negative")
    if not _SHA256_RE.fullmatch(_text(solution_fingerprint_value)):
        raise ValueError("solution_fingerprint_value must be a SHA256 digest")
    if max_sets < 1:
        raise ValueError("max_sets must be positive")
    if max_search_nodes < 1:
        raise ValueError("max_search_nodes must be positive")

    requirements_list = _canonical_required_steps(required_steps)
    requirement_map = {
        item["step_index"]: item for item in requirements_list
    }
    required_indexes = sorted(requirement_map)
    electron_context = dict(electron_inference or {})
    candidate_pool_fingerprint = _candidate_pool_fingerprint(
        solution_fingerprint_value=solution_fingerprint_value,
        chassis_key=chassis_key,
        required_steps=requirements_list,
        candidate_rows=candidate_rows,
        electron_inference=electron_context,
    )
    input_fingerprint = _input_fingerprint(
        candidate_pool_fingerprint,
        max_sets=max_sets,
        max_search_nodes=max_search_nodes,
    )
    parameters = {
        "max_sets": max_sets,
        "max_search_nodes": max_search_nodes,
        "candidate_scope": "top_n_shortlist",
    }

    if not required_indexes:
        result = MainEnzymeSetsResult(
            schema_version=MAIN_ENZYME_SETS_SCHEMA_VERSION,
            algorithm_version=MAIN_ENZYME_SETS_ALGORITHM_VERSION,
            ok=True,
            status="complete",
            selected_solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint=solution_fingerprint_value,
            chassis_key=chassis_key,
            parameters=parameters,
            candidate_pool_fingerprint=candidate_pool_fingerprint,
            input_fingerprint=input_fingerprint,
            required_step_indexes=[],
            minimum_protein_count=None,
            search_complete=True,
            search_nodes=0,
            sets=[],
            uncovered_step_indexes=[],
            blocking_reasons=[],
            warnings=["The selected route has no enzyme-requiring heterologous steps."],
            source_artifacts=dict(source_artifacts or {}),
        )
        return result.model_dump(mode="json")

    proteins, rejected_reasons = _build_protein_pool(
        requirements_list, candidate_rows
    )
    covered_by_pool = {
        step for protein in proteins for step in protein.capable_steps
    }
    missing = sorted(set(required_indexes) - covered_by_pool)
    if missing:
        blocking = [
            "No constructable, non-conflicting Top-N candidate for Step "
            + ", ".join(map(str, missing))
        ]
        if rejected_reasons:
            blocking.extend(sorted(rejected_reasons))
        result = MainEnzymeSetsResult(
            schema_version=MAIN_ENZYME_SETS_SCHEMA_VERSION,
            algorithm_version=MAIN_ENZYME_SETS_ALGORITHM_VERSION,
            ok=False,
            status="infeasible",
            selected_solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint=solution_fingerprint_value,
            chassis_key=chassis_key,
            parameters=parameters,
            candidate_pool_fingerprint=candidate_pool_fingerprint,
            input_fingerprint=input_fingerprint,
            required_step_indexes=required_indexes,
            minimum_protein_count=None,
            search_complete=True,
            search_nodes=0,
            sets=[],
            uncovered_step_indexes=required_indexes,
            blocking_reasons=blocking,
            warnings=[],
            source_artifacts=dict(source_artifacts or {}),
        )
        return result.model_dump(mode="json")

    search = _enumerate_valid_assignments(
        required_indexes,
        proteins,
        max_search_nodes=max_search_nodes,
    )
    protein_by_accession = {item.accession: item for item in proteins}
    set_payloads: list[dict[str, Any]] = []
    for assignment in search.assignments:
        payload = _build_set_payload(
            assignment,
            protein_by_accession,
            requirement_map,
            electron_context,
            candidate_pool_fingerprint,
        )
        if len(set_payloads) < max_sets:
            set_payloads.append(payload)
            set_payloads.sort(key=_set_ranking_key)
        elif _set_ranking_key(payload) < _set_ranking_key(set_payloads[-1]):
            set_payloads[-1] = payload
            set_payloads.sort(key=_set_ranking_key)
    for set_id, payload in enumerate(set_payloads, start=1):
        payload["set_id"] = set_id

    if not search.complete:
        status = "truncated"
        ok = bool(set_payloads)
        warnings = [
            f"Search stopped at max_search_nodes={max_search_nodes}; "
            "reported sets may be incomplete."
        ]
    elif not set_payloads:
        status = "infeasible"
        ok = False
        warnings = []
    elif set_payloads[0]["status"] == "review_required":
        status = "review_required"
        ok = True
        warnings = []
    else:
        status = "complete"
        ok = True
        warnings = (
            ["Additional lower-ranked combinations require review."]
            if any(
                item["status"] == "review_required"
                for item in set_payloads[1:]
            )
            else []
        )

    blocking_reasons: list[str] = []
    if status == "infeasible":
        blocking_reasons.append(
            "No complete combination was found in the selected Top-N shortlist."
        )
    result = MainEnzymeSetsResult(
        schema_version=MAIN_ENZYME_SETS_SCHEMA_VERSION,
        algorithm_version=MAIN_ENZYME_SETS_ALGORITHM_VERSION,
        ok=ok,
        status=status,
        selected_solution_id=solution_id,
        expansion_depth=expansion_depth,
        solution_fingerprint=solution_fingerprint_value,
        chassis_key=chassis_key,
        parameters=parameters,
        candidate_pool_fingerprint=candidate_pool_fingerprint,
        input_fingerprint=input_fingerprint,
        required_step_indexes=required_indexes,
        minimum_protein_count=search.minimum_count if set_payloads else None,
        search_complete=search.complete,
        search_nodes=search.nodes,
        sets=set_payloads,
        uncovered_step_indexes=[] if set_payloads else required_indexes,
        blocking_reasons=blocking_reasons,
        warnings=warnings,
        source_artifacts=dict(source_artifacts or {}),
    )
    return result.model_dump(mode="json")


def _selection_projection(selection: MainEnzymeSelectionResult) -> dict[str, Any]:
    return {
        "schema_version": selection.schema_version,
        "selected_solution_id": selection.selected_solution_id,
        "expansion_depth": selection.expansion_depth,
        "solution_fingerprint": selection.solution_fingerprint,
        "chassis_key": selection.chassis_key,
        "parameters": selection.parameters.model_dump(mode="json"),
        "shortlist_decision_fingerprint": (
            selection.shortlist_decision_fingerprint
        ),
        "candidates_by_step": {
            str(step): [
                candidate.model_dump(mode="json")
                for candidate in sorted(
                    candidates, key=lambda item: item.candidate_rank
                )
            ]
            for step, candidates in sorted(selection.candidates_by_step.items())
        },
        "uncovered_step_indexes": selection.uncovered_step_indexes,
        "direction_rejected_step_indexes": (
            selection.direction_rejected_step_indexes
        ),
        "direction_risk_step_indexes": selection.direction_risk_step_indexes,
    }


def main_enzyme_selection_fingerprint(
    selection: MainEnzymeSelectionResult,
) -> str:
    """Fingerprint a validated candidate selection without absolute paths."""

    return stable_json_hash(_selection_projection(selection))


def _validate_selection_csv(
    selection: MainEnzymeSelectionResult,
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    selected_rows = [
        row
        for row in rows
        if _text(row.get("selection_status")).lower() == "selected"
    ]
    by_key: dict[tuple[int, int, str], dict[str, str]] = {}
    for row in selected_rows:
        if _int(row.get("solution_id")) != selection.selected_solution_id:
            raise ValueError(
                "candidate CSV solution_id does not match canonical selection"
            )
        role = _text(row.get("candidate_role") or row.get("role")).lower()
        if role not in {"main", "catalytic_main"}:
            raise ValueError(
                "candidate CSV contains a selected non-main protein row"
            )
        key = (
            _int(row.get("step_index")),
            _int(row.get("candidate_rank")),
            _text(row.get("accession")).upper(),
        )
        if key in by_key:
            raise ValueError(f"duplicate candidate CSV row: {key}")
        by_key[key] = row

    expected_keys: set[tuple[int, int, str]] = set()
    for step_index, candidates in selection.candidates_by_step.items():
        for candidate in candidates:
            key = (step_index, candidate.candidate_rank, candidate.accession)
            expected_keys.add(key)
            row = by_key.get(key)
            if row is None:
                raise ValueError(f"candidate CSV is missing {key}")
            if _text(row.get("reaction_id")).upper() != candidate.reaction_id:
                raise ValueError(f"candidate CSV reaction mismatch for {key}")
            scalar_checks = (
                (
                    _text(row.get("ec_number")),
                    candidate.ec_number or "",
                    "EC number",
                ),
                (
                    _text(row.get("protein_name")),
                    candidate.protein_name,
                    "protein name",
                ),
                (
                    _text(row.get("organism_name")),
                    candidate.organism_name,
                    "organism name",
                ),
                (
                    _text(row.get("reaction_fit_status")).lower(),
                    candidate.reaction_fit_status,
                    "reaction fit status",
                ),
                (
                    _text(row.get("direction_verdict")).lower(),
                    candidate.direction_verdict.lower(),
                    "direction verdict",
                ),
                (
                    _text(row.get("direction_confidence")).lower(),
                    candidate.direction_confidence.lower(),
                    "direction confidence",
                ),
                (
                    _text(row.get("reaction_confidence")),
                    candidate.reaction_confidence,
                    "reaction confidence",
                ),
            )
            for actual, expected, label in scalar_checks:
                if actual != expected:
                    raise ValueError(
                        f"candidate CSV {label} mismatch for {key}"
                    )
            if _bool(row.get("reviewed")) != candidate.reviewed:
                raise ValueError(f"candidate CSV reviewed mismatch for {key}")
            if candidate.length is None:
                if _text(row.get("length")):
                    raise ValueError(f"candidate CSV length mismatch for {key}")
            elif _int(row.get("length")) != candidate.length:
                raise ValueError(f"candidate CSV length mismatch for {key}")
            if not isclose(
                _float(row.get("score")),
                candidate.protein_score,
                rel_tol=1e-9,
                abs_tol=1e-9,
            ):
                raise ValueError(f"candidate CSV score mismatch for {key}")
            if not isclose(
                _float(row.get("reaction_fit_score")),
                candidate.reaction_fit_score,
                rel_tol=1e-9,
                abs_tol=1e-9,
            ):
                raise ValueError(
                    f"candidate CSV reaction fit score mismatch for {key}"
                )
            list_checks = (
                (
                    _split_values(row.get("retrieval_strategy")),
                    sorted(candidate.retrieval_strategies),
                    "retrieval strategies",
                ),
                (
                    _split_values(row.get("matched_rhea_ids")),
                    sorted(candidate.matched_rhea_ids),
                    "matched Rhea IDs",
                ),
                (
                    _split_values(row.get("matched_ko_ids")),
                    sorted(candidate.matched_ko_ids),
                    "matched KO IDs",
                ),
                (
                    _split_values(row.get("warnings")),
                    sorted(candidate.warnings),
                    "warnings",
                ),
                (
                    _split_values(row.get("reasons")),
                    sorted(candidate.reasons),
                    "reasons",
                ),
            )
            for actual, expected, label in list_checks:
                if actual != expected:
                    raise ValueError(
                        f"candidate CSV {label} mismatch for {key}"
                    )
            if candidate.sequence_version is None:
                if _text(row.get("sequence_version")):
                    raise ValueError(
                        f"candidate CSV sequence version mismatch for {key}"
                    )
            elif _int(row.get("sequence_version")) != candidate.sequence_version:
                raise ValueError(
                    f"candidate CSV sequence version mismatch for {key}"
                )
            if _text(row.get("sequence_sha256")).lower() != candidate.sequence_sha256:
                raise ValueError(f"candidate CSV sequence mismatch for {key}")
            if _text(row.get("sequence")) != candidate.sequence:
                raise ValueError(f"candidate CSV sequence payload mismatch for {key}")
    if set(by_key) != expected_keys:
        extras = sorted(set(by_key) - expected_keys)
        raise ValueError(f"candidate CSV has rows outside canonical selection: {extras}")
    return [by_key[key] for key in sorted(by_key)]


def _failure_result(
    *,
    status: str,
    solution_id: int,
    expansion_depth: int,
    solution_fingerprint_value: str,
    chassis_key: str,
    required_step_indexes: list[int],
    max_sets: int,
    max_search_nodes: int,
    reason: str,
    source_artifacts: Mapping[str, str],
) -> dict[str, Any]:
    candidate_pool_fingerprint = stable_json_hash(
        {
            "algorithm_version": MAIN_ENZYME_SETS_ALGORITHM_VERSION,
            "status": status,
            "solution_fingerprint": solution_fingerprint_value,
            "chassis_key": chassis_key,
            "reason": reason,
        }
    )
    input_fingerprint = _input_fingerprint(
        candidate_pool_fingerprint,
        max_sets=max_sets,
        max_search_nodes=max_search_nodes,
    )
    result = MainEnzymeSetsResult(
        ok=False,
        status=status,
        selected_solution_id=solution_id,
        expansion_depth=expansion_depth,
        solution_fingerprint=solution_fingerprint_value,
        chassis_key=chassis_key,
        parameters={
            "max_sets": max_sets,
            "max_search_nodes": max_search_nodes,
            "candidate_scope": "top_n_shortlist",
        },
        candidate_pool_fingerprint=candidate_pool_fingerprint,
        input_fingerprint=input_fingerprint,
        required_step_indexes=required_step_indexes,
        minimum_protein_count=None,
        search_complete=False,
        search_nodes=0,
        sets=[],
        uncovered_step_indexes=required_step_indexes,
        blocking_reasons=[reason],
        warnings=[],
        source_artifacts=dict(source_artifacts),
    )
    return result.model_dump(mode="json")


def _write_set_outputs(
    output_dir: str | Path,
    result: Mapping[str, Any],
) -> dict[str, str]:
    paths = main_enzyme_set_paths(output_dir)
    write_json_atomic(paths["main_enzyme_sets_json"], dict(result))
    summary_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    for enzyme_set in result.get("sets", []):
        metrics = enzyme_set["metrics"]
        summary_rows.append(
            {
                "set_id": enzyme_set["set_id"],
                "status": enzyme_set["status"],
                "protein_count": enzyme_set["protein_count"],
                "accessions": ";".join(
                    item["accession"] for item in enzyme_set["proteins"]
                ),
                "covered_step_indexes": ";".join(
                    map(str, enzyme_set["covered_step_indexes"])
                ),
                "review_required_step_indexes": ";".join(
                    map(str, enzyme_set["review_required_step_indexes"])
                ),
                "organism_count": metrics["organism_count"],
                "min_reaction_fit_score": metrics["min_reaction_fit_score"],
                "mean_reaction_fit_score": metrics["mean_reaction_fit_score"],
                "min_protein_score": metrics["min_protein_score"],
                "mean_protein_score": metrics["mean_protein_score"],
                "min_host_fit_score": metrics["min_host_fit_score"],
                "mean_host_fit_score": metrics["mean_host_fit_score"],
                "reviewed_fraction": metrics["reviewed_fraction"],
                "reaction_fit_risk_count": metrics[
                    "reaction_fit_risk_count"
                ],
                "direction_risk_count": metrics["direction_risk_count"],
                "low_direction_confidence_count": metrics[
                    "low_direction_confidence_count"
                ],
                "specificity_risk_count": metrics["specificity_risk_count"],
                "exact_specificity_count": metrics[
                    "exact_specificity_count"
                ],
                "carrier_compatibility_status": metrics[
                    "carrier_compatibility_status"
                ],
                "electron_reassessment_status": metrics[
                    "electron_reassessment_status"
                ],
                "set_fingerprint": enzyme_set["set_fingerprint"],
                "warnings": enzyme_set["warnings"],
            }
        )
        reactions_by_accession: dict[str, list[str]] = defaultdict(list)
        for assignment in enzyme_set["step_assignments"]:
            reactions_by_accession[assignment["accession"]].append(
                assignment["reaction_id"]
            )
        for protein in enzyme_set["proteins"]:
            member_rows.append(
                {
                    "set_id": enzyme_set["set_id"],
                    "set_status": enzyme_set["status"],
                    **protein,
                    "cofactors": ";".join(protein["cofactors"]),
                    "capable_step_indexes": ";".join(
                        map(str, protein["capable_step_indexes"])
                    ),
                    "assigned_step_indexes": ";".join(
                        map(str, protein["assigned_step_indexes"])
                    ),
                    "assigned_reaction_ids": ";".join(
                        reactions_by_accession[protein["accession"]]
                    ),
                }
            )
    write_csv(paths["main_enzyme_sets_csv"], summary_rows, _SUMMARY_COLUMNS)
    write_csv(
        paths["main_enzyme_set_members_csv"], member_rows, _MEMBER_COLUMNS
    )
    return {key: rel_or_abs(path) for key, path in paths.items()}


def build_main_enzyme_sets(
    manifest_path: str | Path,
    selection_path: str | Path,
    candidate_csv_path: str | Path,
    output_dir: str | Path,
    *,
    max_sets: int = 20,
    max_search_nodes: int = 1_000_000,
) -> dict[str, Any]:
    """Validate persisted inputs, enumerate sets, and write three artifacts."""

    manifest_path = Path(manifest_path)
    selection_path = Path(selection_path)
    candidate_csv_path = Path(candidate_csv_path)
    manifest = read_manifest(manifest_path)
    solution_id, steps = get_solution_steps(manifest)
    solution = manifest["solution"]
    expansion_depth = _int(solution.get("expansion_depth"))
    route_fingerprint = solution_fingerprint(solution_id, steps)
    requirements = heterologous_requirements(steps)
    required_indexes = sorted(
        _int(item.get("step_index")) for item in requirements
    )
    source_artifacts = {"manifest": rel_or_abs(manifest_path)}

    if not selection_path.exists() or not candidate_csv_path.exists():
        missing = [
            str(path)
            for path in (selection_path, candidate_csv_path)
            if not path.exists()
        ]
        result = _failure_result(
            status="source_unavailable",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key="unknown",
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason="Missing main-enzyme candidate artifact(s): " + ", ".join(missing),
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result

    source_artifacts.update(
        {
            "candidate_selection": rel_or_abs(selection_path),
            "candidate_selection_sha256": file_sha256(selection_path),
            "candidate_details": rel_or_abs(candidate_csv_path),
            "candidate_details_sha256": file_sha256(candidate_csv_path),
        }
    )
    try:
        raw_selection = json.loads(selection_path.read_text(encoding="utf-8"))
        selection = MainEnzymeSelectionResult.model_validate(raw_selection)
    except Exception as exc:
        result = _failure_result(
            status="stale_input",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key="unknown",
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason=f"Invalid main-enzyme selection artifact: {exc}",
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result

    stale_reasons: list[str] = []
    if selection.selected_solution_id != solution_id:
        stale_reasons.append("selected solution ID changed")
    if selection.expansion_depth != expansion_depth:
        stale_reasons.append("selected expansion depth changed")
    if selection.solution_fingerprint != route_fingerprint:
        stale_reasons.append("selected route fingerprint changed")
    if stale_reasons:
        result = _failure_result(
            status="stale_input",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key=selection.chassis_key,
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason="; ".join(stale_reasons) + "; rerun main-enzyme",
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result
    if selection.status == "source_unavailable":
        result = _failure_result(
            status="source_unavailable",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key=selection.chassis_key,
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason="Main-enzyme candidate retrieval was unavailable; rerun main-enzyme.",
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result

    rows = read_csv(candidate_csv_path)
    try:
        shortlist_rows = _validate_selection_csv(selection, rows)
    except ValueError as exc:
        result = _failure_result(
            status="stale_input",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key=selection.chassis_key,
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason=f"Candidate JSON/CSV mismatch: {exc}; rerun main-enzyme",
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result

    detail_fingerprint = shortlist_decision_fingerprint(shortlist_rows)
    if not selection.shortlist_decision_fingerprint:
        result = _failure_result(
            status="stale_input",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key=selection.chassis_key,
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason=(
                "Candidate selection predates shortlist decision binding; "
                "rerun main-enzyme"
            ),
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result
    if selection.shortlist_decision_fingerprint != detail_fingerprint:
        result = _failure_result(
            status="stale_input",
            solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint_value=route_fingerprint,
            chassis_key=selection.chassis_key,
            required_step_indexes=required_indexes,
            max_sets=max_sets,
            max_search_nodes=max_search_nodes,
            reason=(
                "Candidate detail fingerprint does not match selection JSON; "
                "rerun main-enzyme"
            ),
            source_artifacts=source_artifacts,
        )
        result["output_files"] = _write_set_outputs(output_dir, result)
        return result

    source_artifacts["candidate_selection_fingerprint"] = (
        main_enzyme_selection_fingerprint(selection)
    )
    source_artifacts["candidate_detail_fingerprint"] = detail_fingerprint
    result = build_main_enzyme_sets_from_rows(
        solution_id=solution_id,
        expansion_depth=expansion_depth,
        solution_fingerprint_value=route_fingerprint,
        chassis_key=selection.chassis_key,
        required_steps=requirements,
        candidate_rows=shortlist_rows,
        electron_inference=(
            manifest.get("electron_inference")
            if isinstance(manifest.get("electron_inference"), Mapping)
            else {}
        ),
        max_sets=max_sets,
        max_search_nodes=max_search_nodes,
        source_artifacts=source_artifacts,
    )
    output_files = _write_set_outputs(output_dir, result)
    return {**result, "output_files": output_files}


def run_main_enzyme_sets(config: Any) -> dict[str, Any]:
    """Run the independent combination stage with project ``RunConfig``."""

    output_dir = Path(config.project_output_path) / "main_protein_selection"
    result = build_main_enzyme_sets(
        manifest_path=Path(config.manifest_output_path),
        selection_path=output_dir / "main_enzyme_selection.json",
        candidate_csv_path=output_dir / "step_main_enzyme_candidates.csv",
        output_dir=output_dir,
        max_sets=_int(getattr(config, "max_sets", 20), 20),
        max_search_nodes=_int(
            getattr(config, "max_search_nodes", 1_000_000), 1_000_000
        ),
    )
    summary = {
        "运行成功": result["ok"],
        "状态": result["status"],
        "路径编号": result["selected_solution_id"],
        "最少主酶数量": result["minimum_protein_count"],
        "组合数量": len(result["sets"]),
        "搜索是否完整": result["search_complete"],
        "结果文件": result.get("output_files", {}),
        "阻塞原因": result.get("blocking_reasons", []),
        "警告": result.get("warnings", []),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "MAIN_ENZYME_SETS_FILENAME",
    "MAIN_ENZYME_SETS_SUMMARY_FILENAME",
    "MAIN_ENZYME_SET_MEMBERS_FILENAME",
    "build_main_enzyme_sets",
    "build_main_enzyme_sets_from_rows",
    "candidate_pool_fingerprint_from_rows",
    "main_enzyme_selection_fingerprint",
    "main_enzyme_set_paths",
    "run_main_enzyme_sets",
    "shortlist_decision_fingerprint",
]
