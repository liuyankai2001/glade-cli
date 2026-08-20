"""Deterministic whole-cassette recommendation from remote part facts."""

from __future__ import annotations

import hashlib
import heapq
import itertools
import math
from collections import Counter
from collections.abc import Callable
from importlib.metadata import version
from typing import Any

from dnachisel import AvoidPattern, DnaOptimizationProblem, EnforceGCContent

from src.expression_box.config import (
    EXPRESSION_PARTS_CANDIDATE_POOL_MAX,
    EXPRESSION_PARTS_CANDIDATE_POOL_MIN,
    EXPRESSION_PARTS_CANDIDATE_POOL_MULTIPLIER,
    EXPRESSION_SUCCESS_ASSEMBLY_WEIGHT,
    EXPRESSION_SUCCESS_EVIDENCE_WEIGHT,
    EXPRESSION_SUCCESS_MIN_SCORE,
    EXPRESSION_SUCCESS_PROMOTER_WEIGHT,
    EXPRESSION_SUCCESS_TERMINATOR_WEIGHT,
    EXPRESSION_SUCCESS_TRANSLATION_WEIGHT,
    RBS_CONTEXT_CDS_PREFIX_NT,
    RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT,
    RBS_SHORTLIST_PER_STRENGTH,
)
from src.expression_box.ostir_adapter import predict_rbs_context
from src.expression_box.parts_models import (
    ExpressionCds,
    ExpressionPartsCassette,
    ExpressionPartsContext,
    PartCandidate,
    PartsSnapshot,
    RbsPrediction,
)
from src.protein_to_cds.sequence_constraints import (
    DEFAULT_FORBIDDEN_MOTIFS,
    DNA_ALPHABET,
    HOMOPOLYMER_LIMIT,
    LOCAL_GC_MAX,
    LOCAL_GC_MIN,
    LOCAL_GC_WINDOW_NT,
    gc_fraction,
    local_gc_values,
    max_homopolymer_length,
    motif_hits,
    reverse_complement,
)


_CONFIDENCE_SCORE = {"high": 3, "medium": 2, "low": 1, "unknown": 0, "": 0}
_EVIDENCE_SCORE = {"A": 3, "B": 2, "C": 1}
_CONFIDENCE_RELIABILITY = {
    "high": 1.0,
    "medium": 0.85,
    "low": 0.65,
    "unknown": 0.50,
    "": 0.50,
}
_EVIDENCE_RELIABILITY = {"A": 1.0, "B": 0.85, "C": 0.65}
_STRATEGIES = (
    {
        "design_id": 1,
        "strategy": "balanced",
        "name": "均衡表达方案",
        "strength": "medium",
        "target_percentile": 50.0,
        "burden": "moderate",
    },
    {
        "design_id": 2,
        "strategy": "low_burden",
        "name": "低负担方案",
        "strength": "low",
        "target_percentile": 20.0,
        "burden": "low",
    },
    {
        "design_id": 3,
        "strategy": "high_expression",
        "name": "较高表达方案",
        "strength": "high",
        "target_percentile": 80.0,
        "burden": "high",
    },
)


def _quantile(values: list[float], percentile: float) -> float:
    if not values:
        raise ValueError("cannot calculate a quantile from an empty list")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * percentile / 100.0
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _metadata_penalty(part: PartCandidate) -> float:
    evidence = _EVIDENCE_SCORE.get(part.evidence_grade, 0)
    host_penalty = 0.0 if part.host_match_kind == "exact" else 0.10
    confidence_values = (
        part.role_confidence,
        part.sequence_type_confidence,
        part.host_confidence,
        part.strength_confidence,
    )
    confidence_penalty = sum(
        3 - _CONFIDENCE_SCORE.get(value, 0) for value in confidence_values
    ) * 0.02
    return (
        (3 - evidence) * 0.10
        + host_penalty
        + confidence_penalty
        + min(len(part.warnings), 5) * 0.05
    )


def _activity_distance(part: PartCandidate, target: float) -> float:
    if part.activity_percentile is None:
        return 1.0
    return abs(part.activity_percentile - target) / 100.0


def _metadata_rank(part: PartCandidate, target: float) -> tuple[Any, ...]:
    return (
        _activity_distance(part, target) + _metadata_penalty(part),
        -_EVIDENCE_SCORE.get(part.evidence_grade, 0),
        0 if part.host_match_kind == "exact" else 1,
        -_CONFIDENCE_SCORE.get(part.strength_confidence, 0),
        len(part.warnings),
        part.part_id,
    )


def _part_payload(part: PartCandidate) -> dict[str, Any]:
    return {
        "part_id": part.part_id,
        "role": part.role,
        "sequence": part.sequence,
        "sequence_sha256": part.sequence_sha256,
        "length_bp": len(part.sequence),
        "strength": part.strength,
        "regulation": part.regulation,
        "direction": part.direction,
        "host_match_kind": part.host_match_kind,
        "evidence_grade": part.evidence_grade,
        "role_confidence": part.role_confidence,
        "sequence_type_confidence": part.sequence_type_confidence,
        "host_confidence": part.host_confidence,
        "strength_confidence": part.strength_confidence,
        "regulation_confidence": part.regulation_confidence,
        "activity_value": part.activity_value,
        "activity_percentile": part.activity_percentile,
        "activity_dataset": part.activity_dataset,
        "activity_unit": part.activity_unit,
        "activity_context": part.activity_context,
        "activity_source": part.activity_source,
        "source": part.source,
        "evidence": part.evidence,
        "warnings": list(part.warnings),
    }


def _context_sequence(
    cassette: ExpressionPartsCassette,
    cds_index: int,
    rbs: PartCandidate,
) -> tuple[str, int]:
    cds = cassette.cds[cds_index]
    previous = cassette.cds[cds_index - 1] if cds_index > 0 else None
    upstream = (
        previous.sequence[-RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT:]
        if previous is not None
        else ""
    )
    sequence = (
        upstream
        + rbs.sequence
        + cds.sequence[:RBS_CONTEXT_CDS_PREFIX_NT]
    )
    intended_start = len(upstream) + len(rbs.sequence) + 1
    return sequence, intended_start


def _rbs_shortlist(
    candidates: tuple[PartCandidate, ...],
) -> dict[str, list[PartCandidate]]:
    targets = {"low": 20.0, "medium": 50.0, "high": 80.0}
    result: dict[str, list[PartCandidate]] = {}
    for strength, target in targets.items():
        values = [
            item
            for item in candidates
            if item.role == "rbs" and item.strength == strength
            and _part_sequence_safe(item)
        ]
        values.sort(key=lambda item: _metadata_rank(item, target))
        result[strength] = values[:RBS_SHORTLIST_PER_STRENGTH]
        if not result[strength]:
            raise ValueError(f"no usable Milvus RBS candidates for strength {strength}")
    return result


def _part_sequence_safe(part: PartCandidate) -> bool:
    return (
        not motif_hits(part.sequence, DEFAULT_FORBIDDEN_MOTIFS)
        and max_homopolymer_length(part.sequence) < HOMOPOLYMER_LIMIT
    )


def _predict_all_rbs(
    context: ExpressionPartsContext,
    shortlists: dict[str, list[PartCandidate]],
    predictor: Callable[..., RbsPrediction],
) -> tuple[
    dict[tuple[int, str, str], RbsPrediction],
    list[dict[str, str]],
]:
    unique_rbs = {
        item.part_id: item
        for values in shortlists.values()
        for item in values
    }
    predictions: dict[tuple[int, str, str], RbsPrediction] = {}
    failures: list[dict[str, str]] = []
    for cassette in context.cassettes:
        for cds_index, cds in enumerate(cassette.cds):
            for rbs in unique_rbs.values():
                sequence, intended_start = _context_sequence(
                    cassette,
                    cds_index,
                    rbs,
                )
                try:
                    prediction = predictor(
                        sequence=sequence,
                        intended_start_position=intended_start,
                        accession=cds.accession,
                        part_id=rbs.part_id,
                    )
                except Exception as exc:
                    failures.append(
                        {
                            "cassette_index": str(cassette.cassette_index),
                            "accession": cds.accession,
                            "part_id": rbs.part_id,
                            "error_type": type(exc).__name__,
                            "message": str(exc),
                        }
                    )
                    continue
                predictions[(cassette.cassette_index, cds.accession, rbs.part_id)] = (
                    prediction
                )
    return predictions, failures


def _rbs_options(
    *,
    cassette: ExpressionPartsCassette,
    cds: ExpressionCds,
    strategy: dict[str, Any],
    shortlists: dict[str, list[PartCandidate]],
    predictions: dict[tuple[int, str, str], RbsPrediction],
) -> list[tuple[float, PartCandidate, RbsPrediction]]:
    all_predictions = [
        prediction
        for (cassette_index, accession, _), prediction in predictions.items()
        if cassette_index == cassette.cassette_index and accession == cds.accession
    ]
    if not all_predictions:
        return []
    target_log = _quantile(
        [math.log10(item.expression) for item in all_predictions],
        float(strategy["target_percentile"]),
    )
    options: list[tuple[float, PartCandidate, RbsPrediction]] = []
    for part in shortlists[str(strategy["strength"])]:
        prediction = predictions.get(
            (cassette.cassette_index, cds.accession, part.part_id)
        )
        if prediction is None:
            continue
        cost = (
            abs(math.log10(prediction.expression) - target_log)
            + 0.05 * prediction.unintended_start_count
            + _metadata_penalty(part)
        )
        options.append((cost, part, prediction))
    options.sort(
        key=lambda item: (
            item[0],
            item[1].part_id,
        )
    )
    return options[:6]


def _audit_sequence(sequence: str) -> dict[str, Any]:
    normalized = str(sequence or "").strip().upper()
    valid_alphabet = bool(normalized) and set(normalized).issubset(DNA_ALPHABET)
    if not valid_alphabet:
        return {
            "engine": "DNA Chisel",
            "engine_version": version("dnachisel"),
            "gate_status": "FAIL",
            "checks": {"valid_alphabet": False},
            "failed_checks": ["valid_alphabet"],
        }
    constrained_motifs: set[str] = set()
    for motif in DEFAULT_FORBIDDEN_MOTIFS.values():
        constrained_motifs.add(motif)
        constrained_motifs.add(reverse_complement(motif))
    constraints: list[Any] = [
        EnforceGCContent(mini=0.30, maxi=0.70),
        EnforceGCContent(
            mini=LOCAL_GC_MIN,
            maxi=LOCAL_GC_MAX,
            window=LOCAL_GC_WINDOW_NT,
        ),
        *(AvoidPattern(motif) for motif in sorted(constrained_motifs)),
        *(AvoidPattern(base * HOMOPOLYMER_LIMIT) for base in "ACGT"),
    ]
    problem = DnaOptimizationProblem(
        normalized,
        constraints=constraints,
        objectives=[],
        logger=None,
    )
    dnachisel_pass = problem.all_constraints_pass()
    local_values = local_gc_values(normalized)
    hits = motif_hits(normalized, DEFAULT_FORBIDDEN_MOTIFS)
    checks = {
        "valid_alphabet": True,
        "global_gc_pass": 0.30 <= gc_fraction(normalized) <= 0.70,
        "local_gc_pass": all(
            LOCAL_GC_MIN <= value <= LOCAL_GC_MAX for value in local_values
        ),
        "forbidden_motif_pass": not hits,
        "homopolymer_pass": (
            max_homopolymer_length(normalized) < HOMOPOLYMER_LIMIT
        ),
        "dnachisel_constraints_pass": dnachisel_pass,
    }
    return {
        "engine": "DNA Chisel",
        "engine_version": version("dnachisel"),
        "gate_status": "PASS" if all(checks.values()) else "FAIL",
        "sequence_sha256": hashlib.sha256(normalized.encode("utf-8")).hexdigest(),
        "length_nt": len(normalized),
        "gc_percent": round(100.0 * gc_fraction(normalized), 8),
        "local_gc_min_percent": (
            round(100.0 * min(local_values), 8) if local_values else None
        ),
        "local_gc_max_percent": (
            round(100.0 * max(local_values), 8) if local_values else None
        ),
        "forbidden_site_hits": hits,
        "max_homopolymer": max_homopolymer_length(normalized),
        "checks": checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
    }


def _rank_role_candidates(
    candidates: tuple[PartCandidate, ...],
    *,
    role: str,
    strength: str,
    target: float,
    limit: int,
) -> list[tuple[float, PartCandidate]]:
    rows = [
        item
        for item in candidates
        if item.role == role and item.strength == strength
        and _part_sequence_safe(item)
    ]
    ranked = [
        (_activity_distance(item, target) + _metadata_penalty(item), item)
        for item in rows
    ]
    ranked.sort(key=lambda item: (item[0], item[1].part_id))
    return ranked[:limit]


def _candidate_pool_limit(category_quota: int) -> int:
    return min(
        EXPRESSION_PARTS_CANDIDATE_POOL_MAX,
        max(
            EXPRESSION_PARTS_CANDIDATE_POOL_MIN,
            EXPRESSION_PARTS_CANDIDATE_POOL_MULTIPLIER * category_quota,
        ),
    )


def _partial_signature(candidate: dict[str, Any]) -> tuple[str, tuple[str, ...]]:
    return (
        candidate["promoter"].part_id,
        tuple(gene[1].part_id for gene in candidate["genes"]),
    )


def _cassette_signature(payload: dict[str, Any]) -> tuple[str, ...]:
    return (
        str(payload["promoter"]["part_id"]),
        *(str(gene["rbs"]["part_id"]) for gene in payload["genes"]),
        str(payload["terminator"]["part_id"]),
    )


def _design_signature(cassettes: list[dict[str, Any]]) -> tuple[str, ...]:
    signature: list[str] = []
    for cassette in cassettes:
        signature.extend(_cassette_signature(cassette))
    return tuple(signature)


def _cassette_candidates(
    *,
    cassette: ExpressionPartsCassette,
    strategy: dict[str, Any],
    candidates: tuple[PartCandidate, ...],
    shortlists: dict[str, list[PartCandidate]],
    predictions: dict[tuple[int, str, str], RbsPrediction],
    limit: int,
) -> list[dict[str, Any]]:
    promoter_options = _rank_role_candidates(
        candidates,
        role="promoter",
        strength=str(strategy["strength"]),
        target=float(strategy["target_percentile"]),
        limit=6,
    )
    terminator_options = _rank_role_candidates(
        candidates,
        role="terminator",
        strength="high",
        target=95.0,
        limit=8,
    )
    if not promoter_options:
        raise ValueError(
            f"no promoter candidates for strength {strategy['strength']}"
        )
    if not terminator_options:
        raise ValueError("no high-strength terminator candidates")

    beams: list[dict[str, Any]] = []
    for base_cost, promoter in promoter_options:
        beams.append(
            {
                "cost": base_cost,
                "promoter": promoter,
                "genes": [],
                "sequence": promoter.sequence,
                "hashes": {promoter.sequence_sha256},
            }
        )
    for cds in cassette.cds:
        options = _rbs_options(
            cassette=cassette,
            cds=cds,
            strategy=strategy,
            shortlists=shortlists,
            predictions=predictions,
        )
        if not options:
            raise ValueError(
                f"no OSTIR-valid {strategy['strength']} RBS candidates for {cds.accession}"
            )
        expanded: list[dict[str, Any]] = []
        for beam, (base_cost, rbs, prediction) in itertools.product(beams, options):
            repeated = rbs.sequence_sha256 in beam["hashes"]
            expanded.append(
                {
                    **beam,
                    "cost": beam["cost"] + base_cost + (0.25 if repeated else 0.0),
                    "genes": [*beam["genes"], (cds, rbs, prediction)],
                    "sequence": beam["sequence"] + rbs.sequence + cds.sequence,
                    "hashes": {*beam["hashes"], rbs.sequence_sha256},
                }
            )
        expanded.sort(
            key=lambda item: (
                item["cost"],
                item["promoter"].part_id,
                tuple(gene[1].part_id for gene in item["genes"]),
            )
        )
        beams = []
        seen_partial: set[tuple[str, tuple[str, ...]]] = set()
        for item in expanded:
            signature = _partial_signature(item)
            if signature in seen_partial:
                continue
            seen_partial.add(signature)
            beams.append(item)
            if len(beams) >= limit:
                break

    raw_finals: list[tuple[float, dict[str, Any], PartCandidate]] = []
    for beam, (base_cost, terminator) in itertools.product(beams, terminator_options):
        repeated = terminator.sequence_sha256 in beam["hashes"]
        total_cost = beam["cost"] + base_cost + (0.25 if repeated else 0.0)
        raw_finals.append((total_cost, beam, terminator))
    raw_finals.sort(
        key=lambda item: (
            item[0],
            item[1]["promoter"].part_id,
            tuple(gene[1].part_id for gene in item[1]["genes"]),
            item[2].part_id,
        )
    )

    finals: list[dict[str, Any]] = []
    seen_final: set[tuple[str, ...]] = set()
    for cost, selected, terminator in raw_finals:
        raw_signature = (
            selected["promoter"].part_id,
            *(gene[1].part_id for gene in selected["genes"]),
            terminator.part_id,
        )
        if raw_signature in seen_final:
            continue
        seen_final.add(raw_signature)
        audit = _audit_sequence(selected["sequence"] + terminator.sequence)
        if audit["gate_status"] != "PASS":
            continue
        part_hashes = (
            selected["promoter"].sequence_sha256,
            *(gene[1].sequence_sha256 for gene in selected["genes"]),
            terminator.sequence_sha256,
        )
        payload = {
            "cassette_index": cassette.cassette_index,
            "promoter": _part_payload(selected["promoter"]),
            "genes": [
                {
                    "accession": cds.accession,
                    "cds_sequence_sha256": cds.sequence_sha256,
                    "cds_length_nt": len(cds.sequence),
                    "rbs": _part_payload(rbs),
                    "ostir": {
                        "translation_initiation_rate": prediction.expression,
                        "d_g_total": prediction.d_g_total,
                        "intended_start_position": (
                            prediction.intended_start_position
                        ),
                        "unintended_start_count": prediction.unintended_start_count,
                        "context_sha256": prediction.context_sha256,
                    },
                }
                for cds, rbs, prediction in selected["genes"]
            ],
            "terminator": _part_payload(terminator),
            "assembled_sequence": {
                "stored": False,
                "length_nt": audit["length_nt"],
                "sequence_sha256": audit["sequence_sha256"],
            },
            "sequence_audit": audit,
        }
        finals.append(
            {
                "cost": cost,
                "payload": payload,
                "signature": _cassette_signature(payload),
                "part_hashes": part_hashes,
            }
        )
        if len(finals) >= limit:
            break
    if not finals:
        raise ValueError(
            f"no sequence-safe part combination for cassette {cassette.cassette_index}"
        )
    return finals


def _combine_design_candidate(
    beam: dict[str, Any],
    cassette_candidate: dict[str, Any],
) -> dict[str, Any]:
    cross_repeats = sum(
        sequence_hash in beam["hashes"]
        for sequence_hash in cassette_candidate["part_hashes"]
    )
    cassettes = [*beam["cassettes"], cassette_candidate["payload"]]
    return {
        "cost": beam["cost"] + cassette_candidate["cost"] + 0.25 * cross_repeats,
        "cassettes": cassettes,
        "hashes": {*beam["hashes"], *cassette_candidate["part_hashes"]},
        "signature": _design_signature(cassettes),
    }


def _whole_design_candidates(
    *,
    context: ExpressionPartsContext,
    strategy: dict[str, Any],
    candidates: tuple[PartCandidate, ...],
    shortlists: dict[str, list[PartCandidate]],
    predictions: dict[tuple[int, str, str], RbsPrediction],
    category_quota: int,
) -> list[dict[str, Any]]:
    pool_limit = _candidate_pool_limit(category_quota)
    cassette_limit = min(pool_limit, max(32, category_quota * 8))
    beams: list[dict[str, Any]] = [
        {"cost": 0.0, "cassettes": [], "hashes": set(), "signature": ()}
    ]
    for cassette in context.cassettes:
        cassette_candidates = _cassette_candidates(
            cassette=cassette,
            strategy=strategy,
            candidates=candidates,
            shortlists=shortlists,
            predictions=predictions,
            limit=cassette_limit,
        )
        expanded = (
            _combine_design_candidate(beam, cassette_candidate)
            for beam in beams
            for cassette_candidate in cassette_candidates
        )
        beams = heapq.nsmallest(
            pool_limit,
            expanded,
            key=lambda item: (item["cost"], item["signature"]),
        )
        unique: list[dict[str, Any]] = []
        seen: set[tuple[str, ...]] = set()
        for item in beams:
            if item["signature"] in seen:
                continue
            seen.add(item["signature"])
            unique.append(item)
        beams = unique
    return beams


def _central_expression_score(percentile: float | None) -> float:
    if percentile is None:
        return 0.50
    bounded = min(100.0, max(0.0, float(percentile)))
    return max(0.0, 1.0 - abs(bounded - 50.0) / 50.0)


def _empirical_percentile(value: float, reference: tuple[float, ...]) -> float:
    if len(reference) <= 1:
        return 50.0
    target = math.log10(max(float(value), 1e-12))
    logs = sorted(math.log10(max(float(item), 1e-12)) for item in reference)
    lower = sum(item < target for item in logs)
    equal = sum(math.isclose(item, target, rel_tol=1e-12, abs_tol=1e-12) for item in logs)
    average_index = lower + max(equal - 1, 0) / 2.0
    return 100.0 * average_index / (len(logs) - 1)


def _translation_reference(
    predictions: dict[tuple[int, str, str], RbsPrediction],
) -> dict[tuple[int, str], tuple[float, ...]]:
    grouped: dict[tuple[int, str], list[float]] = {}
    for (cassette_index, accession, _), prediction in predictions.items():
        grouped.setdefault((cassette_index, accession), []).append(prediction.expression)
    return {key: tuple(values) for key, values in grouped.items()}


def _selected_parts(cassettes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    parts: list[dict[str, Any]] = []
    for cassette in cassettes:
        parts.append(cassette["promoter"])
        parts.extend(gene["rbs"] for gene in cassette["genes"])
        parts.append(cassette["terminator"])
    return parts


def _part_reliability(part: dict[str, Any]) -> float:
    evidence = _EVIDENCE_RELIABILITY.get(str(part.get("evidence_grade") or ""), 0.50)
    host = 1.0 if part.get("host_match_kind") == "exact" else 0.85
    confidences = [
        _CONFIDENCE_RELIABILITY.get(str(part.get(field) or ""), 0.50)
        for field in (
            "role_confidence",
            "sequence_type_confidence",
            "host_confidence",
            "strength_confidence",
            "regulation_confidence",
        )
    ]
    confidence = sum(confidences) / len(confidences)
    return 0.40 * evidence + 0.35 * host + 0.25 * confidence


def _score_expression_candidate(
    candidate: dict[str, Any],
    translation_reference: dict[tuple[int, str], tuple[float, ...]],
) -> dict[str, Any]:
    cassettes = candidate["cassettes"]
    parts = _selected_parts(cassettes)
    evidence_factor = sum(_part_reliability(part) for part in parts) / len(parts)

    promoter_factors = [
        _central_expression_score(cassette["promoter"].get("activity_percentile"))
        for cassette in cassettes
    ]
    promoter_factor = sum(promoter_factors) / len(promoter_factors)

    translation_factors: list[float] = []
    unintended_start_count = 0
    for cassette in cassettes:
        cassette_index = int(cassette["cassette_index"])
        for gene in cassette["genes"]:
            ostir = gene["ostir"]
            unintended = int(ostir["unintended_start_count"])
            unintended_start_count += unintended
            percentile = _empirical_percentile(
                float(ostir["translation_initiation_rate"]),
                translation_reference[(cassette_index, str(gene["accession"]))],
            )
            robustness = _central_expression_score(percentile)
            translation_factors.append(robustness / (1.0 + 0.50 * unintended))
    translation_factor = sum(translation_factors) / len(translation_factors)

    terminator_factors = []
    for cassette in cassettes:
        percentile = cassette["terminator"].get("activity_percentile")
        terminator_factors.append(
            0.50
            if percentile is None
            else min(1.0, max(0.0, float(percentile)) / 90.0)
        )
    terminator_factor = sum(terminator_factors) / len(terminator_factors)

    hashes = Counter(str(part["sequence_sha256"]) for part in parts)
    repeated_occurrences = sum(count - 1 for count in hashes.values() if count > 1)
    assembly_factor = 1.0 - repeated_occurrences / len(parts)

    components = {
        "evidence_reliability": round(
            EXPRESSION_SUCCESS_EVIDENCE_WEIGHT * evidence_factor, 2
        ),
        "promoter_robustness": round(
            EXPRESSION_SUCCESS_PROMOTER_WEIGHT * promoter_factor, 2
        ),
        "translation_robustness": round(
            EXPRESSION_SUCCESS_TRANSLATION_WEIGHT * translation_factor, 2
        ),
        "terminator_reliability": round(
            EXPRESSION_SUCCESS_TERMINATOR_WEIGHT * terminator_factor, 2
        ),
        "assembly_robustness": round(
            EXPRESSION_SUCCESS_ASSEMBLY_WEIGHT * assembly_factor, 2
        ),
    }
    score = round(sum(components.values()), 2)
    return {
        **candidate,
        "expression_success_score": score,
        "score_components": components,
        "unintended_start_count": unintended_start_count,
        "repeated_regulatory_part_count": repeated_occurrences,
        "passes_success_threshold": score >= EXPRESSION_SUCCESS_MIN_SCORE,
    }


def _success_rank_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -float(candidate["expression_success_score"]),
        -float(candidate["score_components"]["evidence_reliability"]),
        int(candidate["unintended_start_count"]),
        int(candidate["repeated_regulatory_part_count"]),
        candidate["signature"],
    )


def _design_warnings(cassettes: list[dict[str, Any]], burden: str) -> list[str]:
    selected_parts: list[dict[str, Any]] = []
    for cassette in cassettes:
        selected_parts.append(cassette["promoter"])
        selected_parts.extend(gene["rbs"] for gene in cassette["genes"])
        selected_parts.append(cassette["terminator"])
    warnings: list[str] = []
    hashes = Counter(item["sequence_sha256"] for item in selected_parts)
    repeated_ids = sorted(
        {
            item["part_id"]
            for item in selected_parts
            if hashes[item["sequence_sha256"]] > 1
        }
    )
    if repeated_ids:
        warnings.append(
            "the design reuses regulatory-part sequences: "
            + ", ".join(repeated_ids)
        )
    if any(item["host_match_kind"] != "exact" for item in selected_parts):
        warnings.append(
            "some parts use MG1655 lineage-transfer evidence rather than direct strain evidence"
        )
    if burden == "high":
        warnings.append(
            "high-expression settings may impose substantial host burden and require experimental validation"
        )
    return warnings


def generate_expression_parts_designs(
    context: ExpressionPartsContext,
    snapshot: PartsSnapshot,
    *,
    requested_design_count: int,
    predictor: Callable[..., RbsPrediction] = predict_rbs_context,
) -> dict[str, Any]:
    """Generate globally ranked, sequence-safe stable-expression designs."""

    shortlists = _rbs_shortlist(snapshot.candidates)
    predictions, prediction_failures = _predict_all_rbs(
        context,
        shortlists,
        predictor,
    )
    pool_basis = max(1, math.ceil(requested_design_count / len(_STRATEGIES)))
    candidate_pools: dict[str, list[dict[str, Any]]] = {}
    skipped: list[dict[str, Any]] = []
    for strategy in _STRATEGIES:
        strategy_key = str(strategy["strategy"])
        try:
            generated = _whole_design_candidates(
                context=context,
                strategy=strategy,
                candidates=snapshot.candidates,
                shortlists=shortlists,
                predictions=predictions,
                category_quota=pool_basis,
            )
            candidate_pools[strategy_key] = [
                {
                    **candidate,
                    "expression_regime": strategy_key,
                    "expression_target_percentile": strategy["target_percentile"],
                    "estimated_burden": strategy["burden"],
                }
                for candidate in generated
            ]
        except ValueError as exc:
            candidate_pools[strategy_key] = []
            skipped.append(
                {
                    "strategy": strategy_key,
                    "reason": str(exc),
                }
            )

    seen_signatures: set[tuple[str, ...]] = set()
    for strategy in _STRATEGIES:
        strategy_key = str(strategy["strategy"])
        unique_pool: list[dict[str, Any]] = []
        for candidate in candidate_pools[strategy_key]:
            if candidate["signature"] in seen_signatures:
                continue
            seen_signatures.add(candidate["signature"])
            unique_pool.append(candidate)
        candidate_pools[strategy_key] = unique_pool

    capacities = {key: len(value) for key, value in candidate_pools.items()}
    combined_candidates = [
        candidate
        for strategy in _STRATEGIES
        for candidate in candidate_pools[str(strategy["strategy"])]
    ]
    translation_reference = _translation_reference(predictions)
    scored = [
        _score_expression_candidate(candidate, translation_reference)
        for candidate in combined_candidates
    ]
    scored.sort(key=_success_rank_key)
    eligible = [candidate for candidate in scored if candidate["passes_success_threshold"]]
    selected_candidates = eligible[:requested_design_count]

    designs: list[dict[str, Any]] = []
    for design_id, selected in enumerate(selected_candidates, start=1):
        cassettes = selected["cassettes"]
        score = float(selected["expression_success_score"])
        designs.append(
            {
                "design_id": design_id,
                "rank": design_id,
                "strategy": selected["expression_regime"],
                "expression_regime": selected["expression_regime"],
                "name": f"稳定表达候选 {design_id}",
                "recommended": design_id == 1,
                "recommendation_level": (
                    "primary" if design_id == 1 else "qualified_alternative"
                ),
                "confidence": "high" if score >= 85.0 else "medium",
                "expression_target_percentile": selected[
                    "expression_target_percentile"
                ],
                "estimated_burden": selected["estimated_burden"],
                "expression_success_score": score,
                "score_type": "interpretable_heuristic_not_probability",
                "passes_success_threshold": True,
                "score_components": selected["score_components"],
                "unintended_start_count": selected["unintended_start_count"],
                "repeated_regulatory_part_count": selected[
                    "repeated_regulatory_part_count"
                ],
                "design_signature": list(selected["signature"]),
                "cassette_count": len(cassettes),
                "cassettes": cassettes,
                "warnings": _design_warnings(
                    cassettes,
                    str(selected["estimated_burden"]),
                ),
            }
        )

    if len(designs) == requested_design_count:
        status = "complete"
    elif designs:
        status = "partial"
    else:
        status = "failed"
    warnings = [
        "expression success scores are interpretable design heuristics, not experimental success probabilities",
        "recommendations estimate stable expression and do not prove enzyme activity or target-compound production",
        "promoter measurements, vector copy number, culture conditions, and final assembly scars may change expression",
    ]
    if prediction_failures:
        warnings.append(
            f"OSTIR rejected {len(prediction_failures)} candidate-context predictions"
        )
    if len(designs) < requested_design_count:
        warnings.append(
            f"requested {requested_design_count} designs but only {len(designs)} unique "
            f"designs reached the minimum success score {EXPRESSION_SUCCESS_MIN_SCORE:g}"
        )
    selected_scores = [float(design["expression_success_score"]) for design in designs]
    return {
        "status": status,
        "designs": designs,
        "ranking": {
            "method": "stable_expression_success_heuristic.v1",
            "score_type": "interpretable_heuristic_not_probability",
            "requested_design_count": requested_design_count,
            "candidate_pool_counts": capacities,
            "candidate_pool_limit_per_regime": _candidate_pool_limit(pool_basis),
            "hard_gate_pass_count": len(scored),
            "eligible_count": len(eligible),
            "excluded_below_threshold_count": len(scored) - len(eligible),
            "minimum_success_score": EXPRESSION_SUCCESS_MIN_SCORE,
            "selected_score_range": {
                "highest": max(selected_scores) if selected_scores else None,
                "lowest": min(selected_scores) if selected_scores else None,
            },
            "component_weights": {
                "evidence_reliability": EXPRESSION_SUCCESS_EVIDENCE_WEIGHT,
                "promoter_robustness": EXPRESSION_SUCCESS_PROMOTER_WEIGHT,
                "translation_robustness": EXPRESSION_SUCCESS_TRANSLATION_WEIGHT,
                "terminator_reliability": EXPRESSION_SUCCESS_TERMINATOR_WEIGHT,
                "assembly_robustness": EXPRESSION_SUCCESS_ASSEMBLY_WEIGHT,
            },
        },
        "skipped_strategies": skipped,
        "prediction_failure_count": len(prediction_failures),
        "prediction_failures": prediction_failures,
        "warnings": warnings,
    }


__all__ = ["generate_expression_parts_designs"]
