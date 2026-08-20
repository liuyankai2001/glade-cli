"""Transparent intrinsic-expression burden estimates for plasmid matching."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from typing import Any


EXPRESSION_BURDEN_SCHEMA_VERSION = "expression_burden.v1"
EXPRESSION_BURDEN_MODEL_VERSION = "intrinsic_expression_burden.v1.0.0"
REFERENCE_LOAD_UNITS = 0.75
REFERENCE_CDS_LENGTH_NT = 1_500
MIN_LENGTH_FACTOR = 0.25
MAX_LENGTH_FACTOR = 2.0
RBS_PERCENTILE_WEIGHT = 0.40
OSTIR_PERCENTILE_WEIGHT = 0.60
UNINTENDED_START_WEIGHT = 0.05
UNINTENDED_START_CAP = 5
MIN_HIGH_CONFIDENCE_REFERENCE_COUNT = 6


def stable_payload_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def empirical_log_percentile(value: float, reference: tuple[float, ...]) -> float:
    """Return a deterministic mid-rank percentile on a log10 scale."""

    if value <= 0:
        raise ValueError("OSTIR translation initiation rate must be positive")
    positive = tuple(float(item) for item in reference if float(item) > 0)
    if not positive:
        raise ValueError("OSTIR reference distribution is empty")
    if len(positive) == 1:
        return 50.0
    target = math.log10(value)
    logs = sorted(math.log10(item) for item in positive)
    lower = sum(item < target for item in logs)
    equal = sum(
        math.isclose(item, target, rel_tol=1e-12, abs_tol=1e-12)
        for item in logs
    )
    average_index = lower + max(equal - 1, 0) / 2.0
    return 100.0 * average_index / (len(logs) - 1)


def burden_level(score: float) -> str:
    if score < 35.0:
        return "low"
    if score < 65.0:
        return "moderate"
    return "high"


def burden_model_parameters() -> dict[str, Any]:
    return {
        "schema_version": EXPRESSION_BURDEN_SCHEMA_VERSION,
        "model_version": EXPRESSION_BURDEN_MODEL_VERSION,
        "reference_load_units": REFERENCE_LOAD_UNITS,
        "reference_cds_length_nt": REFERENCE_CDS_LENGTH_NT,
        "length_factor_bounds": [MIN_LENGTH_FACTOR, MAX_LENGTH_FACTOR],
        "translation_percentile_weights": {
            "rbs_activity": RBS_PERCENTILE_WEIGHT,
            "ostir_context": OSTIR_PERCENTILE_WEIGHT,
        },
        "unintended_start_weight": UNINTENDED_START_WEIGHT,
        "unintended_start_cap": UNINTENDED_START_CAP,
        "level_thresholds": {
            "low_max_exclusive": 35.0,
            "moderate_max_exclusive": 65.0,
        },
        "minimum_high_confidence_reference_count": (
            MIN_HIGH_CONFIDENCE_REFERENCE_COUNT
        ),
        "score_type": "interpretable_heuristic_not_probability",
    }


def _optional_percentile(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"activity percentile is invalid: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"activity percentile is not finite: {value!r}")
    return min(100.0, max(0.0, parsed))


def _length_factor(length_nt: int) -> float:
    return min(
        MAX_LENGTH_FACTOR,
        max(MIN_LENGTH_FACTOR, length_nt / REFERENCE_CDS_LENGTH_NT),
    )


def _burden_score(raw_load_units: float) -> float:
    if raw_load_units < 0:
        raise ValueError("raw expression load cannot be negative")
    if raw_load_units == 0:
        return 0.0
    return 100.0 * raw_load_units / (raw_load_units + REFERENCE_LOAD_UNITS)


def calculate_expression_burden(
    cassettes: list[dict[str, Any]],
    translation_reference: Mapping[tuple[int, str], tuple[float, ...]],
    *,
    fallback_promoter_percentile: float,
) -> dict[str, Any]:
    """Calculate a complete, auditable burden record for one design."""

    if not cassettes:
        raise ValueError("expression burden requires at least one cassette")
    fallback = _optional_percentile(fallback_promoter_percentile)
    if fallback is None:
        raise ValueError("fallback promoter percentile is required")

    warnings: list[str] = []
    gene_metrics: list[dict[str, Any]] = []
    total_cds_length_nt = 0
    minimum_reference_count: int | None = None
    used_fallback = False
    for cassette in cassettes:
        cassette_index = int(cassette.get("cassette_index") or 0)
        if cassette_index < 1:
            raise ValueError("cassette_index must be positive")
        promoter = cassette.get("promoter")
        if not isinstance(promoter, Mapping):
            raise ValueError("cassette promoter is missing")
        promoter_percentile = _optional_percentile(
            promoter.get("activity_percentile")
        )
        promoter_source = "measured"
        if promoter_percentile is None:
            promoter_percentile = fallback
            promoter_source = "strategy_target_fallback"
            used_fallback = True
            warnings.append(
                f"cassette {cassette_index} promoter percentile used the strategy target fallback"
            )
        genes = cassette.get("genes")
        if not isinstance(genes, list) or not genes:
            raise ValueError(f"cassette {cassette_index} has no genes")
        for gene in genes:
            if not isinstance(gene, Mapping):
                raise ValueError("cassette contains an invalid gene")
            accession = str(gene.get("accession") or "").strip()
            if not accession:
                raise ValueError("gene accession is missing")
            try:
                cds_length_nt = int(gene.get("cds_length_nt") or 0)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid CDS length for {accession}") from exc
            if cds_length_nt < 1:
                raise ValueError(f"invalid CDS length for {accession}")
            rbs = gene.get("rbs")
            ostir = gene.get("ostir")
            if not isinstance(rbs, Mapping) or not isinstance(ostir, Mapping):
                raise ValueError(f"gene {accession} is missing RBS or OSTIR data")
            rbs_percentile = _optional_percentile(rbs.get("activity_percentile"))
            try:
                tir = float(ostir.get("translation_initiation_rate") or 0.0)
                unintended = int(ostir.get("unintended_start_count") or 0)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid OSTIR data for {accession}") from exc
            if tir <= 0 or unintended < 0:
                raise ValueError(f"invalid OSTIR data for {accession}")
            reference = tuple(
                float(item)
                for item in translation_reference.get(
                    (cassette_index, accession), ()
                )
                if float(item) > 0
            )
            if reference:
                ostir_percentile = empirical_log_percentile(tir, reference)
            else:
                ostir_percentile = 50.0
                used_fallback = True
                warnings.append(
                    f"gene {accession} used the neutral OSTIR percentile fallback"
                )
            reference_count = len(reference)
            minimum_reference_count = (
                reference_count
                if minimum_reference_count is None
                else min(minimum_reference_count, reference_count)
            )
            if rbs_percentile is None:
                combined_translation_percentile = ostir_percentile
                translation_source = "ostir_only_fallback"
                used_fallback = True
                warnings.append(
                    f"gene {accession} has no RBS activity percentile; OSTIR alone was used"
                )
            else:
                combined_translation_percentile = (
                    RBS_PERCENTILE_WEIGHT * rbs_percentile
                    + OSTIR_PERCENTILE_WEIGHT * ostir_percentile
                )
                translation_source = "rbs_and_ostir"
            length_factor = _length_factor(cds_length_nt)
            unintended_factor = 1.0 + UNINTENDED_START_WEIGHT * min(
                unintended, UNINTENDED_START_CAP
            )
            load_units = (
                promoter_percentile
                / 100.0
                * combined_translation_percentile
                / 100.0
                * length_factor
                * unintended_factor
            )
            total_cds_length_nt += cds_length_nt
            gene_metrics.append(
                {
                    "cassette_index": cassette_index,
                    "accession": accession,
                    "promoter_activity_percentile": round(
                        promoter_percentile, 8
                    ),
                    "promoter_percentile_source": promoter_source,
                    "rbs_activity_percentile": (
                        round(rbs_percentile, 8)
                        if rbs_percentile is not None
                        else None
                    ),
                    "ostir_translation_initiation_rate": round(tir, 8),
                    "ostir_expression_percentile": round(
                        ostir_percentile, 8
                    ),
                    "ostir_reference_count": reference_count,
                    "combined_translation_percentile": round(
                        combined_translation_percentile, 8
                    ),
                    "translation_percentile_source": translation_source,
                    "cds_length_nt": cds_length_nt,
                    "cds_length_factor": round(length_factor, 8),
                    "unintended_start_count": unintended,
                    "unintended_start_factor": round(unintended_factor, 8),
                    "load_units": round(load_units, 8),
                }
            )
    raw_load_units = round(
        sum(float(item["load_units"]) for item in gene_metrics), 8
    )
    score = round(_burden_score(raw_load_units), 2)
    minimum_count = minimum_reference_count or 0
    if minimum_count < MIN_HIGH_CONFIDENCE_REFERENCE_COUNT:
        warnings.append(
            "at least one gene has fewer than "
            f"{MIN_HIGH_CONFIDENCE_REFERENCE_COUNT} successful OSTIR reference predictions"
        )
    confidence = (
        "high"
        if not used_fallback
        and minimum_count >= MIN_HIGH_CONFIDENCE_REFERENCE_COUNT
        else "medium"
    )
    payload: dict[str, Any] = {
        "schema_version": EXPRESSION_BURDEN_SCHEMA_VERSION,
        "model_version": EXPRESSION_BURDEN_MODEL_VERSION,
        "score_type": "interpretable_heuristic_not_probability",
        "score": score,
        "level": burden_level(score),
        "confidence": confidence,
        "raw_load_units": raw_load_units,
        "reference_load_units": REFERENCE_LOAD_UNITS,
        "gene_count": len(gene_metrics),
        "cassette_count": len(cassettes),
        "total_cds_length_nt": total_cds_length_nt,
        "minimum_ostir_reference_count": minimum_count,
        "gene_metrics": gene_metrics,
        "warnings": list(dict.fromkeys(warnings)),
    }
    payload["fingerprint"] = stable_payload_hash(payload)
    return payload


def validate_expression_burden(
    payload: Mapping[str, Any],
    cassettes: list[dict[str, Any]],
    *,
    fallback_promoter_percentile: float | None = None,
) -> None:
    """Validate stored arithmetic and its direct links to cassette inputs."""

    if payload.get("schema_version") != EXPRESSION_BURDEN_SCHEMA_VERSION:
        raise ValueError("expression burden schema is missing or unsupported")
    if payload.get("model_version") != EXPRESSION_BURDEN_MODEL_VERSION:
        raise ValueError("expression burden model version is unsupported")
    recorded_fingerprint = str(payload.get("fingerprint") or "")
    unsigned = dict(payload)
    unsigned.pop("fingerprint", None)
    if recorded_fingerprint != stable_payload_hash(unsigned):
        raise ValueError("expression burden fingerprint is invalid")
    metrics = payload.get("gene_metrics")
    if not isinstance(metrics, list) or not metrics:
        raise ValueError("expression burden gene_metrics are missing")
    expected_genes: list[tuple[int, Mapping[str, Any], Mapping[str, Any]]] = []
    for cassette in cassettes:
        if not isinstance(cassette, Mapping):
            raise ValueError("invalid cassette in expression burden source")
        cassette_index = int(cassette.get("cassette_index") or 0)
        promoter = cassette.get("promoter")
        genes = cassette.get("genes")
        if not isinstance(promoter, Mapping) or not isinstance(genes, list):
            raise ValueError("invalid cassette in expression burden source")
        for gene in genes:
            if not isinstance(gene, Mapping):
                raise ValueError("invalid gene in expression burden source")
            expected_genes.append((cassette_index, promoter, gene))
    if len(metrics) != len(expected_genes):
        raise ValueError("expression burden gene count does not match cassettes")

    recomputed_loads: list[float] = []
    reference_counts: list[int] = []
    total_length = 0
    used_fallback = False
    for metric, (cassette_index, promoter, gene) in zip(
        metrics, expected_genes, strict=True
    ):
        if not isinstance(metric, Mapping):
            raise ValueError("expression burden contains an invalid gene metric")
        accession = str(gene.get("accession") or "")
        if (
            int(metric.get("cassette_index") or 0) != cassette_index
            or str(metric.get("accession") or "") != accession
        ):
            raise ValueError("expression burden gene order or identity is invalid")
        cds_length = int(gene.get("cds_length_nt") or 0)
        ostir = gene.get("ostir")
        rbs = gene.get("rbs")
        if not isinstance(ostir, Mapping) or not isinstance(rbs, Mapping):
            raise ValueError("expression burden source lacks RBS or OSTIR")
        if int(metric.get("cds_length_nt") or 0) != cds_length:
            raise ValueError("expression burden CDS length is invalid")
        if not math.isclose(
            float(metric.get("ostir_translation_initiation_rate") or 0.0),
            float(ostir.get("translation_initiation_rate") or 0.0),
            rel_tol=1e-9,
            abs_tol=1e-8,
        ):
            raise ValueError("expression burden OSTIR rate is invalid")
        unintended = int(ostir.get("unintended_start_count") or 0)
        if int(metric.get("unintended_start_count") or 0) != unintended:
            raise ValueError("expression burden unintended-start count is invalid")
        promoter_source = str(metric.get("promoter_percentile_source") or "")
        measured_promoter = _optional_percentile(
            promoter.get("activity_percentile")
        )
        promoter_percentile = float(
            metric.get("promoter_activity_percentile")
        )
        if promoter_source == "measured":
            if measured_promoter is None or not math.isclose(
                promoter_percentile, measured_promoter, abs_tol=1e-8
            ):
                raise ValueError("expression burden promoter percentile is invalid")
        elif promoter_source != "strategy_target_fallback":
            raise ValueError("expression burden promoter source is invalid")
        else:
            used_fallback = True
            if fallback_promoter_percentile is not None and not math.isclose(
                promoter_percentile,
                float(fallback_promoter_percentile),
                abs_tol=1e-8,
            ):
                raise ValueError("expression burden promoter fallback is invalid")
        rbs_percentile = _optional_percentile(rbs.get("activity_percentile"))
        recorded_rbs = metric.get("rbs_activity_percentile")
        if rbs_percentile is None:
            if recorded_rbs is not None:
                raise ValueError("expression burden RBS fallback is invalid")
        elif recorded_rbs is None or not math.isclose(
            float(recorded_rbs), rbs_percentile, abs_tol=1e-8
        ):
            raise ValueError("expression burden RBS percentile is invalid")
        ostir_percentile = float(metric.get("ostir_expression_percentile"))
        if not 0.0 <= ostir_percentile <= 100.0:
            raise ValueError("expression burden OSTIR percentile is invalid")
        if rbs_percentile is None:
            expected_translation = ostir_percentile
            expected_source = "ostir_only_fallback"
            used_fallback = True
        else:
            expected_translation = (
                RBS_PERCENTILE_WEIGHT * rbs_percentile
                + OSTIR_PERCENTILE_WEIGHT * ostir_percentile
            )
            expected_source = "rbs_and_ostir"
        if str(metric.get("translation_percentile_source")) != expected_source:
            raise ValueError("expression burden translation source is invalid")
        if not math.isclose(
            float(metric.get("combined_translation_percentile")),
            expected_translation,
            abs_tol=1e-7,
        ):
            raise ValueError("expression burden translation percentile is invalid")
        length_factor = _length_factor(cds_length)
        unintended_factor = 1.0 + UNINTENDED_START_WEIGHT * min(
            unintended, UNINTENDED_START_CAP
        )
        load_units = (
            promoter_percentile
            / 100.0
            * expected_translation
            / 100.0
            * length_factor
            * unintended_factor
        )
        for field, expected in (
            ("cds_length_factor", length_factor),
            ("unintended_start_factor", unintended_factor),
            ("load_units", load_units),
        ):
            if not math.isclose(
                float(metric.get(field)), expected, rel_tol=1e-7, abs_tol=1e-7
            ):
                raise ValueError(f"expression burden {field} is invalid")
        recomputed_loads.append(load_units)
        reference_count = int(metric.get("ostir_reference_count") or 0)
        if reference_count < 0:
            raise ValueError("expression burden OSTIR reference count is invalid")
        reference_counts.append(reference_count)
        if reference_count == 0:
            used_fallback = True
        total_length += cds_length
    raw_load = sum(recomputed_loads)
    score = _burden_score(raw_load)
    if not math.isclose(
        float(payload.get("raw_load_units")), raw_load, abs_tol=1e-7
    ):
        raise ValueError("expression burden raw load is invalid")
    if not math.isclose(float(payload.get("score")), score, abs_tol=0.011):
        raise ValueError("expression burden score is invalid")
    if payload.get("level") != burden_level(score):
        raise ValueError("expression burden level is invalid")
    minimum_reference_count = min(reference_counts, default=0)
    if int(payload.get("minimum_ostir_reference_count") or 0) != minimum_reference_count:
        raise ValueError("expression burden minimum reference count is invalid")
    expected_confidence = (
        "high"
        if not used_fallback
        and minimum_reference_count >= MIN_HIGH_CONFIDENCE_REFERENCE_COUNT
        else "medium"
    )
    if payload.get("confidence") != expected_confidence:
        raise ValueError("expression burden confidence is invalid")
    if not math.isclose(
        float(payload.get("reference_load_units")),
        REFERENCE_LOAD_UNITS,
        abs_tol=1e-12,
    ):
        raise ValueError("expression burden reference load is invalid")
    if (
        int(payload.get("gene_count") or 0) != len(expected_genes)
        or int(payload.get("cassette_count") or 0) != len(cassettes)
        or int(payload.get("total_cds_length_nt") or 0) != total_length
    ):
        raise ValueError("expression burden summary counts are invalid")


def expression_burden_summary(payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": payload.get("schema_version"),
        "model_version": payload.get("model_version"),
        "score_type": payload.get("score_type"),
        "score": payload.get("score"),
        "level": payload.get("level"),
        "confidence": payload.get("confidence"),
        "raw_load_units": payload.get("raw_load_units"),
        "reference_load_units": payload.get("reference_load_units"),
        "gene_count": payload.get("gene_count"),
        "cassette_count": payload.get("cassette_count"),
        "total_cds_length_nt": payload.get("total_cds_length_nt"),
        "minimum_ostir_reference_count": payload.get(
            "minimum_ostir_reference_count"
        ),
        "fingerprint": payload.get("fingerprint"),
        "warnings": list(payload.get("warnings") or []),
    }


def validate_expression_burden_summary(payload: Mapping[str, Any]) -> None:
    if payload.get("schema_version") != EXPRESSION_BURDEN_SCHEMA_VERSION:
        raise ValueError("expression burden summary schema is unsupported")
    if payload.get("model_version") != EXPRESSION_BURDEN_MODEL_VERSION:
        raise ValueError("expression burden summary model is unsupported")
    try:
        score = float(payload.get("score"))
        raw_load = float(payload.get("raw_load_units"))
        reference_load = float(payload.get("reference_load_units"))
        gene_count = int(payload.get("gene_count"))
        cassette_count = int(payload.get("cassette_count"))
        total_length = int(payload.get("total_cds_length_nt"))
        reference_count = int(payload.get("minimum_ostir_reference_count"))
    except (TypeError, ValueError) as exc:
        raise ValueError("expression burden summary contains invalid numbers") from exc
    if (
        not 0.0 <= score <= 100.0
        or raw_load < 0.0
        or not math.isclose(reference_load, REFERENCE_LOAD_UNITS, abs_tol=1e-12)
        or gene_count < 1
        or cassette_count < 1
        or total_length < 1
        or reference_count < 0
    ):
        raise ValueError("expression burden summary values are out of range")
    if not math.isclose(score, _burden_score(raw_load), abs_tol=0.011):
        raise ValueError("expression burden summary score is invalid")
    if payload.get("level") != burden_level(score):
        raise ValueError("expression burden summary level is invalid")
    if payload.get("confidence") not in {"high", "medium"}:
        raise ValueError("expression burden summary confidence is invalid")
    fingerprint = str(payload.get("fingerprint") or "")
    if len(fingerprint) != 64 or set(fingerprint) - set("0123456789abcdef"):
        raise ValueError("expression burden summary fingerprint is invalid")


__all__ = [
    "EXPRESSION_BURDEN_MODEL_VERSION",
    "EXPRESSION_BURDEN_SCHEMA_VERSION",
    "calculate_expression_burden",
    "burden_level",
    "burden_model_parameters",
    "empirical_log_percentile",
    "expression_burden_summary",
    "stable_payload_hash",
    "validate_expression_burden",
    "validate_expression_burden_summary",
]
