"""Pure hard-filtering, scoring, and diversity logic for plasmid templates."""

from __future__ import annotations

import re
import statistics
from collections.abc import Mapping, Sequence
from typing import Any

from src.plasmid_selection.config import (
    COPY_CLASS_SCORES,
    MARKER_ALIASES,
    MARKER_SCORES,
    SUPPORTED_ASSEMBLY_POLICIES,
    SUPPORTED_PRIORITIES,
)
from src.plasmid_selection.models import PlasmidContext


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def normalize_marker_name(value: Any) -> str:
    text = str(value or "").strip().lower()
    collapsed = re.sub(r"[^a-z0-9/]+", "", text)
    if text in MARKER_ALIASES:
        return MARKER_ALIASES[text]
    if collapsed in MARKER_ALIASES:
        return MARKER_ALIASES[collapsed]
    for alias, canonical in sorted(
        MARKER_ALIASES.items(), key=lambda item: len(item[0]), reverse=True
    ):
        if len(alias) > 2 and alias in text:
            return canonical
    return ""


def template_marker(template: Mapping[str, Any]) -> str:
    raw_markers = template.get("selection_markers")
    if isinstance(raw_markers, list):
        for raw in raw_markers:
            if not isinstance(raw, Mapping):
                continue
            for key in ("family", "phenotype", "gene"):
                normalized = normalize_marker_name(raw.get(key))
                if normalized:
                    return normalized
    for value in (
        template.get("bacterial_resistance"),
        template.get("resistance_markers"),
    ):
        normalized = normalize_marker_name(value)
        if normalized:
            return normalized
        if isinstance(value, list):
            for item in value:
                if isinstance(item, Mapping):
                    for nested in item.values():
                        normalized = normalize_marker_name(nested)
                        if normalized:
                            return normalized
                else:
                    normalized = normalize_marker_name(item)
                    if normalized:
                        return normalized
    return ""


def normalize_marker_preferences(
    preferred: str | None,
    excluded: Sequence[str] | None,
) -> tuple[str | None, tuple[str, ...]]:
    preferred_marker: str | None = None
    if preferred:
        preferred_marker = normalize_marker_name(preferred)
        if not preferred_marker:
            raise ValueError(f"不支持的优先抗性标记：{preferred}")
    normalized_excluded: list[str] = []
    for value in excluded or ():
        marker = normalize_marker_name(value)
        if not marker:
            raise ValueError(f"不支持的排除抗性标记：{value}")
        if marker not in normalized_excluded:
            normalized_excluded.append(marker)
    if preferred_marker in normalized_excluded:
        raise ValueError("同一抗性标记不能同时设为优先和排除")
    return preferred_marker, tuple(normalized_excluded)


def _is_nonempty_list(value: Any) -> bool:
    return isinstance(value, list) and bool(value)


def _has_valid_insertion_region(value: Any) -> bool:
    if not isinstance(value, list) or not value:
        return False
    for region in value:
        if not isinstance(region, Mapping):
            continue
        try:
            start = int(region.get("start_bp"))
            end = int(region.get("end_bp"))
        except (TypeError, ValueError):
            continue
        if start >= 0 and end >= start:
            return True
    return False


def hard_filter_reason(
    template: Mapping[str, Any],
    *,
    excluded_markers: Sequence[str] = (),
) -> str | None:
    if template.get("schema_version") != "plasmid_template.v2":
        return "unsupported_schema"
    if template.get("audit_status") != "PASS" or template.get("audit_passed") is not True:
        return "sequence_audit_not_passed"
    if template.get("mg1655_compatible") is not True:
        return "not_mg1655_compatible"
    if str(template.get("topology") or "").strip().lower() != "circular":
        return "not_circular"
    if not str(template.get("plasmid_id") or "").strip():
        return "missing_plasmid_id"
    if not str(template.get("name") or "").strip():
        return "missing_name"
    try:
        length_bp = int(template.get("length_bp") or 0)
    except (TypeError, ValueError):
        length_bp = 0
    if length_bp < 1:
        return "invalid_backbone_length"
    content_hash = str(template.get("sequence_content_sha256") or "").lower()
    canonical_hash = str(template.get("canonical_sequence_sha256") or "").lower()
    source_file_hash = str(template.get("source_file_sha256") or "").lower()
    if (
        not _SHA256_RE.fullmatch(content_hash)
        or not _SHA256_RE.fullmatch(canonical_hash)
        or not _SHA256_RE.fullmatch(source_file_hash)
    ):
        return "missing_sequence_hashes"
    if not str(template.get("source") or "").strip():
        return "missing_sequence_source"
    if not str(template.get("sequence_file") or "").strip():
        return "missing_sequence_source"
    if not (
        str(template.get("source_record_id") or "").strip()
        or str(template.get("sequence_url") or "").strip()
    ):
        return "missing_sequence_source"
    if not isinstance(template.get("source_provenance"), Mapping) or not template.get("source_provenance"):
        return "missing_source_provenance"
    if not _is_nonempty_list(template.get("evidence_refs")):
        return "missing_evidence"
    if not _is_nonempty_list(template.get("origins")):
        return "missing_replication_origin"
    marker = template_marker(template)
    if not marker:
        return "missing_or_unsupported_marker"
    if marker in set(excluded_markers):
        return f"excluded_marker:{marker}"
    if not _is_nonempty_list(template.get("host_compatibility")):
        return "missing_host_compatibility"
    if not _has_valid_insertion_region(template.get("insertion_regions")):
        return "missing_audited_insertion_region"
    if not _is_nonempty_list(template.get("protected_features")):
        return "missing_protected_features"
    dependencies = template.get("replication_dependencies")
    if isinstance(dependencies, list) and dependencies:
        return "unsupported_replication_dependencies"
    if dependencies not in (None, [], ""):
        return "unsupported_replication_dependencies"
    if template.get("has_expression_cargo") is not False:
        return "expression_cargo_present_or_unknown"
    policy = str(template.get("assembly_policy") or "").strip()
    if policy not in SUPPORTED_ASSEMBLY_POLICIES:
        return f"unsupported_assembly_policy:{policy or 'missing'}"
    copy_class = str(template.get("copy_number_class") or "").strip().lower()
    if copy_class not in {"low", "medium", "context_dependent", "high"}:
        return "unsupported_copy_number_class"
    if not str(template.get("replicon_family") or "").strip():
        return "missing_replicon_family"
    return None


def insertion_region_length(template: Mapping[str, Any]) -> int:
    if template.get("assembly_policy") != "replace_seva_cargo_paci_spei":
        return 0
    lengths: list[int] = []
    for region in template.get("insertion_regions") or []:
        if not isinstance(region, Mapping):
            continue
        try:
            start = int(region.get("start_bp"))
            end = int(region.get("end_bp"))
        except (TypeError, ValueError):
            continue
        if start >= 0 and end >= start:
            lengths.append(end - start + 1)
    return max(lengths, default=0)


def estimated_final_length(template: Mapping[str, Any], insert_length_bp: int) -> int:
    backbone_length = int(template.get("length_bp") or 0)
    if template.get("assembly_policy") == "replace_seva_cargo_paci_spei":
        return backbone_length - insertion_region_length(template) + insert_length_bp
    return backbone_length + insert_length_bp


def _size_score(length_bp: int) -> float:
    if length_bp <= 8_000:
        return 10.0
    if length_bp <= 10_000:
        return 8.0
    if length_bp <= 12_000:
        return 5.0
    return 2.0


def _copy_load_score(copy_class: str, priority: str, insert_length_bp: int) -> float:
    score = COPY_CLASS_SCORES[priority][copy_class]
    if insert_length_bp >= 15_000:
        score -= 10.0
    elif insert_length_bp >= 10_000:
        score -= 5.0
    return max(0.0, score)


def _marker_score(marker: str, preferred_marker: str | None) -> float:
    base = MARKER_SCORES[marker]
    if preferred_marker is None:
        return base
    if marker == preferred_marker:
        return 10.0
    return min(base, 4.0)


def score_template(
    template: Mapping[str, Any],
    context: PlasmidContext,
    *,
    priority: str = "stability",
    preferred_marker: str | None = None,
) -> dict[str, Any]:
    if priority not in SUPPORTED_PRIORITIES:
        raise ValueError(f"不支持的质粒优先策略：{priority}")
    marker = template_marker(template)
    copy_class = str(template.get("copy_number_class") or "").lower()
    policy = str(template.get("assembly_policy") or "")
    assembly_score = 25.0 if policy == "insert_into_mcs" else 23.0
    evidence_score = 20.0
    marker_score = _marker_score(marker, preferred_marker)

    pair_scores: list[dict[str, Any]] = []
    for construct in context.constructs:
        final_length = estimated_final_length(template, construct.length_bp)
        breakdown = {
            "copy_load_fit": _copy_load_score(
                copy_class, priority, construct.length_bp
            ),
            "assembly_readiness": assembly_score,
            "source_evidence_completeness": evidence_score,
            "marker_suitability": marker_score,
            "estimated_final_size": _size_score(final_length),
        }
        total = round(sum(breakdown.values()), 2)
        pair_scores.append(
            {
                "parts_design_id": construct.design_id,
                "insert_length_bp": construct.length_bp,
                "estimated_final_length_bp": final_length,
                "score": total,
                "breakdown": breakdown,
            }
        )
    robust_score = min(item["score"] for item in pair_scores)
    median_score = round(statistics.median(item["score"] for item in pair_scores), 2)
    maximum_final_length = max(
        item["estimated_final_length_bp"] for item in pair_scores
    )
    minimum_final_length = min(
        item["estimated_final_length_bp"] for item in pair_scores
    )
    representative = min(
        pair_scores,
        key=lambda item: (item["score"], -item["estimated_final_length_bp"]),
    )["breakdown"]
    rationales = [
        f"{copy_class} copy backbone under the {priority} priority",
        (
            "audited MCS insertion region"
            if policy == "insert_into_mcs"
            else "audited PacI/SpeI SEVA cargo-replacement region"
        ),
        f"{marker} selection marker",
        f"estimated final plasmid size {minimum_final_length}-{maximum_final_length} bp",
        f"one backbone is applicable to all {len(pair_scores)} selected expression designs",
    ]
    warnings = [
        "Backbone capacity is estimated from copy class and final size because the database has no experimentally verified insert-capacity field.",
        "Recommendation score is an interpretable heuristic, not an experimental success probability.",
    ]
    return {
        "plasmid_id": template.get("plasmid_id"),
        "name": template.get("name"),
        "replicon_family": template.get("replicon_family"),
        "copy_number_class": copy_class,
        "marker": marker,
        "assembly_policy": policy,
        "robust_score": robust_score,
        "pair_score_median": median_score,
        "score_breakdown": representative,
        "estimated_final_length_range_bp": {
            "minimum": minimum_final_length,
            "maximum": maximum_final_length,
        },
        "applicable_parts_design_ids": [item.design_id for item in context.constructs],
        "pair_scores": pair_scores,
        "confidence": "medium",
        "rationales": rationales,
        "warnings": warnings,
        "template": dict(template),
    }


def rank_templates(
    templates: Sequence[Mapping[str, Any]],
    context: PlasmidContext,
    *,
    priority: str,
    preferred_marker: str | None,
    excluded_markers: Sequence[str],
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    scored: list[dict[str, Any]] = []
    rejected: list[dict[str, str]] = []
    for template in templates:
        reason = hard_filter_reason(template, excluded_markers=excluded_markers)
        if reason is not None:
            rejected.append(
                {
                    "plasmid_id": str(template.get("plasmid_id") or ""),
                    "name": str(template.get("name") or ""),
                    "stage": "hard_filter",
                    "reason": reason,
                }
            )
            continue
        scored.append(
            score_template(
                template,
                context,
                priority=priority,
                preferred_marker=preferred_marker,
            )
        )
    scored.sort(
        key=lambda item: (
            -float(item["robust_score"]),
            -float(item["pair_score_median"]),
            int(item["estimated_final_length_range_bp"]["maximum"]),
            str(item["plasmid_id"]),
        )
    )
    rejected.sort(key=lambda item: (item["stage"], item["plasmid_id"]))
    return scored, rejected


def select_diverse_candidates(
    ranked: Sequence[dict[str, Any]],
    count: int,
) -> list[dict[str, Any]]:
    """Keep the raw best candidate, then greedily add distinct alternatives."""

    if count < 1 or not ranked:
        return []
    selected: list[dict[str, Any]] = [ranked[0]]
    selected_ids = {str(ranked[0]["plasmid_id"])}
    tuple_keys = {
        (
            str(ranked[0]["replicon_family"]).lower(),
            str(ranked[0]["marker"]).lower(),
            str(ranked[0]["assembly_policy"]).lower(),
        )
    }
    copy_counts = {str(ranked[0]["copy_number_class"]): 1}

    def add_with_rules(*, enforce_copy_cap: bool, enforce_tuple: bool) -> None:
        for item in ranked:
            if len(selected) >= count:
                return
            identity = str(item["plasmid_id"])
            if identity in selected_ids:
                continue
            copy_class = str(item["copy_number_class"])
            key = (
                str(item["replicon_family"]).lower(),
                str(item["marker"]).lower(),
                str(item["assembly_policy"]).lower(),
            )
            if enforce_copy_cap and copy_counts.get(copy_class, 0) >= 2:
                continue
            if enforce_tuple and key in tuple_keys:
                continue
            selected.append(item)
            selected_ids.add(identity)
            tuple_keys.add(key)
            copy_counts[copy_class] = copy_counts.get(copy_class, 0) + 1

    add_with_rules(enforce_copy_cap=True, enforce_tuple=True)
    add_with_rules(enforce_copy_cap=False, enforce_tuple=True)
    add_with_rules(enforce_copy_cap=False, enforce_tuple=False)
    return selected[:count]


__all__ = [
    "estimated_final_length",
    "hard_filter_reason",
    "insertion_region_length",
    "normalize_marker_name",
    "normalize_marker_preferences",
    "rank_templates",
    "score_template",
    "select_diverse_candidates",
    "template_marker",
]
