"""Restriction-ligation plan candidates for one complete expression construct."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from src.final_assemble_plan.common import (
    enzyme_catalog,
    enzyme_cut_positions,
    enzyme_site_start,
    overlaps_spans,
    protected_coordinate_spans,
    region_bounds,
    span_inside_region,
)
from src.final_assemble_plan.models import (
    AssemblyConstruct,
    FinalAssemblyContext,
)


def _size_score(length_bp: int) -> float:
    if length_bp <= 12_000:
        return 5.0
    if length_bp <= 15_000:
        return 3.0
    return 1.0


def _valid_cutters(
    context: FinalAssemblyContext,
    construct: AssemblyConstruct,
) -> list[dict[str, Any]]:
    backbone = context.backbone
    protected_spans = protected_coordinate_spans(
        backbone.protected_features,
        backbone.genbank_features,
    )
    cutters: list[dict[str, Any]] = []
    for item in enzyme_catalog():
        enzyme = item["enzyme"]
        vector_cuts = enzyme_cut_positions(
            enzyme, backbone.sequence, circular=True
        )
        insert_cuts = enzyme_cut_positions(
            enzyme, construct.sequence, circular=False
        )
        if len(vector_cuts) != 1 or insert_cuts:
            continue
        cut_after = vector_cuts[0]
        site_start = enzyme_site_start(
            enzyme, cut_after, backbone.length_bp
        )
        site_end = site_start + int(item["site_length"]) - 1
        if site_end > backbone.length_bp:
            continue
        matching_regions = [
            dict(region)
            for region in backbone.insertion_regions
            if span_inside_region(site_start, site_end, region)
        ]
        if not matching_regions or overlaps_spans(
            site_start, site_end, protected_spans
        ):
            continue
        cutters.append(
            {
                "name": item["name"],
                "recognition_site": item["site"],
                "site_start_bp": site_start,
                "site_end_bp": site_end,
                "cut_after_bp": cut_after,
                "overhang": item["overhang"],
                "priority": item["priority"],
                "insertion_region": matching_regions[0],
            }
        )
    cutters.sort(
        key=lambda item: (
            int(item["site_start_bp"]),
            -int(item["priority"]),
            str(item["name"]),
        )
    )
    return cutters


def _same_region(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    return left.get("insertion_region") == right.get("insertion_region")


def _pair_candidate(
    context: FinalAssemblyContext,
    construct: AssemblyConstruct,
    left: Mapping[str, Any],
    right: Mapping[str, Any],
) -> dict[str, Any] | None:
    if left.get("name") == right.get("name") or not _same_region(left, right):
        return None
    start = int(left["site_start_bp"])
    end = int(right["site_end_bp"])
    if start >= int(right["site_start_bp"]):
        return None
    region = left["insertion_region"]
    if not span_inside_region(start, end, region):
        return None
    protected_spans = protected_coordinate_spans(
        context.backbone.protected_features,
        context.backbone.genbank_features,
    )
    if overlaps_spans(start, end, protected_spans):
        return None
    replaced_length = end - start + 1
    final_length = (
        context.backbone.length_bp - replaced_length + construct.length_bp
        + len(str(left["recognition_site"]))
        + len(str(right["recognition_site"]))
    )
    five_prime_count = sum(
        item.get("overhang") == "five_prime" for item in (left, right)
    )
    junction_quality = 25.0 if five_prime_count == 2 else 22.0 if five_prime_count == 1 else 19.0
    breakdown = {
        "insertion_region_safety": 40.0,
        "junction_quality": junction_quality,
        "experimental_simplicity": 20.0,
        "source_integrity": 10.0,
        "estimated_final_size": _size_score(final_length),
    }
    score = round(sum(breakdown.values()), 2)
    enzymes = [
        {
            "role": "left",
            **{key: left[key] for key in (
                "name",
                "recognition_site",
                "site_start_bp",
                "site_end_bp",
                "cut_after_bp",
                "overhang",
            )},
        },
        {
            "role": "right",
            **{key: right[key] for key in (
                "name",
                "recognition_site",
                "site_start_bp",
                "site_end_bp",
                "cut_after_bp",
                "overhang",
            )},
        },
    ]
    return {
        "parts_design_id": construct.design_id,
        "assembly_method": "restriction",
        "target": {
            "mode": "replace",
            "replace_start_bp": start,
            "replace_end_bp": end,
            "insertion_region": region,
        },
        "backbone_linearization": {
            "mode": "restriction",
            "restriction_enzymes": enzymes,
            "enzyme_summary": f"{left['name']}/{right['name']}",
        },
        "restriction": {
            "left_enzyme": left["name"],
            "right_enzyme": right["name"],
            "left_site": left["recognition_site"],
            "right_site": right["recognition_site"],
            "restriction_site_retention": "retain",
        },
        "gibson": None,
        "estimated_final_length_bp": final_length,
        "score": score,
        "score_breakdown": breakdown,
        "rationales": [
            "Both enzymes cut the circular backbone exactly once.",
            "Neither enzyme cuts this complete expression construct.",
            "Both recognition sites are inside the same audited insertion region.",
            "The replacement span does not overlap mapped protected features.",
        ],
        "warnings": [
            "Double-digest buffer and incubation compatibility require experimental confirmation."
        ],
        "_rank_key": (
            -score,
            -int(left.get("priority") or 0) - int(right.get("priority") or 0),
            end - start,
            str(left["name"]),
            str(right["name"]),
        ),
    }


def _seva_candidate(
    context: FinalAssemblyContext,
    construct: AssemblyConstruct,
) -> dict[str, Any] | None:
    regions = [
        region
        for region in context.backbone.insertion_regions
        if str(region.get("type") or "") == "seva_cargo"
    ]
    if len(regions) != 1:
        return None
    cutters = _valid_cutters(context, construct)
    left = next((item for item in cutters if item["name"] == "PacI"), None)
    right = next((item for item in cutters if item["name"] == "SpeI"), None)
    if left is None or right is None:
        return None
    candidate = _pair_candidate(context, construct, left, right)
    if candidate is None:
        return None
    bounds = region_bounds(regions[0])
    if bounds is None:
        return None
    if (
        candidate["target"]["replace_start_bp"] != bounds[0]
        or candidate["target"]["replace_end_bp"] != bounds[1]
    ):
        return None
    candidate["rationales"].insert(
        0,
        "The selected backbone policy requires audited PacI/SpeI cargo replacement.",
    )
    return candidate


def recommend_restriction_plans(
    context: FinalAssemblyContext,
    construct: AssemblyConstruct,
    *,
    limit: int = 20,
) -> list[dict[str, Any]]:
    """Return safe double-digest candidates for one construct."""

    if context.backbone.assembly_policy == "replace_seva_cargo_paci_spei":
        candidate = _seva_candidate(context, construct)
        return [candidate] if candidate is not None else []
    cutters = _valid_cutters(context, construct)
    candidates: list[dict[str, Any]] = []
    for index, left in enumerate(cutters):
        for right in cutters[index + 1 :]:
            candidate = _pair_candidate(context, construct, left, right)
            if candidate is not None:
                candidates.append(candidate)
    candidates.sort(key=lambda item: item["_rank_key"])
    for item in candidates:
        item.pop("_rank_key", None)
    return candidates[:limit]


__all__ = ["recommend_restriction_plans"]
