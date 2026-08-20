"""Gibson plan candidates for one complete expression construct."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from src.final_assemble_plan.common import (
    arm_quality,
    circular_segment,
    enzyme_catalog,
    enzyme_cut_positions,
    enzyme_site_start,
    overlaps_spans,
    protected_coordinate_spans,
    region_bounds,
    span_inside_region,
)
from src.final_assemble_plan.config import DEFAULT_HOMOLOGY_ARM_LENGTH
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


def _target_arms(
    sequence: str,
    target: Mapping[str, Any],
    arm_length: int,
) -> tuple[str, str]:
    if target.get("mode") == "insert_after":
        insert_after = int(target["insert_after_bp"])
        left_start = insert_after - arm_length + 1
        right_start = insert_after + 1
    else:
        start = int(target["replace_start_bp"])
        end = int(target["replace_end_bp"])
        left_start = start - arm_length
        right_start = end + 1
    return (
        circular_segment(sequence, left_start, arm_length),
        circular_segment(sequence, right_start, arm_length),
    )


def _linearization_enzymes(
    context: FinalAssemblyContext,
) -> list[dict[str, Any]]:
    backbone = context.backbone
    protected_spans = protected_coordinate_spans(
        backbone.protected_features,
        backbone.genbank_features,
    )
    result: list[dict[str, Any]] = []
    for item in enzyme_catalog():
        enzyme = item["enzyme"]
        cuts = enzyme_cut_positions(enzyme, backbone.sequence, circular=True)
        if len(cuts) != 1:
            continue
        cut_after = cuts[0]
        site_start = enzyme_site_start(enzyme, cut_after, backbone.length_bp)
        site_end = site_start + int(item["site_length"]) - 1
        if site_end > backbone.length_bp:
            continue
        regions = [
            dict(region)
            for region in backbone.insertion_regions
            if span_inside_region(site_start, site_end, region)
        ]
        if not regions or overlaps_spans(site_start, site_end, protected_spans):
            continue
        result.append(
            {
                "name": item["name"],
                "recognition_site": item["site"],
                "site_start_bp": site_start,
                "site_end_bp": site_end,
                "cut_after_bp": cut_after,
                "overhang": item["overhang"],
                "priority": item["priority"],
                "insertion_region": regions[0],
            }
        )
    result.sort(
        key=lambda item: (
            -int(item["priority"]),
            int(item["cut_after_bp"]),
            str(item["name"]),
        )
    )
    return result


def _candidate(
    context: FinalAssemblyContext,
    construct: AssemblyConstruct,
    *,
    target: dict[str, Any],
    linearization_mode: str,
    enzymes: list[dict[str, Any]],
    arm_length: int,
) -> dict[str, Any] | None:
    left_sequence, right_sequence = _target_arms(
        context.backbone.sequence,
        target,
        arm_length,
    )
    left_arm = arm_quality(left_sequence, context.backbone.sequence)
    right_arm = arm_quality(right_sequence, context.backbone.sequence)
    if not left_arm["passes"] or not right_arm["passes"]:
        return None
    if target["mode"] == "insert_after":
        final_length = context.backbone.length_bp + construct.length_bp
    else:
        replaced_length = (
            int(target["replace_end_bp"])
            - int(target["replace_start_bp"])
            + 1
        )
        final_length = (
            context.backbone.length_bp - replaced_length + construct.length_bp
        )
    average_arm_quality = (
        float(left_arm["quality"]) + float(right_arm["quality"])
    ) / 2.0
    junction_quality = round(20.0 + 5.0 * average_arm_quality, 2)
    simplicity = 16.0 if linearization_mode == "pcr" else 18.0 if len(enzymes) == 1 else 17.0
    breakdown = {
        "insertion_region_safety": 40.0,
        "junction_quality": junction_quality,
        "experimental_simplicity": simplicity,
        "source_integrity": 10.0,
        "estimated_final_size": _size_score(final_length),
    }
    score = round(sum(breakdown.values()), 2)
    if linearization_mode == "pcr":
        enzyme_summary = "none (PCR linearization)"
        normalized_enzymes: list[dict[str, Any]] = []
    else:
        normalized_enzymes = [
            {
                key: enzyme[key]
                for key in (
                    "name",
                    "recognition_site",
                    "site_start_bp",
                    "site_end_bp",
                    "cut_after_bp",
                    "overhang",
                )
            }
            for enzyme in enzymes
        ]
        enzyme_summary = "/".join(item["name"] for item in enzymes)
    return {
        "parts_design_id": construct.design_id,
        "assembly_method": "gibson",
        "target": target,
        "backbone_linearization": {
            "mode": linearization_mode,
            "restriction_enzymes": normalized_enzymes,
            "enzyme_summary": enzyme_summary,
        },
        "restriction": None,
        "gibson": {
            "homology_arm_length": arm_length,
            "left_homology": left_sequence,
            "right_homology": right_sequence,
            "left_arm_audit": left_arm,
            "right_arm_audit": right_arm,
        },
        "estimated_final_length_bp": final_length,
        "score": score,
        "score_breakdown": breakdown,
        "rationales": [
            "The target is inside an audited insertion region and outside mapped protected features.",
            "Both homology arms are unique in the circular backbone.",
            "Both homology arms pass GC and homopolymer gates.",
            (
                "The backbone is linearized by PCR, so no restriction enzyme is required."
                if linearization_mode == "pcr"
                else f"The backbone is linearized with {enzyme_summary}."
            ),
        ],
        "warnings": [
            "Primer annealing segments and reaction conditions require experimental confirmation."
        ],
        "_rank_key": (
            -score,
            0 if linearization_mode == "restriction" else 1,
            enzyme_summary,
            int(target.get("insert_after_bp") or target.get("replace_start_bp") or 0),
        ),
    }


def _replacement_target(
    context: FinalAssemblyContext,
) -> dict[str, Any] | None:
    regions = [
        dict(region)
        for region in context.backbone.insertion_regions
        if str(region.get("type") or "") == "seva_cargo"
    ]
    if len(regions) != 1:
        return None
    bounds = region_bounds(regions[0])
    if bounds is None:
        return None
    protected = protected_coordinate_spans(
        context.backbone.protected_features,
        context.backbone.genbank_features,
    )
    if overlaps_spans(bounds[0], bounds[1], protected):
        return None
    return {
        "mode": "replace",
        "replace_start_bp": bounds[0],
        "replace_end_bp": bounds[1],
        "insertion_region": regions[0],
    }


def _pcr_insertion_targets(
    context: FinalAssemblyContext,
) -> list[dict[str, Any]]:
    protected = protected_coordinate_spans(
        context.backbone.protected_features,
        context.backbone.genbank_features,
    )
    targets: list[dict[str, Any]] = []
    for region in context.backbone.insertion_regions:
        bounds = region_bounds(region)
        if bounds is None:
            continue
        for cut_after in range(bounds[0], bounds[1] + 1):
            if overlaps_spans(cut_after, cut_after, protected):
                continue
            targets.append(
                {
                    "mode": "insert_after",
                    "insert_after_bp": cut_after,
                    "insertion_region": dict(region),
                }
            )
    return targets


def _seva_enzyme_linearization(
    context: FinalAssemblyContext,
) -> list[dict[str, Any]]:
    cutters = _linearization_enzymes(context)
    paci = next((item for item in cutters if item["name"] == "PacI"), None)
    spei = next((item for item in cutters if item["name"] == "SpeI"), None)
    return [paci, spei] if paci is not None and spei is not None else []


def recommend_gibson_plans(
    context: FinalAssemblyContext,
    construct: AssemblyConstruct,
    *,
    arm_length: int = DEFAULT_HOMOLOGY_ARM_LENGTH,
    limit: int = 20,
) -> list[dict[str, Any]]:
    """Return safe Gibson candidates, including PCR/enzyme linearization."""

    candidates: list[dict[str, Any]] = []
    if context.backbone.assembly_policy == "replace_seva_cargo_paci_spei":
        target = _replacement_target(context)
        if target is None:
            return []
        pcr_candidate = _candidate(
            context,
            construct,
            target=target,
            linearization_mode="pcr",
            enzymes=[],
            arm_length=arm_length,
        )
        if pcr_candidate is not None:
            candidates.append(pcr_candidate)
        enzymes = _seva_enzyme_linearization(context)
        if enzymes:
            enzyme_candidate = _candidate(
                context,
                construct,
                target=target,
                linearization_mode="restriction",
                enzymes=enzymes,
                arm_length=arm_length,
            )
            if enzyme_candidate is not None:
                candidates.append(enzyme_candidate)
    else:
        for target in _pcr_insertion_targets(context):
            candidate = _candidate(
                context,
                construct,
                target=target,
                linearization_mode="pcr",
                enzymes=[],
                arm_length=arm_length,
            )
            if candidate is not None:
                candidates.append(candidate)
        for enzyme in _linearization_enzymes(context):
            target = {
                "mode": "insert_after",
                "insert_after_bp": enzyme["cut_after_bp"],
                "insertion_region": enzyme["insertion_region"],
            }
            candidate = _candidate(
                context,
                construct,
                target=target,
                linearization_mode="restriction",
                enzymes=[enzyme],
                arm_length=arm_length,
            )
            if candidate is not None:
                candidates.append(candidate)
    candidates.sort(key=lambda item: item["_rank_key"])
    unique: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    for item in candidates:
        item.pop("_rank_key", None)
        target = item["target"]
        key = (
            target["mode"],
            target.get("insert_after_bp"),
            target.get("replace_start_bp"),
            target.get("replace_end_bp"),
            item["backbone_linearization"]["mode"],
            item["backbone_linearization"]["enzyme_summary"],
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(item)
        if len(unique) >= limit:
            break
    return unique


__all__ = ["recommend_gibson_plans"]
