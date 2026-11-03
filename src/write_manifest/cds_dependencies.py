"""Preserve user-uploaded parts while invalidating CDS-dependent results."""

from __future__ import annotations

import copy
from collections.abc import Mapping
from typing import Any

from src.expression_box.parts_preparation import (
    bind_cds, make_parts_draft, normalize_parts_draft, stale_parts_designs,
)


CDS_SELECTION_DOWNSTREAM_SECTIONS = (
    "expression_box_selection",
    "expression_cassette_assembly",
    "expression_parts_draft",
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
    "final_design_report",
)


def _retain_uploaded_draft(manifest: Mapping[str, Any]) -> tuple[dict[str, Any], tuple[str, ...]]:
    raw = manifest.get("expression_parts_draft")
    if not isinstance(raw, Mapping) or not (raw.get("source_type") == "user_uploaded" or raw.get("schema_version") == "expression_parts_draft.v1"):
        return {}, CDS_SELECTION_DOWNSTREAM_SECTIONS
    draft = normalize_parts_draft(raw)
    stale_parts_designs(draft["designs"])
    return {"expression_parts_draft": draft}, tuple(
        name for name in CDS_SELECTION_DOWNSTREAM_SECTIONS if name != "expression_parts_draft"
    )


def cds_dependency_update(
    manifest: Mapping[str, Any], selection: Mapping[str, Any],
) -> tuple[dict[str, Any], tuple[str, ...]]:
    """Return retained sections and obsolete sections for one atomic commit.

    Grouping identity and uploaded files are user choices. Prediction and
    assembly records depend on the edited DNA and must be invalidated.
    """
    previous = manifest.get("cds_selection", {})
    original_box = manifest.get("expression_box_selection")
    if (
        not isinstance(previous, Mapping)
        or previous.get('host') != selection.get('host')
        or previous.get("status") != "complete"
        or selection.get("status") != "complete"
        or not isinstance(original_box, Mapping)
        or original_box.get("schema_version") != "expression_box_selection.v2"
        or original_box.get("selection_status") != "user_selected"
    ):
        return _retain_uploaded_draft(manifest)
    old_proteins = previous.get("proteins", [])
    new_proteins = selection.get("proteins", [])
    old_accessions = {row["accession"] for row in old_proteins}
    current = {row["accession"]: row for row in new_proteins}
    if old_accessions != set(current):
        return _retain_uploaded_draft(manifest)
    cassettes = original_box.get("cassettes")
    if (
        not isinstance(cassettes, list) or not cassettes
    ):
        return _retain_uploaded_draft(manifest)
    assigned = []
    for index, cassette in enumerate(cassettes, 1):
        if (
            not isinstance(cassette, Mapping)
            or cassette.get("cassette_index") != index
            or not isinstance(cassette.get("protein_accessions"), list)
            or not cassette["protein_accessions"]
        ):
            return _retain_uploaded_draft(manifest)
        assigned.extend(cassette["protein_accessions"])
    if len(assigned) != len(set(assigned)) or set(assigned) != set(current):
        return _retain_uploaded_draft(manifest)
    box = copy.deepcopy(original_box)
    for cassette in box["cassettes"]:
        cassette["total_cds_length_nt"] = sum(
            current[accession]["optimized_cds"]["length_nt"]
            for accession in cassette["protein_accessions"]
        )
    box.setdefault("summary", {}).update(
        cassette_count=len(cassettes), protein_count=len(current),
        total_cds_length_nt=sum(cassette["total_cds_length_nt"] for cassette in box["cassettes"]),
    )
    retained = {"expression_box_selection": box}
    original_draft = manifest.get("expression_parts_draft")
    if original_draft is not None and not (
        isinstance(original_draft, Mapping) and original_draft.get("source_type") == "recommended"
    ):
        if not isinstance(original_draft, Mapping):
            raise ValueError("上传元件草稿格式无效，CDS 更新未提交")
        draft = normalize_parts_draft(original_draft)
        source = draft["source"]
        if (
            draft["target_compound_id"] != manifest.get("target_compound_id")
        ):
            raise ValueError("元件准备记录与当前分组或 CDS 来源不一致")
        for design in draft["designs"]:
            draft_cassettes = design["cassettes"]
            if len(draft_cassettes) != len(cassettes):
                raise ValueError("元件准备记录与当前表达盒数量不一致")
            for expected, cassette in zip(cassettes, draft_cassettes, strict=True):
                if (
                    cassette.get("cassette_index") != expected["cassette_index"]
                    or cassette.get("protein_accessions") != expected["protein_accessions"]
                ):
                    raise ValueError("元件准备记录与当前表达盒分组不一致")
            bind_cds(draft_cassettes, current)
        stale_parts_designs(draft["designs"])
        draft = make_parts_draft(
            target=draft["target_compound_id"], source=draft["source"],
            source_type=draft["source_type"], designs=draft["designs"],
        )
        retained["expression_parts_draft"] = draft
    return retained, tuple(field for field in CDS_SELECTION_DOWNSTREAM_SECTIONS if field not in retained)
