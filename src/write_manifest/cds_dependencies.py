"""Preserve user-uploaded parts while invalidating CDS-dependent results."""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from typing import Any


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
)


def _fingerprint(content: Any) -> str:
    encoded = json.dumps(content, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


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
        or previous.get("status") != "complete"
        or selection.get("status") != "complete"
        or not isinstance(original_box, Mapping)
        or original_box.get("schema_version") != "expression_box_selection.v1"
        or original_box.get("selection_status") != "user_selected"
    ):
        return {}, CDS_SELECTION_DOWNSTREAM_SECTIONS
    old_proteins = previous.get("proteins", [])
    new_proteins = selection.get("proteins", [])
    old_accessions = {row["accession"] for row in old_proteins}
    current = {row["accession"]: row for row in new_proteins}
    if old_accessions != set(current):
        return {}, CDS_SELECTION_DOWNSTREAM_SECTIONS
    cassettes = original_box.get("cassettes")
    box_fingerprint = original_box.get("selected_design_fingerprint")
    if (
        not isinstance(cassettes, list) or not cassettes
        or not isinstance(box_fingerprint, str) or len(box_fingerprint) != 64
        or any(base not in "0123456789abcdefABCDEF" for base in box_fingerprint)
    ):
        return {}, CDS_SELECTION_DOWNSTREAM_SECTIONS
    assigned = []
    for index, cassette in enumerate(cassettes, 1):
        if (
            not isinstance(cassette, Mapping)
            or cassette.get("cassette_index") != index
            or not isinstance(cassette.get("protein_accessions"), list)
            or not cassette["protein_accessions"]
        ):
            return {}, CDS_SELECTION_DOWNSTREAM_SECTIONS
        assigned.extend(cassette["protein_accessions"])
    if len(assigned) != len(set(assigned)) or set(assigned) != set(current):
        return {}, CDS_SELECTION_DOWNSTREAM_SECTIONS
    old_source = previous.get("source_fingerprint")
    new_source = selection.get("source_fingerprint")
    box_source = original_box.get("source")
    if not old_source or not new_source or not isinstance(box_source, Mapping) or box_source.get("cds_selection_source_fingerprint") != old_source:
        raise ValueError("表达盒分组与当前 CDS 来源不一致，不能保留上传元件")

    box = copy.deepcopy(original_box)
    box["source"]["cds_selection_source_fingerprint"] = new_source
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
    if original_draft is not None:
        if not isinstance(original_draft, Mapping) or original_draft.get("schema_version") != "expression_parts_draft.v1":
            raise ValueError("上传元件草稿格式无效，CDS 更新未提交")
        content = {
            "target_compound_id": original_draft.get("target_compound_id"),
            "source": original_draft.get("source"),
            "cassettes": original_draft.get("cassettes"),
        }
        if original_draft.get("draft_fingerprint") != _fingerprint(content):
            raise ValueError("上传元件草稿 fingerprint 校验失败，CDS 更新未提交")
        source = original_draft.get("source")
        if (
            original_draft.get("target_compound_id") != manifest.get("target_compound_id")
            or not isinstance(source, Mapping)
            or source.get("expression_box_selection_fingerprint") != original_box.get("selected_design_fingerprint")
            or source.get("cds_selection_source_fingerprint") != old_source
        ):
            raise ValueError("上传元件草稿与当前分组或 CDS 来源不一致")
        draft = copy.deepcopy(original_draft)
        draft_cassettes = draft.get("cassettes")
        if not isinstance(draft_cassettes, list) or len(draft_cassettes) != len(cassettes):
            raise ValueError("上传元件草稿与当前表达盒数量不一致")
        for expected, cassette in zip(cassettes, draft_cassettes, strict=True):
            if (
                not isinstance(cassette, Mapping)
                or cassette.get("cassette_index") != expected["cassette_index"]
                or cassette.get("protein_accessions") != expected["protein_accessions"]
            ):
                raise ValueError("上传元件草稿与当前表达盒分组不一致")
            rbs_parts = cassette.get("rbs_by_accession", {})
            if not isinstance(rbs_parts, Mapping) or set(rbs_parts) - set(expected["protein_accessions"]):
                raise ValueError("上传元件草稿包含未知 RBS 蛋白")
            for part in rbs_parts.values():
                if not isinstance(part, dict):
                    raise ValueError("上传 RBS 记录格式无效")
                part.pop("ostir", None)
                part["ostir_status"] = "stale"
        draft["source"]["cds_selection_source_fingerprint"] = new_source
        draft["draft_fingerprint"] = _fingerprint({
            "target_compound_id": draft["target_compound_id"],
            "source": draft["source"], "cassettes": draft_cassettes,
        })
        retained["expression_parts_draft"] = draft
    return retained, tuple(field for field in CDS_SELECTION_DOWNSTREAM_SECTIONS if field not in retained)
