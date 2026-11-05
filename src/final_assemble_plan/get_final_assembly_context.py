"""Manifest adapter for final-assembly planning."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.final_assemble_plan.common import (
    feature_payloads,
    read_genbank,
    resolve_project_file,
)
from src.final_assemble_plan.models import (
    AssemblyBackbone,
    AssemblyConstruct,
    FinalAssemblyContext,
)
from src.plasmid_selection.get_plasmid_context import load_plasmid_context
from src.write_manifest.store import read_design_manifest


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少 {name}")
    return value


def load_final_assembly_context(config: Any) -> FinalAssemblyContext:
    """Load one selected backbone and every selected complete construct."""

    plasmid_context = load_plasmid_context(config)
    manifest = read_design_manifest(plasmid_context.manifest_path)

    return build_final_assembly_context(plasmid_context, manifest.get("plasmid_selection"))


def build_final_assembly_context(plasmid_context: Any, selection: Any) -> FinalAssemblyContext:
    """Validate an explicit selection before its atomic manifest registration."""
    selection = _mapping(selection, "plasmid_selection")
    if selection.get("schema_version") != "plasmid_selection.v3":
        raise ValueError(
            "manifest 中没有新版质粒选择，请重新运行 plasmid --recommend "
            "并使用 write --plasmid N"
        )
    if selection.get("status") != "selected":
        raise ValueError("plasmid_selection.status 必须为 selected")

    vector = _mapping(selection.get("vector"), "plasmid_selection.vector")
    sequence_file = _mapping(
        vector.get("selected_sequence_file"),
        "plasmid_selection.vector.selected_sequence_file",
    )
    backbone_path = resolve_project_file(
        plasmid_context.project_output_path,
        sequence_file.get("path"),
    )
    record, sequence, _ = read_genbank(backbone_path)
    try:
        recorded_length = int(sequence_file.get("length_bp") or 0)
        vector_length = int(vector.get("length_bp") or recorded_length)
    except (TypeError, ValueError) as exc:
        raise ValueError("selected backbone length_bp 无效") from exc
    if recorded_length != len(sequence) or vector_length != len(sequence):
        raise ValueError("selected backbone 长度与 manifest 不一致")
    topology = str(
        vector.get("topology") or record.annotations.get("topology") or ""
    ).strip().lower()
    if topology != "circular":
        raise ValueError("最终组装仅支持 circular 质粒骨架")

    assembly_policy = str(vector.get("assembly_policy") or "").strip()
    if assembly_policy not in {
        "insert_into_mcs",
        "replace_seva_cargo_paci_spei",
    }:
        raise ValueError(f"不支持的 backbone assembly_policy：{assembly_policy}")
    raw_regions = vector.get("insertion_regions")
    regions = tuple(
        dict(item)
        for item in raw_regions
        if isinstance(item, Mapping)
    ) if isinstance(raw_regions, list) else ()
    if not regions:
        raise ValueError("selected backbone 缺少 audited insertion_regions")
    raw_protected = vector.get("protected_features")
    protected = tuple(
        dict(item)
        for item in raw_protected
        if isinstance(item, Mapping)
    ) if isinstance(raw_protected, list) else ()
    if not protected:
        raise ValueError("selected backbone 缺少 protected_features")

    backbone = AssemblyBackbone(
        plasmid_id=str(vector.get("plasmid_id") or ""),
        name=str(vector.get("name") or ""),
        path=backbone_path,
        length_bp=len(sequence),
        sequence=sequence,
        topology=topology,
        assembly_policy=assembly_policy,
        insertion_regions=regions,
        protected_features=protected,
        genbank_features=feature_payloads(record),
    )
    constructs = tuple(
        AssemblyConstruct(
            design_id=item.design_id,
            rank=item.rank,
            path=item.path,
            length_bp=item.length_bp,
            sequence=read_genbank(item.path)[1],
        )
        for item in plasmid_context.constructs
    )
    return FinalAssemblyContext(
        target_compound_id=plasmid_context.target_compound_id,
        manifest_path=plasmid_context.manifest_path,
        project_output_path=plasmid_context.project_output_path,
        constructs=constructs,
        backbone=backbone,
    )


def final_assembly_context_summary(config: Any) -> dict[str, Any]:
    context = load_final_assembly_context(config)
    return {
        "target_compound_id": context.target_compound_id,
        "construct_count": len(context.constructs),
        "parts_design_ids": [item.design_id for item in context.constructs],
        "backbone": {
            "plasmid_id": context.backbone.plasmid_id,
            "name": context.backbone.name,
            "length_bp": context.backbone.length_bp,
            "assembly_policy": context.backbone.assembly_policy,
            "insertion_regions": list(context.backbone.insertion_regions),
        },
    }


__all__ = [
    "build_final_assembly_context",
    "final_assembly_context_summary",
    "load_final_assembly_context",
]
