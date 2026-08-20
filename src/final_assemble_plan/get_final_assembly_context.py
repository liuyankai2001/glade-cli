"""Manifest adapter for final-assembly planning."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.final_assemble_plan.common import (
    feature_payloads,
    read_genbank,
    resolve_project_file,
    sha256_bytes,
    sha256_sequence,
    stable_json_hash,
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
    if int(manifest.get("revision") or 0) != plasmid_context.manifest_revision:
        raise ValueError("manifest revision 在读取最终组装上下文时发生变化，请重试")

    selection = _mapping(manifest.get("plasmid_selection"), "plasmid_selection")
    if selection.get("schema_version") != "plasmid_selection.v2":
        raise ValueError(
            "manifest 中没有新版质粒选择，请重新运行 plasmid --recommend "
            "并使用 write --plasmid N"
        )
    if selection.get("status") != "selected":
        raise ValueError("plasmid_selection.status 必须为 selected")
    selection_fingerprint = str(
        selection.get("selection_fingerprint") or ""
    ).strip().lower()
    if len(selection_fingerprint) != 64:
        raise ValueError("plasmid_selection 缺少有效 selection_fingerprint")

    vector = _mapping(selection.get("vector"), "plasmid_selection.vector")
    sequence_file = _mapping(
        vector.get("selected_sequence_file"),
        "plasmid_selection.vector.selected_sequence_file",
    )
    backbone_path = resolve_project_file(
        plasmid_context.project_output_path,
        sequence_file.get("path"),
    )
    record, sequence, content = read_genbank(backbone_path)
    actual_file_hash = sha256_bytes(content)
    actual_sequence_hash = sha256_sequence(sequence)
    if actual_file_hash != str(sequence_file.get("file_sha256") or ""):
        raise ValueError("selected_backbone.gb 文件哈希与 manifest 不一致")
    if actual_sequence_hash != str(
        sequence_file.get("sequence_content_sha256") or ""
    ):
        raise ValueError("selected_backbone.gb 序列哈希与 manifest 不一致")
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
        file_sha256=actual_file_hash,
        sequence_content_sha256=actual_sequence_hash,
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
            sequence_sha256=item.sequence_sha256,
            file_sha256=item.file_sha256,
        )
        for item in plasmid_context.constructs
    )
    input_fingerprint = stable_json_hash(
        {
            "target_compound_id": plasmid_context.target_compound_id,
            "parts_selection_fingerprint": (
                plasmid_context.parts_selection_fingerprint
            ),
            "assembled_constructs_fingerprint": (
                plasmid_context.assembled_constructs_fingerprint
            ),
            "plasmid_selection_fingerprint": selection_fingerprint,
            "constructs": [item.fingerprint_payload() for item in constructs],
            "backbone": backbone.fingerprint_payload(),
        }
    )
    return FinalAssemblyContext(
        target_compound_id=plasmid_context.target_compound_id,
        manifest_path=plasmid_context.manifest_path,
        project_output_path=plasmid_context.project_output_path,
        manifest_revision=plasmid_context.manifest_revision,
        parts_selection_fingerprint=plasmid_context.parts_selection_fingerprint,
        assembled_constructs_fingerprint=(
            plasmid_context.assembled_constructs_fingerprint
        ),
        plasmid_selection_fingerprint=selection_fingerprint,
        input_fingerprint=input_fingerprint,
        constructs=constructs,
        backbone=backbone,
    )


def final_assembly_context_summary(config: Any) -> dict[str, Any]:
    context = load_final_assembly_context(config)
    return {
        "target_compound_id": context.target_compound_id,
        "manifest_revision": context.manifest_revision,
        "construct_count": len(context.constructs),
        "parts_design_ids": [item.design_id for item in context.constructs],
        "backbone": {
            "plasmid_id": context.backbone.plasmid_id,
            "name": context.backbone.name,
            "length_bp": context.backbone.length_bp,
            "assembly_policy": context.backbone.assembly_policy,
            "insertion_regions": list(context.backbone.insertion_regions),
        },
        "input_fingerprint": context.input_fingerprint,
    }


__all__ = [
    "final_assembly_context_summary",
    "load_final_assembly_context",
]
