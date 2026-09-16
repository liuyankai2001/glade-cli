"""Readiness and source fingerprints for one existing expression construct."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.plasmid_design.catalog import ModuleCatalog
from src.plasmid_design.errors import DesignError
from src.plasmid_selection.get_plasmid_context import (
    load_plasmid_context,
    stable_json_hash,
)
from src.write_manifest.store import read_design_manifest


@dataclass(slots=True)
class ProjectSnapshot:
    manifest: dict
    source: Any
    catalog: ModuleCatalog
    insert_record: SeqRecord
    fingerprint: str


def read_project(config: Any) -> ProjectSnapshot:
    manifest = read_design_manifest(config.manifest_output_path)
    assembled = manifest.get("assembled_expression_constructs")
    if not isinstance(assembled, dict) or assembled.get("status") != "assembled":
        raise DesignError(
            "请先在命令行运行 expression --assemble，完成表达构建后再设计质粒。",
            code="construct_missing",
        )
    constructs = assembled.get("constructs")
    parts = manifest.get("parts_selection", {})
    cds = manifest.get("cds_selection", {})
    if not isinstance(parts, dict) or not isinstance(cds, dict):
        raise DesignError(
            "项目的 CDS 或表达元件登记格式有误，请重新完成表达构建。",
            code="invalid_manifest",
        )
    if (
        not isinstance(constructs, list)
        or len(constructs) != 1
        or len(parts.get("selected_design_ids", [])) != 1
    ):
        raise DesignError(
            "网页首版需要一个已确认的完整表达构建；请在命令行保留一个方案并重新组装。",
            code="construct_count",
        )
    try:
        source = load_plasmid_context(config)
        catalog = ModuleCatalog(config.data_dir)
    except (OSError, ValueError) as exc:
        raise DesignError(
            f"当前表达构建或组件库需更新：{exc}", code="project_not_ready"
        ) from exc
    if not (
        source.host_key.lower().startswith("ecoli")
        or source.host_key.lower() == "escherichia_coli"
    ):
        raise DesignError(
            "当前组件库首版支持大肠杆菌；该项目宿主尚无适用的复制模块。",
            code="unsupported_host",
        )
    if source.manifest_revision != int(manifest.get("revision", 0)):
        raise DesignError("读取期间项目状态已变化，请刷新后重试。", code="stale_source")
    try:
        record = SeqIO.read(source.constructs[0].path, "genbank")
    except (OSError, ValueError) as exc:
        raise DesignError(
            "表达构建 GenBank 无法读取，请重新组装并刷新。", code="project_not_ready"
        ) from exc
    from src.final_assemble_execute.common import sha256_sequence

    if sha256_sequence(str(record.seq)) != source.constructs[0].sequence_sha256:
        raise DesignError("表达构建文件已变化，请重新组装并刷新。", code="stale_source")
    fingerprint = stable_json_hash(
        {
            "context": source.input_fingerprint,
            "catalog": catalog.fingerprint,
            "cds_selection": manifest.get("cds_selection"),
            "parts_selection": parts,
            "assembled_expression_constructs": assembled,
        }
    )
    return ProjectSnapshot(manifest, source, catalog, record, fingerprint)


def authenticate_request(snapshot: ProjectSnapshot, request: dict) -> None:
    if (
        request.get("expected_revision") != snapshot.source.manifest_revision
        or request.get("source_fingerprint") != snapshot.fingerprint
    ):
        raise DesignError(
            "项目已更新，请刷新页面后重新检查与生成。", code="stale_source"
        )
