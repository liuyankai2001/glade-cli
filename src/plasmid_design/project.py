"""Readiness checks for one existing expression construct."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.plasmid_design.catalog import ModuleCatalog
from src.plasmid_design.errors import DesignError
from src.plasmid_selection.get_plasmid_context import load_plasmid_context
from src.write_manifest.store import read_design_manifest


@dataclass(slots=True)
class ProjectSnapshot:
    manifest: dict
    source: Any
    catalog: ModuleCatalog
    insert_record: SeqRecord


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
    try:
        record = SeqIO.read(source.constructs[0].path, "genbank")
    except (OSError, ValueError) as exc:
        raise DesignError(
            "表达构建 GenBank 无法读取，请重新组装并刷新。", code="project_not_ready"
        ) from exc
    return ProjectSnapshot(manifest, source, catalog, record)
