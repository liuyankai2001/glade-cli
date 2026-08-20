"""Read and validate the manifest inputs required for plasmid recommendation."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio import SeqIO

from src.expression_box.expression_burden import (
    validate_expression_burden_summary,
)
from src.pathway_analyze.target_id import validate_target_compound_id
from src.plasmid_selection.models import (
    ExpressionBurden,
    ExpressionConstruct,
    PlasmidContext,
)
from src.write_manifest.store import read_design_manifest


def stable_json_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少 {name}")
    return value


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"manifest {name} 必须是正整数")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"manifest {name} 必须是正整数") from exc
    if parsed < 1:
        raise ValueError(f"manifest {name} 必须是正整数")
    return parsed


def _resolve_construct_path(project_output: Path, value: Any) -> Path:
    text = str(value or "").strip()
    if not text:
        raise ValueError("assembled_expression_constructs 包含空文件路径")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = project_output / path
    path = path.resolve()
    try:
        path.relative_to(project_output.resolve())
    except ValueError as exc:
        raise ValueError(f"表达构建文件不在当前项目输出目录内：{path}") from exc
    return path


def _validate_construct_file(
    path: Path,
    *,
    expected_file_hash: str,
    expected_sequence_hash: str,
    expected_length: int,
) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"缺少表达构建 GenBank 文件：{path}")
    if file_sha256(path) != expected_file_hash:
        raise ValueError(f"表达构建文件哈希不匹配：{path}")
    try:
        record = SeqIO.read(path, "genbank")
    except Exception as exc:
        raise ValueError(f"无法解析表达构建 GenBank 文件：{path}") from exc
    sequence = str(record.seq).upper()
    if len(sequence) != expected_length:
        raise ValueError(f"表达构建长度与 manifest 不一致：{path}")
    if set(sequence) - set("ACGT"):
        raise ValueError(f"表达构建包含非 A/C/G/T 碱基：{path}")
    actual_sequence_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if actual_sequence_hash != expected_sequence_hash:
        raise ValueError(f"表达构建序列哈希与 manifest 不一致：{path}")


def load_plasmid_context(config: Any) -> PlasmidContext:
    """Load current selected constructs; no separate readiness flag is used."""

    target = validate_target_compound_id(config.target_name)
    project_output = Path(config.project_output_path).expanduser().resolve()
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    if str(manifest.get("target_compound_id") or "") != target:
        raise ValueError("design manifest 与当前目标化合物不一致")

    parts = _mapping(manifest.get("parts_selection"), "parts_selection")
    if parts.get("schema_version") != "parts_selection.v2":
        raise ValueError(
            "manifest 的表达元件选择缺少数值负担模型，请依次重新运行："
            "expression --design --parts；write --expression-parts 1:12"
        )
    if parts.get("status") != "selected":
        raise ValueError("manifest 中没有有效的表达元件选择")
    parts_fingerprint = str(parts.get("selection_fingerprint") or "").strip()
    if not parts_fingerprint:
        raise ValueError("parts_selection 缺少 selection_fingerprint")
    raw_ids = parts.get("selected_design_ids")
    if not isinstance(raw_ids, list) or not raw_ids:
        raise ValueError("parts_selection 没有 selected_design_ids")
    selected_ids = tuple(_positive_int(item, "selected_design_ids") for item in raw_ids)
    if len(selected_ids) != len(set(selected_ids)):
        raise ValueError("parts_selection 包含重复方案编号")
    raw_references = parts.get("design_references")
    if not isinstance(raw_references, list) or not raw_references:
        raise ValueError("parts_selection 没有 design_references")
    burden_by_design: dict[int, ExpressionBurden] = {}
    for reference in raw_references:
        if not isinstance(reference, Mapping):
            raise ValueError("parts_selection 包含无效 design_reference")
        design_id = _positive_int(reference.get("design_id"), "design_id")
        burden = reference.get("expression_burden")
        if not isinstance(burden, Mapping):
            raise ValueError(
                f"表达元件方案 {design_id} 缺少数值负担；"
                "请重新运行 write --expression-parts"
            )
        try:
            validate_expression_burden_summary(burden)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"表达元件方案 {design_id} 的负担摘要无效：{exc}"
            ) from exc
        if design_id in burden_by_design:
            raise ValueError("parts_selection 包含重复 design_reference")
        burden_by_design[design_id] = ExpressionBurden(
            schema_version=str(burden["schema_version"]),
            model_version=str(burden["model_version"]),
            score=float(burden["score"]),
            level=str(burden["level"]),
            confidence=str(burden["confidence"]),
            raw_load_units=float(burden["raw_load_units"]),
            reference_load_units=float(burden["reference_load_units"]),
            gene_count=int(burden["gene_count"]),
            cassette_count=int(burden["cassette_count"]),
            total_cds_length_nt=int(burden["total_cds_length_nt"]),
            minimum_ostir_reference_count=int(
                burden["minimum_ostir_reference_count"]
            ),
            fingerprint=str(burden["fingerprint"]),
            warnings=tuple(str(item) for item in burden.get("warnings", [])),
        )
    if set(burden_by_design) != set(selected_ids):
        raise ValueError("parts_selection 负担摘要没有覆盖全部已选方案")

    assembled = _mapping(
        manifest.get("assembled_expression_constructs"),
        "assembled_expression_constructs",
    )
    if (
        assembled.get("schema_version") != "assembled_expression_constructs.v1"
        or assembled.get("status") != "assembled"
    ):
        raise ValueError("manifest 中没有有效的完整串联表达构建")
    if assembled.get("source_parts_selection_fingerprint") != parts_fingerprint:
        raise ValueError("完整表达构建与当前表达元件选择不一致")
    raw_constructs = assembled.get("constructs")
    if not isinstance(raw_constructs, list) or not raw_constructs:
        raise ValueError("assembled_expression_constructs 没有 constructs")

    constructs: list[ExpressionConstruct] = []
    for raw in raw_constructs:
        if not isinstance(raw, Mapping):
            raise ValueError("assembled_expression_constructs 包含无效构建")
        design_id = _positive_int(raw.get("parts_design_id"), "parts_design_id")
        length_bp = _positive_int(raw.get("length_bp"), "length_bp")
        sequence_hash = str(raw.get("sequence_sha256") or "").strip().lower()
        file_hash = str(raw.get("file_sha256") or "").strip().lower()
        if len(sequence_hash) != 64 or len(file_hash) != 64:
            raise ValueError(f"表达构建 {design_id} 缺少有效 SHA-256")
        path = _resolve_construct_path(project_output, raw.get("path"))
        _validate_construct_file(
            path,
            expected_file_hash=file_hash,
            expected_sequence_hash=sequence_hash,
            expected_length=length_bp,
        )
        burden = burden_by_design.get(design_id)
        if burden is None:
            raise ValueError(f"表达构建 {design_id} 缺少负担摘要")
        raw_ranges = raw.get("cassette_ranges")
        raw_ranges = raw_ranges if isinstance(raw_ranges, list) else []
        gene_count = sum(
            len(item.get("protein_accessions") or [])
            for item in raw_ranges
            if isinstance(item, Mapping)
            and isinstance(item.get("protein_accessions"), list)
        )
        cassette_count = _positive_int(raw.get("cassette_count"), "cassette_count")
        if burden.cassette_count != cassette_count or (
            gene_count and burden.gene_count != gene_count
        ):
            raise ValueError(f"表达构建 {design_id} 与负担摘要计数不一致")
        constructs.append(
            ExpressionConstruct(
                design_id=design_id,
                rank=_positive_int(raw.get("rank"), "rank"),
                length_bp=length_bp,
                cassette_count=cassette_count,
                component_count=_positive_int(raw.get("component_count"), "component_count"),
                sequence_sha256=sequence_hash,
                file_sha256=file_hash,
                path=path,
                burden=burden,
            )
        )
    constructs.sort(key=lambda item: (item.rank, item.design_id))
    construct_ids = tuple(item.design_id for item in constructs)
    if set(construct_ids) != set(selected_ids) or len(construct_ids) != len(selected_ids):
        raise ValueError("完整表达构建没有覆盖全部已选表达元件方案")
    recorded_count = _positive_int(assembled.get("design_count"), "design_count")
    if recorded_count != len(constructs):
        raise ValueError("assembled_expression_constructs.design_count 不一致")

    assembled_fingerprint = stable_json_hash(
        {
            "schema_version": assembled.get("schema_version"),
            "source_parts_selection_fingerprint": parts_fingerprint,
            "junction_policy": assembled.get("junction_policy"),
            "linker_sequence": assembled.get("linker_sequence"),
            "constructs": [item.fingerprint_payload() for item in constructs],
        }
    )
    cds = manifest.get("cds_selection")
    cds = cds if isinstance(cds, Mapping) else {}
    host = cds.get("host")
    host = host if isinstance(host, Mapping) else {}
    host_key = str(host.get("chassis_key") or "ecoli_mg1655").strip()
    host_name = str(
        host.get("name") or "Escherichia coli str. K-12 substr. MG1655"
    ).strip()
    input_fingerprint = stable_json_hash(
        {
            "target_compound_id": target,
            "host_key": host_key,
            "host_name": host_name,
            "parts_selection_fingerprint": parts_fingerprint,
            "assembled_constructs_fingerprint": assembled_fingerprint,
            "selected_design_ids": sorted(selected_ids),
        }
    )
    return PlasmidContext(
        target_compound_id=target,
        manifest_path=manifest_path,
        project_output_path=project_output,
        manifest_revision=int(manifest.get("revision") or 0),
        host_key=host_key,
        host_name=host_name,
        parts_selection_fingerprint=parts_fingerprint,
        assembled_constructs_fingerprint=assembled_fingerprint,
        input_fingerprint=input_fingerprint,
        selected_design_ids=selected_ids,
        constructs=tuple(constructs),
    )


def plasmid_context_summary(config: Any) -> dict[str, Any]:
    context = load_plasmid_context(config)
    return {
        "target_compound_id": context.target_compound_id,
        "manifest_revision": context.manifest_revision,
        "host": {"key": context.host_key, "name": context.host_name},
        "selected_design_ids": list(context.selected_design_ids),
        "construct_count": len(context.constructs),
        "insert_length_range_bp": {
            "minimum": context.minimum_insert_length_bp,
            "maximum": context.maximum_insert_length_bp,
        },
        "maximum_cassette_count": context.maximum_cassette_count,
        "expression_burden": {
            "score_range": {
                "minimum": context.minimum_burden_score,
                "maximum": context.maximum_burden_score,
            },
            "levels": sorted({item.burden.level for item in context.constructs}),
            "model_version": context.constructs[0].burden.model_version,
        },
        "input_fingerprint": context.input_fingerprint,
    }


__all__ = [
    "file_sha256",
    "load_plasmid_context",
    "plasmid_context_summary",
    "stable_json_hash",
]
