"""Read and validate the manifest inputs required for plasmid recommendation."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio import SeqIO

from src.pathway_analyze.target_id import validate_target_compound_id
from src.plasmid_selection.models import ExpressionConstruct, PlasmidContext
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
    if parts.get("schema_version") != "parts_selection.v1" or parts.get("status") != "selected":
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
        constructs.append(
            ExpressionConstruct(
                design_id=design_id,
                rank=_positive_int(raw.get("rank"), "rank"),
                length_bp=length_bp,
                cassette_count=_positive_int(raw.get("cassette_count"), "cassette_count"),
                component_count=_positive_int(raw.get("component_count"), "component_count"),
                sequence_sha256=sequence_hash,
                file_sha256=file_hash,
                path=path,
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
        "input_fingerprint": context.input_fingerprint,
    }


__all__ = [
    "file_sha256",
    "load_plasmid_context",
    "plasmid_context_summary",
    "stable_json_hash",
]
