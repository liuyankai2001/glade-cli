"""Import one user-provided promoter into an expression-parts draft."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from collections.abc import Mapping
from io import StringIO
from pathlib import Path
from typing import Any

from Bio import SeqIO

from src.pathway_analyze.target_id import validate_target_compound_id
from src.protein_to_cds.sequence_constraints import (
    DEFAULT_FORBIDDEN_MOTIFS,
    gc_fraction,
    max_homopolymer_length,
    motif_hits,
)
from src.write_manifest.store import read_design_manifest, update_design_manifest


EXPRESSION_PARTS_DRAFT_SCHEMA_VERSION = "expression_parts_draft.v1"
PROMOTER_UPLOAD_DOWNSTREAM_SECTIONS = (
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
)
_SUPPORTED_SUFFIXES = {".txt", ".fa", ".fasta", ".fna"}
_SAFE_PART_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_DNA_ALPHABET = frozenset("ACGT")


def _stable_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} 必须是 JSON 对象")
    return value


def _positive_int(value: Any, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} 必须是正整数")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} 必须是正整数") from exc
    if parsed < 1:
        raise ValueError(f"{field} 必须是正整数")
    return parsed


def _promoter_args(config: Any) -> tuple[int, str]:
    raw = getattr(config, "promoter", None)
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise ValueError("--promoter 后必须提供表达盒编号和文件名")
    cassette_index = _positive_int(raw[0], "表达盒编号")
    filename = str(raw[1] or "").strip()
    if not filename:
        raise ValueError("启动子文件名不能为空")
    path = Path(filename)
    if path.name != filename or path.is_absolute() or "/" in filename or "\\" in filename:
        raise ValueError("--promoter 只接受 inputs/parts 目录下的文件名")
    if path.suffix.lower() not in _SUPPORTED_SUFFIXES:
        raise ValueError("启动子文件仅支持 .txt、.fa、.fasta 或 .fna")
    return cassette_index, filename


def _read_promoter(path: Path) -> tuple[str, str, str]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"无法读取启动子文件：{path}") from exc
    text = text.lstrip("\ufeff")
    if not text.strip():
        raise ValueError(f"启动子文件为空：{path}")

    if path.suffix.lower() == ".txt":
        part_id = path.stem
        sequence = re.sub(r"\s+", "", text).upper()
        input_format = "plain_text"
    else:
        try:
            records = list(SeqIO.parse(StringIO(text), "fasta"))
        except Exception as exc:
            raise ValueError(f"启动子 FASTA 无法解析：{path}") from exc
        if len(records) != 1:
            raise ValueError("启动子 FASTA 必须恰好包含一条序列")
        part_id = str(records[0].id or "").strip()
        sequence = re.sub(r"\s+", "", str(records[0].seq)).upper()
        input_format = "fasta"

    if not part_id or _SAFE_PART_ID.fullmatch(part_id) is None:
        raise ValueError(
            "启动子 ID 只能包含英文字母、数字、下划线、点和连字符"
        )
    if not sequence:
        raise ValueError("启动子 DNA 序列不能为空")
    invalid = sorted(set(sequence) - _DNA_ALPHABET)
    if invalid:
        raise ValueError("启动子 DNA 只能包含 A/C/G/T；发现：" + ", ".join(invalid))
    return part_id, sequence, input_format


def _canonical_fasta(part_id: str, sequence: str) -> bytes:
    lines = [sequence[index : index + 80] for index in range(0, len(sequence), 80)]
    return (f">{part_id} user_uploaded_promoter\n" + "\n".join(lines) + "\n").encode(
        "utf-8"
    )


def _write_bytes_atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _relative(project_root: Path, path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(project_root.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError("启动子快照路径超出当前项目输出目录") from exc


def _current_expression_box(
    manifest: Mapping[str, Any],
) -> tuple[str, str, list[dict[str, Any]]]:
    selection = _mapping(
        manifest.get("expression_box_selection"),
        "expression_box_selection",
    )
    if selection.get("schema_version") != "expression_box_selection.v1":
        raise ValueError("请先写入有效的表达盒方案")
    box_fingerprint = str(selection.get("selected_design_fingerprint") or "").strip()
    if not box_fingerprint:
        raise ValueError("表达盒方案缺少 fingerprint")
    cds_selection = _mapping(manifest.get("cds_selection"), "cds_selection")
    cds_fingerprint = str(cds_selection.get("source_fingerprint") or "").strip()
    if not cds_fingerprint:
        raise ValueError("CDS 选择缺少 source_fingerprint")
    raw_cassettes = selection.get("cassettes")
    if not isinstance(raw_cassettes, list) or not raw_cassettes:
        raise ValueError("表达盒方案没有有效分组")
    cassettes: list[dict[str, Any]] = []
    for expected_index, raw in enumerate(raw_cassettes, start=1):
        item = _mapping(raw, f"expression_box_selection.cassettes[{expected_index - 1}]")
        cassette_index = _positive_int(item.get("cassette_index"), "cassette_index")
        if cassette_index != expected_index:
            raise ValueError("表达盒编号必须从 1 开始连续排列")
        proteins = item.get("protein_accessions")
        if not isinstance(proteins, list) or not proteins:
            raise ValueError(f"表达盒 {cassette_index} 没有蛋白")
        cassettes.append(
            {
                "cassette_index": cassette_index,
                "protein_accessions": [
                    str(accession or "").strip().upper() for accession in proteins
                ],
                "promoter": None,
                "rbs_by_accession": {},
                "terminator": None,
            }
        )
    return box_fingerprint, cds_fingerprint, cassettes


def _merge_existing_draft(
    manifest: Mapping[str, Any],
    cassettes: list[dict[str, Any]],
    *,
    box_fingerprint: str,
    cds_fingerprint: str,
) -> None:
    raw_draft = manifest.get("expression_parts_draft")
    if not isinstance(raw_draft, Mapping):
        return
    if raw_draft.get("schema_version") != EXPRESSION_PARTS_DRAFT_SCHEMA_VERSION:
        return
    source = raw_draft.get("source")
    if not isinstance(source, Mapping) or (
        source.get("expression_box_selection_fingerprint") != box_fingerprint
        or source.get("cds_selection_source_fingerprint") != cds_fingerprint
    ):
        return
    raw_cassettes = raw_draft.get("cassettes")
    if not isinstance(raw_cassettes, list):
        return
    existing = {
        item.get("cassette_index"): item
        for item in raw_cassettes
        if isinstance(item, Mapping)
    }
    for cassette in cassettes:
        previous = existing.get(cassette["cassette_index"])
        if not isinstance(previous, Mapping):
            continue
        if previous.get("protein_accessions") != cassette["protein_accessions"]:
            continue
        promoter = previous.get("promoter")
        cassette["promoter"] = dict(promoter) if isinstance(promoter, Mapping) else None
        raw_rbs = previous.get("rbs_by_accession")
        cassette["rbs_by_accession"] = dict(raw_rbs) if isinstance(raw_rbs, Mapping) else {}
        terminator = previous.get("terminator")
        cassette["terminator"] = (
            dict(terminator) if isinstance(terminator, Mapping) else None
        )


def _assert_part_id_consistent(
    cassettes: list[dict[str, Any]],
    part_id: str,
    sequence_sha256: str,
    *,
    replacing_cassette_index: int,
) -> None:
    for cassette in cassettes:
        parts: list[Any] = [cassette.get("terminator")]
        if cassette.get("cassette_index") != replacing_cassette_index:
            parts.append(cassette.get("promoter"))
        rbs = cassette.get("rbs_by_accession")
        if isinstance(rbs, Mapping):
            parts.extend(rbs.values())
        for raw in parts:
            if not isinstance(raw, Mapping) or raw.get("part_id") != part_id:
                continue
            if raw.get("sequence_sha256") != sequence_sha256:
                raise ValueError(f"元件 ID {part_id} 已对应另一条 DNA 序列")


def _draft_payload(
    *,
    target_compound_id: str,
    box_fingerprint: str,
    cds_fingerprint: str,
    cassettes: list[dict[str, Any]],
) -> dict[str, Any]:
    promoter_count = sum(item.get("promoter") is not None for item in cassettes)
    content = {
        "target_compound_id": target_compound_id,
        "source": {
            "expression_box_selection_fingerprint": box_fingerprint,
            "cds_selection_source_fingerprint": cds_fingerprint,
        },
        "cassettes": cassettes,
    }
    return {
        "schema_version": EXPRESSION_PARTS_DRAFT_SCHEMA_VERSION,
        "status": "partial",
        **content,
        "summary": {
            "cassette_count": len(cassettes),
            "promoter_count": promoter_count,
            "required_promoter_count": len(cassettes),
            "rbs_count": sum(
                len(item.get("rbs_by_accession") or {}) for item in cassettes
            ),
            "terminator_count": sum(
                item.get("terminator") is not None for item in cassettes
            ),
        },
        "draft_fingerprint": _stable_hash(content),
    }


def upload_expression_promoter(config: Any) -> dict[str, Any]:
    """Upload or replace the promoter for one selected expression cassette."""

    cassette_index, filename = _promoter_args(config)
    target_compound_id = validate_target_compound_id(config.target_name)
    inputs_root = Path(config.inputs_dir).expanduser().resolve()
    input_path = (inputs_root / "parts" / filename).resolve()
    parts_root = (inputs_root / "parts").resolve()
    if input_path.parent != parts_root:
        raise ValueError("启动子文件必须直接位于 inputs/parts 目录")
    if not input_path.is_file():
        raise FileNotFoundError(f"未找到启动子文件：{input_path}")
    part_id, sequence, input_format = _read_promoter(input_path)

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if recorded_target != target_compound_id:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target or '空'}，"
            f"与当前输入目标 {target_compound_id} 不一致"
        )
    try:
        current_revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision 必须是整数") from exc
    box_fingerprint, cds_fingerprint, cassettes = _current_expression_box(manifest)
    if cassette_index > len(cassettes):
        raise ValueError(
            f"不存在表达盒 {cassette_index}；可用编号：{list(range(1, len(cassettes) + 1))}"
        )
    _merge_existing_draft(
        manifest,
        cassettes,
        box_fingerprint=box_fingerprint,
        cds_fingerprint=cds_fingerprint,
    )

    sequence_sha256 = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    _assert_part_id_consistent(
        cassettes,
        part_id,
        sequence_sha256,
        replacing_cassette_index=cassette_index,
    )
    fasta_bytes = _canonical_fasta(part_id, sequence)
    snapshot_path = (
        project_root
        / "expression_box"
        / "uploaded_parts"
        / "promoters"
        / f"cassette_{cassette_index:03d}"
        / f"{part_id}.{sequence_sha256[:12]}.fasta"
    ).resolve()
    snapshot_relative = _relative(project_root, snapshot_path)
    promoter = {
        "part_id": part_id,
        "role": "promoter",
        "source": "user_uploaded",
        "source_input_file": f"parts/{filename}",
        "input_format": input_format,
        "sequence_file": {
            "path": snapshot_relative,
            "format": "fasta",
            "file_sha256": _sha256_bytes(fasta_bytes),
        },
        "sequence_sha256": sequence_sha256,
        "length_bp": len(sequence),
        "metrics": {
            "gc_percent": round(100.0 * gc_fraction(sequence), 8),
            "forbidden_site_hits": motif_hits(sequence, DEFAULT_FORBIDDEN_MOTIFS),
            "max_homopolymer": max_homopolymer_length(sequence),
        },
    }
    previous_promoter = cassettes[cassette_index - 1].get("promoter")
    cassettes[cassette_index - 1]["promoter"] = promoter
    payload = _draft_payload(
        target_compound_id=target_compound_id,
        box_fingerprint=box_fingerprint,
        cds_fingerprint=cds_fingerprint,
        cassettes=cassettes,
    )
    current_draft = manifest.get("expression_parts_draft")
    draft_changed = not (
        isinstance(current_draft, Mapping) and dict(current_draft) == payload
    )
    downstream_present = any(
        field in manifest for field in PROMOTER_UPLOAD_DOWNSTREAM_SECTIONS
    )
    manifest_changed = draft_changed or downstream_present

    snapshot_existed = snapshot_path.is_file()
    snapshot_matches = (
        snapshot_existed and _sha256_bytes(snapshot_path.read_bytes()) == _sha256_bytes(fasta_bytes)
    )
    previous_sequence_file = (
        previous_promoter.get("sequence_file")
        if isinstance(previous_promoter, Mapping)
        else None
    )
    snapshot_repair = (
        not snapshot_matches
        and isinstance(previous_sequence_file, Mapping)
        and previous_sequence_file.get("path") == snapshot_relative
    )
    if not snapshot_matches:
        _write_bytes_atomic(snapshot_path, fasta_bytes)
    try:
        if manifest_changed:
            updated = update_design_manifest(
                manifest_path,
                target_compound_id=target_compound_id,
                sections={"expression_parts_draft": payload},
                discard_sections=PROMOTER_UPLOAD_DOWNSTREAM_SECTIONS,
                expected_revision=current_revision,
            )
        else:
            updated = manifest
    except Exception:
        if not snapshot_existed:
            snapshot_path.unlink(missing_ok=True)
        raise

    action = (
        "新增"
        if previous_promoter is None
        else "未变化"
        if not draft_changed
        else "替换"
    )
    return {
        "运行成功": True,
        "目标化合物": target_compound_id,
        "表达盒编号": cassette_index,
        "操作": action,
        "启动子ID": part_id,
        "输入文件": f"parts/{filename}",
        "序列长度_bp": len(sequence),
        "GC百分比": promoter["metrics"]["gc_percent"],
        "序列SHA256": sequence_sha256,
        "快照文件": str(snapshot_path),
        "快照是否写入": not snapshot_matches,
        "快照是否修复": snapshot_repair,
        "已上传启动子数": payload["summary"]["promoter_count"],
        "表达盒总数": payload["summary"]["required_promoter_count"],
        "清单是否更新": manifest_changed,
        "清单版本": updated["revision"],
        "清单文件": str(manifest_path),
    }


__all__ = [
    "EXPRESSION_PARTS_DRAFT_SCHEMA_VERSION",
    "PROMOTER_UPLOAD_DOWNSTREAM_SECTIONS",
    "upload_expression_promoter",
]
