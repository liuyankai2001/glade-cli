"""Import user-provided boundary parts into an expression-parts draft."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Callable, Mapping
from io import StringIO
from pathlib import Path
from typing import Any

from Bio import SeqIO

from src.expression_box.config import (
    RBS_CONTEXT_CDS_PREFIX_NT,
    RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT,
)
from src.expression_box.ostir_adapter import predict_rbs_context
from src.expression_box.parts_manifest_adapter import load_expression_parts_context
from src.expression_box.parts_models import RbsPrediction
from src.pathway_analyze.target_id import validate_target_compound_id
from src.protein_to_cds.sequence_constraints import (
    DEFAULT_FORBIDDEN_MOTIFS,
    gc_fraction,
    max_homopolymer_length,
    motif_hits,
)
from src.write_manifest.store import read_design_manifest, update_design_manifest


EXPRESSION_PARTS_DRAFT_SCHEMA_VERSION = "expression_parts_draft.v1"
EXPRESSION_PART_UPLOAD_DOWNSTREAM_SECTIONS = (
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
)
PROMOTER_UPLOAD_DOWNSTREAM_SECTIONS = EXPRESSION_PART_UPLOAD_DOWNSTREAM_SECTIONS
_SUPPORTED_SUFFIXES = {".txt", ".fa", ".fasta", ".fna"}
_SAFE_PART_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_DNA_ALPHABET = frozenset("ACGT")
_ROLE_LABELS = {
    "promoter": "启动子",
    "rbs": "RBS",
    "terminator": "终止子",
}


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


def _part_args(config: Any, role: str) -> tuple[int, str]:
    label = _ROLE_LABELS[role]
    option = f"--{role}"
    raw = getattr(config, role, None)
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise ValueError(f"{option} 后必须提供表达盒编号和文件名")
    cassette_index = _positive_int(raw[0], "表达盒编号")
    filename = str(raw[1] or "").strip()
    if not filename:
        raise ValueError(f"{label}文件名不能为空")
    path = Path(filename)
    if path.name != filename or path.is_absolute() or "/" in filename or "\\" in filename:
        raise ValueError(f"{option} 只接受 inputs/parts 目录下的文件名")
    if path.suffix.lower() not in _SUPPORTED_SUFFIXES:
        raise ValueError(f"{label}文件仅支持 .txt、.fa、.fasta 或 .fna")
    return cassette_index, filename


def _rbs_args(config: Any) -> tuple[str, str]:
    raw = getattr(config, "rbs", None)
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise ValueError("--rbs 后必须提供蛋白 accession 和文件名")
    accession = str(raw[0] or "").strip().upper()
    if not accession:
        raise ValueError("--rbs 的蛋白 accession 不能为空")
    filename = str(raw[1] or "").strip()
    if not filename:
        raise ValueError("RBS 文件名不能为空")
    path = Path(filename)
    if path.name != filename or path.is_absolute() or "/" in filename or "\\" in filename:
        raise ValueError("--rbs 只接受 inputs/parts 目录下的文件名")
    if path.suffix.lower() not in _SUPPORTED_SUFFIXES:
        raise ValueError("RBS 文件仅支持 .txt、.fa、.fasta 或 .fna")
    return accession, filename


def _read_part(path: Path, role: str) -> tuple[str, str, str]:
    label = _ROLE_LABELS[role]
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"无法读取{label}文件：{path}") from exc
    text = text.lstrip("\ufeff")
    if not text.strip():
        raise ValueError(f"{label}文件为空：{path}")

    if path.suffix.lower() == ".txt":
        part_id = path.stem
        sequence = re.sub(r"\s+", "", text).upper()
        input_format = "plain_text"
    else:
        try:
            records = list(SeqIO.parse(StringIO(text), "fasta"))
        except Exception as exc:
            raise ValueError(f"{label} FASTA 无法解析：{path}") from exc
        if len(records) != 1:
            raise ValueError(f"{label} FASTA 必须恰好包含一条序列")
        part_id = str(records[0].id or "").strip()
        sequence = re.sub(r"\s+", "", str(records[0].seq)).upper()
        input_format = "fasta"

    if not part_id or _SAFE_PART_ID.fullmatch(part_id) is None:
        raise ValueError(
            f"{label} ID 只能包含英文字母、数字、下划线、点和连字符"
        )
    if not sequence:
        raise ValueError(f"{label} DNA 序列不能为空")
    invalid = sorted(set(sequence) - _DNA_ALPHABET)
    if invalid:
        raise ValueError(f"{label} DNA 只能包含 A/C/G/T；发现：" + ", ".join(invalid))
    return part_id, sequence, input_format


def _canonical_fasta(part_id: str, sequence: str, role: str) -> bytes:
    lines = [sequence[index : index + 80] for index in range(0, len(sequence), 80)]
    return (f">{part_id} user_uploaded_{role}\n" + "\n".join(lines) + "\n").encode(
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
        raise ValueError("表达元件快照路径超出当前项目输出目录") from exc


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
    replacing_role: str,
    replacing_accession: str | None = None,
) -> None:
    for cassette in cassettes:
        parts: list[Any] = []
        for role in ("promoter", "terminator"):
            replacing_current = (
                cassette.get("cassette_index") == replacing_cassette_index
                and role == replacing_role
            )
            if not replacing_current:
                parts.append(cassette.get(role))
        rbs = cassette.get("rbs_by_accession")
        if isinstance(rbs, Mapping):
            for accession, raw in rbs.items():
                replacing_current = (
                    cassette.get("cassette_index") == replacing_cassette_index
                    and replacing_role == "rbs"
                    and str(accession).upper() == replacing_accession
                )
                if not replacing_current:
                    parts.append(raw)
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


def _commit_uploaded_part(
    *,
    manifest_path: Path,
    manifest: Mapping[str, Any],
    current_revision: int,
    project_root: Path,
    target_compound_id: str,
    box_fingerprint: str,
    cds_fingerprint: str,
    cassettes: list[dict[str, Any]],
    cassette_index: int,
    role: str,
    filename: str,
    part_id: str,
    sequence: str,
    input_format: str,
    accession: str | None = None,
    extra_fields: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    label = _ROLE_LABELS[role]
    sequence_sha256 = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    _assert_part_id_consistent(
        cassettes,
        part_id,
        sequence_sha256,
        replacing_cassette_index=cassette_index,
        replacing_role=role,
        replacing_accession=accession,
    )
    fasta_bytes = _canonical_fasta(part_id, sequence, role)
    snapshot_dir = (
        project_root
        / "expression_box"
        / "uploaded_parts"
        / ("rbs" if role == "rbs" else f"{role}s")
        / f"cassette_{cassette_index:03d}"
    )
    if accession is not None:
        snapshot_dir = snapshot_dir / accession
    snapshot_path = (
        snapshot_dir / f"{part_id}.{sequence_sha256[:12]}.fasta"
    ).resolve()
    snapshot_relative = _relative(project_root, snapshot_path)
    part = {
        "part_id": part_id,
        "role": role,
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
        **dict(extra_fields or {}),
    }
    if role == "rbs":
        rbs_by_accession = cassettes[cassette_index - 1]["rbs_by_accession"]
        previous_part = rbs_by_accession.get(accession)
        rbs_by_accession[accession] = part
    else:
        previous_part = cassettes[cassette_index - 1].get(role)
        cassettes[cassette_index - 1][role] = part

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
        field in manifest for field in EXPRESSION_PART_UPLOAD_DOWNSTREAM_SECTIONS
    )
    manifest_changed = draft_changed or downstream_present

    snapshot_existed = snapshot_path.is_file()
    snapshot_matches = (
        snapshot_existed
        and _sha256_bytes(snapshot_path.read_bytes()) == _sha256_bytes(fasta_bytes)
    )
    previous_sequence_file = (
        previous_part.get("sequence_file")
        if isinstance(previous_part, Mapping)
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
                discard_sections=EXPRESSION_PART_UPLOAD_DOWNSTREAM_SECTIONS,
                expected_revision=current_revision,
            )
        else:
            updated = dict(manifest)
    except Exception:
        if not snapshot_existed:
            snapshot_path.unlink(missing_ok=True)
        raise

    action = (
        "新增"
        if previous_part is None
        else "未变化"
        if not draft_changed
        else "替换"
    )
    result = {
        "运行成功": True,
        "目标化合物": target_compound_id,
        "表达盒编号": cassette_index,
        "操作": action,
        f"{label}ID": part_id,
        "输入文件": f"parts/{filename}",
        "序列长度_bp": len(sequence),
        "GC百分比": part["metrics"]["gc_percent"],
        "序列SHA256": sequence_sha256,
        "快照文件": str(snapshot_path),
        "快照是否写入": not snapshot_matches,
        "快照是否修复": snapshot_repair,
        "已上传启动子数": payload["summary"]["promoter_count"],
        "已上传RBS数": payload["summary"]["rbs_count"],
        "已上传终止子数": payload["summary"]["terminator_count"],
        "表达盒总数": payload["summary"]["required_promoter_count"],
        "清单是否更新": manifest_changed,
        "清单版本": updated["revision"],
        "清单文件": str(manifest_path),
    }
    if accession is not None:
        result["蛋白ID"] = accession
    ostir = part.get("ostir")
    if isinstance(ostir, Mapping):
        result["翻译起始率"] = ostir.get("translation_initiation_rate")
        result["非预期起始位点数"] = ostir.get("unintended_start_count")
    return result


def _upload_expression_boundary_part(config: Any, role: str) -> dict[str, Any]:
    label = _ROLE_LABELS[role]
    cassette_index, filename = _part_args(config, role)
    target_compound_id = validate_target_compound_id(config.target_name)
    inputs_root = Path(config.inputs_dir).expanduser().resolve()
    input_path = (inputs_root / "parts" / filename).resolve()
    parts_root = (inputs_root / "parts").resolve()
    if input_path.parent != parts_root:
        raise ValueError(f"{label}文件必须直接位于 inputs/parts 目录")
    if not input_path.is_file():
        raise FileNotFoundError(f"未找到{label}文件：{input_path}")
    part_id, sequence, input_format = _read_part(input_path, role)

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

    return _commit_uploaded_part(
        manifest_path=manifest_path,
        manifest=manifest,
        current_revision=current_revision,
        project_root=project_root,
        target_compound_id=target_compound_id,
        box_fingerprint=box_fingerprint,
        cds_fingerprint=cds_fingerprint,
        cassettes=cassettes,
        cassette_index=cassette_index,
        role=role,
        filename=filename,
        part_id=part_id,
        sequence=sequence,
        input_format=input_format,
    )


def _rbs_location(
    cassettes: list[dict[str, Any]],
    accession: str,
) -> tuple[int, int]:
    matches: list[tuple[int, int]] = []
    available: list[str] = []
    for cassette in cassettes:
        cassette_index = int(cassette["cassette_index"])
        for cds_index, value in enumerate(cassette["protein_accessions"]):
            normalized = str(value).upper()
            available.append(normalized)
            if normalized == accession:
                matches.append((cassette_index, cds_index))
    if len(matches) != 1:
        if not matches:
            raise ValueError(
                f"当前表达盒中不存在蛋白 {accession}；可用蛋白：{sorted(available)}"
            )
        raise ValueError(f"蛋白 {accession} 在表达盒中出现多次，无法确定 RBS 位置")
    return matches[0]


def _predict_uploaded_rbs(
    *,
    manifest_path: Path,
    project_root: Path,
    cassette_index: int,
    cds_index: int,
    accession: str,
    part_id: str,
    rbs_sequence: str,
    predictor: Callable[..., RbsPrediction],
) -> dict[str, Any]:
    context = load_expression_parts_context(manifest_path, project_root)
    cassette_matches = [
        cassette
        for cassette in context.cassettes
        if cassette.cassette_index == cassette_index
    ]
    if len(cassette_matches) != 1:
        raise ValueError(f"无法读取表达盒 {cassette_index} 的 CDS 上下文")
    cassette = cassette_matches[0]
    if cds_index >= len(cassette.cds) or cassette.cds[cds_index].accession != accession:
        raise ValueError(f"蛋白 {accession} 的 CDS 顺序与表达盒方案不一致")
    cds = cassette.cds[cds_index]
    previous = cassette.cds[cds_index - 1] if cds_index > 0 else None
    upstream = (
        previous.sequence[-RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT:]
        if previous is not None
        else ""
    )
    prediction_sequence = (
        upstream + rbs_sequence + cds.sequence[:RBS_CONTEXT_CDS_PREFIX_NT]
    )
    intended_start = len(upstream) + len(rbs_sequence) + 1
    prediction = predictor(
        sequence=prediction_sequence,
        intended_start_position=intended_start,
        accession=accession,
        part_id=part_id,
    )
    if not isinstance(prediction, RbsPrediction):
        raise ValueError("OSTIR 预测器没有返回有效的 RbsPrediction")
    expected_context_hash = hashlib.sha256(
        prediction_sequence.encode("utf-8")
    ).hexdigest()
    if prediction.part_id != part_id or prediction.accession != accession:
        raise ValueError("OSTIR 预测结果与当前 RBS 或蛋白不一致")
    if prediction.intended_start_position != intended_start:
        raise ValueError("OSTIR 预测结果的预期起始位置不一致")
    if prediction.context_sha256 != expected_context_hash:
        raise ValueError("OSTIR 预测结果的上下文 SHA-256 不一致")
    if not math.isfinite(prediction.expression) or prediction.expression <= 0:
        raise ValueError("OSTIR 翻译起始率必须是正有限数")
    if not math.isfinite(prediction.d_g_total):
        raise ValueError("OSTIR 总自由能必须是有限数")
    if prediction.unintended_start_count < 0:
        raise ValueError("OSTIR 非预期起始位点数量不能为负数")
    return {
        "translation_initiation_rate": prediction.expression,
        "d_g_total": prediction.d_g_total,
        "intended_start_position": prediction.intended_start_position,
        "unintended_start_count": prediction.unintended_start_count,
        "context_sha256": prediction.context_sha256,
        "context_parameters": {
            "previous_cds_suffix_nt": RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT,
            "current_cds_prefix_nt": RBS_CONTEXT_CDS_PREFIX_NT,
        },
    }


def upload_expression_rbs(
    config: Any,
    *,
    predictor: Callable[..., RbsPrediction] = predict_rbs_context,
) -> dict[str, Any]:
    """Upload one protein-specific RBS and immediately evaluate its CDS context."""

    accession, filename = _rbs_args(config)
    target_compound_id = validate_target_compound_id(config.target_name)
    inputs_root = Path(config.inputs_dir).expanduser().resolve()
    input_path = (inputs_root / "parts" / filename).resolve()
    parts_root = (inputs_root / "parts").resolve()
    if input_path.parent != parts_root:
        raise ValueError("RBS 文件必须直接位于 inputs/parts 目录")
    if not input_path.is_file():
        raise FileNotFoundError(f"未找到 RBS 文件：{input_path}")
    part_id, sequence, input_format = _read_part(input_path, "rbs")

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
    cassette_index, cds_index = _rbs_location(cassettes, accession)
    _merge_existing_draft(
        manifest,
        cassettes,
        box_fingerprint=box_fingerprint,
        cds_fingerprint=cds_fingerprint,
    )
    ostir = _predict_uploaded_rbs(
        manifest_path=manifest_path,
        project_root=project_root,
        cassette_index=cassette_index,
        cds_index=cds_index,
        accession=accession,
        part_id=part_id,
        rbs_sequence=sequence,
        predictor=predictor,
    )
    return _commit_uploaded_part(
        manifest_path=manifest_path,
        manifest=manifest,
        current_revision=current_revision,
        project_root=project_root,
        target_compound_id=target_compound_id,
        box_fingerprint=box_fingerprint,
        cds_fingerprint=cds_fingerprint,
        cassettes=cassettes,
        cassette_index=cassette_index,
        role="rbs",
        filename=filename,
        part_id=part_id,
        sequence=sequence,
        input_format=input_format,
        accession=accession,
        extra_fields={"ostir": ostir},
    )


def upload_expression_promoter(config: Any) -> dict[str, Any]:
    """Upload or replace the promoter for one selected expression cassette."""

    return _upload_expression_boundary_part(config, "promoter")


def upload_expression_terminator(config: Any) -> dict[str, Any]:
    """Upload or replace the terminator for one selected expression cassette."""

    return _upload_expression_boundary_part(config, "terminator")


__all__ = [
    "EXPRESSION_PART_UPLOAD_DOWNSTREAM_SECTIONS",
    "EXPRESSION_PARTS_DRAFT_SCHEMA_VERSION",
    "PROMOTER_UPLOAD_DOWNSTREAM_SECTIONS",
    "upload_expression_promoter",
    "upload_expression_rbs",
    "upload_expression_terminator",
]
