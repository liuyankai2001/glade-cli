"""Validate and publish a matching protein FASTA/UniProt metadata pair."""
from __future__ import annotations

import json
import os
from datetime import UTC, datetime
from functools import cache
from pathlib import Path
from typing import Any

from filelock import FileLock

from src.main_protein_selection.sequence_quality import analyze_protein_sequence
from src.main_protein_selection.uniprot_protein_candidates import (
    _entry_is_reviewed,
    extract_cofactors,
    get_protein_name,
    get_sequence_value,
    get_sequence_version,
)
from src.protein_to_cds.artifacts import ArtifactTransaction, write_bytes_atomic
from src.protein_to_cds.search_protein_sequence import (
    ProteinSequenceError,
    _parse_single_fasta,
    normalize_uniprot_accession,
)

METADATA_SCHEMA = "protein_sequence_metadata.v1"


@cache
def _protein_lock(path: str) -> FileLock:
    return FileLock(path, timeout=30)


def protein_file_lock(folder: Path, accession: str) -> FileLock:
    return _protein_lock(os.path.normcase(str(folder.resolve() / f".{accession}.download.lock")))


def accession_value(value: str) -> str:
    try:
        return normalize_uniprot_accession(value)
    except ProteinSequenceError as exc:
        raise ValueError(str(exc)) from exc


def entry_sequence(entry: dict[str, Any]) -> str:
    quality = analyze_protein_sequence(get_sequence_value(entry))
    if not quality.normalized_sequence:
        raise ValueError("UniProt 氨基酸序列为空")
    if not quality.is_constructable:
        raise ValueError(f"UniProt 氨基酸序列含不支持的字符：{quality.rejection_reason()}")
    declared = entry.get("sequence", {})
    if isinstance(declared, dict) and declared.get("length") is not None:
        if declared["length"] != len(quality.normalized_sequence):
            raise ValueError("UniProt 氨基酸序列与声明长度不一致")
    return quality.normalized_sequence


def entry_metadata(entry: dict[str, Any]) -> dict[str, Any]:
    sequence = entry_sequence(entry)
    accession = accession_value(entry.get("primaryAccession", ""))
    organism = entry.get("organism") or {}
    if not isinstance(organism, dict):
        raise ValueError("UniProt 来源物种信息格式错误")
    audit = entry.get("entryAudit") or {}
    version = audit.get("sequenceVersion") if isinstance(audit, dict) else None
    if version is None:
        version = get_sequence_version(entry)
    elif isinstance(version, bool) or not isinstance(version, int) or version < 1:
        raise ValueError("UniProt 序列版本格式错误")
    return {
        "schema_version": METADATA_SCHEMA,
        "accession": accession,
        "protein_name": get_protein_name(entry),
        "organism_name": organism.get("scientificName", ""),
        "organism_id": organism.get("taxonId"),
        "reviewed": _entry_is_reviewed(entry),
        "cofactors": extract_cofactors(entry),
        "length_aa": len(sequence),
        "sequence_version": version,
        "source_url": f"https://rest.uniprot.org/uniprotkb/{accession}.json",
        "downloaded_at": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "uniprot_entry": entry,
    }


def _paths(folder: Path, accession: str) -> tuple[Path, Path]:
    return folder / f"{accession}.fasta", folder / f"{accession}.json"


def _read_fasta(path: Path, accession: str) -> str:
    primary, _, sequence = _parse_single_fasta(path.read_text(encoding="utf-8"))
    if primary != accession:
        raise ValueError("FASTA 中的编号与文件名不一致")
    return sequence


def _read_metadata(path: Path, accession: str) -> dict[str, Any]:
    metadata = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(metadata, dict) or not isinstance(metadata.get("uniprot_entry"), dict):
        raise ValueError("缺少 UniProt 原始记录")
    expected = entry_metadata(metadata["uniprot_entry"])
    if expected["accession"] != accession:
        raise ValueError("JSON 中的编号与文件名不一致")
    for key, value in expected.items():
        if key != "downloaded_at" and metadata.get(key) != value:
            raise ValueError(f"JSON 元数据字段不一致：{key}")
    if not isinstance(metadata.get("downloaded_at"), str) or not metadata["downloaded_at"]:
        raise ValueError("JSON 缺少下载时间")
    return metadata


def _inspect(folder: Path, accession: str) -> tuple[str | None, dict[str, Any] | None]:
    fasta, metadata_file = _paths(folder, accession)
    try:
        if metadata_file.exists():
            raw = json.loads(metadata_file.read_text(encoding="utf-8-sig"))
            from .local_sequences import LOCAL_METADATA_SCHEMA, read_sequence_pair
            if isinstance(raw, dict) and raw.get("schema_version") == LOCAL_METADATA_SCHEMA:
                metadata, sequence = read_sequence_pair(fasta, accession)
                return sequence, metadata
        sequence = _read_fasta(fasta, accession) if fasta.exists() else None
        metadata = _read_metadata(metadata_file, accession) if metadata_file.exists() else None
        if sequence is not None and metadata is not None:
            if sequence != entry_sequence(metadata["uniprot_entry"]):
                raise ValueError("FASTA 与 JSON 氨基酸序列不一致")
        return sequence, metadata
    except (OSError, ValueError, UnicodeError, ProteinSequenceError, TypeError, AttributeError) as exc:
        raise ValueError(f"本地蛋白文件损坏或内容冲突，请使用 --force 刷新：{accession}（{exc}）") from exc


def read_cached_pair(folder: Path, accession: str) -> dict[str, Any] | None:
    if not folder.exists():
        return None
    with protein_file_lock(folder, accession):
        sequence, metadata = _inspect(folder, accession)
        return metadata if sequence is not None and metadata is not None and metadata.get("schema_version") == METADATA_SCHEMA else None


def publish_entry(folder: Path, entry: dict[str, Any], *, force: bool) -> tuple[dict[str, Any], bool]:
    metadata = entry_metadata(entry)
    accession, sequence = metadata["accession"], entry_sequence(entry)
    folder.mkdir(parents=True, exist_ok=True)
    with protein_file_lock(folder, accession):
        if not force:
            previous_sequence, previous_metadata = _inspect(folder, accession)
            if previous_sequence is not None and previous_metadata is not None and previous_metadata.get("schema_version") == METADATA_SCHEMA:
                return previous_metadata, True
            existing = previous_sequence
            if previous_metadata is not None and existing is None:
                existing = entry_sequence(previous_metadata["uniprot_entry"])
            if existing is not None and existing != sequence:
                raise ValueError(f"本地蛋白序列与 UniProt 内容冲突，请使用 --force 刷新：{accession}")
        fasta, metadata_file = _paths(folder, accession)
        label = " ".join(str(entry.get("uniProtkbId") or metadata["protein_name"]).split())
        fasta_bytes = f">{accession} {label}\n{sequence}\n".encode("utf-8")
        metadata_bytes = (json.dumps(metadata, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
        with ArtifactTransaction((fasta, metadata_file)) as transaction:
            write_bytes_atomic(fasta, fasta_bytes)
            write_bytes_atomic(metadata_file, metadata_bytes)
            transaction.commit()
    return metadata, False
