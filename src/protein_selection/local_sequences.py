"""Offline FASTA records with either UniProt or explicit local annotations."""
from __future__ import annotations

import hashlib
import json
import re
from io import StringIO
from pathlib import Path, PureWindowsPath

from Bio import SeqIO

from src.main_protein_selection.sequence_quality import analyze_protein_sequence
from src.protein_to_cds.search_protein_sequence import ProteinSequenceError
from .sequence_store import METADATA_SCHEMA, _read_fasta, _read_metadata, entry_sequence, protein_file_lock

LOCAL_METADATA_SCHEMA = "protein_sequence_metadata.local.v1"
UNKNOWN_FIELDS = ("protein_name", "organism_name", "organism_id", "reviewed", "cofactors",
                  "sequence_version", "source_url")


def local_filename(value: str) -> str:
    if (not isinstance(value, str) or not value or value != value.strip()
            or "/" in value or "\\" in value or PureWindowsPath(value).drive
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*\.fasta", value, re.IGNORECASE)):
        raise ValueError("-f 只接受 protein_sequences 下的 .fasta 文件名，不接受路径")
    return value


def sequence_identity(path: Path) -> str:
    return Path(local_filename(path.name)).stem.upper()


def read_fasta_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8-sig")
    if not text.lstrip().startswith(">"):
        raise ValueError("蛋白文件不是有效 FASTA")
    try:
        records = list(SeqIO.parse(StringIO(text), "fasta"))
    except Exception as exc:
        raise ValueError("无法解析蛋白 FASTA") from exc
    if len(records) != 1:
        raise ValueError("蛋白 FASTA 必须且只能包含一条序列")
    quality = analyze_protein_sequence(str(records[0].seq))
    if not quality.is_constructable:
        raise ValueError(f"蛋白序列为空或含不支持的氨基酸：{quality.rejection_reason()}")
    return quality.normalized_sequence


def local_metadata(accession: str, sequence: str) -> dict:
    return {"schema_version": LOCAL_METADATA_SCHEMA, "source_type": "local", "accession": accession,
            "length_aa": len(sequence), "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            **dict.fromkeys(UNKNOWN_FIELDS)}


def read_sequence_pair(path: Path, accession: str, *, allow_missing_metadata: bool = False) -> tuple[dict, str]:
    """Missing sidecars are accepted only during an explicit add, never on reads."""
    path = path.expanduser().resolve()
    if sequence_identity(path) != accession:
        raise ValueError("选中蛋白序列文件名与编号不一致")
    with protein_file_lock(path.parent, accession):
        if not path.is_file():
            raise ValueError(f"缺少完整的蛋白 FASTA/JSON：{accession}")
        sequence = read_fasta_sequence(path)
        sidecar = path.with_suffix(".json")
        if not sidecar.exists():
            if not allow_missing_metadata:
                raise ValueError(f"缺少完整的蛋白 FASTA/JSON：{accession}")
            return local_metadata(accession, sequence), sequence
        try:
            metadata = json.loads(sidecar.read_text(encoding="utf-8-sig"))
            if not isinstance(metadata, dict):
                raise ValueError("蛋白 JSON 必须是对象")
            if metadata.get("schema_version") == METADATA_SCHEMA:
                metadata = _read_metadata(sidecar, accession)
                if _read_fasta(path, accession) != entry_sequence(metadata["uniprot_entry"]):
                    raise ValueError("FASTA 与 JSON 氨基酸序列不一致")
            elif metadata.get("schema_version") == LOCAL_METADATA_SCHEMA:
                expected = local_metadata(accession, sequence)
                if any(metadata.get(key) != value or type(metadata.get(key)) is not type(value)
                       for key, value in expected.items()):
                    raise ValueError("本地蛋白 JSON 与序列、编号或未知注释不一致")
                if set(expected) != set(metadata):
                    raise ValueError("本地蛋白 JSON 字段不完整或含不支持字段")
            else:
                raise ValueError("不支持的蛋白注释 JSON 格式")
        except (ValueError, TypeError, AttributeError, ProteinSequenceError) as exc:
            raise ValueError(f"本地蛋白 JSON 损坏或内容冲突：{accession}（{exc}）") from exc
        return metadata, sequence


def runtime_annotations(metadata: dict) -> dict:
    return {**metadata, "protein_name": metadata.get("protein_name") or "",
            "organism_name": metadata.get("organism_name") or "", "cofactors": metadata.get("cofactors") or []}
