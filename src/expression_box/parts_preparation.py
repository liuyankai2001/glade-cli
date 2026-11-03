"""Common, file-backed parts preparation for uploads and recommendations."""

from __future__ import annotations

import copy
import uuid
from collections.abc import Mapping
from io import StringIO
from pathlib import Path
from typing import Any

from Bio import SeqIO


PARTS_DRAFT_SCHEMA = "expression_parts_draft.v3"
PREVIOUS_PARTS_DRAFT_SCHEMA = "expression_parts_draft.v2"
LEGACY_PARTS_DRAFT_SCHEMA = "expression_parts_draft.v1"
SELECTED_PARTS_SCHEMA = "expression_selected_parts.v2"
SELECTED_PARTS_PATH = "expression_box/selected_expression_parts.json"


def draft_content(draft: Mapping[str, Any]) -> dict[str, Any]:
    keys = ("target_compound_id", "source", "source_type", "designs")
    if draft.get("schema_version") == LEGACY_PARTS_DRAFT_SCHEMA:
        keys = ("target_compound_id", "source", "cassettes")
    return {key: draft.get(key) for key in keys}


def make_parts_draft(*, target: str, source: Mapping[str, Any],
                     source_type: str, designs: list[dict[str, Any]]) -> dict[str, Any]:
    cassettes = [cassette for design in designs for cassette in design["cassettes"]]
    promoter_count = sum(c.get("promoter") is not None for c in cassettes)
    terminator_count = sum(c.get("terminator") is not None for c in cassettes)
    rbs_count = sum(len(c.get("rbs_by_accession", {})) for c in cassettes)
    ready = bool(cassettes) and all(
        c.get("promoter") is not None and c.get("terminator") is not None
        and set(c.get("rbs_by_accession", {})) == set(c["protein_accessions"])
        for c in cassettes
    )
    content = {"target_compound_id": target, "source": dict(source),
               "source_type": source_type, "designs": designs}
    return {
        "schema_version": PARTS_DRAFT_SCHEMA,
        "status": "ready" if ready else "partial",
        **content,
        "summary": {"design_count": len(designs), "cassette_count": len(cassettes),
                    "promoter_count": promoter_count,
                    "required_promoter_count": len(cassettes),
                    "rbs_count": rbs_count, "terminator_count": terminator_count},
    }


def _without_validation_digests(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _without_validation_digests(item) for key, item in value.items()
                if not str(key).endswith(("sha256", "fingerprint"))}
    if isinstance(value, list):
        return [_without_validation_digests(item) for item in value]
    return value


def normalize_parts_draft(draft: Mapping[str, Any]) -> dict[str, Any]:
    schema = draft.get("schema_version")
    if schema not in {PARTS_DRAFT_SCHEMA, PREVIOUS_PARTS_DRAFT_SCHEMA, LEGACY_PARTS_DRAFT_SCHEMA}:
        raise ValueError("不支持的 expression_parts_draft schema_version")
    if schema == PREVIOUS_PARTS_DRAFT_SCHEMA and draft.get("source_type") == "recommended":
        raise ValueError("旧版推荐元件草稿已过期，请重新计算并确认表达元件方案")
    source = draft.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("元件准备记录缺少 source")
    if schema == LEGACY_PARTS_DRAFT_SCHEMA:
        designs = [{"design_id": 1, "rank": 1, "name": "用户上传元件方案",
                    "cassettes": copy.deepcopy(draft.get("cassettes"))}]
        source_type = "user_uploaded"
    else:
        designs = copy.deepcopy(draft.get("designs"))
        source_type = draft.get("source_type")
    if source_type not in {"recommended", "user_uploaded"}:
        raise ValueError("元件准备记录的来源无效")
    if schema != PARTS_DRAFT_SCHEMA and source_type != "user_uploaded":
        raise ValueError("旧版推荐元件草稿已过期，请重新计算并确认表达元件方案")
    if not isinstance(designs, list) or not designs:
        raise ValueError("元件准备记录没有方案")
    ids, ranks = set(), set()
    for design in designs:
        if not isinstance(design, dict):
            raise ValueError("元件准备方案无效")
        design_id, rank = design.get("design_id"), design.get("rank")
        if (type(design_id) is not int or design_id < 1 or design_id in ids
                or type(rank) is not int or rank < 1 or rank in ranks):
            raise ValueError("元件准备方案编号或顺序无效")
        ids.add(design_id)
        ranks.add(rank)
        if not isinstance(design.get("cassettes"), list) or not design["cassettes"]:
            raise ValueError("元件准备方案没有表达盒")
        for index, cassette in enumerate(design["cassettes"], start=1):
            if not isinstance(cassette, dict) or not isinstance(cassette.get("protein_accessions"), list):
                raise ValueError("元件准备方案的表达盒无效")
            accessions = cassette["protein_accessions"]
            if (cassette.get("cassette_index") != index or not accessions
                    or any(not isinstance(a, str) or not a for a in accessions)
                    or len(set(accessions)) != len(accessions)):
                raise ValueError("元件准备方案的表达盒分组或顺序无效")
            rbs = cassette.get("rbs_by_accession", {})
            if not isinstance(rbs, Mapping) or set(rbs) - set(cassette["protein_accessions"]):
                raise ValueError("元件准备方案包含未知 RBS 蛋白")
    if source_type == "user_uploaded" and len(designs) != 1:
        raise ValueError("上传元件只能对应一个当前方案")
    if schema != PARTS_DRAFT_SCHEMA:
        source = _without_validation_digests(source)
        designs = _without_validation_digests(designs)
        stale_parts_designs(designs)
    designs.sort(key=lambda d: d["rank"])
    return make_parts_draft(target=draft["target_compound_id"], source=source,
                            source_type=source_type, designs=designs)


def project_file(root: Path, value: Any) -> Path:
    path = (root / str(value or "")).resolve()
    if not str(value or "") or path == root.resolve() or root.resolve() not in path.parents:
        raise ValueError("表达元件文件必须位于当前项目输出目录")
    return path


def read_part_sequence(part: Mapping[str, Any], root: Path,
                       guarded: dict[Path, bytes] | None = None) -> str:
    reference = part.get("sequence_file")
    if not isinstance(reference, Mapping) or reference.get("format") != "fasta":
        raise ValueError("表达元件缺少有效的 FASTA 快照引用")
    path = project_file(root, reference.get("path"))
    content = path.read_bytes()
    try:
        record = SeqIO.read(StringIO(content.decode("utf-8")), "fasta")
    except Exception as exc:
        raise ValueError(f"表达元件快照不是有效 FASTA：{path}") from exc
    sequence = str(record.seq).upper()
    if (not sequence or set(sequence) - set("ACGT") or record.id != part.get("part_id")
            or len(sequence) != part.get("length_bp")):
        raise ValueError(f"表达元件序列与记录不一致：{path}")
    if guarded is not None:
        if path in guarded and guarded[path] != content:
            raise ValueError(f"表达元件快照在读取时发生变化：{path}")
        guarded[path] = content
    return sequence


def snapshot_part(part: Mapping[str, Any], root: Path) -> tuple[dict[str, Any], Path, bytes]:
    """Prepare an immutable snapshot without writing anything."""
    result = copy.deepcopy(dict(part))
    sequence = str(result.pop("sequence")).upper()
    part_id = str(result["part_id"])
    content = (f">{part_id}\n{sequence}\n").encode("utf-8")
    relative = f"expression_box/prepared_parts/{result['role']}/{uuid.uuid4().hex}.fasta"
    result["sequence_file"] = {"path": relative, "format": "fasta"}
    return result, project_file(root, relative), content


def bind_cds(cassettes: list[dict[str, Any]], proteins: Mapping[str, Any]) -> None:
    for cassette in cassettes:
        cassette["optimized_cds_by_accession"] = {
            accession: copy.deepcopy(proteins[accession]["optimized_cds"])
            for accession in cassette["protein_accessions"]
        }


def stale_parts_designs(designs: list[dict[str, Any]]) -> None:
    for design in designs:
        recommendation = design.get("recommendation")
        if isinstance(recommendation, dict):
            recommendation["status"] = "stale"
        for cassette in design["cassettes"]:
            for key, status in (("restriction_site_audit", "restriction_audit_status"),
                                ("homopolymer_audit", "homopolymer_audit_status")):
                if cassette.pop(key, None) is not None:
                    cassette[status] = "stale"
            cassette.pop("sequence_audit", None)
            for part in cassette.get("rbs_by_accession", {}).values():
                if not isinstance(part, dict):
                    raise ValueError("RBS 元件记录无效")
                part.pop("ostir", None)
                part["ostir_status"] = "stale"
