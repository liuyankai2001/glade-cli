"""Validate and commit one or more expression-parts design references."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from src.expression_box.config import (
    EXPRESSION_PARTS_DESIGNS_SCHEMA_VERSION,
    EXPRESSION_SUCCESS_MIN_SCORE,
    PARTS_RECOMMENDATION_ALGORITHM_VERSION,
)
from src.expression_box.parts_manifest_adapter import load_expression_parts_context
from src.expression_box.parts_models import ExpressionPartsContext
from src.expression_box.parts_pipeline import EXPRESSION_PARTS_DESIGNS_FILENAME
from src.pathway_analyze.target_id import validate_target_compound_id
from src.write_manifest.expression_constructs import prepare_expression_constructs
from src.write_manifest.store import read_design_manifest, update_design_manifest


PARTS_SELECTION_SCHEMA_VERSION = "parts_selection.v1"
PARTS_SELECTION_DOWNSTREAM_SECTIONS = (
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "final_assembly_plan",
    "final_assembly",
)
_SELECTION_TOKEN = re.compile(r"^[0-9]+(?::[0-9]+)?$")


def _stable_json_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def parse_expression_parts_design_ids(value: Any) -> tuple[int, ...]:
    """Expand positive IDs and inclusive ``start:end`` ranges."""

    if value is None:
        raise ValueError(
            "未指定表达元件方案，请使用 write --expression-parts 1:12"
        )
    if isinstance(value, bool):
        raise ValueError("表达元件方案编号不能是布尔值")
    if isinstance(value, (str, int)):
        tokens: Sequence[Any] = (value,)
    elif isinstance(value, Sequence):
        tokens = value
    else:
        raise ValueError("表达元件方案必须是编号或 start:end 闭区间")

    expanded: list[int] = []
    seen: set[int] = set()
    for raw_token in tokens:
        if isinstance(raw_token, bool):
            raise ValueError("表达元件方案编号不能是布尔值")
        token = str(raw_token).strip()
        if not token or not _SELECTION_TOKEN.fullmatch(token):
            raise ValueError(
                f"无效的表达元件方案范围 {token!r}；请使用 3 或 1:12"
            )
        if ":" in token:
            start_text, end_text = token.split(":", 1)
            start = int(start_text)
            end = int(end_text)
            if start < 1 or end < 1:
                raise ValueError("表达元件方案编号必须是正整数")
            if start > end:
                raise ValueError(
                    f"表达元件方案范围不能倒序：{token}"
                )
            values = range(start, end + 1)
        else:
            value_id = int(token)
            if value_id < 1:
                raise ValueError("表达元件方案编号必须是正整数")
            values = (value_id,)
        for design_id in values:
            if design_id not in seen:
                seen.add(design_id)
                expanded.append(design_id)
    if not expanded:
        raise ValueError("至少需要选择一个表达元件方案")
    return tuple(expanded)


def _artifact_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "expression_box"
        / EXPRESSION_PARTS_DESIGNS_FILENAME
    )


def _read_artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到表达元件候选方案，请先运行："
            "expression --design --parts -i <输入文件>"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"表达元件候选方案不是有效 JSON：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"表达元件候选方案根节点必须是 JSON 对象：{path}")
    return payload


def _part_sequence_valid(part: Mapping[str, Any], role: str) -> bool:
    if str(part.get("role") or "").strip().lower() != role:
        return False
    sequence = str(part.get("sequence") or "").strip().upper()
    if not sequence or set(sequence) - set("ACGT"):
        return False
    expected_hash = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
    return (
        part.get("sequence_sha256") == expected_hash
        and int(part.get("length_bp") or 0) == len(sequence)
    )


def _validate_design_structure(
    design: Mapping[str, Any],
    context: ExpressionPartsContext,
    minimum_score: float,
) -> None:
    try:
        design_id = int(design.get("design_id"))
        rank = int(design.get("rank"))
        score = float(design.get("expression_success_score"))
    except (TypeError, ValueError) as exc:
        raise ValueError("表达元件方案包含无效的编号、排名或评分") from exc
    if design_id < 1 or rank < 1:
        raise ValueError("表达元件方案编号和排名必须是正整数")
    if design.get("passes_success_threshold") is not True or score < minimum_score:
        raise ValueError(f"表达元件方案 {design_id} 未达到成功评分门槛")

    raw_cassettes = design.get("cassettes")
    if not isinstance(raw_cassettes, list):
        raise ValueError(f"表达元件方案 {design_id} 缺少 cassettes")
    cassette_by_index = {
        int(item.get("cassette_index") or 0): item
        for item in raw_cassettes
        if isinstance(item, Mapping)
    }
    expected_indexes = {cassette.cassette_index for cassette in context.cassettes}
    if set(cassette_by_index) != expected_indexes or len(raw_cassettes) != len(expected_indexes):
        raise ValueError(f"表达元件方案 {design_id} 与当前表达盒分组不一致")

    signature: list[str] = []
    for cassette in context.cassettes:
        raw = cassette_by_index[cassette.cassette_index]
        promoter = raw.get("promoter")
        terminator = raw.get("terminator")
        genes = raw.get("genes")
        audit = raw.get("sequence_audit")
        if not isinstance(promoter, Mapping) or not _part_sequence_valid(promoter, "promoter"):
            raise ValueError(f"表达元件方案 {design_id} 包含无效 promoter")
        if not isinstance(terminator, Mapping) or not _part_sequence_valid(terminator, "terminator"):
            raise ValueError(f"表达元件方案 {design_id} 包含无效 terminator")
        if not isinstance(genes, list) or len(genes) != len(cassette.cds):
            raise ValueError(f"表达元件方案 {design_id} 的基因数量无效")
        if not isinstance(audit, Mapping) or audit.get("gate_status") != "PASS":
            raise ValueError(f"表达元件方案 {design_id} 未通过序列安全检查")

        signature.append(str(promoter.get("part_id") or ""))
        for raw_gene, expected_cds in zip(genes, cassette.cds, strict=True):
            if not isinstance(raw_gene, Mapping):
                raise ValueError(f"表达元件方案 {design_id} 包含无效 gene")
            if str(raw_gene.get("accession") or "") != expected_cds.accession:
                raise ValueError(f"表达元件方案 {design_id} 的 CDS 顺序已经变化")
            if raw_gene.get("cds_sequence_sha256") != expected_cds.sequence_sha256:
                raise ValueError(f"表达元件方案 {design_id} 的 CDS 已经变化")
            if int(raw_gene.get("cds_length_nt") or 0) != len(expected_cds.sequence):
                raise ValueError(f"表达元件方案 {design_id} 的 CDS 长度无效")
            rbs = raw_gene.get("rbs")
            ostir = raw_gene.get("ostir")
            if not isinstance(rbs, Mapping) or not _part_sequence_valid(rbs, "rbs"):
                raise ValueError(f"表达元件方案 {design_id} 包含无效 RBS")
            if not isinstance(ostir, Mapping) or float(
                ostir.get("translation_initiation_rate") or 0.0
            ) <= 0.0:
                raise ValueError(f"表达元件方案 {design_id} 包含无效 OSTIR 结果")
            signature.append(str(rbs.get("part_id") or ""))
        signature.append(str(terminator.get("part_id") or ""))
    if list(design.get("design_signature") or []) != signature:
        raise ValueError(f"表达元件方案 {design_id} 的 design_signature 无效")


def _validate_artifact(
    artifact: Mapping[str, Any],
    *,
    context: ExpressionPartsContext,
) -> list[dict[str, Any]]:
    if artifact.get("schema_version") != EXPRESSION_PARTS_DESIGNS_SCHEMA_VERSION:
        raise ValueError(
            "不支持的表达元件候选方案 schema_version，请重新运行 "
            "expression --design --parts"
        )
    if artifact.get("algorithm_version") != PARTS_RECOMMENDATION_ALGORITHM_VERSION:
        raise ValueError(
            "表达元件推荐算法已经变化，请重新运行 expression --design --parts"
        )
    if artifact.get("status") not in {"complete", "partial"}:
        raise ValueError("表达元件候选方案没有可写入的合格设计")
    if str(artifact.get("target_compound_id") or "") != context.target_compound_id:
        raise ValueError("表达元件候选方案与当前目标化合物不一致")
    source = artifact.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("表达元件候选方案缺少 source")
    if source.get("context_input_fingerprint") != context.input_fingerprint:
        raise ValueError(
            "表达元件候选方案与当前 CDS 或表达盒不一致，请重新运行 "
            "expression --design --parts"
        )
    if (
        source.get("expression_box_selection_fingerprint")
        != context.expression_box_selection_fingerprint
        or source.get("cds_selection_source_fingerprint")
        != context.cds_selection_source_fingerprint
    ):
        raise ValueError("表达元件候选方案的上游来源已经变化")
    ranking = artifact.get("ranking")
    if not isinstance(ranking, Mapping):
        raise ValueError("表达元件候选方案缺少 ranking")
    minimum_score = float(ranking.get("minimum_success_score") or 0.0)
    if minimum_score != EXPRESSION_SUCCESS_MIN_SCORE:
        raise ValueError("表达元件候选方案的成功评分门槛已经变化")
    raw_designs = artifact.get("designs")
    if not isinstance(raw_designs, list) or not raw_designs:
        raise ValueError("表达元件候选方案 designs 为空")
    if int(artifact.get("design_count") or 0) != len(raw_designs):
        raise ValueError("表达元件候选方案 design_count 无效")

    designs: list[dict[str, Any]] = []
    design_ids: set[int] = set()
    ranks: set[int] = set()
    for raw_design in raw_designs:
        if not isinstance(raw_design, dict):
            raise ValueError("表达元件候选方案包含无效 design")
        _validate_design_structure(raw_design, context, minimum_score)
        design_id = int(raw_design["design_id"])
        rank = int(raw_design["rank"])
        if design_id in design_ids or rank in ranks:
            raise ValueError("表达元件候选方案包含重复的 design_id 或 rank")
        design_ids.add(design_id)
        ranks.add(rank)
        designs.append(raw_design)
    return designs


def _select_designs(
    designs: list[dict[str, Any]],
    requested_ids: tuple[int, ...],
) -> list[dict[str, Any]]:
    by_id = {int(item["design_id"]): item for item in designs}
    missing = sorted(set(requested_ids) - set(by_id))
    if missing:
        raise ValueError(
            f"不存在表达元件方案 {missing}；可用方案编号：{sorted(by_id)}"
        )
    return sorted((by_id[item] for item in set(requested_ids)), key=lambda item: int(item["rank"]))


def _selection_payload(
    selected: list[dict[str, Any]],
    artifact: Mapping[str, Any],
    context: ExpressionPartsContext,
) -> dict[str, Any]:
    source = artifact["source"]
    references = [
        {
            "design_id": int(design["design_id"]),
            "rank": int(design["rank"]),
            "expression_success_score": float(design["expression_success_score"]),
            "expression_regime": str(design["expression_regime"]),
            "system_recommended": bool(design.get("recommended")),
            "design_fingerprint": _stable_json_hash(design),
        }
        for design in selected
    ]
    selected_content_fingerprint = _stable_json_hash(
        {
            "schema_version": artifact["schema_version"],
            "algorithm_version": artifact["algorithm_version"],
            "artifact_input_fingerprint": source["input_fingerprint"],
            "designs": selected,
        }
    )
    selection_fingerprint = _stable_json_hash(
        {
            "context_input_fingerprint": context.input_fingerprint,
            "selected_content_fingerprint": selected_content_fingerprint,
            "design_references": references,
        }
    )
    scores = [item["expression_success_score"] for item in references]
    warnings = list(
        dict.fromkeys(
            str(warning)
            for design in selected
            for warning in design.get("warnings", [])
            if str(warning)
        )
    )
    return {
        "schema_version": PARTS_SELECTION_SCHEMA_VERSION,
        "status": "selected",
        "selection_status": "user_selected",
        "selected_design_ids": [item["design_id"] for item in references],
        "primary_design_id": references[0]["design_id"],
        "design_count": len(references),
        "selection_fingerprint": selection_fingerprint,
        "design_references": references,
        "source": {
            "artifact": "expression_box/expression_parts_designs.json",
            "designs_schema_version": artifact["schema_version"],
            "algorithm_version": artifact["algorithm_version"],
            "artifact_status": artifact["status"],
            "artifact_manifest_revision": source.get("manifest_revision"),
            "artifact_input_fingerprint": source["input_fingerprint"],
            "selected_content_fingerprint": selected_content_fingerprint,
            "expression_box_selection_fingerprint": (
                context.expression_box_selection_fingerprint
            ),
            "cds_selection_source_fingerprint": (
                context.cds_selection_source_fingerprint
            ),
            "context_input_fingerprint": context.input_fingerprint,
            "candidate_snapshot_fingerprint": source.get(
                "candidate_snapshot_fingerprint"
            ),
        },
        "summary": {
            "highest_score": max(scores),
            "lowest_score": min(scores),
            "minimum_success_score": EXPRESSION_SUCCESS_MIN_SCORE,
        },
        "warnings": warnings,
    }


def write_expression_parts_selection(config: Any) -> dict[str, Any]:
    """Commit selected designs and their complete concatenated GenBank files."""

    requested_ids = parse_expression_parts_design_ids(
        getattr(config, "expression_parts", None)
    )
    target_compound_id = validate_target_compound_id(config.target_name)
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if not recorded_target:
        raise ValueError("manifest 尚未写入表达盒方案")
    if recorded_target != target_compound_id:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target}，"
            f"与当前输入目标 {target_compound_id} 不一致"
        )

    context = load_expression_parts_context(
        manifest_path,
        Path(config.project_output_path).expanduser().resolve(),
    )
    manifest = read_design_manifest(manifest_path)
    if int(manifest.get("revision", 0)) != context.manifest_revision:
        raise ValueError("manifest revision 在读取表达元件上下文时发生变化，请重试")
    artifact = _read_artifact(_artifact_path(config))
    designs = _validate_artifact(artifact, context=context)
    selected = _select_designs(designs, requested_ids)
    payload = _selection_payload(selected, artifact, context)

    current_selection = manifest.get("parts_selection")
    current_constructs = manifest.get("assembled_expression_constructs")
    transaction = prepare_expression_constructs(
        selected_designs=selected,
        selection_payload=payload,
        context=context,
        current_section=current_constructs,
    )
    manifest_changed = not (
        isinstance(current_selection, Mapping)
        and dict(current_selection) == payload
        and isinstance(current_constructs, Mapping)
        and dict(current_constructs) == transaction.section
    )
    try:
        transaction.install()
        if manifest_changed:
            updated_manifest = update_design_manifest(
                manifest_path,
                target_compound_id=target_compound_id,
                sections={
                    "parts_selection": payload,
                    "assembled_expression_constructs": transaction.section,
                },
                discard_sections=PARTS_SELECTION_DOWNSTREAM_SECTIONS,
                expected_revision=context.manifest_revision,
            )
        else:
            updated_manifest = manifest
    except Exception:
        try:
            transaction.rollback()
        except Exception as rollback_error:
            raise RuntimeError(
                "表达构建文件写入失败，且自动回滚未能完成；"
                "请检查 expression_constructs 目录"
            ) from rollback_error
        raise
    cleanup_warning = transaction.finalize()

    warnings = list(payload["warnings"])
    if cleanup_warning:
        warnings.append(cleanup_warning)

    construct_dir = (
        Path(config.project_output_path).expanduser().resolve()
        / "expression_constructs"
    )

    return {
        "运行成功": True,
        "目标化合物": target_compound_id,
        "表达元件方案编号": payload["selected_design_ids"],
        "主方案编号": payload["primary_design_id"],
        "方案数量": payload["design_count"],
        "完整表达构建数量": transaction.section["design_count"],
        "最高成功评分": payload["summary"]["highest_score"],
        "最低成功评分": payload["summary"]["lowest_score"],
        "GenBank目录": str(construct_dir),
        "GenBank是否复用": not transaction.needs_install,
        "GenBank是否修复": transaction.is_repair,
        "清单是否更新": manifest_changed,
        "清单文件": str(manifest_path),
        "清单版本": updated_manifest["revision"],
        "警告": warnings,
    }


def run_write_expression_parts_selection(config: Any) -> dict[str, Any]:
    """CLI entry point for ``write --expression-parts``."""

    result = write_expression_parts_selection(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "PARTS_SELECTION_DOWNSTREAM_SECTIONS",
    "PARTS_SELECTION_SCHEMA_VERSION",
    "parse_expression_parts_design_ids",
    "run_write_expression_parts_selection",
    "write_expression_parts_selection",
]
