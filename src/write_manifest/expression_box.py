"""Write one user-selected expression-box grouping to the design manifest."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.expression_box.config import (
    EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION,
    GROUPING_ALGORITHM_VERSION,
)
from src.expression_box.manifest_adapter import load_expression_grouping_context
from src.expression_box.models import ExpressionGroupingDesign, ExpressionProtein
from src.expression_box.pipeline import EXPRESSION_BOX_DESIGNS_FILENAME
from src.expression_box.protein_grouping import generate_grouping_designs
from src.pathway_analyze.target_id import validate_target_compound_id
from src.write_manifest.store import read_design_manifest, update_design_manifest


EXPRESSION_BOX_SELECTION_SCHEMA_VERSION = "expression_box_selection.v1"
EXPRESSION_BOX_SELECTION_DOWNSTREAM_SECTIONS = (
    "expression_cassette_assembly",
    "parts_selection",
    "assembled_expression_cassettes",
)


def _stable_json_hash(payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _positive_design_id(config: Any) -> int:
    raw_value = getattr(config, "expression_box", None)
    if raw_value is None:
        raise ValueError(
            "未指定表达盒方案，请使用 write --expression-box N"
        )
    if isinstance(raw_value, bool):
        raise ValueError("表达盒方案编号必须是正整数")
    try:
        design_id = int(raw_value)
    except (TypeError, ValueError) as exc:
        raise ValueError("表达盒方案编号必须是正整数") from exc
    if design_id < 1:
        raise ValueError("表达盒方案编号必须是正整数")
    return design_id


def _designs_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "expression_box"
        / EXPRESSION_BOX_DESIGNS_FILENAME
    )


def _read_designs(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到表达盒候选方案，请先运行："
            "expression --design --box -i <输入文件>"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"表达盒候选方案不是有效 JSON：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"表达盒候选方案根节点必须是 JSON 对象：{path}")
    return payload


def _protein_payload(protein: ExpressionProtein) -> dict[str, Any]:
    return {
        "accession": protein.accession,
        "roles": list(protein.roles),
        "assigned_step_indexes": list(protein.assigned_step_indexes),
        "required_by_main_accessions": list(protein.required_by_main_accessions),
        "optimized_cds_length_nt": protein.optimized_cds_length_nt,
        "optimized_cds_sequence_sha256": protein.optimized_cds_sequence_sha256,
    }


def _design_payload(
    design: ExpressionGroupingDesign,
    design_id: int,
) -> dict[str, Any]:
    cassettes: list[dict[str, Any]] = []
    for cassette_index, cassette in enumerate(design.cassettes, start=1):
        cassettes.append(
            {
                "cassette_index": cassette_index,
                "protein_accessions": [
                    protein.accession for protein in cassette.proteins
                ],
                "protein_count": len(cassette.proteins),
                "main_enzyme_count": cassette.main_enzyme_count,
                "total_cds_length_nt": cassette.total_cds_length_nt,
                "reason": cassette.reason,
                "proteins": [
                    _protein_payload(protein) for protein in cassette.proteins
                ],
            }
        )
    return {
        "design_id": design_id,
        "rank": design_id,
        "strategy": design.strategy,
        "name": design.name,
        "recommended": design.recommended,
        "cassette_count": len(cassettes),
        "protein_count": sum(item["protein_count"] for item in cassettes),
        "total_cds_length_nt": sum(
            item["total_cds_length_nt"] for item in cassettes
        ),
        "cassettes": cassettes,
        "warnings": list(design.warnings),
    }


def _validate_artifact(
    artifact: Mapping[str, Any],
    *,
    target_compound_id: str,
    input_fingerprint: str,
    cds_selection_source_fingerprint: str,
    proteins: tuple[ExpressionProtein, ...],
) -> list[dict[str, Any]]:
    if artifact.get("schema_version") != EXPRESSION_BOX_DESIGNS_SCHEMA_VERSION:
        raise ValueError(
            "不支持的表达盒候选方案 schema_version，请重新运行 "
            "expression --design --box"
        )
    if artifact.get("status") != "complete":
        raise ValueError("表达盒候选方案尚未完成，不能写入 manifest")
    if artifact.get("algorithm_version") != GROUPING_ALGORITHM_VERSION:
        raise ValueError(
            "表达盒候选方案算法版本已经变化，请重新运行 "
            "expression --design --box"
        )
    if str(artifact.get("target_compound_id") or "").strip() != target_compound_id:
        raise ValueError("表达盒候选方案与当前目标化合物不一致")

    source = artifact.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("表达盒候选方案缺少 source")
    if source.get("input_fingerprint") != input_fingerprint:
        raise ValueError(
            "表达盒候选方案与当前 CDS 不一致，请重新运行 "
            "expression --design --box"
        )
    if (
        source.get("cds_selection_source_fingerprint")
        != cds_selection_source_fingerprint
    ):
        raise ValueError(
            "表达盒候选方案的 CDS 来源已经变化，请重新运行 "
            "expression --design --box"
        )

    generation = generate_grouping_designs(proteins)
    expected_designs = [
        _design_payload(design, design_id)
        for design_id, design in enumerate(generation.designs, start=1)
    ]
    raw_designs = artifact.get("designs")
    if raw_designs != expected_designs:
        raise ValueError(
            "表达盒候选方案内容与当前分组算法不一致，结果可能已损坏；"
            "请重新运行 expression --design --box"
        )
    if artifact.get("design_count") != len(expected_designs):
        raise ValueError("表达盒候选方案的 design_count 无效")
    if artifact.get("protein_count") != len(proteins):
        raise ValueError("表达盒候选方案的 protein_count 无效")
    return expected_designs


def _select_design(
    designs: list[dict[str, Any]],
    design_id: int,
) -> dict[str, Any]:
    selected = next(
        (item for item in designs if item["design_id"] == design_id),
        None,
    )
    if selected is not None:
        return selected
    available = [item["design_id"] for item in designs]
    raise ValueError(
        f"不存在表达盒方案 {design_id}；可用方案编号：{available}"
    )


def _selection_payload(
    selected: Mapping[str, Any],
    artifact: Mapping[str, Any],
) -> dict[str, Any]:
    source = artifact["source"]
    design_fingerprint = _stable_json_hash(
        {
            "input_fingerprint": source["input_fingerprint"],
            "design": selected,
        }
    )
    cassettes = [
        {
            "cassette_index": cassette["cassette_index"],
            "protein_accessions": list(cassette["protein_accessions"]),
            "protein_count": cassette["protein_count"],
            "main_enzyme_count": cassette["main_enzyme_count"],
            "total_cds_length_nt": cassette["total_cds_length_nt"],
            "reason": cassette["reason"],
        }
        for cassette in selected["cassettes"]
    ]
    return {
        "schema_version": EXPRESSION_BOX_SELECTION_SCHEMA_VERSION,
        "selection_status": "user_selected",
        "selected_design_id": selected["design_id"],
        "selected_design_fingerprint": design_fingerprint,
        "rank": selected["rank"],
        "strategy": selected["strategy"],
        "name": selected["name"],
        "recommended": selected["recommended"],
        "source": {
            "artifact": "expression_box/expression_box_designs.json",
            "designs_schema_version": artifact["schema_version"],
            "algorithm_version": artifact["algorithm_version"],
            "artifact_manifest_revision": source.get("manifest_revision"),
            "cds_selection_source_fingerprint": (
                source["cds_selection_source_fingerprint"]
            ),
            "input_fingerprint": source["input_fingerprint"],
        },
        "summary": {
            "cassette_count": selected["cassette_count"],
            "protein_count": selected["protein_count"],
            "total_cds_length_nt": selected["total_cds_length_nt"],
        },
        "cassettes": cassettes,
        "warnings": list(selected["warnings"]),
    }


def write_expression_box_selection(config: Any) -> dict[str, Any]:
    """Validate and commit one expression-box grouping to the manifest."""

    design_id = _positive_design_id(config)
    target_compound_id = validate_target_compound_id(config.target_name)
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if not recorded_target:
        raise ValueError("manifest 尚未建立，请先完成 CDS 优化")
    if recorded_target != target_compound_id:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target}，"
            f"与当前输入目标 {target_compound_id} 不一致"
        )

    context = load_expression_grouping_context(manifest_path)
    artifact_path = _designs_path(config)
    artifact = _read_designs(artifact_path)
    designs = _validate_artifact(
        artifact,
        target_compound_id=target_compound_id,
        input_fingerprint=context.input_fingerprint,
        cds_selection_source_fingerprint=(
            context.cds_selection_source_fingerprint
        ),
        proteins=context.proteins,
    )
    selected = _select_design(designs, design_id)
    payload = _selection_payload(selected, artifact)

    current = manifest.get("expression_box_selection")
    unchanged = (
        isinstance(current, Mapping)
        and current.get("schema_version")
        == EXPRESSION_BOX_SELECTION_SCHEMA_VERSION
        and current.get("selected_design_fingerprint")
        == payload["selected_design_fingerprint"]
    )
    if unchanged:
        updated_manifest = manifest
    else:
        updated_manifest = update_design_manifest(
            manifest_path,
            target_compound_id=target_compound_id,
            sections={"expression_box_selection": payload},
            discard_sections=EXPRESSION_BOX_SELECTION_DOWNSTREAM_SECTIONS,
            expected_revision=context.manifest_revision,
        )

    return {
        "运行成功": True,
        "目标化合物": target_compound_id,
        "表达盒方案编号": selected["design_id"],
        "方案名称": selected["name"],
        "方案策略": selected["strategy"],
        "是否系统推荐": selected["recommended"],
        "表达盒数量": selected["cassette_count"],
        "蛋白数量": selected["protein_count"],
        "表达盒分组": [
            {
                "表达盒编号": cassette["cassette_index"],
                "蛋白": list(cassette["protein_accessions"]),
            }
            for cassette in selected["cassettes"]
        ],
        "清单是否更新": not unchanged,
        "清单文件": str(manifest_path),
        "清单版本": updated_manifest["revision"],
    }


def run_write_expression_box_selection(config: Any) -> dict[str, Any]:
    """CLI entry point for ``write --expression-box N``."""

    result = write_expression_box_selection(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "EXPRESSION_BOX_SELECTION_DOWNSTREAM_SECTIONS",
    "EXPRESSION_BOX_SELECTION_SCHEMA_VERSION",
    "run_write_expression_box_selection",
    "write_expression_box_selection",
]
