"""Write the complete auxiliary-protein research result to the manifest."""

from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from src.pathway_analyze.target_id import validate_target_compound_id
from src.protein_selection import (
    AuxiliaryProteinCombinationResult,
    auxiliary_protein_input_fingerprint,
    auxiliary_protein_result_path,
)
from src.protein_selection.manifest_adapter import (
    build_main_enzyme_research_units,
)
from src.protein_selection.reaction_scope import normalize_reaction_ids
from src.write_manifest.store import (
    read_design_manifest,
    update_design_manifest,
)


AUXILIARY_PROTEIN_RESEARCH_MANIFEST_SCHEMA_VERSION = (
    "auxiliary_protein_research_manifest.v1"
)
AUXILIARY_PROTEIN_RESEARCH_DOWNSTREAM_SECTIONS = (
    "auxiliary_protein_selection",
    "enzyme_system_selection",
    "cds_selection",
    "expression_box_selection",
    "expression_cassette_assembly",
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
)


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少 {field_name}")
    return value


def _integer(value: Any, field_name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ValueError(f"manifest {field_name} 必须是整数")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"manifest {field_name} 必须是整数") from exc
    if parsed < minimum:
        raise ValueError(f"manifest {field_name} 不能小于 {minimum}")
    return parsed


def _research_result_path(config: Any) -> Path:
    output_dir = Path(config.project_output_path).expanduser() / "protein_selection"
    return auxiliary_protein_result_path(
        config.manifest_output_path,
        output_dir,
    )


def _read_research_result(path: Path) -> AuxiliaryProteinCombinationResult:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到辅助蛋白研究结果，请先运行："
            "auxiliary-protein -i <输入文件>"
        )
    try:
        return AuxiliaryProteinCombinationResult.model_validate_json(
            path.read_text(encoding="utf-8")
        )
    except (OSError, ValidationError, ValueError) as exc:
        raise ValueError(
            "辅助蛋白研究结果无效或不是当前v3格式，请重新运行："
            "auxiliary-protein -i <输入文件>"
        ) from exc


def _unit_reaction_ids(unit: Mapping[str, Any]) -> list[str]:
    return normalize_reaction_ids(
        [
            str(step.get("reaction_id") or "")
            for step in unit.get("reaction_steps", [])
            if isinstance(step, Mapping)
        ]
    )


def _validate_manifest_binding(
    result: AuxiliaryProteinCombinationResult,
    manifest: Mapping[str, Any],
) -> None:
    solution = _mapping(manifest.get("solution"), "solution")
    selection = _mapping(
        manifest.get("main_enzyme_selection"),
        "main_enzyme_selection",
    )
    source = _mapping(
        selection.get("source"),
        "main_enzyme_selection.source",
    )
    expected = {
        "target_compound_id": str(
            manifest.get("target_compound_id") or ""
        ).upper(),
        "selected_solution_id": _integer(
            solution.get("solution_id"),
            "solution.solution_id",
            minimum=1,
        ),
        "expansion_depth": _integer(
            solution.get("expansion_depth", 0),
            "solution.expansion_depth",
        ),
        "selected_set_id": _integer(
            selection.get("selected_set_id"),
            "main_enzyme_selection.selected_set_id",
            minimum=1,
        ),
        "selected_set_fingerprint": str(
            selection.get("selected_set_fingerprint") or ""
        ).lower(),
        "solution_fingerprint": str(
            source.get("solution_fingerprint") or ""
        ).lower(),
        "chassis_key": str(selection.get("chassis_key") or ""),
    }
    for field_name, expected_value in expected.items():
        actual_value = getattr(result, field_name)
        if actual_value != expected_value:
            raise ValueError(
                f"辅助蛋白研究结果的 {field_name} 已过期；"
                "请重新运行 auxiliary-protein"
            )

    current_revision = _integer(
        manifest.get("revision", 0),
        "revision",
    )
    if result.source_manifest_revision > current_revision:
        raise ValueError("辅助蛋白研究结果来自未来的 manifest revision")

    units = build_main_enzyme_research_units(manifest)
    expected_input_fingerprint = auxiliary_protein_input_fingerprint(
        manifest,
        units,
        result.research_mode,
    )
    if result.input_fingerprint != expected_input_fingerprint:
        raise ValueError(
            "当前路线或主酶组合已经变化，请重新运行 auxiliary-protein"
        )

    if len(result.main_enzyme_results) != len(units):
        raise ValueError("辅助蛋白研究结果的主酶数量已过期")
    for item, unit in zip(result.main_enzyme_results, units, strict=True):
        if item.accession != str(unit["accession"]).upper():
            raise ValueError("辅助蛋白研究结果的主酶 accession 已过期")
        if item.sequence_sha256 != str(unit["sequence_sha256"]).lower():
            raise ValueError("辅助蛋白研究结果的主酶序列已经变化")
        if item.assigned_step_indexes != list(unit["assigned_step_indexes"]):
            raise ValueError("辅助蛋白研究结果的主酶步骤分配已经变化")
        if item.reaction_ids != _unit_reaction_ids(unit):
            raise ValueError("辅助蛋白研究结果的反应范围已经变化")


def _confirmed_auxiliary_payload(protein: Any) -> dict[str, Any]:
    return {
        "accession": protein.uniprot_id,
        "protein_name": protein.protein_name,
        "role": protein.role,
        "necessity": protein.necessity,
        "organism_name": protein.organism_name,
        "taxon_id": protein.taxon_id,
        "availability": protein.availability,
        "confidence": protein.confidence,
        "reason": protein.reason,
        "evidence_citation_ids": list(protein.evidence_citation_ids),
        "dependency_synthesis_ids": list(protein.dependency_synthesis_ids),
        "curated_assertion_ids": list(protein.curated_assertion_ids),
    }


def _candidate_auxiliary_payload(protein: Any) -> dict[str, Any]:
    return {
        "accession": protein.uniprot_id,
        "protein_name": protein.protein_name,
        "role": protein.role,
        "proposed_necessity": protein.proposed_necessity,
        "reason": protein.reason,
        "evidence_citation_ids": list(protein.evidence_citation_ids),
        "unresolved_reasons": list(protein.unresolved_reasons),
    }


def _manifest_payload(
    result: AuxiliaryProteinCombinationResult,
) -> dict[str, Any]:
    main_enzymes: list[dict[str, Any]] = []
    for item in result.main_enzyme_results:
        research = item.research_result
        main_enzymes.append({
            "accession": item.accession,
            "sequence_sha256": item.sequence_sha256,
            "assigned_step_indexes": item.assigned_step_indexes,
            "reaction_ids": item.reaction_ids,
            "research_status": item.status,
            "reaction_match": (
                research.reaction_match if research is not None else None
            ),
            "research_outcome": (
                research.outcome if research is not None else None
            ),
            "auxiliary_requirement_assessment": (
                item.auxiliary_requirement_assessment.model_dump(mode="json")
            ),
            "confirmed_auxiliary_proteins": (
                [
                    _confirmed_auxiliary_payload(protein)
                    for protein in research.auxiliary_proteins
                ]
                if research is not None
                else []
            ),
            "candidate_auxiliary_proteins": (
                [
                    _candidate_auxiliary_payload(protein)
                    for protein in research.candidate_auxiliary_proteins
                ]
                if research is not None
                else []
            ),
            "error": item.error,
        })

    return {
        "schema_version": AUXILIARY_PROTEIN_RESEARCH_MANIFEST_SCHEMA_VERSION,
        "research_status": result.status,
        "can_advance": result.can_advance,
        "selected_solution_id": result.selected_solution_id,
        "expansion_depth": result.expansion_depth,
        "selected_set_id": result.selected_set_id,
        "selected_set_fingerprint": result.selected_set_fingerprint,
        "chassis_key": result.chassis_key,
        "source": {
            "artifact": "protein_selection/auxiliary_protein_research.json",
            "artifact_schema_version": result.schema_version,
            "pipeline_version": result.pipeline_version,
            "source_manifest_revision": result.source_manifest_revision,
            "solution_fingerprint": result.solution_fingerprint,
            "input_fingerprint": result.input_fingerprint,
            "result_fingerprint": result.result_fingerprint,
            "research_mode": result.research_mode,
        },
        "main_enzymes": main_enzymes,
        "required_auxiliary_protein_accessions": (
            result.required_auxiliary_protein_accessions
        ),
        "recommended_auxiliary_protein_accessions": (
            result.recommended_auxiliary_protein_accessions
        ),
        "host_available_auxiliary_protein_accessions": (
            result.host_available_auxiliary_protein_accessions
        ),
        "auxiliary_proteins_to_introduce": (
            result.auxiliary_proteins_to_introduce
        ),
        "candidate_auxiliary_protein_accessions": (
            result.candidate_auxiliary_protein_accessions
        ),
        "complete_protein_list": result.complete_protein_list,
        "blocking_reasons": result.blocking_reasons,
        "warnings": result.warnings,
    }


def _scientific_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    normalized = dict(payload)
    source = normalized.get("source")
    if isinstance(source, Mapping):
        normalized["source"] = {
            key: value
            for key, value in source.items()
            if key != "source_manifest_revision"
        }
    return normalized


def _classification_lists(
    result: AuxiliaryProteinCombinationResult,
) -> dict[str, list[str]]:
    return {
        classification: [
            item.accession
            for item in result.main_enzyme_results
            if item.auxiliary_requirement_assessment.classification
            == classification
        ]
        for classification in (
            "confirmed_required",
            "possibly_required",
            "possibly_not_required",
            "confirmed_not_required",
            "unknown",
        )
    }


def write_auxiliary_protein_research(config: Any) -> dict[str, Any]:
    """Validate and commit all auxiliary-protein research conclusions."""

    target_compound = validate_target_compound_id(config.target_name)
    manifest_path = Path(config.manifest_output_path).expanduser()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if not recorded_target:
        raise ValueError(
            "manifest 尚未写入主酶组合，请先运行 write --main-enzyme-set N"
        )
    if recorded_target != target_compound:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target}，"
            f"与当前输入目标 {target_compound} 不一致"
        )
    if not isinstance(manifest.get("main_enzyme_selection"), Mapping):
        raise ValueError(
            "manifest 尚未写入主酶组合，请先运行 write --main-enzyme-set N"
        )

    result_path = _research_result_path(config)
    result = _read_research_result(result_path)
    _validate_manifest_binding(result, manifest)
    payload = _manifest_payload(result)

    current_payload = manifest.get("protein_selection")
    unchanged = (
        isinstance(current_payload, Mapping)
        and current_payload.get("schema_version")
        == AUXILIARY_PROTEIN_RESEARCH_MANIFEST_SCHEMA_VERSION
        and _scientific_payload(current_payload) == _scientific_payload(payload)
    )
    if unchanged:
        updated_manifest = manifest
    else:
        updated_manifest = update_design_manifest(
            manifest_path,
            target_compound_id=target_compound,
            sections={"protein_selection": payload},
            discard_sections=AUXILIARY_PROTEIN_RESEARCH_DOWNSTREAM_SECTIONS,
        )

    classifications = _classification_lists(result)
    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "主酶组合编号": result.selected_set_id,
        "研究状态": result.status,
        "可以进入后续流程": result.can_advance,
        "确认需要辅助蛋白的主酶": classifications["confirmed_required"],
        "可能需要辅助蛋白的主酶": classifications["possibly_required"],
        "可能不需要辅助蛋白的主酶": classifications[
            "possibly_not_required"
        ],
        "确认不需要辅助蛋白的主酶": classifications[
            "confirmed_not_required"
        ],
        "辅助蛋白需求未知的主酶": classifications["unknown"],
        "必需辅助蛋白": result.required_auxiliary_protein_accessions,
        "推荐辅助蛋白": result.recommended_auxiliary_protein_accessions,
        "候选辅助蛋白": result.candidate_auxiliary_protein_accessions,
        "需要导入的辅助蛋白": result.auxiliary_proteins_to_introduce,
        "阻断原因": result.blocking_reasons,
        "风险警告": result.warnings,
        "清单文件": str(manifest_path.resolve(strict=False)),
        "清单版本": updated_manifest["revision"],
    }


def run_write_auxiliary_protein_research(config: Any) -> dict[str, Any]:
    """CLI entry point for ``write --auxiliary-protein``."""

    result = write_auxiliary_protein_research(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "AUXILIARY_PROTEIN_RESEARCH_MANIFEST_SCHEMA_VERSION",
    "run_write_auxiliary_protein_research",
    "write_auxiliary_protein_research",
]
