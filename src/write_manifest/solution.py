from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from src.pathway_analyze.gem_validation import validation_depth_output_dir
from src.pathway_analyze.kegg_gap_analyze import (
    ELECTRON_INFERENCE_VERSION,
    gap_depth_output_dir,
)
from src.pathway_analyze.target_id import validate_target_compound_id
from src.write_manifest.store import update_design_manifest


SOLUTION_SUMMARY_FIELDS = (
    "target_compound_id",
    "target_compound_name",
    "total_steps",
    "heterologous_steps",
    "heterologous_reaction_ids",
    "heterologous_ko_ids",
    "heterologous_enzyme_ecs",
    "reaction_resolution_status",
    "normalization_event_count",
    "normalization_events",
    "blocking_reaction_count",
    "blocking_reaction_ids",
    "eligible_for_recommendation",
    "reachable_anchor_compounds",
    "reachable_anchor_labels",
)

RETROPATH_SOLUTION_SUMMARY_FIELDS = (
    "solution_source",
    "retropath_candidate_rank",
    "retropath_candidate_id",
    "retropath_combination_id",
    "prediction_review_required",
    "materialization_id",
    "validation_status",
    "stoichiometry_status",
    "gem_status",
    "validation_issue",
    "promotion_id",
    "combination_truncated",
    "upstream_enumeration_truncated",
    "candidate_top_k_truncated",
)

SOLUTION_STEP_FIELDS = (
    "step_index",
    "gap_step_index",
    "status",
    "reaction_id",
    "reaction_name",
    "reaction_comment",
    "equation",
    "direction",
    "produced_compound_id",
    "produced_compound_name",
    "precursor_compound_ids",
    "precursor_compound_labels",
    "ko_ids",
    "module_ids",
    "enzyme_ecs",
    "rhea_ids",
    "oxygen_required",
    "thermo_direction",
    "screening_rule_hits",
    "source_reaction_ids",
    "resolution_action",
    "resolution_evidence",
    "auxiliary_requirements",
)

RETROPATH_SOLUTION_STEP_FIELDS = (
    "step_source",
    "locked_enzyme_ecs",
    "ec_status",
    "enzyme_search_eligible",
    "retropath_step_id",
    "retropath_hypothesis_id",
    "retropath_rule_id",
    "source_mnxr_id",
    "source_ec_numbers",
    "source_uniprot_ids",
    "source_rhea_ids",
    "exact_kegg_reaction_ids",
    "exact_rhea_ids",
    "formal_mapping_exact",
    "reaction_signature_sha256",
    "full_reaction_smiles",
    "core_reaction_smiles",
    "stoichiometry_terms",
    "prediction_provenance",
    "prediction_review_required",
    "depends_on_step_ids",
    "validation_status",
    "stoichiometry_status",
    "gem_status",
)

VALIDATION_FIELDS = (
    "validation_status",
    "fba_status",
    "fba_product_flux",
    "fba_growth_flux",
    "pfba_status",
    "pfba_product_flux",
    "pfba_growth_flux",
    "fva_status",
    "target_fva_minimum",
    "target_fva_maximum",
    "route_reaction_count",
    "active_route_reaction_count",
    "required_route_reaction_count",
    "blocked_route_reaction_count",
    "active_route_reaction_ids",
    "required_route_reaction_ids",
    "blocked_route_reaction_ids",
    "cofactor_mode",
    "cofactor_relaxed",
    "opened_generic_compound_ids",
    "issues",
)

RETROPATH_VALIDATION_FIELDS = (
    "retropath_candidate_rank",
    "retropath_candidate_id",
    "retropath_combination_id",
    "stoichiometry_hypothesis_ids",
    "combination_truncated",
)

ELECTRON_SUMMARY_FIELDS = (
    "max_electron_risk_level",
    "max_electron_risk_score",
    "electron_system_status",
    "electron_balance_status",
    "requires_external_electron_regeneration",
    "requires_carrier_compatibility_check",
    "electron_carrier_ids",
    "electron_requirement_classes",
    "electron_carrier_net_changes",
    "balanced_electron_carrier_pairs",
    "unbalanced_electron_carrier_pairs",
    "unresolved_electron_carrier_ids",
    "annotation_only_electron_requirements",
    "required_auxiliary_roles",
    "auxiliary_requirement_status",
)

ELECTRON_RISK_STEP_FIELDS = (
    "step_index",
    "gap_step_index",
    "reaction_id",
    "electron_carrier_ids",
    "electron_requirement_classes",
    "electron_risk_level",
    "electron_risk_score",
    "electron_risk_evidence",
    "electron_carrier_net_changes",
    "required_auxiliary_roles",
    "auxiliary_requirement_status",
    "auxiliary_requirements",
)

SOLUTION_DOWNSTREAM_SECTIONS = (
    "main_enzyme_selection",
    "protein_selection",
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


def _read_csv_rows(path: Path, *, required: bool = True) -> list[dict[str, str]]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"缺少文件：{path}")
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _convert_value(value: Any) -> Any:
    text = str(value or "").strip()
    if not text:
        return None
    if text.lower() == "true":
        return True
    if text.lower() == "false":
        return False
    try:
        number = float(text)
    except ValueError:
        return text
    return int(number) if number.is_integer() else number


def _normalize_row(row: dict[str, str]) -> dict[str, Any]:
    return {key: _convert_value(value) for key, value in row.items()}


def _keep_fields(row: dict[str, Any], fields: tuple[str, ...]) -> dict[str, Any]:
    return {key: row.get(key) for key in fields}


def _row_solution_id(row: dict[str, str]) -> int:
    try:
        return int(str(row.get("solution_id", "")).strip())
    except ValueError as exc:
        raise ValueError(
            f"CSV 中存在无效的 solution_id：{row.get('solution_id')!r}"
        ) from exc


def _select_solution_summary(path: Path, solution_id: int) -> dict[str, Any]:
    rows = _read_csv_rows(path)
    for row in rows:
        if _row_solution_id(row) == solution_id:
            return _normalize_row(row)
    available = sorted({_row_solution_id(row) for row in rows})
    raise ValueError(f"未找到 solution {solution_id}；可用编号：{available}")


def _select_solution_steps(path: Path, solution_id: int) -> list[dict[str, Any]]:
    rows = [row for row in _read_csv_rows(path) if _row_solution_id(row) == solution_id]
    if not rows:
        raise ValueError(f"solution {solution_id} 没有反应步骤：{path}")
    rows.sort(key=lambda row: int(str(row.get("step_index", "0")).strip()))
    return [_normalize_row(row) for row in rows]


def _split_values(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def _parse_auxiliary_requirements(value: Any) -> list[dict[str, Any]]:
    text = str(value or "").strip()
    if not text:
        return []
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError("路线辅助蛋白需求字段不是有效 JSON") from exc
    if not isinstance(payload, list) or any(
        not isinstance(item, dict) for item in payload
    ):
        raise ValueError("路线辅助蛋白需求字段必须是对象列表")
    required_fields = {
        "role",
        "necessity",
        "confidence",
        "selection_status",
        "carrier_ids",
        "evidence",
    }
    for item in payload:
        if required_fields.difference(item):
            raise ValueError("路线辅助蛋白需求字段缺少必要信息")
    return payload


def _parse_json_value(
    value: Any,
    *,
    field_name: str,
    expected_type: type,
) -> Any:
    text = str(value or "").strip()
    if not text:
        return expected_type()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"路线 {field_name} 字段不是有效 JSON") from exc
    if not isinstance(payload, expected_type):
        raise ValueError(f"路线 {field_name} 字段类型无效")
    return payload


def _integer_step_index(value: Any) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"无效的 step_index：{value!r}") from exc


def _order_steps_forward(
    rows: list[dict[str, Any]],
    target_compound: str,
) -> list[dict[str, Any]]:
    """按底盘前体到目标产物的依赖顺序排列路线步骤。"""

    row_by_product: dict[str, dict[str, Any]] = {}
    for row in rows:
        product_id = str(row.get("produced_compound_id") or "").strip()
        if not product_id:
            raise ValueError(
                f"路线步骤缺少 produced_compound_id：step {row.get('step_index')}"
            )
        if product_id in row_by_product:
            raise ValueError(f"路线中存在重复产物，无法稳定编号：{product_id}")
        row_by_product[product_id] = row

    ordered_rows: list[dict[str, Any]] = []
    visited_products: set[str] = set()
    active_products: set[str] = set()

    def visit(product_id: str) -> None:
        if product_id in visited_products:
            return
        if product_id in active_products:
            raise ValueError(f"路线步骤存在循环依赖：{product_id}")
        row = row_by_product.get(product_id)
        if row is None:
            return
        active_products.add(product_id)
        for precursor_id in _split_values(row.get("precursor_compound_ids")):
            visit(precursor_id)
        active_products.remove(product_id)
        visited_products.add(product_id)
        ordered_rows.append(row)

    visit(target_compound)
    for row in reversed(rows):
        visit(str(row.get("produced_compound_id") or "").strip())
    if len(ordered_rows) != len(rows):
        raise ValueError("路线步骤无法完整转换为正向顺序")
    return ordered_rows


def _renumber_steps_forward(
    rows: list[dict[str, Any]],
    target_compound: str,
) -> tuple[list[dict[str, Any]], dict[int, int]]:
    """生成正向步骤，并返回 gap 原始编号到正向编号的映射。"""

    forward_rows: list[dict[str, Any]] = []
    gap_to_forward: dict[int, int] = {}
    for forward_index, row in enumerate(
        _order_steps_forward(rows, target_compound),
        start=1,
    ):
        gap_step_index = _integer_step_index(row.get("step_index"))
        if gap_step_index in gap_to_forward:
            raise ValueError(f"路线中存在重复 step_index：{gap_step_index}")
        gap_to_forward[gap_step_index] = forward_index
        forward_row = dict(row)
        forward_row["gap_step_index"] = gap_step_index
        forward_row["step_index"] = forward_index
        forward_row["auxiliary_requirements"] = (
            _parse_auxiliary_requirements(
                forward_row.get("auxiliary_requirements_json")
            )
        )
        if str(forward_row.get("step_source") or "").strip() == "retropath":
            forward_row["stoichiometry_terms"] = _parse_json_value(
                forward_row.get("stoichiometry_terms_json"),
                field_name="stoichiometry_terms",
                expected_type=list,
            )
            forward_row["prediction_provenance"] = _parse_json_value(
                forward_row.get("prediction_provenance_json"),
                field_name="prediction_provenance",
                expected_type=dict,
            )
            hypothesis_id = str(
                forward_row.get("retropath_hypothesis_id") or ""
            ).strip()
            step_id = str(
                forward_row.get("retropath_step_id")
                or forward_row["prediction_provenance"].get("step_id")
                or ""
            ).strip()
            reaction_id = str(forward_row.get("reaction_id") or "").strip()
            provenance_hypothesis = str(
                forward_row["prediction_provenance"].get("hypothesis_id") or ""
            ).strip()
            identity_ok = bool(step_id) and reaction_id in {
                step_id,
                hypothesis_id,
            }
            if hypothesis_id:
                identity_ok = identity_ok and provenance_hypothesis == hypothesis_id
            else:
                identity_ok = identity_ok and reaction_id == step_id
            if not identity_ok:
                raise ValueError("RetroPath 路线步骤的计量假设身份不一致")
        forward_rows.append(forward_row)
    return forward_rows, gap_to_forward


def _remap_electron_steps(
    rows: list[dict[str, Any]],
    gap_to_forward: dict[int, int],
) -> list[dict[str, Any]]:
    remapped_rows: list[dict[str, Any]] = []
    for row in rows:
        gap_step_index = _integer_step_index(row.get("step_index"))
        if gap_step_index not in gap_to_forward:
            raise ValueError(
                f"电子风险步骤引用了不存在的 gap step：{gap_step_index}"
            )
        remapped_row = dict(row)
        remapped_row["gap_step_index"] = gap_step_index
        remapped_row["step_index"] = gap_to_forward[gap_step_index]
        remapped_row["auxiliary_requirements"] = (
            _parse_auxiliary_requirements(
                remapped_row.get("auxiliary_requirements_json")
            )
        )
        remapped_rows.append(remapped_row)
    remapped_rows.sort(key=lambda row: int(row["step_index"]))
    return remapped_rows


def _select_optional_solution_row(path: Path, solution_id: int) -> dict[str, Any]:
    for row in _read_csv_rows(path, required=False):
        if _row_solution_id(row) == solution_id:
            return _normalize_row(row)
    return {}


def _select_optional_solution_rows(
    path: Path,
    solution_id: int,
) -> list[dict[str, Any]]:
    rows = [
        row
        for row in _read_csv_rows(path, required=False)
        if _row_solution_id(row) == solution_id
    ]
    rows.sort(key=lambda row: int(str(row.get("step_index", "0")).strip()))
    return [_normalize_row(row) for row in rows]


def _parse_solution_ids(value: Any) -> tuple[int, ...]:
    parts = str(value or "").replace(",", ";").split(";")
    try:
        return tuple(int(part.strip()) for part in parts if part.strip())
    except ValueError as exc:
        raise ValueError(f"通量验证文件中的 solution_ids 无效：{value!r}") from exc


def _select_passed_validation(
    path: Path,
    solution_id: int,
    expansion_depth: int,
) -> dict[str, Any]:
    validation_command = (
        f"validate -i <输入文件> -s {solution_id} -m per -d {expansion_depth}"
    )
    if not path.is_file():
        raise FileNotFoundError(
            f"solution {solution_id} 尚未进行独立 GEM 验证；请先执行："
            f"{validation_command}"
        )
    matches = [
        row
        for row in _read_csv_rows(path)
        if str(row.get("validation_mode", "")).strip() == "per-solution"
        and _parse_solution_ids(row.get("solution_ids")) == (solution_id,)
    ]
    if not matches:
        raise ValueError(
            f"solution {solution_id} 没有独立通量验证结果；请先执行 "
            f"{validation_command}"
        )

    validation = _normalize_row(matches[-1])
    status = str(validation.get("validation_status") or "")
    if not status.startswith("PASS_"):
        raise ValueError(
            f"solution {solution_id} 未通过通量验证，当前状态：{status or '未知'}"
        )
    return validation


def _select_retropath_validation(
    promotion: dict[str, Any],
    solution_id: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    mappings = [
        item
        for item in promotion.get("solution_mappings", [])
        if isinstance(item, dict)
        and int(item.get("solution_id") or 0) == solution_id
    ]
    if len(mappings) != 1:
        raise ValueError(f"RetroPath solution {solution_id} 的晋升身份不唯一")
    mapping = mappings[0]
    p8_manifest_path = Path(
        str(
            promotion.get("inputs", {}).get("p8_validation_manifest")
            or ""
        )
    ).expanduser().resolve()
    p8_manifest = json.loads(p8_manifest_path.read_text(encoding="utf-8-sig"))
    summary_record = p8_manifest.get("artifacts", {}).get(
        "gem_validation_summary.csv"
    )
    if not isinstance(summary_record, dict):
        raise ValueError("P8 validation summary 记录缺失")
    summary_path = Path(str(summary_record.get("path") or "")).expanduser().resolve()
    matches = [
        _normalize_row(row)
        for row in _read_csv_rows(summary_path)
        if str(row.get("candidate_id") or "").strip()
        == str(mapping.get("candidate_id") or "").strip()
        and str(row.get("combination_id") or "").strip()
        == str(mapping.get("combination_id") or "").strip()
        and str(row.get("validation_status") or "").strip()
        == "PASS_STRICT_ROUTE_FLUX"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"RetroPath solution {solution_id} 没有唯一的 P8 严格通过记录"
        )
    row = matches[0]
    normalized = {
        "validation_status": row.get("validation_status"),
        "fba_status": row.get("fba_status"),
        "fba_product_flux": row.get("target_flux"),
        "fba_growth_flux": row.get("growth_flux"),
        "pfba_status": row.get("pfba_status"),
        "pfba_product_flux": row.get("pfba_target_flux"),
        "pfba_growth_flux": row.get("growth_flux"),
        "fva_status": row.get("fva_status"),
        "target_fva_minimum": None,
        "target_fva_maximum": None,
        "route_reaction_count": row.get("route_step_count"),
        "active_route_reaction_count": row.get("active_route_step_count"),
        "required_route_reaction_count": row.get("route_step_count"),
        "blocked_route_reaction_count": len(
            _split_values(row.get("blocked_route_step_ids"))
        ),
        "active_route_reaction_ids": None,
        "required_route_reaction_ids": None,
        "blocked_route_reaction_ids": row.get("blocked_route_step_ids"),
        "cofactor_mode": "strict",
        "cofactor_relaxed": False,
        "opened_generic_compound_ids": None,
        "issues": row.get("issues"),
        "retropath_candidate_rank": mapping.get("candidate_rank"),
        "retropath_candidate_id": mapping.get("candidate_id"),
        "retropath_combination_id": mapping.get("combination_id"),
        "stoichiometry_hypothesis_ids": mapping.get(
            "stoichiometry_hypothesis_ids", []
        ),
        "combination_truncated": row.get("combination_truncated"),
    }
    return normalized, mapping


def _select_retropath_materialization_validation(
    materialization: dict[str, Any],
    solution_id: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    mappings = [
        item
        for item in materialization.get("solution_mappings", [])
        if isinstance(item, dict)
        and int(item.get("solution_id") or 0) == solution_id
    ]
    if len(mappings) != 1:
        raise ValueError(f"RetroPath solution {solution_id} 的物化身份不唯一")
    mapping = mappings[0]
    state = mapping.get("validation")
    if not isinstance(state, dict):
        state = {
            "status": "not_run",
            "stoichiometry_status": "core_only",
            "gem_status": "not_run",
            "combination_id": "",
            "stoichiometry_hypothesis_ids": [],
            "issues": [],
        }
    normalized: dict[str, Any] = {
        "validation_status": state.get("status", "not_run"),
        "fba_status": state.get("gem_status", "not_run"),
        "cofactor_mode": (
            "strict" if state.get("stoichiometry_status") == "completed" else "not_run"
        ),
        "cofactor_relaxed": False,
        "issues": "; ".join(str(item) for item in state.get("issues", [])),
        "retropath_candidate_rank": mapping.get("candidate_rank"),
        "retropath_candidate_id": mapping.get("candidate_id"),
        "retropath_combination_id": state.get("combination_id", ""),
        "stoichiometry_hypothesis_ids": state.get(
            "stoichiometry_hypothesis_ids", []
        ),
        "combination_truncated": False,
    }
    if state.get("status") != "passed":
        return normalized, mapping

    overlay = materialization.get("validation_overlay")
    p8_manifest_value = state.get("manifest_path")
    if not p8_manifest_value and isinstance(overlay, dict):
        p8_manifest_value = overlay.get("path")
    if not p8_manifest_value:
        raise ValueError("RetroPath solution 声明验证通过但缺少 P8 overlay")
    p8_manifest_path = Path(str(p8_manifest_value)).expanduser().resolve()
    p8_manifest = json.loads(p8_manifest_path.read_text(encoding="utf-8-sig"))
    summary_record = p8_manifest.get("artifacts", {}).get(
        "gem_validation_summary.csv"
    )
    if not isinstance(summary_record, dict):
        raise ValueError("P8 validation summary 记录缺失")
    summary_path = Path(str(summary_record.get("path") or "")).expanduser().resolve()
    matches = [
        _normalize_row(row)
        for row in _read_csv_rows(summary_path)
        if str(row.get("candidate_id") or "").strip()
        == str(mapping.get("candidate_id") or "").strip()
        and str(row.get("combination_id") or "").strip()
        == str(state.get("combination_id") or "").strip()
        and str(row.get("validation_status") or "").strip()
        == "PASS_STRICT_ROUTE_FLUX"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"RetroPath solution {solution_id} 没有唯一的 P8 严格通过记录"
        )
    row = matches[0]
    normalized.update(
        {
            "validation_status": row.get("validation_status"),
            "fba_status": row.get("fba_status"),
            "fba_product_flux": row.get("target_flux"),
            "fba_growth_flux": row.get("growth_flux"),
            "pfba_status": row.get("pfba_status"),
            "pfba_product_flux": row.get("pfba_target_flux"),
            "pfba_growth_flux": row.get("growth_flux"),
            "fva_status": row.get("fva_status"),
            "route_reaction_count": row.get("route_step_count"),
            "active_route_reaction_count": row.get("active_route_step_count"),
            "required_route_reaction_count": row.get("route_step_count"),
            "blocked_route_reaction_count": len(
                _split_values(row.get("blocked_route_step_ids"))
            ),
            "blocked_route_reaction_ids": row.get("blocked_route_step_ids"),
            "issues": row.get("issues"),
            "combination_truncated": row.get("combination_truncated"),
        }
    )
    return normalized, mapping


def write_solution(config: Any) -> dict[str, Any]:
    """将一条 KEGG 或可选验证的 RetroPath 路线写入 design manifest。"""

    target_compound = validate_target_compound_id(config.target_name)
    expansion_depth = int(getattr(config, "depth", 0))
    if expansion_depth < 0:
        raise ValueError(
            "depth 必须大于等于 0 (must be greater than or equal to 0)"
        )
    configured_solution_id = getattr(config, "solution", None)
    if configured_solution_id is None:
        raise ValueError("未指定 solution，请使用 write --solution N")
    solution_id = int(configured_solution_id)
    if solution_id < 1:
        raise ValueError("solution 必须是正整数")

    gap_root_dir = Path(config.gap_output_path).expanduser().resolve()
    gap_dir = gap_depth_output_dir(gap_root_dir, expansion_depth)
    validation_dir = validation_depth_output_dir(gap_root_dir, expansion_depth)
    manifest_path = Path(config.manifest_output_path).expanduser()
    solutions_path = gap_dir / "solutions.csv"
    if not solutions_path.is_file() and (gap_dir / "run_config.json").is_file():
        try:
            early_run_config = json.loads(
                (gap_dir / "run_config.json").read_text(encoding="utf-8-sig")
            )
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"无法读取 gap 运行配置：{gap_dir / 'run_config.json'}"
            ) from exc
        if (
            not isinstance(early_run_config, dict)
            or early_run_config.get("electron_inference_version")
            != ELECTRON_INFERENCE_VERSION
        ):
            raise ValueError(
                "gap 结果使用旧版电子推断，请重新运行："
                f"gap -i <输入文件> -d {expansion_depth}"
            )
    summary = _select_solution_summary(solutions_path, solution_id)
    solution_source = str(summary.get("solution_source") or "kegg").strip().lower()
    promotion: dict[str, Any] | None = None
    materialization: dict[str, Any] | None = None
    promotion_mapping: dict[str, Any] | None = None
    if solution_source == "retropath":
        from src.pathway_analyze.retropath_materialization import (
            MATERIALIZATION_MANIFEST_FILE_NAME,
            verify_retropath_solution_materialization,
        )

        if (gap_dir / "retropath" / MATERIALIZATION_MANIFEST_FILE_NAME).is_file():
            materialization = verify_retropath_solution_materialization(
                gap_dir=gap_dir,
                target_compound=target_compound,
                expansion_depth=expansion_depth,
                solution_id=solution_id,
            )
        else:
            from src.pathway_analyze.retropath_promotion import (
                verify_retropath_solution_promotion,
            )

            promotion = verify_retropath_solution_promotion(
                gap_dir=gap_dir,
                target_compound=target_compound,
                expansion_depth=expansion_depth,
                solution_id=solution_id,
            )
    else:
        run_config_path = gap_dir / "run_config.json"
        try:
            run_config = json.loads(run_config_path.read_text(encoding="utf-8-sig"))
        except (FileNotFoundError, json.JSONDecodeError) as exc:
            raise ValueError(f"无法读取 gap 运行配置：{run_config_path}") from exc
        if (
            not isinstance(run_config, dict)
            or run_config.get("electron_inference_version")
            != ELECTRON_INFERENCE_VERSION
        ):
            raise ValueError(
                "gap 结果使用旧版电子推断，请重新运行："
                f"gap -i <输入文件> -d {expansion_depth}"
            )
    if summary.get("target_compound_id") != target_compound:
        raise ValueError(
            f"solution {solution_id} 的目标化合物是 {summary.get('target_compound_id')}，"
            f"与当前目标 {target_compound} 不一致"
        )
    if int(summary.get("blocking_reaction_count") or 0) > 0:
        raise ValueError(
            f"solution {solution_id} 含有阻断反应："
            f"{summary.get('blocking_reaction_ids') or '未知'}"
        )
    if summary.get("eligible_for_recommendation") is not True:
        raise ValueError(f"solution {solution_id} 不可作为推荐路径")

    gap_steps = _select_solution_steps(
        gap_dir / "all_solution_steps.csv",
        solution_id,
    )
    steps, gap_to_forward = _renumber_steps_forward(gap_steps, target_compound)
    if materialization is not None:
        validation, promotion_mapping = _select_retropath_materialization_validation(
            materialization,
            solution_id,
        )
    elif promotion is not None:
        validation, promotion_mapping = _select_retropath_validation(
            promotion,
            solution_id,
        )
    else:
        validation = _select_passed_validation(
            validation_dir / "gem_validation_summary.csv",
            solution_id,
            expansion_depth,
        )
    electron_summary = _select_optional_solution_row(
        gap_dir / "solution_electron_summary.csv",
        solution_id,
    )
    required_electron_fields = {
        "electron_system_status",
        "electron_balance_status",
        "requires_external_electron_regeneration",
        "requires_carrier_compatibility_check",
    }
    if (
        not electron_summary
        or required_electron_fields.difference(electron_summary)
    ):
        raise ValueError(
            "gap 结果使用旧版电子推断，请重新运行："
            f"gap -i <输入文件> -d {expansion_depth}"
        )
    electron_steps = _remap_electron_steps(
        _select_optional_solution_rows(
            gap_dir / "route_electron_requirements.csv",
            solution_id,
        ),
        gap_to_forward,
    )

    summary_payload = _keep_fields(summary, SOLUTION_SUMMARY_FIELDS)
    validation_payload = _keep_fields(validation, VALIDATION_FIELDS)
    step_payloads = [
        _keep_fields(step, SOLUTION_STEP_FIELDS) for step in steps
    ]
    if materialization is not None or promotion is not None:
        summary_payload.update(
            _keep_fields(summary, RETROPATH_SOLUTION_SUMMARY_FIELDS)
        )
        validation_payload.update(
            _keep_fields(validation, RETROPATH_VALIDATION_FIELDS)
        )
        step_payloads = [
            {
                **_keep_fields(step, SOLUTION_STEP_FIELDS),
                **_keep_fields(step, RETROPATH_SOLUTION_STEP_FIELDS),
            }
            for step in steps
        ]
    solution_payload = {
        "gap_dir": str(gap_dir.resolve()),
        "expansion_depth": expansion_depth,
        "solution_id": solution_id,
        "summary": summary_payload,
        "steps": step_payloads,
        "validation": validation_payload,
    }
    if (materialization is not None or promotion is not None) and promotion_mapping is not None:
        commit = materialization if materialization is not None else promotion
        assert commit is not None
        state = promotion_mapping.get("validation", {})
        solution_payload["prediction"] = {
            "source": "retropath",
            "review_status": "pending",
            "review_required": True,
            "materialization_id": commit.get(
                "materialization_id", commit.get("promotion_id")
            ),
            "materialization_schema_version": commit.get("schema_version"),
            # Compatibility aliases retained for existing main-enzyme/write
            # consumers.  New code should prefer the materialization names.
            "promotion_id": commit.get(
                "materialization_id", commit.get("promotion_id")
            ),
            "promotion_schema_version": commit.get("schema_version"),
            "candidate_rank": promotion_mapping.get("candidate_rank"),
            "candidate_id": promotion_mapping.get("candidate_id"),
            "combination_id": state.get(
                "combination_id", promotion_mapping.get("combination_id")
            ),
            "stoichiometry_hypothesis_ids": state.get(
                "stoichiometry_hypothesis_ids",
                promotion_mapping.get("stoichiometry_hypothesis_ids", []),
            ),
            "validation_status": state.get(
                "status", validation.get("validation_status", "not_run")
            ),
            "stoichiometry_status": state.get("stoichiometry_status", "completed"),
            "gem_status": state.get("gem_status", "passed"),
            "validation_issues": state.get("issues", []),
            "materialization_manifest": str(
                (
                    gap_dir
                    / "retropath"
                    / (
                        "solution_materialization.json"
                        if materialization is not None
                        else "formal_solution_promotion.json"
                    )
                ).resolve()
            ),
            "promotion_manifest": str(
                (
                    gap_dir
                    / "retropath"
                    / (
                        "solution_materialization.json"
                        if materialization is not None
                        else "formal_solution_promotion.json"
                    )
                ).resolve()
            ),
            "p8_validation_sha256": (
                state.get("manifest_sha256")
                or (commit.get("validation_overlay") or {}).get("sha256")
                if materialization is not None
                else commit.get("inputs", {}).get("p8_validation_sha256")
            ),
            "rr02_sha256": (
                commit.get("rr02_sha256")
                or commit.get("inputs", {}).get("rr02_sha256")
                or (commit.get("inputs", {}).get("rr02") or {}).get("sha256")
            ),
        }
    electron_payload = {
        "summary": (
            _keep_fields(electron_summary, ELECTRON_SUMMARY_FIELDS)
            if electron_summary
            else {}
        ),
        "risk_steps": [
            _keep_fields(step, ELECTRON_RISK_STEP_FIELDS) for step in electron_steps
        ],
    }

    manifest = update_design_manifest(
        manifest_path,
        target_compound_id=target_compound,
        sections={
            "solution": solution_payload,
            "electron_inference": electron_payload,
        },
        discard_sections=SOLUTION_DOWNSTREAM_SECTIONS,
    )

    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "路径编号": solution_id,
        "扩展深度": expansion_depth,
        "步骤数量": len(steps),
        "通量验证状态": validation.get("validation_status"),
        "辅因子模式": validation.get("cofactor_mode"),
        "电子系统状态": electron_summary.get("electron_system_status"),
        "电子载体平衡状态": electron_summary.get("electron_balance_status"),
        "是否需要额外电子再生系统": electron_summary.get(
            "requires_external_electron_regeneration"
        ),
        "是否需要确认真实载体兼容性": electron_summary.get(
            "requires_carrier_compatibility_check"
        ),
        "清单文件": str(manifest_path.resolve()),
        "清单版本": manifest["revision"],
    }


def run_write_solution(config: Any) -> dict[str, Any]:
    """``write --solution N`` 的命令行入口。"""

    result = write_solution(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["run_write_solution", "write_solution"]
