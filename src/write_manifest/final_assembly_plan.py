"""Validate and commit the complete final-assembly plan bundle."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.final_assemble_plan.common import stable_json_hash
from src.final_assemble_plan.config import (
    ASSEMBLY_PLAN_ALGORITHM_VERSION,
    ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME,
    ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION,
    FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION,
    SUPPORTED_METHODS,
)
from src.final_assemble_plan.get_final_assembly_context import (
    load_final_assembly_context,
)
from src.write_manifest.store import read_design_manifest, update_design_manifest


FINAL_ASSEMBLY_PLAN_DOWNSTREAM_SECTIONS = (
    "final_assembly",
    "final_design_report",
)


def _artifact_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "final_assemble_plan"
        / ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME
    )


def _read_artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到最终组装计划，请先运行：assembly --plan -i <输入文件>"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"最终组装计划不是有效 JSON：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("最终组装计划根节点必须是 JSON 对象")
    return payload


def _validate_method_payload(plan: Mapping[str, Any]) -> None:
    method = str(plan.get("assembly_method") or "")
    if method not in SUPPORTED_METHODS:
        raise ValueError(f"组装计划包含不支持的方法：{method}")
    target = plan.get("target")
    linearization = plan.get("backbone_linearization")
    if not isinstance(target, Mapping) or not isinstance(linearization, Mapping):
        raise ValueError("组装计划缺少 target 或 backbone_linearization")
    backbone = plan.get("backbone")
    try:
        backbone_length = int(
            backbone.get("length_bp") if isinstance(backbone, Mapping) else 0
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("组装计划包含无效 backbone length_bp") from exc
    if backbone_length < 1:
        raise ValueError("组装计划包含无效 backbone length_bp")
    if target.get("mode") == "insert_after":
        try:
            insert_after = int(target.get("insert_after_bp"))
        except (TypeError, ValueError) as exc:
            raise ValueError("insert_after 计划包含无效坐标") from exc
        if not 1 <= insert_after <= backbone_length:
            raise ValueError("insert_after 计划包含无效坐标")
    elif target.get("mode") == "replace":
        try:
            start = int(target.get("replace_start_bp"))
            end = int(target.get("replace_end_bp"))
        except (TypeError, ValueError) as exc:
            raise ValueError("replace 计划包含无效坐标") from exc
        if start < 1 or end < start or end > backbone_length:
            raise ValueError("replace 计划包含无效坐标")
    else:
        raise ValueError("组装计划 target.mode 无效")
    enzymes = linearization.get("restriction_enzymes")
    if not isinstance(enzymes, list):
        raise ValueError("组装计划缺少 restriction_enzymes 列表")
    mode = linearization.get("mode")
    if mode == "pcr":
        if enzymes or linearization.get("enzyme_summary") != "none (PCR linearization)":
            raise ValueError("PCR 线性化必须明确记录无限制酶")
    elif mode == "restriction":
        if not enzymes:
            raise ValueError("限制酶线性化必须记录酶")
        for enzyme in enzymes:
            if not isinstance(enzyme, Mapping):
                raise ValueError("限制酶线性化包含不完整的酶记录")
            if not enzyme.get("name") or not enzyme.get("recognition_site"):
                raise ValueError("限制酶线性化包含不完整的酶记录")
            if enzyme.get("site_start_bp") is None or enzyme.get("cut_after_bp") is None:
                raise ValueError("限制酶线性化包含不完整的酶记录")
    else:
        raise ValueError("backbone_linearization.mode 无效")
    if method == "restriction":
        restriction = plan.get("restriction")
        if not isinstance(restriction, Mapping) or not all(
            restriction.get(field)
            for field in ("left_enzyme", "right_enzyme", "left_site", "right_site")
        ):
            raise ValueError("restriction 计划缺少左右酶信息")
    else:
        gibson = plan.get("gibson")
        if not isinstance(gibson, Mapping) or not all(
            gibson.get(field)
            for field in ("homology_arm_length", "left_homology", "right_homology")
        ):
            raise ValueError("Gibson 计划缺少同源臂")


def _validate_artifact(
    artifact: Mapping[str, Any],
    *,
    context: Any,
) -> list[dict[str, Any]]:
    if (
        artifact.get("schema_version")
        != ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION
        or artifact.get("algorithm_version") != ASSEMBLY_PLAN_ALGORITHM_VERSION
    ):
        raise ValueError("最终组装计划版本已过期，请重新运行 assembly --plan")
    if artifact.get("status") != "complete":
        raise ValueError("最终组装计划不完整，不能写入 manifest")
    if artifact.get("target_compound_id") != context.target_compound_id:
        raise ValueError("最终组装计划与当前目标化合物不一致")
    source = artifact.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("最终组装计划缺少 source")
    expected_source = {
        "context_input_fingerprint": context.input_fingerprint,
        "parts_selection_fingerprint": context.parts_selection_fingerprint,
        "assembled_constructs_fingerprint": context.assembled_constructs_fingerprint,
        "plasmid_selection_fingerprint": context.plasmid_selection_fingerprint,
    }
    for key, value in expected_source.items():
        if source.get(key) != value:
            raise ValueError(
                "最终组装计划与当前表达构建或质粒选择不一致，请重新运行 assembly --plan"
            )
    raw_plans = artifact.get("plans")
    if not isinstance(raw_plans, list) or not raw_plans:
        raise ValueError("最终组装计划 plans 为空")
    plans = [dict(item) for item in raw_plans if isinstance(item, Mapping)]
    if len(plans) != len(raw_plans):
        raise ValueError("最终组装计划包含无效条目")
    expected_ids = {item.design_id for item in context.constructs}
    plan_ids = {int(item.get("parts_design_id") or 0) for item in plans}
    if (
        len(plans) != len(context.constructs)
        or plan_ids != expected_ids
        or int(artifact.get("design_count") or 0) != len(context.constructs)
    ):
        raise ValueError("最终组装计划没有覆盖全部已选表达构建")
    construct_by_id = {item.design_id: item for item in context.constructs}
    for plan in plans:
        design_id = int(plan["parts_design_id"])
        construct = construct_by_id[design_id]
        insert = plan.get("insert")
        backbone = plan.get("backbone")
        if not isinstance(insert, Mapping) or not isinstance(backbone, Mapping):
            raise ValueError("最终组装计划缺少 insert 或 backbone")
        if (
            insert.get("file_sha256") != construct.file_sha256
            or insert.get("sequence_sha256") != construct.sequence_sha256
            or backbone.get("file_sha256") != context.backbone.file_sha256
            or backbone.get("sequence_content_sha256")
            != context.backbone.sequence_content_sha256
        ):
            raise ValueError("最终组装计划文件指纹与当前输入不一致")
        _validate_method_payload(plan)
        recorded = str(plan.get("plan_fingerprint") or "")
        unsigned = dict(plan)
        unsigned.pop("plan_fingerprint", None)
        if recorded != stable_json_hash(unsigned):
            raise ValueError(f"design {design_id} 的 plan_fingerprint 无效")
    if artifact.get("plan_set_fingerprint") != stable_json_hash(plans):
        raise ValueError("最终组装计划集合指纹无效")
    requested_method = artifact.get("requested_method")
    if requested_method is not None and any(
        item.get("assembly_method") != requested_method for item in plans
    ):
        raise ValueError("用户指定 method 未应用到全部 design")
    return sorted(plans, key=lambda item: int(item["parts_design_id"]))


def _selection_payload(
    artifact: Mapping[str, Any],
    plans: list[dict[str, Any]],
    *,
    context: Any,
) -> dict[str, Any]:
    base: dict[str, Any] = {
        "schema_version": FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION,
        "status": "selected",
        "selection_mode": artifact.get("selection_mode"),
        "requested_method": artifact.get("requested_method"),
        "design_count": len(plans),
        "design_plans": plans,
        "source": {
            "artifact": (
                "final_assemble_plan/"
                + ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME
            ),
            "artifact_schema_version": artifact.get("schema_version"),
            "algorithm_version": artifact.get("algorithm_version"),
            "plan_set_fingerprint": artifact.get("plan_set_fingerprint"),
            "request_fingerprint": artifact.get("source", {}).get(
                "request_fingerprint"
            ),
            "context_input_fingerprint": context.input_fingerprint,
            "parts_selection_fingerprint": context.parts_selection_fingerprint,
            "assembled_constructs_fingerprint": context.assembled_constructs_fingerprint,
            "plasmid_selection_fingerprint": context.plasmid_selection_fingerprint,
        },
        "warnings": list(artifact.get("warnings") or []),
    }
    base["selection_fingerprint"] = stable_json_hash(base)
    return base


def write_final_assembly_plan(config: Any) -> dict[str, Any]:
    if not bool(getattr(config, "assembly_plan", False)):
        raise ValueError("请使用 write --assembly-plan 接受当前完整组装计划")
    context = load_final_assembly_context(config)
    artifact_path = _artifact_path(config)
    artifact = _read_artifact(artifact_path)
    plans = _validate_artifact(artifact, context=context)
    payload = _selection_payload(artifact, plans, context=context)
    manifest = read_design_manifest(context.manifest_path)
    current = manifest.get("final_assembly_plan")
    unchanged = (
        isinstance(current, Mapping)
        and current.get("schema_version") == FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION
        and current.get("selection_fingerprint")
        == payload["selection_fingerprint"]
    )
    if unchanged:
        updated = manifest
    else:
        updated = update_design_manifest(
            context.manifest_path,
            target_compound_id=context.target_compound_id,
            sections={"final_assembly_plan": payload},
            discard_sections=FINAL_ASSEMBLY_PLAN_DOWNSTREAM_SECTIONS,
            expected_revision=context.manifest_revision,
        )
    return {
        "ok": True,
        "status": "selected",
        "target_compound_id": context.target_compound_id,
        "selection_mode": payload["selection_mode"],
        "requested_method": payload["requested_method"],
        "design_count": len(plans),
        "method_counts": dict(
            sorted(
                {
                    method: sum(
                        item["assembly_method"] == method for item in plans
                    )
                    for method in SUPPORTED_METHODS
                    if any(item["assembly_method"] == method for item in plans)
                }.items()
            )
        ),
        "manifest_path": str(context.manifest_path.resolve()),
        "manifest_revision": updated["revision"],
        "manifest_modified": not unchanged,
    }


def run_write_final_assembly_plan(config: Any) -> dict[str, Any]:
    result = write_final_assembly_plan(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "FINAL_ASSEMBLY_PLAN_DOWNSTREAM_SECTIONS",
    "run_write_final_assembly_plan",
    "write_final_assembly_plan",
]
