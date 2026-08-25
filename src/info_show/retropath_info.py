"""Read-only, integrity-checked views of isolated RetroPath candidates."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from src.main_protein_selection.retropath_enzyme_selection import (
    ENZYME_CANDIDATE_COLUMNS,
    ENZYME_REQUIREMENTS_FILE_NAME,
    ENZYME_SELECTION_DIR_NAME,
    ENZYME_SELECTION_FILE_NAME,
    RETROPATH_ENZYME_SELECTION_SCHEMA,
    SELENZYME_EVIDENCE_FILE_NAME,
    STEP_ENZYME_AUDIT_FILE_NAME,
    STEP_ENZYME_CANDIDATES_FILE_NAME,
)
from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEP_COLUMNS,
    CANDIDATE_STEPS_FILE_NAME,
    REJECTED_ROUTE_COLUMNS,
    REJECTED_ROUTES_FILE_NAME,
)
from src.pathway_analyze.retropath_gem_validation import (
    RETROPATH_GEM_VALIDATION_SCHEMA,
    VALIDATION_FLUX_FILE_NAME,
    VALIDATION_MANIFEST_FILE_NAME,
    VALIDATION_SUMMARY_FILE_NAME,
)
from src.pathway_analyze.retropath_pipeline import (
    PIPELINE_RESULT_FILE_NAME,
    RETROPATH_PIPELINE_SCHEMA,
)
from src.pathway_analyze.target_id import validate_target_compound_id

STATUS_NAMES = {
    "retropath_candidates_found": "已找到 RetroPath 候选路线",
    "retropath_no_scope": "未找到命中可信 sink 的完整候选路线",
    "retropath_source_in_sink": "目标化合物已属于所选 sink",
    "retropath_input_invalid": "RetroPath 输入无效",
    "retropath_expansion_missing": "指定深度的底盘扩展缺失或已过期",
    "retropath_rules_missing": "RetroRules 规则文件缺失",
    "retropath_service_unavailable": "本地 RetroPath 服务不可用",
    "retropath_timeout": "RetroPath 运行超时",
    "retropath_execution_failed": "RetroPath 执行失败",
    "retropath_parse_failed": "RetroPath 结果解析失败",
    "retropath_merge_failed": "RetroPath 候选拼接失败",
    "retropath_configuration_invalid": "RetroPath 配置无效",
}

REJECTION_REASON_NAMES = {
    "candidate_dependency_cycle": "候选步骤依赖形成环路",
    "candidate_enzyme_limit_exceeded": "候选新增酶数量超过上限",
    "candidate_limit_exceeded": "候选数量超过 Top-K 上限",
    "candidate_merge_invalid": "候选路线无法安全拼接",
    "candidate_step_limit_exceeded": "候选总步骤数超过上限",
    "cycle_detected": "逆合成网络包含环路",
    "depth_exceeded": "逆合成分支超过最大深度",
    "enumeration_limit_reached": "路径枚举达到搜索上限",
    "expansion_witness_missing": "缺少从 A0 到 sink 的 KEGG 扩展证据",
    "no_heavy_atom_auxiliary_fragment": "发现无重原子的辅助片段",
    "non_structural_auxiliary_fragment": "发现未结构化的辅助片段",
    "reaction_direction_invalid": "反应方向或迭代层次不一致",
    "reaction_noop": "预测反应没有产生结构变化",
    "rule_evidence_missing": "预测步骤缺少对应 RetroRules 证据",
    "sink_depth_mismatch": "sink 深度与扩展记录不一致",
    "sink_identity_mismatch": "sink 结构身份与 P2 输入不一致",
    "sink_not_in_expansion_bundle": "sink 不在所选累计可达集合中",
    "structure_invalid": "预测结构无法有效解析",
    "top_k_pruned": "路线在 Top-K 排序中被裁剪",
    "unresolved_non_sink_leaf": "分支末端既不是 sink，也无法继续展开",
}

VALIDATION_STATUS_NAMES = {
    "raw": "原始预测候选，尚未验证",
    "validated": "已验证",
    "rejected": "验证未通过",
    "promoted": "已晋升为正式候选",
}

STEP_SOURCE_NAMES = {
    "kegg_expansion": "KEGG 扩展前缀",
    "retropath": "RetroPath 预测步骤",
}

STEP_STATUS_NAMES = {
    "endogenous": "内源反应",
    "heterologous": "异源 KEGG 反应",
    "predicted": "预测反应",
}

DIRECTION_NAMES = {
    "left_to_right": "从左到右",
    "right_to_left": "从右到左",
}

ORIENTATION_NAMES = {
    "biosynthetic": "合成方向",
    "retrosynthetic": "逆合成方向",
}

SCORE_SEMANTICS_NAMES = {
    "lower_is_better": "分数越低越好",
    "higher_is_better": "分数越高越好",
}


@dataclass(frozen=True)
class _RetroPathViewContext:
    target_compound: str
    depth: int
    output_dir: Path
    pipeline: Mapping[str, Any]
    candidate_routes: tuple[Mapping[str, str], ...] = tuple()
    candidate_steps: tuple[Mapping[str, str], ...] = tuple()
    rejected_routes: tuple[Mapping[str, str], ...] = tuple()


@dataclass(frozen=True)
class _P8Overlay:
    state: str
    candidate_statuses: Mapping[int, str]
    passing_combinations: Mapping[int, tuple[str, ...]]
    flux_by_candidate_step: Mapping[tuple[int, str], float]
    warning: str = ""


@dataclass(frozen=True)
class _P9Overlay:
    state: str
    status: str
    recommended_combination_id: str
    combinations: tuple[Mapping[str, Any], ...]
    candidates_by_step: Mapping[int, tuple[Mapping[str, str], ...]]
    warning: str = ""


def _rerun_hint(depth: int) -> str:
    return f"请重新运行 gap --input <输入文件> --depth {depth} --retropath"


def _read_json(path: Path, depth: int) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            f"未找到 RetroPath 运行结果：{path}；{_rerun_hint(depth)}"
        )
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"RetroPath 运行结果不是有效 JSON：{path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"RetroPath 运行结果根节点必须是对象：{path}")
    return value


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _required_mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"RetroPath 结果字段 {field_name} 必须是对象")
    return value


def _as_int(value: Any, field_name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ValueError(f"RetroPath 结果字段 {field_name} 必须是整数")
    try:
        normalized = int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"RetroPath 结果字段 {field_name} 必须是整数") from exc
    if normalized < minimum:
        raise ValueError(
            f"RetroPath 结果字段 {field_name} 必须大于等于 {minimum}"
        )
    return normalized


def _optional_float(value: Any, field_name: str) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        normalized = float(text)
    except ValueError as exc:
        raise ValueError(f"RetroPath 结果字段 {field_name} 必须是数字") from exc
    if not math.isfinite(normalized):
        raise ValueError(f"RetroPath 结果字段 {field_name} 必须是有限数字")
    return normalized


def _as_bool(value: Any, field_name: str) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ValueError(f"RetroPath 结果字段 {field_name} 必须是布尔值")


def _optional_bool(value: Any, field_name: str) -> bool | None:
    if str(value or "").strip() == "":
        return None
    return _as_bool(value, field_name)


def _split(value: Any) -> list[str]:
    return [item.strip() for item in str(value or "").split(";") if item.strip()]


def _read_csv(
    path: Path,
    expected_columns: Sequence[str],
    depth: int,
) -> tuple[Mapping[str, str], ...]:
    if not path.is_file():
        raise FileNotFoundError(
            f"缺少 RetroPath 候选文件：{path}；{_rerun_hint(depth)}"
        )
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            observed_columns = tuple(reader.fieldnames or tuple())
            if observed_columns != tuple(expected_columns):
                raise ValueError(
                    f"RetroPath 候选文件表头不兼容：{path}；"
                    f"期望 {list(expected_columns)}，实际 {list(observed_columns)}"
                )
            return tuple(dict(row) for row in reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"无法读取 RetroPath 候选文件：{path}") from exc


def _verify_artifact(
    *,
    output_dir: Path,
    artifacts: Mapping[str, Any],
    artifact_key: str,
    file_name: str,
    depth: int,
) -> Path:
    record = _required_mapping(artifacts.get(artifact_key), f"artifacts.{artifact_key}")
    expected_path = (output_dir / file_name).resolve()
    recorded_path_text = str(record.get("path") or "").strip()
    if not recorded_path_text:
        raise ValueError(f"RetroPath 结果缺少 artifacts.{artifact_key}.path")
    if Path(recorded_path_text).expanduser().resolve() != expected_path:
        raise ValueError(
            f"RetroPath 结果记录的 {artifact_key} 路径与当前 depth 不一致"
        )
    recorded_hash = str(record.get("sha256") or "").strip().lower()
    if len(recorded_hash) != 64 or any(
        character not in "0123456789abcdef" for character in recorded_hash
    ):
        raise ValueError(f"RetroPath 结果缺少有效的 {artifact_key} SHA-256")
    if not expected_path.is_file():
        raise FileNotFoundError(
            f"缺少 RetroPath 候选文件：{expected_path}；{_rerun_hint(depth)}"
        )
    observed_hash = _sha256_file(expected_path)
    if observed_hash != recorded_hash:
        raise ValueError(
            f"RetroPath 候选文件校验失败，文件可能已变化：{expected_path}；"
            f"{_rerun_hint(depth)}"
        )
    return expected_path


def _validate_candidate_relationships(
    *,
    target_compound: str,
    routes: Sequence[Mapping[str, str]],
    steps: Sequence[Mapping[str, str]],
) -> None:
    route_ids: dict[str, Mapping[str, str]] = {}
    ranks: list[int] = []
    for route in routes:
        rank = _as_int(route.get("candidate_rank"), "candidate_rank", minimum=1)
        ranks.append(rank)
        candidate_id = str(route.get("candidate_id") or "").strip()
        if not candidate_id or candidate_id in route_ids:
            raise ValueError("RetroPath 候选 ID 为空或重复")
        recorded_target = str(route.get("target_compound_id") or "").strip().upper()
        if recorded_target != target_compound:
            raise ValueError(
                f"RetroPath 候选目标不一致：{recorded_target} != {target_compound}"
            )
        route_ids[candidate_id] = route
    if ranks != list(range(1, len(routes) + 1)):
        raise ValueError("RetroPath candidate_rank 必须从 1 连续递增")

    steps_by_candidate: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for step in steps:
        candidate_id = str(step.get("candidate_id") or "").strip()
        if candidate_id not in route_ids:
            raise ValueError(f"RetroPath 步骤引用未知候选 ID：{candidate_id}")
        steps_by_candidate[candidate_id].append(step)

    for candidate_id, route in route_ids.items():
        candidate_steps = sorted(
            steps_by_candidate.get(candidate_id, []),
            key=lambda row: _as_int(row.get("step_index"), "step_index", minimum=1),
        )
        total_steps = _as_int(route.get("total_steps"), "total_steps", minimum=1)
        if len(candidate_steps) != total_steps:
            raise ValueError(
                f"RetroPath 候选 {candidate_id} 的步骤数与摘要不一致"
            )
        indexes = [
            _as_int(row.get("step_index"), "step_index", minimum=1)
            for row in candidate_steps
        ]
        if indexes != list(range(1, total_steps + 1)):
            raise ValueError(
                f"RetroPath 候选 {candidate_id} 的 step_index 不连续"
            )
        step_index_by_id: dict[str, int] = {}
        for index, row in zip(indexes, candidate_steps):
            step_id = str(row.get("step_id") or "").strip()
            if not step_id or step_id in step_index_by_id:
                raise ValueError(
                    f"RetroPath 候选 {candidate_id} 的步骤 ID 为空或重复"
                )
            step_index_by_id[step_id] = index
        for index, row in zip(indexes, candidate_steps):
            for dependency_id in _split(row.get("depends_on_step_ids")):
                dependency_index = step_index_by_id.get(dependency_id)
                if dependency_index is None or dependency_index >= index:
                    raise ValueError(
                        f"RetroPath 候选 {candidate_id} 的步骤依赖不是有效拓扑顺序"
                    )


def _load_context(config: Any) -> _RetroPathViewContext:
    target_compound = validate_target_compound_id(config.target_name)
    depth = _as_int(getattr(config, "depth", 0), "depth")
    output_dir = (
        gap_depth_output_dir(Path(config.gap_output_path), depth) / "retropath"
    ).expanduser().resolve()
    pipeline = _read_json(output_dir / PIPELINE_RESULT_FILE_NAME, depth)
    if pipeline.get("schema_version") != RETROPATH_PIPELINE_SCHEMA:
        raise ValueError(
            "RetroPath 运行结果版本不兼容；"
            f"{_rerun_hint(depth)}"
        )
    if pipeline.get("retropath_requested") is not True or (
        pipeline.get("search_engine") != "retropath"
    ):
        raise ValueError("指定文件不是 RetroPath pipeline 运行结果")
    recorded_target = str(pipeline.get("target_compound") or "").strip().upper()
    if recorded_target != target_compound:
        raise ValueError(
            f"RetroPath 运行目标与当前输入不一致："
            f"{recorded_target} != {target_compound}"
        )
    recorded_depth = _as_int(pipeline.get("expansion_depth"), "expansion_depth")
    if recorded_depth != depth:
        raise ValueError(
            f"RetroPath 运行深度与请求深度不一致：{recorded_depth} != {depth}"
        )
    recorded_output_dir = str(pipeline.get("output_dir") or "").strip()
    if recorded_output_dir and (
        Path(recorded_output_dir).expanduser().resolve() != output_dir
    ):
        raise ValueError("RetroPath 运行结果的输出目录与当前项目不一致")
    if not isinstance(pipeline.get("ok"), bool):
        raise ValueError("RetroPath 运行结果字段 ok 必须是布尔值")
    status = str(pipeline.get("status") or "").strip()
    if status not in STATUS_NAMES:
        raise ValueError(f"RetroPath 运行结果包含未知状态：{status!r}")
    successful_statuses = {
        "retropath_candidates_found",
        "retropath_no_scope",
        "retropath_source_in_sink",
    }
    if (pipeline["ok"] is True) != (status in successful_statuses):
        raise ValueError("RetroPath 运行结果的 ok 与 status 不一致")
    expected_sink_source = (
        "chassis_A0" if depth == 0 else f"cumulative_expansion_A{depth}"
    )
    if pipeline.get("sink_source") != expected_sink_source:
        raise ValueError("RetroPath 运行结果的 sink_source 与 depth 不一致")
    if pipeline["ok"] is False:
        return _RetroPathViewContext(
            target_compound=target_compound,
            depth=depth,
            output_dir=output_dir,
            pipeline=pipeline,
        )

    artifacts = _required_mapping(pipeline.get("artifacts"), "artifacts")
    route_path = _verify_artifact(
        output_dir=output_dir,
        artifacts=artifacts,
        artifact_key="candidate_routes",
        file_name=CANDIDATE_ROUTES_FILE_NAME,
        depth=depth,
    )
    step_path = _verify_artifact(
        output_dir=output_dir,
        artifacts=artifacts,
        artifact_key="candidate_steps",
        file_name=CANDIDATE_STEPS_FILE_NAME,
        depth=depth,
    )
    rejected_path = _verify_artifact(
        output_dir=output_dir,
        artifacts=artifacts,
        artifact_key="rejected_routes",
        file_name=REJECTED_ROUTES_FILE_NAME,
        depth=depth,
    )
    candidate_routes = _read_csv(route_path, CANDIDATE_ROUTE_COLUMNS, depth)
    candidate_steps = _read_csv(step_path, CANDIDATE_STEP_COLUMNS, depth)
    rejected_routes = _read_csv(rejected_path, REJECTED_ROUTE_COLUMNS, depth)
    candidate_count = _as_int(pipeline.get("candidate_count"), "candidate_count")
    rejection_count = _as_int(pipeline.get("rejection_count"), "rejection_count")
    if candidate_count != len(candidate_routes):
        raise ValueError("RetroPath pipeline 候选数量与 candidate_routes.csv 不一致")
    if rejection_count != len(rejected_routes):
        raise ValueError("RetroPath pipeline 拒绝数量与 rejected_routes.csv 不一致")
    if (status == "retropath_candidates_found") != bool(candidate_routes):
        raise ValueError("RetroPath 运行状态与候选文件是否包含路线不一致")
    _validate_candidate_relationships(
        target_compound=target_compound,
        routes=candidate_routes,
        steps=candidate_steps,
    )
    return _RetroPathViewContext(
        target_compound=target_compound,
        depth=depth,
        output_dir=output_dir,
        pipeline=pipeline,
        candidate_routes=candidate_routes,
        candidate_steps=candidate_steps,
        rejected_routes=rejected_routes,
    )


def _sink_boundaries(row: Mapping[str, str]) -> list[dict[str, Any]]:
    sink_ids = _split(row.get("sink_kegg_ids"))
    depths: dict[str, int] = {}
    for item in _split(row.get("sink_depths")):
        compound_id, separator, raw_depth = item.rpartition(":")
        if not separator or not compound_id:
            raise ValueError(f"RetroPath sink_depths 格式无效：{item!r}")
        depths[compound_id] = _as_int(raw_depth, "sink_depths")
    if set(sink_ids) != set(depths):
        raise ValueError("RetroPath sink_kegg_ids 与 sink_depths 不一致")
    return [
        {"KEGG化合物ID": compound_id, "首次可达深度": depths[compound_id]}
        for compound_id in sink_ids
    ]


def _route_summary(row: Mapping[str, str]) -> dict[str, Any]:
    score_semantics = str(row.get("score_semantics") or "").strip()
    validation_status = str(row.get("validation_status") or "").strip()
    return {
        "候选排名": _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1),
        "候选ID": str(row.get("candidate_id") or "").strip(),
        "原始逆合成路径ID": str(
            row.get("source_retrosynthetic_path_id") or ""
        ).strip(),
        "命中边界化合物": _sink_boundaries(row),
        "KEGG前缀反应": _split(row.get("kegg_prefix_reaction_ids")),
        "RetroPath预测步骤ID": _split(row.get("retropath_step_ids")),
        "RetroPath反应选项": _split(row.get("retropath_reaction_option_ids")),
        "KEGG前缀步骤数": _as_int(row.get("kegg_prefix_steps"), "kegg_prefix_steps"),
        "RetroPath预测步骤数": _as_int(row.get("retropath_steps"), "retropath_steps"),
        "总步骤数": _as_int(row.get("total_steps"), "total_steps", minimum=1),
        "最大Sink深度": _as_int(row.get("maximum_sink_depth"), "maximum_sink_depth"),
        "最小规则特异性": _as_int(
            row.get("minimum_rule_specificity"),
            "minimum_rule_specificity",
        ),
        "最差规则分数": _optional_float(row.get("worst_rule_score"), "worst_rule_score"),
        "分数含义": SCORE_SEMANTICS_NAMES.get(score_semantics, score_semantics),
        "包含辅助片段": _as_bool(
            row.get("contains_auxiliary_fragments"),
            "contains_auxiliary_fragments",
        ),
        "路线来源": str(row.get("route_source") or "").strip(),
        "包含预测步骤": _as_bool(
            row.get("contains_predicted_steps"),
            "contains_predicted_steps",
        ),
        "验证状态": VALIDATION_STATUS_NAMES.get(validation_status, validation_status),
        "需要人工复核": _as_bool(row.get("review_required"), "review_required"),
        "上游枚举被截断": _as_bool(
            row.get("upstream_enumeration_truncated"),
            "upstream_enumeration_truncated",
        ),
        "候选TopK被截断": _as_bool(
            row.get("candidate_top_k_truncated"),
            "candidate_top_k_truncated",
        ),
    }


def _stoichiometry(value: Any, field_name: str) -> list[dict[str, Any]]:
    try:
        parsed = json.loads(str(value or ""))
    except json.JSONDecodeError as exc:
        raise ValueError(f"RetroPath {field_name} 不是有效 JSON") from exc
    if not isinstance(parsed, list):
        raise ValueError(f"RetroPath {field_name} 必须是列表")
    result: list[dict[str, Any]] = []
    for item in parsed:
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError(f"RetroPath {field_name} 项格式无效")
        compound_id = str(item[0] or "").strip()
        coefficient = _optional_float(item[1], field_name)
        if not compound_id or coefficient is None:
            raise ValueError(f"RetroPath {field_name} 项格式无效")
        result.append({"化合物ID": compound_id, "计量系数": coefficient})
    return result


def _step_display(row: Mapping[str, str]) -> dict[str, Any]:
    step_source = str(row.get("step_source") or "").strip()
    status = str(row.get("status") or "").strip()
    direction = str(row.get("direction") or "").strip()
    orientation = str(row.get("orientation") or "").strip()
    score_semantics = str(row.get("score_semantics") or "").strip()
    balance_status = str(row.get("balance_status") or "").strip()
    cofactor_status = str(row.get("cofactor_reconstruction_status") or "").strip()
    warnings: list[str] = []
    if step_source == "retropath":
        warnings.append("该步骤为 RetroPath 预测反应，需结合 P8/P9 结果和人工复核")
    if balance_status != "balanced":
        warnings.append(f"反应平衡状态：{balance_status or '未知'}")
    if cofactor_status != "complete":
        warnings.append(f"辅因子恢复状态：{cofactor_status or '未知'}")
    return {
        "步骤编号": _as_int(row.get("step_index"), "step_index", minimum=1),
        "步骤ID": str(row.get("step_id") or "").strip(),
        "步骤来源": STEP_SOURCE_NAMES.get(step_source, step_source),
        "反应类型": STEP_STATUS_NAMES.get(status, status),
        "方向语义": ORIENTATION_NAMES.get(orientation, orientation),
        "方向": DIRECTION_NAMES.get(direction, direction),
        "反应选项ID": _split(row.get("reaction_option_ids")),
        "Reaction SMILES": str(row.get("reaction_smiles") or "").strip(),
        "底物计量": _stoichiometry(
            row.get("substrate_stoichiometry_json"),
            "substrate_stoichiometry_json",
        ),
        "产物计量": _stoichiometry(
            row.get("product_stoichiometry_json"),
            "product_stoichiometry_json",
        ),
        "依赖步骤ID": _split(row.get("depends_on_step_ids")),
        "原始Transformation ID": _split(row.get("source_transformation_ids")),
        "Sink锚点": _split(row.get("sink_anchor_kegg_ids")),
        "扩展深度": _as_int(row.get("expansion_depth"), "expansion_depth"),
        "是否内源": _optional_bool(row.get("is_endogenous"), "is_endogenous"),
        "RetroRules规则ID": _split(row.get("rule_ids")),
        "来源反应ID": _split(row.get("source_reaction_ids")),
        "来源EC编号": _split(row.get("source_ec_numbers")),
        "来源UniProt": _split(row.get("source_uniprot_ids")),
        "最小规则特异性": (
            None
            if str(row.get("minimum_rule_specificity") or "").strip() == ""
            else _as_int(
                row.get("minimum_rule_specificity"),
                "minimum_rule_specificity",
            )
        ),
        "最差规则分数": _optional_float(row.get("worst_rule_score"), "worst_rule_score"),
        "分数含义": SCORE_SEMANTICS_NAMES.get(score_semantics, score_semantics),
        "平衡状态": balance_status,
        "辅因子恢复状态": cofactor_status,
        "风险提示": warnings,
    }


def _rejection_summary(
    rows: Sequence[Mapping[str, str]],
) -> list[dict[str, Any]]:
    counts: Counter[tuple[str, str]] = Counter()
    examples: dict[tuple[str, str], str] = {}
    for row in rows:
        stage = str(row.get("source_stage") or "").strip()
        code = str(row.get("reason_code") or "").strip()
        if not stage or not code:
            raise ValueError("RetroPath 拒绝记录缺少 source_stage 或 reason_code")
        key = (stage, code)
        counts[key] += 1
        detail = str(row.get("reason_detail") or "").strip()
        if detail and key not in examples:
            examples[key] = detail
    return [
        {
            "阶段": stage,
            "原因代码": code,
            "原因说明": REJECTION_REASON_NAMES.get(code, "未识别的稳定原因代码"),
            "数量": counts[(stage, code)],
            "示例详情": examples.get((stage, code), ""),
        }
        for stage, code in sorted(counts)
    ]


def _unique(values: Sequence[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _success_warnings(context: _RetroPathViewContext) -> list[str]:
    warnings: list[str] = []
    if context.candidate_routes:
        warnings.append(
            "候选含 RetroPath 预测步骤，尚未完成计量、GEM 和酶证据验证，"
            "不可直接进入正式设计"
        )
    if any(
        _as_bool(row.get("contains_auxiliary_fragments"), "contains_auxiliary_fragments")
        for row in context.candidate_routes
    ):
        warnings.append("部分候选包含辅助片段，共底物或辅因子仍需恢复")
    if any(
        _as_bool(
            row.get("upstream_enumeration_truncated"),
            "upstream_enumeration_truncated",
        )
        for row in context.candidate_routes
    ) or _as_bool(
        context.pipeline.get("upstream_enumeration_truncated"),
        "upstream_enumeration_truncated",
    ):
        warnings.append("RetroPath 完整路径枚举达到上限，结果并非穷尽")
    if any(
        _as_bool(
            row.get("candidate_top_k_truncated"),
            "candidate_top_k_truncated",
        )
        for row in context.candidate_routes
    ) or _as_bool(
        context.pipeline.get("candidate_top_k_truncated"),
        "candidate_top_k_truncated",
    ):
        warnings.append("候选经过 Top-K 裁剪，未展示全部可拼接路线")
    status = str(context.pipeline.get("status") or "").strip()
    if status == "retropath_source_in_sink":
        warnings.append("目标化合物已属于所选 A0/AN sink，不需要 RetroPath 预测步骤")
    elif status == "retropath_no_scope":
        warnings.append("本次运行没有得到命中可信 sink 的完整 RetroPath 候选")
    return _unique(warnings)


def _p8_overlay(context: _RetroPathViewContext) -> _P8Overlay:
    validation_dir = context.output_dir / "gem_validation"
    manifest_path = validation_dir / VALIDATION_MANIFEST_FILE_NAME
    if not manifest_path.is_file():
        return _P8Overlay(
            "not_run",
            {},
            {},
            {},
            "尚未运行 P8 严格计量和 GEM 验证",
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        if (
            not isinstance(manifest, dict)
            or manifest.get("schema_version") != RETROPATH_GEM_VALIDATION_SCHEMA
            or manifest.get("target_compound") != context.target_compound
            or _as_int(manifest.get("expansion_depth"), "P8 expansion_depth")
            != context.depth
        ):
            raise ValueError("P8 manifest identity mismatch")
        inputs = _required_mapping(manifest.get("inputs"), "P8 inputs")
        pipeline_record = _required_mapping(
            inputs.get("pipeline_result"),
            "P8 inputs.pipeline_result",
        )
        if pipeline_record.get("sha256") != _sha256_file(
            context.output_dir / PIPELINE_RESULT_FILE_NAME
        ):
            raise ValueError("P8 input pipeline hash is stale")
        pipeline_artifacts = _required_mapping(
            context.pipeline.get("artifacts"),
            "artifacts",
        )
        if inputs.get("candidate_routes_sha256") != _required_mapping(
            pipeline_artifacts.get("candidate_routes"),
            "artifacts.candidate_routes",
        ).get("sha256") or inputs.get("candidate_steps_sha256") != _required_mapping(
            pipeline_artifacts.get("candidate_steps"),
            "artifacts.candidate_steps",
        ).get("sha256"):
            raise ValueError("P8 candidate input hashes are stale")
        artifacts = _required_mapping(manifest.get("artifacts"), "P8 artifacts")

        def artifact_path(name: str) -> Path:
            record = _required_mapping(artifacts.get(name), f"P8 artifacts.{name}")
            expected = (validation_dir / name).resolve()
            if Path(str(record.get("path") or "")).expanduser().resolve() != expected:
                raise ValueError(f"P8 artifact path is stale: {name}")
            if not expected.is_file() or record.get("sha256") != _sha256_file(expected):
                raise ValueError(f"P8 artifact checksum mismatch: {name}")
            return expected

        summary_path = artifact_path(VALIDATION_SUMMARY_FILE_NAME)
        flux_path = artifact_path(VALIDATION_FLUX_FILE_NAME)
        with summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
            summary_rows = tuple(dict(row) for row in csv.DictReader(handle))
        with flux_path.open("r", encoding="utf-8-sig", newline="") as handle:
            flux_rows = tuple(dict(row) for row in csv.DictReader(handle))
        raw_statuses = _required_mapping(
            manifest.get("candidate_statuses"),
            "P8 candidate_statuses",
        )
        statuses = {
            _as_int(rank, "P8 candidate rank", minimum=1): str(status)
            for rank, status in raw_statuses.items()
        }
        passing: dict[int, tuple[str, ...]] = {}
        for rank in statuses:
            passing[rank] = tuple(
                sorted(
                    {
                        str(row.get("combination_id") or "")
                        for row in summary_rows
                        if _as_int(
                            row.get("candidate_rank"),
                            "P8 candidate_rank",
                            minimum=1,
                        )
                        == rank
                        and row.get("validation_status")
                        == "PASS_STRICT_ROUTE_FLUX"
                        and str(row.get("combination_id") or "")
                    }
                )
            )
        flux_by_step: dict[tuple[int, str], float] = {}
        for rank, combinations in passing.items():
            if not combinations:
                continue
            selected_combination = combinations[0]
            for row in flux_rows:
                if (
                    _as_int(
                        row.get("candidate_rank"),
                        "P8 flux candidate_rank",
                        minimum=1,
                    )
                    == rank
                    and row.get("combination_id") == selected_combination
                ):
                    flux_by_step[(rank, str(row.get("step_id") or ""))] = float(
                        row.get("directed_fba_flux") or "nan"
                    )
        return _P8Overlay("current", statuses, passing, flux_by_step)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        return _P8Overlay(
            "stale",
            {},
            {},
            {},
            f"P8 验证结果已过期或损坏，未应用其状态：{exc}",
        )


def _p9_overlay(
    context: _RetroPathViewContext,
    rank: int,
    candidate_id: str,
) -> _P9Overlay:
    output_dir = (
        context.output_dir
        / ENZYME_SELECTION_DIR_NAME
        / f"candidate_{rank}"
    )
    selection_path = output_dir / ENZYME_SELECTION_FILE_NAME
    if not selection_path.is_file():
        return _P9Overlay(
            "not_run",
            "not_run",
            "",
            tuple(),
            {},
            "尚未运行 P9 RetroPath 主酶候选检索",
        )
    try:
        selection = json.loads(selection_path.read_text(encoding="utf-8-sig"))
        if (
            not isinstance(selection, dict)
            or selection.get("schema_version") != RETROPATH_ENZYME_SELECTION_SCHEMA
            or selection.get("target_compound") != context.target_compound
            or _as_int(selection.get("expansion_depth"), "P9 expansion_depth")
            != context.depth
            or _as_int(selection.get("candidate_rank"), "P9 candidate_rank", minimum=1)
            != rank
            or selection.get("candidate_id") != candidate_id
        ):
            raise ValueError("P9 manifest identity mismatch")
        inputs = _required_mapping(selection.get("inputs"), "P9 inputs")
        pipeline_record = _required_mapping(
            inputs.get("pipeline_result"),
            "P9 inputs.pipeline_result",
        )
        if pipeline_record.get("sha256") != _sha256_file(
            context.output_dir / PIPELINE_RESULT_FILE_NAME
        ):
            raise ValueError("P9 pipeline input hash is stale")
        p8_record = _required_mapping(
            inputs.get("p8_validation_manifest"),
            "P9 inputs.p8_validation_manifest",
        )
        expected_p8 = (
            context.output_dir / "gem_validation" / VALIDATION_MANIFEST_FILE_NAME
        ).resolve()
        if (
            Path(str(p8_record.get("path") or "")).expanduser().resolve()
            != expected_p8
            or not expected_p8.is_file()
            or p8_record.get("sha256") != _sha256_file(expected_p8)
        ):
            raise ValueError("P9 P8-validation input hash is stale")
        pipeline_artifacts = _required_mapping(
            context.pipeline.get("artifacts"),
            "pipeline artifacts",
        )
        if inputs.get("candidate_routes_sha256") != _required_mapping(
            pipeline_artifacts.get("candidate_routes"),
            "pipeline artifacts.candidate_routes",
        ).get("sha256") or inputs.get("candidate_steps_sha256") != _required_mapping(
            pipeline_artifacts.get("candidate_steps"),
            "pipeline artifacts.candidate_steps",
        ).get("sha256"):
            raise ValueError("P9 candidate input hashes are stale")
        artifacts = _required_mapping(selection.get("artifacts"), "P9 artifacts")
        verified_paths: dict[str, Path] = {}
        for name in (
            ENZYME_REQUIREMENTS_FILE_NAME,
            STEP_ENZYME_CANDIDATES_FILE_NAME,
            STEP_ENZYME_AUDIT_FILE_NAME,
            SELENZYME_EVIDENCE_FILE_NAME,
        ):
            record = _required_mapping(
                artifacts.get(name),
                f"P9 artifacts.{name}",
            )
            expected = (output_dir / name).resolve()
            if (
                Path(str(record.get("path") or "")).expanduser().resolve()
                != expected
                or not expected.is_file()
                or record.get("sha256") != _sha256_file(expected)
            ):
                raise ValueError(f"P9 artifact checksum mismatch: {name}")
            verified_paths[name] = expected
        candidate_path = verified_paths[STEP_ENZYME_CANDIDATES_FILE_NAME]
        with candidate_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != ENZYME_CANDIDATE_COLUMNS:
                raise ValueError("P9 candidate artifact schema mismatch")
            candidate_rows = tuple(dict(row) for row in reader)
        by_step: dict[int, list[Mapping[str, str]]] = defaultdict(list)
        for row in candidate_rows:
            if (
                _as_int(row.get("candidate_rank"), "P9 row candidate_rank", minimum=1)
                != rank
                or row.get("candidate_id") != candidate_id
                or row.get("selection_status") != "selected"
            ):
                raise ValueError("P9 candidate row binding mismatch")
            by_step[_as_int(row.get("step_index"), "P9 step_index", minimum=1)].append(
                row
            )
        combinations_value = selection.get("combinations")
        if not isinstance(combinations_value, list) or any(
            not isinstance(item, Mapping) for item in combinations_value
        ):
            raise ValueError("P9 combinations are invalid")
        combinations = tuple(dict(item) for item in combinations_value)
        combination_ids = {
            str(item.get("combination_id") or "") for item in combinations
        }
        if any(row.get("combination_id") not in combination_ids for row in candidate_rows):
            raise ValueError("P9 candidate references an unknown combination")
        status = str(selection.get("status") or "").strip()
        if status not in {
            "ready_for_review",
            "partial_no_candidate",
            "source_unavailable",
        }:
            raise ValueError("P9 status is invalid")
        if (
            selection.get("review_required") is not True
            or selection.get("formal_promotion_allowed") is not False
        ):
            raise ValueError("P9 promotion/review gate is invalid")
        return _P9Overlay(
            "current",
            status,
            str(selection.get("recommended_combination_id") or ""),
            combinations,
            {
                step_index: tuple(sorted(
                    rows,
                    key=lambda row: (
                        _as_int(
                            row.get("combination_rank"),
                            "P9 combination_rank",
                            minimum=1,
                        ),
                        _as_int(
                            row.get("protein_candidate_rank"),
                            "P9 protein_candidate_rank",
                            minimum=1,
                        ),
                    ),
                ))
                for step_index, rows in by_step.items()
            },
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        return _P9Overlay(
            "stale",
            "stale",
            "",
            tuple(),
            {},
            f"P9 主酶候选结果已过期或损坏，未应用：{exc}",
        )


def _p9_candidate_display(row: Mapping[str, str]) -> dict[str, Any]:
    return {
        "计量组合排名": _as_int(
            row.get("combination_rank"),
            "P9 combination_rank",
            minimum=1,
        ),
        "计量组合ID": str(row.get("combination_id") or ""),
        "计量假设ID": str(row.get("hypothesis_id") or ""),
        "蛋白候选排名": _as_int(
            row.get("protein_candidate_rank"),
            "P9 protein_candidate_rank",
            minimum=1,
        ),
        "UniProt": str(row.get("accession") or ""),
        "蛋白名称": str(row.get("protein_name") or ""),
        "来源物种": str(row.get("organism_name") or ""),
        "Reviewed": _as_bool(row.get("reviewed"), "P9 reviewed"),
        "证据等级": str(row.get("evidence_tier") or ""),
        "反应匹配状态": str(row.get("fit_status") or ""),
        "必须人工复核": _as_bool(
            row.get("manual_review_required"),
            "P9 manual_review_required",
        ),
        "Reaction similarity": str(row.get("reaction_similarity") or ""),
        "sim_RF": str(row.get("sim_rf") or ""),
        "匹配反应": str(row.get("matched_reaction_id") or ""),
        "EC": _split(row.get("ec_numbers")),
        "Rhea": _split(row.get("matched_rhea_ids")),
        "警告": _split(str(row.get("warnings") or "").replace(" | ", ";")),
    }


def get_retropath_info(config: Any) -> dict[str, Any]:
    """Return a compact Chinese summary of one P6 RetroPath run."""

    context = _load_context(config)
    pipeline = context.pipeline
    status = str(pipeline.get("status") or "").strip()
    common: dict[str, Any] = {
        "运行成功": bool(pipeline["ok"]),
        "搜索引擎": "RetroPath",
        "目标化合物ID": context.target_compound,
        "扩展深度": context.depth,
        "Sink来源": str(pipeline.get("sink_source") or "").strip(),
        "运行状态": STATUS_NAMES.get(status, status),
        "原始状态代码": status,
    }
    if pipeline["ok"] is False:
        detail = str(pipeline.get("detail") or "").strip()
        stage = str(pipeline.get("stage") or "").strip()
        return {
            **common,
            "失败阶段": stage,
            "失败详情": detail,
            "候选路线数": 0,
            "警告": [detail] if detail else [],
        }

    input_summary = _required_mapping(pipeline.get("input_summary"), "input_summary")
    overlay = _p8_overlay(context)
    route_summaries = [_route_summary(row) for row in context.candidate_routes]
    p9_warnings: list[str] = []
    for route_row, route_summary in zip(context.candidate_routes, route_summaries):
        rank = int(route_summary["候选排名"])
        route_summary["P8严格GEM状态"] = overlay.candidate_statuses.get(
            rank,
            "未验证",
        )
        route_summary["P8可行计量假设数"] = len(
            overlay.passing_combinations.get(rank, tuple())
        )
        p9 = _p9_overlay(
            context,
            rank,
            str(route_row.get("candidate_id") or "").strip(),
        )
        route_summary["P9主酶候选状态"] = p9.status
        route_summary["P9候选计量组合数"] = len(p9.combinations)
        if p9.state == "stale" and p9.warning:
            p9_warnings.append(p9.warning)
    return {
        **common,
        "RetroPath任务ID": str(pipeline.get("job_id") or "").strip(),
        "服务终态": str(pipeline.get("service_status") or "").strip(),
        "使用缓存": _as_bool(pipeline.get("cache_hit"), "cache_hit"),
        "Scope文件存在": _as_bool(pipeline.get("scope_present"), "scope_present"),
        "累计可达化合物数": _as_int(
            input_summary.get("reachable_compound_count"),
            "input_summary.reachable_compound_count",
        ),
        "有效Sink结构数": _as_int(
            input_summary.get("sink_structure_count"),
            "input_summary.sink_structure_count",
        ),
        "输入结构拒绝数": _as_int(
            input_summary.get("rejected_compound_count"),
            "input_summary.rejected_compound_count",
        ),
        "命中Sink数": _as_int(pipeline.get("sink_match_count"), "sink_match_count"),
        "完整逆向路径数": _as_int(
            pipeline.get("complete_path_count"),
            "complete_path_count",
        ),
        "候选路线数": len(context.candidate_routes),
        "拒绝路线数": len(context.rejected_routes),
        "P8验证状态": overlay.state,
        "候选路线摘要": route_summaries,
        "拒绝原因统计": _rejection_summary(context.rejected_routes),
        "警告": _unique([
            *_success_warnings(context),
            overlay.warning,
            *p9_warnings,
        ]),
    }


def get_retropath_candidate_info(config: Any) -> dict[str, Any]:
    """Return one candidate DAG, or one topologically ordered candidate step."""

    context = _load_context(config)
    if context.pipeline["ok"] is False:
        raise ValueError(
            "RetroPath 运行未成功，不能查看候选；请先使用 info --retropath "
            "查看失败阶段和详情"
        )
    raw_rank = getattr(config, "retropath_candidate", None)
    if raw_rank is None:
        raise ValueError("未指定候选排名，请使用 --retropath-candidate N")
    rank = _as_int(raw_rank, "retropath_candidate", minimum=1)
    route = next(
        (
            row
            for row in context.candidate_routes
            if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1)
            == rank
        ),
        None,
    )
    if route is None:
        available = [
            _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1)
            for row in context.candidate_routes
        ]
        raise ValueError(f"不存在 RetroPath 候选 {rank}；可用候选排名：{available}")
    candidate_id = str(route.get("candidate_id") or "").strip()
    candidate_steps = sorted(
        (
            row
            for row in context.candidate_steps
            if str(row.get("candidate_id") or "").strip() == candidate_id
        ),
        key=lambda row: _as_int(row.get("step_index"), "step_index", minimum=1),
    )
    common = {
        "运行成功": True,
        "搜索引擎": "RetroPath",
        "目标化合物ID": context.target_compound,
        "扩展深度": context.depth,
        "候选排名": rank,
        "候选ID": candidate_id,
    }
    overlay = _p8_overlay(context)
    common["P8验证状态"] = overlay.state
    common["P8严格GEM状态"] = overlay.candidate_statuses.get(rank, "未验证")
    common["P8可行计量假设"] = list(
        overlay.passing_combinations.get(rank, tuple())
    )
    p9 = _p9_overlay(context, rank, candidate_id)
    common["P9结果状态"] = p9.state
    common["P9主酶候选状态"] = p9.status
    common["P9推荐计量组合"] = p9.recommended_combination_id
    common["P9计量组合摘要"] = list(p9.combinations)
    raw_step = getattr(config, "step", None)
    if raw_step is not None:
        step_index = _as_int(raw_step, "step", minimum=1)
        if step_index > len(candidate_steps):
            raise ValueError(
                f"RetroPath 候选 {rank} 中没有步骤 {step_index}；"
                f"可用步骤：{list(range(1, len(candidate_steps) + 1))}"
            )
        step_detail = _step_display(candidate_steps[step_index - 1])
        step_detail["P8严格验证定向通量"] = overlay.flux_by_candidate_step.get(
            (rank, str(candidate_steps[step_index - 1].get("step_id") or ""))
        )
        step_detail["P9主酶候选"] = [
            _p9_candidate_display(row)
            for row in p9.candidates_by_step.get(step_index, tuple())
        ]
        return {
            **common,
            "步骤编号": step_index,
            "步骤详情": step_detail,
        }

    route_summary = _route_summary(route)
    displayed_steps = [_step_display(row) for row in candidate_steps]
    for row, displayed in zip(candidate_steps, displayed_steps):
        step_index = _as_int(row.get("step_index"), "step_index", minimum=1)
        displayed["P8严格验证定向通量"] = overlay.flux_by_candidate_step.get(
            (rank, str(row.get("step_id") or ""))
        )
        displayed["P9主酶候选数"] = len(
            p9.candidates_by_step.get(step_index, tuple())
        )
    if overlay.state != "current":
        validation_warning = (
            "该候选尚未完成有效的 P8 计量和严格 GEM 验证，不可进入正式设计"
        )
    elif p9.state != "current":
        validation_warning = (
            "该候选已有 P8 结果，但尚未完成有效的 P9 主酶候选检索"
        )
    else:
        validation_warning = (
            "P8/P9 仅提供计算可行性和候选酶证据，RP2 步骤仍需人工复核"
        )
    return {
        **common,
        "候选概览": route_summary,
        "KEGG前缀步骤编号": [
            step["步骤编号"]
            for row, step in zip(candidate_steps, displayed_steps)
            if str(row.get("step_source") or "").strip() == "kegg_expansion"
        ],
        "RetroPath预测步骤编号": [
            step["步骤编号"]
            for row, step in zip(candidate_steps, displayed_steps)
            if str(row.get("step_source") or "").strip() == "retropath"
        ],
        "反应DAG步骤": displayed_steps,
        "风险提示": _unique(
            [
                validation_warning,
                *(
                    ["该候选包含辅助片段，共底物或辅因子仍需恢复"]
                    if route_summary["包含辅助片段"]
                    else []
                ),
                *[
                    warning
                    for step in displayed_steps
                    for warning in step["风险提示"]
                ],
                overlay.warning,
                p9.warning,
                "RP2 主酶候选即使结构相似度为 1 也仍需人工复核，P9 不允许正式晋升",
            ]
        ),
    }


def run_retropath_info(config: Any) -> dict[str, Any]:
    result = get_retropath_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def run_retropath_candidate_info(config: Any) -> dict[str, Any]:
    result = get_retropath_candidate_info(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "get_retropath_candidate_info",
    "get_retropath_info",
    "run_retropath_candidate_info",
    "run_retropath_info",
]
