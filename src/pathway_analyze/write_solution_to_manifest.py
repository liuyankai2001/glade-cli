from __future__ import annotations

import csv
import json
import re
import tempfile
from pathlib import Path
from typing import Any

from src.pathway_analyze.gem_validation import validation_depth_output_dir
from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir


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

SOLUTION_STEP_FIELDS = (
    "step_index",
    "status",
    "reaction_id",
    "reaction_name",
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

ELECTRON_SUMMARY_FIELDS = (
    "max_electron_risk_level",
    "max_electron_risk_score",
    "electron_system_status",
    "requires_downstream_electron_design",
)

ELECTRON_RISK_STEP_FIELDS = (
    "step_index",
    "reaction_id",
    "electron_carrier_ids",
    "electron_requirement_classes",
    "electron_risk_level",
    "electron_risk_score",
    "electron_risk_evidence",
    "requires_downstream_electron_design",
)

KEGG_COMPOUND_ID_PATTERN = re.compile(r"^C\d{5}$")


def _validate_target_compound_id(value: Any) -> str:
    target_compound = str(value or "").strip()
    if not KEGG_COMPOUND_ID_PATTERN.fullmatch(target_compound):
        raise ValueError(
            f"target_name 必须是 Cxxxxx 格式的 KEGG Compound ID，当前值："
            f"{target_compound!r}"
        )
    return target_compound


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
        raise ValueError(f"CSV 中存在无效的 solution_id：{row.get('solution_id')!r}") from exc


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


def _select_optional_solution_row(path: Path, solution_id: int) -> dict[str, Any]:
    for row in _read_csv_rows(path, required=False):
        if _row_solution_id(row) == solution_id:
            return _normalize_row(row)
    return {}


def _select_optional_solution_rows(path: Path, solution_id: int) -> list[dict[str, Any]]:
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


def _select_passed_validation(path: Path, solution_id: int) -> dict[str, Any]:
    matches = [
        row
        for row in _read_csv_rows(path)
        if str(row.get("validation_mode", "")).strip() == "per-solution"
        and _parse_solution_ids(row.get("solution_ids")) == (solution_id,)
    ]
    if not matches:
        raise ValueError(
            f"solution {solution_id} 没有独立通量验证结果；请先执行 "
            f"validate -s {solution_id} -m per"
        )

    validation = _normalize_row(matches[-1])
    status = str(validation.get("validation_status") or "")
    if not status.startswith("PASS_"):
        raise ValueError(
            f"solution {solution_id} 未通过通量验证，当前状态：{status or '未知'}"
        )
    return validation


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"schema_version": "design_manifest.v1", "revision": 0}
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict):
        raise ValueError(f"manifest 根节点必须是 JSON 对象：{path}")
    manifest.setdefault("schema_version", "design_manifest.v1")
    manifest.setdefault("revision", 0)
    return manifest


def _write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            delete=False,
            suffix=".tmp",
        ) as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def select_solution(config: Any) -> dict[str, Any]:
    """选择一条已通过独立 GEM 验证的路径，并写入项目 manifest。"""

    target_compound = _validate_target_compound_id(config.target_name)
    expansion_depth = int(getattr(config, "depth", 0))
    if expansion_depth < 0:
        raise ValueError("depth must be greater than or equal to 0")
    configured_solution_id = getattr(config, "solution_id", None)
    if configured_solution_id is None:
        raise ValueError("未指定 solution_id，请将 select 命令的 --solution 写入 RunConfig")
    solution_id = int(configured_solution_id)
    if solution_id < 1:
        raise ValueError("solution_id 必须是正整数")

    gap_root_dir = Path(config.gap_output_path).expanduser().resolve()
    gap_dir = gap_depth_output_dir(gap_root_dir, expansion_depth)
    validation_dir = validation_depth_output_dir(gap_root_dir, expansion_depth)
    manifest_path = Path(config.manifest_output_path)

    summary = _select_solution_summary(gap_dir / "solutions.csv", solution_id)
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

    steps = _select_solution_steps(gap_dir / "all_solution_steps.csv", solution_id)
    validation = _select_passed_validation(
        validation_dir / "gem_validation_summary.csv",
        solution_id,
    )
    electron_summary = _select_optional_solution_row(
        gap_dir / "solution_electron_summary.csv",
        solution_id,
    )
    electron_steps = _select_optional_solution_rows(
        gap_dir / "route_electron_requirements.csv",
        solution_id,
    )

    solution_payload = {
        "gap_dir": str(gap_dir.resolve()),
        "expansion_depth": expansion_depth,
        "solution_id": solution_id,
        "summary": _keep_fields(summary, SOLUTION_SUMMARY_FIELDS),
        "steps": [_keep_fields(step, SOLUTION_STEP_FIELDS) for step in steps],
        "validation": _keep_fields(validation, VALIDATION_FIELDS),
    }
    electron_payload = {
        "summary": (
            _keep_fields(electron_summary, ELECTRON_SUMMARY_FIELDS)
            if electron_summary
            else {}
        ),
        "risk_steps": [
            _keep_fields(step, ELECTRON_RISK_STEP_FIELDS)
            for step in electron_steps
        ],
    }

    manifest = _read_manifest(manifest_path)
    try:
        current_revision = int(manifest["revision"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"manifest revision 无效：{manifest.get('revision')!r}") from exc

    manifest["target_compound_id"] = target_compound
    manifest["solution"] = solution_payload
    manifest["electron_inference"] = electron_payload
    manifest["revision"] = current_revision + 1
    _write_json_atomic(manifest_path, manifest)

    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "路径编号": solution_id,
        "扩展深度": expansion_depth,
        "步骤数量": len(steps),
        "通量验证状态": validation.get("validation_status"),
        "辅因子模式": validation.get("cofactor_mode"),
        "电子系统状态": electron_summary.get("electron_system_status"),
        "是否需要后续电子传递设计": electron_summary.get(
            "requires_downstream_electron_design"
        ),
        "清单文件": str(manifest_path.resolve()),
        "清单版本": manifest["revision"],
    }


def run_select_solution(config: Any) -> dict[str, Any]:
    """argparse 命令入口；参数解析和 RunConfig 构造由 main.py 负责。"""

    result = select_solution(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result
