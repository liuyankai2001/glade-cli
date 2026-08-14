from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any, Literal

from src.pathway_analyze.target_id import validate_target_compound_id
from src.runtime.monitor import monitor
from src.tools.common.manifest import MANIFEST_STAGE_ORDER
from src.tools.common.session_paths import (
    chassis_reachability_dir,
    design_manifest_file,
    outputs_dir as resolve_outputs_dir,
    producible_kegg_compounds_file,
    session_dir as resolve_session_dir,
)

try:
    from langchain.tools import tool
except ModuleNotFoundError:
    def tool(func):
        func.name = func.__name__
        func.invoke = lambda payload: func(**payload)
        return func


GAP_DIR_PATTERN = re.compile(r"^kegg_gap_(C\d{5})$")
GEM_VALIDATION_DIRNAME = "gem_validation"


def _count_csv_rows(path: Path) -> int | None:
    if not path.exists() or not path.is_file():
        return None
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return sum(1 for _ in csv.DictReader(handle))
    except Exception:
        return None


def _file_summary(path: Path, *, count_rows: bool = False) -> dict[str, Any]:
    exists = path.exists() and path.is_file()
    return {
        "path": str(path.resolve()),
        "exists": exists,
        "size_bytes": path.stat().st_size if exists else 0,
        "row_count": _count_csv_rows(path) if count_rows and exists else None,
    }


def _read_json_object(path: Path) -> dict[str, Any]:
    if not path.exists() or not path.is_file():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return data


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists() or not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except Exception:
        return []


def _normalize_target(target_compound_id: str | None) -> str | None:
    if target_compound_id is None:
        return None
    if not str(target_compound_id).strip():
        return None
    return validate_target_compound_id(target_compound_id)


def _solution_ids(solutions_path: Path) -> list[int]:
    ids: list[int] = []
    for row in _read_csv_rows(solutions_path):
        value = row.get("solution_id")
        try:
            ids.append(int(float(str(value).strip())))
        except (TypeError, ValueError):
            continue
    return sorted(set(ids))


def _brief_run_config(path: Path) -> dict[str, Any]:
    data = _read_json_object(path)
    if not data:
        return {}
    keys = (
        "target",
        "resolved_target_compound",
        "search_mode_used",
        "did_fallback",
        "screening_rule_version",
        "max_total_steps",
        "max_new_enzymes",
        "max_solutions",
        "module_filter_mode",
    )
    return {key: data.get(key) for key in keys if key in data}


def _gap_target_from_dir(gap_dir: Path) -> str | None:
    match = GAP_DIR_PATTERN.fullmatch(gap_dir.name)
    return match.group(1) if match else None


def _gap_output_summary(gap_dir: Path, detail: Literal["compact", "full"]) -> dict[str, Any]:
    solutions_path = gap_dir / "solutions.csv"
    steps_path = gap_dir / "all_solution_steps.csv"
    run_config_path = gap_dir / "run_config.json"
    validation_dir = gap_dir / GEM_VALIDATION_DIRNAME
    validation_summary_path = validation_dir / "gem_validation_summary.csv"
    validation_fluxes_path = validation_dir / "gem_validation_route_fluxes.csv"
    electron_summary_path = gap_dir / "solution_electron_summary.csv"
    electron_requirements_path = gap_dir / "route_electron_requirements.csv"

    solutions_summary = _file_summary(solutions_path, count_rows=True)
    steps_summary = _file_summary(steps_path, count_rows=True)
    validation_summary = _file_summary(validation_summary_path, count_rows=True)
    validation_fluxes = _file_summary(validation_fluxes_path, count_rows=True)
    electron_summary = _file_summary(electron_summary_path, count_rows=True)
    electron_requirements = _file_summary(electron_requirements_path, count_rows=True)
    has_validation = validation_summary["exists"] or validation_fluxes["exists"]

    result: dict[str, Any] = {
        "target_compound_id": _gap_target_from_dir(gap_dir),
        "gap_dir": str(gap_dir.resolve()),
        "exists": gap_dir.exists() and gap_dir.is_dir(),
        "solutions_csv": solutions_summary,
        "all_solution_steps_csv": steps_summary,
        "solution_count": solutions_summary["row_count"] or 0,
        "step_count": steps_summary["row_count"] or 0,
        "has_electron_inference": electron_summary["exists"] or electron_requirements["exists"],
        "electron_summary_csv": electron_summary,
        "route_electron_requirements_csv": electron_requirements,
        "has_gem_validation": has_validation,
    }

    if detail == "full":
        result.update(
            {
                "solution_ids": _solution_ids(solutions_path),
                "run_config_json": _file_summary(run_config_path),
                "run_config": _brief_run_config(run_config_path),
                "gem_validation": {
                    "dir": str(validation_dir.resolve()),
                    "summary_csv": validation_summary,
                    "route_fluxes_csv": validation_fluxes,
                },
            }
        )

    return result


def _gap_dirs(output_root: Path, target_compound_id: str | None) -> list[Path]:
    if target_compound_id:
        return [output_root / f"kegg_gap_{target_compound_id}"]
    if not output_root.exists() or not output_root.is_dir():
        return []
    return sorted(
        path
        for path in output_root.glob("kegg_gap_C*")
        if path.is_dir() and _gap_target_from_dir(path)
    )


def _selected_solution(manifest: dict[str, Any]) -> dict[str, Any] | None:
    solution = manifest.get("solution")
    if not isinstance(solution, dict):
        return None

    summary = solution.get("summary") if isinstance(solution.get("summary"), dict) else {}
    steps = solution.get("steps") if isinstance(solution.get("steps"), list) else []
    return {
        "solution_id": solution.get("solution_id"),
        "target_compound_id": summary.get("target_compound_id"),
        "target_compound_name": summary.get("target_compound_name"),
        "step_count": len(steps),
        "gap_dir": solution.get("gap_dir"),
        "has_electron_inference": isinstance(manifest.get("electron_inference"), dict),
    }


def _downstream_fields_present(manifest: dict[str, Any]) -> list[str]:
    downstream_fields = MANIFEST_STAGE_ORDER[1:]
    return [field for field in downstream_fields if field in manifest]


def _current_stage(
    *,
    has_chassis_reachability: bool,
    has_gap_candidates: bool,
    has_selected_solution: bool,
    has_gem_validation: bool,
) -> str:
    if has_gem_validation:
        return "gem_validation_available"
    if has_selected_solution:
        return "solution_selected"
    if has_gap_candidates:
        return "gap_candidates_ready"
    if has_chassis_reachability:
        return "chassis_reachability_ready"
    return "empty"


def _next_actions(stage: str, target_compound_id: str | None, has_gap_outputs: bool) -> list[str]:
    if target_compound_id and not has_gap_outputs:
        return [
            f"No gap output found for {target_compound_id}; run analyze_chassis_metabolites if needed, then kegg_gap_analyze.",
        ]
    if stage == "empty":
        return [
            "Run analyze_chassis_metabolites before starting KEGG gap analysis.",
            "Confirm the target is a standard KEGG Compound ID in Cxxxxx format.",
        ]
    if stage == "chassis_reachability_ready":
        return [
            "Run kegg_gap_analyze with the confirmed target KEGG Compound ID.",
        ]
    if stage == "gap_candidates_ready":
        return [
            "Use list_all_solutions to compare candidate routes.",
            "Use list_solution_steps for detailed reaction evidence.",
            "Only call write_solution_to_manifest after the user explicitly selects a Solution.",
        ]
    if stage == "solution_selected":
        return [
            "Use the selected solution from design_manifest.json; do not reselect a route unless the user asks.",
            "Run gem_validate_relaxed_or_strict_l1 if the user requests GEM validation.",
            "Hand off to enzyme_system_selection_agent for protein selection after route validation or confirmation.",
        ]
    return [
        "Read existing GEM validation files before deciding whether to rerun validation.",
        "If validation is acceptable, continue downstream without rerunning KEGG gap analysis.",
    ]


@tool
def get_pathway_analysis_context(
    target_compound_id: str | None = None,
    detail: Literal["compact", "full"] = "compact",
) -> dict[str, Any]:
    """
    读取当前通路分析阶段上下文。

    调用时机：继续当前 session、查看已有通路/gap/validation 结果、验证前或避免重复运行分析前。
    返回：manifest、底盘可达性、kegg_gap_Cxxxxx 输出、已选 solution 和下一步建议。
    限制：只读；只支持新版 outputs/kegg_gap_Cxxxxx 和 manifest 顶层 solution，不迁移旧字段。
    """

    tool_name = "get_pathway_analysis_context"
    monitor.report_start(tool_name, {"target_compound_id": target_compound_id, "detail": detail})
    try:
        if detail not in {"compact", "full"}:
            raise ValueError("detail must be compact or full")

        normalized_target = _normalize_target(target_compound_id)
        session_root = resolve_session_dir()
        output_root = resolve_outputs_dir()
        manifest_path = design_manifest_file()
        manifest = _read_json_object(manifest_path)

        reachability_dir = chassis_reachability_dir()
        producible_path = producible_kegg_compounds_file()
        reachability_summary_path = reachability_dir / "analyze_chassis_metabolites_summary.csv"
        producible_summary = _file_summary(producible_path, count_rows=True)
        reachability_summary = _file_summary(reachability_summary_path, count_rows=True)
        has_chassis_reachability = producible_summary["exists"]

        selected_solution = _selected_solution(manifest)
        gap_outputs = [
            _gap_output_summary(gap_dir, detail)
            for gap_dir in _gap_dirs(output_root, normalized_target)
        ]
        has_gap_candidates = any(
            item.get("solution_count", 0) > 0 or item.get("step_count", 0) > 0
            for item in gap_outputs
        )
        has_gem_validation = any(item.get("has_gem_validation") for item in gap_outputs)
        stage = _current_stage(
            has_chassis_reachability=has_chassis_reachability,
            has_gap_candidates=has_gap_candidates,
            has_selected_solution=selected_solution is not None,
            has_gem_validation=has_gem_validation,
        )

        solution = manifest.get("solution") if isinstance(manifest.get("solution"), dict) else {}
        solution_summary = solution.get("summary") if isinstance(solution.get("summary"), dict) else {}
        result: dict[str, Any] = {
            "ok": True,
            "detail": detail,
            "session": {
                "session_dir": str(session_root.resolve()),
                "outputs_dir": str(output_root.resolve()),
            },
            "manifest": {
                "path": str(manifest_path.resolve()),
                "exists": manifest_path.exists() and manifest_path.is_file(),
                "revision": manifest.get("revision"),
                "has_solution": selected_solution is not None,
                "downstream_fields_present": _downstream_fields_present(manifest),
            },
            "current_stage": stage,
            "chassis_reachability": {
                "dir": str(reachability_dir.resolve()),
                "producible_kegg_compounds_csv": producible_summary,
                "summary_csv": reachability_summary,
            },
            "selected_solution": selected_solution,
            "gap_outputs": gap_outputs,
            "next_actions": _next_actions(stage, normalized_target, bool(gap_outputs and has_gap_candidates)),
        }

        if detail == "full" and selected_solution is not None:
            result["selected_solution_summary"] = {
                key: solution_summary.get(key)
                for key in (
                    "target_compound_id",
                    "target_compound_name",
                    "total_steps",
                    "heterologous_steps",
                    "heterologous_reaction_ids",
                    "heterologous_enzyme_ecs",
                    "route_total_nadph_burden",
                    "route_total_sam_burden",
                    "route_total_coa_burden",
                    "oxygen_required_steps",
                    "thermo_disfavored_steps",
                    "max_electron_risk_level",
                    "max_electron_risk_score",
                    "electron_system_status",
                    "requires_downstream_electron_design",
                )
                if key in solution_summary
            }
            if isinstance(manifest.get("electron_inference"), dict):
                result["selected_electron_inference"] = manifest["electron_inference"]

        monitor.report_end(
            tool_name,
            {
                "current_stage": stage,
                "gap_output_count": len(gap_outputs),
                "has_solution": selected_solution is not None,
            },
        )
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise
