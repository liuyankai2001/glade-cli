from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
from langchain.tools import tool

from src.pathway_analyze.target_id import validate_target_compound_id
from src.config.pathway_config import (
    FLUX_GATE_PANEL_C_FALLBACK_CONFIG,
    FLUX_GATE_PRIMARY_SEARCH_CONFIG,
)
from src.runtime.monitor import monitor
from src.tools.common.session_paths import kegg_gap_dir as resolve_kegg_gap_dir
from src.tools.common.session_paths import outputs_dir as resolve_outputs_dir
from src.tools.pathway_analyse_tools.gem_validation import (
    DEFAULT_FLUX_THRESHOLD,
    gem_validate_relaxed_or_strict_l1,
    safe_float,
)
from src.tools.pathway_analyse_tools.kegg_gap_analyze import kegg_gap_analyze


# 使用新字典保留模块内可替换性，避免调用方修改全局配置实例。
PRIMARY_SEARCH_PARAMS = FLUX_GATE_PRIMARY_SEARCH_CONFIG.as_dict()
PANEL_C_FALLBACK_SEARCH_PARAMS = FLUX_GATE_PANEL_C_FALLBACK_CONFIG.as_dict()

PASS_PREFIX = "PASS"


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_ready(item) for item in value]
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return _json_ready(value.item())
        except (TypeError, ValueError):
            pass
    return value


def _read_csv(path: str | Path) -> pd.DataFrame:
    csv_path = Path(path)
    if not csv_path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        # A zero-solution gap search currently materializes its aggregate CSVs
        # as UTF-8 BOM-only files.  That is a valid fail-closed result, not an
        # exceptional pipeline failure.
        return pd.DataFrame()


def _parse_first_solution_id(value: Any, route_scope: Any = "") -> int | None:
    tokens = []
    for part in str(value or "").replace(",", ";").split(";"):
        part = part.strip()
        if part:
            tokens.append(part)
    if not tokens and str(route_scope or "").startswith("solution_"):
        tokens.append(str(route_scope).replace("solution_", "", 1))
    if not tokens:
        return None
    try:
        return int(float(tokens[0]))
    except ValueError:
        return None


def _pass_rows(
    summary_df: pd.DataFrame,
    solutions_df: pd.DataFrame | None = None,
    electron_df: pd.DataFrame | None = None,
    flux_threshold: float = DEFAULT_FLUX_THRESHOLD,
) -> pd.DataFrame:
    if summary_df.empty or "validation_status" not in summary_df.columns:
        return pd.DataFrame()
    rows = summary_df.copy()
    rows["solution_id_for_rank"] = rows.apply(
        lambda row: _parse_first_solution_id(row.get("solution_ids", ""), row.get("route_scope", "")),
        axis=1,
    )
    rows["fba_product_flux_numeric"] = rows["fba_product_flux"].apply(
        lambda value: safe_float(value, 0.0)
    ) if "fba_product_flux" in rows.columns else 0.0
    solutions_df = solutions_df if solutions_df is not None else pd.DataFrame()
    electron_df = electron_df if electron_df is not None else pd.DataFrame()
    if not solutions_df.empty and "solution_id" in solutions_df.columns:
        solution_columns = [
            column
            for column in (
                "solution_id",
                "reaction_resolution_status",
                "blocking_reaction_count",
                "eligible_for_recommendation",
                "heterologous_steps",
                "total_steps",
            )
            if column in solutions_df.columns
        ]
        rows = rows.merge(
            solutions_df[solution_columns],
            how="left",
            left_on="solution_id_for_rank",
            right_on="solution_id",
        )
    if not electron_df.empty and {
        "solution_id", "max_electron_risk_score"
    }.issubset(electron_df.columns):
        rows = rows.merge(
            electron_df[["solution_id", "max_electron_risk_score"]],
            how="left",
            left_on="solution_id_for_rank",
            right_on="solution_id",
            suffixes=("", "_electron"),
        )
    blocking_count = (
        pd.to_numeric(rows["blocking_reaction_count"], errors="coerce").fillna(0)
        if "blocking_reaction_count" in rows.columns
        else pd.Series(0, index=rows.index)
    )
    eligible_values = rows.get("eligible_for_recommendation", True)
    if not isinstance(eligible_values, pd.Series):
        eligible_mask = pd.Series(bool(eligible_values), index=rows.index)
    else:
        eligible_mask = eligible_values.astype(str).str.lower().isin(
            {"true", "1", "yes"}
        )
    passed = rows[
        rows["validation_status"].astype(str).str.startswith(PASS_PREFIX)
        & (rows["fba_product_flux_numeric"] > flux_threshold)
        & rows["solution_id_for_rank"].notna()
        & (blocking_count == 0)
        & eligible_mask
    ].copy()
    if passed.empty:
        return passed
    for column, default in (
        ("heterologous_steps", 10**9),
        ("total_steps", 10**9),
        ("max_electron_risk_score", 10**9),
    ):
        if column not in passed.columns:
            passed[column] = default
        passed[column] = pd.to_numeric(passed[column], errors="coerce").fillna(default)
    return passed.sort_values(
        [
            "heterologous_steps",
            "total_steps",
            "max_electron_risk_score",
            "fba_product_flux_numeric",
            "solution_id_for_rank",
        ],
        ascending=[True, True, True, False, True],
    )


def _row_for_solution(df: pd.DataFrame, solution_id: int | None) -> dict[str, Any]:
    if solution_id is None or df.empty or "solution_id" not in df.columns:
        return {}
    rows = df[df["solution_id"].astype(str) == str(solution_id)]
    if rows.empty:
        return {}
    return rows.iloc[0].where(pd.notnull(rows.iloc[0]), None).to_dict()


def _stage_summary(
    *,
    stage_name: str,
    search_result: dict[str, Any],
    validation_result: dict[str, Any] | None,
    flux_threshold: float,
) -> dict[str, Any]:
    target_compound = str(search_result.get("target_compound") or "")
    output_dir = Path(str(search_result.get("output_dir") or resolve_kegg_gap_dir(target_compound, resolve_outputs_dir())))
    solutions_df = _read_csv(output_dir / "solutions.csv")
    electron_df = _read_csv(output_dir / "solution_electron_summary.csv")
    summary_csv = ""
    validation_dir = ""
    if validation_result:
        summary_csv = str(validation_result.get("validation_summary_csv") or "")
        validation_dir = str(validation_result.get("validation_dir") or "")
    summary_df = _read_csv(summary_csv) if summary_csv else pd.DataFrame()

    pass_df = _pass_rows(
        summary_df,
        solutions_df,
        electron_df,
        flux_threshold,
    )
    selected_solution_id: int | None = None
    selected_validation: dict[str, Any] = {}
    if not pass_df.empty:
        selected_solution_id = int(pass_df.iloc[0]["solution_id_for_rank"])
        selected_validation = pass_df.iloc[0].where(pd.notnull(pass_df.iloc[0]), None).to_dict()

    solution_count = int(search_result.get("solution_count") or len(solutions_df))
    recommendation_status = "flux_pass" if selected_solution_id is not None else (
        "route_found_no_flux" if solution_count > 0 else "not_found"
    )
    solution_row = _row_for_solution(solutions_df, selected_solution_id)
    electron_row = _row_for_solution(electron_df, selected_solution_id)

    return {
        "stage_name": stage_name,
        "recommendation_status": recommendation_status,
        "target_compound": target_compound,
        "output_dir": str(output_dir.resolve()),
        "solutions_file": str((output_dir / "solutions.csv").resolve()),
        "all_solution_steps_file": str((output_dir / "all_solution_steps.csv").resolve()),
        "validation_dir": validation_dir,
        "validation_summary_csv": summary_csv,
        "solution_count": solution_count,
        "flux_pass_solution_count": int(len(pass_df)),
        "recommended_solution_id": selected_solution_id,
        "validation_status": selected_validation.get("validation_status"),
        "fba_product_flux": selected_validation.get("fba_product_flux"),
        "pfba_product_flux": selected_validation.get("pfba_product_flux"),
        "fva_minimum": selected_validation.get("target_fva_minimum"),
        "fva_maximum": selected_validation.get("target_fva_maximum"),
        "electron_system_status": electron_row.get("electron_system_status"),
        "max_electron_risk_level": electron_row.get("max_electron_risk_level"),
        "max_electron_risk_score": electron_row.get("max_electron_risk_score"),
        "requires_downstream_electron_design": electron_row.get("requires_downstream_electron_design"),
        "total_steps": solution_row.get("total_steps"),
        "heterologous_steps": solution_row.get("heterologous_steps"),
        "search_mode_used": search_result.get("search_mode_used"),
        "electron_avoidance_mode": search_result.get("electron_avoidance_mode"),
        "electron_avoidance_fallback": search_result.get("electron_avoidance_fallback"),
        "reaction_resolution_mode": search_result.get("reaction_resolution_mode"),
        "rejected_reaction_route_count": int(
            search_result.get("rejected_reaction_route_count") or 0
        ),
    }


def _run_search_stage(
    *,
    target: str,
    stage_name: str,
    search_params: dict[str, Any],
    flux_threshold: float,
    validation_mode: str,
    fva_fraction: float,
    processes: int,
) -> dict[str, Any]:
    search_result = kegg_gap_analyze.invoke({"target": target, **search_params})
    validation_result: dict[str, Any] | None = None
    if int(search_result.get("solution_count") or 0) > 0:
        validation_result = gem_validate_relaxed_or_strict_l1.invoke(
            {
                "target": search_result.get("target_compound") or target,
                "mode": validation_mode,
                "cofactor_mode": "strict_l1",
                "flux_threshold": flux_threshold,
                "fva_fraction": fva_fraction,
                "processes": processes,
            }
        )
    return _stage_summary(
        stage_name=stage_name,
        search_result=search_result,
        validation_result=validation_result,
        flux_threshold=flux_threshold,
    )


def _write_audit(final_stage: dict[str, Any], audit_payload: dict[str, Any]) -> str:
    output_dir = Path(str(final_stage["output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = output_dir / "search_flux_gate_summary.json"
    audit_path.write_text(json.dumps(_json_ready(audit_payload), ensure_ascii=False, indent=2), encoding="utf-8")
    return str(audit_path.resolve())


@tool
def kegg_gap_search_with_flux_gate(
    target: str,
    flux_threshold: float = DEFAULT_FLUX_THRESHOLD,
    fva_fraction: float = 0.99,
    processes: int = 1,
) -> dict[str, Any]:
    """
    先用轻量 KEGG 搜索找路线，再用 strict_l1 GEM 通量门控筛选路线；若路线不存在或无通量，
    自动切换到 panel c 验证过的扩参 fallback 搜索。

    调用时机：新目标产物的默认通路搜索入口。
    限制：不写 design_manifest.json；用户确认 solution 后仍需调用 write_solution_to_manifest。
    """

    tool_name = "kegg_gap_search_with_flux_gate"
    monitor.report_start(tool_name, {"target": target})
    try:
        target = validate_target_compound_id(target)
        monitor.report_running(tool_name, "正在运行第一轮轻量 prefer 搜索和 strict_l1 通量验证...", progress=0.1)
        primary = _run_search_stage(
            target=target,
            stage_name="primary_prefer",
            search_params=PRIMARY_SEARCH_PARAMS,
            flux_threshold=flux_threshold,
            validation_mode="per-solution",
            fva_fraction=fva_fraction,
            processes=processes,
        )

        stages = [primary]
        final_stage = primary
        fallback_triggered = primary["recommendation_status"] != "flux_pass"

        if fallback_triggered:
            monitor.report_running(
                tool_name,
                "第一轮没有通过 strict_l1 flux 的路线，正在运行 panel c 扩参 fallback...",
                progress=0.55,
            )
            fallback = _run_search_stage(
                target=primary.get("target_compound") or target,
                stage_name="panel_c_expanded_fallback",
                search_params=PANEL_C_FALLBACK_SEARCH_PARAMS,
                flux_threshold=flux_threshold,
                validation_mode="per-solution",
                fva_fraction=fva_fraction,
                processes=processes,
            )
            stages.append(fallback)
            final_stage = fallback

        audit_payload = {
            "target": target,
            "target_compound": final_stage.get("target_compound"),
            "fallback_triggered": fallback_triggered,
            "final_stage_used": final_stage.get("stage_name"),
            "recommendation_status": final_stage.get("recommendation_status"),
            "recommended_solution_id": final_stage.get("recommended_solution_id"),
            "primary_search_params": PRIMARY_SEARCH_PARAMS,
            "panel_c_fallback_search_params": PANEL_C_FALLBACK_SEARCH_PARAMS,
            "stages": stages,
        }
        audit_path = _write_audit(final_stage, audit_payload)

        result = {
            "target_compound": final_stage.get("target_compound"),
            "recommended_solution_id": final_stage.get("recommended_solution_id"),
            "recommendation_status": final_stage.get("recommendation_status"),
            "stage_used": final_stage.get("stage_name"),
            "fallback_triggered": fallback_triggered,
            "validation_status": final_stage.get("validation_status"),
            "fba_product_flux": final_stage.get("fba_product_flux"),
            "pfba_product_flux": final_stage.get("pfba_product_flux"),
            "fva_minimum": final_stage.get("fva_minimum"),
            "fva_maximum": final_stage.get("fva_maximum"),
            "electron_system_status": final_stage.get("electron_system_status"),
            "max_electron_risk_level": final_stage.get("max_electron_risk_level"),
            "max_electron_risk_score": final_stage.get("max_electron_risk_score"),
            "requires_downstream_electron_design": final_stage.get("requires_downstream_electron_design"),
            "solution_count": final_stage.get("solution_count"),
            "flux_pass_solution_count": final_stage.get("flux_pass_solution_count"),
            "reaction_resolution_mode": final_stage.get("reaction_resolution_mode"),
            "rejected_reaction_route_count": final_stage.get(
                "rejected_reaction_route_count"
            ),
            "gap_dir": final_stage.get("output_dir"),
            "solutions_file": final_stage.get("solutions_file"),
            "all_solution_steps_file": final_stage.get("all_solution_steps_file"),
            "rejected_reaction_routes_file": str(
                (Path(str(final_stage.get("output_dir"))) / "rejected_reaction_routes.csv").resolve()
            ),
            "validation_dir": final_stage.get("validation_dir"),
            "validation_summary_csv": final_stage.get("validation_summary_csv"),
            "search_flux_gate_summary_json": audit_path,
        }
        result = _json_ready(result)
        monitor.report_end(tool_name, result)
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise
