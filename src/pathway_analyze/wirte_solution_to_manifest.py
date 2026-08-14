from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path
from typing import Any
from langchain.tools import tool
from src.pathway_analyze.target_id import validate_target_compound_id
from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import design_manifest_file, kegg_gap_dir as resolve_kegg_gap_dir


SOLUTION_SUMMARY_FIELDS = (
    "target_compound_id",
    "target_compound_name",
    "reaction_resolution_status",
    "normalization_event_count",
    "normalization_events",
    "blocking_reaction_count",
    "blocking_reaction_ids",
    "eligible_for_recommendation",
)

SOLUTION_STEP_FIELDS = (
    "step_index",
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
    "source_reaction_ids",
    "resolution_action",
    "resolution_evidence",
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


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing CSV file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_optional_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "schema_version": "design_manifest.v1",
            "revision": 0,
        }
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Manifest root must be a JSON object: {path}")
    data.setdefault("schema_version", "design_manifest.v1")
    data.setdefault("revision", 0)
    return data


def _write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
        suffix=".tmp",
    ) as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        tmp_path = Path(handle.name)
    tmp_path.replace(path)


def _to_int(value: Any) -> int:
    return int(float(str(value).strip()))


def _optional_number(value: str) -> int | float | str | bool | None:
    text = str(value or "").strip()
    if text == "":
        return None
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        number = float(text)
    except ValueError:
        return text
    if number.is_integer():
        return int(number)
    return number


def _normalize_row(row: dict[str, str]) -> dict[str, Any]:
    return {key: _optional_number(value) for key, value in row.items()}


def _keep_fields(row: dict[str, Any], fields: tuple[str, ...]) -> dict[str, Any]:
    """只保留阶段交接所需字段，详细证据继续留在源 CSV。"""

    return {key: row.get(key) for key in fields}


def _select_solution_summary(solutions_path: Path, solution_id: int) -> dict[str, Any]:
    rows = _read_csv_rows(solutions_path)
    for row in rows:
        if _to_int(row.get("solution_id", "")) == solution_id:
            return _normalize_row(row)
    available = [row.get("solution_id", "") for row in rows]
    raise ValueError(f"solution_id {solution_id} not found in {solutions_path}. Available: {available}")


def _select_solution_steps(steps_path: Path, solution_id: int) -> list[dict[str, Any]]:
    rows = [
        row
        for row in _read_csv_rows(steps_path)
        if _to_int(row.get("solution_id", "")) == solution_id
    ]
    if not rows:
        raise ValueError(f"solution_id {solution_id} has no steps in {steps_path}")
    rows.sort(key=lambda row: _to_int(row.get("step_index", "0")))
    return [_normalize_row(row) for row in rows]


def _select_optional_solution_row(path: Path, solution_id: int) -> dict[str, Any]:
    for row in _read_optional_csv_rows(path):
        if _to_int(row.get("solution_id", "")) == solution_id:
            return _normalize_row(row)
    return {}


def _select_optional_solution_rows(path: Path, solution_id: int) -> list[dict[str, Any]]:
    rows = [
        row
        for row in _read_optional_csv_rows(path)
        if _to_int(row.get("solution_id", "")) == solution_id
    ]
    rows.sort(key=lambda row: _to_int(row.get("step_index", "0")))
    return [_normalize_row(row) for row in rows]


@tool
def write_solution_to_manifest(
    target_compound_id: str,
    solution_id: int,
    expected_revision: int | None = None,
) -> dict[str, Any]:
    """
    把用户确认的通路候选方案写入 design_manifest.json。
    
    调用时机：kegg_gap_analyze/list_solution_steps 后，用户确认某个 solution。
    输入：target_compound_id、solution_id 和可选的 expected_revision。
    返回：ok、写入的 solution 摘要、manifest_path 和 revision。
    限制：只写 solution 和 electron_inference；不选择蛋白、不做 CDS 或表达盒设计。
    """

    tool_name = "write_solution_to_manifest"
    monitor.report_start(tool_name, {"target_compound_id": target_compound_id, "solution_id": solution_id})
    try:
        target_compound_id = validate_target_compound_id(target_compound_id)
        manifest_file = design_manifest_file()
        gap_dir_path = resolve_kegg_gap_dir(target_compound_id)
        solutions_path = gap_dir_path / "solutions.csv"
        steps_path = gap_dir_path / "all_solution_steps.csv"
        electron_summary_path = gap_dir_path / "solution_electron_summary.csv"
        electron_requirements_path = gap_dir_path / "route_electron_requirements.csv"

        selected_solution_id = int(solution_id)
        summary = _select_solution_summary(solutions_path, selected_solution_id)
        if int(summary.get("blocking_reaction_count") or 0) > 0 or (
            summary.get("eligible_for_recommendation") is False
        ):
            raise ValueError(
                f"solution_id {selected_solution_id} is blocked by reaction resolution: "
                f"{summary.get('blocking_reaction_ids') or 'unknown reaction'}"
            )
        steps = _select_solution_steps(steps_path, selected_solution_id)
        electron_summary = _select_optional_solution_row(electron_summary_path, selected_solution_id)
        electron_steps = _select_optional_solution_rows(electron_requirements_path, selected_solution_id)

        solution_payload = {
            "gap_dir": str(gap_dir_path.resolve()),
            "solution_id": selected_solution_id,
            "summary": _keep_fields(summary, SOLUTION_SUMMARY_FIELDS),
            "steps": [
                _keep_fields(step, SOLUTION_STEP_FIELDS)
                for step in steps
            ],
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

        manifest = _read_manifest(manifest_file)
        current_revision = int(manifest.get("revision", 0))
        if expected_revision is not None and int(expected_revision) != current_revision:
            raise ValueError(
                f"manifest revision mismatch: expected {expected_revision}, current {current_revision}"
        )

        manifest["solution"] = solution_payload
        clear_downstream_fields(manifest, "solution")
        manifest["electron_inference"] = electron_payload
        clear_downstream_fields(manifest, "electron_inference")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_file, manifest)

        result = {
            "manifest_json": str(manifest_file.resolve()),
            "solution_id": selected_solution_id,
            "step_count": len(steps),
            "target_compound_id": summary.get("target_compound_id"),
            "target_compound_name": summary.get("target_compound_name"),
            "electron_inference_available": bool(electron_summary),
            "max_electron_risk_level": electron_summary.get("max_electron_risk_level"),
            "electron_system_status": electron_summary.get("electron_system_status"),
            "requires_downstream_electron_design": electron_summary.get("requires_downstream_electron_design"),
            "revision": manifest["revision"],
        }
        monitor.report_end(tool_name, result)
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise

if __name__ == "__main__":

    target = 'C05432'
    print(write_solution_to_manifest.invoke({"target_compound_id": target, "solution_id": 1}))
