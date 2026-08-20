
from __future__ import annotations

import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import design_manifest_file


class WriteExpressionCassettePlanArgs(BaseModel):
    plan: list[list[str]] = Field(
        description="表达盒组装方式。第一维代表第几个表达盒，第二维为该表达盒包含的优化 CDS cds_id 列表"
    )


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"design_manifest.json root must be a JSON object: {path}")
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


def _selected_cds_ids(manifest: dict[str, Any]) -> list[str]:
    cds_selection = manifest.get("cds_selection")
    if not isinstance(cds_selection, dict):
        return []

    selected_cds = cds_selection.get("selected_cds")
    if selected_cds is None:
        return []
    if not isinstance(selected_cds, list):
        raise ValueError('design_manifest.json field "cds_selection.selected_cds" must be a list')

    cds_ids = []
    for item in selected_cds:
        if isinstance(item, dict):
            cds_id = str(item.get("cds_id") or "").strip()
            if cds_id:
                cds_ids.append(cds_id)
    return cds_ids


def _normalize_plan(plan: list[list[str]]) -> list[list[str]]:
    if not isinstance(plan, list) or not plan:
        raise ValueError("plan must be a non-empty two-dimensional list")

    normalized_plan = []
    for cassette_index, cassette in enumerate(plan, start=1):
        if not isinstance(cassette, list) or not cassette:
            raise ValueError(f"plan cassette {cassette_index} must be a non-empty list")

        cds_ids = []
        for cds_id in cassette:
            text = str(cds_id or "").strip()
            if not text:
                raise ValueError(f"plan cassette {cassette_index} contains empty cds_id")
            cds_ids.append(text)
        normalized_plan.append(cds_ids)

    return normalized_plan


def _validate_plan_cds_ids(plan: list[list[str]], selected_cds_ids: list[str]) -> None:
    selected_set = set(selected_cds_ids)
    assigned_ids = [cds_id for cassette in plan for cds_id in cassette]
    assigned_set = set(assigned_ids)

    duplicated_ids = sorted({
        cds_id
        for cds_id in assigned_ids
        if assigned_ids.count(cds_id) > 1
    })
    if duplicated_ids:
        raise ValueError(f"plan contains duplicated cds_id: {duplicated_ids}")

    unknown_ids = [cds_id for cds_id in assigned_ids if cds_id not in selected_set]
    if unknown_ids:
        raise ValueError(f"plan contains cds_id not found in selected_cds: {unknown_ids}")

    missing_ids = [cds_id for cds_id in selected_cds_ids if cds_id not in assigned_set]
    if missing_ids:
        raise ValueError(f"plan does not include all selected_cds cds_id: {missing_ids}")


def _cassette_items(plan: list[list[str]]) -> list[dict[str, Any]]:
    return [
        {
            "cassette_index": index,
            "cds_count": len(cds_ids),
            "cds_ids": cds_ids,
        }
        for index, cds_ids in enumerate(plan, start=1)
    ]


@tool(args_schema=WriteExpressionCassettePlanArgs)
def write_expression_cassette_plan(plan: list[list[str]]) -> str:
    """
    写入用户确认的表达盒分组计划。
    
    调用时机：已选 CDS 后，用户确认每个 cassette 包含哪些 cds_id。
    输入：plan，格式为 [[cds_id, ...], ...]。
    返回：ok、cassette_count、cassette_items、manifest_path 和 revision。
    限制：不选择 parts、不生成 GenBank；后续调用 recommend_expression_parts/select_parts。
    """

    tool_name = "write_expression_cassette_plan"
    monitor.report_start(tool_name, {"cassette_count": len(plan) if isinstance(plan, list) else None})
    try:
        manifest_path = design_manifest_file()
        manifest = _read_manifest(manifest_path)

        selected_cds_ids = _selected_cds_ids(manifest)
        if not selected_cds_ids:
            raise ValueError("design_manifest.json has no selected CDS in cds_selection.selected_cds")

        normalized_plan = _normalize_plan(plan)
        _validate_plan_cds_ids(normalized_plan, selected_cds_ids)

        current_revision = int(manifest.get("revision", 0))
        cassette_items = _cassette_items(normalized_plan)
        manifest["expression_cassette_assembly"] = {
            "status": "planned",
            "source": "expression_cassette_plan",
            "planned_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "cassette_count": len(cassette_items),
            "total_cds_count": sum(item["cds_count"] for item in cassette_items),
            "cassettes": cassette_items,
        }
        clear_downstream_fields(manifest, "expression_cassette_assembly")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_path, manifest)

        monitor.report_end(tool_name, {"cassette_count": len(cassette_items), "revision": manifest["revision"]})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": str(manifest_path),
                "revision": manifest["revision"],
                "cassette_count": len(cassette_items),
                "total_cds_count": sum(item["cds_count"] for item in cassette_items),
                "cassettes": cassette_items,
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(write_expression_cassette_plan.invoke({
        "plan": [["P21685_optimized"]],
    }))
