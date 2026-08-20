from __future__ import annotations

import json
import sys
from pathlib import Path

from langchain.tools import tool
from pydantic import BaseModel
from src.runtime.monitor import monitor

if __package__:
    from .cds_manifest_context import (
        build_cds_manifest_context,
        compact_validation_result,
        json_error,
    )
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from src.tools.expression_cassette_assembly_tools.cds_manifest_context import (
        build_cds_manifest_context,
        compact_validation_result,
        json_error,
    )


class CdsCountValidationArgs(BaseModel):
    pass


@tool(args_schema=CdsCountValidationArgs)
def cds_count_validation(
) -> str:
    """
    校验已选 CDS 数量是否满足表达盒规划要求。
    
    调用时机：write_expression_cassette_plan 前，或用户询问能否开始表达盒设计时。
    返回：ok、valid、selected_cds_count、recommended_accession_count、缺失/多余 accessions 和建议。
    限制：只读；不选择 CDS，不写表达盒计划。
    """

    tool_name = "cds_count_validation"
    monitor.report_start(tool_name)
    try:
        context = build_cds_manifest_context()
        result = compact_validation_result(context)
        monitor.report_end(tool_name, {"valid": result.get("valid")})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    session_dir = (
        Path(__file__).resolve().parents[3]
        / "agent_workspace"
        / "users"
        / "admin"
        / "sessions"
        / "0e09e4680bb24cb4a364f3d4c6c316a5"
    )
    manifest_path = session_dir / "outputs" / "design_manifest.json"
    print(cds_count_validation.invoke({
        "session_dir": str(session_dir),
        "manifest_path": str(manifest_path),
    }))
