from __future__ import annotations

import json
import sys
from pathlib import Path

from langchain.tools import tool
from pydantic import BaseModel
from src.runtime.monitor import monitor

if __package__:
    from .cds_manifest_context import build_cds_manifest_context, json_error
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from src.tools.expression_cassette_assembly_tools.cds_manifest_context import (
        build_cds_manifest_context,
        json_error,
    )


class GetCdsListArgs(BaseModel):
    pass


@tool(args_schema=GetCdsListArgs)
def get_cds_list(
) -> str:
    """
    读取当前 cds_selection，供表达盒规划使用。
    
    调用时机：用户要查看已选 CDS、规划 expression cassette 或校验 CDS 数量。
    返回：ok、valid、selected_cds、推荐 accessions、缺失项和 cassette 规划提示。
    限制：只读；不写 manifest，不选择 CDS，不组装表达盒。
    """

    tool_name = "get_cds_list"
    monitor.report_start(tool_name)
    try:
        result = build_cds_manifest_context()
        result["next_actions"] = [
            "根据 valid 判断 CDS 信息是否已经全部获取。",
            "根据 cds_count 和 assembly_plan_hints 选择每个表达盒包含几个 CDS。",
            "确认表达方案后，再选择 promoter、RBS、terminator 等表达元件。",
        ]
        monitor.report_end(tool_name, {"cds_count": result.get("cds_count"), "valid": result.get("valid")})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(get_cds_list.invoke({}))
