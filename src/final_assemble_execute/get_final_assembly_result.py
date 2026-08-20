from __future__ import annotations

import json

from langchain.tools import tool
from pydantic import BaseModel

from src.runtime.monitor import monitor
from src.tools.final_assemble_execute_tools.common import (
    final_assembly_summary,
    json_error,
    manifest_file,
    outputs_dir,
    read_json,
    relpath_or_abs,
    session_root,
)


class GetFinalAssemblyResultArgs(BaseModel):
    pass


@tool(args_schema=GetFinalAssemblyResultArgs)
def get_final_assembly_result(
) -> str:
    """
    读取已生成的 final_assembly 结果摘要。
    
    调用时机：用户要求查看最终组装产物、导出文件路径或执行结果状态。
    返回：ok、final_assembly 摘要、GenBank/FASTA/report 路径和 warnings。
    限制：只读；不执行组装，不修改计划或 manifest。
    """

    tool_name = "get_final_assembly_result"
    monitor.report_start(tool_name)
    try:
        root = session_root(None, None, None)
        outputs = outputs_dir(root, None)
        manifest = manifest_file(root, outputs, None)
        data = read_json(manifest)
        summary = final_assembly_summary(data, outputs)
        monitor.report_end(tool_name, {"available": summary.get("available")})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": relpath_or_abs(root, manifest),
                "revision": data.get("revision"),
                "final_assembly": summary,
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return json_error(type(exc).__name__, message=str(exc))
