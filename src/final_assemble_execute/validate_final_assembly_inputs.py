from __future__ import annotations

import json

from langchain.tools import tool
from pydantic import BaseModel

from src.runtime.monitor import monitor
from src.tools.final_assemble_execute_tools.common import (
    json_error,
    prepare_context,
    relpath_or_abs,
)


class ValidateFinalAssemblyInputsArgs(BaseModel):
    pass


@tool(args_schema=ValidateFinalAssemblyInputsArgs)
def validate_final_assembly_inputs(
) -> str:
    """
    校验当前 final_assembly_plan 是否具备执行最终理论组装的输入。
    
    调用时机：execute_final_assembly 前或用户要求检查最终组装输入时。
    返回：ok、ready_for_execution、vector/insert/strategy 摘要和 warnings。
    限制：只读；不生成文件、不修改 manifest；失败时提示回到 final_assemble_plan_agent 修正。
    """

    tool_name = "validate_final_assembly_inputs"
    monitor.report_start(tool_name)
    try:
        context = prepare_context(None, None, None)
        vector = context["vector"]
        insert = context["insert"]
        target = context["target"]
        plan = context["plan"]

        result = {
            "ok": True,
            "ready_for_execution": True,
            "manifest_path": relpath_or_abs(context["session_root"], context["manifest_file"]),
            "construct_name": plan.get("construct_name", "final_construct"),
            "assembly_method": context["strategy"].get("assembly_method"),
            "vector": {
                "path": relpath_or_abs(context["outputs"], vector["path"]),
                "length_bp": len(vector["sequence"]),
                "topology": vector.get("topology", ""),
            },
            "insert": {
                "length_bp": len(insert["sequence"]),
                "cassette_count": len(insert["cassette_files"]),
                "cassette_files": insert["cassette_files"],
            },
            "target": {key: value for key, value in target.items() if key != "final_sequence"},
            "expected_final_length_bp": len(target["final_sequence"]),
            "warnings": context["warnings"],
            "next_action": "调用 execute_final_assembly 生成最终质粒 GenBank/FASTA 和报告，并写入 final_assembly。",
        }
        monitor.report_end(tool_name, {"ready_for_execution": True, "assembly_method": context["strategy"].get("assembly_method")})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return json_error(type(exc).__name__, message=str(exc), ready_for_execution=False)
