from __future__ import annotations

import json
from datetime import datetime
from langchain.tools import tool
from pydantic import BaseModel, Field

from src.runtime.monitor import monitor

from src.tools.common.manifest import clear_downstream_fields
from src.tools.final_assemble_execute_tools.common import (
    as_dict,
    file_summary,
    final_assembly_dir as resolve_final_assembly_dir,
    final_features,
    json_error,
    prepare_context,
    relpath_or_abs,
    sanitize_name,
    sha256_file,
    sha256_text,
    validate_written_outputs,
    write_fasta,
    write_genbank,
    write_json_atomic,
)


class ExecuteFinalAssemblyArgs(BaseModel):
    overwrite: bool = Field(default=True, description="是否覆盖同名 final construct 输出文件")


@tool(args_schema=ExecuteFinalAssemblyArgs)
def execute_final_assembly(
    overwrite: bool = True,
) -> str:
    """
    按已写入的 final_assembly_plan 生成理论最终质粒文件。
    
    调用时机：final_assembly_plan.status 为 planned，且用户要求执行最终组装。
    返回：ok、final_assembly 摘要、GenBank/FASTA/manifest/report 路径和 warnings。
    限制：不制定或修改组装计划；缺少计划或输入不完整时停止并提示回到 plan 阶段。
    """

    tool_name = "execute_final_assembly"
    monitor.report_start(tool_name, {"overwrite": overwrite})
    try:
        monitor.report_running(tool_name, "正在读取最终组装上下文...", progress=0.15)
        context = prepare_context(None, None, None)
        outputs = context["outputs"]
        manifest = context["manifest"]
        plan = context["plan"]
        vector = context["vector"]
        insert = context["insert"]
        target = context["target"]

        construct_name = sanitize_name(plan.get("construct_name", "final_construct"))
        assembly_dir = resolve_final_assembly_dir(outputs, None)
        assembly_dir.mkdir(parents=True, exist_ok=True)

        genbank_path = assembly_dir / f"{construct_name}.gb"
        fasta_path = assembly_dir / f"{construct_name}.fasta"
        report_path = assembly_dir / f"{construct_name}_assembly_report.json"
        if not overwrite:
            existing = [path for path in [genbank_path, fasta_path, report_path] if path.exists()]
            if existing:
                raise FileExistsError(f"final assembly output files already exist: {existing}")

        monitor.report_running(tool_name, "正在拼接最终 construct 序列并生成注释...", progress=0.45)
        final_sequence = target["final_sequence"]
        features, feature_warnings = final_features(vector, insert, target)
        warnings_list = context["warnings"] + feature_warnings
        topology = "circular" if str(vector.get("topology") or "").lower() == "circular" else "linear"

        monitor.report_running(tool_name, "正在写入 GenBank、FASTA 和报告文件...", progress=0.75)
        write_genbank(genbank_path, construct_name, final_sequence, topology, features)
        write_fasta(fasta_path, construct_name, final_sequence)
        output_validation = validate_written_outputs(genbank_path, fasta_path, final_sequence)
        warnings_list.extend(output_validation["genbank_parse_warnings"])
        warnings_list.extend(output_validation["fasta_parse_warnings"])

        report = {
            "ok": True,
            "source": "execute_final_assembly",
            "assembled_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "construct_name": construct_name,
            "assembly_method": context["strategy"].get("assembly_method"),
            "vector": {
                "path": relpath_or_abs(outputs, vector["path"]),
                "length_bp": len(vector["sequence"]),
                "topology": vector.get("topology", ""),
            },
            "insert": {
                "length_bp": len(insert["sequence"]),
                "cassette_files": insert["cassette_files"],
                "linker_sequence": insert["linker_sequence"],
            },
            "target": {key: value for key, value in target.items() if key != "final_sequence"},
            "final_construct": {
                "length_bp": len(final_sequence),
                "topology": topology,
                "genbank": relpath_or_abs(outputs, genbank_path),
                "fasta": relpath_or_abs(outputs, fasta_path),
                "sequence_sha256": sha256_text(final_sequence),
            },
            "components": features,
            "validation": output_validation,
            "warnings": warnings_list,
        }
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

        current_revision = int(manifest.get("revision", 0))
        manifest["final_assembly"] = {
            "status": "assembled",
            "source": "execute_final_assembly",
            "assembled_at": report["assembled_at"],
            "construct_name": construct_name,
            "strategy": context["strategy"],
            "target": report["target"],
            "vector": report["vector"],
            "insert": report["insert"],
            "construct_files": {
                "genbank": file_summary(outputs, genbank_path, sequence=final_sequence, file_format="genbank"),
                "fasta": file_summary(outputs, fasta_path, sequence=final_sequence, file_format="fasta"),
                "report": {
                    "path": relpath_or_abs(outputs, report_path),
                    "sha256": sha256_file(report_path),
                    "format": "json",
                },
            },
            "length_bp": len(final_sequence),
            "topology": topology,
            "component_count": len(features),
            "components": features,
            "validation": {
                "vector_genbank_validated": True,
                "insert_genbank_validated": True,
                "final_files_written": True,
                "sequence_sha256": sha256_text(final_sequence),
                **output_validation,
            },
            "warnings": warnings_list,
        }
        clear_downstream_fields(manifest, "final_assembly")
        manifest["revision"] = current_revision + 1
        write_json_atomic(context["manifest_file"], manifest)

        monitor.report_end(
            tool_name,
            {
                "construct_name": construct_name,
                "length_bp": len(final_sequence),
                "final_assembly_dir": str(assembly_dir),
            },
        )
        return json.dumps(
            {
                "ok": True,
                "manifest_path": relpath_or_abs(context["session_root"], context["manifest_file"]),
                "revision": manifest["revision"],
                "construct_name": construct_name,
                "assembly_method": context["strategy"].get("assembly_method"),
                "length_bp": len(final_sequence),
                "topology": topology,
                "construct_files": as_dict(manifest["final_assembly"].get("construct_files")),
                "warning_count": len(warnings_list),
                "warnings": warnings_list,
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return json_error(type(exc).__name__, message=str(exc))
