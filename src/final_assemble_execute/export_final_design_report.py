from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

from langchain.chat_models import init_chat_model
from langchain.tools import tool
from langchain_core.messages import HumanMessage, SystemMessage
from pydantic import BaseModel, Field

from src.runtime.monitor import monitor
from src.runtime.context import resolve_llm_settings
from src.tools.final_assemble_execute_tools.common import (
    as_dict,
    as_list,
    final_assembly_dir as resolve_final_assembly_dir,
    json_error,
    manifest_file,
    outputs_dir,
    read_json,
    relpath_or_abs,
    resolve_file_path,
    sanitize_name,
    session_root,
    sha256_file,
    write_json_atomic,
)

def _limit_items(items: list[Any], limit: int) -> list[Any]:
    return items[:limit] if len(items) > limit else items


def _extract_route_steps(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    steps = as_list(as_dict(manifest.get("solution")).get("steps"))
    result = []
    for step in _limit_items([item for item in steps if isinstance(item, dict)], 30):
        result.append({
            "step_index": step.get("step_index"),
            "reaction_id": step.get("reaction_id"),
            "reaction_name": step.get("reaction_name"),
            "produced_compound_id": step.get("produced_compound_id"),
            "produced_compound_name": step.get("produced_compound_name"),
            "enzyme_ecs": step.get("enzyme_ecs"),
            "ko_ids": step.get("ko_ids"),
            "oxygen_required": step.get("oxygen_required"),
            "direction": step.get("direction"),
            "status": step.get("status"),
        })
    return result


def _extract_selected_cds(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    selected = as_list(as_dict(manifest.get("cds_selection")).get("selected_cds"))
    result = []
    for item in selected:
        if not isinstance(item, dict):
            continue
        protein = as_dict(item.get("protein"))
        optimized = as_dict(item.get("optimized_cds"))
        result.append({
            "cds_id": item.get("cds_id"),
            "step_index": item.get("step_index"),
            "reaction_id": item.get("reaction_id"),
            "ec_number": item.get("ec_number"),
            "protein": {
                "accession": protein.get("accession"),
                "protein_name": protein.get("protein_name"),
                "organism_name": protein.get("organism_name"),
                "annotation_score": protein.get("annotation_score"),
            },
            "optimized_cds": {
                "length_nt": optimized.get("length_nt"),
                "gc_percent": optimized.get("gc_percent"),
                "start_codon": optimized.get("start_codon"),
                "stop_codon": optimized.get("stop_codon"),
                "length_multiple_of_3": optimized.get("length_multiple_of_3"),
                "sequence_file": as_dict(optimized.get("sequence_file")),
            },
        })
    return result


def _extract_parts_selection(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    selected = as_list(as_dict(manifest.get("parts_selection")).get("selected_parts"))
    result = []
    for item in selected:
        if not isinstance(item, dict):
            continue
        result.append({
            "cassette_index": item.get("cassette_index"),
            "role": item.get("role"),
            "part_id": item.get("part_id"),
            "part_name": item.get("part_name") or item.get("name"),
            "cds_id": item.get("cds_id"),
            "length_bp": as_dict(item.get("sequence_file")).get("length_bp"),
        })
    return result


def _extract_cassette_files(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    cassette_files = as_list(as_dict(manifest.get("assembled_expression_cassettes")).get("cassette_files"))
    result = []
    for item in cassette_files:
        if not isinstance(item, dict):
            continue
        components = []
        for component in as_list(item.get("components")):
            if not isinstance(component, dict):
                continue
            components.append({
                "type": component.get("type"),
                "label": component.get("label"),
                "source_id": component.get("source_id"),
                "cds_id": component.get("cds_id"),
                "protein_name": component.get("protein_name"),
                "start_bp": component.get("start_bp"),
                "end_bp": component.get("end_bp"),
            })
        result.append({
            "cassette_index": item.get("cassette_index"),
            "path": item.get("path"),
            "length_bp": item.get("length_bp"),
            "component_count": item.get("component_count"),
            "components": components,
        })
    return result


def _file_status(outputs: Path, file_info: dict[str, Any]) -> dict[str, Any]:
    path_value = file_info.get("path")
    if not path_value:
        return {"path": "", "exists": False, "sha256_matches": False}

    path = resolve_file_path(outputs, path_value)
    exists = path.exists()
    expected_sha = str(file_info.get("sha256") or "")
    actual_sha = sha256_file(path) if exists else ""
    return {
        "path": relpath_or_abs(outputs, path),
        "exists": exists,
        "sha256": actual_sha or expected_sha,
        "sha256_matches": bool(exists and (not expected_sha or actual_sha == expected_sha)),
        "format": file_info.get("format"),
        "length_bp": file_info.get("length_bp"),
    }


def _read_assembly_report(outputs: Path, final_assembly: dict[str, Any]) -> tuple[Path | None, dict[str, Any]]:
    report_info = as_dict(as_dict(final_assembly.get("construct_files")).get("report"))
    path_value = report_info.get("path")
    if not path_value:
        return None, {}
    report_path = resolve_file_path(outputs, path_value)
    if not report_path.exists():
        return report_path, {}
    return report_path, read_json(report_path)


def _clean_legacy_markdown_registration(
    *,
    manifest: dict[str, Any],
    outputs: Path,
    assembly_report_path: Path | None,
    assembly_report: dict[str, Any],
) -> list[str]:
    changes = []
    final_assembly = as_dict(manifest.get("final_assembly"))
    construct_files = as_dict(final_assembly.get("construct_files"))

    if "markdown_report" in construct_files:
        construct_files.pop("markdown_report", None)
        changes.append("removed legacy final_assembly.construct_files.markdown_report")

    if assembly_report_path and assembly_report_path.exists() and "markdown_report" in assembly_report:
        assembly_report.pop("markdown_report", None)
        assembly_report_path.write_text(
            json.dumps(assembly_report, ensure_ascii=False, indent=2, default=str) + "\n",
            encoding="utf-8",
        )
        report_info = as_dict(construct_files.get("report"))
        if report_info.get("path"):
            report_info["sha256"] = sha256_file(assembly_report_path)
            construct_files["report"] = report_info
        changes.append("removed legacy assembly_report.markdown_report")

    if changes:
        final_assembly["construct_files"] = construct_files
        manifest["final_assembly"] = final_assembly

    return changes


def _markdown_report_payload(
    *,
    manifest: dict[str, Any],
    assembly_report: dict[str, Any],
    outputs: Path,
    assembly_report_path: Path | None,
) -> dict[str, Any]:
    solution = as_dict(manifest.get("solution"))
    plasmid_selection = as_dict(manifest.get("plasmid_selection"))
    vector = as_dict(plasmid_selection.get("vector"))
    final_assembly = as_dict(manifest.get("final_assembly"))
    construct_files = as_dict(final_assembly.get("construct_files"))

    return {
        "manifest_revision": manifest.get("revision"),
        "solution_summary": as_dict(solution.get("summary")),
        "route_steps": _extract_route_steps(manifest),
        "protein_selection": {
            "chassis_key": as_dict(manifest.get("protein_selection")).get("chassis_key"),
            "recommended_design": as_dict(as_dict(manifest.get("protein_selection")).get("recommended_design")),
        },
        "selected_cds": _extract_selected_cds(manifest),
        "expression_cassette_assembly": as_dict(manifest.get("expression_cassette_assembly")),
        "parts_selection": _extract_parts_selection(manifest),
        "assembled_expression_cassettes": _extract_cassette_files(manifest),
        "plasmid_selection": {
            "source_kind": plasmid_selection.get("source_kind"),
            "vector": {
                "plasmid_id": vector.get("plasmid_id"),
                "addgene_id": vector.get("addgene_id"),
                "name": vector.get("name"),
                "description": vector.get("description"),
                "vector_type": vector.get("vector_type"),
                "bacterial_resistance": vector.get("bacterial_resistance"),
                "copy_number": vector.get("copy_number"),
                "length_bp": vector.get("length_bp"),
                "topology": vector.get("topology"),
                "features": as_dict(vector.get("features")),
            },
        },
        "final_assembly_plan": as_dict(manifest.get("final_assembly_plan")),
        "final_assembly": {
            "assembled_at": final_assembly.get("assembled_at"),
            "construct_name": final_assembly.get("construct_name"),
            "assembly_method": as_dict(final_assembly.get("strategy")).get("assembly_method"),
            "target": final_assembly.get("target"),
            "vector": final_assembly.get("vector"),
            "insert": final_assembly.get("insert"),
            "length_bp": final_assembly.get("length_bp"),
            "topology": final_assembly.get("topology"),
            "component_count": final_assembly.get("component_count"),
            "components": _limit_items(as_list(final_assembly.get("components")), 120),
            "validation": final_assembly.get("validation"),
            "warnings": final_assembly.get("warnings"),
        },
        "assembly_report": {
            "path": relpath_or_abs(outputs, assembly_report_path) if assembly_report_path else "",
            "assembled_at": assembly_report.get("assembled_at"),
            "construct_name": assembly_report.get("construct_name"),
            "assembly_method": assembly_report.get("assembly_method"),
            "final_construct": assembly_report.get("final_construct"),
            "validation": assembly_report.get("validation"),
            "warnings": assembly_report.get("warnings"),
        },
        "output_files": {
            key: _file_status(outputs, as_dict(value))
            for key, value in construct_files.items()
        },
    }


def _extract_message_text(message: Any) -> str:
    content = getattr(message, "content", message)
    if isinstance(content, str):
        return content.strip()
    if isinstance(content, list):
        parts = []
        for item in content:
            if isinstance(item, str):
                parts.append(item)
            elif isinstance(item, dict) and item.get("type") == "text":
                parts.append(str(item.get("text", "")))
        return "".join(parts).strip()
    return str(content).strip()


def _strip_markdown_fence(markdown: str) -> str:
    text = markdown.strip()
    if not text.startswith("```"):
        return text
    lines = text.splitlines()
    if len(lines) >= 2 and lines[0].startswith("```") and lines[-1].strip() == "```":
        return "\n".join(lines[1:-1]).strip()
    return text


def _deepseek_compat_kwargs() -> dict[str, Any]:
    settings = resolve_llm_settings()
    text = " ".join([
        settings.get("provider", ""),
        settings.get("base_url", ""),
        settings.get("model", ""),
    ]).lower()
    if "deepseek" not in text:
        return {}
    return {"extra_body": {"thinking": {"type": "disabled"}}}


def _model_provider() -> str | None:
    settings = resolve_llm_settings()
    provider = settings.get("provider")
    if provider:
        return provider
    text = " ".join([settings.get("base_url", ""), settings.get("model", "")]).lower()
    if "deepseek" in text:
        return "deepseek"
    return None


def _fallback_markdown_report(payload: dict[str, Any], warning: str) -> str:
    solution = as_dict(payload.get("solution_summary"))
    final = as_dict(payload.get("final_assembly"))
    vector = as_dict(as_dict(payload.get("plasmid_selection")).get("vector"))
    selected_cds = as_list(payload.get("selected_cds"))
    output_files = as_dict(payload.get("output_files"))
    warnings_list = as_list(final.get("warnings"))

    lines = [
        f"# {final.get('construct_name') or '最终设计报告'}",
        "",
        "> 大模型报告生成失败，本文件为系统基于设计清单整理的结构化兜底版本。",
        "",
        "## 设计摘要",
        f"- 目标产物：{solution.get('target_compound_name', '')} ({solution.get('target_compound_id', '')})",
        f"- 最终构建：{final.get('construct_name', '')}",
        f"- 组装方式：{final.get('assembly_method', '')}",
        f"- 最终长度：{final.get('length_bp', '')} bp",
        f"- 拓扑结构：{final.get('topology', '')}",
        "",
        "## 载体与插入片段",
        f"- 载体：{vector.get('name', '')}，长度 {vector.get('length_bp', '')} bp，拓扑 {vector.get('topology', '')}",
        f"- 抗性：{vector.get('bacterial_resistance', '')}",
        f"- 插入片段长度：{as_dict(final.get('insert')).get('length_bp', '')} bp",
        "",
        "## 蛋白与 CDS",
    ]
    for cds in selected_cds:
        protein = as_dict(cds.get("protein"))
        optimized = as_dict(cds.get("optimized_cds"))
        lines.append(
            f"- {cds.get('cds_id')}: {protein.get('protein_name')} ({protein.get('accession')}), "
            f"{optimized.get('length_nt')} nt, GC {optimized.get('gc_percent')}%"
        )

    lines.extend([
        "",
        "## 输出文件",
        f"- GenBank：`{as_dict(output_files.get('genbank')).get('path', '')}`",
        f"- FASTA：`{as_dict(output_files.get('fasta')).get('path', '')}`",
        f"- 组装 JSON 报告：`{as_dict(output_files.get('report')).get('path', '')}`",
        "",
        "## 验证与风险提示",
    ])
    if warnings_list:
        lines.extend(f"- {item}" for item in warnings_list)
    else:
        lines.append("- 当前 final_assembly 未记录额外 warning。")
    lines.extend([
        f"- 报告生成 warning：{warning}",
        "",
        "## 后续建议",
        "- 该文件为理论构建结果，后续仍需人工复核酶切位点、阅读框、调控元件方向和实验可行性。",
    ])
    return "\n".join(lines) + "\n"


def generate_chinese_markdown_report(payload: dict[str, Any]) -> tuple[str, str, str | None]:
    system_prompt = """
你是 GLADE 合成生物学助手的中文技术报告撰写器。
请严格基于用户提供的 design_manifest 摘要和 final assembly 事实，生成一份中文 Markdown 报告。
要求：
- 只输出 Markdown，不要输出解释性前后缀。
- 不要编造未提供的实验结果、转化结果、产量、测序结果或湿实验验证。
- 不要输出完整 DNA 序列。
- 使用清晰标题、短段落和表格。
- 必须包含：设计摘要、目标产物与底盘/通路、蛋白与 CDS、表达盒与 parts、质粒载体、最终组装方案、输出文件、验证与风险提示、后续实验建议。
- 报告面向合成生物学研发人员，中文表达专业但易读。
""".strip()
    user_prompt = (
        "请根据下面 JSON 内容生成中文版 Markdown 最终设计报告。\n\n"
        "```json\n"
        f"{json.dumps(payload, ensure_ascii=False, indent=2, default=str)}\n"
        "```"
    )

    try:
        settings = resolve_llm_settings()
        llm = init_chat_model(
            model=settings["model"],
            model_provider=_model_provider(),
            api_key=settings.get("api_key"),
            base_url=settings.get("base_url"),
        )
        response = llm.invoke(
            [SystemMessage(content=system_prompt), HumanMessage(content=user_prompt)],
            **_deepseek_compat_kwargs(),
        )
        markdown = _strip_markdown_fence(_extract_message_text(response))
        if not markdown:
            raise ValueError("LLM returned empty Markdown report")
        return markdown.strip() + "\n", "llm", None
    except Exception as exc:
        warning = f"LLM Markdown report generation failed; fallback report was written: {exc}"
        return _fallback_markdown_report(payload, warning), "fallback", warning


class ExportFinalDesignReportArgs(BaseModel):
    overwrite: bool = Field(default=True, description="是否覆盖同名中文 Markdown 最终设计报告")


@tool(args_schema=ExportFinalDesignReportArgs)
def export_final_design_report(
    overwrite: bool = True,
) -> str:
    """
    根据当前 design_manifest 导出最终设计报告。
    
    调用时机：final_assembly 已生成或用户要求导出完整设计报告。
    返回：ok、报告文件路径、章节摘要和 warnings。
    限制：不修改设计，不执行组装；报告内容来自 manifest 和已有输出。
    """

    tool_name = "export_final_design_report"
    monitor.report_start(tool_name, {"overwrite": overwrite})
    try:
        monitor.report_running(tool_name, "正在读取最终组装结果和设计清单...", progress=0.18)
        root = session_root(None, None, None)
        outputs = outputs_dir(root, None)
        manifest_path = manifest_file(root, outputs, None)
        manifest = read_json(manifest_path)

        final_assembly = as_dict(manifest.get("final_assembly"))
        if final_assembly.get("status") != "assembled":
            raise ValueError("design_manifest.json 中还没有已完成的 final_assembly，请先执行最终组装。")

        assembly_report_path, assembly_report = _read_assembly_report(outputs, final_assembly)
        legacy_changes = _clean_legacy_markdown_registration(
            manifest=manifest,
            outputs=outputs,
            assembly_report_path=assembly_report_path,
            assembly_report=assembly_report,
        )

        construct_name = sanitize_name(
            final_assembly.get("construct_name")
            or assembly_report.get("construct_name")
            or "final_construct"
        )
        report_dir = resolve_final_assembly_dir(outputs, None)
        report_dir.mkdir(parents=True, exist_ok=True)
        report_path = report_dir / f"{construct_name}_中文设计报告.md"
        if report_path.exists() and not overwrite:
            raise FileExistsError(f"final design report already exists: {report_path}")

        monitor.report_running(tool_name, "正在调用大模型生成中文版 Markdown 报告...", progress=0.55)
        payload = _markdown_report_payload(
            manifest=manifest,
            assembly_report=assembly_report,
            outputs=outputs,
            assembly_report_path=assembly_report_path,
        )
        markdown, generated_by, generation_warning = generate_chinese_markdown_report(payload)
        report_path.write_text(markdown, encoding="utf-8")

        monitor.report_running(tool_name, "正在登记报告文件到设计清单...", progress=0.86)
        current_revision = int(manifest.get("revision", 0))
        warnings_list = []
        if generation_warning:
            warnings_list.append(generation_warning)

        construct_files = as_dict(final_assembly.get("construct_files"))
        genbank_info = as_dict(construct_files.get("genbank"))
        report_info = as_dict(construct_files.get("report"))
        manifest["final_design_report"] = {
            "status": "exported",
            "source": "export_final_design_report",
            "exported_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "construct_name": construct_name,
            "language": "zh-CN",
            "report_file": {
                "path": relpath_or_abs(outputs, report_path),
                "sha256": sha256_file(report_path),
                "format": "markdown",
                "language": "zh-CN",
                "generated_by": generated_by,
            },
            "based_on": {
                "manifest_revision": current_revision,
                "final_assembly_status": final_assembly.get("status"),
                "final_genbank_sha256": genbank_info.get("sha256"),
                "assembly_report_path": report_info.get("path"),
                "assembly_report_sha256": report_info.get("sha256"),
            },
            "legacy_cleanup": legacy_changes,
            "warnings": warnings_list,
        }
        manifest["revision"] = current_revision + 1
        write_json_atomic(manifest_path, manifest)

        monitor.report_end(
            tool_name,
            {
                "construct_name": construct_name,
                "report_file": relpath_or_abs(outputs, report_path),
                "generated_by": generated_by,
            },
        )
        return json.dumps(
            {
                "ok": True,
                "manifest_path": relpath_or_abs(root, manifest_path),
                "revision": manifest["revision"],
                "construct_name": construct_name,
                "markdown_report": relpath_or_abs(outputs, report_path),
                "generated_by": generated_by,
                "warning_count": len(warnings_list),
                "warnings": warnings_list,
                "legacy_cleanup": legacy_changes,
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return json_error(type(exc).__name__, message=str(exc))
