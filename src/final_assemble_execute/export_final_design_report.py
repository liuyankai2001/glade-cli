"""Generate a Chinese portfolio report for final in-silico assemblies."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from typing import Any

from langchain_core.messages import HumanMessage, SystemMessage

from src.config.config import ROOT
from src.protein_selection.config import build_chat_model, load_model_settings


ReportGenerator = Callable[[dict[str, Any]], str]


def build_report_payload(
    *,
    manifest: Mapping[str, Any],
    final_assembly: Mapping[str, Any],
    run_summary: Mapping[str, Any],
) -> dict[str, Any]:
    solution = manifest.get("solution")
    solution = solution if isinstance(solution, Mapping) else {}
    plasmid = manifest.get("plasmid_selection")
    plasmid = plasmid if isinstance(plasmid, Mapping) else {}
    vector = plasmid.get("vector")
    vector = vector if isinstance(vector, Mapping) else {}
    parts = manifest.get("parts_selection")
    parts = parts if isinstance(parts, Mapping) else {}
    constructs = final_assembly.get("constructs")
    constructs = constructs if isinstance(constructs, list) else []
    failures = final_assembly.get("failures")
    failures = failures if isinstance(failures, list) else []
    return {
        "target_compound_id": manifest.get("target_compound_id"),
        "solution_summary": solution.get("summary"),
        "plasmid": {
            "plasmid_id": vector.get("plasmid_id"),
            "name": vector.get("name"),
            "length_bp": vector.get("length_bp"),
            "copy_number_class": vector.get("copy_number_class"),
            "marker": vector.get("marker"),
            "assembly_policy": vector.get("assembly_policy"),
        },
        "selected_parts_design_ids": parts.get("selected_design_ids"),
        "execution": {
            "status": final_assembly.get("status"),
            "planned_design_count": final_assembly.get("planned_design_count"),
            "succeeded_count": final_assembly.get("succeeded_count"),
            "failed_count": final_assembly.get("failed_count"),
            "method_counts": run_summary.get("method_counts"),
            "output_dir": final_assembly.get("output_dir"),
            "run_summary": final_assembly.get("run_summary_file"),
        },
        "constructs": [
            {
                "parts_design_id": item.get("parts_design_id"),
                "assembly_method": item.get("assembly_method"),
                "enzymes": item.get("enzyme_summary"),
                "target": item.get("target"),
                "length_bp": item.get("length_bp"),
                "sequence_sha256": item.get("sequence_sha256"),
                "files": item.get("files"),
                "validation": item.get("validation"),
                "warnings": item.get("warnings"),
            }
            for item in constructs
            if isinstance(item, Mapping)
        ],
        "failures": [dict(item) for item in failures if isinstance(item, Mapping)],
        "warnings": final_assembly.get("warnings"),
    }


def _strip_markdown_fence(value: str) -> str:
    text = value.strip()
    if not text.startswith("```"):
        return text
    lines = text.splitlines()
    if len(lines) >= 2 and lines[-1].strip() == "```":
        return "\n".join(lines[1:-1]).strip()
    return text


def _message_text(message: Any) -> str:
    content = getattr(message, "content", message)
    if isinstance(content, str):
        return content.strip()
    if isinstance(content, list):
        return "".join(
            str(item.get("text") or "")
            for item in content
            if isinstance(item, Mapping) and item.get("type") == "text"
        ).strip()
    return str(content).strip()


def _llm_report(payload: dict[str, Any]) -> str:
    settings = load_model_settings(ROOT / ".env")
    model = build_chat_model(settings, max_tokens=4096, timeout_seconds=90)
    system = """
你是 GLADE 合成生物学设计报告撰写器。请仅基于输入 JSON 生成中文 Markdown。
必须明确这是理论计算组装，不代表湿实验成功。不要编造产量、转化、测序或表达结果，
不要输出完整 DNA 序列。报告应包含：设计摘要、质粒骨架、12 个方案表格、失败方案、
输出文件、计算验证和实验复核提示。只输出 Markdown。
""".strip()
    user = (
        "请根据下面事实生成最终设计报告：\n\n```json\n"
        + json.dumps(payload, ensure_ascii=False, indent=2, default=str)
        + "\n```"
    )
    kwargs: dict[str, Any] = {}
    identity = " ".join(
        (settings.model, settings.base_url, settings.provider)
    ).lower()
    if "deepseek" in identity:
        kwargs["extra_body"] = {"thinking": {"type": "disabled"}}
    response = model.invoke(
        [SystemMessage(content=system), HumanMessage(content=user)],
        **kwargs,
    )
    markdown = _strip_markdown_fence(_message_text(response))
    if not markdown:
        raise ValueError("LLM returned an empty report")
    return markdown.rstrip() + "\n"


def _fallback_report(payload: Mapping[str, Any], warning: str) -> str:
    execution = payload.get("execution")
    execution = execution if isinstance(execution, Mapping) else {}
    plasmid = payload.get("plasmid")
    plasmid = plasmid if isinstance(plasmid, Mapping) else {}
    constructs = payload.get("constructs")
    constructs = constructs if isinstance(constructs, list) else []
    failures = payload.get("failures")
    failures = failures if isinstance(failures, list) else []
    lines = [
        "# GLADE 最终理论组装报告",
        "",
        "> 本报告描述计算生成的质粒设计，不代表湿实验组装、表达或产物合成已经成功。",
        "",
        "## 运行摘要",
        "",
        f"- 目标化合物：{payload.get('target_compound_id', '')}",
        f"- 状态：{execution.get('status', '')}",
        f"- 成功方案：{execution.get('succeeded_count', 0)}",
        f"- 失败方案：{execution.get('failed_count', 0)}",
        f"- 输出目录：`{execution.get('output_dir', '')}`",
        "",
        "## 质粒骨架",
        "",
        f"- 名称：{plasmid.get('name', '')}",
        f"- ID：{plasmid.get('plasmid_id', '')}",
        f"- 原始长度：{plasmid.get('length_bp', '')} bp",
        f"- 拷贝数类型：{plasmid.get('copy_number_class', '')}",
        f"- 筛选标记：{plasmid.get('marker', '')}",
        "",
        "## 成功生成的方案",
        "",
        "| Design | 方法 | 酶/线性化 | 最终长度 | GenBank |",
        "|---:|---|---|---:|---|",
    ]
    for item in constructs:
        if not isinstance(item, Mapping):
            continue
        files = item.get("files")
        files = files if isinstance(files, Mapping) else {}
        genbank = files.get("genbank")
        genbank = genbank if isinstance(genbank, Mapping) else {}
        lines.append(
            f"| {item.get('parts_design_id', '')} | "
            f"{item.get('assembly_method', '')} | "
            f"{item.get('enzymes', '')} | "
            f"{item.get('length_bp', '')} | "
            f"`{genbank.get('path', '')}` |"
        )
    lines.extend(["", "## 失败方案", ""])
    if failures:
        for item in failures:
            if isinstance(item, Mapping):
                lines.append(
                    f"- Design {item.get('parts_design_id', '')}: "
                    f"{item.get('error_type', '')}: {item.get('message', '')}"
                )
    else:
        lines.append("- 无。")
    lines.extend(
        [
            "",
            "## 计算验证与实验提示",
            "",
            "- GenBank 和 FASTA 已重新解析并与理论序列逐碱基比较。",
            "- 限制酶、缓冲液、引物、反应条件和宿主适配仍需实验人员复核。",
            f"- LLM 报告生成提示：{warning}",
            "",
        ]
    )
    return "\n".join(lines)


def generate_final_design_report(
    payload: dict[str, Any],
    *,
    reporter: ReportGenerator | None = None,
) -> tuple[str, str, str | None]:
    try:
        markdown = reporter(payload) if reporter is not None else _llm_report(payload)
        markdown = _strip_markdown_fence(str(markdown))
        if not markdown:
            raise ValueError("report generator returned empty Markdown")
        return markdown.rstrip() + "\n", ("injected" if reporter else "llm"), None
    except Exception as exc:
        warning = f"LLM report generation failed; deterministic fallback used: {exc}"
        return _fallback_report(payload, warning), "fallback", warning


__all__ = [
    "ReportGenerator",
    "build_report_payload",
    "generate_final_design_report",
]
