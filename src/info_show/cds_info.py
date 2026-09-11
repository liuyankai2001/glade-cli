"""Read-only, separate views of raw and optimized CDS results."""

from __future__ import annotations

import math
import unicodedata
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.write_manifest.store import read_design_manifest


_HEADERS = ("蛋白 ID", "CDS 长度（nt）", "GC（%）", "CAI", "禁止位点数", "修改密码子数量")
_STATUS_LABELS = {"complete": "全部完成", "partial": "部分完成", "failed": "全部失败"}


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} 必须是 JSON 对象")
    return value


def _number(value: Any, *, integer: bool = False) -> int | float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if not math.isfinite(value) or value < 0:
        return None
    if integer:
        return int(value) if value == int(value) else None
    return value


def _directory(root: Path, path_value: Any) -> str | None:
    if not isinstance(path_value, str) or not path_value.strip():
        return None
    path = (root / path_value).resolve()
    if not path.is_relative_to(root):
        raise ValueError("CDS 文件路径超出当前项目输出目录")
    return str(path.parent)


def _selected_site_count(metrics: Mapping[str, Any]) -> tuple[bool, int | None]:
    # Legacy forbidden_site_count includes built-in defaults, not a user selection.
    hits = metrics.get("user_forbidden_site_hits")
    if not isinstance(hits, Mapping) or not hits:
        return False, None
    counts = [_number(value, integer=True) for value in hits.values()]
    if any(value is None for value in counts):
        return True, None
    return True, sum(counts)


def get_cds_info(config: Any, *, accession: str | None = None) -> dict[str, Any]:
    """Read one stage only, without falling back to the other stage."""
    show_raw = bool(getattr(config, "raw", False))
    manifest = read_design_manifest(config.manifest_output_path)
    target = str(manifest.get("target_compound_id") or "").strip()
    if target and target != config.target_name:
        raise ValueError(f"manifest 目标化合物为 {target}，与当前输入目标 {config.target_name} 不一致")
    result: dict[str, Any] = {
        "运行成功": True,
        "状态": "当前无结果",
        "视图": "raw" if show_raw else "optimized",
        "总数": 0,
        "优化成功数": 0,
        "直接使用数": 0,
        "失败数": 0,
        "CDS列表": [],
        "CDS文件目录": [],
        "失败信息": [],
        "未优化数": 0,
        "跳过数": 0,
    }
    if "cds_selection" not in manifest:
        result["提示"] = (
            "当前无原始 CDS 结果，请先运行 protein-to-cds。" if show_raw else
            "当前无优化后 CDS，请先运行 protein-to-cds 生成原始序列，再运行 optimize。"
        )
        return result
    if not target:
        raise ValueError("manifest 缺少 target_compound_id")
    selection = _mapping(manifest["cds_selection"], "cds_selection")
    if selection.get("schema_version") != "protein_to_cds.selection.v2":
        raise ValueError("不支持的 cds_selection schema_version，请重新运行 protein-to-cds")
    status = selection.get("status")
    if status not in _STATUS_LABELS:
        raise ValueError("cds_selection.status 无效")
    proteins = selection.get("proteins")
    failures = selection.get("failures")
    if not isinstance(proteins, list) or not isinstance(failures, list):
        raise ValueError("cds_selection.proteins 和 failures 必须是列表")

    root = Path(config.project_output_path).expanduser().resolve()
    directories: set[str] = set()
    seen: set[str] = set()
    for raw in proteins + failures:
        item = _mapping(raw, "CDS 记录")
        record_id = item.get("accession")
        if not isinstance(record_id, str) or not record_id.strip():
            raise ValueError("CDS 记录缺少 accession")
        if record_id in seen:
            raise ValueError(f"CDS 记录包含重复蛋白 ID：{record_id}")
        seen.add(record_id)

    for item in sorted(proteins, key=lambda row: row["accession"]):
        if accession is not None and item["accession"].upper() != accession.upper():
            continue
        optimized = _mapping(item.get("optimized_cds"), "optimized_cds")
        skipped = optimized.get("optimization_skipped") is True
        if skipped:
            result["跳过数"] += 1
            continue
        raw = item.get("raw_cds")
        if show_raw:
            if not isinstance(raw, Mapping) or not raw.get("path"):
                result["跳过数"] += 1
                continue
            selected = raw
        else:
            directory = _directory(root, optimized.get("path"))
            mode = optimized.get("processing_mode")
            is_optimized = (
                mode in {None, "codon_transformer_and_repair", "dna_chisel_gc_only"}
                and directory == str(root / "protein_to_cds" / "optimized_cds")
            )
            if not is_optimized:
                result["未优化数"] += 1
                continue
            selected = optimized
        metrics = _mapping(optimized.get("metrics", {}), "optimized_cds.metrics")
        metric_key = "raw" if show_raw else "final"
        final = _mapping(metrics.get(metric_key, {}), f"optimized_cds.metrics.{metric_key}")
        changes = _mapping(metrics.get("changes", {}), "optimized_cds.metrics.changes")
        sites_configured, site_count = _selected_site_count(final)
        result["CDS列表"].append({
            "蛋白ID": item["accession"],
            "直接使用": False,
            "长度_nt": _number(selected.get("length_nt"), integer=True),
            "GC百分比": _number(final.get("gc_percent")),
            "CAI": _number(final.get("cai")),
            "位点已配置": sites_configured,
            "禁止位点数": site_count,
            "修改密码子数量": None if show_raw else _number(changes.get("codon_change_count"), integer=True),
        })
        result["优化成功数"] += 1
        directory = _directory(root, selected.get("path"))
        if directory:
            directories.add(directory)

    result["失败信息"] = [
        {
            "蛋白ID": item["accession"],
            "错误类型": str(item.get("error_type") or ""),
            "原因": str(item.get("message") or "未记录失败原因"),
        }
        for item in sorted(failures, key=lambda row: row["accession"])
        if show_raw and (accession is None or item["accession"].upper() == accession.upper())
    ]
    result["状态"] = _STATUS_LABELS[status] if show_raw else "已优化"
    result["总数"] = len(result["CDS列表"])
    result["失败数"] = len(result["失败信息"])
    result["CDS文件目录"] = sorted(directories)
    if not result["CDS列表"] and not result["失败信息"]:
        result["提示"] = (
            "当前无原始 CDS 结果，请先运行 protein-to-cds。" if show_raw else
            "当前无优化后 CDS，请先运行 optimize；可用 info --cds --raw 查看原始结果。"
        )
    return result


def _cell(value: Any) -> str:
    return " ".join(str(value).split())


def _width(value: str) -> int:
    return sum(2 if unicodedata.east_asian_width(char) in "WF" else 1 for char in value)


def _metric(value: int | float | None, decimals: int | None = None) -> str:
    if value is None:
        return "未评估"
    return str(value) if decimals is None else f"{value:.{decimals}f}"


def format_cds_info(result: Mapping[str, Any]) -> str:
    """Render final metrics in six columns, or raw metrics in five columns."""
    if "提示" in result:
        return str(result["提示"])
    show_raw = result.get("视图") == "raw"
    summary = (
        f"原始 CDS（CodonTransformer，未经修正）：共 {result['总数']} 条"
        if show_raw else f"优化后 CDS：共 {result['总数']} 条｜未优化 {result['未优化数']} 条"
    )
    if result["失败数"]:
        summary += f"｜生成失败 {result['失败数']} 条"
    if result["跳过数"]:
        summary += f"｜跳过无对应阶段结果 {result['跳过数']} 条"
    lines = [summary]
    headers = _HEADERS[:5] if show_raw else _HEADERS
    rows = [
        [
            _cell(item["蛋白ID"]),
            _metric(item["长度_nt"]),
            _metric(item["GC百分比"], 2),
            _metric(item["CAI"], 4),
            _metric(item["禁止位点数"]) if item["位点已配置"] else "未配置",
            _metric(item["修改密码子数量"]),
        ][:len(headers)]
        for item in result["CDS列表"]
    ]
    if rows:
        widths = [max(_width(row[i]) for row in [headers, *rows]) for i in range(len(headers))]

        def render(row: Any) -> str:
            return " | ".join(cell + " " * (width - _width(cell)) for cell, width in zip(row, widths)).rstrip()

        lines.extend(["", render(headers), "-+-".join("-" * width for width in widths)])
        lines.extend(render(row) for row in rows)
    if result["CDS文件目录"]:
        lines.append("")
        lines.extend(f"CDS 文件目录：{directory}" for directory in result["CDS文件目录"])
    if result["失败信息"]:
        lines.extend(["", "失败原因："])
        for failure in result["失败信息"]:
            error_type = f"{_cell(failure['错误类型'])}：" if failure["错误类型"] else ""
            lines.append(f"- {_cell(failure['蛋白ID'])}：{error_type}{_cell(failure['原因'])}")
    return "\n".join(lines)


def run_cds_info(config: Any) -> dict[str, Any]:
    result = get_cds_info(config)
    print(format_cds_info(result))
    return result


__all__ = ["get_cds_info", "format_cds_info", "run_cds_info"]
