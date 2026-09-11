"""Read-only, compact view of the final protein-to-CDS results."""

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


def get_cds_info(config: Any) -> dict[str, Any]:
    """Read recorded final metrics without rerunning optimization or auditing files."""
    manifest = read_design_manifest(config.manifest_output_path)
    target = str(manifest.get("target_compound_id") or "").strip()
    if target and target != config.target_name:
        raise ValueError(f"manifest 目标化合物为 {target}，与当前输入目标 {config.target_name} 不一致")
    result: dict[str, Any] = {
        "运行成功": True,
        "状态": "当前无结果",
        "总数": 0,
        "优化成功数": 0,
        "直接使用数": 0,
        "失败数": 0,
        "CDS列表": [],
        "CDS文件目录": [],
        "失败信息": [],
    }
    if "cds_selection" not in manifest:
        result["提示"] = "当前无 CDS 结果，请先运行 protein-to-cds。"
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
        accession = item.get("accession")
        if not isinstance(accession, str) or not accession.strip():
            raise ValueError("CDS 记录缺少 accession")
        if accession in seen:
            raise ValueError(f"CDS 记录包含重复蛋白 ID：{accession}")
        seen.add(accession)

    for item in sorted(proteins, key=lambda row: row["accession"]):
        optimized = _mapping(item.get("optimized_cds"), "optimized_cds")
        skipped = optimized.get("optimization_skipped") is True
        metrics = _mapping(optimized.get("metrics", {}), "optimized_cds.metrics")
        final = _mapping(metrics.get("final", {}), "optimized_cds.metrics.final")
        changes = _mapping(metrics.get("changes", {}), "optimized_cds.metrics.changes")
        result["CDS列表"].append({
            "蛋白ID": item["accession"],
            "直接使用": skipped,
            "长度_nt": _number(optimized.get("length_nt"), integer=True),
            "GC百分比": _number(final.get("gc_percent")),
            "CAI": _number(final.get("cai")),
            "禁止位点数": _number(final.get("forbidden_site_count"), integer=True),
            "修改密码子数量": _number(changes.get("codon_change_count"), integer=True),
        })
        result["直接使用数" if skipped else "优化成功数"] += 1
        directory = _directory(root, optimized.get("path"))
        if directory:
            directories.add(directory)

    result["失败信息"] = [
        {
            "蛋白ID": item["accession"],
            "错误类型": str(item.get("error_type") or ""),
            "原因": str(item.get("message") or "未记录失败原因"),
        }
        for item in sorted(failures, key=lambda row: row["accession"])
    ]
    result["状态"] = _STATUS_LABELS[status]
    result["总数"] = len(proteins) + len(failures)
    result["失败数"] = len(failures)
    result["CDS文件目录"] = sorted(directories)
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
    """Render the agreed six columns with Chinese terminal-width alignment."""
    if "提示" in result:
        return str(result["提示"])
    lines = [
        f"CDS：{result['状态']}｜共 {result['总数']} 条"
        f"｜优化成功 {result['优化成功数']} 条"
        f"｜直接使用 {result['直接使用数']} 条｜失败 {result['失败数']} 条"
    ]
    rows = [
        [
            _cell(item["蛋白ID"]) + ("（直接使用）" if item["直接使用"] else ""),
            _metric(item["长度_nt"]),
            _metric(item["GC百分比"], 2),
            _metric(item["CAI"], 4),
            _metric(item["禁止位点数"]),
            _metric(item["修改密码子数量"]),
        ]
        for item in result["CDS列表"]
    ]
    if rows:
        widths = [max(_width(row[i]) for row in [_HEADERS, *rows]) for i in range(6)]

        def render(row: Any) -> str:
            return " | ".join(cell + " " * (width - _width(cell)) for cell, width in zip(row, widths)).rstrip()

        lines.extend(["", render(_HEADERS), "-+-".join("-" * width for width in widths)])
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
