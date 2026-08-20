from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel

from src.runtime.monitor import monitor
from src.tools.common.session_paths import design_manifest_file, outputs_dir as resolve_outputs_dir, session_dir as resolve_session_dir


class GetPlasmidContextArgs(BaseModel):
    pass


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _session_root(session_dir: str | None, output_dir: str | None, manifest_path: str | None) -> Path:
    if session_dir:
        return Path(session_dir).resolve()
    if output_dir:
        output_path = Path(output_dir).resolve()
        if output_path.name == "outputs":
            return output_path.parent.resolve()
        return output_path.parent.resolve()
    if manifest_path:
        path = Path(manifest_path).resolve()
        if path.name == "design_manifest.json" and path.parent.name == "outputs":
            return path.parent.parent.resolve()
    return resolve_session_dir()


def _output_dir(session_root: Path, output_dir: str | None) -> Path:
    if output_dir:
        return Path(output_dir).resolve()
    return resolve_outputs_dir()


def _manifest_file(session_root: Path, output_dir: Path, manifest_path: str | None) -> Path:
    if manifest_path:
        path = Path(manifest_path)
        if not path.is_absolute():
            path = session_root / path
        return path.resolve()
    return design_manifest_file()


def _relpath_or_abs(base_dir: Path, path_value: Any) -> str:
    if not path_value:
        return ""
    path = Path(str(path_value))
    if not path.is_absolute():
        path = base_dir / path
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"design_manifest.json root must be a JSON object: {path}")
    return data


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _assembled_expression_cassettes(manifest: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    assembled = _as_dict(manifest.get("assembled_expression_cassettes"))
    cassette_files = [
        item
        for item in _as_list(assembled.get("cassette_files"))
        if isinstance(item, dict)
    ]
    cassette_files.sort(key=lambda item: int(item.get("cassette_index") or 0))

    files = []
    for item in cassette_files:
        components = [
            component
            for component in _as_list(item.get("components"))
            if isinstance(component, dict)
        ]
        files.append({
            "cassette_index": item.get("cassette_index"),
            "format": item.get("format"),
            "path": _relpath_or_abs(output_dir, item.get("path")),
            "length_bp": item.get("length_bp"),
            "component_count": item.get("component_count"),
            "components": [
                {
                    "type": component.get("type"),
                    "source_id": component.get("source_id"),
                    "cds_id": component.get("cds_id"),
                    "protein_name": component.get("protein_name"),
                    "start_bp": component.get("start_bp"),
                    "end_bp": component.get("end_bp"),
                }
                for component in components
            ],
        })

    return {
        "available": bool(cassette_files),
        "status": assembled.get("status"),
        "source": assembled.get("source"),
        "assembled_at": assembled.get("assembled_at"),
        "cassette_count": len(files),
        "total_length_bp": sum(
            int(item.get("length_bp") or 0)
            for item in files
        ),
        "cassette_files": files,
    }


def _cassette_plan(manifest: dict[str, Any]) -> dict[str, Any]:
    assembly = _as_dict(manifest.get("expression_cassette_assembly"))
    cassettes = [
        {
            "cassette_index": item.get("cassette_index"),
            "cds_count": item.get("cds_count"),
            "cds_ids": item.get("cds_ids", []),
        }
        for item in _as_list(assembly.get("cassettes"))
        if isinstance(item, dict)
    ]
    cassettes.sort(key=lambda item: int(item.get("cassette_index") or 0))
    return {
        "available": bool(cassettes),
        "status": assembly.get("status"),
        "cassette_count": assembly.get("cassette_count") or len(cassettes),
        "total_cds_count": assembly.get("total_cds_count"),
        "cassettes": cassettes,
    }


def _parts_summary(manifest: dict[str, Any]) -> dict[str, Any]:
    parts_selection = _as_dict(manifest.get("parts_selection"))
    selected_parts = [
        part
        for part in _as_list(parts_selection.get("selected_parts"))
        if isinstance(part, dict)
    ]
    by_role: dict[str, list[dict[str, Any]]] = {
        "promoter": [],
        "rbs": [],
        "terminator": [],
    }
    for part in selected_parts:
        role = str(part.get("role") or "").strip().lower()
        if role not in by_role:
            by_role[role] = []
        by_role[role].append({
            "cassette_index": part.get("cassette_index"),
            "cds_id": part.get("cds_id"),
            "part_id": part.get("part_id"),
            "length_bp": _as_dict(part.get("sequence_file")).get("length_bp"),
        })

    return {
        "available": bool(selected_parts),
        "status": parts_selection.get("status"),
        "selected_part_count": len(selected_parts),
        "promoters": by_role.get("promoter", []),
        "rbs": by_role.get("rbs", []),
        "terminators": by_role.get("terminator", []),
        "other_parts": {
            role: values
            for role, values in by_role.items()
            if role not in {"promoter", "rbs", "terminator"}
        },
    }


def _host_context(manifest: dict[str, Any]) -> dict[str, Any]:
    cds_selection = _as_dict(manifest.get("cds_selection"))
    host = _as_dict(cds_selection.get("host"))
    protein_selection = _as_dict(manifest.get("protein_selection"))
    return {
        "codon_optimization_host": {
            "name": host.get("name"),
            "codon_transformer_organism_id": host.get("codon_transformer_organism_id"),
        },
        "chassis_key": protein_selection.get("chassis_key"),
    }


def _plasmid_selection(manifest: dict[str, Any]) -> dict[str, Any]:
    plasmid_selection = manifest.get("plasmid_selection")
    if not isinstance(plasmid_selection, dict):
        return {"available": False}
    return {
        "available": True,
        "status": plasmid_selection.get("status"),
        "source": plasmid_selection.get("source"),
        "selected_at": plasmid_selection.get("selected_at"),
        "vector": plasmid_selection.get("vector", {}),
        "warnings": plasmid_selection.get("warnings", []),
    }


def _vector_requirements(assembled: dict[str, Any], parts: dict[str, Any]) -> dict[str, Any]:
    total_length = int(assembled.get("total_length_bp") or 0)
    cassette_count = int(assembled.get("cassette_count") or 0)
    warnings = []
    notes = []

    if assembled.get("available"):
        notes.append("当前表达盒已包含 promoter/RBS/CDS/terminator，优先选择可承载表达盒插入片段的 vector backbone。")
    if parts.get("available"):
        notes.append("载体通常不需要额外提供表达调控元件，除非用户明确希望使用载体自带 promoter。")
    if total_length >= 10000:
        warnings.append("表达盒总长度较大，建议优先考虑中/低拷贝或更稳定的载体骨架。")
    elif total_length >= 5000:
        notes.append("表达盒长度中等，建议比较高拷贝表达压力与中拷贝稳定性。")
    if cassette_count > 1:
        notes.append("存在多个表达盒，推荐关注载体稳定性、插入位点和装配策略。")

    return {
        "insert_total_length_bp": total_length,
        "cassette_count": cassette_count,
        "vector_should_provide": [
            "replication_origin",
            "selection_marker",
            "cloning_or_insertion_site",
            "host_compatible_backbone",
        ],
        "vector_usually_not_required_to_provide": [
            "promoter",
            "rbs",
            "terminator",
        ] if assembled.get("available") else [],
        "compatibility_notes": notes,
        "warnings": warnings,
    }


@tool(args_schema=GetPlasmidContextArgs)
def get_plasmid_context(
) -> str:
    """
    读取质粒选择阶段所需的表达盒、parts 和宿主上下文。
    
    调用时机：推荐或写入 vector/backbone 选择前。
    返回：ok、ready_for_plasmid_selection、missing_requirements、assembled 表达盒摘要、vector_requirements 和已有选择。
    限制：只读；不推荐模板、不写 plasmid_selection。
    """

    tool_name = "get_plasmid_context"
    monitor.report_start(tool_name)
    try:
        session_root = _session_root(None, None, None)
        outputs = _output_dir(session_root, None)
        manifest_file = _manifest_file(session_root, outputs, None)
        manifest = _read_json(manifest_file)

        assembled = _assembled_expression_cassettes(manifest, outputs)
        cassette_plan = _cassette_plan(manifest)
        parts = _parts_summary(manifest)
        plasmid_selection = _plasmid_selection(manifest)
        ready = (
            assembled.get("available") is True
            and assembled.get("status") == "assembled"
            and assembled.get("cassette_count", 0) > 0
        )

        if ready:
            message = "表达盒已组装完成，可以开始质粒/vector 选择。"
        else:
            message = "表达盒尚未组装完成，请先完成表达盒组装后再进行质粒选择。"

        result = {
            "ok": True,
            "ready_for_plasmid_selection": ready,
            "message": message,
            "manifest_path": _relpath_or_abs(session_root, manifest_file),
            "output_dir": _relpath_or_abs(session_root, outputs),
            "revision": manifest.get("revision"),
            "host_context": _host_context(manifest),
            "expression_cassette_assembly": cassette_plan,
            "parts_selection": parts,
            "assembled_expression_cassettes": assembled,
            "existing_plasmid_selection": plasmid_selection,
            "vector_requirements": _vector_requirements(assembled, parts),
            "next_actions": (
                [
                    "如果用户已有 vector 序列，调用 annotate_plasmid_sequence 注释复制起点、抗性标记等元件。",
                    "如果用户需要推荐质粒，调用 recommend_plasmid_templates 推荐 3-5 个质粒模板候选。",
                    "用户确认采用某个 vector 后，再调用 write_plasmid_selection 写入 plasmid_selection。",
                ]
                if ready
                else [
                    "先回到 expression_cassette_assembly_agent 完成 assembled_expression_cassettes。",
                ]
            ),
        }
        monitor.report_end(tool_name, {"ready_for_plasmid_selection": ready})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(get_plasmid_context.invoke({}))
