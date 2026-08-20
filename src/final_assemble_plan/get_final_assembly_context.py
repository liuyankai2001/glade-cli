from __future__ import annotations
import json
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel

from src.runtime.monitor import monitor
from src.tools.common.session_paths import (
    design_manifest_file,
    outputs_dir as resolve_outputs_dir,
    session_dir as resolve_session_dir,
)


class GetFinalAssemblyContextArgs(BaseModel):
    pass


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _session_root(
    session_dir: str | None,
    output_dir: str | None,
    manifest_path: str | None,
) -> Path:
    if session_dir:
        return Path(session_dir).resolve()

    if output_dir:
        output_path = Path(output_dir).resolve()
        return output_path.parent.resolve() if output_path.name == "outputs" else output_path.parent.resolve()

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
            session_relative = session_root / path
            output_relative = output_dir / path
            path = output_relative if output_relative.exists() else session_relative
        return path.resolve()

    return design_manifest_file()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"design_manifest.json root must be a JSON object: {path}")
    return data


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


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


def _resolve_file_path(output_dir: Path, path_value: Any) -> Path | None:
    if not path_value:
        return None

    path = Path(str(path_value))
    if path.is_absolute():
        return path.resolve()

    output_relative = (output_dir / path).resolve()
    if output_relative.exists():
        return output_relative

    session_relative = (output_dir.parent / path).resolve()
    if session_relative.exists():
        return session_relative

    return output_relative


def _file_summary(output_dir: Path, file_info: dict[str, Any]) -> dict[str, Any]:
    path_value = file_info.get("path")
    resolved_path = _resolve_file_path(output_dir, path_value)
    exists = bool(resolved_path and resolved_path.exists())

    return {
        "path": _relpath_or_abs(output_dir, resolved_path or path_value),
        "exists": exists,
        "sha256": file_info.get("sha256", ""),
        "sequence_sha256": file_info.get("sequence_sha256", ""),
        "length_bp": file_info.get("length_bp"),
        "format": file_info.get("format", ""),
    }


def _component_summary(component: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    sequence_file = _as_dict(component.get("sequence_file"))
    return {
        "type": component.get("type"),
        "source_id": component.get("source_id"),
        "cds_id": component.get("cds_id"),
        "protein_name": component.get("protein_name"),
        "start_bp": component.get("start_bp"),
        "end_bp": component.get("end_bp"),
        "sequence_file": _file_summary(output_dir, sequence_file) if sequence_file else {},
    }


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
        file_summary = _file_summary(output_dir, item)
        files.append({
            "cassette_index": item.get("cassette_index"),
            **file_summary,
            "component_count": item.get("component_count") or len(components),
            "components": [
                _component_summary(component, output_dir)
                for component in components
            ],
        })

    return {
        "available": bool(cassette_files),
        "status": assembled.get("status"),
        "source": assembled.get("source"),
        "assembled_at": assembled.get("assembled_at"),
        "cassette_count": len(files),
        "total_length_bp": sum(int(item.get("length_bp") or 0) for item in files),
        "all_files_exist": all(bool(item.get("exists")) for item in files) if files else False,
        "cassette_files": files,
    }


def _plasmid_selection(manifest: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    plasmid_selection = _as_dict(manifest.get("plasmid_selection"))
    if not plasmid_selection:
        return {"available": False}

    vector = _as_dict(plasmid_selection.get("vector"))
    sequence_file = _as_dict(vector.get("sequence_file"))
    annotation_files = _as_dict(vector.get("annotation_files"))

    vector_summary = {
        "plasmid_id": vector.get("plasmid_id", ""),
        "addgene_id": vector.get("addgene_id", ""),
        "name": vector.get("name", ""),
        "description": vector.get("description", ""),
        "vector_type": vector.get("vector_type", ""),
        "bacterial_resistance": vector.get("bacterial_resistance", ""),
        "growth_strain": vector.get("growth_strain", ""),
        "copy_number": vector.get("copy_number", ""),
        "copy_number_class": vector.get("copy_number_class", ""),
        "replicon_family": vector.get("replicon_family", ""),
        "cargo_type": vector.get("cargo_type", ""),
        "assembly_policy": vector.get("assembly_policy", ""),
        "requires_cargo_replacement": vector.get("requires_cargo_replacement"),
        "has_expression_cargo": vector.get("has_expression_cargo"),
        "host_compatibility": _as_list(vector.get("host_compatibility")),
        "replication_dependencies": _as_list(vector.get("replication_dependencies")),
        "regulatory_contexts": _as_list(vector.get("regulatory_contexts")),
        "insertion_regions": _as_list(vector.get("insertion_regions")),
        "module_boundaries": _as_dict(vector.get("module_boundaries")),
        "audit": _as_dict(vector.get("audit")),
        "cloning_method": vector.get("cloning_method", ""),
        "backbone": vector.get("backbone", ""),
        "length_bp": vector.get("length_bp"),
        "topology": vector.get("topology", ""),
        "sequence_file": _file_summary(output_dir, sequence_file) if sequence_file else {},
        "features": _as_dict(vector.get("features")),
        "annotation_files": {
            key: _file_summary(output_dir, {"path": value})
            for key, value in annotation_files.items()
            if value
        },
    }

    return {
        "available": True,
        "status": plasmid_selection.get("status"),
        "source": plasmid_selection.get("source"),
        "source_kind": plasmid_selection.get("source_kind"),
        "selected_at": plasmid_selection.get("selected_at"),
        "vector": vector_summary,
        "recommendation": _as_dict(plasmid_selection.get("recommendation")),
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


def _existing_final_assembly(manifest: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    final_assembly = _as_dict(manifest.get("final_assembly"))
    if not final_assembly:
        return {"available": False}

    construct_files = _as_dict(final_assembly.get("construct_files"))
    return {
        "available": True,
        "status": final_assembly.get("status"),
        "source": final_assembly.get("source"),
        "assembled_at": final_assembly.get("assembled_at"),
        "strategy": final_assembly.get("strategy", {}),
        "construct_files": {
            key: _file_summary(output_dir, value if isinstance(value, dict) else {"path": value})
            for key, value in construct_files.items()
        },
        "length_bp": final_assembly.get("length_bp"),
        "warnings": _as_list(final_assembly.get("warnings")),
    }


def _missing_requirements(assembled: dict[str, Any], plasmid: dict[str, Any]) -> list[str]:
    missing = []

    if not assembled.get("available"):
        missing.append("missing assembled_expression_cassettes")
    elif assembled.get("status") != "assembled":
        missing.append('assembled_expression_cassettes.status is not "assembled"')
    elif int(assembled.get("cassette_count") or 0) <= 0:
        missing.append("assembled_expression_cassettes.cassette_files is empty")
    elif not assembled.get("all_files_exist"):
        missing.append("one or more expression cassette files do not exist")

    vector = _as_dict(plasmid.get("vector"))
    vector_sequence = _as_dict(vector.get("sequence_file"))
    if not plasmid.get("available"):
        missing.append("missing plasmid_selection")
    elif plasmid.get("status") != "selected":
        missing.append('plasmid_selection.status is not "selected"')
    elif not vector:
        missing.append("plasmid_selection.vector is missing")
    elif not vector_sequence.get("path"):
        missing.append("plasmid_selection.vector.sequence_file.path is missing")
    elif not vector_sequence.get("exists"):
        missing.append("plasmid_selection.vector.sequence_file.path does not exist")

    return missing


def _final_assembly_requirements(assembled: dict[str, Any], plasmid: dict[str, Any]) -> dict[str, Any]:
    insert_total_length = int(assembled.get("total_length_bp") or 0)
    cassette_count = int(assembled.get("cassette_count") or 0)
    vector = _as_dict(plasmid.get("vector"))
    vector_length = int(vector.get("length_bp") or _as_dict(vector.get("sequence_file")).get("length_bp") or 0)
    estimated_final_length = vector_length + insert_total_length if vector_length and insert_total_length else None
    assembly_policy = str(vector.get("assembly_policy") or "")
    insertion_regions = [
        region
        for region in _as_list(vector.get("insertion_regions"))
        if isinstance(region, dict)
    ]
    audited_region = insertion_regions[0] if insertion_regions else {}

    warnings = []
    notes = []
    if cassette_count > 1:
        notes.append("存在多个表达盒，最终组装需要确认表达盒顺序和连接方式。")
    if insert_total_length >= 10000:
        warnings.append("表达盒插入片段较长，最终载体稳定性和拷贝数需要人工复核。")
    if not vector.get("topology"):
        warnings.append("质粒 topology 未登记，最终组装前建议确认 circular/linear。")
    if assembly_policy == "replace_seva_cargo_paci_spei":
        notes.append("该载体必须整体替换已审计的 PacI–SpeI cargo，不能保留 lacZα 克隆筛选上下文。")
        if not audited_region:
            warnings.append("载体要求替换 SEVA cargo，但缺少已审计 insertion_region。")
    elif not _as_dict(vector.get("features")).get("mcs_features"):
        notes.append("未发现已登记的 MCS/insertion site 信息，后续装配计划可能需要用户指定插入位置。")
    if vector.get("has_expression_cargo") is True:
        warnings.append("载体自身含表达 cargo；与完整上游表达盒组合前必须移除或人工复核。")
    if assembly_policy == "context_specific_review":
        warnings.append("该载体的插入策略未完全标准化，不能自动确认装配方案。")

    requires_user_confirmation = [
        item
        for item in [
            "cassette_order" if cassette_count > 1 else "",
            (
                "assembly_method_for_audited_cargo_replacement"
                if assembly_policy == "replace_seva_cargo_paci_spei"
                else "insertion_site_or_cloning_strategy"
            ),
            "whether_to_replace_existing_insert",
        ]
        if item
    ]

    return {
        "insert_total_length_bp": insert_total_length,
        "cassette_count": cassette_count,
        "vector_length_bp": vector_length,
        "estimated_final_length_bp": estimated_final_length,
        "assembly_policy": assembly_policy,
        "audited_insertion_region": audited_region,
        "recommended_strategy": (
            "replace the complete SEVA PacI-SpeI cargo with the assembled expression cassette"
            if assembly_policy == "replace_seva_cargo_paci_spei"
            else "insert assembled expression cassette files into the selected vector backbone"
        ),
        "requires_user_confirmation": requires_user_confirmation,
        "compatibility_notes": [note for note in notes if note],
        "warnings": warnings,
    }


def _next_actions(ready: bool, missing: list[str]) -> list[str]:
    if ready:
        return [
            "向用户确认最终组装策略：表达盒顺序、插入位点/克隆方式，以及是否替换 vector 原有 insert。",
            "用户确认策略后，调用 plan_final_assembly 写入 final_assembly_plan。",
            "final_assembly_plan 写入完成后，再交给 final_assemble_execute_agent 执行最终组装。"

        ]

    actions = []
    if any("assembled_expression_cassettes" in item or "cassette" in item for item in missing):
        actions.append("先回到 expression_cassette_assembly_agent 完成表达盒 GenBank 组装。")
    if any("plasmid_selection" in item or "vector" in item for item in missing):
        actions.append("先回到 plasmid_selection_agent 明确选择可用 vector，并登记 sequence_file。")
    return actions or ["补齐缺失的上游设计信息后，再重新读取最终组装上下文。"]


@tool(args_schema=GetFinalAssemblyContextArgs)
def get_final_assembly_context(
) -> str:
    """
    读取最终组装计划所需的表达盒和质粒上下文。
    
    调用时机：制定 final_assembly_plan 或推荐 Gibson/restriction 前。
    返回：ok、ready_for_final_assembly、missing_requirements、表达盒/vector 摘要、requirements 和 next_actions。
    限制：只读；不写计划、不生成最终质粒、不重新选择上游设计。
    """

    tool_name = "get_final_assembly_context"
    monitor.report_start(tool_name)
    try:
        session_root = _session_root(None, None, None)
        outputs = _output_dir(session_root, None)
        manifest_file = _manifest_file(session_root, outputs, None)
        manifest = _read_json(manifest_file)

        assembled = _assembled_expression_cassettes(manifest, outputs)
        plasmid = _plasmid_selection(manifest, outputs)
        missing = _missing_requirements(assembled, plasmid)
        ready = len(missing) == 0

        message = (
            "表达盒和质粒选择均已完成，可以开始最终组装。"
            if ready
            else "最终组装所需信息尚未齐全，请先补齐缺失的表达盒或质粒选择信息。"
        )

        result = {
            "ok": True,
            "ready_for_final_assembly": ready,
            "message": message,
            "manifest_path": _relpath_or_abs(session_root, manifest_file),
            "output_dir": _relpath_or_abs(session_root, outputs),
            "revision": manifest.get("revision"),
            "host_context": _host_context(manifest),
            "assembled_expression_cassettes": assembled,
            "plasmid_selection": plasmid,
            "final_assembly_requirements": _final_assembly_requirements(assembled, plasmid),
            "existing_final_assembly": _existing_final_assembly(manifest, outputs),
            "missing_requirements": missing,
            "next_actions": _next_actions(ready, missing),
        }
        monitor.report_end(tool_name, {"ready_for_final_assembly": ready})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(get_final_assembly_context.invoke({}))
