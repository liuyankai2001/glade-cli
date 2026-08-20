from __future__ import annotations

import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Literal

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.config.service_config import get_plasmid_collection_name
from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import (
    design_manifest_file,
    outputs_dir as resolve_outputs_dir,
    plasmid_annotation_dir as resolve_plasmid_annotation_dir,
    session_dir as resolve_session_dir,
)

# 与推荐和详情工具共享同一个默认 collection。
DEFAULT_COLLECTION_NAME = get_plasmid_collection_name()

if __package__:
    from .get_plasmid_template_detail import get_plasmid_template_detail
else:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from src.tools.plasmid_selection_tools.get_plasmid_template_detail import (
        get_plasmid_template_detail,
    )


class RecommendationArgs(BaseModel):
    score: float | None = Field(default=None, description="推荐分数，可来自 recommend_plasmid_templates")
    rationale: list[str] = Field(default_factory=list, description="推荐理由")
    warnings: list[str] = Field(default_factory=list, description="风险提示")


class WritePlasmidSelectionArgs(BaseModel):
    source_type: Literal["template", "annotated_upload"] = Field(
        description="template 表示 Milvus 推荐模板；annotated_upload 表示用户上传并已注释的质粒。"
    )
    expected_revision: int | None = Field(default=None, description="可选，用于防止覆盖旧版本")

    plasmid_id: str | None = Field(default=None, description="template 来源的 plasmid_id，例如 addgene:50005")
    addgene_id: str | None = Field(default=None, description="template 来源的 Addgene ID")
    name: str | None = Field(default=None, description="template 来源的质粒名称，或 upload 来源的自定义名称")
    collection_name: str = Field(default=DEFAULT_COLLECTION_NAME, description="Milvus collection 名称")

    annotation_result: dict[str, Any] | str | None = Field(
        default=None,
        description="annotated_upload 来源时，传入 annotate_plasmid_sequence 返回的 JSON 对象或 JSON 字符串。",
    )
    recommendation: RecommendationArgs | None = Field(default=None, description="可选推荐信息")


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
        return output_path.parent if output_path.name == "outputs" else output_path.parent
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


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "schema_version": "design_manifest.v1",
            "revision": 0,
        }
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return data


def _write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
        suffix=".tmp",
    ) as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        tmp_path = Path(handle.name)
    tmp_path.replace(path)


def _relpath_or_abs(base_dir: Path, path_value: Any) -> str:
    if not path_value:
        return ""
    path = Path(str(path_value))
    if not path.is_absolute():
        project_relative = (Path.cwd() / path).resolve()
        path = project_relative if project_relative.exists() else base_dir / path
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _load_tool_json(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        data = json.loads(value)
        if isinstance(data, dict):
            return data
    raise ValueError("expected JSON object or JSON string")


def _load_annotation_result(
    *,
    annotation_result: dict[str, Any] | str | None,
) -> dict[str, Any]:
    if annotation_result is not None:
        return _load_tool_json(annotation_result)

    annotation_dir = resolve_plasmid_annotation_dir()
    candidates = sorted(
        annotation_dir.glob("*_annotation_result.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise ValueError(
            "missing annotation_result and no annotation_result JSON found in current session plasmid_annotation output"
        )
    return _read_json(candidates[0].resolve())


def _template_detail(
    *,
    plasmid_id: str | None,
    addgene_id: str | None,
    name: str | None,
    collection_name: str,
) -> dict[str, Any]:
    detail = _load_tool_json(get_plasmid_template_detail.invoke({
        "plasmid_id": plasmid_id,
        "addgene_id": addgene_id,
        "name": name,
        "collection_name": collection_name,
    }))
    if not detail.get("ok"):
        raise ValueError(detail.get("message") or detail.get("error") or "get_plasmid_template_detail failed")
    if not detail.get("found"):
        raise ValueError("plasmid template not found")
    plasmid = detail.get("plasmid")
    if not isinstance(plasmid, dict):
        raise ValueError("plasmid detail missing plasmid object")
    return plasmid


def _vector_from_template(plasmid: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    if plasmid.get("schema_version") == "plasmid_template.v2":
        if plasmid.get("audit_status") != "PASS" or plasmid.get("audit_passed") is not True:
            raise ValueError("plasmid template did not pass the v2 sequence audit")
        if plasmid.get("mg1655_compatible") is not True:
            raise ValueError("plasmid template is not marked compatible with E. coli K-12 MG1655")

    sequence_file_path = _relpath_or_abs(output_dir, plasmid.get("sequence_file"))
    length_bp = int(plasmid.get("length_bp") or 0)
    return {
        "plasmid_id": plasmid.get("plasmid_id", ""),
        "addgene_id": plasmid.get("addgene_id", ""),
        "name": plasmid.get("name", ""),
        "description": plasmid.get("description", ""),
        "vector_type": plasmid.get("vector_type", ""),
        "bacterial_resistance": plasmid.get("bacterial_resistance", ""),
        "growth_strain": plasmid.get("growth_strain", ""),
        "copy_number": plasmid.get("copy_number", ""),
        "copy_number_class": plasmid.get("copy_number_class", ""),
        "replicon_family": plasmid.get("replicon_family", ""),
        "cargo_type": plasmid.get("cargo_type", ""),
        "assembly_policy": plasmid.get("assembly_policy", ""),
        "requires_cargo_replacement": plasmid.get("requires_cargo_replacement"),
        "has_expression_cargo": plasmid.get("has_expression_cargo"),
        "host_compatibility": _as_list(plasmid.get("host_compatibility")),
        "replication_dependencies": _as_list(plasmid.get("replication_dependencies")),
        "regulatory_contexts": _as_list(plasmid.get("regulatory_contexts")),
        "insertion_regions": _as_list(plasmid.get("insertion_regions")),
        "module_boundaries": _as_dict(plasmid.get("module_boundaries")),
        "cloning_method": plasmid.get("cloning_method", ""),
        "backbone": plasmid.get("backbone", ""),
        "insert_name": plasmid.get("insert_name", ""),
        "insert_species": plasmid.get("insert_species", ""),
        "length_bp": length_bp,
        "topology": plasmid.get("topology", ""),
        "sequence_file": {
            "path": sequence_file_path,
            "sha256": plasmid.get("sequence_sha256", ""),
            "sequence_content_sha256": plasmid.get("sequence_content_sha256", ""),
            "canonical_sequence_sha256": plasmid.get("canonical_sequence_sha256", ""),
            "length_bp": length_bp,
            "format": "genbank" if sequence_file_path.lower().endswith((".gb", ".gbk", ".gbj")) else "",
        },
        "source_url": plasmid.get("source_url", ""),
        "sequence_url": plasmid.get("sequence_url", ""),
        "features": {
            "origins": _as_list(plasmid.get("origins")),
            "resistance_markers": _as_list(plasmid.get("resistance_markers")),
            "promoters": _as_list(plasmid.get("promoters")),
            "terminators": _as_list(plasmid.get("terminators")),
            "mcs_features": _as_list(plasmid.get("mcs_features")),
            "selection_markers": _as_list(plasmid.get("selection_markers")),
            "protected_features": _as_list(plasmid.get("protected_features")),
        },
        "audit": {
            "schema_version": plasmid.get("schema_version", ""),
            "status": plasmid.get("audit_status", ""),
            "version": plasmid.get("audit_version", ""),
            "passed": plasmid.get("audit_passed"),
            "normalization_applied": plasmid.get("normalization_applied"),
            "checks": _as_dict(plasmid.get("audit_checks")),
            "normalization_events": _as_list(plasmid.get("normalization_events")),
        },
        "provenance": {
            "source_record_id": plasmid.get("source_record_id", ""),
            "source_record_version": plasmid.get("source_record_version", ""),
            "source": _as_dict(plasmid.get("source_provenance")),
            "evidence_refs": _as_list(plasmid.get("evidence_refs")),
        },
        "annotation_files": {},
    }


def _vector_from_annotation(annotation: dict[str, Any], output_dir: Path, override_name: str | None) -> dict[str, Any]:
    if not annotation.get("ok"):
        raise ValueError(annotation.get("message") or annotation.get("error") or "annotation_result is not ok")

    input_sequence = _as_dict(annotation.get("input_sequence_file"))
    annotation_files = _as_dict(annotation.get("annotation_files"))
    detected = _as_dict(annotation.get("detected_features"))
    vector_id = str(annotation.get("vector_id") or override_name or "uploaded_vector")
    length_bp = int(annotation.get("length_bp") or 0)
    input_path = input_sequence.get("path", "")
    return {
        "plasmid_id": f"upload:{vector_id}",
        "addgene_id": "",
        "name": override_name or vector_id,
        "description": "User uploaded plasmid sequence annotated by annotate_plasmid_sequence.",
        "vector_type": "",
        "bacterial_resistance": ", ".join(str(item) for item in _as_list(detected.get("resistance_markers"))),
        "growth_strain": "",
        "copy_number": "",
        "cloning_method": "",
        "backbone": "",
        "insert_name": "",
        "insert_species": "",
        "length_bp": length_bp,
        "gc_percent": annotation.get("gc_percent"),
        "topology": annotation.get("topology", ""),
        "sequence_file": {
            "path": _relpath_or_abs(output_dir, input_path),
            "sha256": input_sequence.get("sha256", ""),
            "length_bp": length_bp,
            "format": "genbank" if str(input_path).lower().endswith((".gb", ".gbk", ".gbj")) else "fasta",
        },
        "source_url": "",
        "sequence_url": "",
        "features": {
            "origins": _as_list(detected.get("origins")),
            "resistance_markers": _as_list(detected.get("resistance_markers")),
            "promoters": _as_list(detected.get("promoters")),
            "terminators": _as_list(detected.get("terminators")),
            "mcs_features": _as_list(detected.get("mcs")),
            "reporters": _as_list(detected.get("reporters")),
            "tags": _as_list(detected.get("tags")),
        },
        "annotation_files": {
            "genbank": _relpath_or_abs(output_dir, annotation_files.get("genbank")),
            "html_map": _relpath_or_abs(output_dir, annotation_files.get("html_map")),
            "features_csv": _relpath_or_abs(output_dir, annotation_files.get("features_csv")),
            "high_confidence_features_csv": _relpath_or_abs(
                output_dir,
                annotation_files.get("high_confidence_features_csv"),
            ),
        },
        "annotation_summary": {
            "database": annotation.get("database", ""),
            "feature_count": annotation.get("feature_count"),
            "high_confidence_feature_count": annotation.get("high_confidence_feature_count"),
            "warnings": _as_list(annotation.get("warnings")),
        },
    }


def _recommendation_payload(recommendation: RecommendationArgs | None) -> dict[str, Any]:
    if recommendation is None:
        return {
            "score": None,
            "rationale": [],
            "warnings": [],
        }
    return {
        "score": recommendation.score,
        "rationale": recommendation.rationale,
        "warnings": recommendation.warnings,
    }


@tool(args_schema=WritePlasmidSelectionArgs)
def write_plasmid_selection(
    source_type: Literal["template", "annotated_upload"],
    expected_revision: int | None = None,
    plasmid_id: str | None = None,
    addgene_id: str | None = None,
    name: str | None = None,
    collection_name: str = DEFAULT_COLLECTION_NAME,
    annotation_result: dict[str, Any] | str | None = None,
    recommendation: RecommendationArgs | None = None,
) -> str:
    """
    把用户确认的 vector/backbone 选择写入 design_manifest.json。
    
    调用时机：用户确认模板库候选或自定义质粒注释结果后。
    输入：source、template/annotation 信息、vector 名称和推荐依据。
    返回：ok、plasmid_selection 摘要、manifest_path 和 revision。
    限制：不推荐模板、不注释序列；只登记已确认 vector。
    """

    tool_name = "write_plasmid_selection"
    monitor.report_start(tool_name, {"source_type": source_type, "plasmid_id": plasmid_id, "name": name})
    try:
        root = _session_root(None, None, None)
        outputs = _output_dir(root, None)
        manifest_file = _manifest_file(root, outputs, None)
        manifest = _read_json(manifest_file)

        current_revision = int(manifest.get("revision", 0))
        if expected_revision is not None and expected_revision != current_revision:
            raise ValueError(
                f"manifest revision mismatch: expected {expected_revision}, current {current_revision}"
            )

        if source_type == "template":
            plasmid = _template_detail(
                plasmid_id=plasmid_id,
                addgene_id=addgene_id,
                name=name,
                collection_name=collection_name,
            )
            vector = _vector_from_template(plasmid, outputs)
            source_kind = "milvus_template"
        elif source_type == "annotated_upload":
            annotation = _load_annotation_result(
                annotation_result=annotation_result,
            )
            vector = _vector_from_annotation(annotation, outputs, name)
            source_kind = "annotated_upload"
        else:
            raise ValueError(f"unsupported source_type: {source_type}")

        manifest["plasmid_selection"] = {
            "status": "selected",
            "source": "write_plasmid_selection",
            "source_kind": source_kind,
            "selected_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "vector": vector,
            "recommendation": _recommendation_payload(recommendation),
        }
        clear_downstream_fields(manifest, "plasmid_selection")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_file, manifest)

        monitor.report_end(tool_name, {"source_kind": source_kind, "revision": manifest["revision"]})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": _relpath_or_abs(root, manifest_file),
                "revision": manifest["revision"],
                "source_kind": source_kind,
                "vector": {
                    "plasmid_id": vector.get("plasmid_id"),
                    "name": vector.get("name"),
                    "sequence_file": vector.get("sequence_file"),
                    "length_bp": vector.get("length_bp"),
                },
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        try:
            monitor.report_error("write_plasmid_selection", exc)
        except Exception:
            pass
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(write_plasmid_selection.invoke({
        "source_type": "template",
        "plasmid_id": "addgene:50005",
    }))
