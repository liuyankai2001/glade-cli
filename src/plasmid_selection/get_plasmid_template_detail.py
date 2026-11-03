from __future__ import annotations

import json
import re
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.plasmid_selection.config import PLASMID_COLLECTION
from src.plasmid_selection.milvus_plasmids import _client_from_environment
from src.plasmid_selection.models import without_audit_hashes


# 兼容工具参数默认值；环境变量与默认库名集中在 service_config.py。
COLLECTION_NAME = PLASMID_COLLECTION

OUTPUT_FIELDS = [
    "plasmid_id",
    "addgene_id",
    "name",
    "description",
    "depositor",
    "article",
    "pubmed_id",
    "vector_type",
    "bacterial_resistance",
    "growth_strain",
    "copy_number",
    "cloning_method",
    "backbone",
    "insert_name",
    "insert_species",
    "sequence_file",
    "length_bp",
    "topology",
    "features",
    "origins",
    "resistance_markers",
    "promoters",
    "terminators",
    "mcs_features",
    "source",
    "source_url",
    "sequence_url",
    "text",
    "created_at",
    "updated_at",
    "schema_version",
    "source_record_id",
    "source_record_version",
    "replicon_family",
    "copy_number_class",
    "cargo_type",
    "assembly_policy",
    "audit_status",
    "audit_version",
    "audit_passed",
    "mg1655_compatible",
    "has_expression_cargo",
    "requires_cargo_replacement",
    "normalization_applied",
    "source_provenance",
    "host_compatibility",
    "replication_dependencies",
    "selection_markers",
    "regulatory_contexts",
    "insertion_regions",
    "module_boundaries",
    "protected_features",
    "evidence_refs",
    "audit_checks",
    "normalization_events",
]


class GetPlasmidTemplateDetailArgs(BaseModel):
    plasmid_id: str | None = Field(
        default=None,
        description="质粒模板主键，例如 addgene:50005。优先使用该字段。",
    )
    addgene_id: str | None = Field(default=None, description="Addgene 编号，例如 50005。")
    name: str | None = Field(default=None, description="质粒名称，例如 pUC19。")
    collection_name: str = Field(default=COLLECTION_NAME, description="Milvus collection 名称")


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _milvus_client():
    return _client_from_environment()


def _escape_expr_string(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _name_variants(name: str) -> list[str]:
    text = name.strip()
    if not text:
        return []
    variants = [text]
    collapsed = re.sub(r"[\s_-]+", "", text)
    if collapsed and collapsed not in variants:
        variants.append(collapsed)
    lower = text.lower()
    if lower not in variants:
        variants.append(lower)
    return variants


def _build_expr(plasmid_id: str | None, addgene_id: str | None, name: str | None) -> str:
    clauses = []
    if plasmid_id:
        clauses.append(f'plasmid_id == "{_escape_expr_string(plasmid_id.strip())}"')
    if addgene_id:
        clauses.append(f'addgene_id == "{_escape_expr_string(addgene_id.strip())}"')
        addgene_plasmid_id = f"addgene:{addgene_id.strip()}"
        clauses.append(f'plasmid_id == "{_escape_expr_string(addgene_plasmid_id)}"')
    if name:
        for variant in _name_variants(name):
            clauses.append(f'name == "{_escape_expr_string(variant)}"')
    if not clauses:
        raise ValueError("missing plasmid_id, addgene_id, or name")
    return " or ".join(f"({clause})" for clause in clauses)


def _query_collection(collection_name: str, expr: str) -> list[dict[str, Any]]:
    client = _milvus_client()
    if not client.has_collection(collection_name):
        raise ValueError(f"Milvus collection not found: {collection_name}")
    client.load_collection(collection_name)
    rows = client.query(
        collection_name=collection_name,
        filter=expr,
        output_fields=OUTPUT_FIELDS,
        limit=20,
    )
    cleaned = []
    for row in rows:
        values = dict(row)
        for key, value in values.items():
            if isinstance(value, str) and value.lstrip().startswith(("{", "[")):
                try:
                    values[key] = json.loads(value)
                except json.JSONDecodeError:
                    pass
        cleaned.append(without_audit_hashes(values))
    return cleaned


def _summarize_candidates(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "plasmid_id": row.get("plasmid_id", ""),
            "addgene_id": row.get("addgene_id", ""),
            "name": row.get("name", ""),
            "description": row.get("description", ""),
            "bacterial_resistance": row.get("bacterial_resistance", ""),
            "copy_number": row.get("copy_number", ""),
            "copy_number_class": row.get("copy_number_class", ""),
            "replicon_family": row.get("replicon_family", ""),
            "cargo_type": row.get("cargo_type", ""),
            "assembly_policy": row.get("assembly_policy", ""),
            "audit_status": row.get("audit_status", ""),
            "source_url": row.get("source_url", ""),
        }
        for row in rows
    ]


@tool(args_schema=GetPlasmidTemplateDetailArgs)
def get_plasmid_template_detail(
    plasmid_id: str | None = None,
    addgene_id: str | None = None,
    name: str | None = None,
    collection_name: str = COLLECTION_NAME,
) -> str:
    """
    按 plasmid_id、Addgene ID 或名称查询质粒模板详情。
    
    调用时机：recommend_plasmid_templates 返回候选后，需要查看某个模板细节。
    返回：ok、模板元数据、序列/注释路径、features 摘要和候选匹配。
    限制：只读；不写 plasmid_selection；找不到唯一模板时返回候选供用户确认。
    """

    try:
        expr = _build_expr(plasmid_id, addgene_id, name)
        rows = _query_collection(collection_name, expr)
        if not rows:
            result = {
                "ok": True,
                "found": False,
                "message": "未找到匹配的质粒模板。",
                "query": {
                    "plasmid_id": plasmid_id,
                    "addgene_id": addgene_id,
                    "name": name,
                },
                "collection_name": collection_name,
            }
            return json.dumps(result, ensure_ascii=False, default=str)

        exact_rows = rows
        if plasmid_id:
            exact_rows = [row for row in rows if row.get("plasmid_id") == plasmid_id] or rows
        elif addgene_id:
            exact_rows = [row for row in rows if str(row.get("addgene_id") or "") == str(addgene_id)] or rows
        elif name:
            exact_rows = [row for row in rows if str(row.get("name") or "").lower() == name.lower()] or rows

        plasmid = exact_rows[0]
        result = {
            "ok": True,
            "found": True,
            "collection_name": collection_name,
            "plasmid": plasmid,
            "candidate_count": len(rows),
            "other_candidates": _summarize_candidates(rows[1:]),
            "next_actions": [
                "如果用户认可该质粒，可继续调用 write_plasmid_selection 写入 plasmid_selection。",
                "如果需要检查用户上传或本地序列，可调用 annotate_plasmid_sequence。",
            ],
        }
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(get_plasmid_template_detail.invoke({"plasmid_id": "addgene:50005"}))
