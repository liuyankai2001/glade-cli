from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.config.global_config import BGE_M3_MODEL_PATH
from src.config.service_config import (
    get_milvus_config,
    get_plasmid_collection_name,
    get_plasmid_embedding_model_path,
)
from src.runtime.monitor import monitor

if __package__:
    from .get_plasmid_context import get_plasmid_context
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from src.tools.plasmid_selection_tools.get_plasmid_context import get_plasmid_context


# 兼容工具参数默认值；环境变量与字段说明集中在 service_config.py。
COLLECTION_NAME = get_plasmid_collection_name()
DEFAULT_MODEL_PATH = BGE_M3_MODEL_PATH

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
    "sequence_sha256",
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
    "source_file_sha256",
    "sequence_content_sha256",
    "canonical_sequence_sha256",
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

_EMBEDDING_MODEL: Any | None = None


class RecommendPlasmidTemplatesArgs(BaseModel):
    query: str | None = Field(
        default=None,
        description="用户额外偏好，例如低拷贝、Kan抗性、E. coli、广宿主等。",
    )
    top_k: int = Field(default=5, ge=1, le=20, description="最终返回候选数量")
    search_k: int = Field(default=20, ge=1, le=100, description="Milvus 初始召回数量")
    preferred_resistance: str | None = Field(default=None, description="偏好的抗性标记")
    excluded_resistance: list[str] = Field(default_factory=list, description="需要排除的抗性标记")
    preferred_copy_number: str | None = Field(default=None, description="偏好的拷贝数，例如 Low Copy")
    collection_name: str = Field(default=COLLECTION_NAME, description="Milvus collection 名称")
    allow_without_context: bool = Field(
        default=False,
        description="没有可用 design_manifest 上下文时，是否允许只按 query 检索质粒模板。",
    )


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _import_sentence_transformers():
    try:
        from sentence_transformers import SentenceTransformer
    except Exception as exc:
        raise RuntimeError(
            "sentence-transformers is not importable. If the error mentions torchvision, "
            "remove or fix torchvision because bge-m3 text embedding does not require it."
        ) from exc
    return SentenceTransformer


def _import_milvus_client():
    try:
        from pymilvus import MilvusClient
    except ImportError as exc:
        raise RuntimeError("pymilvus is not installed or not importable.") from exc
    return MilvusClient


def _milvus_uri() -> str:
    return get_milvus_config().uri


def _milvus_client():
    milvus_client = _import_milvus_client()
    kwargs = get_milvus_config().client_kwargs()
    kwargs["uri"] = _milvus_uri()
    return milvus_client(**kwargs)


def _embedding_model():
    global _EMBEDDING_MODEL
    if _EMBEDDING_MODEL is None:
        model_path = get_plasmid_embedding_model_path()
        sentence_transformer = _import_sentence_transformers()
        _EMBEDDING_MODEL = sentence_transformer(str(model_path), trust_remote_code=True)
    return _EMBEDDING_MODEL


def _embed_query(text: str) -> list[float]:
    model = _embedding_model()
    embedding = model.encode(text, normalize_embeddings=True)
    return embedding.tolist()


def _tool_json(tool_result: Any) -> dict[str, Any]:
    if isinstance(tool_result, dict):
        return tool_result
    if isinstance(tool_result, str):
        data = json.loads(tool_result)
        if not isinstance(data, dict):
            raise ValueError("tool result JSON must be an object")
        return data
    raise ValueError(f"unsupported tool result type: {type(tool_result).__name__}")


def _load_context(
) -> dict[str, Any]:
    return _tool_json(get_plasmid_context.invoke({}))


def _empty_context() -> dict[str, Any]:
    return {
        "ok": True,
        "ready_for_plasmid_selection": True,
        "message": "未读取 design_manifest，仅根据用户 query 检索质粒模板。",
        "host_context": {},
        "parts_selection": {"available": False},
        "assembled_expression_cassettes": {
            "available": False,
            "cassette_count": 0,
            "total_length_bp": 0,
            "cassette_files": [],
        },
        "vector_requirements": {
            "insert_total_length_bp": 0,
            "cassette_count": 0,
            "compatibility_notes": [],
            "warnings": [],
        },
    }


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _lower_text(*values: Any) -> str:
    parts = []
    for value in values:
        if isinstance(value, list):
            parts.extend(str(item) for item in value)
        else:
            parts.append(str(value or ""))
    return " ".join(parts).lower()


def _marker_text(plasmid: dict[str, Any]) -> str:
    normalized = _as_list(plasmid.get("selection_markers"))
    return _lower_text(
        plasmid.get("bacterial_resistance"),
        _as_list(plasmid.get("resistance_markers")),
        [
            " ".join(
                str(marker.get(key) or "")
                for key in ("phenotype", "family", "gene")
            )
            for marker in normalized
            if isinstance(marker, dict)
        ],
    )


def _is_mg1655_context(context: dict[str, Any]) -> bool:
    host = _as_dict(context.get("host_context"))
    host_name = _lower_text(
        host.get("chassis_key"),
        _as_dict(host.get("codon_optimization_host")).get("name"),
    )
    return any(
        token in host_name
        for token in ("mg1655", "escherichia coli", "e. coli", "ecoli")
    )


def _hard_filter_reason(
    plasmid: dict[str, Any],
    *,
    context: dict[str, Any],
    excluded_resistance: list[str],
    require_v2_audit: bool,
) -> str | None:
    if require_v2_audit:
        if plasmid.get("schema_version") != "plasmid_template.v2":
            return "record is not a plasmid_template.v2 card"
        if plasmid.get("audit_status") != "PASS" or plasmid.get("audit_passed") is not True:
            return "sequence audit did not pass"
    if _is_mg1655_context(context) and plasmid.get("mg1655_compatible") is False:
        return "replicon is not compatible with E. coli K-12 MG1655"

    resistance_text = _marker_text(plasmid)
    for marker in excluded_resistance:
        marker_text = marker.strip().lower()
        if marker_text and marker_text in resistance_text:
            return f"excluded resistance marker: {marker}"

    assembled = _as_dict(context.get("assembled_expression_cassettes"))
    parts = _as_dict(context.get("parts_selection"))
    expression_already_regulated = bool(parts.get("available") or assembled.get("available"))
    if expression_already_regulated and plasmid.get("has_expression_cargo") is True:
        return "vector expression cargo conflicts with a complete upstream expression cassette"

    if plasmid.get("assembly_policy") == "replace_seva_cargo_paci_spei":
        if not _as_list(plasmid.get("insertion_regions")):
            return "SEVA cargo replacement policy is missing an audited insertion region"
    return None


def _build_query_text(context: dict[str, Any], user_query: str | None) -> str:
    assembled = _as_dict(context.get("assembled_expression_cassettes"))
    requirements = _as_dict(context.get("vector_requirements"))
    host = _as_dict(context.get("host_context"))
    parts = _as_dict(context.get("parts_selection"))

    insert_length = requirements.get("insert_total_length_bp") or assembled.get("total_length_bp") or 0
    cassette_count = requirements.get("cassette_count") or assembled.get("cassette_count") or 0
    chassis = host.get("chassis_key") or _as_dict(host.get("codon_optimization_host")).get("name") or ""
    already_has_parts = parts.get("available") or assembled.get("available")

    lines = [
        "Recommend a plasmid vector template for synthetic biology construct assembly.",
        f"Insert total length: {insert_length} bp.",
        f"Expression cassette count: {cassette_count}.",
        f"Host or chassis: {chassis}.",
        "The vector should provide replication origin, selectable marker, and cloning or insertion sites.",
    ]
    if already_has_parts:
        lines.append(
            "The expression cassette already contains promoter, RBS, CDS, and terminator; prefer backbone compatibility and avoid regulatory conflicts."
        )
    notes = _as_list(requirements.get("compatibility_notes"))
    if notes:
        lines.append("Compatibility notes: " + " ".join(str(note) for note in notes))
    warnings = _as_list(requirements.get("warnings"))
    if warnings:
        lines.append("Risk notes: " + " ".join(str(warning) for warning in warnings))
    if user_query:
        lines.append(f"User preference: {user_query}")
    return " ".join(lines)


def _hit_entity(hit: Any) -> tuple[dict[str, Any], float]:
    if isinstance(hit, dict):
        entity = hit.get("entity") or hit
        distance = hit.get("distance", hit.get("score", 0.0))
        return dict(entity), float(distance or 0.0)
    entity = getattr(hit, "entity", {}) or {}
    distance = getattr(hit, "distance", getattr(hit, "score", 0.0))
    return dict(entity), float(distance or 0.0)


def _milvus_search(collection_name: str, query_embedding: list[float], search_k: int) -> list[tuple[dict[str, Any], float]]:
    client = _milvus_client()
    if not client.has_collection(collection_name):
        raise ValueError(f"Milvus collection not found: {collection_name}")
    client.load_collection(collection_name)
    result = client.search(
        collection_name=collection_name,
        data=[query_embedding],
        anns_field="embedding",
        limit=search_k,
        output_fields=OUTPUT_FIELDS,
        search_params={"metric_type": "COSINE"},
    )
    if not result:
        return []
    return [_hit_entity(hit) for hit in result[0]]


def _rule_score(
    plasmid: dict[str, Any],
    *,
    context: dict[str, Any],
    preferred_resistance: str | None,
    excluded_resistance: list[str],
    preferred_copy_number: str | None,
) -> tuple[float, list[str], list[str]]:
    requirements = _as_dict(context.get("vector_requirements"))
    assembled = _as_dict(context.get("assembled_expression_cassettes"))
    parts = _as_dict(context.get("parts_selection"))
    insert_length = int(requirements.get("insert_total_length_bp") or assembled.get("total_length_bp") or 0)

    score = 0.0
    rationale: list[str] = []
    warnings: list[str] = []

    sequence_file = str(plasmid.get("sequence_file") or "").strip()
    origins = _as_list(plasmid.get("origins"))
    resistance_markers = _as_list(plasmid.get("resistance_markers"))
    promoters = _as_list(plasmid.get("promoters"))
    mcs_features = _as_list(plasmid.get("mcs_features"))
    bacterial_resistance = str(plasmid.get("bacterial_resistance") or "")
    copy_number = str(
        plasmid.get("copy_number_class")
        or plasmid.get("copy_number")
        or ""
    )
    vector_type = str(plasmid.get("vector_type") or "")
    cargo_type = str(plasmid.get("cargo_type") or "")
    assembly_policy = str(plasmid.get("assembly_policy") or "")

    if sequence_file:
        score += 10
        rationale.append("有本地 GenBank 序列文件，可继续注释和装配检查。")
    if origins:
        score += 10
        rationale.append("已识别复制起点。")
    if resistance_markers or bacterial_resistance:
        score += 10
        rationale.append("有筛选标记。")
    if mcs_features:
        score += 8
        rationale.append("有 MCS/克隆位点，便于插入表达盒。")

    lower_copy = copy_number.lower()
    if insert_length >= 5000:
        if "low" in lower_copy:
            score += 15
            rationale.append("表达盒较长，低拷贝载体更利于稳定性。")
        elif "high" in lower_copy:
            score -= 10
            warnings.append("表达盒较长，高拷贝载体可能增加代谢负担。")
    elif 0 < insert_length < 3000 and "high" in lower_copy:
        score += 5
        rationale.append("插入片段较短，高拷贝载体适合快速克隆验证。")

    if preferred_copy_number and preferred_copy_number.lower() in lower_copy:
        score += 10
        rationale.append(f"符合偏好的拷贝数：{preferred_copy_number}。")

    resistance_text = _marker_text(plasmid)
    if preferred_resistance and preferred_resistance.lower() in resistance_text:
        score += 12
        rationale.append(f"符合偏好的抗性标记：{preferred_resistance}。")

    expression_already_regulated = bool(parts.get("available") or assembled.get("available"))
    if expression_already_regulated and plasmid.get("has_expression_cargo") is True:
        score -= 100
        warnings.append("载体含表达 cargo，与上游完整表达盒冲突。")
    elif expression_already_regulated and promoters:
        score -= 3
        warnings.append("载体含非表达 cargo 的局部 promoter，需按调控角色检查读穿风险。")
    if "expression" in vector_type.lower() and expression_already_regulated:
        score -= 2
        warnings.append("该模板标记为表达载体，若表达盒已自带调控元件需复核装配策略。")

    if cargo_type == "mcs_default" and assembly_policy == "insert_into_mcs":
        score += 10
        rationale.append("载体提供标准 MCS，可按已审计插入区装配完整表达盒。")
    elif assembly_policy == "replace_seva_cargo_paci_spei":
        score += 8
        rationale.append("可整体替换 PacI–SpeI cargo，移除 lacZα 克隆筛选上下文。")
    elif assembly_policy == "context_specific_review":
        score -= 15
        warnings.append("该旧载体需要人工确认插入区或原有调控上下文。")

    if plasmid.get("audit_passed") is True:
        score += 12
        rationale.append("完整序列与标准化生物学字段已通过审计。")

    if not rationale:
        rationale.append("Milvus 语义检索与当前需求匹配。")
    return score, rationale, warnings


def _candidate_payload(
    plasmid: dict[str, Any],
    *,
    rank: int,
    semantic_score: float,
    rule_score: float,
    rationale: list[str],
    warnings: list[str],
) -> dict[str, Any]:
    return {
        "rank": rank,
        "plasmid_id": plasmid.get("plasmid_id", ""),
        "addgene_id": plasmid.get("addgene_id", ""),
        "name": plasmid.get("name", ""),
        "description": plasmid.get("description", ""),
        "final_score": round(semantic_score * 60 + rule_score, 3),
        "semantic_score": round(semantic_score, 4),
        "rule_score": round(rule_score, 3),
        "vector_type": plasmid.get("vector_type", ""),
        "bacterial_resistance": plasmid.get("bacterial_resistance", ""),
        "copy_number": plasmid.get("copy_number", ""),
        "copy_number_class": plasmid.get("copy_number_class", ""),
        "replicon_family": plasmid.get("replicon_family", ""),
        "cargo_type": plasmid.get("cargo_type", ""),
        "assembly_policy": plasmid.get("assembly_policy", ""),
        "audit_status": plasmid.get("audit_status", ""),
        "audit_passed": plasmid.get("audit_passed"),
        "mg1655_compatible": plasmid.get("mg1655_compatible"),
        "has_expression_cargo": plasmid.get("has_expression_cargo"),
        "requires_cargo_replacement": plasmid.get("requires_cargo_replacement"),
        "backbone": plasmid.get("backbone", ""),
        "insert_name": plasmid.get("insert_name", ""),
        "sequence_file": plasmid.get("sequence_file", ""),
        "sequence_sha256": plasmid.get("sequence_sha256", ""),
        "length_bp": plasmid.get("length_bp"),
        "topology": plasmid.get("topology", ""),
        "origins": _as_list(plasmid.get("origins")),
        "resistance_markers": _as_list(plasmid.get("resistance_markers")),
        "promoters": _as_list(plasmid.get("promoters")),
        "terminators": _as_list(plasmid.get("terminators")),
        "mcs_features": _as_list(plasmid.get("mcs_features")),
        "selection_markers": _as_list(plasmid.get("selection_markers")),
        "regulatory_contexts": _as_list(plasmid.get("regulatory_contexts")),
        "insertion_regions": _as_list(plasmid.get("insertion_regions")),
        "host_compatibility": _as_list(plasmid.get("host_compatibility")),
        "replication_dependencies": _as_list(plasmid.get("replication_dependencies")),
        "evidence_refs": _as_list(plasmid.get("evidence_refs")),
        "source_record_id": plasmid.get("source_record_id", ""),
        "source_record_version": plasmid.get("source_record_version", ""),
        "source_url": plasmid.get("source_url", ""),
        "sequence_url": plasmid.get("sequence_url", ""),
        "rationale": rationale,
        "warnings": warnings,
    }


@tool(args_schema=RecommendPlasmidTemplatesArgs)
def recommend_plasmid_templates(
    query: str | None = None,
    top_k: int = 3,
    search_k: int = 20,
    preferred_resistance: str | None = None,
    excluded_resistance: list[str] | None = None,
    preferred_copy_number: str | None = None,
    collection_name: str = COLLECTION_NAME,
    allow_without_context: bool = False,
) -> str:
    """
    根据当前表达盒和宿主需求推荐质粒模板。
    
    调用时机：plasmid_selection 尚未确定，用户要求推荐 vector/backbone。
    输入：query、top_k、过滤偏好。
    返回：ok、context 摘要、ranked recommendations、get_detail/write_selection 参数和 warnings。
    限制：只读；用户确认模板后再调用 write_plasmid_selection。
    """

    tool_name = "recommend_plasmid_templates"
    monitor.report_start(tool_name, {"query": query, "top_k": top_k})
    try:
        excluded = excluded_resistance or []
        try:
            context = _load_context()
        except Exception:
            if not allow_without_context:
                raise
            context = _empty_context()

        if not context.get("ok") and allow_without_context:
            context = _empty_context()
        if not context.get("ok"):
            return json.dumps(context, ensure_ascii=False, default=str)
        if not context.get("ready_for_plasmid_selection") and not allow_without_context:
            result = {
                "ok": True,
                "ready_for_plasmid_selection": False,
                "message": context.get("message") or "表达盒尚未准备好，暂不推荐质粒模板。",
                "context": context,
                "candidates": [],
            }
            return json.dumps(result, ensure_ascii=False, default=str)

        query_text = _build_query_text(context, query)
        monitor.report_running(tool_name, "正在生成质粒模板检索向量...", progress=0.25)
        query_embedding = _embed_query(query_text)
        monitor.report_running(
            tool_name,
            f"正在检索 Milvus {collection_name}...",
            progress=0.55,
        )
        hits = _milvus_search(collection_name, query_embedding, max(search_k, top_k))

        ranked = []
        rejected: list[dict[str, str]] = []
        require_v2_audit = collection_name != "plasmid_templates"
        for plasmid, semantic_score in hits:
            filter_reason = _hard_filter_reason(
                plasmid,
                context=context,
                excluded_resistance=excluded,
                require_v2_audit=require_v2_audit,
            )
            if filter_reason:
                rejected.append(
                    {
                        "plasmid_id": str(plasmid.get("plasmid_id") or ""),
                        "name": str(plasmid.get("name") or ""),
                        "reason": filter_reason,
                    }
                )
                continue
            rule_score, rationale, warnings = _rule_score(
                plasmid,
                context=context,
                preferred_resistance=preferred_resistance,
                excluded_resistance=excluded,
                preferred_copy_number=preferred_copy_number,
            )
            ranked.append((
                semantic_score * 60 + rule_score,
                plasmid,
                semantic_score,
                rule_score,
                rationale,
                warnings,
            ))

        ranked.sort(key=lambda item: item[0], reverse=True)
        candidates = [
            _candidate_payload(
                plasmid,
                rank=index + 1,
                semantic_score=semantic_score,
                rule_score=rule_score,
                rationale=rationale,
                warnings=warnings,
            )
            for index, (_, plasmid, semantic_score, rule_score, rationale, warnings)
            in enumerate(ranked[:top_k])
        ]

        requirements = _as_dict(context.get("vector_requirements"))
        assembled = _as_dict(context.get("assembled_expression_cassettes"))
        result = {
            "ok": True,
            "ready_for_plasmid_selection": bool(context.get("ready_for_plasmid_selection")),
            "message": f"已基于 Milvus {collection_name} 和当前表达盒上下文生成质粒候选。",
            "collection_name": collection_name,
            "query_text": query_text,
            "insert_total_length_bp": requirements.get("insert_total_length_bp")
            or assembled.get("total_length_bp")
            or 0,
            "cassette_count": requirements.get("cassette_count")
            or assembled.get("cassette_count")
            or 0,
            "candidate_count": len(candidates),
            "candidates": candidates,
            "rejected_count": len(rejected),
            "rejected": rejected,
            "next_actions": [
                "向用户展示候选质粒、推荐理由和风险提示。",
                "如果用户想看某个质粒详情，可查询该 plasmid_id 的完整 Milvus 字段。",
                "只有用户明确选择/采用某个质粒后，才调用写入 plasmid_selection 的工具。",
            ],
        }
        monitor.report_end(tool_name, {"candidate_count": len(candidates)})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        try:
            monitor.report_error("recommend_plasmid_templates", exc)
        except Exception:
            pass
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":

    print(recommend_plasmid_templates.invoke({"allow_without_context": True, "query": "low copy E. coli plasmid"}))
