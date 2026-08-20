from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Literal

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.config.expression_config import (
    EXPRESSION_RETRIEVAL_CONFIG,
    get_expression_embedding_model_path,
)
from src.config.global_config import BGE_M3_MODEL_PATH
from src.config.service_config import get_milvus_config
from src.runtime.monitor import monitor
from src.tools.common.session_paths import design_manifest_file, outputs_dir as resolve_outputs_dir, parts_dir as resolve_parts_dir
from src.tools.expression_cassette_assembly_tools.context_aware_rbs import (
    RegulatoryTarget,
    RelativeExpressionConstraint,
    context_aware_reselection,
)
from src.tools.expression_cassette_assembly_tools.context_aware_expression import (
    CassetteExpressionTarget,
    optimize_expression_regulatory_parts,
)
from src.tools.expression_cassette_assembly_tools.expression_host_context import (
    ExpressionHostContextError,
    ExpressionHostResolution,
    NoCompatibleExpressionPartsError,
    candidate_host_match_kind,
    resolve_expression_host,
)
from src.tools.expression_cassette_assembly_tools.part_identifiers import (
    normalize_part_id,
    part_sequence_path,
)


ALLOWED_ROLES = {"promoter", "rbs", "terminator"}
ALLOWED_STRENGTHS = {"low", "medium", "high", "unknown", "any"}
ALLOWED_REGULATIONS = {"constitutive", "inducible", "repressible", "unknown", "any"}
# 兼容已有测试、工具 schema 和外部导入；唯一配置源位于 expression_config.py。
COLLECTION_NAME = EXPRESSION_RETRIEVAL_CONFIG.collection_name
DEFAULT_SEARCH_K = EXPRESSION_RETRIEVAL_CONFIG.default_search_k
V3_EMBEDDING_DIM = EXPRESSION_RETRIEVAL_CONFIG.embedding_dim
DEFAULT_MODEL_PATH = BGE_M3_MODEL_PATH
CONFIDENCE_SCORES = {"high": 6, "medium": 3, "low": 1, "unknown": -2, "": -2}
UNSAFE_TERMINATOR_IDS = {"BBa_B1002", "BBa_B1006", "BBa_B1010"}
OUTPUT_FIELDS = [
    "doc_id",
    "part_id",
    "role",
    "role_confidence",
    "sequence_type",
    "sequence_type_confidence",
    "name",
    "description",
    "source",
    "evidence",
    "sequence",
    "sequence_available",
    "length_bp",
    "sequence_sha256",
    "gc_content",
    "host",
    "host_confidence",
    "strength",
    "strength_confidence",
    "regulation",
    "regulation_confidence",
    "warnings",
    "notes",
    "registry_metadata",
    "rag_text",
    "activity_value",
    "activity_sd",
    "activity_percentile",
    "activity_dataset",
    "activity_unit",
    "activity_context",
    "activity_source",
    "evidence_grade",
    "direction",
    "regulator_required",
    "quantitative_metadata",
]
V3_REQUIRED_FIELDS = {
    "part_id",
    "role",
    "sequence_type",
    "sequence",
    "sequence_available",
    "length_bp",
    "host",
    "strength",
    "regulation",
    "activity_value",
    "activity_percentile",
    "activity_dataset",
    "activity_unit",
    "evidence_grade",
    "direction",
    "regulator_required",
    "quantitative_metadata",
    "embedding",
}

_EMBEDDING_MODEL: Any | None = None


class CandidateExpressionPart(BaseModel):
    part_id: str = Field(description="BioBrick 或来源限定元件 ID，例如 BBa_J23106、KOSURI2013:P001")
    role: Literal["promoter", "rbs", "terminator"] = Field(description="元件角色")
    source: str = Field(default="user", description="候选来源，例如 web_search、milvus、user")
    strength: str = Field(default="unknown", description="表达强度：low、medium、high、unknown")
    regulation: str = Field(default="unknown", description="调控类型：constitutive、inducible、repressible、unknown")
    host: list[str] = Field(default_factory=list, description="适用宿主，例如 E. coli")
    evidence: str = Field(default="", description="证据链接或来源说明")
    name: str = Field(default="", description="可选；get_parts_facts 返回的元件名称")
    part_short_desc: str = Field(default="", description="可选；get_parts_facts 返回的元件简短描述")
    sequence_type: str = Field(default="", description="可选；get_parts_facts 返回的 sequence_type")
    notes: str = Field(default="", description="可选备注")
    activity_value: float | None = None
    activity_sd: float | None = None
    activity_percentile: float | None = None
    activity_dataset: str = ""
    activity_unit: str = ""
    activity_context: str = ""
    activity_source: str = ""
    evidence_grade: str = "C"
    direction: str = "unknown"
    regulator_required: str = "unknown"
    quantitative_metadata: dict[str, Any] = Field(default_factory=dict)


class RecommendExpressionPartsArgs(BaseModel):
    candidate_parts: list[CandidateExpressionPart] | None = Field(
        default=None,
        description="可选；用户提供的候选表达元件",
    )
    host: str | None = Field(
        default=None,
        description="可选宿主一致性断言；省略时只读继承 cds_selection.host，不覆盖上游宿主",
    )
    expression_strength: str = Field(
        default="medium",
        description="期望表达强度：low、medium、high 或 any",
    )
    promoter_regulation: str = Field(
        default="any",
        description="期望 promoter 调控类型：constitutive、inducible、repressible 或 any",
    )
    strategy: str = Field(
        default="balanced",
        description="推荐策略：balanced、high_expression、low_burden 或 custom",
    )
    alternatives_count: int = Field(default=0, ge=0, le=5, description="每类元件返回多少个备选")
    detail: Literal["compact", "full"] = Field(
        default="compact",
        description="返回详细程度：compact 精简输出；full 返回完整候选对象",
    )
    collection_name: str = Field(default=COLLECTION_NAME, description="Milvus expression_parts collection 名称")
    search_k: int = Field(
        default=DEFAULT_SEARCH_K,
        ge=1,
        le=100,
        description="每类 role 的 Milvus 初始召回数量；正式系统默认每类 Top-10",
    )
    context_rbs_limit: int = Field(
        default=24,
        ge=3,
        le=100,
        description="进入逐 CDS OSTIR 计算的分层 RBS 候选上限；从完整召回池按 low/medium/high 轮转抽取",
    )
    use_milvus: bool = Field(default=True, description="是否优先使用 Milvus expression_parts RAG 召回候选")
    selection_mode: Literal["metadata_only", "context_aware", "expression_aware"] = Field(
        default="metadata_only",
        description="metadata_only 保持原推荐；context_aware 仅重选 RBS；expression_aware 联合选择 promoter、逐 CDS RBS 和 terminator",
    )
    expression_targets: list[CassetteExpressionTarget] | None = Field(
        default=None,
        description="expression_aware 模式下每个 cassette 的 promoter/terminator 生物学意图",
    )
    regulatory_targets: list[RegulatoryTarget] | None = Field(
        default=None,
        description="context_aware 模式下每个 CDS 的绝对 TIR 区间或候选档位目标",
    )
    relative_constraints: list[RelativeExpressionConstraint] | None = Field(
        default=None,
        description="可选；同一设计内 CDS 之间的平衡或相对表达约束",
    )
    target_set_id: str = Field(default="", description="预注册目标集 ID")
    target_set_sha256: str = Field(default="", description="可选；预注册目标集内容哈希")
    strict_context_coverage: bool = Field(
        default=True,
        description="context_aware 模式下是否要求每个计划 CDS 都有且仅有一个目标",
    )


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _stable_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _read_manifest(path: Path) -> dict[str, Any]:
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


def _clean_text(value: Any) -> str:
    text = str(value or "").replace("\r", " ").replace("\n", " ")
    return re.sub(r"\s+", " ", text).strip()


def _normalize_part_id(part_id: str) -> str:
    return normalize_part_id(part_id)


def _normalize_strength(value: str) -> str:
    normalized = str(value or "unknown").strip().lower()
    return normalized if normalized in ALLOWED_STRENGTHS else "unknown"


def _normalize_regulation(value: str) -> str:
    normalized = str(value or "unknown").strip().lower()
    return normalized if normalized in ALLOWED_REGULATIONS else "unknown"


def _sequence_file(parts_dir: Path, part_id: str) -> Path:
    return part_sequence_path(parts_dir, part_id)


def _read_sequence_if_available(parts_dir: Path, part_id: str) -> tuple[bool, int | None]:
    path = _sequence_file(parts_dir, part_id)
    if not path.exists():
        return False, None
    text = path.read_text(encoding="utf-8").strip()
    sequence = "".join(
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()
    return bool(sequence), len(sequence) if sequence else None


def _write_sequence_if_available(parts_dir: Path, part_id: str, sequence: Any) -> None:
    normalized = "".join(str(sequence or "").split()).upper()
    if not normalized:
        return
    parts_dir.mkdir(parents=True, exist_ok=True)
    _sequence_file(parts_dir, part_id).write_text(normalized, encoding="utf-8")


def _cassette_map(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    assembly = manifest.get("expression_cassette_assembly")
    if not isinstance(assembly, dict):
        raise ValueError('design_manifest.json missing "expression_cassette_assembly"')

    cassettes = assembly.get("cassettes")
    if not isinstance(cassettes, list) or not cassettes:
        raise ValueError('design_manifest.json field "expression_cassette_assembly.cassettes" must be a non-empty list')

    normalized = []
    for cassette in cassettes:
        if not isinstance(cassette, dict):
            continue
        cassette_index = int(cassette.get("cassette_index") or 0)
        cds_ids = [str(cds_id).strip() for cds_id in _as_list(cassette.get("cds_ids")) if str(cds_id).strip()]
        if cassette_index > 0 and cds_ids:
            normalized.append({"cassette_index": cassette_index, "cds_ids": cds_ids})

    if not normalized:
        raise ValueError("expression_cassette_assembly has no valid cassettes")
    return normalized


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
        model_path = get_expression_embedding_model_path()
        sentence_transformer = _import_sentence_transformers()
        _EMBEDDING_MODEL = sentence_transformer(str(model_path), trust_remote_code=True)
    return _EMBEDDING_MODEL


def _embed_query(text: str) -> list[float]:
    model = _embedding_model()
    embedding = model.encode(text, normalize_embeddings=True)
    return embedding.tolist()


def _hit_entity(hit: Any) -> tuple[dict[str, Any], float]:
    if isinstance(hit, dict):
        entity = hit.get("entity") or hit
        distance = hit.get("distance", hit.get("score", 0.0))
        return dict(entity), float(distance or 0.0)
    entity = getattr(hit, "entity", {}) or {}
    distance = getattr(hit, "distance", getattr(hit, "score", 0.0))
    return dict(entity), float(distance or 0.0)


def _milvus_filter(
    role: str,
    host_labels: tuple[str, ...] | None = None,
) -> str:
    if not host_labels:
        return f'role == "{role}" and sequence_available == true'
    encoded_labels = json.dumps(list(host_labels), ensure_ascii=False)
    return (
        f'role == "{role}" and sequence_available == true '
        f"and ARRAY_CONTAINS_ANY(host, {encoded_labels})"
    )


def _build_role_query(
    *,
    role: str,
    host: str,
    preferred_strength: str,
    promoter_regulation: str,
    strategy: str,
) -> str:
    parts = [
        f"Recommend a {role} expression part for a synthetic biology expression cassette.",
        f"Host or chassis: {host}.",
        f"Expression strategy: {strategy}.",
        f"Desired expression strength: {preferred_strength}.",
    ]
    if role == "promoter":
        parts.append(f"Desired promoter regulation: {promoter_regulation}.")
    if role == "terminator":
        parts.append("Prefer robust common terminators with reliable Registry evidence; avoid weak, failed, or experimental short artificial terminators.")
    return " ".join(parts)


def _validate_v3_collection_schema(
    collection_name: str,
    description: dict[str, Any],
) -> None:
    """Fail closed when the formal V3 collection does not expose its required schema."""

    if collection_name != COLLECTION_NAME:
        return

    fields = {
        str(field.get("name") or ""): field
        for field in description.get("fields", [])
        if field.get("name")
    }
    missing = sorted(V3_REQUIRED_FIELDS - set(fields))
    if missing:
        raise ValueError(
            f"{COLLECTION_NAME} schema missing required fields: {missing}"
        )

    embedding_params = fields["embedding"].get("params") or {}
    embedding_dim = int(embedding_params.get("dim") or 0)
    if embedding_dim != V3_EMBEDDING_DIM:
        raise ValueError(
            f"{COLLECTION_NAME} embedding dimension mismatch: "
            f"expected {V3_EMBEDDING_DIM}, found {embedding_dim}"
        )


def _milvus_search_role(
    *,
    collection_name: str,
    role: str,
    query_embedding: list[float],
    search_k: int,
    host_labels: tuple[str, ...] | None = None,
) -> list[tuple[dict[str, Any], float]]:
    client = _milvus_client()
    if not client.has_collection(collection_name):
        raise ValueError(f"Milvus collection not found: {collection_name}")
    client.load_collection(collection_name)
    description = client.describe_collection(collection_name)
    _validate_v3_collection_schema(collection_name, description)
    schema_fields = {
        str(field.get("name") or "")
        for field in description.get("fields", [])
    }
    output_fields = [field for field in OUTPUT_FIELDS if field in schema_fields]
    result = client.search(
        collection_name=collection_name,
        data=[query_embedding],
        anns_field="embedding",
        limit=search_k,
        filter=_milvus_filter(role, host_labels),
        output_fields=output_fields,
        search_params={"metric_type": "COSINE"},
    )
    if not result:
        return []
    return [_hit_entity(hit) for hit in result[0]]


def _milvus_record_to_candidate(record: dict[str, Any], semantic_score: float) -> dict[str, Any]:
    description = _clean_text(record.get("description"))
    notes = _as_list(record.get("notes"))
    warnings = _as_list(record.get("warnings"))
    return {
        "part_id": _normalize_part_id(str(record.get("part_id") or "")),
        "role": str(record.get("role") or "").lower(),
        "source": _clean_text(record.get("source") or "milvus"),
        "strength": _normalize_strength(str(record.get("strength") or "unknown")),
        "strength_confidence": _clean_text(record.get("strength_confidence") or "unknown").lower(),
        "regulation": _normalize_regulation(str(record.get("regulation") or "unknown")),
        "regulation_confidence": _clean_text(record.get("regulation_confidence") or "unknown").lower(),
        "host": [str(item).strip() for item in _as_list(record.get("host")) if str(item).strip()],
        "host_confidence": _clean_text(record.get("host_confidence") or "unknown").lower(),
        "evidence": _clean_text(record.get("evidence")),
        "name": _clean_text(record.get("name")),
        "part_short_desc": description,
        "description": description,
        "sequence_type": _clean_text(record.get("sequence_type")).lower(),
        "sequence_type_confidence": _clean_text(record.get("sequence_type_confidence") or "unknown").lower(),
        "sequence": _clean_text(record.get("sequence")).upper(),
        "sequence_available": bool(record.get("sequence_available")),
        "sequence_length_bp": int(record.get("length_bp") or 0) or None,
        "registry_metadata": _as_dict(record.get("registry_metadata")),
        "warnings": [str(item) for item in warnings if str(item).strip()],
        "notes": "; ".join(str(item) for item in notes if str(item).strip()),
        "semantic_score": semantic_score,
        "candidate_source": "milvus",
        "activity_value": float(record.get("activity_value", -1.0) or -1.0),
        "activity_sd": float(record.get("activity_sd", -1.0) or -1.0),
        "activity_percentile": float(record.get("activity_percentile", -1.0) or -1.0),
        "activity_dataset": _clean_text(record.get("activity_dataset")),
        "activity_unit": _clean_text(record.get("activity_unit")),
        "activity_context": _clean_text(record.get("activity_context")),
        "activity_source": _clean_text(record.get("activity_source")),
        "evidence_grade": _clean_text(record.get("evidence_grade") or "C").upper(),
        "direction": _clean_text(record.get("direction") or "unknown").lower(),
        "regulator_required": _clean_text(record.get("regulator_required") or "unknown").lower(),
        "quantitative_metadata": _as_dict(record.get("quantitative_metadata")),
    }


def _retrieve_milvus_candidates(
    *,
    collection_name: str,
    search_k: int,
    host_resolution: ExpressionHostResolution,
    preferred_strength: str,
    promoter_regulation: str,
    strategy: str,
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for role in sorted(ALLOWED_ROLES):
        query_text = _build_role_query(
            role=role,
            host=host_resolution.query_name,
            preferred_strength=preferred_strength,
            promoter_regulation=promoter_regulation,
            strategy=strategy,
        )
        query_embedding = _embed_query(query_text)
        role_hits = _milvus_search_role(
            collection_name=collection_name,
            role=role,
            query_embedding=query_embedding,
            search_k=search_k,
            host_labels=host_resolution.milvus_host_labels,
        )
        if not role_hits:
            raise NoCompatibleExpressionPartsError(
                f"expression_parts_v3 returned no host-compatible {role} candidates "
                f"for {host_resolution.codon_transformer_name or host_resolution.query_name}"
            )
        for record, semantic_score in role_hits:
            try:
                candidates.append(_milvus_record_to_candidate(record, semantic_score))
            except Exception:
                continue
    return candidates


def _candidate_to_dict(candidate: CandidateExpressionPart) -> dict[str, Any]:
    payload = candidate.model_dump()
    payload["part_id"] = _normalize_part_id(payload["part_id"])
    payload["role"] = str(payload["role"]).lower()
    payload["strength"] = _normalize_strength(payload.get("strength", "unknown"))
    payload["regulation"] = _normalize_regulation(payload.get("regulation", "unknown"))
    return payload


def _merge_candidates(
    *,
    milvus_parts: list[dict[str, Any]],
    candidate_parts: list[CandidateExpressionPart] | None,
) -> list[dict[str, Any]]:
    merged: dict[tuple[str, str], dict[str, Any]] = {}

    for candidate in candidate_parts or []:
        part = _candidate_to_dict(candidate)
        key = (part["role"], part["part_id"].upper())
        existing = merged.get(key, {})
        merged[key] = {
            **existing,
            **{key_: value for key_, value in part.items() if value not in (None, "", [])},
            "source": part.get("source") or existing.get("source") or "user",
        }

    for part in milvus_parts:
        key = (str(part.get("role") or "").lower(), str(part.get("part_id") or "").upper())
        existing = merged.get(key, {})
        merged[key] = {
            **existing,
            **{key_: value for key_, value in part.items() if value not in (None, "", [])},
            "source": part.get("source") or existing.get("source") or "milvus",
        }

    return list(merged.values())


def _host_match_kind(
    part: dict[str, Any],
    host_resolution: ExpressionHostResolution,
) -> str:
    return candidate_host_match_kind(
        _as_list(part.get("host")),
        host_resolution,
    )


def _filter_host_compatible_candidates(
    candidates: list[dict[str, Any]],
    host_resolution: ExpressionHostResolution,
) -> tuple[list[dict[str, Any]], int]:
    compatible: list[dict[str, Any]] = []
    excluded = 0
    for candidate in candidates:
        match_kind = _host_match_kind(candidate, host_resolution)
        if match_kind == "none":
            excluded += 1
            continue
        compatible.append({**candidate, "_host_match_kind": match_kind})
    return compatible, excluded


def _confidence_score(value: Any) -> int:
    return CONFIDENCE_SCORES.get(str(value or "unknown").strip().lower(), -2)


def _registry_metadata(part: dict[str, Any]) -> dict[str, Any]:
    return _as_dict(part.get("registry_metadata"))


def _registry_parameter_text(part: dict[str, Any], name: str) -> str:
    parameters = _as_dict(_registry_metadata(part).get("parameters"))
    return str(parameters.get(name) or "").strip()


def _sequence_type_matches_role(part: dict[str, Any], role: str) -> bool:
    sequence_type = str(part.get("sequence_type") or "").strip().lower()
    if not sequence_type:
        return True
    if role == "terminator":
        return "terminator" in sequence_type
    return sequence_type == role


def _terminator_is_unsafe_default(part: dict[str, Any], sequence_length: int | None) -> bool:
    part_id = str(part.get("part_id") or "")
    if part_id in UNSAFE_TERMINATOR_IDS:
        return True
    metadata = _registry_metadata(part)
    result = str(metadata.get("part_results") or "").strip().lower()
    text = " ".join(
        str(value or "")
        for value in (
            part.get("description"),
            part.get("part_short_desc"),
            part.get("notes"),
            part.get("name"),
        )
    ).lower()
    if result == "fails" or "should not work" in text or "terninator" in text:
        return True
    if "artificial" in text and sequence_length is not None and sequence_length < 45:
        return True
    return False


def _write_candidate_sequences(parts_dir: Path, candidates: list[dict[str, Any]]) -> list[str]:
    written = []
    for candidate in candidates:
        sequence = candidate.get("sequence")
        if not sequence:
            continue
        part_id = _normalize_part_id(str(candidate.get("part_id") or ""))
        _write_sequence_if_available(parts_dir, part_id, sequence)
        written.append(part_id)
    return sorted(set(written))


def _score_part(
    part: dict[str, Any],
    *,
    host_resolution: ExpressionHostResolution,
    preferred_strength: str,
    promoter_regulation: str,
    parts_dir: Path,
) -> dict[str, Any]:
    role = str(part.get("role") or "").lower()
    part_id = _normalize_part_id(str(part.get("part_id") or ""))
    strength = _normalize_strength(str(part.get("strength") or "unknown"))
    regulation = _normalize_regulation(str(part.get("regulation") or "unknown"))
    sequence_available, sequence_length = _read_sequence_if_available(parts_dir, part_id)
    if not sequence_length and part.get("sequence_length_bp"):
        sequence_length = int(part.get("sequence_length_bp") or 0) or None

    score = 0
    reasons = []
    warnings = [str(item) for item in _as_list(part.get("warnings")) if str(item).strip()]

    if role in ALLOWED_ROLES:
        score += 25
        reasons.append("角色明确匹配")

    if _sequence_type_matches_role(part, role):
        score += 12
        if part.get("sequence_type"):
            reasons.append("sequence_type 与 role 匹配")
    else:
        score -= 80
        warnings.append(f"sequence_type {part.get('sequence_type')} 与 role {role} 不匹配")

    if sequence_available:
        score += 40
        reasons.append("本地已存在序列文件，可直接进入 select_parts")
    else:
        warnings.append("本地缺少序列文件，需先调用 get_parts_facts")

    host_match_kind = str(
        part.get("_host_match_kind")
        or _host_match_kind(part, host_resolution)
    )
    if host_match_kind == "exact":
        score += 20
        reasons.append(
            "宿主精确匹配 "
            f"{host_resolution.codon_transformer_name or host_resolution.query_name}"
        )
    elif host_match_kind == "lineage_transfer_inference":
        score += 8
        reasons.append("宿主通过显式 lineage-transfer inference 兼容")
        warnings.append(
            "宿主兼容性来自 lineage-transfer inference，不等同于目标菌株直接测量"
        )
    else:
        raise NoCompatibleExpressionPartsError(
            f"candidate {part_id} is not compatible with the upstream host"
        )

    score += _confidence_score(part.get("host_confidence"))

    preferred = _normalize_strength(preferred_strength)
    if role == "terminator":
        if strength == "high":
            score += 14
            reasons.append("terminator 强终止优先")
        elif strength == "medium":
            score += 6
            reasons.append("terminator 中等终止强度可作为备选")
        elif strength == "unknown":
            score -= 3
            warnings.append("terminator 强度未知")
    elif preferred == "any":
        score += 5
    elif strength == preferred:
        score += 15
        reasons.append(f"强度符合偏好 {preferred}")
    elif strength == "unknown":
        score -= 8
        warnings.append("强度未知")
    else:
        score -= 5
    score += _confidence_score(part.get("strength_confidence"))

    desired_regulation = _normalize_regulation(promoter_regulation)
    if role == "promoter":
        if desired_regulation == "any":
            score += 5
        elif regulation == desired_regulation:
            score += 15
            reasons.append(f"调控类型符合偏好 {desired_regulation}")
        elif regulation == "unknown":
            score -= 8
            warnings.append("promoter 调控类型未知")
        else:
            score -= 10
        score += _confidence_score(part.get("regulation_confidence"))

    source = str(part.get("source") or "").lower()
    if source in {"web_search", "igem", "igem registry", "registry", "milvus"}:
        score += 8
        reasons.append("候选来自外部检索或 Registry")
    if part.get("evidence"):
        score += 5
        reasons.append("包含证据来源")

    metadata = _registry_metadata(part)
    release_status = str(metadata.get("release_status") or "").lower()
    part_results = str(metadata.get("part_results") or "").lower()
    if "released hq" in release_status:
        score += 6
        reasons.append("Registry release_status 为 Released HQ")
    if part_results == "works":
        score += 8
        reasons.append("Registry part_results 为 Works")
    elif part_results == "fails":
        score -= 60
        warnings.append("Registry part_results 为 Fails")

    if role == "terminator":
        forward_efficiency = _registry_parameter_text(part, "forward_efficiency")
        if forward_efficiency and forward_efficiency.upper() != "NA":
            score += 5
            reasons.append("包含 terminator forward_efficiency 证据")
        if _terminator_is_unsafe_default(part, sequence_length):
            score -= 120
            warnings.append("默认推荐排除短人工、失败或不可靠 terminator")

    return {
        **part,
        "part_id": part_id,
        "role": role,
        "strength": strength,
        "regulation": regulation,
        "score": score,
        "sequence_available": sequence_available,
        "sequence_length_bp": sequence_length,
        "reasons": reasons,
        "warnings": warnings,
        "sequence_file": str(_sequence_file(parts_dir, part_id)),
    }


def _rank_by_role(scored_parts: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    ranked: dict[str, list[dict[str, Any]]] = {role: [] for role in sorted(ALLOWED_ROLES)}
    for part in scored_parts:
        role = str(part.get("role") or "").lower()
        if role in ranked:
            ranked[role].append(part)

    for role, items in ranked.items():
        ranked[role] = sorted(
            items,
            key=lambda item: (
                not bool(item.get("sequence_available")),
                -int(item.get("score") or 0),
                str(item.get("part_id") or ""),
            ),
        )
    return ranked


def _choose_recommended(ranked: dict[str, list[dict[str, Any]]]) -> dict[str, dict[str, Any]]:
    selected = {}
    for role in ALLOWED_ROLES:
        if not ranked.get(role):
            raise ValueError(f"no candidate parts available for role: {role}")
        selected[role] = ranked[role][0]
    return selected


def _selection_item(
    *,
    cassette_index: int,
    role: str,
    part: dict[str, Any],
    cds_id: str | None = None,
) -> dict[str, Any]:
    item = {
        "cassette_index": cassette_index,
        "role": role,
        "part_id": part["part_id"],
    }
    if cds_id:
        item["cds_id"] = cds_id
    for key in ("name", "part_short_desc", "sequence_type"):
        if part.get(key):
            item[key] = part[key]
    return item


def _build_selections(cassettes: list[dict[str, Any]], selected: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    selections = []
    for cassette in cassettes:
        cassette_index = int(cassette["cassette_index"])
        cds_ids = list(cassette["cds_ids"])
        selections.append(
            _selection_item(
                cassette_index=cassette_index,
                role="promoter",
                part=selected["promoter"],
            )
        )
        for cds_id in cds_ids:
            selections.append(
                _selection_item(
                    cassette_index=cassette_index,
                    role="rbs",
                    part=selected["rbs"],
                    cds_id=cds_id,
                )
            )
        selections.append(
            _selection_item(
                cassette_index=cassette_index,
                role="terminator",
                part=selected["terminator"],
            )
        )
    return selections


def _eligible_context_rbs(ranked: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    eligible = []
    for part in ranked.get("rbs", []):
        metadata = _registry_metadata(part)
        if not part.get("sequence_available"):
            continue
        if not _sequence_type_matches_role(part, "rbs"):
            continue
        if str(metadata.get("part_results") or "").strip().lower() == "fails":
            continue
        try:
            sequence = _sequence_file(Path(part["sequence_file"]).parent, part["part_id"]).read_text(
                encoding="utf-8"
            )
            normalized = "".join(sequence.split()).upper()
        except Exception:
            continue
        if normalized and not (set(normalized) - set("ACGT")):
            eligible.append(part)
    return eligible


def _shortlist_context_rbs(
    eligible: list[dict[str, Any]],
    limit: int,
) -> list[dict[str, Any]]:
    """Bound OSTIR work while preserving strength-tier coverage and ranking order."""

    if len(eligible) <= limit:
        return eligible
    grouped: dict[str, list[dict[str, Any]]] = {
        "medium": [], "low": [], "high": [], "unknown": [],
    }
    for part in eligible:
        tier = _normalize_strength(str(part.get("strength") or "unknown"))
        grouped[tier if tier in grouped else "unknown"].append(part)
    selected: list[dict[str, Any]] = []
    offsets = {tier: 0 for tier in grouped}
    tier_order = ("medium", "low", "high", "unknown")
    while len(selected) < limit:
        added = False
        for tier in tier_order:
            offset = offsets[tier]
            if offset < len(grouped[tier]):
                selected.append(grouped[tier][offset])
                offsets[tier] += 1
                added = True
                if len(selected) == limit:
                    break
        if not added:
            break
    return selected


def _build_context_selections(
    cassettes: list[dict[str, Any]],
    selected: dict[str, dict[str, Any]],
    selected_by_cds: dict[str, dict[str, Any]],
    rbs_by_id: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    selections = []
    for cassette in cassettes:
        cassette_index = int(cassette["cassette_index"])
        selections.append(
            _selection_item(cassette_index=cassette_index, role="promoter", part=selected["promoter"])
        )
        for cds_id in cassette["cds_ids"]:
            prediction = selected_by_cds[str(cds_id)]
            part_id = str(prediction["part_id"])
            selections.append(
                _selection_item(
                    cassette_index=cassette_index,
                    role="rbs",
                    part=rbs_by_id[part_id],
                    cds_id=str(cds_id),
                )
            )
        selections.append(
            _selection_item(cassette_index=cassette_index, role="terminator", part=selected["terminator"])
        )
    return selections


def _build_expression_selections(
    cassettes: list[dict[str, Any]],
    promoter_by_cassette: dict[int, dict[str, Any]],
    terminator_by_cassette: dict[int, dict[str, Any]],
    selected_by_cds: dict[str, dict[str, Any]],
    rbs_by_id: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    selections: list[dict[str, Any]] = []
    for cassette in cassettes:
        cassette_index = int(cassette["cassette_index"])
        selections.append(
            _selection_item(
                cassette_index=cassette_index,
                role="promoter",
                part=promoter_by_cassette[cassette_index],
            )
        )
        for cds_id in cassette["cds_ids"]:
            prediction = selected_by_cds[str(cds_id)]
            part_id = str(prediction["part_id"])
            selections.append(
                _selection_item(
                    cassette_index=cassette_index,
                    role="rbs",
                    part=rbs_by_id[part_id],
                    cds_id=str(cds_id),
                )
            )
        selections.append(
            _selection_item(
                cassette_index=cassette_index,
                role="terminator",
                part=terminator_by_cassette[cassette_index],
            )
        )
    return selections


def _json_safe_context(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_safe_context(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe_context(item) for item in value]
    return value


def _compact_part(part: dict[str, Any]) -> dict[str, Any]:
    payload = {
        "part_id": part.get("part_id"),
        "role": part.get("role"),
        "name": part.get("name", ""),
        "source": part.get("source", ""),
        "evidence": part.get("evidence", ""),
        "sequence_type": part.get("sequence_type", ""),
        "strength": part.get("strength"),
        "strength_confidence": part.get("strength_confidence", ""),
        "regulation": part.get("regulation"),
        "regulation_confidence": part.get("regulation_confidence", ""),
        "score": part.get("score"),
        "sequence_available": part.get("sequence_available"),
        "sequence_length_bp": part.get("sequence_length_bp"),
        "activity_value": part.get("activity_value"),
        "activity_percentile": part.get("activity_percentile"),
        "activity_dataset": part.get("activity_dataset", ""),
        "activity_unit": part.get("activity_unit", ""),
        "evidence_grade": part.get("evidence_grade", ""),
        "direction": part.get("direction", ""),
    }
    warnings = part.get("warnings") or []
    if warnings:
        payload["warnings"] = warnings
    return payload


def _compact_alternative(part: dict[str, Any]) -> dict[str, Any]:
    payload = {
        "part_id": part.get("part_id"),
        "score": part.get("score"),
        "sequence_available": part.get("sequence_available"),
        "sequence_length_bp": part.get("sequence_length_bp"),
    }
    name = part.get("name")
    if name:
        payload["name"] = name
    return payload


def _compact_selection(selection: dict[str, Any]) -> dict[str, Any]:
    payload = {
        "cassette_index": selection.get("cassette_index"),
        "role": selection.get("role"),
        "part_id": selection.get("part_id"),
    }
    cds_id = selection.get("cds_id")
    if cds_id:
        payload["cds_id"] = cds_id
    return payload


@tool(args_schema=RecommendExpressionPartsArgs)
def recommend_expression_parts(
    candidate_parts: list[CandidateExpressionPart] | None = None,
    host: str | None = None,
    expression_strength: str = "medium",
    promoter_regulation: str = "any",
    strategy: str = "balanced",
    alternatives_count: int = 0,
    detail: Literal["compact", "full"] = "compact",
    collection_name: str = COLLECTION_NAME,
    search_k: int = DEFAULT_SEARCH_K,
    context_rbs_limit: int = 24,
    use_milvus: bool = True,
    selection_mode: Literal["metadata_only", "context_aware", "expression_aware"] = "metadata_only",
    expression_targets: list[CassetteExpressionTarget] | None = None,
    regulatory_targets: list[RegulatoryTarget] | None = None,
    relative_constraints: list[RelativeExpressionConstraint] | None = None,
    target_set_id: str = "",
    target_set_sha256: str = "",
    strict_context_coverage: bool = True,
) -> str:
    """
    为当前表达盒计划推荐 promoter、RBS 和 terminator 组合。

    调用时机：CDS 和 expression_cassette_assembly 已确定，但 parts_selection 尚未写入。
    输入：可选宿主一致性断言、表达强度、调控偏好、Milvus collection 和可选候选 parts；可显式启用 context_aware。
    返回：默认精简推荐组合、select_parts_args、缺失序列和 warnings；detail=full 返回完整候选对象。
    限制：只读；用户确认后才调用 select_parts 写入 manifest。
    """

    tool_name = "recommend_expression_parts"
    monitor.report_start(tool_name, {"host": host, "strategy": strategy, "selection_mode": selection_mode})
    try:
        output_path = resolve_outputs_dir()
        manifest_path = design_manifest_file()
        monitor.report_running(tool_name, "正在读取 design_manifest 和候选表达元件...", progress=0.2)
        manifest = _read_manifest(manifest_path)
        resolved_parts_dir = resolve_parts_dir()
        host_resolution = resolve_expression_host(manifest, host)
        host_resolution_payload = host_resolution.to_dict()

        cassettes = _cassette_map(manifest)
        strategy_text = str(strategy or "balanced").strip().lower()
        strength = _normalize_strength(expression_strength)
        if strategy_text == "high_expression" and strength == "medium":
            strength = "high"
        elif strategy_text == "low_burden" and strength == "medium":
            strength = "low"

        retrieval = {
            "source": "none",
            "collection_name": collection_name,
            "search_k": search_k,
            "use_milvus": use_milvus,
            "fallback_used": False,
            "milvus_candidate_count": 0,
            "written_sequence_parts": [],
            "host_filter": list(host_resolution.milvus_host_labels),
            "host_incompatible_candidate_count": 0,
            "warnings": [],
        }
        milvus_candidates: list[dict[str, Any]] = []
        if use_milvus:
            try:
                monitor.report_running(
                    tool_name,
                    f"正在从 Milvus {collection_name} 按每类 Top-{search_k} 召回候选...",
                    progress=0.35,
                )
                milvus_candidates = _retrieve_milvus_candidates(
                    collection_name=collection_name,
                    search_k=max(search_k, 1),
                    host_resolution=host_resolution,
                    preferred_strength=strength,
                    promoter_regulation=promoter_regulation,
                    strategy=strategy_text,
                )
                retrieval["source"] = "milvus" if milvus_candidates else "none"
                retrieval["milvus_candidate_count"] = len(milvus_candidates)
                retrieval["written_sequence_parts"] = _write_candidate_sequences(resolved_parts_dir, milvus_candidates)
            except ExpressionHostContextError:
                raise
            except Exception as exc:
                raise RuntimeError(
                    f"formal Milvus retrieval from {collection_name} failed; "
                    "refusing fallback. Set use_milvus=false only for an explicitly "
                    f"candidate-only workflow. Cause: {type(exc).__name__}: {exc}"
                ) from exc
            if not milvus_candidates:
                raise ValueError(
                    f"formal Milvus retrieval from {collection_name} returned no candidates; "
                    "refusing fallback"
                )

        candidates = _merge_candidates(
            milvus_parts=milvus_candidates,
            candidate_parts=candidate_parts,
        )
        candidates, incompatible_count = _filter_host_compatible_candidates(
            candidates,
            host_resolution,
        )
        retrieval["host_incompatible_candidate_count"] = incompatible_count
        if not candidates:
            raise NoCompatibleExpressionPartsError(
                "no expression part candidates are compatible with cds_selection.host"
            )
        missing_host_compatible_roles = sorted(
            ALLOWED_ROLES
            - {
                str(candidate.get("role") or "").lower()
                for candidate in candidates
            }
        )
        if missing_host_compatible_roles:
            raise NoCompatibleExpressionPartsError(
                "no host-compatible candidates for roles "
                f"{missing_host_compatible_roles}; upstream host="
                f"{host_resolution.codon_transformer_name or host_resolution.query_name}"
            )

        if milvus_candidates:
            retrieval["source"] = "milvus"
        elif candidate_parts:
            retrieval["source"] = "candidate_parts"

        monitor.report_running(tool_name, "正在为 promoter/RBS/terminator 评分排序...", progress=0.55)
        scored = [
            _score_part(
                candidate,
                host_resolution=host_resolution,
                preferred_strength=strength,
                promoter_regulation=promoter_regulation,
                parts_dir=resolved_parts_dir,
            )
            for candidate in candidates
        ]
        ranked = _rank_by_role(scored)
        selected = _choose_recommended(ranked)
        context_analysis: dict[str, Any] | None = None
        expression_analysis: dict[str, Any] | None = None
        rbs_by_cds: list[dict[str, Any]] = []
        promoter_by_cassette = {
            int(cassette["cassette_index"]): selected["promoter"] for cassette in cassettes
        }
        terminator_by_cassette = {
            int(cassette["cassette_index"]): selected["terminator"] for cassette in cassettes
        }
        if selection_mode in {"context_aware", "expression_aware"}:
            targets = list(regulatory_targets or [])
            if not targets:
                raise ValueError(f"{selection_mode} selection requires regulatory_targets")
            eligible_rbs_unbounded = _eligible_context_rbs(ranked)
            if not eligible_rbs_unbounded:
                raise ValueError(f"{selection_mode} selection found no eligible RBS candidates with DNA sequences")
            eligible_rbs = _shortlist_context_rbs(eligible_rbs_unbounded, context_rbs_limit)
            expression_rows = list(expression_targets or [])
            if selection_mode == "expression_aware":
                if not expression_rows:
                    raise ValueError("expression_aware selection requires expression_targets")
                regulatory_result = optimize_expression_regulatory_parts(
                    cassettes=cassettes,
                    promoter_candidates=ranked["promoter"],
                    terminator_candidates=ranked["terminator"],
                    expression_targets=expression_rows,
                )
                promoter_by_cassette = regulatory_result["promoter_by_cassette"]
                terminator_by_cassette = regulatory_result["terminator_by_cassette"]
                expression_analysis = _json_safe_context(
                    {
                        "target_set_id": target_set_id,
                        "promoter_evaluations": regulatory_result["promoter_evaluations"],
                        "terminator_evaluations": regulatory_result["terminator_evaluations"],
                        "candidate_evaluations": regulatory_result["candidate_evaluations"],
                        "objective": regulatory_result["objective"],
                        "optimizer": regulatory_result["optimizer"],
                        "selected_promoter_by_cassette": [
                            {
                                "cassette_index": cassette_index,
                                **_compact_part(part),
                            }
                            for cassette_index, part in sorted(promoter_by_cassette.items())
                        ],
                        "selected_terminator_by_cassette": [
                            {
                                "cassette_index": cassette_index,
                                **_compact_part(part),
                            }
                            for cassette_index, part in sorted(terminator_by_cassette.items())
                        ],
                    }
                )
            monitor.report_running(
                tool_name,
                f"正在预测 {sum(len(item['cds_ids']) for item in cassettes)} 个 CDS 的候选 RBS 上下文 TIR...",
                progress=0.68,
            )
            raw_context = context_aware_reselection(
                manifest=manifest,
                output_dir=output_path,
                cassettes=cassettes,
                promoter=promoter_by_cassette,
                terminator=terminator_by_cassette,
                baseline_rbs=selected["rbs"],
                rbs_candidates=eligible_rbs,
                regulatory_targets=targets,
                relative_constraints=list(relative_constraints or []),
                strict_context_coverage=strict_context_coverage,
            )
            if selection_mode == "expression_aware":
                target_payload = [row.model_dump(mode="json") for row in targets]
                relative_payload = [row.model_dump(mode="json") for row in relative_constraints or []]
                expression_payload = [row.model_dump(mode="json") for row in expression_rows]
                actual_target_hash = _stable_sha256(
                    {
                        "expression_targets": expression_payload,
                        "regulatory_targets": target_payload,
                        "relative_constraints": relative_payload,
                    }
                )
            else:
                actual_target_hash = str(raw_context["target_set_sha256"])
            if target_set_sha256 and target_set_sha256 != actual_target_hash:
                raise ValueError(
                    f"target_set_sha256 mismatch: supplied={target_set_sha256}, computed={actual_target_hash}"
                )
            rbs_by_id = {str(part["part_id"]): part for part in eligible_rbs}
            if selection_mode == "expression_aware":
                selections = _build_expression_selections(
                    cassettes,
                    promoter_by_cassette,
                    terminator_by_cassette,
                    raw_context["selected_by_cds"],
                    rbs_by_id,
                )
                assert expression_analysis is not None
                expression_analysis["target_set_sha256"] = actual_target_hash
            else:
                selections = _build_context_selections(
                    cassettes,
                    selected,
                    raw_context["selected_by_cds"],
                    rbs_by_id,
                )
            rbs_by_cds = [
                {
                    "cassette_index": int(row["cassette_index"]),
                    "cds_id": cds_id,
                    "part_id": str(row["part_id"]),
                    "predicted_tir": round(float(row["tir"]), 6),
                    "candidate_percentile": round(float(row["percentile"]), 6),
                }
                for cds_id, row in sorted(
                    raw_context["selected_by_cds"].items(),
                    key=lambda item: (
                        int(item[1]["cassette_index"]),
                        int(item[1]["cds_order"]),
                        item[0],
                    ),
                )
            ]
            context_analysis = _json_safe_context(
                {
                    "target_set_id": target_set_id,
                    "target_set_sha256": actual_target_hash,
                    "candidate_count": len(eligible_rbs),
                    "candidate_universe_count": len(eligible_rbs_unbounded),
                    "candidate_shortlist_policy": "rank-preserving round-robin across medium/low/high/unknown tiers",
                    "candidate_shortlist_limit": context_rbs_limit,
                    "placement_count": len(raw_context["selected_by_cds"]),
                    "candidate_prediction_count": len(raw_context["candidate_predictions"]),
                    "candidate_predictions": raw_context["candidate_predictions"],
                    "target_evaluations": raw_context["target_evaluations"],
                    "relative_evaluations": raw_context["relative_evaluations"],
                    "objective": raw_context["objective"],
                    "optimizer": raw_context["optimizer"],
                    "selected_rbs_by_cds": rbs_by_cds,
                    "predictor": {
                        "name": "OSTIR",
                        "parameters": {
                            "threads": 1,
                            "decimal_places": 6,
                            "circular": False,
                        },
                    },
                }
            )
        else:
            selections = _build_selections(cassettes, selected)

        selected_part_lookup = {
            (str(part.get("role") or ""), str(part.get("part_id") or "")): part
            for part in scored
        }
        required_part_fetches = sorted({
            item["part_id"]
            for item in selections
            if not selected_part_lookup.get((item["role"], item["part_id"]), {}).get("sequence_available")
        })
        all_selected_have_sequences = not required_part_fetches

        alternatives = {
            role: items[: max(0, int(alternatives_count))]
            for role, items in ranked.items()
        }
        warnings = []
        warnings.extend(retrieval["warnings"])
        if required_part_fetches:
            warnings.append(
                "推荐组合中存在尚未保存序列的 part；请先调用 get_parts_facts，再调用 select_parts。"
            )
        selected_actual_parts = list(promoter_by_cassette.values()) + list(terminator_by_cassette.values())
        selected_actual_parts.append(selected["rbs"])
        for part in {str(row["part_id"]): row for row in selected_actual_parts}.values():
            warnings.extend(
                [f"{part['role']}:{part['part_id']} - {warning}" for warning in part.get("warnings", [])]
            )
        if selection_mode in {"context_aware", "expression_aware"} and any(
            target.analysis_tier == "exploratory" for target in regulatory_targets or []
        ):
            warnings.append("部分调控目标为 exploratory；不得将预测达标表述为实验表达或通路活性。")
        if selection_mode == "expression_aware" and any(
            row.promoter.analysis_tier == "exploratory" or row.terminator.analysis_tier == "exploratory"
            for row in expression_targets or []
        ):
            warnings.append("部分 promoter/terminator 意图为 exploratory；主分析仅统计 A/B 级约束。")

        first_cassette_index = min(promoter_by_cassette)
        primary_promoter = promoter_by_cassette[first_cassette_index]
        primary_terminator = terminator_by_cassette[first_cassette_index]
        promoter_by_cassette_payload = [
            {"cassette_index": index, **_compact_part(part)}
            for index, part in sorted(promoter_by_cassette.items())
        ]
        terminator_by_cassette_payload = [
            {"cassette_index": index, **_compact_part(part)}
            for index, part in sorted(terminator_by_cassette.items())
        ]

        full_result = {
            "ok": True,
            "manifest_path": str(manifest_path),
            "parts_dir": str(resolved_parts_dir),
            "strategy": strategy_text,
            "selection_mode": selection_mode,
            "host": host_resolution.codon_transformer_name,
            "host_resolution": host_resolution_payload,
            "expression_strength": strength,
            "promoter_regulation": _normalize_regulation(promoter_regulation),
            "retrieval": retrieval,
            "cassette_count": len(cassettes),
            "recommended": {
                "promoter": primary_promoter,
                "promoter_by_cassette": promoter_by_cassette_payload,
                "rbs": selected["rbs"],
                "rbs_by_cds": rbs_by_cds,
                "terminator": primary_terminator,
                "terminator_by_cassette": terminator_by_cassette_payload,
                "selections": selections,
                "all_selected_have_sequences": all_selected_have_sequences,
                "required_part_fetches": required_part_fetches,
                "rationale": [
                    "按 expression_cassette_assembly 为每个 cassette 分配 1 个 promoter 和 1 个 terminator。",
                    "按每个 cds_id 分配 1 个 RBS。",
                    "优先选择本地已有序列、宿主匹配、强度/调控偏好匹配的候选。",
                    "context_aware 仅逐 CDS 重选 RBS；expression_aware 同时按预注册意图选择 promoter 和 terminator。",
                ],
            },
            "context_analysis": context_analysis,
            "expression_analysis": expression_analysis,
            "alternatives": alternatives,
            "select_parts_args": {
                "selections": selections,
            },
            "next_actions": [
                "如果 required_part_fetches 非空，先对这些 part_id 调用 get_parts_facts。",
                "确认采用推荐组合后，调用 select_parts 写入 design_manifest.json 的 parts_selection 字段。",
                "parts_selection 写入后，再调用 assemble_expression_cassette 生成表达盒 GenBank 文件。",
            ],
            "warnings": warnings,
        }
        compact_result = {
            "ok": True,
            "strategy": strategy_text,
            "selection_mode": selection_mode,
            "host": host_resolution.codon_transformer_name,
            "host_resolution": host_resolution_payload,
            "expression_strength": strength,
            "promoter_regulation": _normalize_regulation(promoter_regulation),
            "retrieval": retrieval,
            "cassette_count": len(cassettes),
            "recommended": {
                "promoter": _compact_part(primary_promoter),
                "promoter_by_cassette": promoter_by_cassette_payload,
                "rbs": _compact_part(selected["rbs"]),
                "rbs_by_cds": rbs_by_cds,
                "terminator": _compact_part(primary_terminator),
                "terminator_by_cassette": terminator_by_cassette_payload,
                "all_selected_have_sequences": all_selected_have_sequences,
                "required_part_fetches": required_part_fetches,
            },
            "select_parts_args": {
                "selections": [_compact_selection(item) for item in selections],
            },
            "warnings": warnings,
        }
        if context_analysis is not None:
            compact_result["context_analysis"] = {
                key: context_analysis[key]
                for key in (
                    "target_set_id",
                    "target_set_sha256",
                    "candidate_count",
                    "placement_count",
                    "candidate_prediction_count",
                    "target_evaluations",
                    "relative_evaluations",
                    "objective",
                    "optimizer",
                    "selected_rbs_by_cds",
                    "predictor",
                )
            }
        if expression_analysis is not None:
            compact_result["expression_analysis"] = {
                key: expression_analysis[key]
                for key in (
                    "target_set_id",
                    "target_set_sha256",
                    "promoter_evaluations",
                    "terminator_evaluations",
                    "objective",
                    "optimizer",
                    "selected_promoter_by_cassette",
                    "selected_terminator_by_cassette",
                )
            }
        if any(alternatives.values()):
            compact_result["alternatives"] = {
                role: [_compact_alternative(item) for item in items]
                for role, items in alternatives.items()
                if items
            }

        result = full_result if detail == "full" else compact_result
        monitor.report_end(tool_name, {"cassette_count": len(cassettes), "required_part_fetches": required_part_fetches})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(recommend_expression_parts.invoke({}))
