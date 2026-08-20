"""Validate and commit one recommended plasmid backbone."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.plasmid_selection.config import (
    PLASMID_CANDIDATES_FILENAME,
    PLASMID_CANDIDATES_SCHEMA_VERSION,
    PLASMID_RECOMMENDATION_ALGORITHM_VERSION,
    PLASMID_SELECTION_SCHEMA_VERSION,
)
from src.plasmid_selection.get_plasmid_context import (
    load_plasmid_context,
    stable_json_hash,
)
from src.plasmid_selection.sequence_fetch import validate_genbank_bytes
from src.write_manifest.store import read_design_manifest, update_design_manifest


PLASMID_SELECTION_DOWNSTREAM_SECTIONS = (
    "final_assembly_plan",
    "final_assembly",
)


def _artifact_path(config: Any) -> Path:
    return (
        Path(config.project_output_path).expanduser().resolve()
        / "plasmid_selection"
        / PLASMID_CANDIDATES_FILENAME
    )


def _read_artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(
            "未找到质粒候选，请先运行：plasmid --recommend -i <输入文件>"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"质粒候选文件不是有效 JSON：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"质粒候选文件根节点必须是 JSON 对象：{path}")
    return payload


def _positive_rank(config: Any) -> int:
    raw = getattr(config, "plasmid", None)
    if isinstance(raw, bool) or raw is None:
        raise ValueError("未指定质粒候选，请使用 write --plasmid N")
    try:
        rank = int(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError("--plasmid 必须是正整数") from exc
    if rank < 1:
        raise ValueError("--plasmid 必须是正整数")
    return rank


def _candidate_set_fingerprint(candidates: list[Any]) -> str:
    return stable_json_hash(
        [
            {
                key: value
                for key, value in candidate.items()
                if key not in {"system_recommended"}
            }
            for candidate in candidates
            if isinstance(candidate, Mapping)
        ]
    )


def _resolve_project_file(project_output: Path, value: Any) -> Path:
    text = str(value or "").strip()
    if not text:
        raise ValueError("候选质粒缺少本地 GenBank 路径")
    path = Path(text)
    if not path.is_absolute():
        path = project_output / path
    path = path.resolve()
    try:
        path.relative_to(project_output.resolve())
    except ValueError as exc:
        raise ValueError(f"候选质粒文件不在当前项目输出目录内：{path}") from exc
    return path


def _validate_artifact(
    artifact: Mapping[str, Any],
    *,
    context: Any,
) -> list[dict[str, Any]]:
    if artifact.get("schema_version") != PLASMID_CANDIDATES_SCHEMA_VERSION:
        raise ValueError("质粒候选 schema 已过期，请重新运行 plasmid --recommend")
    if (
        artifact.get("algorithm_version")
        != PLASMID_RECOMMENDATION_ALGORITHM_VERSION
    ):
        raise ValueError("质粒推荐算法已变化，请重新运行 plasmid --recommend")
    if artifact.get("status") not in {"complete", "partial"}:
        raise ValueError("质粒推荐结果中没有可写入候选")
    if artifact.get("target_compound_id") != context.target_compound_id:
        raise ValueError("质粒候选与当前目标化合物不一致")
    source = artifact.get("source")
    if not isinstance(source, Mapping):
        raise ValueError("质粒候选缺少 source")
    expected = {
        "context_input_fingerprint": context.input_fingerprint,
        "parts_selection_fingerprint": context.parts_selection_fingerprint,
        "assembled_constructs_fingerprint": context.assembled_constructs_fingerprint,
    }
    for key, value in expected.items():
        if source.get(key) != value:
            raise ValueError(
                "质粒候选与当前表达构建不一致，请重新运行 plasmid --recommend"
            )
    raw_candidates = artifact.get("candidates")
    if not isinstance(raw_candidates, list) or not raw_candidates:
        raise ValueError("质粒候选文件中没有 candidates")
    candidates = [dict(item) for item in raw_candidates if isinstance(item, Mapping)]
    if len(candidates) != len(raw_candidates):
        raise ValueError("质粒候选包含无效条目")
    ranks = [item.get("rank") for item in candidates]
    if ranks != list(range(1, len(candidates) + 1)):
        raise ValueError("质粒候选排名不连续")
    primary_count = sum(item.get("system_recommended") is True for item in candidates)
    if primary_count != 1 or candidates[0].get("system_recommended") is not True:
        raise ValueError("质粒候选必须且只能有一个首选方案")
    if artifact.get("candidate_set_fingerprint") != _candidate_set_fingerprint(candidates):
        raise ValueError("质粒候选内容指纹无效")
    return candidates


def _validate_candidate_file(
    candidate: Mapping[str, Any],
    *,
    project_output: Path,
) -> tuple[Path, bytes, dict[str, Any]]:
    sequence_file = candidate.get("local_sequence_file")
    template = candidate.get("template")
    if not isinstance(sequence_file, Mapping) or not isinstance(template, Mapping):
        raise ValueError("候选质粒缺少序列或模板元数据")
    path = _resolve_project_file(project_output, sequence_file.get("path"))
    if not path.is_file():
        raise FileNotFoundError(
            f"候选质粒 GenBank 文件缺失：{path}；请重新运行 plasmid --recommend"
        )
    content = path.read_bytes()
    actual_file_hash = hashlib.sha256(content).hexdigest()
    if actual_file_hash != sequence_file.get("file_sha256"):
        raise ValueError(
            f"候选质粒文件哈希无效：{path}；请重新运行 plasmid --recommend"
        )
    validated = validate_genbank_bytes(
        content,
        template,
        source_download_url=str(sequence_file.get("source_download_url") or "local"),
    )
    if (
        validated.sequence_content_sha256
        != sequence_file.get("sequence_content_sha256")
        or validated.canonical_sequence_sha256
        != sequence_file.get("canonical_sequence_sha256")
    ):
        raise ValueError("候选质粒序列审计值与候选文件不一致")
    return path, content, dict(sequence_file)


def _install_bytes_atomic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _selected_file_is_intact(
    path: Path,
    expected_file_hash: str,
    template: Mapping[str, Any],
) -> bool:
    if not path.is_file():
        return False
    try:
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != expected_file_hash:
            return False
        validate_genbank_bytes(content, template, source_download_url="selected")
    except Exception:
        return False
    return True


def _selection_payload(
    candidate: Mapping[str, Any],
    *,
    context: Any,
    artifact_path: Path,
    artifact: Mapping[str, Any],
    sequence_file: Mapping[str, Any],
) -> dict[str, Any]:
    template = candidate["template"]
    selected_relative_path = "plasmid_selection/selected_backbone.gb"
    vector = {
        "plasmid_id": template.get("plasmid_id"),
        "name": template.get("name"),
        "description": template.get("description"),
        "vector_type": template.get("vector_type"),
        "length_bp": template.get("length_bp"),
        "topology": template.get("topology"),
        "replicon_family": template.get("replicon_family"),
        "copy_number": template.get("copy_number"),
        "copy_number_class": template.get("copy_number_class"),
        "marker": candidate.get("marker"),
        "bacterial_resistance": template.get("bacterial_resistance"),
        "assembly_policy": template.get("assembly_policy"),
        "cargo_type": template.get("cargo_type"),
        "requires_cargo_replacement": template.get("requires_cargo_replacement"),
        "host_compatibility": template.get("host_compatibility", []),
        "replication_dependencies": template.get("replication_dependencies", []),
        "origins": template.get("origins", []),
        "selection_markers": template.get("selection_markers", []),
        "insertion_regions": template.get("insertion_regions", []),
        "protected_features": template.get("protected_features", []),
        "selected_sequence_file": {
            "path": selected_relative_path,
            "format": "genbank",
            "length_bp": sequence_file.get("length_bp"),
            "file_sha256": sequence_file.get("file_sha256"),
            "sequence_content_sha256": sequence_file.get(
                "sequence_content_sha256"
            ),
            "canonical_sequence_sha256": sequence_file.get(
                "canonical_sequence_sha256"
            ),
        },
        "source": {
            "source": template.get("source"),
            "source_url": template.get("source_url"),
            "sequence_url": template.get("sequence_url"),
            "source_record_id": template.get("source_record_id"),
            "source_record_version": template.get("source_record_version"),
            "source_provenance": template.get("source_provenance", {}),
            "evidence_refs": template.get("evidence_refs", []),
            "audit": {
                "schema_version": template.get("schema_version"),
                "status": template.get("audit_status"),
                "version": template.get("audit_version"),
                "passed": template.get("audit_passed"),
                "checks": template.get("audit_checks", {}),
            },
        },
    }
    recommendation = {
        "candidate_id": candidate.get("candidate_id"),
        "rank": candidate.get("rank"),
        "system_recommended": candidate.get("system_recommended"),
        "score_type": candidate.get("score_type"),
        "robust_score": candidate.get("robust_score"),
        "pair_score_median": candidate.get("pair_score_median"),
        "score_breakdown": candidate.get("score_breakdown"),
        "copy_load_fit_details": candidate.get("copy_load_fit_details"),
        "expression_burden_score_range": candidate.get(
            "expression_burden_score_range"
        ),
        "expression_burden_levels": candidate.get(
            "expression_burden_levels", []
        ),
        "expression_burden_model_version": candidate.get(
            "expression_burden_model_version"
        ),
        "confidence": candidate.get("confidence"),
        "rationales": candidate.get("rationales", []),
        "warnings": candidate.get("warnings", []),
    }
    base = {
        "schema_version": PLASMID_SELECTION_SCHEMA_VERSION,
        "status": "selected",
        "selection_status": "user_selected",
        "source": "write_plasmid_selection",
        "backbone_policy": "one_backbone_for_all_selected_constructs",
        "covered_parts_design_ids": list(context.selected_design_ids),
        "covered_design_count": len(context.selected_design_ids),
        "vector": vector,
        "recommendation": recommendation,
        "provenance": {
            "candidate_artifact": str(
                artifact_path.relative_to(context.project_output_path).as_posix()
            ),
            "candidate_artifact_schema_version": artifact.get("schema_version"),
            "candidate_algorithm_version": artifact.get("algorithm_version"),
            "candidate_set_fingerprint": artifact.get("candidate_set_fingerprint"),
            "request_fingerprint": artifact.get("source", {}).get(
                "request_fingerprint"
            ),
            "candidate_snapshot_fingerprint": artifact.get("source", {}).get(
                "candidate_snapshot_fingerprint"
            ),
            "context_input_fingerprint": context.input_fingerprint,
            "parts_selection_fingerprint": context.parts_selection_fingerprint,
            "assembled_constructs_fingerprint": context.assembled_constructs_fingerprint,
            "expression_burden_model_version": context.constructs[
                0
            ].burden.model_version,
        },
        "warnings": list(candidate.get("warnings") or []),
    }
    base["selection_fingerprint"] = stable_json_hash(base)
    return base


def write_plasmid_selection(config: Any) -> dict[str, Any]:
    """Write one ranked candidate and install its exact GenBank backbone."""

    rank = _positive_rank(config)
    context = load_plasmid_context(config)
    artifact_path = _artifact_path(config)
    artifact = _read_artifact(artifact_path)
    candidates = _validate_artifact(artifact, context=context)
    if rank > len(candidates):
        raise ValueError(
            f"质粒候选 {rank} 不存在；可用编号：1-{len(candidates)}"
        )
    candidate = candidates[rank - 1]
    _, candidate_content, sequence_file = _validate_candidate_file(
        candidate,
        project_output=context.project_output_path,
    )
    payload = _selection_payload(
        candidate,
        context=context,
        artifact_path=artifact_path,
        artifact=artifact,
        sequence_file=sequence_file,
    )
    destination = (
        context.project_output_path
        / "plasmid_selection"
        / "selected_backbone.gb"
    )
    manifest = read_design_manifest(context.manifest_path)
    current = manifest.get("plasmid_selection")
    unchanged = (
        isinstance(current, Mapping)
        and current.get("schema_version") == PLASMID_SELECTION_SCHEMA_VERSION
        and current.get("selection_fingerprint")
        == payload["selection_fingerprint"]
    )
    intact = _selected_file_is_intact(
        destination,
        str(sequence_file.get("file_sha256") or ""),
        candidate["template"],
    )
    repaired = False
    if unchanged:
        if not intact:
            _install_bytes_atomic(destination, candidate_content)
            repaired = True
        updated_manifest = manifest
    else:
        old_content = destination.read_bytes() if destination.is_file() else None
        _install_bytes_atomic(destination, candidate_content)
        try:
            updated_manifest = update_design_manifest(
                context.manifest_path,
                target_compound_id=context.target_compound_id,
                sections={"plasmid_selection": payload},
                discard_sections=PLASMID_SELECTION_DOWNSTREAM_SECTIONS,
                expected_revision=context.manifest_revision,
            )
        except Exception:
            if old_content is None:
                destination.unlink(missing_ok=True)
            else:
                _install_bytes_atomic(destination, old_content)
            raise

    return {
        "ok": True,
        "status": "selected",
        "target_compound_id": context.target_compound_id,
        "candidate_rank": rank,
        "plasmid_id": candidate.get("plasmid_id"),
        "name": candidate.get("name"),
        "score": candidate.get("robust_score"),
        "covered_design_count": len(context.selected_design_ids),
        "selected_backbone": str(destination.resolve()),
        "manifest_path": str(context.manifest_path.resolve()),
        "manifest_revision": updated_manifest["revision"],
        "manifest_modified": not unchanged,
        "selected_file_repaired": repaired,
    }


def run_write_plasmid_selection(config: Any) -> dict[str, Any]:
    result = write_plasmid_selection(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "PLASMID_SELECTION_DOWNSTREAM_SECTIONS",
    "run_write_plasmid_selection",
    "write_plasmid_selection",
]
