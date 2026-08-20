"""End-to-end plasmid backbone recommendation pipeline."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Callable, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from src.plasmid_selection.config import (
    COPY_CLASS_SCORES,
    DEFAULT_CANDIDATE_COUNT,
    MARKER_SCORES,
    MAX_CANDIDATE_COUNT,
    MIN_CANDIDATE_COUNT,
    PLASMID_CANDIDATES_FILENAME,
    PLASMID_CANDIDATES_SCHEMA_VERSION,
    PLASMID_RECOMMENDATION_ALGORITHM_VERSION,
    SUPPORTED_PRIORITIES,
)
from src.plasmid_selection.get_plasmid_context import (
    load_plasmid_context,
    stable_json_hash,
)
from src.plasmid_selection.milvus_plasmids import fetch_plasmid_snapshot
from src.plasmid_selection.models import DownloadedSequence, PlasmidContext
from src.plasmid_selection.recommend_plasmid_templates import (
    normalize_marker_preferences,
    rank_templates,
    select_diverse_candidates,
)
from src.plasmid_selection.sequence_fetch import (
    download_and_validate_template,
    validate_genbank_bytes,
)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _parameters(
    *,
    count: int,
    priority: str,
    preferred_marker: str | None,
    excluded_markers: Sequence[str],
) -> dict[str, Any]:
    return {
        "requested_candidate_count": count,
        "priority": priority,
        "preferred_resistance": preferred_marker,
        "excluded_resistance": list(excluded_markers),
        "score_weights": {
            "copy_load_fit": 35,
            "assembly_readiness": 25,
            "source_evidence_completeness": 20,
            "marker_suitability": 10,
            "estimated_final_size": 10,
        },
        "copy_class_scores": COPY_CLASS_SCORES[priority],
        "default_marker_scores": MARKER_SCORES,
        "robust_aggregation": "minimum_pair_score_across_all_constructs",
        "diversity": {
            "primary_is_raw_best": True,
            "maximum_per_copy_class_before_relaxation": 2,
            "unique_tuple_before_relaxation": [
                "replicon_family",
                "marker",
                "assembly_policy",
            ],
        },
    }


def _request_fingerprint(
    context: PlasmidContext,
    *,
    snapshot_fingerprint: str,
    schema_fingerprint: str,
    parameters: Mapping[str, Any],
) -> str:
    return stable_json_hash(
        {
            "algorithm_version": PLASMID_RECOMMENDATION_ALGORITHM_VERSION,
            "context_input_fingerprint": context.input_fingerprint,
            "candidate_snapshot_fingerprint": snapshot_fingerprint,
            "milvus_schema_fingerprint": schema_fingerprint,
            "parameters": parameters,
        }
    )


def _candidate_set_fingerprint(candidates: Sequence[Mapping[str, Any]]) -> str:
    return stable_json_hash(
        [
            {
                key: value
                for key, value in candidate.items()
                if key not in {"system_recommended"}
            }
            for candidate in candidates
        ]
    )


def _resolve_local_candidate_path(
    project_output: Path,
    path_value: Any,
) -> Path | None:
    text = str(path_value or "").strip()
    if not text:
        return None
    path = Path(text)
    if not path.is_absolute():
        path = project_output / path
    path = path.resolve()
    try:
        path.relative_to(project_output.resolve())
    except ValueError:
        return None
    return path


def _cached_payload(
    artifact_path: Path,
    *,
    context: PlasmidContext,
    request_fingerprint: str,
) -> dict[str, Any] | None:
    if not artifact_path.is_file():
        return None
    try:
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    source = payload.get("source")
    candidates = payload.get("candidates")
    if (
        payload.get("schema_version") != PLASMID_CANDIDATES_SCHEMA_VERSION
        or payload.get("algorithm_version")
        != PLASMID_RECOMMENDATION_ALGORITHM_VERSION
        or payload.get("status") not in {"complete", "partial"}
        or not isinstance(source, Mapping)
        or source.get("request_fingerprint") != request_fingerprint
        or not isinstance(candidates, list)
        or not candidates
        or payload.get("candidate_count") != len(candidates)
    ):
        return None
    if (
        any(not isinstance(candidate, Mapping) for candidate in candidates)
        or payload.get("candidate_set_fingerprint")
        != _candidate_set_fingerprint(candidates)
        or [candidate.get("rank") for candidate in candidates]
        != list(range(1, len(candidates) + 1))
        or candidates[0].get("system_recommended") is not True
        or sum(
            candidate.get("system_recommended") is True
            for candidate in candidates
        )
        != 1
    ):
        return None
    for candidate in candidates:
        sequence_file = candidate.get("local_sequence_file")
        template = candidate.get("template")
        if not isinstance(sequence_file, Mapping) or not isinstance(template, Mapping):
            return None
        path = _resolve_local_candidate_path(
            context.project_output_path,
            sequence_file.get("path"),
        )
        if path is None or not path.is_file():
            return None
        try:
            content = path.read_bytes()
        except OSError:
            return None
        if _sha256_bytes(content) != sequence_file.get("file_sha256"):
            return None
        try:
            validated = validate_genbank_bytes(
                content,
                template,
                source_download_url=str(
                    sequence_file.get("source_download_url") or "cached"
                ),
            )
        except Exception:
            return None
        if (
            validated.sequence_content_sha256
            != sequence_file.get("sequence_content_sha256")
            or validated.canonical_sequence_sha256
            != sequence_file.get("canonical_sequence_sha256")
        ):
            return None
    return payload


def _strict_diverse_count(ranked: Sequence[dict[str, Any]], count: int) -> int:
    if not ranked:
        return 0
    selected = [ranked[0]]
    copy_counts = {str(ranked[0]["copy_number_class"]): 1}
    keys = {
        (
            str(ranked[0]["replicon_family"]).lower(),
            str(ranked[0]["marker"]).lower(),
            str(ranked[0]["assembly_policy"]).lower(),
        )
    }
    for item in ranked[1:]:
        copy_class = str(item["copy_number_class"])
        key = (
            str(item["replicon_family"]).lower(),
            str(item["marker"]).lower(),
            str(item["assembly_policy"]).lower(),
        )
        if copy_counts.get(copy_class, 0) >= 2 or key in keys:
            continue
        selected.append(item)
        copy_counts[copy_class] = copy_counts.get(copy_class, 0) + 1
        keys.add(key)
        if len(selected) >= count:
            break
    return len(selected)


def _candidate_payload(
    scored: Mapping[str, Any],
    downloaded: DownloadedSequence,
    *,
    rank: int,
    context: PlasmidContext,
) -> dict[str, Any]:
    payload = dict(scored)
    payload.update(
        {
            "candidate_id": f"candidate_{rank:03d}",
            "rank": rank,
            "system_recommended": rank == 1,
            "score_type": "interpretable_heuristic_not_probability",
            "local_sequence_file": {
                "path": (
                    f"plasmid_selection/candidates/candidate_{rank:03d}.gb"
                ),
                "format": "genbank",
                "length_bp": downloaded.length_bp,
                "file_sha256": downloaded.file_sha256,
                "sequence_content_sha256": downloaded.sequence_content_sha256,
                "canonical_sequence_sha256": downloaded.canonical_sequence_sha256,
                "source_download_url": downloaded.source_download_url,
            },
        }
    )
    return payload


def _commit_artifact(
    *,
    plasmid_dir: Path,
    artifact: Mapping[str, Any],
    selected: Sequence[dict[str, Any]],
    downloaded_by_id: Mapping[str, DownloadedSequence],
) -> None:
    project_output = plasmid_dir.parent
    project_output.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=".plasmid_selection.stage.", dir=project_output)
    )
    backup = Path(
        tempfile.mkdtemp(prefix=".plasmid_selection.backup.", dir=project_output)
    )
    target_candidates = plasmid_dir / "candidates"
    target_artifact = plasmid_dir / PLASMID_CANDIDATES_FILENAME
    backup_candidates = backup / "candidates"
    backup_artifact = backup / PLASMID_CANDIDATES_FILENAME
    moved_old_candidates = False
    moved_old_artifact = False
    installed_candidates = False
    installed_artifact = False
    try:
        stage_candidates = staging / "candidates"
        stage_candidates.mkdir(parents=True)
        for rank, item in enumerate(selected, start=1):
            downloaded = downloaded_by_id[str(item["plasmid_id"])]
            path = stage_candidates / f"candidate_{rank:03d}.gb"
            path.write_bytes(downloaded.content)
            if _sha256_bytes(path.read_bytes()) != downloaded.file_sha256:
                raise OSError(f"failed to verify staged candidate file: {path}")
        stage_artifact = staging / PLASMID_CANDIDATES_FILENAME
        _write_json(stage_artifact, artifact)

        plasmid_dir.mkdir(parents=True, exist_ok=True)
        if target_candidates.exists():
            os.replace(target_candidates, backup_candidates)
            moved_old_candidates = True
        if target_artifact.exists():
            os.replace(target_artifact, backup_artifact)
            moved_old_artifact = True
        os.replace(stage_candidates, target_candidates)
        installed_candidates = True
        os.replace(stage_artifact, target_artifact)
        installed_artifact = True
    except Exception:
        if installed_artifact and target_artifact.exists():
            target_artifact.unlink()
        if installed_candidates and target_candidates.exists():
            shutil.rmtree(target_candidates)
        if moved_old_artifact and backup_artifact.exists():
            os.replace(backup_artifact, target_artifact)
        if moved_old_candidates and backup_candidates.exists():
            os.replace(backup_candidates, target_candidates)
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)
        shutil.rmtree(backup, ignore_errors=True)


def plasmid_recommendation_summary(
    payload: Mapping[str, Any],
    artifact_path: Path,
    *,
    reused_existing: bool,
) -> dict[str, Any]:
    candidates = payload.get("candidates")
    candidates = candidates if isinstance(candidates, list) else []
    return {
        "ok": payload.get("status") == "complete",
        "status": payload.get("status"),
        "target_compound_id": payload.get("target_compound_id"),
        "candidate_count": len(candidates),
        "requested_candidate_count": (
            payload.get("ranking_parameters", {}).get("requested_candidate_count")
            if isinstance(payload.get("ranking_parameters"), Mapping)
            else None
        ),
        "primary_recommendation": (
            {
                "rank": candidates[0].get("rank"),
                "plasmid_id": candidates[0].get("plasmid_id"),
                "name": candidates[0].get("name"),
                "score": candidates[0].get("robust_score"),
                "copy_number_class": candidates[0].get("copy_number_class"),
                "marker": candidates[0].get("marker"),
            }
            if candidates
            else None
        ),
        "candidates": [
            {
                "rank": item.get("rank"),
                "name": item.get("name"),
                "score": item.get("robust_score"),
                "copy": item.get("copy_number_class"),
                "marker": item.get("marker"),
                "assembly": item.get("assembly_policy"),
            }
            for item in candidates
        ],
        "rejected_count": payload.get("rejected_count", 0),
        "warnings": payload.get("warnings", []),
        "output_path": str(artifact_path.resolve()),
        "reused_existing": reused_existing,
        "manifest_modified": False,
    }


def run_plasmid_recommendation(
    config: Any,
    *,
    client: Any | None = None,
    downloader: Callable[[Mapping[str, Any]], DownloadedSequence] = (
        download_and_validate_template
    ),
) -> dict[str, Any]:
    """Recommend a diverse set of one-backbone-for-all-designs candidates."""

    raw_count = getattr(config, "n_candidates", DEFAULT_CANDIDATE_COUNT)
    count = DEFAULT_CANDIDATE_COUNT if raw_count is None else int(raw_count)
    if not MIN_CANDIDATE_COUNT <= count <= MAX_CANDIDATE_COUNT:
        raise ValueError(
            f"--n-candidates 必须在 {MIN_CANDIDATE_COUNT}-{MAX_CANDIDATE_COUNT} 之间"
        )
    priority = str(getattr(config, "priority", "stability") or "stability")
    if priority not in SUPPORTED_PRIORITIES:
        raise ValueError(f"不支持的 --priority：{priority}")
    preferred_marker, excluded_markers = normalize_marker_preferences(
        getattr(config, "preferred_resistance", None),
        getattr(config, "exclude_resistance", None),
    )
    context = load_plasmid_context(config)
    snapshot = fetch_plasmid_snapshot(client=client)
    parameters = _parameters(
        count=count,
        priority=priority,
        preferred_marker=preferred_marker,
        excluded_markers=excluded_markers,
    )
    request_fingerprint = _request_fingerprint(
        context,
        snapshot_fingerprint=snapshot.candidate_fingerprint,
        schema_fingerprint=snapshot.schema_fingerprint,
        parameters=parameters,
    )
    plasmid_dir = context.project_output_path / "plasmid_selection"
    artifact_path = plasmid_dir / PLASMID_CANDIDATES_FILENAME
    cached = _cached_payload(
        artifact_path,
        context=context,
        request_fingerprint=request_fingerprint,
    )
    if cached is not None:
        return plasmid_recommendation_summary(
            cached, artifact_path, reused_existing=True
        )

    ranked, rejected = rank_templates(
        snapshot.candidates,
        context,
        priority=priority,
        preferred_marker=preferred_marker,
        excluded_markers=excluded_markers,
    )
    validated: list[dict[str, Any]] = []
    downloaded_by_id: dict[str, DownloadedSequence] = {}
    for item in ranked:
        template = item["template"]
        try:
            downloaded = downloader(template)
        except Exception as exc:
            rejected.append(
                {
                    "plasmid_id": str(item.get("plasmid_id") or ""),
                    "name": str(item.get("name") or ""),
                    "stage": "sequence_download",
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
            continue
        identity = str(item["plasmid_id"])
        downloaded_by_id[identity] = downloaded
        validated.append(item)
        if _strict_diverse_count(validated, count) >= count:
            break

    selected = select_diverse_candidates(validated, count)
    candidates = [
        _candidate_payload(
            item,
            downloaded_by_id[str(item["plasmid_id"])],
            rank=rank,
            context=context,
        )
        for rank, item in enumerate(selected, start=1)
    ]
    status = (
        "complete"
        if len(candidates) == count
        else "partial"
        if candidates
        else "failed"
    )
    warnings = [
        "Plasmid insert capacity is estimated because plasmid_templates_v2 has no experimentally verified capacity field.",
        "Scores are interpretable heuristics and are not experimental success probabilities.",
        "One selected backbone will later be combined separately with every selected expression construct.",
    ]
    if status == "partial":
        warnings.append(
            f"Only {len(candidates)} of {count} requested downloadable candidates were available."
        )
    if status == "failed":
        warnings.append("No candidate passed both metadata and sequence validation.")
    candidate_set_fingerprint = _candidate_set_fingerprint(candidates)
    payload: dict[str, Any] = {
        "schema_version": PLASMID_CANDIDATES_SCHEMA_VERSION,
        "algorithm_version": PLASMID_RECOMMENDATION_ALGORITHM_VERSION,
        "status": status,
        "score_type": "interpretable_heuristic_not_probability",
        "target_compound_id": context.target_compound_id,
        "generated_at": _utc_now(),
        "source": {
            "manifest_path": str(context.manifest_path),
            "manifest_revision": context.manifest_revision,
            "context_input_fingerprint": context.input_fingerprint,
            "parts_selection_fingerprint": context.parts_selection_fingerprint,
            "assembled_constructs_fingerprint": context.assembled_constructs_fingerprint,
            "milvus_collection": snapshot.collection_name,
            "milvus_row_count": snapshot.collection_row_count,
            "milvus_schema_fingerprint": snapshot.schema_fingerprint,
            "candidate_snapshot_fingerprint": snapshot.candidate_fingerprint,
            "request_fingerprint": request_fingerprint,
        },
        "requirements": {
            "host": {"key": context.host_key, "name": context.host_name},
            "selected_parts_design_ids": list(context.selected_design_ids),
            "construct_count": len(context.constructs),
            "insert_length_range_bp": {
                "minimum": context.minimum_insert_length_bp,
                "maximum": context.maximum_insert_length_bp,
            },
            "maximum_cassette_count": context.maximum_cassette_count,
            "backbone_policy": "one_backbone_applicable_to_all_selected_constructs",
        },
        "ranking_parameters": parameters,
        "candidate_count": len(candidates),
        "candidate_set_fingerprint": candidate_set_fingerprint,
        "candidates": candidates,
        "rejected_count": len(rejected),
        "rejected": sorted(
            rejected,
            key=lambda item: (
                str(item.get("stage") or ""),
                str(item.get("plasmid_id") or ""),
            ),
        ),
        "warnings": warnings,
    }
    _commit_artifact(
        plasmid_dir=plasmid_dir,
        artifact=payload,
        selected=selected,
        downloaded_by_id=downloaded_by_id,
    )
    return plasmid_recommendation_summary(
        payload, artifact_path, reused_existing=False
    )


__all__ = [
    "PLASMID_CANDIDATES_FILENAME",
    "plasmid_recommendation_summary",
    "run_plasmid_recommendation",
]
