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
from filelock import FileLock
from src.write_manifest.store import manifest_update_lock, update_design_manifest

from src.plasmid_selection.config import (
    BURDEN_COPY_FIT_SCORES,
    DEFAULT_CANDIDATE_COUNT,
    MARKER_SCORES,
    MAX_CANDIDATE_COUNT,
    MIN_CANDIDATE_COUNT,
    PLASMID_CANDIDATES_FILENAME,
    PLASMID_CANDIDATES_SCHEMA_VERSION,
    PLASMID_RECOMMENDATION_ALGORITHM_VERSION,
    PRIORITY_COPY_ADJUSTMENTS,
    SUPPORTED_PRIORITIES,
)
from src.plasmid_selection.get_plasmid_context import (
    load_plasmid_context,
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
    burden_model_version: str,
) -> dict[str, Any]:
    return {
        "requested_candidate_count": count,
        "priority": priority,
        "preferred_resistance": preferred_marker,
        "excluded_resistance": list(excluded_markers),
        "expression_burden_model_version": burden_model_version,
        "score_weights": {
            "copy_load_fit": 35,
            "assembly_readiness": 25,
            "source_evidence_completeness": 20,
            "marker_suitability": 10,
            "estimated_final_size": 10,
        },
        "dynamic_copy_fit_scores": BURDEN_COPY_FIT_SCORES,
        "priority_copy_adjustments": PRIORITY_COPY_ADJUSTMENTS,
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


def _raw_sequence(plasmid_dir: Path, template: Mapping[str, Any], downloader):
    """Reuse only downloaded source bytes, never a ranked recommendation."""
    locator = [template.get(k) for k in ("source", "source_record_id", "source_record_version", "sequence_url", "source_url")]
    key = hashlib.sha256(json.dumps(locator, ensure_ascii=False).encode("utf-8")).hexdigest()
    cache_path = plasmid_dir / "raw_sequences" / (key + ".gb")
    if cache_path.is_file():
        try:
            return validate_genbank_bytes(cache_path.read_bytes(), template, source_download_url="raw_cache")
        except (OSError, ValueError, RuntimeError):
            pass
    downloaded = downloader(template)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=cache_path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(downloaded.content)
            handle.flush()
            os.fsync(handle.fileno())
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, cache_path)
    finally:
        temporary.unlink(missing_ok=True)
    return downloaded


def _seed_old_raw_sequences(plasmid_dir: Path, artifact_path: Path) -> None:
    """Recover source bytes from old reports without reusing their ranking."""
    try:
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        candidates = payload.get("candidates", [])
        for candidate in candidates:
            template = candidate.get("template")
            info = candidate.get("local_sequence_file")
            if not isinstance(template, Mapping) or not isinstance(info, Mapping):
                continue
            path = Path(str(info.get("path") or ""))
            if not path.is_absolute():
                path = plasmid_dir.parent / path
            path = path.resolve()
            path.relative_to((plasmid_dir / "candidates").resolve())
            content = path.read_bytes()
            downloaded = validate_genbank_bytes(
                content, template, source_download_url=str(info.get("source_download_url") or "old_raw_file")
            )
            _raw_sequence(plasmid_dir, template, lambda _: downloaded)
    except (OSError, ValueError, RuntimeError, TypeError, AttributeError):
        # A malformed historical artifact is never a reason to reuse it.
        return


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
                "burden_range": candidates[0].get(
                    "expression_burden_score_range"
                ),
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
                "burden_levels": item.get("expression_burden_levels"),
            }
            for item in candidates
        ],
        "rejected_count": payload.get("rejected_count", 0),
        "warnings": payload.get("warnings", []),
        "output_path": str(artifact_path.resolve()),
        "reused_existing": reused_existing,
        "manifest_modified": True,
    }


def _run_plasmid_recommendation(
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
    parameters = _parameters(
        count=count,
        priority=priority,
        preferred_marker=preferred_marker,
        excluded_markers=excluded_markers,
        burden_model_version=context.constructs[0].burden.model_version,
    )
    plasmid_dir = context.project_output_path / "plasmid_selection"
    artifact_path = plasmid_dir / PLASMID_CANDIDATES_FILENAME
    with manifest_update_lock(context.manifest_path):
        obsolete = ("plasmid_selection", "final_assembly_plan", "final_assembly", "final_design_report")
        update_design_manifest(context.manifest_path, target_compound_id=context.target_compound_id,
                               sections={}, discard_sections=obsolete)
        # Move computed reports under the same lock as current-result invalidation.
        if artifact_path.exists():
            _seed_old_raw_sequences(plasmid_dir, artifact_path)
            from uuid import uuid4
            history = plasmid_dir / "history" / uuid4().hex
            history.mkdir(parents=True)
            os.replace(artifact_path, history / artifact_path.name)
            candidates_dir = plasmid_dir / "candidates"
            if candidates_dir.exists():
                os.replace(candidates_dir, history / "candidates")
        context = load_plasmid_context(config)

    snapshot = fetch_plasmid_snapshot(client=client)
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
            downloaded = _raw_sequence(plasmid_dir, template, downloader)
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
        "Expression burden dynamically determines copy-number fit; --priority contributes at most a three-point adjustment.",
        "Scores are interpretable heuristics and are not experimental success probabilities.",
        "One selected backbone will later be combined separately with every selected expression construct.",
    ]
    if status == "partial":
        warnings.append(
            f"Only {len(candidates)} of {count} requested downloadable candidates were available."
        )
    if status == "failed":
        warnings.append("No candidate passed both metadata and sequence validation.")
    payload: dict[str, Any] = {
        "schema_version": PLASMID_CANDIDATES_SCHEMA_VERSION,
        "algorithm_version": PLASMID_RECOMMENDATION_ALGORITHM_VERSION,
        "status": status,
        "score_type": "interpretable_heuristic_not_probability",
        "target_compound_id": context.target_compound_id,
        "generated_at": _utc_now(),
        "source": {
            "manifest_path": str(context.manifest_path),
            "milvus_collection": snapshot.collection_name,
            "milvus_row_count": snapshot.collection_row_count,
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
            "expression_burden": {
                "model_version": context.constructs[0].burden.model_version,
                "score_range": {
                    "minimum": context.minimum_burden_score,
                    "maximum": context.maximum_burden_score,
                },
                "levels": sorted(
                    {item.burden.level for item in context.constructs}
                ),
                "confidence_levels": sorted(
                    {item.burden.confidence for item in context.constructs}
                ),
            },
            "backbone_policy": "one_backbone_applicable_to_all_selected_constructs",
        },
        "ranking_parameters": parameters,
        "candidate_count": len(candidates),
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


def run_plasmid_recommendation(config: Any, *, client: Any | None = None,
                               downloader: Callable[[Mapping[str, Any]], DownloadedSequence] = download_and_validate_template) -> dict[str, Any]:
    directory = Path(config.project_output_path) / "plasmid_selection"
    directory.mkdir(parents=True, exist_ok=True)
    with FileLock(str(directory / ".recommendation.lock"), timeout=0):
        return _run_plasmid_recommendation(config, client=client, downloader=downloader)


__all__ = [
    "PLASMID_CANDIDATES_FILENAME",
    "plasmid_recommendation_summary",
    "run_plasmid_recommendation",
]
