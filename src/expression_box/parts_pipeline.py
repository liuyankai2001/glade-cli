"""Manifest-driven expression-parts recommendation pipeline."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from src.expression_box.config import (
    DEFAULT_EXPRESSION_PARTS_DESIGN_COUNT,
    EXPRESSION_PARTS_CANDIDATE_POOL_MAX,
    EXPRESSION_PARTS_CANDIDATE_POOL_MIN,
    EXPRESSION_PARTS_CANDIDATE_POOL_MULTIPLIER,
    EXPRESSION_PARTS_DESIGNS_SCHEMA_VERSION,
    EXPRESSION_SUCCESS_ASSEMBLY_WEIGHT,
    EXPRESSION_SUCCESS_EVIDENCE_WEIGHT,
    EXPRESSION_SUCCESS_MIN_SCORE,
    EXPRESSION_SUCCESS_PROMOTER_WEIGHT,
    EXPRESSION_SUCCESS_TERMINATOR_WEIGHT,
    EXPRESSION_SUCCESS_TRANSLATION_WEIGHT,
    MAX_EXPRESSION_PARTS_DESIGN_COUNT,
    MIN_EXPRESSION_PARTS_DESIGN_COUNT,
    PARTS_RECOMMENDATION_ALGORITHM_VERSION,
    RBS_CONTEXT_CDS_PREFIX_NT,
    RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT,
    RBS_SHORTLIST_PER_STRENGTH,
)
from src.expression_box.milvus_parts import fetch_expression_part_candidates
from src.expression_box.ostir_adapter import backend_versions
from src.expression_box.parts_manifest_adapter import load_expression_parts_context
from src.expression_box.recommend_expression_parts import (
    generate_expression_parts_designs,
)


EXPRESSION_PARTS_DESIGNS_FILENAME = "expression_parts_designs.json"


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _stable_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _cached_payload(path: Path, input_fingerprint: str) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    source = payload.get("source") if isinstance(payload, dict) else None
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != EXPRESSION_PARTS_DESIGNS_SCHEMA_VERSION
        or payload.get("status") not in {"complete", "partial", "failed"}
        or not isinstance(source, dict)
        or source.get("input_fingerprint") != input_fingerprint
        or not isinstance(payload.get("designs"), list)
        or (payload.get("status") != "failed" and not payload["designs"])
        or (payload.get("status") == "failed" and payload["designs"])
    ):
        return None
    return payload


def _summary(
    payload: dict[str, Any],
    output_path: Path,
    *,
    reused_existing: bool,
) -> dict[str, Any]:
    designs = payload["designs"]
    ranking = payload["ranking"]
    primary = next(
        (design for design in designs if bool(design.get("recommended"))),
        None,
    )
    high_confidence_count = sum(
        float(design["expression_success_score"]) >= 85.0 for design in designs
    )
    return {
        "ok": payload["status"] == "complete",
        "status": payload["status"],
        "target_compound_id": payload["target_compound_id"],
        "requested_design_count": ranking["requested_design_count"],
        "design_count": payload["design_count"],
        "minimum_success_score": ranking["minimum_success_score"],
        "score_range": ranking["selected_score_range"],
        "score_tiers": {
            "high_confidence_85_or_above": high_confidence_count,
            "qualified_70_to_84_99": len(designs) - high_confidence_count,
        },
        "primary_recommendation": (
            {
                "design_id": primary["design_id"],
                "score": primary["expression_success_score"],
                "expression_regime": primary["expression_regime"],
            }
            if primary is not None
            else None
        ),
        "skipped_strategies": payload["skipped_strategies"],
        "warnings": payload["warnings"],
        "output_path": str(output_path),
        "reused_existing": reused_existing,
        "manifest_modified": False,
    }


def run_expression_parts_design(
    config: Any,
    *,
    milvus_client: Any | None = None,
    predictor: Any | None = None,
) -> dict[str, Any]:
    """Retrieve remote facts, score contexts, and persist ranked designs."""

    requested_design_count = int(
        getattr(config, "n_designs", DEFAULT_EXPRESSION_PARTS_DESIGN_COUNT)
    )
    if not (
        MIN_EXPRESSION_PARTS_DESIGN_COUNT
        <= requested_design_count
        <= MAX_EXPRESSION_PARTS_DESIGN_COUNT
    ):
        raise ValueError(
            "n_designs must be between "
            f"{MIN_EXPRESSION_PARTS_DESIGN_COUNT} and "
            f"{MAX_EXPRESSION_PARTS_DESIGN_COUNT}"
        )

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    context = load_expression_parts_context(manifest_path, project_root)
    snapshot = fetch_expression_part_candidates(
        context,
        client=milvus_client,
    )
    backends = backend_versions()
    backends["dnachisel"] = __import__("dnachisel").__version__
    backends["pymilvus"] = __import__("pymilvus").__version__
    input_fingerprint = _stable_hash(
        {
            "algorithm_version": PARTS_RECOMMENDATION_ALGORITHM_VERSION,
            "context_input_fingerprint": context.input_fingerprint,
            "collection_name": snapshot.collection_name,
            "collection_schema_fingerprint": snapshot.schema_fingerprint,
            "candidate_snapshot_fingerprint": snapshot.candidate_fingerprint,
            "backends": backends,
            "parameters": {
                "requested_design_count": requested_design_count,
                "rbs_shortlist_per_strength": RBS_SHORTLIST_PER_STRENGTH,
                "rbs_context_cds_prefix_nt": RBS_CONTEXT_CDS_PREFIX_NT,
                "rbs_context_previous_cds_suffix_nt": (
                    RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT
                ),
                "candidate_pool_min": EXPRESSION_PARTS_CANDIDATE_POOL_MIN,
                "candidate_pool_multiplier": (
                    EXPRESSION_PARTS_CANDIDATE_POOL_MULTIPLIER
                ),
                "candidate_pool_max": EXPRESSION_PARTS_CANDIDATE_POOL_MAX,
                "minimum_success_score": EXPRESSION_SUCCESS_MIN_SCORE,
                "success_score_weights": {
                    "evidence": EXPRESSION_SUCCESS_EVIDENCE_WEIGHT,
                    "promoter": EXPRESSION_SUCCESS_PROMOTER_WEIGHT,
                    "translation": EXPRESSION_SUCCESS_TRANSLATION_WEIGHT,
                    "terminator": EXPRESSION_SUCCESS_TERMINATOR_WEIGHT,
                    "assembly": EXPRESSION_SUCCESS_ASSEMBLY_WEIGHT,
                },
            },
        }
    )
    output_path = (
        project_root
        / "expression_box"
        / EXPRESSION_PARTS_DESIGNS_FILENAME
    )
    cached = _cached_payload(output_path, input_fingerprint)
    if cached is not None:
        try:
            return _summary(cached, output_path, reused_existing=True)
        except (KeyError, TypeError, ValueError):
            pass

    generation_kwargs: dict[str, Any] = {}
    if predictor is not None:
        generation_kwargs["predictor"] = predictor
    generation = generate_expression_parts_designs(
        context,
        snapshot,
        requested_design_count=requested_design_count,
        **generation_kwargs,
    )
    payload = {
        "schema_version": EXPRESSION_PARTS_DESIGNS_SCHEMA_VERSION,
        "status": generation["status"],
        "generated_at": _utc_now(),
        "algorithm_version": PARTS_RECOMMENDATION_ALGORITHM_VERSION,
        "source": {
            "manifest_path": str(manifest_path),
            "manifest_revision": context.manifest_revision,
            "expression_box_selection_fingerprint": (
                context.expression_box_selection_fingerprint
            ),
            "cds_selection_source_fingerprint": (
                context.cds_selection_source_fingerprint
            ),
            "context_input_fingerprint": context.input_fingerprint,
            "milvus_collection": snapshot.collection_name,
            "milvus_schema_fingerprint": snapshot.schema_fingerprint,
            "candidate_snapshot_fingerprint": snapshot.candidate_fingerprint,
            "input_fingerprint": input_fingerprint,
        },
        "target_compound_id": context.target_compound_id,
        "host": {
            "name": context.host_name,
            "key": context.host_key,
            "milvus_labels": list(context.host_labels),
        },
        "backends": backends,
        "retrieval": {
            "source": "remote_milvus",
            "collection": snapshot.collection_name,
            "collection_row_count": snapshot.collection_row_count,
            "raw_candidate_counts": snapshot.raw_counts,
            "accepted_candidate_counts": snapshot.accepted_counts,
            "rejected_candidate_counts": snapshot.rejected_counts,
            "vector_search_used": False,
            "embedding_model_used": False,
        },
        "parameters": {
            "requested_design_count": requested_design_count,
            "rbs_shortlist_per_strength": RBS_SHORTLIST_PER_STRENGTH,
            "rbs_context_cds_prefix_nt": RBS_CONTEXT_CDS_PREFIX_NT,
            "rbs_context_previous_cds_suffix_nt": (
                RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT
            ),
            "candidate_pool_min": EXPRESSION_PARTS_CANDIDATE_POOL_MIN,
            "candidate_pool_multiplier": (
                EXPRESSION_PARTS_CANDIDATE_POOL_MULTIPLIER
            ),
            "candidate_pool_max": EXPRESSION_PARTS_CANDIDATE_POOL_MAX,
            "minimum_success_score": EXPRESSION_SUCCESS_MIN_SCORE,
            "success_score_weights": {
                "evidence": EXPRESSION_SUCCESS_EVIDENCE_WEIGHT,
                "promoter": EXPRESSION_SUCCESS_PROMOTER_WEIGHT,
                "translation": EXPRESSION_SUCCESS_TRANSLATION_WEIGHT,
                "terminator": EXPRESSION_SUCCESS_TERMINATOR_WEIGHT,
                "assembly": EXPRESSION_SUCCESS_ASSEMBLY_WEIGHT,
            },
        },
        "ranking": generation["ranking"],
        "design_count": len(generation["designs"]),
        "designs": generation["designs"],
        "skipped_strategies": generation["skipped_strategies"],
        "prediction_failure_count": generation["prediction_failure_count"],
        "prediction_failures": generation["prediction_failures"],
        "warnings": generation["warnings"],
        "manifest_modified": False,
    }
    _write_json_atomic(output_path, payload)
    return _summary(payload, output_path, reused_existing=False)


__all__ = [
    "EXPRESSION_PARTS_DESIGNS_FILENAME",
    "run_expression_parts_design",
]
