"""Generate one recommended assembly plan for every selected design."""

from __future__ import annotations

import json
import os
import tempfile
from collections import Counter
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from filelock import FileLock
from src.write_manifest.store import manifest_update_lock, update_design_manifest

from src.final_assemble_plan.common import (
    relative_project_path,
)
from src.final_assemble_plan.config import (
    ASSEMBLY_PLAN_ALGORITHM_VERSION,
    ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME,
    ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION,
    DEFAULT_HOMOLOGY_ARM_LENGTH,
    PLAN_SCORE_WEIGHTS,
    SUPPORTED_METHODS,
)
from src.final_assemble_plan.get_final_assembly_context import (
    load_final_assembly_context,
)
from src.final_assemble_plan.recommend_gibson_insertion_sites import (
    recommend_gibson_plans,
)
from src.final_assemble_plan.recommend_restriction_sites import (
    recommend_restriction_plans,
)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
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


def _parameters(requested_method: str | None) -> dict[str, Any]:
    return {
        "selection_mode": (
            "user_method" if requested_method is not None else "auto"
        ),
        "requested_method": requested_method,
        "supported_methods": list(SUPPORTED_METHODS),
        "homology_arm_length": DEFAULT_HOMOLOGY_ARM_LENGTH,
        "restriction_site_retention": "retain",
        "score_weights": PLAN_SCORE_WEIGHTS,
        "gibson_linearization_modes": ["restriction", "pcr"],
    }


def _candidate_rank(candidate: Mapping[str, Any]) -> tuple[Any, ...]:
    method = str(candidate.get("assembly_method") or "")
    linearization = candidate.get("backbone_linearization")
    linearization = linearization if isinstance(linearization, Mapping) else {}
    return (
        -float(candidate.get("score") or 0.0),
        0 if method == "restriction" else 1,
        0 if linearization.get("mode") == "restriction" else 1,
        str(linearization.get("enzyme_summary") or ""),
    )


def _candidates_for_design(
    context: Any,
    construct: Any,
    requested_method: str | None,
) -> tuple[list[dict[str, Any]], list[str]]:
    methods = (requested_method,) if requested_method else SUPPORTED_METHODS
    candidates: list[dict[str, Any]] = []
    attempted: list[str] = []
    for method in methods:
        attempted.append(str(method))
        if method == "restriction":
            candidates.extend(
                recommend_restriction_plans(context, construct)
            )
        elif method == "gibson":
            candidates.extend(
                recommend_gibson_plans(context, construct)
            )
    candidates.sort(key=_candidate_rank)
    return candidates, attempted


def _decorate_plan(
    context: Any,
    construct: Any,
    candidate: Mapping[str, Any],
    *,
    requested_method: str | None,
) -> dict[str, Any]:
    plan = dict(candidate)
    plan.update(
        {
            "selection_source": (
                "user_method" if requested_method else "system_recommended"
            ),
            "insert": {
                "path": relative_project_path(
                    context.project_output_path, construct.path
                ),
                "format": "genbank",
                "length_bp": construct.length_bp,
            },
            "backbone": {
                "plasmid_id": context.backbone.plasmid_id,
                "name": context.backbone.name,
                "path": relative_project_path(
                    context.project_output_path, context.backbone.path
                ),
                "length_bp": context.backbone.length_bp,
                "topology": context.backbone.topology,
                "assembly_policy": context.backbone.assembly_policy,
            },
        }
    )
    return plan


def assembly_plan_summary(
    payload: Mapping[str, Any],
    output_path: Path,
    *,
    reused_existing: bool,
) -> dict[str, Any]:
    plans = payload.get("plans")
    plans = plans if isinstance(plans, list) else []
    failures = payload.get("failures")
    failures = failures if isinstance(failures, list) else []
    return {
        "ok": payload.get("status") == "complete",
        "status": payload.get("status"),
        "target_compound_id": payload.get("target_compound_id"),
        "selection_mode": payload.get("selection_mode"),
        "requested_method": payload.get("requested_method"),
        "design_count": payload.get("design_count"),
        "planned_design_count": len(plans),
        "failed_design_count": len(failures),
        "plans": [
            {
                "parts_design_id": plan.get("parts_design_id"),
                "method": plan.get("assembly_method"),
                "score": plan.get("score"),
                "target": plan.get("target"),
                "linearization": (
                    plan.get("backbone_linearization", {}).get("mode")
                    if isinstance(plan.get("backbone_linearization"), Mapping)
                    else None
                ),
                "enzymes": (
                    plan.get("backbone_linearization", {}).get(
                        "enzyme_summary"
                    )
                    if isinstance(plan.get("backbone_linearization"), Mapping)
                    else None
                ),
            }
            for plan in plans
        ],
        "failures": failures,
        "warnings": payload.get("warnings", []),
        "output_path": str(output_path.resolve()),
        "reused_existing": reused_existing,
        "manifest_modified": True,
    }


def _run_final_assembly_plan(config: Any) -> dict[str, Any]:
    """Generate one best plan per complete expression construct."""

    raw_method = getattr(config, "assembly_method", None)
    requested_method = str(raw_method).strip().lower() if raw_method else None
    if requested_method is not None and requested_method not in SUPPORTED_METHODS:
        raise ValueError(
            f"--method 必须是以下之一：{', '.join(SUPPORTED_METHODS)}"
        )
    context = load_final_assembly_context(config)
    parameters = _parameters(requested_method)
    output_path = (context.project_output_path / "final_assemble_plan" / ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME)
    with manifest_update_lock(context.manifest_path):
        obsolete = ("final_assembly_plan", "final_assembly", "final_design_report")
        update_design_manifest(context.manifest_path, target_compound_id=context.target_compound_id,
                               sections={}, discard_sections=obsolete)
        if output_path.exists():
            from uuid import uuid4
            history = output_path.parent / "history" / uuid4().hex
            history.mkdir(parents=True)
            os.replace(output_path, history / output_path.name)
        context = load_final_assembly_context(config)

    plans: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for construct in context.constructs:
        candidates, attempted = _candidates_for_design(
            context,
            construct,
            requested_method,
        )
        if not candidates:
            failures.append(
                {
                    "parts_design_id": construct.design_id,
                    "attempted_methods": attempted,
                    "reason": (
                        "no safe insertion target, enzyme pair, or homology-arm plan passed all gates"
                    ),
                }
            )
            continue
        plans.append(
            _decorate_plan(
                context,
                construct,
                candidates[0],
                requested_method=requested_method,
            )
        )
    plans.sort(key=lambda item: int(item["parts_design_id"]))
    failures.sort(key=lambda item: int(item["parts_design_id"]))
    status = (
        "complete"
        if len(plans) == len(context.constructs)
        else "partial"
        if plans
        else "failed"
    )
    method_counts = Counter(str(item["assembly_method"]) for item in plans)
    warnings = [
        "Assembly-plan scores are transparent engineering heuristics, not experimental success probabilities.",
        "Restriction digest conditions, primer annealing segments, and reaction conditions require experimental confirmation.",
        "This command writes no final GenBank or FASTA files.",
    ]
    if failures:
        warnings.append(
            f"{len(failures)} design(s) have no complete plan; write --assembly-plan is disabled."
        )
    payload: dict[str, Any] = {
        "schema_version": ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION,
        "algorithm_version": ASSEMBLY_PLAN_ALGORITHM_VERSION,
        "status": status,
        "generated_at": _utc_now(),
        "target_compound_id": context.target_compound_id,
        "selection_mode": (
            "user_method" if requested_method else "system_recommended"
        ),
        "requested_method": requested_method,
        "source": {
            "manifest_path": str(context.manifest_path),
        },
        "parameters": parameters,
        "design_count": len(context.constructs),
        "planned_design_count": len(plans),
        "failed_design_count": len(failures),
        "method_counts": dict(sorted(method_counts.items())),
        "plans": plans,
        "failures": failures,
        "warnings": warnings,
        "manifest_modified": True,
    }
    _write_json_atomic(output_path, payload)
    return assembly_plan_summary(
        payload, output_path, reused_existing=False
    )


def run_final_assembly_plan(config: Any) -> dict[str, Any]:
    directory = Path(config.project_output_path) / "final_assemble_plan"
    directory.mkdir(parents=True, exist_ok=True)
    with FileLock(str(directory / ".plan.lock"), timeout=0):
        return _run_final_assembly_plan(config)


__all__ = [
    "assembly_plan_summary",
    "run_final_assembly_plan",
]
