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

from src.final_assemble_plan.common import (
    relative_project_path,
    stable_json_hash,
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


def _request_fingerprint(
    context_fingerprint: str,
    parameters: Mapping[str, Any],
) -> str:
    return stable_json_hash(
        {
            "algorithm_version": ASSEMBLY_PLAN_ALGORITHM_VERSION,
            "context_input_fingerprint": context_fingerprint,
            "parameters": parameters,
        }
    )


def _plan_set_fingerprint(plans: list[Mapping[str, Any]]) -> str:
    return stable_json_hash(plans)


def _cached_payload(
    path: Path,
    *,
    request_fingerprint: str,
) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    source = payload.get("source")
    plans = payload.get("plans")
    if (
        payload.get("schema_version")
        != ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION
        or payload.get("algorithm_version") != ASSEMBLY_PLAN_ALGORITHM_VERSION
        or payload.get("status") not in {"complete", "partial"}
        or not isinstance(source, Mapping)
        or source.get("request_fingerprint") != request_fingerprint
        or not isinstance(plans, list)
        or not plans
        or int(payload.get("planned_design_count") or 0) != len(plans)
        or payload.get("plan_set_fingerprint")
        != _plan_set_fingerprint(plans)
    ):
        return None
    identities = [item.get("parts_design_id") for item in plans if isinstance(item, Mapping)]
    if len(identities) != len(plans) or len(identities) != len(set(identities)):
        return None
    for plan in plans:
        if plan.get("assembly_method") not in SUPPORTED_METHODS:
            return None
        recorded = str(plan.get("plan_fingerprint") or "")
        unsigned = dict(plan)
        unsigned.pop("plan_fingerprint", None)
        if recorded != stable_json_hash(unsigned):
            return None
    return payload


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
                "file_sha256": construct.file_sha256,
                "sequence_sha256": construct.sequence_sha256,
            },
            "backbone": {
                "plasmid_id": context.backbone.plasmid_id,
                "name": context.backbone.name,
                "path": relative_project_path(
                    context.project_output_path, context.backbone.path
                ),
                "length_bp": context.backbone.length_bp,
                "file_sha256": context.backbone.file_sha256,
                "sequence_content_sha256": (
                    context.backbone.sequence_content_sha256
                ),
                "topology": context.backbone.topology,
                "assembly_policy": context.backbone.assembly_policy,
            },
        }
    )
    plan["plan_fingerprint"] = stable_json_hash(plan)
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
        "manifest_modified": False,
    }


def run_final_assembly_plan(config: Any) -> dict[str, Any]:
    """Generate one best plan per complete expression construct."""

    raw_method = getattr(config, "assembly_method", None)
    requested_method = str(raw_method).strip().lower() if raw_method else None
    if requested_method is not None and requested_method not in SUPPORTED_METHODS:
        raise ValueError(
            f"--method 必须是以下之一：{', '.join(SUPPORTED_METHODS)}"
        )
    context = load_final_assembly_context(config)
    parameters = _parameters(requested_method)
    request_fingerprint = _request_fingerprint(
        context.input_fingerprint,
        parameters,
    )
    output_path = (
        context.project_output_path
        / "final_assemble_plan"
        / ASSEMBLY_PLAN_RECOMMENDATIONS_FILENAME
    )
    cached = _cached_payload(
        output_path,
        request_fingerprint=request_fingerprint,
    )
    if cached is not None:
        return assembly_plan_summary(
            cached, output_path, reused_existing=True
        )

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
            "manifest_revision": context.manifest_revision,
            "context_input_fingerprint": context.input_fingerprint,
            "parts_selection_fingerprint": context.parts_selection_fingerprint,
            "assembled_constructs_fingerprint": (
                context.assembled_constructs_fingerprint
            ),
            "plasmid_selection_fingerprint": (
                context.plasmid_selection_fingerprint
            ),
            "request_fingerprint": request_fingerprint,
        },
        "parameters": parameters,
        "design_count": len(context.constructs),
        "planned_design_count": len(plans),
        "failed_design_count": len(failures),
        "method_counts": dict(sorted(method_counts.items())),
        "plan_set_fingerprint": _plan_set_fingerprint(plans),
        "plans": plans,
        "failures": failures,
        "warnings": warnings,
        "manifest_modified": False,
    }
    _write_json_atomic(output_path, payload)
    return assembly_plan_summary(
        payload, output_path, reused_existing=False
    )


__all__ = [
    "assembly_plan_summary",
    "run_final_assembly_plan",
]
