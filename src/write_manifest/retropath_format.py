"""RetroPath-specific projection and restoration for the shared compact format."""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any
import warnings

from src.write_manifest.solution_format import _present, _relative, _resolve, _true, _values

_META_ALIASES = {
    "step_id": "retropath_step_id", "hypothesis_id": "retropath_hypothesis_id",
    "rule_ids": "retropath_rule_id", "source_mnxr_id": "source_mnxr_id",
    "source_reaction_ids": "source_reaction_ids", "source_ec_numbers": "source_ec_numbers",
    "source_uniprot_ids": "source_uniprot_ids", "source_rhea_ids": "source_rhea_ids",
    "exact_kegg_reaction_ids": "exact_kegg_reaction_ids", "exact_rhea_ids": "exact_rhea_ids",
    "formal_mapping_exact": "formal_mapping_exact", "reaction_signature_sha256": "reaction_signature_sha256",
    "full_reaction_smiles": "full_reaction_smiles", "core_reaction_smiles": "core_reaction_smiles",
}
_LIST_META = frozenset({"rule_ids", "source_reaction_ids", "source_ec_numbers", "source_uniprot_ids",
                        "source_rhea_ids", "exact_kegg_reaction_ids", "exact_rhea_ids"})
_TRUNCATION = ("combination_truncated", "upstream_enumeration_truncated", "candidate_top_k_truncated")
_STATUS_FIELDS = ("validation_status", "stoichiometry_status", "gem_status")


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _pick(name: str, values: list[Any], normalize=None) -> Any:
    present = [value for value in values if _present(value)]
    if not present:
        return None
    normalized = [normalize(value) if normalize else value for value in present]
    first = normalized[0]
    for value in normalized[1:]:
        equal = set(first) == set(value) if isinstance(first, list) and isinstance(value, list) else first == value
        if not equal:
            raise ValueError(f"RetroPath identity/binding conflict for {name}: {first!r} != {value!r}")
    return deepcopy(first)


def _stage(status: Any) -> str:
    value = str(status or "not_run")
    if value.startswith("PASS_"):
        return "passed"
    if value.startswith("FAIL_"):
        return "failed"
    return value


def _prediction(solution: Mapping[str, Any], path: Path) -> dict[str, Any]:
    old = _mapping(solution.get("prediction"))
    summary, validation = _mapping(solution.get("summary")), _mapping(solution.get("validation"))
    result = {}
    for key in ("candidate_rank", "candidate_id", "combination_id"):
        value = _pick(key, [old.get(key), summary.get("retropath_" + key), validation.get("retropath_" + key)],
                      int if key == "candidate_rank" else str)
        if value is not None:
            result[key] = value
    for suffix in ("id", "schema_version", "manifest"):
        value = _pick("materialization_" + suffix,
            [old.get("materialization_" + suffix), old.get("promotion_" + suffix),
             summary.get("materialization_" + suffix), summary.get("promotion_" + suffix)],
            (lambda value: _relative(value, path)) if suffix == "manifest" else str)
        if value is not None:
            result["materialization_" + suffix] = value
    for key in ("review_status", "review_required"):
        value = old.get(key)
        if value is None and key == "review_required":
            value = summary.get("prediction_review_required")
        if value is not None:
            result[key] = _true(value) if key == "review_required" else value
    for key in _TRUNCATION:
        value = _pick(key, [old.get(key), validation.get(key), summary.get(key)], _true)
        if value is not None:
            result[key] = value
    return result


def _validation(solution: Mapping[str, Any], core: dict[str, Any]) -> dict[str, Any]:
    old, summary = _mapping(solution.get("validation")), _mapping(solution.get("summary"))
    prediction = _mapping(solution.get("prediction"))
    result = deepcopy(core)
    if not result.get("validation_status"):
        value = prediction.get("validation_status") or summary.get("validation_status")
        if value:
            result["validation_status"] = value
    for key in ("stoichiometry_status", "gem_status"):
        # Solver/overlay outcomes supersede the earlier prediction snapshot;
        # these are states at different times, not conflicting identities.
        value = next((item for item in (old.get(key), prediction.get(key), summary.get(key)) if _present(item)), None)
        if value is not None:
            result[key] = value
    for key in ("cofactor_mode", "opened_generic_compound_ids"):
        if not _present(result.get(key)) and _present(prediction.get(key)):
            result[key] = _values(prediction[key]) if key == "opened_generic_compound_ids" else prediction[key]
    if not result.get("issues") and prediction.get("validation_issues"):
        result["issues"] = deepcopy(prediction["validation_issues"])
    return result


def _validation_view(validation: Mapping[str, Any]) -> dict[str, Any]:
    result = {key: deepcopy(validation[key]) for key in ("stoichiometry_status", "gem_status", "cofactor_mode") if key in validation}
    result["validation_status"] = _stage(validation.get("validation_status"))
    result["cofactor_relaxed"] = str(validation.get("cofactor_mode") or "").startswith("relaxed")
    result["opened_generic_compound_ids"] = _values(validation.get("opened_generic_compound_ids"))
    issues = validation.get("issues")
    result["validation_issues"] = deepcopy(issues) if isinstance(issues, list) else [issues] if issues else []
    return result


def _step(raw: Mapping[str, Any], core: dict[str, Any], prediction: Mapping[str, Any], validation: Mapping[str, Any], *, compact_input: bool) -> dict[str, Any]:
    result = deepcopy(core)
    source = str(raw.get("step_source") or "")
    result["step_source"] = source
    for key in ("ec_status", "enzyme_search_eligible"):
        if _present(raw.get(key)):
            result[key] = _true(raw[key]) if key == "enzyme_search_eligible" else raw[key]
    # Equal locks are recoverable from annotation; explicit empty/different locks are not.
    ecs = _values(raw.get("enzyme_ecs"))
    locks = _values(raw.get("locked_enzyme_ecs", ecs if compact_input else None))
    if locks != ecs:
        result["locked_enzyme_ecs"] = locks
    dependencies = _values(raw.get("depends_on_step_ids"))
    if dependencies:
        result["depends_on_step_ids"] = dependencies
    if "stoichiometry_terms" in raw and raw["stoichiometry_terms"] is not None:
        result["stoichiometry_terms"] = deepcopy(raw["stoichiometry_terms"])
    overrides = {}
    for key in _STATUS_FIELDS:
        value = _mapping(raw.get("validation")).get(key, raw.get(key))
        global_value = _stage(validation.get(key)) if key == "validation_status" else validation.get(key)
        if _present(value) and value != global_value:
            overrides[key] = value
    if overrides:
        result["validation"] = overrides
    review = _true(raw.get("prediction_review_required", prediction.get("review_required", True) if compact_input else False))
    if review != prediction.get("review_required", True):
        result["prediction_review_required"] = review
    old = _mapping(raw.get("prediction_provenance"))
    if source != "retropath":
        # Expansion witnesses can carry additional source reaction identities.
        ids = _values(raw.get("source_reaction_ids"))
        if ids and ids != [str(raw.get("reaction_id") or "")]:
            result["source_reaction_ids"] = ids
        return result
    meta = deepcopy(dict(old))
    for key in ("candidate_rank", "candidate_id", "combination_id"):
        value = _pick(key, [prediction.get(key), old.get(key)], int if key == "candidate_rank" else str)
        meta.pop(key, None)
        if value is not None and not _present(prediction.get(key)):
            meta[key] = value
    for suffix in ("id", "schema_version", "manifest"):
        key = "materialization_" + suffix
        value = _pick(key, [old.get(key), old.get("promotion_" + suffix)])
        meta.pop("promotion_" + suffix, None)
        meta.pop(key, None)
        if value is not None and value != prediction.get(key):
            meta[key] = value
    for key, alias in _META_ALIASES.items():
        fallback = old.get("rule_id") if key == "rule_ids" else None
        normalize = _values if key in _LIST_META else _true if key == "formal_mapping_exact" else str
        value = _pick(key, [raw.get(alias), old.get(key), fallback], normalize)
        meta.pop(key, None)
        if value is not None and (_present(value) or key == "formal_mapping_exact"):
            meta[key] = value
    meta.pop("rule_id", None)
    sid, hid = meta.get("step_id"), meta.get("hypothesis_id")
    if not sid or str(raw.get("reaction_id") or "") not in {sid, hid}:
        raise ValueError("RetroPath step/hypothesis identity is inconsistent")
    for key in (*_STATUS_FIELDS, "review_status"):
        global_value = prediction.get(key)
        if key in _STATUS_FIELDS:
            global_value = _stage(validation.get(key)) if key == "validation_status" else validation.get(key)
        # Provenance records the state at prediction creation. A later solver
        # outcome belongs to the current step, not this historical annotation.
        value = old.get(key)
        meta.pop(key, None)
        if _present(value) and value != global_value:
            meta[key] = value
    meta.pop("prediction_review_required", None)
    result["prediction_provenance"] = meta
    return result


def compact_retropath_solution(solution: Mapping[str, Any], core: dict[str, Any], path: Path) -> dict[str, Any]:
    result = deepcopy(core)
    prediction = _prediction(solution, path)
    validation = _validation(solution, result.get("validation", {}))
    previous = _mapping(solution.get("prediction"))
    origin = deepcopy(dict(_mapping(previous.get("origin_validation"))))
    for key, current in _validation_view(validation).items():
        if key in previous and previous[key] is not None and previous[key] != current:
            origin[key] = deepcopy(previous[key])
    if origin:
        prediction["origin_validation"] = origin
    result.update(solution_source="retropath", prediction=prediction, validation=validation)
    compact_input = isinstance(solution.get("compound_names"), Mapping)
    result["steps"] = [_step(raw, step, prediction, validation, compact_input=compact_input)
        for raw, step in zip(solution.get("steps", []), result.get("steps", []))]
    return result


def restore_retropath_runtime(solution: dict[str, Any], path: Path) -> None:
    prediction, validation = solution["prediction"], solution.get("validation", {})
    prediction["source"] = "retropath"
    for suffix in ("id", "schema_version", "manifest"):
        key = "materialization_" + suffix
        if key in prediction:
            if suffix == "manifest":
                prediction[key] = str(_resolve(prediction[key], path))
            prediction["promotion_" + suffix] = prediction[key]
    current_validation = _validation_view(validation)
    prediction.update(deepcopy(current_validation))
    prediction.update(deepcopy(dict(_mapping(prediction.get("origin_validation")))))
    validation["cofactor_relaxed"] = current_validation["cofactor_relaxed"]
    hypotheses = []
    for step in solution.get("steps", []):
        step["locked_enzyme_ecs"] = ";".join(_values(step.get("locked_enzyme_ecs", step.get("enzyme_ecs"))))
        step["prediction_review_required"] = step.get("prediction_review_required", prediction.get("review_required", True))
        for key in _STATUS_FIELDS:
            step[key] = _mapping(step.get("validation")).get(key, current_validation.get(key))
        if step.get("step_source") != "retropath":
            continue
        meta = step["prediction_provenance"]
        for key in ("candidate_rank", "candidate_id", "combination_id", "materialization_id", "materialization_schema_version"):
            if key not in meta and key in prediction:
                meta[key] = deepcopy(prediction[key])
        if meta.get("materialization_id"):
            meta["promotion_id"] = meta["materialization_id"]
        for key in (*_STATUS_FIELDS, "review_status"):
            meta.setdefault(key, current_validation.get(key) if key in _STATUS_FIELDS else prediction.get(key))
        meta["prediction_review_required"] = step["prediction_review_required"]
        for key, alias in _META_ALIASES.items():
            value = meta.get(key)
            if key in _LIST_META:
                value = ";".join(_values(value))
            step[alias] = deepcopy(value) if value is not None else False if key == "formal_mapping_exact" else ""
        step["resolution_action"] = "retropath_prediction"
        if meta.get("hypothesis_id"):
            hypotheses.append(meta["hypothesis_id"])
    prediction["stoichiometry_hypothesis_ids"] = list(dict.fromkeys(hypotheses))
    for key in ("candidate_rank", "candidate_id", "combination_id"):
        validation["retropath_" + key] = prediction.get(key)
    validation["stoichiometry_hypothesis_ids"] = prediction["stoichiometry_hypothesis_ids"]
    for key in _TRUNCATION:
        if key in prediction:
            solution["summary"][key] = prediction[key]
    if "combination_truncated" in prediction:
        validation["combination_truncated"] = prediction["combination_truncated"]
    summary = solution["summary"]
    summary.update(solution_source="retropath", prediction_review_required=prediction.get("review_required", True))
    summary["requires_review"] = summary["requires_review"] or bool(summary["prediction_review_required"])
    summary["reaction_resolution_status"] = "predicted_core_stoichiometry" if validation.get("stoichiometry_status") == "core_only" else "predicted_stoichiometry"
    evidence = prediction.get("materialization_manifest")
    if evidence and not Path(evidence).is_file():
        prediction["evidence_available"] = False
        warnings.warn("RetroPath 来源明细不可用 (prediction evidence unavailable): " + str(evidence), RuntimeWarning, stacklevel=2)
