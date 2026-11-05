"""Compact route storage with a legacy-compatible, derived runtime view."""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import os
from pathlib import Path
import re
from typing import Any
import warnings

_FLAGS = {
    "unclear_reaction": re.compile(r"unclear[\s_-]+reaction", re.I),
    "multi_step_reaction": re.compile(r"multi[\s_-]*step[\s_-]+reaction", re.I),
    "incomplete_reaction": re.compile(r"incomplete[\s_-]+reaction", re.I),
}
_FLAG_LABELS = {
    "unclear_reaction": "unclear reaction",
    "multi_step_reaction": "multi-step reaction",
    "incomplete_reaction": "incomplete reaction",
}
_STEP_FIELDS = ("step_index", "reaction_id", "reaction_name", "status", "equation", "direction", "produced_compound_id")
_STEP_LISTS = ("precursor_compound_ids", "enzyme_ecs", "ko_ids", "rhea_ids")
_VALIDATION_FIELDS = ("validation_status", "cofactor_mode", "fba_product_flux", "fba_growth_flux", "issues")
_REPORT_FIELDS = ("summary_file", "route_fluxes_file", "manifest_file")
_SPONTANEOUS = re.compile(r"spontaneous|non[\s-]*enzymatic", re.I)


def _values(value: Any) -> list[str]:
    raw = value if isinstance(value, (list, tuple)) else re.split(r"[;|]", str(value or ""))
    return list(dict.fromkeys(str(item).strip() for item in raw if str(item).strip()))


def _present(value: Any) -> bool:
    return value is not None and value != "" and value != [] and value != {}


def _true(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _retropath(solution: Mapping[str, Any]) -> bool:
    summary = solution.get("summary") or {}
    prediction = solution.get("prediction") or {}
    return (
        str(solution.get("solution_source") or "").lower() == "retropath"
        or isinstance(summary, Mapping) and str(summary.get("solution_source") or "").lower() == "retropath"
        or isinstance(prediction, Mapping) and str(prediction.get("source") or "").lower() == "retropath"
        or any(str(step.get("step_source") or "").lower() == "retropath"
               for step in solution.get("steps", []) if isinstance(step, Mapping))
    )


def _resolve(value: Any, manifest_path: Path) -> Path:
    path = Path(str(value)).expanduser()
    return (path if path.is_absolute() else manifest_path.resolve().parent / path).resolve()


def _relative(value: Any, manifest_path: Path) -> str:
    path = _resolve(value, manifest_path)
    try:
        return Path(os.path.relpath(path, manifest_path.resolve().parent)).as_posix()
    except ValueError:  # Legacy projects may reference another Windows drive.
        return str(path)


def _labels(value: Any) -> list[tuple[str, str]]:
    result = []
    for label in _values(value):
        match = re.fullmatch(r"([CDG]\d{5})\s*\((.*)\)", label)
        if match:
            result.append((match[1], match[2]))
    return result


def _flags(step: Mapping[str, Any]) -> list[str]:
    annotation = str(step.get("reaction_comment") or "")
    explicit = _values(step.get("reaction_flags"))
    return [key for key, pattern in _FLAGS.items()
            if key in explicit or _true(step.get(key)) or pattern.search(annotation)]


def _comment(step: Mapping[str, Any]) -> str:
    text = str(step.get("reaction_comment") or "")
    action = str(step.get("resolution_action") or "")
    if _SPONTANEOUS.search(action) and not _SPONTANEOUS.search(text):
        text = "; ".join(filter(None, (text, action)))
    meaningful = []
    for clause in re.split(r"[;|\r\n]+", text):
        remainder = clause
        for pattern in _FLAGS.values():
            remainder = pattern.sub("", remainder)
        # Labels can stand alone or be concatenated. Preserve complete clauses
        # when they contain an explanation, even if they mention a label.
        if remainder.strip(" ,.:;|-\t"):
            meaningful.append(clause.strip())
    return "; ".join(meaningful)


def compact_kegg_solution(
    solution: Mapping[str, Any], *, manifest_path: str | Path,
    target_compound_id: str, validation_files: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Project a selected pure KEGG route; never compact a RetroPath route."""
    if _retropath(solution):
        return deepcopy(dict(solution))
    return _compact_solution_core(solution, manifest_path=manifest_path,
        target_compound_id=target_compound_id, validation_files=validation_files)


def _compact_solution_core(
    solution: Mapping[str, Any], *, manifest_path: str | Path,
    target_compound_id: str, validation_files: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    path = Path(manifest_path)
    summary = solution.get("summary") or {}
    names = {str(key): str(value) for key, value in (solution.get("compound_names") or {}).items() if _present(value)}
    for key, name in _labels(summary.get("reachable_anchor_labels")):
        names.setdefault(key, name)
    for step in solution.get("steps", []):
        for key, name in _labels(step.get("precursor_compound_labels")):
            names.setdefault(key, name)
        if _present(step.get("produced_compound_name")):
            names.setdefault(str(step["produced_compound_id"]), str(step["produced_compound_name"]))
    if _present(summary.get("target_compound_name")):
        names.setdefault(target_compound_id, str(summary["target_compound_name"]))
    steps = []
    for raw in solution.get("steps", []):
        step = {key: deepcopy(raw[key]) for key in _STEP_FIELDS if _present(raw.get(key))}
        step.update({key: _values(raw[key]) for key in _STEP_LISTS if _values(raw.get(key))})
        flags, comment = _flags(raw), _comment(raw)
        if flags:
            step["reaction_flags"] = flags
        if comment:
            step["reaction_comment"] = comment
        steps.append(step)
    raw_validation = solution.get("validation") or {}
    validation = {key: deepcopy(raw_validation[key]) for key in _VALIDATION_FIELDS if _present(raw_validation.get(key))}
    opened = _values(raw_validation.get("opened_generic_compound_ids"))
    if opened:
        validation["opened_generic_compound_ids"] = opened
    for key in _REPORT_FIELDS:
        if validation_files is not None:
            value = validation_files.get(key)
            if value and _resolve(value, path).is_file():
                validation[key] = _relative(value, path)
        elif raw_validation.get(key):
            # Existing evidence references survive a later missing file.
            validation[key] = _relative(raw_validation[key], path)
    result = {key: solution[key] for key in ("solution_id", "expansion_depth") if key in solution}
    if solution.get("gap_dir"):
        result["gap_dir"] = _relative(solution["gap_dir"], path)
    anchors = _values(solution.get("reachable_anchor_compounds", summary.get("reachable_anchor_compounds")))
    if anchors:
        result["reachable_anchor_compounds"] = anchors
    result.update(compound_names=names, steps=steps)
    if validation:
        result["validation"] = validation
    return result


def compact_solution(
    solution: Mapping[str, Any], *, manifest_path: str | Path,
    target_compound_id: str, validation_files: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Share the compact core while keeping RetroPath-specific meanings."""
    core = _compact_solution_core(solution, manifest_path=manifest_path,
        target_compound_id=target_compound_id, validation_files=validation_files)
    if _retropath(solution):
        from src.write_manifest.retropath_format import compact_retropath_solution
        return compact_retropath_solution(solution, core, Path(manifest_path))
    return core


def _compact(solution: Any) -> bool:
    return isinstance(solution, Mapping) and isinstance(solution.get("compound_names"), Mapping)


def manifest_for_storage(manifest: Mapping[str, Any], manifest_path: str | Path) -> dict[str, Any]:
    """Recompact only routes already using a new compact format."""
    result = deepcopy(dict(manifest))
    if _compact(result.get("solution")):
        result["solution"] = compact_solution(result["solution"], manifest_path=manifest_path,
            target_compound_id=str(result.get("target_compound_id") or ""))
    return result


def manifest_for_runtime(manifest: Mapping[str, Any], manifest_path: str | Path) -> dict[str, Any]:
    """Restore names, labels and computed summaries without rewriting files."""
    result = deepcopy(dict(manifest))
    solution = result.get("solution")
    if not _compact(solution):
        return result
    path = Path(manifest_path)
    names = solution["compound_names"]
    if solution.get("gap_dir"):
        solution["gap_dir"] = str(_resolve(solution["gap_dir"], path))
    steps = solution.get("steps", [])
    for step in steps:
        flags = _flags(step)
        step.update({key: key in flags for key in _FLAGS})
        step["requires_review"] = bool(flags)
        step["produced_compound_name"] = names.get(step.get("produced_compound_id"), "")
        precursors = _values(step.get("precursor_compound_ids"))
        step["precursor_compound_labels"] = "; ".join(
            f"{key} ({names[key]})" if names.get(key) else key for key in precursors)
        for key in _STEP_LISTS:
            step[key] = ";".join(_values(step.get(key)))
        step["reaction_comment"] = "; ".join(filter(None, (
            str(step.get("reaction_comment") or ""),
            *[_FLAG_LABELS[key] for key in flags],
        ))) or None
        step.setdefault("source_reaction_ids", step.get("reaction_id", ""))
        step.setdefault("resolution_action", "none")
    external = [step for step in steps if str(step.get("status") or "").lower() == "heterologous"]
    marked = [str(step["reaction_id"]) for step in steps if step.get("requires_review")]
    anchors = _values(solution.get("reachable_anchor_compounds"))
    solution["summary"] = {
        "target_compound_id": result.get("target_compound_id"),
        "target_compound_name": names.get(result.get("target_compound_id"), ""),
        "total_steps": len(steps), "heterologous_steps": len(external),
        "heterologous_reaction_ids": ";".join(str(step["reaction_id"]) for step in external),
        "heterologous_ko_ids": ";".join(dict.fromkeys(key for step in external for key in _values(step.get("ko_ids")))),
        "heterologous_enzyme_ecs": ";".join(dict.fromkeys(key for step in external for key in _values(step.get("enzyme_ecs")))),
        "reaction_resolution_status": "marked_for_review" if marked else "resolved",
        "requires_review": bool(marked), "marked_reaction_ids": ";".join(marked),
        "reachable_anchor_compounds": ";".join(anchors),
        "reachable_anchor_labels": "; ".join(f"{key} ({names[key]})" if names.get(key) else key for key in anchors),
    }
    validation = solution.get("validation")
    if isinstance(validation, dict):
        if validation.get("cofactor_mode"):
            validation["cofactor_relaxed"] = validation["cofactor_mode"] == "relaxed"
        references = [validation[key] for key in _REPORT_FIELDS if validation.get(key)]
        missing = [str(value) for value in references if not _resolve(value, path).is_file()]
        if references:
            validation["details_available"] = not missing
        if missing:
            message = "验证明细不可用 (validation details unavailable): " + "; ".join(missing)
            validation["details_issues"] = [message]
            warnings.warn(message, RuntimeWarning, stacklevel=2)
    if _retropath(solution):
        from src.write_manifest.retropath_format import restore_retropath_runtime
        restore_retropath_runtime(solution, path)
    return result
