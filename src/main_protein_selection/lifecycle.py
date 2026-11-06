"""Invalidate current calculation references while retaining uploaded drafts."""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from functools import wraps
from inspect import signature

from filelock import FileLock, Timeout
from typing import Any



def retained_uploaded_draft(manifest: Mapping[str, Any]) -> dict[str, Any]:
    raw = manifest.get("expression_parts_draft")
    if not isinstance(raw, Mapping):
        return {}
    if raw.get("schema_version") != "expression_parts_draft.v1" and raw.get("source_type") != "user_uploaded":
        return {}
    from src.expression_box.parts_preparation import normalize_parts_draft, stale_parts_designs
    draft = normalize_parts_draft(raw)
    if draft.get("source_type") != "user_uploaded":
        return {}
    stale_parts_designs(draft["designs"])
    for design in draft["designs"]:
        for cassette in design["cassettes"]:
            cassette.pop("optimized_cds_by_accession", None)
    return {"expression_parts_draft": draft}


def invalidate_current(path: str | Path, sections: tuple[str, ...]) -> dict[str, Any]:
    from src.write_manifest.store import read_design_manifest, update_design_manifest
    manifest = read_design_manifest(path)
    if not any(name in manifest for name in sections):
        return manifest
    retained = retained_uploaded_draft(manifest)
    if "main_enzyme_selection" in sections:
        from src.write_manifest.reaction_selection import retained_manual_selection
        retained.update(retained_manual_selection(manifest))
    return update_design_manifest(
        path, target_compound_id=str(manifest.get("target_compound_id") or ""),
        sections=retained, discard_sections=sections,
    )


def serialized_main_enzyme_run(function):
    """Keep candidate and set publication mutually exclusive with manifest writes."""
    @wraps(function)
    def wrapped(*args, **kwargs):
        arguments = signature(function).bind(*args, **kwargs).arguments
        output = Path(arguments["output_dir"]).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        from src.write_manifest.store import manifest_update_lock
        try:
            with FileLock(str(output / ".main_enzyme_calculation.lock"), timeout=0):
                with manifest_update_lock(arguments["manifest_path"]):
                    return function(*args, **kwargs)
        except Timeout as exc:
            raise ValueError("主酶计算正在运行，请稍后重试") from exc
    return wrapped
