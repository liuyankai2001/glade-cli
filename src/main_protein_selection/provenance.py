"""Small deterministic digests used by the enzyme-selection pipeline."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def stable_json_hash(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def solution_fingerprint(solution_id: int | str, steps: list[dict[str, Any]]) -> str:
    projected = []
    contains_prediction = any(
        str(step.get("step_source") or "").strip() == "retropath"
        for step in steps
        if isinstance(step, dict)
    )
    for step in sorted(steps, key=lambda item: int(item.get("step_index") or 0)):
        if not isinstance(step, dict):
            continue
        row = {
            "step_index": int(step.get("step_index") or 0),
            "status": str(step.get("status") or ""),
            "reaction_id": str(step.get("reaction_id") or ""),
            "reaction_name": str(step.get("reaction_name") or ""),
            "reaction_comment": str(step.get("reaction_comment") or ""),
            "equation": str(step.get("equation") or ""),
            "direction": str(step.get("direction") or ""),
            "produced_compound_id": str(step.get("produced_compound_id") or ""),
            "produced_compound_name": str(step.get("produced_compound_name") or ""),
            "precursor_compound_ids": step.get("precursor_compound_ids"),
            "ko_ids": step.get("ko_ids"),
            "module_ids": step.get("module_ids"),
            "enzyme_ecs": step.get("enzyme_ecs"),
            "rhea_ids": step.get("rhea_ids"),
            "source_reaction_ids": step.get("source_reaction_ids"),
            "resolution_action": str(step.get("resolution_action") or ""),
            "resolution_evidence": step.get("resolution_evidence"),
        }
        if contains_prediction:
            row.update({
                "step_source": str(step.get("step_source") or ""),
                "retropath_step_id": str(step.get("retropath_step_id") or ""),
                "retropath_hypothesis_id": str(
                    step.get("retropath_hypothesis_id") or ""
                ),
                "retropath_rule_id": str(step.get("retropath_rule_id") or ""),
                "source_mnxr_id": str(step.get("source_mnxr_id") or ""),
                "source_ec_numbers": step.get("source_ec_numbers"),
                "source_uniprot_ids": step.get("source_uniprot_ids"),
                "source_rhea_ids": step.get("source_rhea_ids"),
                "exact_kegg_reaction_ids": step.get(
                    "exact_kegg_reaction_ids"
                ),
                "exact_rhea_ids": step.get("exact_rhea_ids"),
                "formal_mapping_exact": step.get("formal_mapping_exact"),
                "reaction_signature_sha256": str(
                    step.get("reaction_signature_sha256") or ""
                ),
                "full_reaction_smiles": str(
                    step.get("full_reaction_smiles") or ""
                ),
                "core_reaction_smiles": str(
                    step.get("core_reaction_smiles") or ""
                ),
                "stoichiometry_terms": step.get("stoichiometry_terms"),
                "prediction_provenance": step.get("prediction_provenance"),
            })
        projected.append(row)
    return stable_json_hash({"solution_id": int(solution_id), "steps": projected})


__all__ = ["file_sha256", "solution_fingerprint", "stable_json_hash"]
