"""Validate manually protected reaction groups before invalidating artifacts."""
from __future__ import annotations

from copy import deepcopy
import hashlib

from src.write_manifest.reaction_selection import (
    is_reaction_selection, read_protein_pair, retained_manual_selection,
)
from .models import PreselectedReaction


def scan_preselected_reactions(manifest: dict) -> dict:
    current = manifest.get("main_enzyme_selection")
    if not is_reaction_selection(current):
        return {}
    manual_groups = {reaction: group for reaction, group in current["reaction_assignments"].items()
                     if any(protein.get("selection_source") == "manual" for protein in group["proteins"])}
    if not manual_groups:
        return {}
    if current.get("selected_solution_id") != manifest.get("solution", {}).get("solution_id"):
        raise ValueError("人工已选反应对应的路线已变化，请核对 manifest")
    retained = retained_manual_selection(manifest, whole_reactions=True)["main_enzyme_selection"]
    snapshots, sequences = {}, {}
    for reaction, group in retained["reaction_assignments"].items():
        snapshot = {"step_indexes": list(group["step_indexes"]), "proteins": []}
        for protein in group["proteins"]:
            _, sequence = read_protein_pair(protein)
            accession = protein["accession"]
            if sequences.setdefault(accession, sequence) != sequence:
                raise ValueError(f"同一已选蛋白编号对应不同序列：{accession}")
            snapshot["proteins"].append({**deepcopy(protein),
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()})
        PreselectedReaction.model_validate(snapshot)
        snapshots[reaction] = snapshot
    return snapshots
