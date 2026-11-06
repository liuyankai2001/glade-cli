"""Read and validate candidate artifacts before changing selected proteins."""
from __future__ import annotations

import csv
import json
from collections.abc import Mapping
from pathlib import Path

from .models import MainEnzymeSelectionResult
from .preselected_sets import validate_preselected_snapshot
from .sequence_quality import analyze_protein_sequence
from src.write_manifest.main_enzyme_format import resolve_sequence_file


def read_candidate_inputs(selection_path: Path, csv_path: Path, manifest: dict, genome_path):
    """Check the route, fixed snapshot, genome and canonical JSON/CSV identities."""
    source = MainEnzymeSelectionResult.model_validate(
        json.loads(selection_path.read_text(encoding="utf-8-sig")))
    selection = manifest.get("main_enzyme_selection")
    groups = selection.get("reaction_assignments", {}) if isinstance(selection, Mapping) else {}
    paths = [p.sequence_file for group in source.preselected_reactions.values() for p in group.proteins]
    paths.extend(p["sequence_file"] for group in groups.values()
                 if any(p.get("selection_source") == "manual" for p in group["proteins"]) for p in group["proteins"])
    for value in paths:
        path = resolve_sequence_file(value)
        if not path.is_file() or not path.with_suffix(".json").is_file():
            raise ValueError("缺少完整的已选蛋白 FASTA/JSON，请重新运行 main-enzyme")
    records = validate_preselected_snapshot(source, manifest, genome_path)
    with csv_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or not {"step_index", "reaction_id", "accession", "candidate_rank", "sequence", "selection_status"}.issubset(fields):
            raise ValueError("候选 CSV 缺少必要字段，请重新运行 main-enzyme")
        rows = list(reader)

    def identity(row):
        return (int(row["step_index"]), row["reaction_id"], row["accession"], int(row["candidate_rank"]),
                analyze_protein_sequence(row["sequence"]).normalized_sequence)

    expected = {(c.step_index, c.reaction_id, c.accession, c.candidate_rank,
                 analyze_protein_sequence(c.sequence).normalized_sequence)
                for candidates in source.candidates_by_step.values() for c in candidates}
    if (len(rows) != len(expected) or {identity(row) for row in rows} != expected
            or any(row.get("selection_status", "selected") != "selected" for row in rows)):
        raise ValueError("候选 JSON 与 CSV 不一致，请重新运行 main-enzyme")
    return source, fields, rows, records
