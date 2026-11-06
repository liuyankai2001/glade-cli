"""Update a valid candidate snapshot after an offline manual selection."""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

from src.main_protein_selection.common import get_solution_steps, heterologous_requirements, merge_step_candidates
from src.main_protein_selection.genome_ec_filter import validate_genome_ec_filter
from src.main_protein_selection.models import MainEnzymeSelectionResult
from src.main_protein_selection.preselected_reactions import scan_preselected_reactions
from src.main_protein_selection.preselected_sets import validate_preselected_snapshot
from src.main_protein_selection.sequence_quality import analyze_protein_sequence
from src.protein_to_cds.artifacts import write_bytes_atomic


def _csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or "accession" not in fields:
            raise ValueError("候选 CSV 缺少必要字段")
        return fields, list(reader)


def _csv_bytes(fields, rows):
    text = io.StringIO(newline="")
    writer = csv.DictWriter(text, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value
                         for key, value in row.items()})
    return text.getvalue().encode("utf-8-sig")


def _without_fixed(value, indexes, reactions):
    if isinstance(value, list):
        return [_without_fixed(item, indexes, reactions) for item in value
                if not isinstance(item, dict) or not (
                    str(item.get("step_index", "")) in {str(index) for index in indexes}
                    or item.get("reaction_id") in reactions)]
    if isinstance(value, dict):
        return {key: _without_fixed(item, indexes, reactions) for key, item in value.items()}
    return value


def _source_failures(value, inherited=None):
    """Find attributable service failures without inventing a successful retry."""
    failures, unscoped = set(), False
    if isinstance(value, list):
        children = value
    elif isinstance(value, dict):
        index = value.get("step_index", inherited)
        if (value.get("source_status") == "source_unavailable" or value.get("status") == "source_unavailable"
                or value.get("source_unavailable") is True
                or value.get("reason_code") == "protein_evidence_source_unavailable"):
            if str(index).isdigit() and int(index) > 0:
                failures.add(int(index))
            else:
                unscoped = True
        inherited = index
        children = value.values()
    else:
        return failures, unscoped
    for child in children:
        child_failures, child_unscoped = _source_failures(child, inherited)
        failures.update(child_failures)
        unscoped = unscoped or child_unscoped
    return failures, unscoped


def candidate_updates(directory: Path, before: dict, after: dict, genome_path) -> dict[Path, bytes] | None:
    """Return a complete publication, or None for missing/stale candidate inputs."""
    source_path = directory / "main_enzyme_selection.json"
    try:
        source = MainEnzymeSelectionResult.model_validate(json.loads(source_path.read_text(encoding="utf-8-sig")))
        validate_preselected_snapshot(source, before, genome_path)
        fields, rows = _csv(directory / "step_main_enzyme_candidates.csv")
        def identity(row):
            return (int(row["step_index"]), row["reaction_id"], row["accession"], int(row["candidate_rank"]),
                    analyze_protein_sequence(row["sequence"]).normalized_sequence)
        expected = {(c.step_index, c.reaction_id, c.accession, c.candidate_rank,
                     analyze_protein_sequence(c.sequence).normalized_sequence)
                    for candidates in source.candidates_by_step.values() for c in candidates}
        if len(rows) != len(expected) or {identity(row) for row in rows} != expected:
            return None
        groups = scan_preselected_reactions(after)
        fixed = {index for group in groups.values() for index in group["step_indexes"]}
        _, steps = get_solution_steps(after)
        requirements = heterologous_requirements(steps)
        genome = validate_genome_ec_filter(source.genome_ec_filter, requirements, source.chassis_key, genome_path)
        external = set(genome.external_required_step_indexes) if genome else {s["step_index"] for s in requirements}
        pending = external - fixed
        evidence = {}
        failed, unscoped = set(), False
        for value in source.evidence_files.values():
            path = Path(value)
            if path.suffix != ".json" or path.resolve().parent != directory.resolve() or not path.exists() or path == source_path:
                continue
            data = json.loads(path.read_text(encoding="utf-8-sig"))
            evidence[path] = data
            indexes, unknown = _source_failures(data)
            failed.update(indexes)
            unscoped = unscoped or unknown
        payload = source.model_dump(mode="json")
        payload.update(schema_version="main_enzyme_selection.v8", preselected_reactions=groups,
                       candidates_by_step={index: candidates for index, candidates in source.candidates_by_step.items() if index not in fixed},
                       uncovered_step_indexes=sorted(pending - {index for index, candidates in source.candidates_by_step.items() if candidates}),
                       direction_rejected_step_indexes=sorted(set(source.direction_rejected_step_indexes) - fixed),
                       direction_risk_step_indexes=sorted(set(source.direction_risk_step_indexes) - fixed))
        if not pending or (source.status == "source_unavailable" and failed and failed.issubset(fixed) and not unscoped):
            payload.update(ok=True, status="complete")
        source = MainEnzymeSelectionResult.model_validate(payload)
        updates = {source_path: (json.dumps(source.model_dump(mode="json"), ensure_ascii=False, indent=2) + "\n").encode("utf-8")}
        remaining = [row for row in rows if int(row["step_index"]) not in fixed]
        updates[directory / "step_main_enzyme_candidates.csv"] = _csv_bytes(fields, remaining)
        audit = directory / "step_main_enzyme_candidate_audit.csv"
        if audit.exists():
            audit_fields, audit_rows = _csv(audit)
            updates[audit] = _csv_bytes(audit_fields, [r for r in audit_rows if int(r["step_index"]) not in fixed])
        merged = directory / "main_enzyme_candidates.csv"
        if merged.exists():
            merged_fields, _ = _csv(merged)
            updates[merged] = _csv_bytes(merged_fields, merge_step_candidates(remaining))
        for path, data in evidence.items():
            filtered = _without_fixed(data, fixed, set(groups))
            if filtered != data:
                updates[path] = (json.dumps(filtered, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
        return updates
    except (OSError, ValueError, KeyError, TypeError):
        return None


def publish_candidate_updates(updates: dict[Path, bytes]) -> None:
    for path, content in updates.items():
        if not path.exists() or path.read_bytes() != content:
            write_bytes_atomic(path, content)
