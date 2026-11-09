"""Offline search and publication with complete preselected reaction groups."""
from __future__ import annotations

import hashlib
import json
from .transmembrane import annotation

from .models import MainEnzymeSelectionResult, PreselectedReaction
from .preselected_set_models import PreselectedEnzymeSetsResult, metadata_metrics
from .preselected_reactions import scan_preselected_reactions
from .genome_ec_filter import validate_genome_ec_filter
from .common import get_solution_steps, heterologous_requirements, read_manifest, read_csv, rel_or_abs
from src.write_manifest.reaction_selection import read_protein_pair
from src.write_manifest.main_enzyme_format import resolve_sequence_file


def validate_preselected_snapshot(source, manifest, genome_path):
    solution_id, steps = get_solution_steps(manifest)
    if source.selected_solution_id != solution_id or source.expansion_depth != int(manifest["solution"].get("expansion_depth") or 0):
        raise ValueError("路线已变化，请重新运行 main-enzyme")
    current = scan_preselected_reactions(manifest)
    def signature(groups):
        return {reaction: (tuple(group["step_indexes"]), [
            (protein["accession"], resolve_sequence_file(protein["sequence_file"]), protein["selection_source"], protein["sequence_sha256"])
            for protein in group["proteins"]]) for reaction, group in groups.items()}
    saved = {reaction: group.model_dump(mode="json") for reaction, group in source.preselected_reactions.items()}
    if signature(saved) != signature(current):
        raise ValueError("已选反应快照已变化，请重新运行 main-enzyme")
    requirements = heterologous_requirements(steps)
    live_filter = validate_genome_ec_filter(source.genome_ec_filter, requirements, source.chassis_key, genome_path)
    route_reactions = {int(step["step_index"]): str(step.get("reaction_id") or "") for step in steps}
    for index, candidates in getattr(source, "candidates_by_step", {}).items():
        if any(candidate.reaction_id != route_reactions.get(index) for candidate in candidates):
            raise ValueError("候选反应与当前路线不一致，请重新运行 main-enzyme")
    fixed_indexes = {index for group in source.preselected_reactions.values() for index in group.step_indexes}
    if hasattr(source, "required_step_indexes"):
        external = (set(live_filter.external_required_step_indexes) if live_filter is not None else
                    {int(step["step_index"]) for step in requirements})
        if set(source.required_step_indexes) != external | fixed_indexes:
            raise ValueError("组合所需步骤与当前路线不一致，请重新运行 main-enzyme 并生成组合")
        for selected in source.sets:
            if any(item.reaction_id != route_reactions.get(item.step_index) for item in selected.step_assignments):
                raise ValueError("组合反应与当前路线不一致，请重新运行 main-enzyme")
    records = {}
    for group in saved.values():
        for protein in group["proteins"]:
            metadata, sequence = read_protein_pair(protein)
            if hashlib.sha256(sequence.encode()).hexdigest() != protein["sequence_sha256"]:
                raise ValueError(f"已选蛋白序列快照已变化：{protein['accession']}")
            if protein["accession"] in records and records[protein["accession"]]["sequence"] != sequence:
                raise ValueError("同编号固定蛋白序列冲突")
            records.setdefault(protein["accession"], {**metadata, "sequence": sequence,
                                                      "transmembrane": annotation(protein.get("transmembrane"))})
    return records


def build_preselected_sets_from_rows(*, solution_id, expansion_depth, chassis_key, required_steps, candidate_rows,
                                    max_sets, max_search_nodes, source_artifacts, genome_ec_filter,
                                    preselected_reactions, preselected_proteins, preselected_review_step_indexes):
    """Pure v8 search over fixed, validated metadata and automatic shortlist rows."""
    import importlib
    engine = importlib.import_module("src.main_protein_selection.build_main_enzyme_sets")
    if solution_id < 1 or expansion_depth < 0 or max_sets < 1 or max_search_nodes < 1:
        raise ValueError("invalid preselected search parameters")
    groups = {reaction: PreselectedReaction.model_validate(group).model_dump(mode="json")
              for reaction, group in preselected_reactions.items()}
    original = engine._canonical_required_steps(required_steps)
    original_map = {step["step_index"]: step for step in original}
    if genome_ec_filter is not None:
        if sorted(original_map) != genome_ec_filter.required_step_indexes or chassis_key != genome_ec_filter.chassis_key:
            raise ValueError("Genome filter does not match the original requirements")
        external = set(genome_ec_filter.external_required_step_indexes)
    else:
        external = set(original_map)
    fixed_steps, fixed_assignments = set(), {}
    for reaction, group in groups.items():
        for index in group["step_indexes"]:
            if index in fixed_steps or index in original_map and original_map[index]["reaction_id"] != reaction:
                raise ValueError("fixed reactions overlap or differ from route steps")
            fixed_steps.add(index)
        for protein in group["proteins"]:
            accession = protein["accession"]
            record = preselected_proteins[accession]
            if hashlib.sha256(record["sequence"].encode()).hexdigest() != protein["sequence_sha256"]:
                raise ValueError("fixed sequence differs from its snapshot")
            fixed_assignments.setdefault(accession, set()).update(group["step_indexes"])
    reviews = sorted(set(preselected_review_step_indexes))
    if not set(reviews).issubset(fixed_steps):
        raise ValueError("fixed review indexes exceed fixed steps")
    required = sorted(external | fixed_steps)
    automatic = sorted(external - fixed_steps)
    base = dict(selected_solution_id=solution_id, expansion_depth=expansion_depth, chassis_key=chassis_key,
        parameters={"max_sets": max_sets, "max_search_nodes": max_search_nodes, "candidate_scope": "top_n_shortlist"},
        required_step_indexes=required, preselected_reactions=groups, preselected_review_step_indexes=reviews,
        genome_ec_filter=genome_ec_filter, source_artifacts=dict(source_artifacts or {}))
    def result(**kwargs):
        return PreselectedEnzymeSetsResult(**base, **kwargs).model_dump(mode="json")
    if not automatic:
        return result(ok=True, status="already_selected", minimum_protein_count=len(fixed_assignments),
            search_complete=True, search_nodes=0, sets=[], uncovered_step_indexes=[])
    pending = [original_map[index] for index in automatic]
    filtered, conflicts = [], []
    for row in candidate_rows:
        if engine._int(row.get("step_index")) not in automatic:
            continue
        accession = engine._text(row.get("accession")).upper()
        sequence = engine.analyze_protein_sequence(row.get("sequence")).normalized_sequence
        if accession in fixed_assignments and sequence != preselected_proteins[accession]["sequence"]:
            conflicts.append(f"{accession}: 自动候选序列与固定蛋白序列冲突，已排除该候选")
        else:
            filtered.append(row)
    pool, rejections = engine._build_protein_pool(pending, filtered)
    pool_steps = {index for protein in pool for index in protein.capable_steps}
    if set(automatic) - pool_steps:
        return result(ok=False, status="infeasible", search_complete=True, search_nodes=0, sets=[],
            minimum_protein_count=None, uncovered_step_indexes=automatic, warnings=sorted(set(conflicts)),
            blocking_reasons=["No constructable, non-conflicting Top-N candidate for Step " +
                ", ".join(map(str, sorted(set(automatic) - pool_steps))), *sorted(rejections)])
    by_accession = {protein.accession: protein for protein in pool}
    requirements = {step["step_index"]: step for step in pending}
    def payload(assignment):
        selected = engine._build_set_payload(assignment, by_accession, requirements)
        proteins = {protein["accession"]: protein for protein in selected["proteins"]}
        for accession, indexes in fixed_assignments.items():
            record = preselected_proteins[accession]
            protein = proteins.setdefault(accession, {"accession": accession, "assigned_step_indexes": [], "capable_step_indexes": []})
            protein.update({key: record[key] for key in ("protein_name", "organism_name", "reviewed", "cofactors")})
            protein["transmembrane"] = annotation(record.get("transmembrane"))
            protein["assigned_step_indexes"] = sorted(set(protein["assigned_step_indexes"]) | indexes)
            protein["capable_step_indexes"] = sorted(set(protein["capable_step_indexes"]) | indexes)
        selected["proteins"] = [proteins[accession] for accession in sorted(proteins)]
        selected["protein_count"] = len(proteins)
        selected["covered_step_indexes"] = required
        selected["review_required_step_indexes"] = sorted(set(selected["review_required_step_indexes"]) | set(reviews))
        selected["status"] = "review_required" if selected["review_required_step_indexes"] or selected["warnings"] else "complete"
        selected["metrics"].update(metadata_metrics(selected["proteins"]))
        selected["reasons"] = ["covers every required heterologous step", f"uses {len(proteins)} distinct protein(s)",
                               *selected["reasons"][2:], "includes every preselected protein"]
        return selected
    search = engine._enumerate_valid_assignments(automatic, pool, max_search_nodes=max_search_nodes,
        fixed_accessions=frozenset(fixed_assignments),
        assignment_quality=lambda assignment: (*engine._set_ranking_key(payload(assignment)), assignment))
    selected_sets = sorted((payload(assignment) for assignment in search.assignments), key=engine._set_ranking_key)[:max_sets]
    for index, selected in enumerate(selected_sets, 1):
        selected["set_id"] = index
    warnings = sorted(set(conflicts))
    if not search.complete:
        status = "truncated"
        warnings.append(f"Search stopped at max_search_nodes={max_search_nodes}; reported sets may be incomplete.")
    elif not selected_sets:
        status = "infeasible"
    else:
        status = selected_sets[0]["status"]
    return result(ok=bool(selected_sets), status=status, minimum_protein_count=search.minimum_count if selected_sets else None,
        search_complete=search.complete, search_nodes=search.nodes, sets=selected_sets,
        uncovered_step_indexes=[] if selected_sets else automatic, warnings=warnings,
        blocking_reasons=[] if selected_sets else ["No complete combination was found in the selected Top-N shortlist."])


def build_persisted_preselected_sets(manifest_path, selection_path, candidate_csv_path, output_dir, *, max_sets, max_search_nodes, genome_path):
    import importlib
    engine = importlib.import_module("src.main_protein_selection.build_main_enzyme_sets")
    selection = MainEnzymeSelectionResult.model_validate(json.loads(selection_path.read_text(encoding="utf-8-sig")))
    manifest = read_manifest(manifest_path)
    records = validate_preselected_snapshot(selection, manifest, genome_path)
    rows = read_csv(candidate_csv_path)
    fixed = {index for group in selection.preselected_reactions.values() for index in group.step_indexes}
    result = build_preselected_sets_from_rows(solution_id=selection.selected_solution_id, expansion_depth=selection.expansion_depth,
        chassis_key=selection.chassis_key, required_steps=heterologous_requirements(manifest["solution"]["steps"]),
        candidate_rows=rows, max_sets=max_sets, max_search_nodes=max_search_nodes,
        source_artifacts={"manifest": rel_or_abs(manifest_path), "candidate_selection": rel_or_abs(selection_path),
                          "candidate_details": rel_or_abs(candidate_csv_path)},
        genome_ec_filter=selection.genome_ec_filter, preselected_reactions=selection.preselected_reactions,
        preselected_proteins=records, preselected_review_step_indexes=sorted(set(manifest["main_enzyme_selection"]["review_required_step_indexes"]) & fixed))
    if selection.status == "source_unavailable":
        automatic = sorted(set(result["required_step_indexes"]) - fixed)
        result.update(ok=False, status="source_unavailable", sets=[], minimum_protein_count=None,
                      search_complete=False, search_nodes=0, uncovered_step_indexes=automatic,
                      blocking_reasons=["Main-enzyme candidate retrieval was unavailable; rerun main-enzyme."])
        result = PreselectedEnzymeSetsResult.model_validate(result).model_dump(mode="json")
    validate_preselected_snapshot(selection, read_manifest(manifest_path), genome_path)
    from src.protein_to_cds.artifacts import ArtifactTransaction
    with ArtifactTransaction((manifest_path, *engine.main_enzyme_set_paths(output_dir).values())) as transaction:
        if result["status"] != "already_selected":
            from .lifecycle import invalidate_current
            from src.write_manifest.main_enzyme import MAIN_ENZYME_DOWNSTREAM_SECTIONS
            invalidate_current(manifest_path, ("main_enzyme_selection", *MAIN_ENZYME_DOWNSTREAM_SECTIONS), preserve_manual_reactions=True)
        result["output_files"] = engine._write_set_outputs(output_dir, result)
        transaction.commit()
    return result
