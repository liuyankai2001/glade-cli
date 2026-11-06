"""Reaction-oriented selected proteins, independent of a ranked set ID."""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any

from .main_enzyme_format import resolve_sequence_file

FIELDS = ("selected_solution_id", "reaction_assignments", "review_required_step_indexes")
PROTEIN_FIELDS = ("accession", "sequence_file", "selection_source")


def is_reaction_selection(value: Any) -> bool:
    return isinstance(value, Mapping) and "reaction_assignments" in value


def compact_reaction_selection(selection: Mapping[str, Any]) -> dict[str, Any]:
    try:
        return {
            "selected_solution_id": selection["selected_solution_id"],
            "reaction_assignments": {
                reaction: {"step_indexes": list(item["step_indexes"]), "proteins": [
                    {key: protein[key] for key in PROTEIN_FIELDS} for protein in item["proteins"]]}
                for reaction, item in selection["reaction_assignments"].items()},
            "review_required_step_indexes": list(selection["review_required_step_indexes"]),
        }
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError("按反应组织的主酶选择格式无效：reaction_assignments") from exc


def route_reactions(steps: list[dict]) -> dict[str, list[int]]:
    reactions = {}
    for step in steps:
        reactions.setdefault(str(step["reaction_id"]), []).append(int(step["step_index"]))
    return {reaction: sorted(set(indexes)) for reaction, indexes in reactions.items()}


def from_step_selection(selection: Mapping[str, Any], steps: list[dict], *, paths=None) -> dict[str, Any]:
    """Project both supported step formats; old selections were automatic."""
    reactions = route_reactions(steps)
    result = {"selected_solution_id": selection["selected_solution_id"], "reaction_assignments": {},
              "review_required_step_indexes": sorted(set(selection.get("review_required_step_indexes", [
                  step for item in selection.get("unresolved_reviews", []) for step in item["step_indexes"]])))}
    for assignment in sorted(selection["step_assignments"], key=lambda item: item["step_index"]):
        reaction, accession = assignment["reaction_id"], assignment["accession"]
        if reaction not in reactions or assignment["step_index"] not in reactions[reaction]:
            raise ValueError("主酶反应与当前路线不一致，请重新运行 main-enzyme")
        group = result["reaction_assignments"].setdefault(reaction, {"step_indexes": reactions[reaction], "proteins": []})
        if not any(protein["accession"] == accession for protein in group["proteins"]):
            group["proteins"].append({"accession": accession,
                "sequence_file": assignment.get("sequence_file") or (paths or {}).get(accession, ""),
                "selection_source": "automatic"})
    return result


def _validate_structure(selection: Mapping[str, Any], manifest: Mapping[str, Any], *, route_binding: bool = True) -> None:
    from src.main_protein_selection.common import get_solution_steps
    solution_id, steps = get_solution_steps(dict(manifest))
    if (type(selection.get("selected_solution_id")) is not int or selection["selected_solution_id"] < 1
            or route_binding and selection["selected_solution_id"] != solution_id):
        raise ValueError("主酶选择对应的路线已变化，请重新运行 main-enzyme")
    route = route_reactions(steps)
    assignments = selection.get("reaction_assignments")
    if not isinstance(assignments, Mapping):
        raise ValueError("reaction_assignments 必须是按反应编号组织的对象")
    assigned = set()
    for reaction, group in assignments.items():
        if not isinstance(group, Mapping) or not isinstance(reaction, str) or route_binding and reaction not in route:
            raise ValueError("主酶反应不属于当前路线，请重新运行 main-enzyme")
        indexes = group.get("step_indexes")
        if (not isinstance(indexes, list) or not indexes or any(type(index) is not int or index < 1 for index in indexes)
                or indexes != sorted(set(indexes)) or route_binding and indexes != route[reaction]):
            raise ValueError("反应步骤与当前路线不一致，请重新运行 main-enzyme")
        proteins = group.get("proteins")
        if not isinstance(proteins, list) or not proteins:
            raise ValueError("每个已选反应必须包含蛋白")
        for protein in proteins:
            if (not isinstance(protein, Mapping)
                    or any(not isinstance(protein.get(key), str) or not protein[key].strip() for key in PROTEIN_FIELDS)
                    or protein["selection_source"] not in {"automatic", "manual"}):
                raise ValueError("反应中的蛋白记录格式无效")
            if resolve_sequence_file(protein["sequence_file"]).name != f"{protein['accession']}.fasta":
                raise ValueError("选中蛋白序列文件名与编号不一致")
        assigned.update(indexes)
    reviews = selection.get("review_required_step_indexes")
    if (not isinstance(reviews, list) or any(type(index) is not int for index in reviews)
            or reviews != sorted(set(reviews)) or not set(reviews).issubset(assigned)):
        raise ValueError("待复核步骤与已选反应不一致")


def read_protein_pair(protein: Mapping[str, Any]):
    from src.protein_selection.sequence_store import entry_sequence, read_cached_pair
    path = resolve_sequence_file(protein["sequence_file"])
    if path.name != f"{protein['accession']}.fasta":
        raise ValueError("选中蛋白序列文件名与编号不一致")
    metadata = read_cached_pair(path.parent, protein["accession"])
    if metadata is None:
        raise ValueError(f"缺少完整的蛋白 FASTA/JSON：{protein['accession']}")
    return metadata, entry_sequence(metadata["uniprot_entry"])


def merge_manual_selection(automatic: dict[str, Any], current: Any, manifest: Mapping[str, Any]) -> dict[str, Any]:
    result = deepcopy(automatic)
    if is_reaction_selection(current) and current.get("selected_solution_id") == result["selected_solution_id"]:
        _validate_structure(current, manifest)
        for reaction, group in current["reaction_assignments"].items():
            manual = [deepcopy(protein) for protein in group["proteins"] if protein["selection_source"] == "manual"]
            if not manual:
                continue
            target = result["reaction_assignments"].setdefault(reaction, {"step_indexes": list(group["step_indexes"]), "proteins": []})
            manual_ids = {protein["accession"] for protein in manual}
            # Validate colliding automatic/manual records before preferring manual.
            sequences = {}
            for protein in target["proteins"] + manual:
                _, sequence = read_protein_pair(protein)
                if sequences.setdefault(protein["accession"], sequence) != sequence:
                    raise ValueError(f"同一蛋白编号对应不同序列：{protein['accession']}")
            target["proteins"] = [protein for protein in target["proteins"] if protein["accession"] not in manual_ids] + manual
            result["review_required_step_indexes"] = sorted(set(result["review_required_step_indexes"]) | set(group["step_indexes"]))
    # Validate globally as one accession also shares one CDS across reactions.
    sequences = {}
    for group in result["reaction_assignments"].values():
        unique = {}
        for protein in group["proteins"]:
            _, sequence = read_protein_pair(protein)
            accession = protein["accession"]
            if sequences.setdefault(accession, sequence) != sequence:
                raise ValueError(f"同一蛋白编号对应不同序列：{accession}")
            unique.setdefault(accession, protein)
        group["proteins"] = list(unique.values())
    _validate_structure(result, manifest)
    return result


def retained_manual_selection(manifest: Mapping[str, Any]) -> dict[str, Any]:
    current = manifest.get("main_enzyme_selection")
    if not is_reaction_selection(current):
        return {}
    if current.get("selected_solution_id") != manifest.get("solution", {}).get("solution_id"):
        return {}
    result = compact_reaction_selection(current)
    result["reaction_assignments"] = {
        reaction: {"step_indexes": group["step_indexes"], "proteins": [
            protein for protein in group["proteins"] if protein["selection_source"] == "manual"]}
        for reaction, group in result["reaction_assignments"].items()
        if any(protein["selection_source"] == "manual" for protein in group["proteins"])}
    if not result["reaction_assignments"]:
        return {}
    result["review_required_step_indexes"] = sorted({
        index for group in result["reaction_assignments"].values() for index in group["step_indexes"]})
    _validate_structure(result, manifest)
    return {"main_enzyme_selection": result}


def _candidate_context(manifest: Mapping[str, Any], manifest_path: Path):
    from src.main_protein_selection.common import get_solution_steps, heterologous_requirements
    from src.main_protein_selection.models import MainEnzymeSelectionResult
    path = manifest_path.resolve().parent / "main_protein_selection/main_enzyme_selection.json"
    if not path.is_file():
        raise ValueError("缺少当前底盘与候选上下文，请重新运行 main-enzyme")
    context = MainEnzymeSelectionResult.model_validate(json.loads(path.read_text(encoding="utf-8-sig")))
    solution_id, steps = get_solution_steps(dict(manifest))
    if (context.selected_solution_id != solution_id
            or context.expansion_depth != int(manifest["solution"].get("expansion_depth") or 0)
            or context.genome_ec_filter is not None and context.genome_ec_filter.required_step_indexes != sorted(
                int(step["step_index"]) for step in heterologous_requirements(steps))):
        raise ValueError("主酶上下文与当前路线不一致，请重新运行 main-enzyme")
    return context


def restore_reaction_selection(manifest: dict[str, Any], manifest_path: Path) -> dict[str, Any]:
    from src.main_protein_selection.common import get_solution_steps, heterologous_requirements
    from src.main_protein_selection.sequence_quality import analyze_protein_sequence
    from .main_enzyme import MAIN_ENZYME_MANIFEST_SCHEMA_VERSION

    raw = manifest["main_enzyme_selection"]
    compact = compact_reaction_selection(raw)
    selection = {**deepcopy(compact), "_reaction_storage": True,
                 "schema_version": MAIN_ENZYME_MANIFEST_SCHEMA_VERSION, "proteins": [], "warnings": []}
    manifest["main_enzyme_selection"] = selection
    errors = []
    if "schema_version" in raw and not raw.get("_reaction_storage"):
        errors.append("主酶选择结果格式已过期，请重新运行 main-enzyme 并写入")
    try:
        _validate_structure(compact, manifest)
    except (ValueError, TypeError, KeyError) as exc:
        errors.append(str(exc))
        try:
            _validate_structure(compact, manifest, route_binding=False)
        except (ValueError, TypeError, KeyError):
            selection["_compact_details_error"] = "；".join(errors)
            return manifest
    _, steps = get_solution_steps(manifest)
    selection["expansion_depth"] = int(manifest["solution"].get("expansion_depth") or 0)
    context = None
    try:
        context = _candidate_context(manifest, manifest_path)
        selection["chassis_key"] = context.chassis_key
        if context.genome_ec_filter is not None:
            selection["genome_ec_filter"] = context.genome_ec_filter.model_dump(mode="json")
    except (OSError, ValueError, TypeError, KeyError) as exc:
        selection["_selection_context_error"] = f"主酶上下文不可用，请重新运行 main-enzyme：{exc}"
        selection["warnings"].append(selection["_selection_context_error"])
    required = sorted(int(step["step_index"]) for step in heterologous_requirements(steps))
    supported = context.genome_ec_filter.supported_step_indexes if context and context.genome_ec_filter else []
    assigned = {index for group in compact["reaction_assignments"].values() for index in group["step_indexes"]}
    covered = sorted((assigned | set(supported)) & set(required))
    selection["coverage"] = {"complete": covered == required, "required_step_indexes": required,
        "covered_step_indexes": covered, "genome_supported_step_indexes": supported,
        "external_required_step_indexes": sorted(set(required) - set(supported))}
    reviews = compact["review_required_step_indexes"]
    selection["selection_status"] = "user_selected_pending_review" if reviews else "user_selected"
    selection["unresolved_reviews"] = ([{"review_type": "reaction_review", "status": "pending",
        "step_indexes": reviews, "reason": "部分已选蛋白仍需复核，详见主酶分析记录"}] if reviews else [])
    proteins, sequences = {}, {}
    selection["step_assignments"] = []
    for reaction, group in compact["reaction_assignments"].items():
        for selected in group["proteins"]:
            accession = selected["accession"]
            record = proteins.setdefault(accession, {"accession": accession,
                "assigned_step_indexes": [], "sequence_file": str(resolve_sequence_file(selected["sequence_file"]))})
            record["assigned_step_indexes"] = sorted(set(record["assigned_step_indexes"]) | set(group["step_indexes"]))
            selection["step_assignments"].extend({"step_index": index, "reaction_id": reaction, **selected}
                                                 for index in group["step_indexes"])
            try:
                metadata, sequence = read_protein_pair(selected)
                if sequences.setdefault(accession, sequence) != sequence:
                    raise ValueError(f"同一蛋白编号对应不同序列：{accession}")
                for key in ("protein_name", "organism_name", "reviewed", "cofactors"):
                    record[key] = metadata[key]
                record["sequence_sha256"] = hashlib.sha256(sequence.encode()).hexdigest()
                if selected["selection_source"] == "automatic" and context is not None:
                    matches = [candidate for index in group["step_indexes"] for candidate in context.candidates_by_step.get(index, [])
                               if candidate.accession == accession and candidate.reaction_id == reaction]
                    if not matches:
                        raise ValueError(f"当前候选中缺少已选蛋白：{accession}；请重新运行 main-enzyme")
                    for candidate in matches:
                        quality = analyze_protein_sequence(candidate.sequence)
                        if not quality.is_constructable or quality.normalized_sequence != sequence:
                            raise ValueError(f"本地蛋白序列与候选序列不一致：{accession}")
            except (OSError, ValueError) as exc:
                errors.append(str(exc))
    selection["proteins"] = [proteins[key] for key in sorted(proteins)]
    if errors:
        selection["_compact_details_error"] = "；".join(dict.fromkeys(errors))
        selection["warnings"].extend(errors)
    return manifest
