"""Compact main-enzyme storage and an offline, legacy-compatible runtime view."""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any

from src.config.config import ROOT

FIELDS = ("selected_solution_id", "selected_set_id", "selection_status",
          "step_assignments", "review_required_step_indexes")
ASSIGNMENT_FIELDS = ("step_index", "reaction_id", "accession", "sequence_file")


def is_compact_selection(selection: Any) -> bool:
    return isinstance(selection, Mapping) and (
        selection.get("_compact_storage") is True
        or "schema_version" not in selection and all(key in selection for key in FIELDS)
    )


def compact_selection(selection: Mapping[str, Any]) -> dict[str, Any]:
    result = {key: deepcopy(selection[key]) for key in FIELDS}
    result["step_assignments"] = [
        {key: deepcopy(item[key]) for key in ASSIGNMENT_FIELDS}
        for item in selection["step_assignments"]
    ]
    return result


def selection_for_storage(manifest: dict[str, Any]) -> dict[str, Any]:
    selection = manifest.get("main_enzyme_selection")
    if is_compact_selection(selection):
        manifest["main_enzyme_selection"] = compact_selection(selection)
    return manifest


def resolve_sequence_file(value: str) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else ROOT / path).resolve()


def sequence_file_value(path: Path) -> str:
    resolved = path.expanduser().resolve()
    return resolved.relative_to(ROOT.resolve()).as_posix() if resolved.is_relative_to(ROOT.resolve()) else str(resolved)


def _bindings(assignments: list) -> list[tuple[int, str, str]]:
    return [(item["step_index"], item["reaction_id"], item["accession"]) for item in assignments]


def selected_sequences(directory: Path, result, selected) -> dict[str, str]:
    """Read the exact sequences used by the current candidate calculation."""
    from src.main_protein_selection.models import MainEnzymeSelectionResult
    from src.main_protein_selection.sequence_quality import analyze_protein_sequence

    path = directory / "main_enzyme_selection.json"
    if not path.is_file():
        raise ValueError("缺少当前候选序列，请重新运行 main-enzyme")
    candidates = MainEnzymeSelectionResult.model_validate(json.loads(path.read_text(encoding="utf-8-sig")))
    if (candidates.selected_solution_id != result.selected_solution_id
            or candidates.expansion_depth != result.expansion_depth
            or candidates.chassis_key != result.chassis_key
            or candidates.genome_ec_filter != result.genome_ec_filter):
        raise ValueError("候选结果与主酶组合不一致，请重新运行 main-enzyme、main-enzyme-sets")
    sequences = {}
    for assignment in selected.step_assignments:
        matches = [item for item in candidates.candidates_by_step.get(assignment.step_index, [])
                   if item.accession == assignment.accession and item.reaction_id == assignment.reaction_id]
        if len(matches) != 1:
            raise ValueError(f"找不到选中蛋白的唯一候选序列：{assignment.accession}；请重新运行 main-enzyme")
        quality = analyze_protein_sequence(matches[0].sequence)
        if not quality.normalized_sequence or not quality.is_constructable:
            raise ValueError(f"选中蛋白候选序列无效：{assignment.accession}；请重新运行 main-enzyme")
        previous = sequences.setdefault(assignment.accession, quality.normalized_sequence)
        if previous != quality.normalized_sequence:
            raise ValueError(f"同一蛋白在不同步骤的候选序列不一致：{assignment.accession}")
    return sequences


def _analysis_payload(manifest: Mapping[str, Any], compact: Mapping[str, Any], directory: Path):
    from src.main_protein_selection.common import get_solution_steps, heterologous_requirements
    from src.write_manifest.main_enzyme import _manifest_payload, _read_sets, _select_set

    result = _read_sets(directory / "main_enzyme_sets.json")
    selected = _select_set(result, compact["selected_set_id"])
    solution_id, steps = get_solution_steps(dict(manifest))
    depth = int(manifest["solution"].get("expansion_depth") or 0)
    required = {int(item["step_index"]): str(item.get("reaction_id") or "")
                for item in heterologous_requirements(steps)}
    external = set(result.genome_ec_filter.external_required_step_indexes) if result.genome_ec_filter else set(required)
    if (solution_id != compact["selected_solution_id"] or solution_id != result.selected_solution_id
            or depth != result.expansion_depth or set(result.required_step_indexes) != external
            or result.genome_ec_filter is not None and set(result.genome_ec_filter.required_step_indexes) != set(required)
            or any(required.get(item.step_index) != item.reaction_id for item in selected.step_assignments)
            or _bindings(compact["step_assignments"]) != _bindings([item.model_dump() for item in selected.step_assignments])):
        raise ValueError("路线或蛋白分配已变化，请重新运行 main-enzyme、main-enzyme-sets 和 write")
    sequences = selected_sequences(directory, result, selected)
    payload = _manifest_payload(result, selected, manifest["solution"].get("prediction"))
    reviews = sorted({step for item in payload["unresolved_reviews"] for step in item["step_indexes"]})
    if payload["selection_status"] != compact["selection_status"] or reviews != compact["review_required_step_indexes"]:
        raise ValueError("主酶复核信息已变化，请重新运行 main-enzyme-sets 和 write")
    return payload, sequences


def _validate_compact(compact: Mapping[str, Any]) -> None:
    for key in ("selected_solution_id", "selected_set_id"):
        if type(compact[key]) is not int or compact[key] < 1:
            raise ValueError(f"主酶清单字段无效：{key}")
    if compact["selection_status"] not in {"user_selected", "user_selected_pending_review"}:
        raise ValueError("主酶组合尚未选中")
    assignments = compact["step_assignments"]
    if not isinstance(assignments, list) or not assignments:
        raise ValueError("主酶步骤分配不能为空")
    steps = []
    paths = {}
    for item in assignments:
        if (not isinstance(item, Mapping) or type(item.get("step_index")) is not int or item["step_index"] < 1
                or any(not isinstance(item.get(key), str) or not item[key].strip()
                       for key in ("reaction_id", "accession", "sequence_file"))):
            raise ValueError("主酶步骤分配格式无效")
        steps.append(item["step_index"])
        path = resolve_sequence_file(item["sequence_file"])
        if path.name != f"{item['accession']}.fasta":
            raise ValueError("选中蛋白的序列文件名与编号不一致")
        previous = paths.setdefault(item["accession"], path)
        if previous != path:
            raise ValueError("同一蛋白绑定了不同序列文件")
    reviews = compact["review_required_step_indexes"]
    if steps != sorted(set(steps)) or not isinstance(reviews, list) or any(type(step) is not int for step in reviews):
        raise ValueError("主酶步骤编号格式无效")
    if reviews != sorted(set(reviews)) or not set(reviews).issubset(steps):
        raise ValueError("待复核步骤与主酶步骤不一致")


def restore_selection(manifest: dict[str, Any], manifest_path: Path) -> dict[str, Any]:
    """Missing evidence must not prevent reading a route or running recovery."""
    compact = manifest.get("main_enzyme_selection")
    if not is_compact_selection(compact):
        return manifest
    compact = compact_selection(compact)
    selection = {**deepcopy(compact), "_compact_storage": True, "proteins": [], "warnings": []}
    manifest["main_enzyme_selection"] = selection
    try:
        _validate_compact(compact)
    except (ValueError, TypeError, KeyError) as exc:
        selection["_compact_details_error"] = str(exc)
        return manifest

    errors = []
    payload, sequences = {}, {}
    try:
        payload, sequences = _analysis_payload(manifest, compact, manifest_path.resolve().parent / "main_protein_selection")
    except (OSError, ValueError, TypeError, KeyError) as exc:
        errors.append(f"主酶分析明细不可用，请重新运行 main-enzyme：{exc}")
    selection.update(payload)
    selection.update({key: value for key, value in compact.items() if key != "step_assignments"})
    full_assignments = {item["step_index"]: item for item in payload.get("step_assignments", [])}
    selection["step_assignments"] = [{**full_assignments.get(item["step_index"], {}), **item}
                                      for item in compact["step_assignments"]]

    from src.protein_selection.sequence_store import entry_sequence, read_cached_pair
    full_proteins = {item["accession"]: item for item in payload.get("proteins", [])}
    proteins = {}
    for assignment in compact["step_assignments"]:
        accession = assignment["accession"]
        if accession in proteins:
            proteins[accession]["assigned_step_indexes"].append(assignment["step_index"])
            continue
        path = resolve_sequence_file(assignment["sequence_file"])
        protein = {**full_proteins.get(accession, {}), "accession": accession,
                   "assigned_step_indexes": [assignment["step_index"]], "sequence_file": str(path)}
        if accession in sequences:
            protein["sequence_sha256"] = hashlib.sha256(sequences[accession].encode()).hexdigest()
        proteins[accession] = protein
        try:
            metadata = read_cached_pair(path.parent, accession)
            if metadata is None:
                raise ValueError(f"缺少完整的蛋白 FASTA/JSON：{accession}")
            for key in ("protein_name", "organism_name", "reviewed", "cofactors"):
                protein[key] = metadata[key]
            if accession in sequences and entry_sequence(metadata["uniprot_entry"]) != sequences[accession]:
                raise ValueError(f"本地蛋白序列与候选序列不一致：{accession}")
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    selection["proteins"] = [proteins[key] for key in sorted(proteins)]
    if errors:
        selection["_compact_details_error"] = "；".join(errors)
        selection["warnings"] = [*selection.get("warnings", []), *errors]
    return manifest
