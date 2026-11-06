"""Remove one reaction/protein association without deleting sequence files."""
from __future__ import annotations

from collections.abc import Mapping
from contextlib import ExitStack
from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any

from src.config.config import INPUTS_DIR
from src.main_protein_selection.candidate_inputs import read_candidate_inputs
from src.main_protein_selection.common import get_solution_steps
from src.main_protein_selection.lifecycle import serialized_main_enzyme_run
from src.protein_to_cds.artifacts import ArtifactTransaction, write_bytes_atomic
from src.write_manifest.main_enzyme_format import resolve_sequence_file, sequence_file_value
from src.write_manifest.reaction_selection import (
    _validate_structure, compact_reaction_selection, from_step_selection, is_reaction_selection,
    read_protein_pair, route_reactions,
)
from .candidate_sync import candidate_removal_updates, publish_candidate_updates
from .sequence_store import protein_file_lock


def _selection(raw, steps, folder):
    old = raw.get("main_enzyme_selection")
    if old is None:
        return None, False
    if not isinstance(old, Mapping):
        raise ValueError("已有主酶选择必须是对象")
    legacy = not is_reaction_selection(old)
    try:
        paths = {item["accession"]: sequence_file_value(folder / f"{item['accession']}.fasta")
                 for item in old.get("step_assignments", [])} if legacy else None
        result = from_step_selection(old, steps, paths=paths) if legacy else compact_reaction_selection(old)
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError("已有主酶选择格式无效，无法删除蛋白") from exc
    _validate_structure(result, raw)
    for group in result["reaction_assignments"].values():
        ids = [protein["accession"] for protein in group["proteins"]]
        if len(ids) != len(set(ids)):
            raise ValueError("同一反应下存在重复蛋白编号，请核对 manifest")
    return result, legacy


def remove_protein_from_reaction(config: Any, *, reaction_id: str, accession: str) -> dict:
    """Remove one selected association, offline, accepting IDs rather than paths."""
    if not isinstance(reaction_id, str) or not re.fullmatch(r"R\d{5}", reaction_id.strip().upper()):
        raise ValueError("反应编号必须为 R 加五位数字，例如 R02061")
    if not isinstance(accession, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", accession.strip()):
        raise ValueError("蛋白编号无效，只接受 UniProt 或本地编号，不接受路径")
    identifier = accession.strip().upper()
    return _remove(
        manifest_path=Path(config.manifest_output_path).expanduser().resolve(),
        output_dir=Path(config.project_output_path).expanduser().resolve() / "main_protein_selection",
        folder=Path(getattr(config, "inputs_dir", INPUTS_DIR)).expanduser().resolve() / "protein_sequences",
        reaction=reaction_id.strip().upper(), accession=identifier,
        genome_path=getattr(config, "genome_file", None), target_compound=getattr(config, "target_name", None))


@serialized_main_enzyme_run
def _remove(*, manifest_path, output_dir, folder, reaction, accession, genome_path, target_compound):
    if not manifest_path.is_file():
        raise ValueError("尚未选择路线，请先写入 solution")
    raw = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    if not isinstance(raw, dict) or raw.get("schema_version", "design_manifest.v1") != "design_manifest.v1":
        raise ValueError("manifest 格式无效，无法删除蛋白")
    if target_compound and raw.get("target_compound_id") != str(target_compound).strip().upper():
        raise ValueError("当前输入的目标化合物与 manifest 不一致")
    try:
        solution_id, steps = get_solution_steps(raw)
        route = route_reactions(steps)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("尚未选择有效路线，请先写入 solution") from exc
    if reaction not in route:
        raise ValueError(f"反应 {reaction} 不属于当前已选路线")
    selection, legacy = _selection(raw, steps, folder)
    group = selection["reaction_assignments"].get(reaction) if selection else None
    proteins = group["proteins"] if group else []
    changed = any(protein["accession"] == accession for protein in proteins)
    response = {"ok": True, "changed": changed, "selected_solution_id": solution_id, "reaction_id": reaction,
                "step_indexes": route[reaction], "accession": accession, "remaining_protein_count": len(proteins),
                "reaction_removed": False, "manifest_path": str(manifest_path), "candidates_updated": False,
                "combinations_cleared": False, "requires_candidate_regeneration": False}
    if not changed:
        return {**response, "message": "该反应未选择此蛋白，无变化"}
    before = {**deepcopy(raw), "main_enzyme_selection": deepcopy(selection)}
    group["proteins"] = [protein for protein in proteins if protein["accession"] != accession]
    response["remaining_protein_count"] = len(group["proteins"])
    if not group["proteins"]:
        del selection["reaction_assignments"][reaction]
        selection["review_required_step_indexes"] = sorted(
            set(selection["review_required_step_indexes"]) - set(route[reaction]))
        response["reaction_removed"] = True
    _validate_structure(selection, raw)
    after = {**deepcopy(raw), "main_enzyme_selection": selection}
    # Lock existing directories only: deleting an association must not need its lost files.
    locks = {(resolve_sequence_file(p["sequence_file"]).parent, p["accession"])
             for g in before["main_enzyme_selection"]["reaction_assignments"].values() for p in g["proteins"]
             if resolve_sequence_file(p["sequence_file"]).parent.is_dir()}
    with ExitStack() as stack:
        for directory, identifier in sorted(locks, key=lambda item: (str(item[0]), item[1])):
            stack.enter_context(protein_file_lock(directory, identifier))
        if legacy:
            sequences = {}
            for g in selection["reaction_assignments"].values():
                for protein in g["proteins"]:
                    _, sequence = read_protein_pair(protein)
                    if sequences.setdefault(protein["accession"], sequence) != sequence:
                        raise ValueError(f"同一蛋白编号对应不同序列：{protein['accession']}")
        updates = None
        was_fixed = any(p["selection_source"] == "manual" for p in proteins)
        still_fixed = any(p["selection_source"] == "manual" for p in group["proteins"])
        if was_fixed and not still_fixed:
            response["requires_candidate_regeneration"] = True
        elif still_fixed:
            updates = candidate_removal_updates(output_dir, before, after, genome_path)
            response["requires_candidate_regeneration"] = updates is None
        else:
            try:
                source, _, _, _ = read_candidate_inputs(output_dir / "main_enzyme_selection.json",
                                                       output_dir / "step_main_enzyme_candidates.csv", before, genome_path)
                response["requires_candidate_regeneration"] = not source.ok
            except (OSError, ValueError, KeyError, TypeError):
                response["requires_candidate_regeneration"] = True
        from src.main_protein_selection.build_main_enzyme_sets import main_enzyme_set_paths
        sets = tuple(main_enzyme_set_paths(output_dir).values())
        with ArtifactTransaction((manifest_path, *(updates or {}), *sets)) as transaction:
            if updates is not None:
                publish_candidate_updates(updates)
                response["candidates_updated"] = True
            for artifact in sets:
                if artifact.exists():
                    artifact.unlink()
                    response["combinations_cleared"] = True
            # Patch the raw section so legacy solution/downstream sections remain untouched.
            write_bytes_atomic(manifest_path, (json.dumps(after, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
            transaction.commit()
    response["message"] = ("已删除蛋白关联；请重新运行 main-enzyme 后再生成组合"
                           if response["requires_candidate_regeneration"] else "已删除蛋白关联，可重新生成主酶组合")
    return response
