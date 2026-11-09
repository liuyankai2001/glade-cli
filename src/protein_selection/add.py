"""Append an offline selected protein to one reaction of the current route."""
from __future__ import annotations

from contextlib import ExitStack
from collections.abc import Mapping
from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any

from src.config.config import INPUTS_DIR
from src.main_protein_selection.common import get_solution_steps
from src.main_protein_selection.transmembrane import extract_annotation
from src.main_protein_selection.lifecycle import serialized_main_enzyme_run
from src.protein_to_cds.artifacts import ArtifactTransaction, write_bytes_atomic
from src.write_manifest.main_enzyme_format import resolve_sequence_file, sequence_file_value
from src.write_manifest.reaction_selection import (
    _validate_structure, compact_reaction_selection, from_step_selection, is_reaction_selection,
    read_protein_pair, route_reactions,
)
from src.write_manifest.store import read_design_manifest, update_design_manifest
from .candidate_sync import candidate_updates, publish_candidate_updates
from .local_sequences import local_filename, read_sequence_pair, sequence_identity
from .sequence_store import protein_file_lock


def _selection(manifest, steps, folder):
    old = manifest.get("main_enzyme_selection")
    if old is None:
        return {"selected_solution_id": manifest["solution"]["solution_id"], "reaction_assignments": {},
                "review_required_step_indexes": []}
    if not isinstance(old, Mapping):
        raise ValueError("已有主酶选择必须是对象")
    try:
        if is_reaction_selection(old):
            result = compact_reaction_selection(old)
        else:
            paths = {item["accession"]: sequence_file_value(folder / f"{item['accession']}.fasta")
                     for item in old.get("step_assignments", [])}
            result = from_step_selection(old, steps, paths=paths)
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError("已有主酶选择格式无效，无法恢复本地蛋白") from exc
    _validate_structure(result, manifest)
    for group in result["reaction_assignments"].values():
        identifiers = [protein["accession"] for protein in group["proteins"]]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("同一反应下存在重复蛋白编号，请核对 manifest")
    return result


def add_protein_to_reaction(config: Any, *, reaction_id: str, sequence_file: str) -> dict:
    """Use a basename from inputs/protein_sequences, without querying any service."""
    name = local_filename(sequence_file)
    if not isinstance(reaction_id, str) or not re.fullmatch(r"R\d{5}", reaction_id.strip().upper()):
        raise ValueError("反应编号必须为 R 加五位数字，例如 R02061")
    reaction = reaction_id.strip().upper()
    folder = Path(getattr(config, "inputs_dir", INPUTS_DIR)).expanduser().resolve() / "protein_sequences"
    path = folder / name
    if path.resolve().parent != folder.resolve():
        raise ValueError("蛋白序列文件必须位于 inputs/protein_sequences 中")
    return _add(manifest_path=Path(config.manifest_output_path).expanduser().resolve(),
                output_dir=Path(config.project_output_path).expanduser().resolve() / "main_protein_selection",
                reaction=reaction, path=path, genome_path=getattr(config, "genome_file", None),
                target_compound=getattr(config, "target_name", None))


@serialized_main_enzyme_run
def _add(*, manifest_path, output_dir, reaction, path, genome_path, target_compound):
    manifest = read_design_manifest(manifest_path)
    if target_compound and manifest.get("target_compound_id") != str(target_compound).strip().upper():
        raise ValueError("当前输入的目标化合物与 manifest 不一致")
    solution_id, steps = get_solution_steps(manifest)
    route = route_reactions(steps)
    if reaction not in route:
        raise ValueError(f"反应 {reaction} 不属于当前已选路线")
    selection = _selection(manifest, steps, path.parent)
    accession = sequence_identity(path)
    selected = [p for group in selection["reaction_assignments"].values() for p in group["proteins"]]
    locks = {(path.parent.resolve(), accession)} | {
        (resolve_sequence_file(p["sequence_file"]).parent, p["accession"]) for p in selected}
    with ExitStack() as stack:
        for folder, identifier in sorted(locks, key=lambda item: (str(item[0]), item[1])):
            stack.enter_context(protein_file_lock(folder, identifier))
        metadata, sequence = read_sequence_pair(path, accession, allow_missing_metadata=True)
        canonical = sequence_file_value(path)
        sequences = {}
        existing_paths = {}
        for protein in selected:
            _, existing = read_protein_pair(protein)
            if sequences.setdefault(protein["accession"], existing) != existing:
                raise ValueError(f"同一蛋白编号对应不同序列：{protein['accession']}")
            existing_paths.setdefault(protein["accession"], protein["sequence_file"])
        if accession in sequences:
            if sequences[accession] != sequence:
                raise ValueError(f"同一蛋白编号对应不同序列：{accession}")
            canonical = existing_paths[accession]
        target = selection["reaction_assignments"].setdefault(reaction, {"step_indexes": route[reaction], "proteins": []})
        previous = next((p for p in target["proteins"] if p["accession"] == accession), None)
        if previous is not None:
            canonical = previous["sequence_file"]
        changed = previous is None or previous["selection_source"] != "manual"
        response = {"ok": True, "changed": changed, "selected_solution_id": solution_id, "reaction_id": reaction,
                    "step_indexes": route[reaction], "accession": accession, "sequence_file": canonical,
                    "manifest_path": str(manifest_path), "candidates_updated": False, "combinations_cleared": False}
        if not changed:
            return response
        if previous is None:
            target["proteins"].append({"accession": accession, "sequence_file": canonical, "selection_source": "manual",
                                       "transmembrane": extract_annotation(metadata.get("uniprot_entry"))})
        else:
            previous["selection_source"] = "manual"
        selection["review_required_step_indexes"] = sorted(set(selection["review_required_step_indexes"]) | set(route[reaction]))
        _validate_structure(selection, manifest)
        from src.main_protein_selection.build_main_enzyme_sets import main_enzyme_set_paths
        sets = tuple(main_enzyme_set_paths(output_dir).values())
        sidecar = path.with_suffix(".json")
        files = [manifest_path, sidecar, *sets, *output_dir.glob("*.json"), *output_dir.glob("*.csv")]
        original_fasta = path.read_bytes()
        with ArtifactTransaction(dict.fromkeys(files)) as transaction:
            if not sidecar.exists():
                write_bytes_atomic(sidecar, (json.dumps(metadata, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
            after = {**deepcopy(manifest), "main_enzyme_selection": selection}
            updates = candidate_updates(output_dir, manifest, after, genome_path)
            if updates is not None:
                publish_candidate_updates(updates)
                response["candidates_updated"] = True
            for artifact in sets:
                if artifact.exists():
                    artifact.unlink()
                    response["combinations_cleared"] = True
            if path.read_bytes() != original_fasta:
                raise ValueError("添加期间蛋白序列文件已变化，请重试")
            # Repeat pair validation immediately before the manifest commit.
            read_sequence_pair(path, accession)
            for protein in selected:
                if read_protein_pair(protein)[1] != sequences[protein["accession"]]:
                    raise ValueError("添加期间已有蛋白序列已变化，请重试")
            update_design_manifest(manifest_path, target_compound_id=str(manifest.get("target_compound_id") or ""),
                                   sections={"main_enzyme_selection": selection})
            transaction.commit()
        response["message"] = ("已添加人工蛋白，可继续生成主酶组合" if response["candidates_updated"] else
                               "已添加人工蛋白；请运行 main-enzyme 生成当前候选后再生成组合")
        return response
