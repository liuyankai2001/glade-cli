"""Replace one selected reaction/protein association in a single transaction."""
from __future__ import annotations

from contextlib import ExitStack
from copy import deepcopy
import json
from pathlib import Path
import re
from typing import Any

from src.config.config import INPUTS_DIR
from src.main_protein_selection.common import get_solution_steps
from src.main_protein_selection.lifecycle import serialized_main_enzyme_run
from src.protein_to_cds.artifacts import ArtifactTransaction, write_bytes_atomic
from src.write_manifest.main_enzyme_format import resolve_sequence_file, sequence_file_value
from src.write_manifest.reaction_selection import _validate_structure, read_protein_pair, route_reactions
from .candidate_sync import candidate_updates, candidate_removal_updates, publish_candidate_updates
from .local_sequences import local_filename, read_sequence_pair, sequence_identity
from .remove import _selection
from .sequence_store import protein_file_lock


def replace_protein_in_reaction(config: Any, *, reaction_id: str, old_accession: str, sequence_file: str) -> dict:
    """Replace an existing ID with a FASTA basename from inputs/protein_sequences."""
    name = local_filename(sequence_file)
    if not isinstance(reaction_id, str) or not re.fullmatch(r"R\d{5}", reaction_id.strip().upper()):
        raise ValueError("反应编号必须为 R 加五位数字，例如 R02061")
    if not isinstance(old_accession, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", old_accession.strip()):
        raise ValueError("旧蛋白编号无效，只接受 UniProt 或本地编号，不接受路径")
    folder = Path(getattr(config, "inputs_dir", INPUTS_DIR)).expanduser().resolve() / "protein_sequences"
    path = folder / name
    if path.resolve().parent != folder.resolve():
        raise ValueError("蛋白序列文件必须位于 inputs/protein_sequences 中")
    return _replace(
        manifest_path=Path(config.manifest_output_path).expanduser().resolve(),
        output_dir=Path(config.project_output_path).expanduser().resolve() / "main_protein_selection",
        reaction=reaction_id.strip().upper(), old_accession=old_accession.strip().upper(), path=path,
        genome_path=getattr(config, "genome_file", None), target_compound=getattr(config, "target_name", None))


@serialized_main_enzyme_run
def _replace(*, manifest_path, output_dir, reaction, old_accession, path, genome_path, target_compound):
    if not manifest_path.is_file():
        raise ValueError("尚未选择路线，请先写入 solution")
    raw = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    if not isinstance(raw, dict) or raw.get("schema_version", "design_manifest.v1") != "design_manifest.v1":
        raise ValueError("manifest 格式无效，无法替换蛋白")
    if target_compound and raw.get("target_compound_id") != str(target_compound).strip().upper():
        raise ValueError("当前输入的目标化合物与 manifest 不一致")
    try:
        solution_id, steps = get_solution_steps(raw)
        route = route_reactions(steps)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("尚未选择有效路线，请先写入 solution") from exc
    if reaction not in route:
        raise ValueError(f"反应 {reaction} 不属于当前已选路线")
    selection, _ = _selection(raw, steps, path.parent, action="替换")
    group = selection["reaction_assignments"].get(reaction) if selection else None
    old = next((p for p in group["proteins"] if p["accession"] == old_accession), None) if group else None
    if old is None:
        raise ValueError(f"反应 {reaction} 未选择旧蛋白 {old_accession}，无法替换")
    if not path.is_file():
        raise ValueError(f"缺少蛋白 FASTA：{path.name}")
    accession = sequence_identity(path)
    before = {**deepcopy(raw), "main_enzyme_selection": deepcopy(selection)}
    selected = [p for g in selection["reaction_assignments"].values() for p in g["proteins"]]
    # A removed, different ID need not have readable files. Every survivor must.
    retained = [p for p in selected if p is not old or accession == old_accession]
    locks = {(path.parent.resolve(), accession)} | {
        (resolve_sequence_file(p["sequence_file"]).parent, p["accession"]) for p in selected
        if resolve_sequence_file(p["sequence_file"]).parent.is_dir()}
    with ExitStack() as stack:
        for directory, identifier in sorted(locks, key=lambda item: (str(item[0]), item[1])):
            stack.enter_context(protein_file_lock(directory, identifier))
        metadata, sequence = read_sequence_pair(path, accession, allow_missing_metadata=True)
        original_fasta = path.read_bytes()
        sequences, existing_paths = {}, {}
        for protein in retained:
            _, existing = read_protein_pair(protein)
            if sequences.setdefault(protein["accession"], existing) != existing:
                raise ValueError(f"同一蛋白编号对应不同序列：{protein['accession']}")
            existing_paths.setdefault(protein["accession"], protein["sequence_file"])
        if accession in sequences and sequences[accession] != sequence:
            raise ValueError(f"同一蛋白编号对应不同序列：{accession}")
        existing_new = next((p for p in group["proteins"] if p["accession"] == accession), None)
        canonical = (existing_new["sequence_file"] if existing_new else
                     existing_paths.get(accession, sequence_file_value(path)))
        changed = old_accession != accession or old["selection_source"] != "manual"
        response = {"ok": True, "changed": changed, "selected_solution_id": solution_id,
                    "reaction_id": reaction, "step_indexes": route[reaction], "old_accession": old_accession,
                    "accession": accession, "sequence_file": canonical, "manifest_path": str(manifest_path),
                    "remaining_protein_count": len(group["proteins"]), "candidates_updated": False,
                    "combinations_cleared": False, "requires_candidate_regeneration": False}
        if not changed:
            return {**response, "message": "已是相同的人工蛋白，无变化"}
        if existing_new is not None and old_accession != accession:
            group["proteins"] = [p for p in group["proteins"] if p is not old]
            existing_new["selection_source"] = "manual"
        else:
            index = group["proteins"].index(old)
            group["proteins"][index] = {"accession": accession, "sequence_file": canonical, "selection_source": "manual"}
        response["remaining_protein_count"] = len(group["proteins"])
        selection["review_required_step_indexes"] = sorted(set(selection["review_required_step_indexes"]) | set(route[reaction]))
        _validate_structure(selection, raw)
        after = {**deepcopy(raw), "main_enzyme_selection": selection}
        from src.main_protein_selection.build_main_enzyme_sets import main_enzyme_set_paths
        sets = tuple(main_enzyme_set_paths(output_dir).values())
        sidecar = path.with_suffix(".json")
        files = dict.fromkeys((manifest_path, sidecar, *sets, *output_dir.glob("*.json"), *output_dir.glob("*.csv")))
        with ArtifactTransaction(files) as transaction:
            if not sidecar.exists():
                write_bytes_atomic(sidecar, (json.dumps(metadata, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
            was_fixed = any(p["selection_source"] == "manual" for p in
                            before["main_enzyme_selection"]["reaction_assignments"][reaction]["proteins"])
            sync = candidate_removal_updates if was_fixed else candidate_updates
            updates = sync(output_dir, before, after, genome_path)
            response["requires_candidate_regeneration"] = updates is None
            if updates is not None:
                publish_candidate_updates(updates)
                response["candidates_updated"] = True
                response["requires_candidate_regeneration"] = not json.loads(
                    updates[output_dir / "main_enzyme_selection.json"])["ok"]
            for artifact in sets:
                if artifact.exists():
                    artifact.unlink()
                    response["combinations_cleared"] = True
            if path.read_bytes() != original_fasta:
                raise ValueError("替换期间新蛋白序列文件已变化，请重试")
            read_sequence_pair(path, accession)
            for protein in retained:
                if read_protein_pair(protein)[1] != sequences[protein["accession"]]:
                    raise ValueError("替换期间保留的蛋白序列已变化，请重试")
            write_bytes_atomic(manifest_path, (json.dumps(after, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
            transaction.commit()
    response["message"] = ("已替换蛋白；请重新运行 main-enzyme 后再生成组合"
                           if response["requires_candidate_regeneration"] else "已替换蛋白，可重新生成主酶组合")
    return response
