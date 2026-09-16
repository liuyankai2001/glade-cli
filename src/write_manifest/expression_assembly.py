"""Assemble either prepared source through one transactional entry point."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio.Seq import Seq

from src.expression_box.config import RBS_CONTEXT_CDS_PREFIX_NT, RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT
from src.expression_box.expression_burden import (
    calculate_expression_burden, expression_burden_summary, validate_expression_burden,
)
from src.expression_box.ostir_adapter import predict_rbs_context
from src.expression_box.parts_manifest_adapter import load_expression_parts_context
from src.expression_box.parts_preparation import (
    SELECTED_PARTS_PATH, SELECTED_PARTS_SCHEMA, bind_cds, make_parts_draft,
    normalize_parts_draft, project_file, read_part_sequence, stable_hash,
)
from src.expression_box.sequence_audit import audit_expression_sequence
from src.pathway_analyze.target_id import validate_target_compound_id
from src.protein_to_cds.artifacts import ArtifactTransaction, write_bytes_atomic
from src.write_manifest.expression_constructs import prepare_expression_constructs
from src.write_manifest.expression_parts_draft import _predict_uploaded_rbs
from src.write_manifest.store import read_design_manifest, update_design_manifest


ASSEMBLY_DOWNSTREAM = (
    "assembled_expression_cassettes", "plasmid_selection", "final_assembly_plan", "final_assembly",
)


def _prediction_current(ostir, rbs_sequence, cassette, index):
    if not isinstance(ostir, Mapping):
        return False
    upstream = cassette.cds[index - 1].sequence[-RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT:] if index else ""
    sequence = upstream + rbs_sequence + cassette.cds[index].sequence[:RBS_CONTEXT_CDS_PREFIX_NT]
    parameters = {"previous_cds_suffix_nt": RBS_CONTEXT_PREVIOUS_CDS_SUFFIX_NT,
                  "current_cds_prefix_nt": RBS_CONTEXT_CDS_PREFIX_NT}
    try:
        expression = float(ostir["translation_initiation_rate"])
        energy = float(ostir["d_g_total"])
        return (
            math.isfinite(expression) and expression > 0 and math.isfinite(energy)
            and type(ostir["unintended_start_count"]) is int and ostir["unintended_start_count"] >= 0
            and ostir.get("intended_start_position") == len(upstream) + len(rbs_sequence) + 1
            and ostir.get("context_sha256") == hashlib.sha256(sequence.encode("utf-8")).hexdigest()
            and ostir.get("context_parameters", parameters) == parameters
        )
    except (KeyError, TypeError, ValueError):
        return False


def _normalize_design(prepared, context, source_type, guarded, predictor):
    design = {"design_id": prepared["design_id"], "rank": prepared["rank"],
              "name": prepared.get("name", ""), "source_type": source_type, "cassettes": []}
    raw_cassettes = prepared["cassettes"]
    if len(raw_cassettes) != len(context.cassettes):
        raise ValueError("元件方案与当前表达盒数量不一致")
    prediction_changed = False
    for raw, cassette in zip(raw_cassettes, context.cassettes, strict=True):
        accessions = [cds.accession for cds in cassette.cds]
        if raw.get("cassette_index") != cassette.cassette_index or raw.get("protein_accessions") != accessions:
            raise ValueError("元件方案与当前表达盒分组或 CDS 顺序不一致")
        missing = []
        for role in ("promoter", "terminator"):
            if not isinstance(raw.get(role), Mapping):
                missing.append(role)
        for accession in accessions:
            if not isinstance(raw.get("rbs_by_accession", {}).get(accession), Mapping):
                missing.append(f"RBS:{accession}")
        if missing:
            raise ValueError(f"表达盒 {cassette.cassette_index} 缺少元件：{', '.join(missing)}")

        def loaded(part, role):
            if part.get("role") != role:
                raise ValueError(f"表达盒 {cassette.cassette_index} 元件角色不匹配：{role}")
            return {**copy.deepcopy(dict(part)), "sequence": read_part_sequence(part, context.project_root, guarded)}

        promoter = loaded(raw["promoter"], "promoter")
        terminator = loaded(raw["terminator"], "terminator")
        genes, chunks = [], [promoter["sequence"]]
        for index, cds in enumerate(cassette.cds):
            reference = raw.get("optimized_cds_by_accession", {}).get(cds.accession)
            if reference is not None and (
                reference.get("sequence_sha256") != cds.sequence_sha256
                or reference.get("length_nt") != len(cds.sequence)
                or project_file(context.project_root, reference.get("path")) != cds.path.resolve()
            ):
                raise ValueError(f"{cds.accession} 的当前 CDS 引用已变化，请更新元件准备记录")
            part = raw["rbs_by_accession"][cds.accession]
            rbs = loaded(part, "rbs")
            ostir = part.get("ostir")
            if part.get("ostir_status") == "stale" or not _prediction_current(ostir, rbs["sequence"], cassette, index):
                ostir = _predict_uploaded_rbs(
                    manifest_path=context.manifest_path, project_root=context.project_root,
                    cassette_index=cassette.cassette_index, cds_index=index,
                    accession=cds.accession, part_id=part["part_id"], rbs_sequence=rbs["sequence"],
                    predictor=predictor,
                )
                part["ostir"] = ostir
                part.pop("ostir_status", None)
                prediction_changed = True
            rbs.pop("ostir", None)
            rbs.pop("ostir_status", None)
            genes.append({"accession": cds.accession, "cds_sequence_sha256": cds.sequence_sha256,
                          "cds_length_nt": len(cds.sequence), "optimized_cds": copy.deepcopy(reference),
                          "rbs": rbs, "ostir": copy.deepcopy(ostir)})
            chunks.extend((rbs["sequence"], cds.sequence))
        chunks.append(terminator["sequence"])
        sequence = "".join(chunks)
        audit = audit_expression_sequence(sequence, context.restriction_enzymes, context.homopolymer_max)
        # The complete construct builder reports conflicts with global coordinates
        # and feature labels, including cassette junctions.
        raw["sequence_audit"] = audit
        raw["restriction_site_audit"] = copy.deepcopy(audit["restriction_site_audit"])
        raw["homopolymer_audit"] = copy.deepcopy(audit["homopolymer_audit"])
        raw.pop("restriction_audit_status", None)
        raw.pop("homopolymer_audit_status", None)
        design["cassettes"].append({"cassette_index": cassette.cassette_index,
                                   "promoter": promoter, "terminator": terminator, "genes": genes,
                                   "sequence_audit": audit,
                                   "assembled_sequence": {"stored": False, "length_nt": len(sequence),
                                                          "sequence_sha256": audit["sequence_sha256"]}})

    recommendation = prepared.get("recommendation", {})
    fallback = 50.0
    burden = None
    if recommendation.get("status") == "current" and not prediction_changed:
        fallback = float(recommendation.get("expression_target_percentile", 50.0))
        burden = copy.deepcopy(recommendation.get("expression_burden"))
        if burden is not None:
            validate_expression_burden(burden, design["cassettes"], fallback_promoter_percentile=fallback)
    if burden is None:
        burden = calculate_expression_burden(
            design["cassettes"], {}, fallback_promoter_percentile=fallback,
            fallback_promoter_source="neutral_percentile_fallback",
        )
    design.update(expression_burden=burden, estimated_burden=burden["level"],
                  expression_target_percentile=fallback, warnings=list(burden.get("warnings", [])))
    if recommendation:
        design["recommendation"] = copy.deepcopy(recommendation)
        design["strategy"] = recommendation.get("strategy", "")
        design["expression_regime"] = recommendation.get("expression_regime", "")
        if recommendation.get("status") == "current" and not prediction_changed:
            design["expression_success_score"] = recommendation.get("expression_success_score")
        else:
            recommendation["status"] = "stale"
            design["recommendation"]["status"] = "stale"
    return design


def _selection_payload(designs, draft, context, snapshot_bytes):
    references = [{"design_id": design["design_id"], "rank": design["rank"],
                   "expression_success_score": design.get("expression_success_score"),
                   "expression_regime": design.get("expression_regime", ""),
                   "expression_burden": expression_burden_summary(design["expression_burden"]),
                   "system_recommended": draft["source_type"] == "recommended",
                   "design_fingerprint": stable_hash(design)} for design in designs]
    content = {"context_input_fingerprint": context.input_fingerprint,
               "draft_fingerprint": draft["draft_fingerprint"], "design_references": references}
    return {"schema_version": "parts_selection.v2", "status": "selected",
            "selection_status": "user_selected", "source_type": draft["source_type"],
            "selected_design_ids": [d["design_id"] for d in designs],
            "primary_design_id": designs[0]["design_id"], "design_count": len(designs),
            "selection_fingerprint": stable_hash(content), "design_references": references,
            "source": {**draft["source"], "artifact": SELECTED_PARTS_PATH,
                       "artifact_file_sha256": hashlib.sha256(snapshot_bytes).hexdigest(),
                       "designs_schema_version": SELECTED_PARTS_SCHEMA,
                       "context_input_fingerprint": context.input_fingerprint,
                       "draft_fingerprint": draft["draft_fingerprint"]},
            "warnings": list(dict.fromkeys(w for d in designs for w in d.get("warnings", [])))}


def assemble_expression_constructs(config: Any, *, predictor=predict_rbs_context) -> dict[str, Any]:
    root = Path(config.project_output_path).expanduser().resolve()
    snapshot_path = project_file(root, SELECTED_PARTS_PATH)
    with ArtifactTransaction([snapshot_path], lock_path=root / "protein_to_cds" / ".gc_optimization.lock") as files:
        path = Path(config.manifest_output_path).expanduser().resolve()
        original_manifest = path.read_bytes()
        manifest = read_design_manifest(path)
        target = validate_target_compound_id(config.target_name)
        if manifest.get("target_compound_id") != target:
            raise ValueError("manifest 与当前目标化合物不一致")
        context = load_expression_parts_context(path, root)
        if context.manifest_revision != int(manifest["revision"]):
            raise ValueError("manifest revision 已变化，请重试")
        for cassette in context.cassettes:
            for cds in cassette.cds:
                try:
                    Seq(cds.sequence).translate(table=11, cds=True)
                except Exception as exc:
                    raise ValueError(f"{cds.accession} 的 CDS 不是有效完整编码序列：{exc}") from exc
        raw_draft = manifest.get("expression_parts_draft")
        if not isinstance(raw_draft, Mapping):
            raise ValueError("请先上传表达元件，或使用 write --expression-parts 确认推荐方案")
        draft = normalize_parts_draft(raw_draft)
        if draft["target_compound_id"] != target or (
            draft["source"].get("expression_box_selection_fingerprint") != context.expression_box_selection_fingerprint
            or draft["source"].get("cds_selection_source_fingerprint") != context.cds_selection_source_fingerprint
        ):
            raise ValueError("元件准备记录与当前表达盒分组或 CDS 来源不一致")
        guarded = {cds.path: cds.path.read_bytes() for c in context.cassettes for cds in c.cds}
        proteins = {row["accession"]: row for row in manifest["cds_selection"]["proteins"]}
        for design in draft["designs"]:
            if len(design["cassettes"]) != len(context.cassettes):
                raise ValueError("元件方案与当前表达盒数量不一致")
            for raw, cassette in zip(design["cassettes"], context.cassettes, strict=True):
                if raw["protein_accessions"] != [cds.accession for cds in cassette.cds]:
                    raise ValueError("元件方案与当前表达盒分组或 CDS 顺序不一致")
                if "optimized_cds_by_accession" not in raw:
                    bind_cds([raw], proteins)
        designs = [_normalize_design(d, context, draft["source_type"], guarded, predictor) for d in draft["designs"]]
        for design in draft["designs"]:
            bind_cds(design["cassettes"], proteins)
        draft = make_parts_draft(target=target, source=draft["source"], source_type=draft["source_type"], designs=draft["designs"])
        snapshot = {"schema_version": SELECTED_PARTS_SCHEMA, "target_compound_id": target,
                    "source": {**draft["source"], "context_input_fingerprint": context.input_fingerprint},
                    "source_type": draft["source_type"], "designs": designs}
        snapshot_bytes = (json.dumps(snapshot, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode("utf-8")
        selection = _selection_payload(designs, draft, context, snapshot_bytes)
        transaction = prepare_expression_constructs(selected_designs=designs, selection_payload=selection,
                                                    context=context, current_section=manifest.get("assembled_expression_constructs"))
        changed = (manifest.get("parts_selection") != selection
                   or manifest.get("assembled_expression_constructs") != transaction.section
                   or manifest.get("expression_parts_draft") != draft)
        try:
            if path.read_bytes() != original_manifest or any(p.read_bytes() != content for p, content in guarded.items()):
                raise ValueError("manifest、CDS 或元件快照在构建时发生变化，请重试")
            if load_expression_parts_context(path, root).input_fingerprint != context.input_fingerprint:
                raise ValueError("当前 CDS 或表达盒上下文在构建时发生变化，请重试")
            transaction.install()
            if not snapshot_path.is_file() or snapshot_path.read_bytes() != snapshot_bytes:
                write_bytes_atomic(snapshot_path, snapshot_bytes)
            updated = update_design_manifest(
                path, target_compound_id=target,
                sections={"expression_parts_draft": draft, "parts_selection": selection,
                          "assembled_expression_constructs": transaction.section},
                discard_sections=ASSEMBLY_DOWNSTREAM, expected_revision=context.manifest_revision,
            ) if changed else manifest
            files.commit()
        except Exception:
            transaction.rollback()
            raise
        cleanup_warning = transaction.finalize()
        return {"运行成功": True, "目标化合物": target, "元件来源": draft["source_type"],
                "表达元件方案编号": selection["selected_design_ids"], "完整表达构建数量": len(designs),
                "GenBank目录": str(root / "expression_constructs"), "GenBank是否复用": not transaction.needs_install,
                "GenBank是否修复": transaction.is_repair, "清单是否更新": changed,
                "清单文件": str(path), "清单版本": updated["revision"],
                "警告": selection["warnings"] + ([cleanup_warning] if cleanup_warning else [])}
