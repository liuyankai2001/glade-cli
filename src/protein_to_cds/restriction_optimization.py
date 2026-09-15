"""Atomically remove user-selected enzyme sites from every current CDS."""

from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path
from typing import Any

from src.protein_to_cds.gc_constraints import count_bounds, gc_settings_audit, local_gc_summary, saved_gc_settings
from src.protein_to_cds.gc_optimization import (
    RESTRICTION_MODE, RESTRICTION_REPORT_SCHEMA, _adjust_gc, _changes,
    _load_edit_input, _path,
)
from src.protein_to_cds.restriction_sites import normalize_enzymes
from src.protein_to_cds.sequence_constraints import sha256_text
from src.write_manifest.cds_dependencies import cds_dependency_update
from src.write_manifest.cds_optimization import commit_cds_optimization
from src.write_manifest.store import read_design_manifest


def optimize_restriction_sites(config: Any) -> dict[str, Any]:
    values = getattr(config, "enzyme", None)
    if not values:
        raise ValueError("请通过 --enzyme 指定一个或多个限制酶")
    if any(getattr(config, name, None) is not None for name in ("cds", "gc_min", "gc_max", "window")):
        raise ValueError("--enzyme 批量模式不能与 --cds、GC 或窗口参数混用")
    if isinstance(values, str):
        values = [values]
    clear = len(values) == 1 and isinstance(values[0], str) and values[0].strip().lower() == "none"
    enzymes = [] if clear else normalize_enzymes(values)
    root = Path(config.project_output_path).expanduser().resolve()
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest.get("target_compound_id") != config.target_name:
        raise ValueError("manifest 目标与当前输入不一致")
    selection = manifest.get("cds_selection")
    if not isinstance(selection, dict) or selection.get("schema_version") != "protein_to_cds.selection.v2" or selection.get("status") != "complete":
        raise ValueError("请先运行 protein-to-cds，得到完整 CDS 批次")
    proteins = selection.get("proteins")
    if not isinstance(proteins, list) or not proteins:
        raise ValueError("当前没有可编辑的 CDS")
    accessions = [row.get("accession") for row in proteins]
    if any(not isinstance(name, str) or not name for name in accessions) or len(set(name.upper() for name in accessions)) != len(accessions):
        raise ValueError("当前 CDS 批次的蛋白编号无效或重复")
    guards = {manifest_path: manifest_bytes}
    updated = copy.deepcopy(selection)
    updated["restriction_enzymes"] = enzymes
    files = {}
    results = []
    for item, updated_item in zip(proteins, updated["proteins"], strict=True):
        source = _load_edit_input(root, selection, item, guards)
        accession = item["accession"]
        settings = saved_gc_settings(source.previous, source.source_report, len(source.raw))
        raw_audit = source.audit(source.raw, enzymes)
        input_audit = source.audit(source.source_sequence, enzymes)
        _, input_local, input_pass = gc_settings_audit(source.source_sequence, settings)
        output = _path(root, f"protein_to_cds/optimized_cds/{accession}.fasta")
        report_path = _path(root, f"protein_to_cds/reports/{accession}.restriction_optimization.json")
        guards.setdefault(report_path, report_path.read_bytes() if report_path.exists() else None)
        request = {
            "restriction_enzymes": enzymes, "gc_settings": settings,
            "raw_sequence_sha256": sha256_text(source.raw), "protein_sequence_sha256": sha256_text(source.protein),
            "dnachisel_version": importlib.metadata.version("dnachisel"), "algorithm_version": RESTRICTION_REPORT_SCHEMA,
        }
        reused = (
            selection.get("restriction_enzymes") == enzymes
            and source.previous.get("restriction_enzymes") == enzymes
            and input_pass and input_audit["restriction_site_audit"]["passed"]
            and source.source_report.get("restriction_enzymes") == enzymes
            and source.source_path == output and "raw_cds" not in item
            and "raw" not in source.previous.get("metrics", {})
        )
        if reused:
            results.append({"accession": accession, "reused_existing": True, "output_path": str(output)})
            continue
        lower, upper = count_bounds(settings["global"]["range_percent"], len(source.raw)) if "global" in settings else (0, len(source.raw))
        seed = int(sha256_text(json.dumps(request, sort_keys=True) + sha256_text(source.source_sequence))[:8], 16)
        final = _adjust_gc(source.source_sequence, source.protein, lower, upper, seed,
                           local=settings.get("local"), enzymes=enzymes)
        if len(final) != len(source.raw) or final[:3] != source.raw[:3] or final[-3:] != source.raw[-3:]:
            raise ValueError(f"{accession} 优化结果改变了长度或起止密码子，整批未保存")
        final_audit = source.audit(final, enzymes)
        global_audit, final_local, passed = gc_settings_audit(final, settings)
        if not passed or not final_audit["restriction_site_audit"]["passed"]:
            raise ValueError(f"{accession} 优化结果未满足 GC 或限制酶约束，整批未保存")
        report = {
            "schema_version": RESTRICTION_REPORT_SCHEMA, "processing_mode": RESTRICTION_MODE, "status": "PASS",
            "request": request, "deterministic_seed": seed,
            "source": {"raw_path": source.raw_path.relative_to(root).as_posix(), "raw_sequence_sha256": sha256_text(source.raw),
                       "input_path": source.source_path.relative_to(root).as_posix(), "input_sequence_sha256": sha256_text(source.source_sequence),
                       "cds_selection_source_fingerprint": selection.get("source_fingerprint")},
            "raw_cds": source.raw_meta, "additional_forbidden_motifs": source.extras,
            "restriction_enzymes": enzymes, "gc_settings": settings,
            "enforced_checks": ["encoding_identity", "start_stop_unchanged"] + [f"{scope}_gc" for scope in settings] + (["restriction_sites"] if enzymes else []),
            "raw": raw_audit, "input": input_audit, "final": final_audit,
            "changes": _changes(source.raw, final), "changes_from_input": _changes(source.source_sequence, final),
            "restriction_validation": {"input": input_audit["restriction_site_audit"], "final": final_audit["restriction_site_audit"]},
        }
        if global_audit is not None:
            report["gc_validation"] = global_audit
        if final_local is not None:
            report["local_gc_validation"] = {"input": input_local, "final": final_local}
        report_bytes = (json.dumps(report, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
        fasta = (f">{accession} restriction_optimized\n" + "\n".join(final[i:i + 80] for i in range(0, len(final), 80)) + "\n").encode("utf-8")
        files.update({output: fasta, report_path: report_bytes})
        metadata = {
            "path": output.relative_to(root).as_posix(), "file_sha256": hashlib.sha256(fasta).hexdigest(),
            "sequence_sha256": sha256_text(final), "length_nt": len(final),
            "report": {"path": report_path.relative_to(root).as_posix(), "file_sha256": hashlib.sha256(report_bytes).hexdigest()},
            "processing_mode": RESTRICTION_MODE, "constraint_repair_applied": True, "optimization_skipped": False,
            "quality_checks_enforced": False, "enforced_checks": report["enforced_checks"],
            "restriction_enzymes": enzymes, "gc_settings": settings,
            "metrics": {"final": final_audit, "changes": report["changes"]},
        }
        if enzymes:
            metadata["metrics"]["restriction_sites"] = {
                "input_count": input_audit["restriction_site_audit"]["site_count"],
                "final_count": final_audit["restriction_site_audit"]["site_count"],
            }
        if "global" in settings:
            metadata["gc_range_percent"] = settings["global"]["range_percent"]
        if final_local is not None:
            metadata["metrics"]["local_gc"] = local_gc_summary(input_local, final_local)
        updated_item["optimized_cds"] = metadata
        updated_item.pop("raw_cds", None)
        results.append({"accession": accession, "reused_existing": False, "output_path": str(output)})
    if files or selection.get("restriction_enzymes") != enzymes:
        retained, discard = cds_dependency_update(manifest, updated)
        commit_cds_optimization(manifest_path=manifest_path, project_root=root, target=config.target_name,
                                revision=int(manifest.get("revision", 0)), selection=updated,
                                discard_sections=discard, dependent_sections=retained, files=files, guards=guards)
    return {"restriction_enzymes": enzymes, "results": results,
            "reused_existing": all(row["reused_existing"] for row in results)}


def run_restriction_optimization(config: Any) -> dict[str, Any]:
    from src.info_show.cds_info import format_cds_info, get_cds_info
    try:
        result = optimize_restriction_sites(config)
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"限制酶位点消除失败：{exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    view = get_cds_info(config)
    print(format_cds_info(view))
    return {**view, **result}


def run_cds_optimization(config: Any) -> dict[str, Any]:
    if getattr(config, "enzyme", None) is not None:
        return run_restriction_optimization(config)
    if getattr(config, "gc_min", None) is None or getattr(config, "gc_max", None) is None:
        print("GC 模式必须同时提供 --gc-min 和 --gc-max", file=sys.stderr)
        raise SystemExit(2)
    from src.protein_to_cds.gc_optimization import run_gc_optimization
    return run_gc_optimization(config)
