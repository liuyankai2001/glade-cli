"""Incremental global/local CDS GC editing; raw sequences stay immutable."""

from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import io
import json
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from dnachisel import AvoidChanges, DnaOptimizationProblem, EnforceGCContent, EnforceTranslation

from src.protein_to_cds.sequence_constraints import (
    LEGACY_FORBIDDEN_MOTIFS, _OPTIMIZER_LOCK, assess_generated_cds,
    assess_uploaded_cds, uploaded_cds_protein, sha256_text,
)
from src.protein_to_cds.restriction_sites import (
    enzyme_constraints, normalize_enzymes, restriction_site_audit, with_restriction_audit,
)
from src.protein_to_cds.homopolymers import (
    homopolymer_audit, homopolymer_constraints, homopolymer_summary,
    saved_homopolymer_max, with_homopolymer_audit,
)
from src.write_manifest.cds_dependencies import cds_dependency_update
from src.write_manifest.cds_optimization import commit_cds_optimization
from src.write_manifest.store import read_design_manifest
from src.protein_to_cds.artifacts import raw_cds_metadata
from src.protein_to_cds.gc_constraints import (
    count_bounds, effective_gc_settings, gc_settings_audit, local_gc_summary,
    validate_gc_range, validate_window, window_gc_audit,
)

GC_MODE = "dna_chisel_gc_only"
GC_REPORT_SCHEMA = "protein_to_cds.gc_optimization.v2"
GC_REPORT_SCHEMAS = {GC_REPORT_SCHEMA, "protein_to_cds.gc_optimization.v1"}
RESTRICTION_MODE = "dna_chisel_restriction_sites"
RESTRICTION_REPORT_SCHEMA = "protein_to_cds.restriction_optimization.v1"
HOMOPOLYMER_MODE = "dna_chisel_homopolymer"
HOMOPOLYMER_REPORT_SCHEMA = "protein_to_cds.homopolymer_optimization.v1"


def _path(root: Path, value: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError("缺少 CDS 或报告文件路径")
    path = (root / value).resolve()
    if not path.is_relative_to(root):
        raise ValueError("CDS 文件路径超出当前项目输出目录")
    return path


def _read_record(root: Path, metadata: dict, accession: str, guards: dict) -> tuple[Path, str]:
    path = _path(root, metadata.get("path"))
    content = path.read_bytes()
    guards[path] = content
    if metadata.get("file_sha256") and hashlib.sha256(content).hexdigest() != metadata["file_sha256"]:
        raise ValueError(f"文件内容与来源记录不一致：{path.name}")
    record = SeqIO.read(io.StringIO(content.decode("utf-8")), "fasta")
    if record.id.upper() != accession:
        raise ValueError(f"FASTA 编号与 {accession} 不一致：{path.name}")
    sequence = str(record.seq).upper()
    if metadata.get("sequence_sha256") and sha256_text(sequence) != metadata["sequence_sha256"]:
        raise ValueError(f"序列与来源记录不一致：{path.name}")
    length = metadata.get("length_nt", metadata.get("length_aa"))
    if length is not None and len(sequence) != length:
        raise ValueError(f"序列长度与来源记录不一致：{path.name}")
    return path, sequence


def _read_report(root: Path, metadata: dict, guards: dict) -> dict:
    path = _path(root, metadata.get("path"))
    content = path.read_bytes()
    guards[path] = content
    if metadata.get("file_sha256") and hashlib.sha256(content).hexdigest() != metadata["file_sha256"]:
        raise ValueError(f"报告内容与来源记录不一致：{path.name}")
    report = json.loads(content)
    if not isinstance(report, dict):
        raise ValueError("CDS 报告必须是 JSON 对象")
    return report


def _changes(before: str, after: str) -> dict[str, int | float]:
    bases = sum(a != b for a, b in zip(before, after, strict=True))
    codons = sum(before[i:i + 3] != after[i:i + 3] for i in range(0, len(before), 3))
    return {
        "nucleotide_change_count": bases, "nucleotide_change_fraction": bases / len(before),
        "codon_change_count": codons, "codon_change_fraction": codons / (len(before) // 3),
    }


def _gc_count(sequence: str) -> int:
    return sequence.count("G") + sequence.count("C")


@dataclass
class CdsEditInput:
    item: dict
    previous: dict
    source_report: dict
    raw_meta: dict
    raw_path: Path
    raw: str
    protein: str
    organism_id: int
    extras: list[str]
    source_path: Path
    source_sequence: str
    uploaded: bool

    def audit(self, sequence: str, enzymes: list[str], maximum: int | None = None) -> dict:
        assess = assess_uploaded_cds if self.uploaded else assess_generated_cds
        return with_homopolymer_audit(
            with_restriction_audit(assess(sequence, self.protein, self.organism_id, self.extras), sequence, enzymes),
            sequence, maximum,
        )


def _load_edit_input(root: Path, selection: dict, item: dict, guards: dict) -> CdsEditInput:
    accession = item["accession"]
    previous = item.get("optimized_cds", {})
    if not previous.get("report", {}).get("path"):
        raise ValueError(f"{accession} 缺少原始副本或来源报告，请重新运行 protein-to-cds")
    report = _read_report(root, previous["report"], guards)
    raw_meta = raw_cds_metadata(item, report)
    if not isinstance(raw_meta, dict) or not raw_meta.get("sequence_sha256"):
        raise ValueError(f"{accession} 缺少可信 raw CDS，请重新运行 protein-to-cds")
    raw_path, raw = _read_record(root, raw_meta, accession, guards)
    if not raw_path.is_relative_to(root / "protein_to_cds/raw_cds"):
        raise ValueError("原始 CDS 必须位于当前项目的 raw_cds 目录")
    uploaded = item.get("sequence_input", {}).get("type") == "cds" or previous.get("source_type") == "user_uploaded"
    if uploaded:
        protein = uploaded_cds_protein(raw)
    else:
        _, protein = _read_record(root, item.get("protein_sequence", {}), accession, guards)
    organism_id = int(selection.get("host", {}).get("codon_transformer_organism_id", 0))
    if not organism_id:
        raise ValueError("CDS 结果缺少宿主信息")
    extras = report.get("additional_forbidden_motifs")
    if extras is None:
        extras = [v for k, v in report.get("constraints", {}).get("forbidden_motifs", {}).items() if k not in LEGACY_FORBIDDEN_MOTIFS]
    if not isinstance(extras, list) or any(not isinstance(value, str) for value in extras):
        raise ValueError("来源报告中的 motif 配置无效")
    if report.get("status") not in {"PASS", "IMPORTED"} or report.get("raw", {}).get("sequence_sha256") != sha256_text(raw):
        raise ValueError("CDS 来源报告与当前 raw 序列不一致")
    current_path = _path(root, previous.get("path"))
    if not current_path.is_file():
        raise ValueError("已登记的优化文件缺失，请重新运行 protein-to-cds")
    source_path, source_sequence = _read_record(root, previous, accession, guards)
    mode = previous.get("processing_mode")
    if not source_path.is_relative_to(root / "protein_to_cds/optimized_cds"):
        if mode != "codon_transformer_only" or source_path != raw_path:
            raise ValueError("当前 CDS 不在 optimized_cds 目录，请重新运行 protein-to-cds")
    output = _path(root, f"protein_to_cds/optimized_cds/{accession}.fasta")
    guards.setdefault(output, output.read_bytes() if output.exists() else None)
    if output.exists() and source_path != output:
        raise ValueError("已有优化文件来源不明，不能覆盖；请检查文件与 manifest")
    if report.get("final", {}).get("sequence_sha256") != sha256_text(source_sequence):
        raise ValueError("当前优化文件与来源报告不一致，不能继续编辑")
    edit_reports = {GC_MODE: (GC_REPORT_SCHEMAS, "gc_optimization"),
                    RESTRICTION_MODE: ({RESTRICTION_REPORT_SCHEMA}, "restriction_optimization"),
                    HOMOPOLYMER_MODE: ({HOMOPOLYMER_REPORT_SCHEMA}, "homopolymer_optimization")}
    if mode in edit_reports:
        schemas, suffix = edit_reports[mode]
        if (
            report.get("schema_version") not in schemas
            or _path(root, previous["report"]["path"]) != root / f"protein_to_cds/reports/{accession}.{suffix}.json"
            or report.get("source", {}).get("raw_sequence_sha256") != sha256_text(raw)
            or report.get("source", {}).get("cds_selection_source_fingerprint") != selection.get("source_fingerprint")
        ):
            raise ValueError("已有优化文件与当前 raw 来源不一致，不能继续编辑")
    elif mode not in {None, "codon_transformer_only", "codon_transformer_and_repair", "user_uploaded_cds"}:
        raise ValueError("不支持当前 CDS 的处理方式")
    if len(source_sequence) != len(raw) or source_sequence[:3] != raw[:3] or source_sequence[-3:] != raw[-3:]:
        raise ValueError("已有优化序列的长度或起止密码子与 raw 不一致")
    return CdsEditInput(item, previous, report, raw_meta, raw_path, raw, protein,
                        organism_id, extras, source_path, source_sequence, uploaded)


def _adjust_gc(
    sequence: str, protein: str, lower: int, upper: int, seed: int,
    *, local: dict[str, Any] | None = None, enzymes: list[str] | None = None,
    homopolymer_max: int | None = None,
) -> str:
    enzymes = normalize_enzymes(enzymes)
    if (lower <= _gc_count(sequence) <= upper
        and (local is None or window_gc_audit(sequence, local)["violation_count"] == 0)
        and restriction_site_audit(sequence, enzymes)["passed"]
        and (homopolymer_max is None or homopolymer_audit(sequence, homopolymer_max)["passed"])):
        return sequence
    # Fixed start/stop codons plus synonymous choices bound the possible GC count.
    codons = CodonTable.unambiguous_dna_by_id[11].forward_table
    possible = [[_gc_count(c) for c, aa in codons.items() if aa == residue] for residue in protein[1:]]
    fixed = _gc_count(sequence[:3] + sequence[-3:])
    if fixed + sum(min(values) for values in possible) > upper or fixed + sum(max(values) for values in possible) < lower:
        raise ValueError("指定 GC 范围无法在保持蛋白及起止密码子不变的条件下达到")
    length = len(sequence)
    constraints = [
        EnforceTranslation(genetic_table="Bacterial", start_codon="keep"),
        AvoidChanges(location=(0, 3)), AvoidChanges(location=(length - 3, length)),
        EnforceGCContent(mini=lower / length, maxi=upper / length),
        *enzyme_constraints(enzymes),
        *homopolymer_constraints(homopolymer_max, length),
    ]
    if local is not None:
        window = local["window_nt"]
        local_lower, local_upper = count_bounds(local["range_percent"], window)
        constraints.append(EnforceGCContent(mini=local_lower / window, maxi=local_upper / window, window=window))
    with _OPTIMIZER_LOCK:
        numpy_state, python_state = np.random.get_state(), random.getstate()
        try:
            np.random.seed(seed)
            random.seed(seed)
            problem = DnaOptimizationProblem(
                sequence,
                constraints=constraints,
                objectives=[AvoidChanges()], logger=None,
            )
            problem.resolve_constraints(final_check=True)
            problem.optimize()
            return str(problem.sequence).upper()
        except Exception as exc:
            raise ValueError(f"CDS 优化未找到满足全部约束的序列：{exc}") from exc
        finally:
            np.random.set_state(numpy_state)
            random.setstate(python_state)


def optimize_cds_gc(config: Any) -> dict[str, Any]:
    low, high = validate_gc_range(getattr(config, "gc_min", None), getattr(config, "gc_max", None))
    accession = str(getattr(config, "cds", "") or "").strip().upper()
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9_.-]*[A-Z0-9]|[A-Z0-9]", accession):
        raise ValueError("--cds 必须是一个有效蛋白编号")
    root = Path(config.project_output_path).expanduser().resolve()
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    if manifest.get("target_compound_id") != config.target_name:
        raise ValueError("manifest 目标与当前输入不一致，或尚未生成 CDS")
    selection = manifest.get("cds_selection")
    if not isinstance(selection, dict) or selection.get("schema_version") != "protein_to_cds.selection.v2":
        raise ValueError("请先运行 protein-to-cds 生成原始 CDS")
    matches = [item for item in selection.get("proteins", []) if item.get("accession", "").upper() == accession]
    if len(matches) != 1:
        raise ValueError(f"当前 CDS 结果中未找到唯一编号 {accession}")
    item = matches[0]
    guards: dict[Path, bytes | None] = {manifest_path: manifest_path.read_bytes()}
    source = _load_edit_input(root, selection, item, guards)
    previous, source_report = source.previous, source.source_report
    raw_meta, raw_path, raw, protein = source.raw_meta, source.raw_path, source.raw, source.protein
    extras = source.extras
    enzymes = normalize_enzymes(selection.get("restriction_enzymes"))
    maximum = saved_homopolymer_max(selection)
    raw_audit = source.audit(raw, enzymes, maximum)
    output = _path(root, f"protein_to_cds/optimized_cds/{accession}.fasta")
    report_path = _path(root, f"protein_to_cds/reports/{accession}.gc_optimization.json")
    guards.setdefault(output, output.read_bytes() if output.exists() else None)
    guards.setdefault(report_path, report_path.read_bytes() if report_path.exists() else None)
    source_path, source_sequence = source.source_path, source.source_sequence
    input_audit = source.audit(source_sequence, enzymes, maximum)
    length = len(raw)
    window = getattr(config, "window", None)
    if window is not None:
        window = validate_window(window, length)
    settings = effective_gc_settings(previous, source_report, low, high, window, length)
    lower, upper = count_bounds(settings["global"]["range_percent"], length) if "global" in settings else (0, length)
    _, input_local, input_pass = gc_settings_audit(source_sequence, settings)
    request = {
        "gc_min": settings["global" if window is None else "local"]["range_percent"][0],
        "gc_max": settings["global" if window is None else "local"]["range_percent"][1],
        "scope": "global" if window is None else "local", "gc_settings": settings,
        "raw_sequence_sha256": sha256_text(raw), "protein_sequence_sha256": sha256_text(protein),
        "dnachisel_version": importlib.metadata.version("dnachisel"),
        "algorithm_version": GC_REPORT_SCHEMA,
    }
    if enzymes:
        request["restriction_enzymes"] = enzymes
    if maximum is not None:
        request["homopolymer_max"] = maximum
    if (
        output.exists() and source_report.get("request") == request
        and input_pass and input_audit["restriction_site_audit"]["passed"]
        and (maximum is None or input_audit["homopolymer_audit"]["passed"])
        and source_report.get("raw", {}).get("user_forbidden_site_hits") == raw_audit["user_forbidden_site_hits"]
        and source_report.get("final", {}).get("user_forbidden_site_hits") == input_audit["user_forbidden_site_hits"]
        and source_path == output
        and not any("raw_cds" in row or "raw" in row.get("optimized_cds", {}).get("metrics", {}) for row in selection["proteins"])
    ):
        return {"accession": accession, "reused_existing": True, "output_path": str(output)}
    seed = int(sha256_text(json.dumps(request, sort_keys=True) + sha256_text(source_sequence))[:8], 16)
    options = {}
    if "local" in settings:
        options["local"] = settings["local"]
    if enzymes:
        options["enzymes"] = enzymes
    if maximum is not None:
        options["homopolymer_max"] = maximum
    final = _adjust_gc(source_sequence, protein, lower, upper, seed, **options)
    if len(final) != length or final[:3] != raw[:3] or final[-3:] != raw[-3:]:
        raise ValueError("优化结果改变了序列长度或起止密码子")
    final_audit = source.audit(final, enzymes, maximum)
    if maximum is not None and not final_audit["homopolymer_audit"]["passed"]:
        raise ValueError("优化结果仍包含超过用户阈值的同聚物，未保存")
    if not final_audit["restriction_site_audit"]["passed"]:
        raise ValueError("优化结果仍包含用户禁止的限制酶识别位点，未保存")
    global_validation, final_local, final_pass = gc_settings_audit(final, settings)
    if not final_pass:
        raise ValueError("优化结果的精确 GC 含量未达到全部整体/局部范围，未保存")
    changes = _changes(raw, final)
    output_relative = output.relative_to(root).as_posix()
    baseline = copy.deepcopy(previous.get("metrics", {}).get("raw") or source_report.get("raw") or raw_audit)
    baseline["user_forbidden_site_hits"] = raw_audit["user_forbidden_site_hits"]
    baseline["restriction_site_audit"] = raw_audit["restriction_site_audit"]
    if maximum is not None:
        baseline["homopolymer_audit"] = raw_audit["homopolymer_audit"]
        baseline["checks"] = {**baseline.get("checks", {}), "homopolymer_pass": raw_audit["checks"]["homopolymer_pass"]}
        baseline["failed_checks"] = [name for name, passed in baseline["checks"].items() if not passed]
        baseline["gate_status"] = "FAIL" if baseline["failed_checks"] else "PASS"
    report = {
        "schema_version": GC_REPORT_SCHEMA, "processing_mode": GC_MODE, "status": "PASS",
        "request": request, "deterministic_seed": seed,
        "source": {
            "raw_path": raw_path.relative_to(root).as_posix(), "raw_sequence_sha256": sha256_text(raw),
            "cds_selection_source_fingerprint": selection.get("source_fingerprint"),
            "input_path": source_path.relative_to(root).as_posix(),
            "input_sequence_sha256": sha256_text(source_sequence),
        },
        "additional_forbidden_motifs": extras,
        "gc_settings": settings,
        "restriction_enzymes": enzymes,
        "enforced_checks": ["encoding_identity", "start_stop_unchanged"] + [f"{scope}_gc" for scope in settings] + (["restriction_sites"] if enzymes else []),
        "raw_cds": raw_meta,
        "raw": baseline, "input": input_audit, "final": final_audit,
        "changes": changes, "changes_from_input": _changes(source_sequence, final),
        "restriction_validation": {"input": input_audit["restriction_site_audit"], "final": final_audit["restriction_site_audit"]},
    }
    if global_validation is not None:
        report["gc_validation"] = global_validation
    if final_local is not None:
        report["local_gc_validation"] = {"input": input_local, "final": final_local}
    if maximum is not None:
        report["homopolymer_max"] = maximum
        report["enforced_checks"].append("homopolymers")
        report["homopolymer_validation"] = {"input": input_audit["homopolymer_audit"], "final": final_audit["homopolymer_audit"]}
    report_bytes = (json.dumps(report, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    fasta = (f">{accession} gc_optimized\n" + "\n".join(final[i:i + 80] for i in range(0, length, 80)) + "\n").encode("utf-8")
    updated = copy.deepcopy(selection)
    if "restriction_enzymes" in selection:
        updated["restriction_enzymes"] = enzymes
    updated_item = next(row for row in updated["proteins"] if row["accession"].upper() == accession)
    updated_item["optimized_cds"] = {
        "path": output_relative, "file_sha256": hashlib.sha256(fasta).hexdigest(),
        "sequence_sha256": sha256_text(final), "length_nt": length,
        "report": {"path": report_path.relative_to(root).as_posix(), "file_sha256": hashlib.sha256(report_bytes).hexdigest()},
        "processing_mode": GC_MODE, "constraint_repair_applied": True,
        "optimization_skipped": False, "enforced_checks": report["enforced_checks"],
        "quality_checks_enforced": False, "gc_settings": settings,
        "restriction_enzymes": enzymes,
        "metrics": {"final": final_audit, "changes": changes},
    }
    if "global" in settings:
        updated_item["optimized_cds"]["gc_range_percent"] = settings["global"]["range_percent"]
    if final_local is not None:
        updated_item["optimized_cds"]["metrics"]["local_gc"] = local_gc_summary(input_local, final_local)
    if maximum is not None:
        updated_item["optimized_cds"]["homopolymer_max"] = maximum
        updated_item["optimized_cds"]["metrics"]["homopolymers"] = homopolymer_summary(input_audit["homopolymer_audit"], final_audit["homopolymer_audit"])
    if enzymes:
        updated_item["optimized_cds"]["metrics"]["restriction_sites"] = {
            "input_count": input_audit["restriction_site_audit"]["site_count"],
            "final_count": final_audit["restriction_site_audit"]["site_count"],
        }
    files = {output: fasta, report_path: report_bytes}
    for row in updated["proteins"]:
        current = row["optimized_cds"]
        # Moving legacy raw metadata out of the manifest must also preserve
        # other proteins' baselines and give every current reference a copy.
        if row is not updated_item:
            reference = current.get("report", {})
            if "raw_cds" in row and reference.get("path"):
                legacy_report = _read_report(root, reference, guards)
                if "raw_cds" not in legacy_report:
                    legacy_report["raw_cds"] = row["raw_cds"]
                    content = (json.dumps(legacy_report, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
                    files[_path(root, reference["path"])] = content
                    reference["file_sha256"] = hashlib.sha256(content).hexdigest()
            path = _path(root, current.get("path"))
            if not path.is_relative_to(root / "protein_to_cds/optimized_cds"):
                if current.get("optimization_skipped") is not True and (
                    current.get("processing_mode") != "codon_transformer_only"
                    or path != root / "protein_to_cds/raw_cds" / f"{row['accession']}.raw.fasta"
                ):
                    raise ValueError("旧 CDS 当前引用来源不明，请重新运行 protein-to-cds")
                _read_record(root, current, row["accession"], guards)
                working = _path(root, f"protein_to_cds/optimized_cds/{row['accession']}.fasta")
                if working.exists():
                    raise ValueError("已有优化文件来源不明，不能覆盖；请检查文件与 manifest")
                guards[working] = None
                files[working] = guards[path]
                current["path"] = working.relative_to(root).as_posix()
                current["file_sha256"] = hashlib.sha256(files[working]).hexdigest()
        row.pop("raw_cds", None)
        current.get("metrics", {}).pop("raw", None)
        current.pop("source_raw_sequence_sha256", None)
    retained_sections, discard_sections = cds_dependency_update(manifest, updated)
    commit_cds_optimization(
        manifest_path=manifest_path, project_root=root, target=config.target_name,
        revision=int(manifest.get("revision", 0)), selection=updated,
        discard_sections=discard_sections, dependent_sections=retained_sections,
        files=files, guards=guards,
    )
    return {"accession": accession, "reused_existing": False, "output_path": str(output)}


def run_gc_optimization(config: Any) -> dict[str, Any]:
    # Import after core modules are initialized to avoid info/writer import cycles.
    from src.info_show.cds_info import format_cds_info, get_cds_info

    try:
        result = optimize_cds_gc(config)
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"CDS 优化失败：{exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    view = get_cds_info(config, accession=result["accession"])
    print(format_cds_info(view))
    return {**view, "复用已有结果": result["reused_existing"]}


__all__ = ["GC_MODE", "validate_gc_range", "optimize_cds_gc", "run_gc_optimization"]
