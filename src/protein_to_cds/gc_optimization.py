"""Explicit, incremental whole-CDS GC editing; raw sequences stay immutable."""

from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import io
import json
import random
import re
import sys
from decimal import Decimal, InvalidOperation, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path
from typing import Any

import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable
from dnachisel import AvoidChanges, DnaOptimizationProblem, EnforceGCContent, EnforceTranslation

from src.protein_to_cds.sequence_constraints import (
    DEFAULT_FORBIDDEN_MOTIFS, _OPTIMIZER_LOCK, assess_generated_cds, sha256_text,
)
from src.protein_to_cds.write_protein_to_manifest import CDS_SELECTION_DOWNSTREAM_SECTIONS
from src.write_manifest.cds_optimization import commit_cds_optimization
from src.write_manifest.store import read_design_manifest

GC_MODE = "dna_chisel_gc_only"
GC_REPORT_SCHEMA = "protein_to_cds.gc_optimization.v1"


def validate_gc_range(minimum: Any, maximum: Any) -> tuple[Decimal, Decimal]:
    try:
        low, high = Decimal(str(minimum)), Decimal(str(maximum))
    except InvalidOperation as exc:
        raise ValueError("--gc-min 和 --gc-max 必须是百分比数值") from exc
    if not low.is_finite() or not high.is_finite() or not Decimal(0) <= low <= high <= Decimal(100):
        raise ValueError("GC 范围必须满足 0 <= gc-min <= gc-max <= 100")
    return low, high


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


def _adjust_gc(sequence: str, protein: str, lower: int, upper: int, seed: int) -> str:
    if lower <= _gc_count(sequence) <= upper:
        return sequence
    # Fixed start/stop codons plus synonymous choices bound the possible GC count.
    codons = CodonTable.unambiguous_dna_by_id[11].forward_table
    possible = [[_gc_count(c) for c, aa in codons.items() if aa == residue] for residue in protein[1:]]
    fixed = _gc_count(sequence[:3] + sequence[-3:])
    if fixed + sum(min(values) for values in possible) > upper or fixed + sum(max(values) for values in possible) < lower:
        raise ValueError("指定 GC 范围无法在保持蛋白及起止密码子不变的条件下达到")
    length = len(sequence)
    with _OPTIMIZER_LOCK:
        numpy_state, python_state = np.random.get_state(), random.getstate()
        try:
            np.random.seed(seed)
            random.seed(seed)
            problem = DnaOptimizationProblem(
                sequence,
                constraints=[
                    EnforceTranslation(genetic_table="Bacterial", start_codon="keep"),
                    AvoidChanges(location=(0, 3)),
                    AvoidChanges(location=(length - 3, length)),
                    EnforceGCContent(mini=lower / length, maxi=upper / length),
                ],
                objectives=[AvoidChanges()], logger=None,
            )
            problem.resolve_constraints(final_check=True)
            problem.optimize()
            return str(problem.sequence).upper()
        except Exception as exc:
            raise ValueError(f"整体 GC 优化未找到满足约束的序列：{exc}") from exc
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
    raw_meta = item.get("raw_cds")
    if not isinstance(raw_meta, dict) or not raw_meta.get("sequence_sha256"):
        raise ValueError(f"{accession} 没有 CodonTransformer 原始 CDS；直接上传 CDS 暂不支持此命令")
    guards: dict[Path, bytes | None] = {manifest_path: manifest_path.read_bytes()}
    raw_path, raw = _read_record(root, raw_meta, accession, guards)
    if not raw_path.is_relative_to(root / "protein_to_cds/raw_cds"):
        raise ValueError("原始 CDS 必须位于当前项目的 raw_cds 目录")
    _, protein = _read_record(root, item.get("protein_sequence", {}), accession, guards)
    organism_id = int(selection.get("host", {}).get("codon_transformer_organism_id", 0))
    if not organism_id:
        raise ValueError("CDS 结果缺少宿主信息")
    previous = item.get("optimized_cds", {})
    source_report = _read_report(root, previous.get("report", {}), guards)
    extras = source_report.get("additional_forbidden_motifs")
    if extras is None:
        extras = [v for k, v in source_report.get("constraints", {}).get("forbidden_motifs", {}).items() if k not in DEFAULT_FORBIDDEN_MOTIFS]
    if not isinstance(extras, list) or any(not isinstance(value, str) for value in extras):
        raise ValueError("来源报告中的 motif 配置无效")
    raw_audit = assess_generated_cds(raw, protein, organism_id, extras)
    if source_report.get("status") != "PASS" or source_report.get("raw", {}).get("sequence_sha256") != sha256_text(raw):
        raise ValueError("CDS 来源报告与当前 raw 序列不一致")
    output = _path(root, f"protein_to_cds/optimized_cds/{accession}.fasta")
    report_path = _path(root, f"protein_to_cds/reports/{accession}.gc_optimization.json")
    guards.setdefault(output, output.read_bytes() if output.exists() else None)
    guards.setdefault(report_path, report_path.read_bytes() if report_path.exists() else None)
    source_path, source_sequence = raw_path, raw
    if output.exists():
        if previous.get("processing_mode") != GC_MODE or _path(root, previous.get("path")) != output:
            raise ValueError("已有优化文件来源不明或不是当前版本，不能覆盖；请检查文件与 manifest")
        source_path, source_sequence = _read_record(root, previous, accession, guards)
        if (
            source_report.get("schema_version") != GC_REPORT_SCHEMA
            or _path(root, previous.get("report", {}).get("path")) != report_path
            or source_report.get("source", {}).get("raw_sequence_sha256") != sha256_text(raw)
            or source_report.get("source", {}).get("cds_selection_source_fingerprint") != selection.get("source_fingerprint")
            or source_report.get("final", {}).get("sequence_sha256") != sha256_text(source_sequence)
        ):
            raise ValueError("已有优化文件与当前 raw 来源不一致，不能继续编辑")
    elif report_path.exists() or previous.get("processing_mode") == GC_MODE:
        raise ValueError("已登记的优化文件缺失或优化报告孤立，不能自动从 raw 重建")
    input_audit = assess_generated_cds(source_sequence, protein, organism_id, extras)
    if len(source_sequence) != len(raw) or source_sequence[:3] != raw[:3] or source_sequence[-3:] != raw[-3:]:
        raise ValueError("已有优化序列的长度或起止密码子与 raw 不一致")
    length = len(raw)
    lower = int((low * length / 100).to_integral_value(rounding=ROUND_CEILING))
    upper = int((high * length / 100).to_integral_value(rounding=ROUND_FLOOR))
    if lower > upper:
        raise ValueError("指定 GC 范围对当前序列长度没有可取的整数 GC 数量")
    request = {
        "gc_min": str(low.normalize()), "gc_max": str(high.normalize()),
        "raw_sequence_sha256": sha256_text(raw), "protein_sequence_sha256": sha256_text(protein),
        "dnachisel_version": importlib.metadata.version("dnachisel"),
        "algorithm_version": GC_REPORT_SCHEMA,
    }
    if (
        output.exists() and source_report.get("request") == request
        and lower <= _gc_count(source_sequence) <= upper
        and source_report.get("raw", {}).get("user_forbidden_site_hits") == raw_audit["user_forbidden_site_hits"]
        and source_report.get("final", {}).get("user_forbidden_site_hits") == input_audit["user_forbidden_site_hits"]
    ):
        return {"accession": accession, "reused_existing": True, "output_path": str(output)}
    seed = int(sha256_text(json.dumps(request, sort_keys=True) + sha256_text(source_sequence))[:8], 16)
    final = _adjust_gc(source_sequence, protein, lower, upper, seed)
    if len(final) != length or final[:3] != raw[:3] or final[-3:] != raw[-3:]:
        raise ValueError("优化结果改变了序列长度或起止密码子")
    final_audit = assess_generated_cds(final, protein, organism_id, extras)
    if not lower <= _gc_count(final) <= upper:
        raise ValueError("优化结果的精确 GC 含量未达到用户范围，未保存")
    changes = _changes(raw, final)
    output_relative = output.relative_to(root).as_posix()
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
        "enforced_checks": ["encoding_identity", "start_stop_unchanged", "global_gc"],
        "gc_validation": {"minimum_count": lower, "maximum_count": upper, "observed_count": _gc_count(final)},
        "raw": raw_audit, "input": input_audit, "final": final_audit,
        "changes": changes, "changes_from_input": _changes(source_sequence, final),
    }
    report_bytes = (json.dumps(report, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    fasta = (f">{accession} gc_optimized\n" + "\n".join(final[i:i + 80] for i in range(0, length, 80)) + "\n").encode("utf-8")
    updated = copy.deepcopy(selection)
    updated_item = next(row for row in updated["proteins"] if row["accession"].upper() == accession)
    # The raw metrics are a stable baseline for --raw, even after repeated edits.
    baseline = copy.deepcopy(previous.get("metrics", {}).get("raw") or raw_audit)
    baseline["user_forbidden_site_hits"] = raw_audit["user_forbidden_site_hits"]
    updated_item["optimized_cds"] = {
        "path": output_relative, "file_sha256": hashlib.sha256(fasta).hexdigest(),
        "sequence_sha256": sha256_text(final), "length_nt": length,
        "report": {"path": report_path.relative_to(root).as_posix(), "file_sha256": hashlib.sha256(report_bytes).hexdigest()},
        "processing_mode": GC_MODE, "constraint_repair_applied": True,
        "optimization_skipped": False, "enforced_checks": report["enforced_checks"],
        "quality_checks_enforced": False, "gc_range_percent": [str(low), str(high)],
        "source_raw_sequence_sha256": sha256_text(raw),
        "metrics": {"raw": baseline, "final": final_audit, "changes": changes},
    }
    commit_cds_optimization(
        manifest_path=manifest_path, project_root=root, target=config.target_name,
        revision=int(manifest.get("revision", 0)), selection=updated,
        discard_sections=CDS_SELECTION_DOWNSTREAM_SECTIONS,
        files={output: fasta, report_path: report_bytes}, guards=guards,
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
