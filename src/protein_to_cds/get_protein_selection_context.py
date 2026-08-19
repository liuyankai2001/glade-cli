from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Literal

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.runtime.monitor import monitor
from src.tools.common.session_paths import design_manifest_file, session_dir as resolve_session_dir
from src.tools.enzyme_system_selection_tools.protein_selection_gate import (
    evaluate_protein_selection_gate,
)


class GetProteinSelectionContextArgs(BaseModel):
    top_candidates_per_step: int = Field(default=3, ge=0, le=50, description="每个 step 返回的候选蛋白数量")
    detail: Literal["compact", "full"] = Field(
        default="compact",
        description="返回详细程度：compact 精简输出；full 返回完整候选和 CDS 上下文",
    )


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _session_root(session_dir: str | None, manifest_path: str | None) -> Path:
    if session_dir:
        return Path(session_dir).resolve()
    if manifest_path:
        path = Path(manifest_path).resolve()
        if path.name == "design_manifest.json" and path.parent.name == "outputs":
            return path.parent.parent.resolve()
    return resolve_session_dir()


def _safe_path(session_root: Path, path_value: str | Path) -> Path:
    path = Path(path_value)
    if not path.is_absolute():
        path = session_root / path
    path = path.resolve()
    session_root = session_root.resolve()
    if path != session_root and session_root not in path.parents:
        raise ValueError(f"path escapes session_dir: {path}")
    return path


def _relpath(session_root: Path, path: Path) -> str:
    return path.resolve().relative_to(session_root.resolve()).as_posix()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return data


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        return list(csv.DictReader(file))


def _split_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return _unique([str(item) for item in value])
    text = str(value).strip()
    if not text:
        return []
    return _unique([
        part.strip()
        for chunk in text.split("|")
        for part in chunk.split(";")
        if part.strip()
    ])


def _unique(values: list[str]) -> list[str]:
    seen = set()
    result = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            result.append(text)
    return result


def _number(value: Any) -> int | float | str | bool | None:
    text = str(value or "").strip()
    if not text:
        return None
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        number = float(text)
    except ValueError:
        return text
    return int(number) if number.is_integer() else number


def _to_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def _candidate_file(manifest: dict[str, Any], session_root: Path, key: str) -> Path | None:
    protein_selection = manifest.get("protein_selection")
    if not isinstance(protein_selection, dict):
        return None
    files = protein_selection.get("files")
    if not isinstance(files, dict):
        return None
    value = files.get(key)
    if not value:
        return None
    return _safe_path(session_root, value)


def _pathway_steps_file(manifest: dict[str, Any], session_root: Path) -> Path | None:
    pathway_selection = manifest.get("pathway_selection")
    if isinstance(pathway_selection, dict):
        evidence_files = pathway_selection.get("evidence_files")
        if isinstance(evidence_files, dict):
            steps_csv = evidence_files.get("steps_csv")
            if isinstance(steps_csv, dict) and steps_csv.get("path"):
                return _safe_path(session_root, steps_csv["path"])

    solution = manifest.get("solution")
    if isinstance(solution, dict):
        files = solution.get("files")
        if isinstance(files, dict):
            path_value = files.get("all_solution_steps_csv") or files.get("steps_csv")
            if path_value:
                return _safe_path(session_root, path_value)
    return None


def _selected_solution_id(manifest: dict[str, Any]) -> int | None:
    pathway_selection = manifest.get("pathway_selection")
    if isinstance(pathway_selection, dict):
        value = pathway_selection.get("selected_solution_id")
        if value is not None:
            return _to_int(value)

    solution = manifest.get("solution")
    if isinstance(solution, dict):
        value = solution.get("solution_id")
        if value is not None:
            return _to_int(value)

    protein_selection = manifest.get("protein_selection")
    if isinstance(protein_selection, dict):
        value = protein_selection.get("selected_solution_id")
        if value is not None:
            return _to_int(value)
    return None


def _pathway_steps(manifest: dict[str, Any], session_root: Path, solution_id: int | None) -> tuple[list[dict[str, Any]], str | None]:
    steps_file = _pathway_steps_file(manifest, session_root)
    if steps_file is not None:
        rows = _read_csv(steps_file)
        if solution_id is not None:
            rows = [row for row in rows if _to_int(row.get("solution_id")) == solution_id]
        return [{key: _number(value) for key, value in row.items()} for row in rows], _relpath(session_root, steps_file)

    solution = manifest.get("solution")
    if isinstance(solution, dict) and isinstance(solution.get("steps"), list):
        return [step for step in solution["steps"] if isinstance(step, dict)], None

    return [], None


def _candidate_summary(row: dict[str, str], recommended_accessions: set[str]) -> dict[str, Any]:
    accession = str(row.get("accession", "")).strip()
    payload = {
        "accession": accession,
        "entry_name": row.get("entry_name", ""),
        "protein_name": row.get("protein_name", ""),
        "organism_name": row.get("organism_name", ""),
        "organism_id": _number(row.get("organism_id")),
        "reviewed": _number(row.get("reviewed")),
        "length": _number(row.get("length")),
        "role": row.get("role", "main") or "main",
        "accessory_type": row.get("accessory_type", ""),
        "accessory_for_step_index": _number(row.get("accessory_for_step_index")),
        "ec_number": row.get("ec_number", ""),
        "score": _number(row.get("score")),
        "gene_names": _split_list(row.get("gene_names")),
        "warnings": _split_list(row.get("warnings")),
        "recommended_for_route": accession in recommended_accessions,
    }
    return {key: value for key, value in payload.items() if value not in (None, "", [])}


def _compact_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    payload = {
        "accession": candidate.get("accession"),
        "protein_name": candidate.get("protein_name", ""),
        "organism_name": candidate.get("organism_name", ""),
        "role": candidate.get("role", "main"),
        "accessory_type": candidate.get("accessory_type", ""),
        "reviewed": candidate.get("reviewed"),
        "length": candidate.get("length"),
        "score": candidate.get("score"),
        "recommended_for_route": candidate.get("recommended_for_route"),
    }
    warnings = candidate.get("warnings") or []
    if warnings:
        payload["warnings"] = warnings
    return payload


def _compact_step(step: dict[str, Any]) -> dict[str, Any]:
    return {
        "step_index": step.get("step_index"),
        "reaction_id": step.get("reaction_id", ""),
        "reaction_name": step.get("reaction_name", ""),
        "ec_numbers": step.get("ec_numbers", []),
        "candidate_count": step.get("candidate_count", 0),
        "recommended_accessions": step.get("recommended_accessions", []),
        "candidates": [
            _compact_candidate(candidate)
            for candidate in step.get("candidates", [])
            if isinstance(candidate, dict)
        ],
    }


def _compact_selected_cds(item: dict[str, Any]) -> dict[str, Any]:
    protein = item.get("protein") if isinstance(item.get("protein"), dict) else {}
    optimized = item.get("optimized_cds") if isinstance(item.get("optimized_cds"), dict) else {}
    sequence_file = optimized.get("sequence_file") if isinstance(optimized.get("sequence_file"), dict) else {}
    payload = {
        "cds_id": item.get("cds_id"),
        "step_index": item.get("step_index"),
        "reaction_id": item.get("reaction_id"),
        "ec_number": item.get("ec_number"),
        "accession": protein.get("accession"),
        "protein_name": protein.get("protein_name", ""),
        "organism_name": protein.get("organism_name", ""),
        "sequence_file": sequence_file.get("path"),
        "length_nt": optimized.get("length_nt"),
        "gc_percent": optimized.get("gc_percent"),
    }
    return {key: value for key, value in payload.items() if value not in (None, "", [])}


def _recommended_accessions(manifest: dict[str, Any]) -> list[str]:
    protein_selection = manifest.get("protein_selection")
    if not isinstance(protein_selection, dict):
        return []
    design = protein_selection.get("recommended_design")
    if not isinstance(design, dict):
        return []
    return _split_list(design.get("selected_accessions"))


def _recommended_proteins(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    protein_selection = manifest.get("protein_selection")
    if not isinstance(protein_selection, dict):
        return []
    design = protein_selection.get("recommended_design")
    if not isinstance(design, dict):
        return []
    proteins = design.get("selected_proteins")
    if isinstance(proteins, list):
        return [item for item in proteins if isinstance(item, dict)]
    return []


def _protein_selection_gate(
    manifest: dict[str, Any],
    heterologous_steps: list[dict[str, Any]],
) -> dict[str, Any]:
    required = {
        _to_int(step.get("step_index"))
        for step in heterologous_steps
        if _to_int(step.get("step_index"))
    }
    # Figure-3's immutable v2.5 benchmark uses a deliberately synthetic object
    # without a route solution. Preserve its historical diagnostic semantics,
    # while every real manifest (which has ``solution``) uses the v2.6 gate and
    # rejects v2.5 for CDS hand-off.
    if (
        not isinstance(manifest.get("solution"), dict)
        and not str(manifest.get("schema_version") or "").strip()
    ):
        selection = manifest.get("protein_selection")
        design = selection.get("recommended_design") if isinstance(selection, dict) else None
        if not isinstance(design, dict):
            return {
                "can_proceed_to_cds": False,
                "status": "missing_protein_selection" if not isinstance(selection, dict) else "missing_recommended_design",
                "verified_step_indexes": [],
                "unverified_step_indexes": sorted(required),
                "errors": ["synthetic context lacks a complete protein selection"],
            }
        schema = str(selection.get("schema_version") or design.get("schema_version") or "")
        verified = {
            _to_int(value)
            for value in _split_list(design.get("verified_step_indexes"))
            if _to_int(value)
        }
        blocking = design.get("blocking_issues", [])
        blocking = blocking if isinstance(blocking, list) else []
        accessory_requirements = selection.get("accessory_requirements", [])
        accessory_requirements = accessory_requirements if isinstance(accessory_requirements, list) else []
        selected_proteins = design.get("selected_proteins", [])
        selected_proteins = selected_proteins if isinstance(selected_proteins, list) else []
        string_required = bool(accessory_requirements) or any(
            isinstance(item, dict) and item.get("role") == "required_main_component"
            for item in selected_proteins
        )
        string_summary = design.get("string_evidence_summary", {})
        string_summary = string_summary if isinstance(string_summary, dict) else {}
        string_ok = not string_required or (
            string_summary.get("completion_attempted") is True
            and string_summary.get("all_required_systems_verified") is True
        )
        can_proceed = (
            schema == "protein_selection.v2.5"
            and str(design.get("status") or "") in {"complete", "complete_with_risks"}
            and verified == required
            and not blocking
            and string_ok
            and bool(_split_list(design.get("selected_accessions")))
        )
        return {
            "can_proceed_to_cds": can_proceed,
            "status": str(design.get("status") or "") if schema == "protein_selection.v2.5" else "requires_reselection",
            "schema_version": schema,
            "verified_step_indexes": sorted(verified),
            "unverified_step_indexes": sorted(required - verified),
            "blocking_issue_count": len(blocking),
            "string_gate_ok": string_ok,
            "legacy_synthetic_context": True,
            "errors": [] if can_proceed else ["legacy synthetic benchmark gate is blocked"],
        }

    result = evaluate_protein_selection_gate(manifest)
    verified = {
        _to_int(value)
        for value in result.get("verified_step_indexes", [])
        if _to_int(value)
    }
    string_errors = {
        "STRING completion was not run for network-required components",
        "network-required component evidence is incomplete",
        "selected network-required accessory system is not auto-verified",
    }
    result["unverified_step_indexes"] = sorted(required - verified)
    result["string_gate_ok"] = not bool(string_errors & set(result.get("errors", [])))
    return result


def _route_summary(manifest: dict[str, Any]) -> dict[str, Any]:
    pathway_selection = manifest.get("pathway_selection")
    if isinstance(pathway_selection, dict):
        return {
            "source": pathway_selection.get("source"),
            "target_compound_id": pathway_selection.get("target_compound_id"),
            "selected_solution_id": pathway_selection.get("selected_solution_id"),
            "gap_dir": pathway_selection.get("gap_dir"),
        }

    solution = manifest.get("solution")
    if isinstance(solution, dict):
        summary = solution.get("summary") if isinstance(solution.get("summary"), dict) else {}
        return {
            "source": solution.get("source"),
            "target_compound_id": summary.get("target_compound_id"),
            "selected_solution_id": solution.get("solution_id"),
            "gap_dir": solution.get("gap_dir"),
        }
    return {}


@tool(args_schema=GetProteinSelectionContextArgs)
def get_protein_selection_context(
    top_candidates_per_step: int = 3,
    detail: Literal["compact", "full"] = "compact",
) -> str:
    """
    读取当前 protein_selection 和候选蛋白上下文，供选择序列和密码子优化使用。
    
    调用时机：用户要求根据已选路线继续做 CDS、查看蛋白候选或开始优化。
    返回：默认精简的 solution 摘要、推荐 accession、步骤候选摘要和已选 CDS 摘要；detail=full 返回完整上下文。
    限制：只读；不搜索序列、不写 cds_selection。
    """

    tool_name = "get_protein_selection_context"
    monitor.report_start(tool_name, {"top_candidates_per_step": top_candidates_per_step})
    try:
        session_root = _session_root(None, None)
        manifest_file = design_manifest_file()
        manifest = _read_json(manifest_file)
        solution_id = _selected_solution_id(manifest)
        steps, steps_path = _pathway_steps(manifest, session_root, solution_id)

        heterologous_steps = [
            step
            for step in steps
            if str(step.get("status", "")).strip().lower() == "heterologous"
        ]
        heterologous_steps.sort(key=lambda step: _to_int(step.get("step_index")))

        step_candidates_path = _candidate_file(manifest, session_root, "step_protein_candidates_csv")
        candidate_rows = _read_csv(step_candidates_path) if step_candidates_path else []
        candidates_by_step: dict[int, list[dict[str, str]]] = {}
        for row in candidate_rows:
            step_index = _to_int(row.get("step_index"))
            if step_index:
                candidates_by_step.setdefault(step_index, []).append(row)

        gate = _protein_selection_gate(manifest, heterologous_steps)
        recommended = _recommended_accessions(manifest) if gate["can_proceed_to_cds"] else []
        recommended_proteins = _recommended_proteins(manifest) if gate["can_proceed_to_cds"] else []
        recommended_set = set(recommended)
        selected_cds = []
        cds_selection = manifest.get("cds_selection")
        if isinstance(cds_selection, dict) and isinstance(cds_selection.get("selected_cds"), list):
            selected_cds = [
                item
                for item in cds_selection["selected_cds"]
                if isinstance(item, dict)
            ]

        step_contexts = []
        missing_candidate_step_indexes = []
        for step in heterologous_steps:
            step_index = _to_int(step.get("step_index"))
            rows = candidates_by_step.get(step_index, [])
            rows.sort(key=lambda row: float(row.get("score") or 0.0), reverse=True)
            if not rows:
                missing_candidate_step_indexes.append(step_index)

            candidates = [
                _candidate_summary(row, recommended_set)
                for row in rows[:top_candidates_per_step]
            ]
            step_contexts.append({
                "step_index": step_index,
                "reaction_id": step.get("reaction_id", ""),
                "reaction_name": step.get("reaction_name", ""),
                "produced_compound_id": step.get("produced_compound_id", ""),
                "produced_compound_name": step.get("produced_compound_name", ""),
                "ec_numbers": _split_list(step.get("enzyme_ecs")),
                "ko_ids": _split_list(step.get("ko_ids")),
                "candidate_count": len(rows),
                "recommended_accessions": [
                    accession
                    for accession in recommended
                    if any(str(row.get("accession", "")).strip() == accession for row in rows)
                ],
                "candidates": candidates,
            })

        protein_selection = manifest.get("protein_selection")
        protein_selection_summary = protein_selection if isinstance(protein_selection, dict) else {}

        full_result = {
            "ok": True,
            "manifest_path": _relpath(session_root, manifest_file),
            "revision": manifest.get("revision"),
            "solution_summary": _route_summary(manifest),
            "protein_selection": {
                "available": isinstance(protein_selection, dict),
                "source": protein_selection_summary.get("source"),
                "selected_solution_id": protein_selection_summary.get("selected_solution_id"),
                "chassis_key": protein_selection_summary.get("chassis_key"),
                "recommended_accessions": recommended,
                "recommended_proteins": recommended_proteins,
                "recommended_design": protein_selection_summary.get("recommended_design", {}),
                "gate": gate,
            },
            "evidence_files": {
                "steps_csv": steps_path,
                "step_protein_candidates_csv": (
                    _relpath(session_root, step_candidates_path)
                    if step_candidates_path and step_candidates_path.exists()
                    else None
                ),
            },
            "heterologous_step_count": len(step_contexts),
            "missing_candidate_step_indexes": missing_candidate_step_indexes,
            "steps": step_contexts,
            "selected_cds": selected_cds,
            "next_actions": [
                "Use recommended_accessions or candidates[*].accession with search_protein_sequence.",
                "Then run codon_optimization and write_selected_optimized_cds_to_manifest for user-selected CDS.",
            ],
        }
        compact_result = {
            "ok": True,
            "revision": manifest.get("revision"),
            "solution_summary": _route_summary(manifest),
            "protein_selection": {
                "available": isinstance(protein_selection, dict),
                "source": protein_selection_summary.get("source"),
                "selected_solution_id": protein_selection_summary.get("selected_solution_id"),
                "chassis_key": protein_selection_summary.get("chassis_key"),
                "recommended_accessions": recommended,
                "recommended_proteins": recommended_proteins,
                "recommended_design": protein_selection_summary.get("recommended_design", {}),
                "gate": gate,
            },
            "evidence_files": {
                "steps_csv": steps_path,
                "step_protein_candidates_csv": (
                    _relpath(session_root, step_candidates_path)
                    if step_candidates_path and step_candidates_path.exists()
                    else None
                ),
            },
            "heterologous_step_count": len(step_contexts),
            "missing_candidate_step_indexes": missing_candidate_step_indexes,
            "steps": [_compact_step(step) for step in step_contexts],
            "selected_cds": [
                _compact_selected_cds(item)
                for item in selected_cds
                if isinstance(item, dict)
            ],
        }
        result = full_result if detail == "full" else compact_result
        monitor.report_end(tool_name, {"heterologous_step_count": len(step_contexts)})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))
