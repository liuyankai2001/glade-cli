"""Publish KEGG validation evidence without reselecting the route."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from functools import cache
import hashlib
import json
import math
from numbers import Real
import os
from pathlib import Path
from typing import Any

from filelock import FileLock, Timeout
import pandas as pd

from src.protein_to_cds.artifacts import ArtifactTransaction, write_bytes_atomic
from .solution_format import _resolve, _retropath, _values, compact_kegg_solution
from .store import manifest_update_lock


@cache
def _calculation_lock(path: str):
    return FileLock(path, timeout=0)


@contextmanager
def kegg_validation_lock(gap_dir: Path):
    directory = gap_dir / "gem_validation"
    directory.mkdir(parents=True, exist_ok=True)
    try:
        with _calculation_lock(os.path.normcase(str(directory.resolve() / ".kegg_validation.lock"))):
            yield
    except Timeout as exc:
        raise ValueError("当前项目和深度正在验证或写入路线，请稍后重试") from exc


def _raw_manifest(path: Path) -> dict:
    result = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(result, dict) or result.get("schema_version", "design_manifest.v1") != "design_manifest.v1":
        raise ValueError("manifest 格式无效，无法同步验证结果")
    return result


def _route_identity(solution: dict, manifest_path: Path):
    return (int(solution.get("solution_id") or 0), int(solution.get("expansion_depth") or 0),
            str(_resolve(solution.get("gap_dir") or "", manifest_path)),
            [(int(step["step_index"]), str(step.get("reaction_id") or ""), str(step.get("status") or "").lower(),
              str(step.get("direction") or ""), " ".join(str(step.get("equation") or "").split()),
              str(step.get("produced_compound_id") or ""), tuple(sorted(_values(step.get("precursor_compound_ids")))))
             for step in sorted(solution.get("steps", []), key=lambda item: int(item["step_index"]))])


def _hashes(paths):
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


@dataclass(frozen=True)
class ValidationSnapshot:
    manifest_path: Path
    target: str
    solution_id: int
    route_identity: Any
    inputs: dict[str, str]


def capture_validation_selection(config, *, gap_dir, selected_ids, mode, steps_frame, steps_sha256):
    if mode == "pooled":
        return None, "pooled 验证仅生成报告，不回写单路线状态"
    value = getattr(config, "manifest_output_path", None)
    if not value or not Path(value).is_file():
        return None, "尚无已选路线，仅生成验证报告"
    path = Path(value).expanduser().resolve()
    with manifest_update_lock(path):
        raw = _raw_manifest(path)
    solution = raw.get("solution")
    if not isinstance(solution, dict):
        return None, "尚无已选路线，仅生成验证报告"
    if _retropath(solution):
        return None, "本次不自动回写 RetroPath 路线"
    target = str(config.target_name).strip().upper()
    solution_id = int(solution.get("solution_id") or 0)
    if raw.get("target_compound_id") != target or solution_id not in selected_ids:
        return None, "本次未验证当前已选路线，仅生成验证报告"
    from .solution import _normalize_row, _renumber_steps_forward
    selected = steps_frame[steps_frame["solution_id"] == solution_id]
    rows = [_normalize_row(row) for row in selected.astype(object).where(pd.notna(selected), None).to_dict("records")]
    ordered, _ = _renumber_steps_forward(rows, target)
    expected = {"solution_id": solution_id, "expansion_depth": int(getattr(config, "depth", 0)),
                "gap_dir": str(gap_dir), "steps": ordered}
    identity = _route_identity(solution, path)
    if identity != _route_identity(expected, path):
        return None, "当前路线、深度或 gap 目录与验证输入不一致，结果不回写"
    inputs = [Path(config.model_path).expanduser().resolve(), Path(config.medium_path).expanduser().resolve(),
              gap_dir / "all_solution_steps.csv"]
    hashes = _hashes(inputs)
    hashes[str(gap_dir / "all_solution_steps.csv")] = steps_sha256
    return ValidationSnapshot(path, target, solution_id, identity, hashes), ""


def _matching_rows(frame, solution_id):
    if not {"validation_mode", "solution_ids", "validation_status"}.issubset(frame.columns):
        return frame.iloc[:0]
    from .solution import _parse_solution_ids
    indexes = []
    for index, row in frame.iterrows():
        try:
            if row["validation_mode"] == "per-solution" and _parse_solution_ids(row["solution_ids"]) == (solution_id,):
                indexes.append(index)
        except ValueError:
            continue
    return frame.loc[indexes]


def validation_evidence(summary, fluxes, *, solution_id, output_dir):
    """Create stable, selected-route reports from one in-memory result row."""
    matches = _matching_rows(summary, solution_id)
    if len(matches) != 1:
        raise ValueError("找不到已选路线的唯一单路线验证记录")
    row = matches.iloc[0].to_dict()
    if not str(row.get("validation_status") or "").strip():
        raise ValueError("单路线验证记录缺少验证状态")
    destination = output_dir / f"selected_solution_{solution_id}"
    summary_path = destination / "gem_validation_summary.csv"
    updates = {summary_path: matches.to_csv(index=False).encode("utf-8-sig")}
    files = {"summary_file": summary_path}
    if fluxes is not None:
        scoped = fluxes
        if not fluxes.empty:
            if "route_scope" not in fluxes.columns:
                scoped = None
            else:
                scope = str(row.get("route_scope") or f"solution_{solution_id}")
                scoped = fluxes[fluxes["route_scope"] == scope]
                if "validation_mode" in scoped.columns:
                    scoped = scoped[scoped["validation_mode"] == "per-solution"]
        if scoped is not None:
            if not len(scoped.columns):
                scoped = pd.DataFrame(columns=["route_scope", "validation_mode"])
            path = destination / "gem_validation_route_fluxes.csv"
            updates[path] = scoped.to_csv(index=False).encode("utf-8-sig")
            files["route_fluxes_file"] = path
    return row, files, updates


def project_validation(row, files, *, manifest_path, target):
    clean = {}
    for key, value in row.items():
        if isinstance(value, Real):
            value = float(value) if math.isfinite(float(value)) else None
        elif value is pd.NA:
            value = None
        clean[key] = value
    for key in ("fba_product_flux", "fba_growth_flux"):
        try:
            value = float(clean.get(key))
            clean[key] = value if math.isfinite(value) else None
        except (ValueError, TypeError):
            clean[key] = None
    from .solution import _KEGG_VALIDATION_FAILED_WARNING
    status = str(clean.get("validation_status") or "")
    if status and status != "not_run" and not status.startswith("PASS_"):
        issue = str(clean.get("issues") or "").strip()
        if _KEGG_VALIDATION_FAILED_WARNING not in issue:
            clean["issues"] = " | ".join(filter(None, (issue, _KEGG_VALIDATION_FAILED_WARNING)))
    clean.update({key: str(path) for key, path in files.items()})
    return compact_kegg_solution({"validation": clean}, manifest_path=manifest_path,
                                 target_compound_id=target).get("validation", {})


def publish_evidence(updates):
    for path, content in updates.items():
        if not path.exists() or path.read_bytes() != content:
            write_bytes_atomic(path, content)


def sync_validation_selection(snapshot, *, skip_message, summary, fluxes, output_dir):
    if snapshot is None:
        return {"manifest_updated": False, "manifest_update_message": skip_message}
    path = snapshot.manifest_path
    try:
        with manifest_update_lock(path):
            if not path.exists():
                return {"manifest_updated": False, "manifest_update_message": "验证期间清单已移除，结果未回写"}
            raw = _raw_manifest(path)
            solution = raw.get("solution")
            try:
                stable = (raw.get("target_compound_id") == snapshot.target and isinstance(solution, dict)
                          and not _retropath(solution) and _route_identity(solution, path) == snapshot.route_identity
                          and _hashes(map(Path, snapshot.inputs)) == snapshot.inputs)
            except (OSError, ValueError, KeyError, TypeError):
                stable = False
            if not stable:
                return {"manifest_updated": False, "manifest_update_message": "验证期间路线或模型、培养基输入已变化，报告已生成但结果未回写"}
            row, files, updates = validation_evidence(summary, fluxes, solution_id=snapshot.solution_id, output_dir=output_dir)
            validation = project_validation(row, files, manifest_path=path, target=snapshot.target)
            changed = solution.get("validation") != validation
            with ArtifactTransaction((path, *updates)) as transaction:
                publish_evidence(updates)
                if changed:
                    solution["validation"] = validation
                    from . import store
                    store._write_json_atomic(path, raw)
                transaction.commit()
            return {"manifest_updated": changed, "manifest_update_message": "已自动更新当前路线的验证结果" if changed else "验证结果相同，清单未重写"}
    except (OSError, ValueError) as exc:
        raise ValueError(f"验证报告已生成，但清单回写失败，原清单已保留：{exc}") from exc


def existing_validation_evidence(directory, solution_id):
    """Freeze an existing report when selecting a route after validation."""
    path = directory / "gem_validation_summary.csv"
    if not path.is_file():
        return None
    try:
        summary = pd.read_csv(path, keep_default_na=False)
        matches = _matching_rows(summary, solution_id).tail(1)
        if matches.empty or not str(matches.iloc[0]["validation_status"]).strip():
            return None
        flux_path = directory / "gem_validation_route_fluxes.csv"
        fluxes = pd.read_csv(flux_path, keep_default_na=False) if flux_path.is_file() else None
        return validation_evidence(matches, fluxes, solution_id=solution_id, output_dir=directory)
    except (OSError, ValueError, KeyError, pd.errors.EmptyDataError):
        return None
