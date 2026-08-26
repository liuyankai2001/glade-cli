"""Standalone P11.1 RetroPath benchmark runner.

Usage::

    python -m src.pathway_analyze.retropath_benchmark validate --cases CASES.json
    python -m src.pathway_analyze.retropath_benchmark run --cases CASES.json
    python -m src.pathway_analyze.retropath_benchmark report --run RUN_DIR

The module intentionally does not register a regular GLADE CLI command.  It
isolates all runtime outputs and caches below the requested benchmark run
directory while reusing the production P2--P10 functions.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sqlite3
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

from src.config.run_config import RunConfig
from src.main_protein_selection.select_main_enzymes import (
    run_main_protein_selection,
)
from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir
from src.pathway_analyze.retropath_benchmark_models import (
    BENCHMARK_RUN_SCHEMA,
    BENCHMARK_TASK_SCHEMA,
    BenchmarkCase,
    BenchmarkDataset,
    canonical_sha256,
    load_benchmark_dataset,
    repository_root,
    sha256_file,
)
from src.pathway_analyze.retropath_benchmark_report import (
    generate_benchmark_report,
    score_core_artifacts,
    score_enzyme_artifacts,
)
from src.pathway_analyze.retropath_client import (
    RetroPathClientError,
    RetroPathHttpClient,
)
from src.pathway_analyze.retropath_gem_validation import (
    validate_retropath_candidates,
)
from src.pathway_analyze.retropath_pipeline import (
    PIPELINE_RESULT_FILE_NAME,
    RetroPathPipelineError,
    run_retropath_pipeline,
)
from src.write_manifest.solution import write_solution


DEFAULT_CASES_PATH = (
    repository_root()
    / "docs"
    / "retropath_benchmark"
    / "fixtures"
    / "p11_1_hidden_reactions.json"
)
DEFAULT_OUTPUT_ROOT = (
    repository_root() / "tests" / ".retropath_benchmark_runtime"
)
RUN_MANIFEST_FILE_NAME = "benchmark_run_manifest.json"
TASK_RESULT_FILE_NAME = "task_result.json"
BENCHMARK_PROFILES = ("controlled", "full_a0")
TERMINAL_TASK_STATUSES = frozenset({"completed", "failed", "skipped"})
_BIOLOGICAL_DEFAULT_FIELDS = (
    "max_steps",
    "topx",
    "dmin",
    "dmax",
    "max_candidates",
    "top_k",
    "run_fva",
)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat()


def _atomic_write_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(text)
            temporary = Path(handle.name)
        temporary.replace(path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return sha256_file(path)


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> str:
    return _atomic_write_text(
        path,
        json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        + "\n",
    )


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid benchmark JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"benchmark JSON must be an object: {path}")
    return value


def _git_provenance(root: Path) -> dict[str, Any]:
    def command(*args: str) -> str:
        try:
            result = subprocess.run(
                ["git", *args],
                cwd=root,
                check=True,
                capture_output=True,
                text=True,
                encoding="utf-8",
            )
        except (OSError, subprocess.CalledProcessError):
            return ""
        return result.stdout.strip()

    status = command("status", "--porcelain")
    return {
        "commit": command("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status_sha256": canonical_sha256(status.splitlines()),
    }


@lru_cache(maxsize=1)
def _source_tree_record() -> dict[str, Any]:
    root = repository_root()
    files: list[Path] = []
    for relative in (
        "src/pathway_analyze",
        "src/main_protein_selection",
        "src/write_manifest",
    ):
        files.extend(
            path
            for path in (root / relative).glob("*.py")
            if not path.name.startswith("retropath_benchmark")
        )
    records = {
        str(path.resolve().relative_to(root)): sha256_file(path)
        for path in sorted(set(files))
        if path.is_file()
    }
    return {
        "sha256": canonical_sha256(records),
        "files": records,
    }


@lru_cache(maxsize=1)
def _runner_source_record() -> dict[str, Any]:
    root = repository_root()
    files = sorted(Path(__file__).parent.glob("retropath_benchmark*.py"))
    records = {
        str(path.resolve().relative_to(root)): sha256_file(path)
        for path in files
        if path.is_file()
    }
    return {"sha256": canonical_sha256(records), "files": records}


def validate_gold_against_mnxref(dataset: BenchmarkDataset) -> dict[str, Any]:
    """Verify that case identities and reaction orientation exist in the pinned index."""

    index = dataset.resources["mnxref_index"].path
    connection = sqlite3.connect(f"file:{index.as_posix()}?mode=ro", uri=True)
    try:
        audits: list[dict[str, Any]] = []
        for case in dataset.cases:
            target_mapping = connection.execute(
                "SELECT 1 FROM chemical_xrefs WHERE mnxm_id = ? AND lower(xref) = ?",
                (case.target_mnxm_id, f"kegg:{case.target_kegg_id.lower()}"),
            ).fetchone()
            if target_mapping is None:
                raise ValueError(
                    f"{case.case_id} target KEGG/MNXM mapping is absent from MNXref"
                )
            gold_mnxr = case.gold_mnxr_ids
            if len(gold_mnxr) != 1:
                raise ValueError("P11.1 pilot currently requires one-step gold routes")
            reaction = connection.execute(
                """
                SELECT balanced, transport, parse_status
                FROM reactions WHERE mnxr_id = ?
                """,
                (gold_mnxr[0],),
            ).fetchone()
            if reaction is None or not reaction[0] or reaction[1] or reaction[2] != "ok":
                raise ValueError(
                    f"{case.case_id} gold MNXref reaction is not balanced/non-transport/parseable"
                )
            reaction_mapping = connection.execute(
                "SELECT 1 FROM reaction_xrefs WHERE mnxr_id = ? AND lower(xref) = ?",
                (gold_mnxr[0], f"kegg:{case.gold_reaction_ids[0].lower()}"),
            ).fetchone()
            if reaction_mapping is None:
                raise ValueError(
                    f"{case.case_id} gold KEGG/MNXref reaction mapping is absent"
                )
            target_sides = {
                row[0]
                for row in connection.execute(
                    """
                    SELECT side FROM reaction_terms
                    WHERE mnxr_id = ? AND mnxm_id = ?
                    """,
                    (gold_mnxr[0], case.target_mnxm_id),
                )
            }
            if len(target_sides) != 1:
                raise ValueError(
                    f"{case.case_id} target does not occur on exactly one reaction side"
                )
            target_side = next(iter(target_sides))
            opposite = "right" if target_side == "left" else "left"
            opposite_mnxm = {
                row[0]
                for row in connection.execute(
                    "SELECT mnxm_id FROM reaction_terms WHERE mnxr_id = ? AND side = ?",
                    (gold_mnxr[0], opposite),
                )
            }
            if not set(case.controlled_sink_mnxm_ids) <= opposite_mnxm:
                raise ValueError(
                    f"{case.case_id} controlled MNXM sinks are not on the gold precursor side"
                )
            for sink in case.controlled_sink_kegg_ids:
                mapping = connection.execute(
                    """
                    SELECT mnxm_id FROM chemical_xrefs
                    WHERE lower(xref) = ?
                    """,
                    (f"kegg:{sink.lower()}",),
                ).fetchall()
                mapped = {row[0] for row in mapping}
                if not mapped.intersection(case.controlled_sink_mnxm_ids):
                    raise ValueError(
                        f"{case.case_id} sink {sink} does not match a declared precursor MNXM"
                    )
            rule_count = connection.execute(
                """
                SELECT COUNT(*) FROM rule_templates
                WHERE mnxr_id = ? AND main_mnxm_id = ?
                """,
                (gold_mnxr[0], case.target_mnxm_id),
            ).fetchone()[0]
            if rule_count < 1:
                raise ValueError(
                    f"{case.case_id} has no pinned RR02 rule for the target identity"
                )
            audits.append(
                {
                    "case_id": case.case_id,
                    "target_side": target_side,
                    "gold_precursor_side": opposite,
                    "rule_template_count": int(rule_count),
                    "controlled_sink_count": len(case.controlled_sink_kegg_ids),
                }
            )
    finally:
        connection.close()
    return {
        "ok": True,
        "benchmark_id": dataset.benchmark_id,
        "case_count": len(dataset.cases),
        "cases": audits,
    }


def _service_health(
    service_url: str,
    *,
    request_timeout_seconds: float = 30.0,
    get_attempts: int = 3,
) -> dict[str, Any]:
    with RetroPathHttpClient(
        service_url,
        request_timeout_seconds=request_timeout_seconds,
        get_attempts=get_attempts,
    ) as client:
        return client.health().to_dict()


_RUNTIME_IDENTITY_FIELDS = (
    "service_version",
    "wrapper_version",
    "wrapper_reported_version",
    "workflow_version",
    "knime_version",
    "rdkit_plugin_version",
    "rules_version",
    "rules_sha256",
    "worker_concurrency",
)


def _runtime_identity(health: Mapping[str, Any]) -> dict[str, Any]:
    """Exclude dynamic readiness/queue fields from the runtime version lock."""

    return {field: health.get(field) for field in _RUNTIME_IDENTITY_FIELDS}


def _biological_defaults(value: Mapping[str, Any]) -> dict[str, Any]:
    return {field: value.get(field) for field in _BIOLOGICAL_DEFAULT_FIELDS}


def _resources_compatible(
    manifest_resources: Any,
    dataset: BenchmarkDataset,
) -> bool:
    if not isinstance(manifest_resources, Mapping):
        return False
    return all(
        isinstance(manifest_resources.get(name), Mapping)
        and manifest_resources[name].get("sha256")
        in {resource.sha256, *resource.compatible_sha256}
        for name, resource in dataset.resources.items()
    )


def _render_controlled_a0(case: BenchmarkCase) -> str:
    output = []
    columns = ("source", "met_id", "met_name", "compartment", "kegg_id")
    buffer = tempfile.SpooledTemporaryFile(mode="w+", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(buffer, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for index, compound_id in enumerate(case.controlled_sink_kegg_ids, start=1):
            writer.writerow(
                {
                    "source": "benchmark_controlled",
                    "met_id": f"benchmark_sink_{index}",
                    "met_name": compound_id,
                    "compartment": "c",
                    "kegg_id": compound_id,
                }
            )
        buffer.seek(0)
        output.append(buffer.read())
    finally:
        buffer.close()
    return "".join(output)


def _prepare_a0(
    dataset: BenchmarkDataset,
    case: BenchmarkCase,
    profile: str,
    path: Path,
) -> str:
    if profile == "controlled":
        return _atomic_write_text(path, _render_controlled_a0(case))
    if profile != "full_a0":
        raise ValueError(f"unsupported benchmark profile: {profile}")
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(dataset.resources["a0_snapshot"].path, path)
    return sha256_file(path)


def _benchmark_config(
    dataset: BenchmarkDataset,
    case: BenchmarkCase,
    task_dir: Path,
    shared_cache: Path,
    *,
    service_url: str,
) -> RunConfig:
    config = RunConfig(target_name=case.target_kegg_id)
    config.project_output_path = task_dir
    config.chassis_output_path = task_dir / "chassis_result"
    config.chassis_producible_csv = (
        config.chassis_output_path / "producible_kegg_compounds.csv"
    )
    config.chassis_metabolites_summary_csv = (
        config.chassis_output_path / "analyze_chassis_metabolites_summary.csv"
    )
    config.gap_output_path = task_dir / f"kegg_gap_{case.target_kegg_id}"
    config.validation_output_path = config.gap_output_path / "gem_validation"
    config.manifest_output_path = task_dir / "design_manifest.json"
    config.cache_dir = shared_cache
    config.model_path = dataset.resources["model"].path
    config.medium_path = dataset.resources["medium"].path
    config.data_dir = repository_root() / "data"
    config.retropath_rules_path = dataset.resources["rr02"].path
    config.retropath_service_url = service_url
    config.retropath = True
    config.depth = 0
    config.retropath_max_steps = dataset.defaults.max_steps
    config.retropath_topx = dataset.defaults.topx
    config.retropath_dmin = dataset.defaults.dmin
    config.retropath_dmax = dataset.defaults.dmax
    config.retropath_max_candidates = dataset.defaults.max_candidates
    config.retropath_candidates = []
    config.validation_mode = "per"
    config.validation_cofactor_mode = "strict"
    config.validation_skip_fva = not dataset.defaults.run_fva
    config.retropath_request_timeout_seconds = (
        dataset.defaults.request_timeout_seconds
    )
    config.retropath_get_attempts = dataset.defaults.get_attempts
    config.retropath_wait_timeout_seconds = dataset.defaults.wait_timeout_seconds
    return config


def _task_fingerprint(
    dataset: BenchmarkDataset,
    case: BenchmarkCase,
    profile: str,
    *,
    service_health: Mapping[str, Any],
    with_enzymes: bool,
) -> str:
    return canonical_sha256(
        {
            "schema_version": BENCHMARK_TASK_SCHEMA,
            "dataset_sha256": dataset.dataset_sha256,
            "case": case.to_dict(),
            "profile": profile,
            "defaults": dataset.defaults.to_dict(),
            "resources": {
                name: resource.sha256
                for name, resource in sorted(dataset.resources.items())
            },
            "service_identity": _runtime_identity(service_health),
            "with_enzymes": with_enzymes and profile == "full_a0",
            "source_tree_sha256": _source_tree_record()["sha256"],
        }
    )


def _artifact_record(path: Path, *, root: Path) -> dict[str, str]:
    return {
        "path": str(path.resolve().relative_to(root.resolve())),
        "sha256": sha256_file(path),
    }


def _existing_task_valid(
    record: Mapping[str, Any],
    *,
    run_dir: Path,
    fingerprint: str,
    retry_failed: bool = False,
) -> bool:
    relative = Path(str(record.get("result_path") or ""))
    result_path = relative if relative.is_absolute() else run_dir / relative
    if not result_path.is_file() or sha256_file(result_path) != record.get("sha256"):
        return False
    try:
        result = _read_json(result_path)
    except ValueError:
        return False
    status = str(result.get("status") or "")
    if status not in TERMINAL_TASK_STATUSES:
        return False
    if retry_failed and status != "completed":
        return False
    if result.get("task_fingerprint") != fingerprint:
        return False
    artifacts = result.get("artifacts", {})
    if not isinstance(artifacts, Mapping):
        return False
    for artifact in artifacts.values():
        if not isinstance(artifact, Mapping):
            return False
        path = Path(str(artifact.get("path") or ""))
        path = path if path.is_absolute() else run_dir / path
        if not path.is_file() or sha256_file(path) != artifact.get("sha256"):
            return False
    return True


def _task_from_record(
    record: Mapping[str, Any],
    *,
    run_dir: Path,
) -> dict[str, Any]:
    relative = Path(str(record.get("result_path") or ""))
    result_path = relative if relative.is_absolute() else run_dir / relative
    return _read_json(result_path)


def _is_infrastructure_failure(task: Mapping[str, Any] | None) -> bool:
    if not isinstance(task, Mapping) or task.get("status") != "failed":
        return False
    error = task.get("error", {})
    return bool(
        task.get("outcome") == "infrastructure_error"
        or (
            isinstance(error, Mapping)
            and error.get("stage") == "infrastructure_error"
        )
    )


def _classify_pipeline_error(error: RetroPathPipelineError) -> tuple[str, str]:
    if error.status in {
        "retropath_service_unavailable",
        "retropath_timeout",
        "retropath_execution_failed",
    }:
        return "infrastructure_error", error.status
    if error.status in {
        "retropath_input_invalid",
        "retropath_expansion_missing",
        "retropath_rules_missing",
        "retropath_configuration_invalid",
    }:
        return "data_error", error.status
    return "pipeline_error", error.status


def _outcome(metrics: Mapping[str, Any]) -> str:
    if metrics.get("formal_exact_gold_rank") is not None:
        return "formal_exact_recovered"
    if metrics.get("strict_gem_gold_rank") is not None:
        return "strict_gem_nonexact"
    if metrics.get("balanced_gold_rank") is not None:
        return "balanced_gold_only"
    if metrics.get("gold_template_rank") is not None:
        return "gold_template_only"
    if int(metrics.get("candidate_count") or 0) == 0:
        return "no_candidate"
    return "gold_not_recovered"


def _run_task(
    dataset: BenchmarkDataset,
    case: BenchmarkCase,
    profile: str,
    *,
    run_dir: Path,
    service_url: str,
    service_health: Mapping[str, Any],
    with_enzymes: bool,
    force_retropath: bool = False,
) -> dict[str, Any]:
    task_dir = run_dir / "cases" / case.case_id / profile
    task_dir.mkdir(parents=True, exist_ok=True)
    shared_cache = run_dir / "shared_cache"
    config = _benchmark_config(
        dataset,
        case,
        task_dir,
        shared_cache,
        service_url=service_url,
    )
    config.retropath_force = bool(force_retropath)
    a0_sha = _prepare_a0(
        dataset,
        case,
        profile,
        Path(config.chassis_producible_csv),
    )
    fingerprint = _task_fingerprint(
        dataset,
        case,
        profile,
        service_health=service_health,
        with_enzymes=with_enzymes,
    )
    started_at = _utc_now()
    started = time.perf_counter()
    result: dict[str, Any] = {
        "schema_version": BENCHMARK_TASK_SCHEMA,
        "case_id": case.case_id,
        "ec_class": case.ec_class,
        "profile": profile,
        "target_kegg_id": case.target_kegg_id,
        "task_fingerprint": fingerprint,
        "started_at": started_at,
        "completed_at": None,
        "status": "failed",
        "outcome": "not_run",
        "runtime_seconds": 0.0,
        "a0_sha256": a0_sha,
        "search_input": case.search_dict(),
        "gold_standard": {
            "reaction_ids": list(case.gold_reaction_ids),
            "mnxr_ids": list(case.gold_mnxr_ids),
            "ec_numbers": list(case.gold_ec_numbers),
            "uniprot_accessions": list(case.gold_uniprot_accessions),
        },
        "parameters": dataset.defaults.to_dict(),
        "artifacts": {},
        "metrics": {},
        "error": {},
    }
    try:
        pipeline = run_retropath_pipeline(config)
        gap_dir = gap_depth_output_dir(config.gap_output_path, 0)
        retropath_dir = gap_dir / "retropath"
        pipeline_path = retropath_dir / PIPELINE_RESULT_FILE_NAME
        result["artifacts"][PIPELINE_RESULT_FILE_NAME] = _artifact_record(
            pipeline_path, root=run_dir
        )
        validation_dir: Path | None = None
        validation_error: Exception | None = None
        if pipeline.get("status") == "retropath_candidates_found":
            try:
                validation = validate_retropath_candidates(config)
                validation_dir = Path(validation["validation_dir"])
            except (ValueError, RuntimeError) as exc:
                # A candidate-specific reconstruction or strict validation
                # failure is an evaluable algorithm outcome, not an
                # infrastructure outage.  Preserve the P2--P5 metrics and the
                # error instead of dropping the case from the evaluable
                # denominator.
                validation_error = exc
                possible_validation = retropath_dir / "gem_validation"
                if possible_validation.is_dir():
                    validation_dir = possible_validation
            if validation_dir is not None:
                for name in (
                    "validation_manifest.json",
                    "stoichiometry_hypotheses.csv",
                    "stoichiometry_terms.csv",
                    "gem_validation_summary.csv",
                    "gem_validation_route_fluxes.csv",
                    "rejected_hypotheses.csv",
                ):
                    path = validation_dir / name
                    if path.is_file():
                        result["artifacts"][name] = _artifact_record(path, root=run_dir)
                materialization = retropath_dir / "solution_materialization.json"
                if materialization.is_file():
                    result["artifacts"][materialization.name] = _artifact_record(
                        materialization, root=run_dir
                    )
                for name in (
                    "solutions.csv",
                    "all_solution_steps.csv",
                    "solution_electron_summary.csv",
                    "route_electron_requirements.csv",
                ):
                    path = gap_dir / name
                    if path.is_file():
                        result["artifacts"][name] = _artifact_record(path, root=run_dir)
        metrics = score_core_artifacts(
            case,
            pipeline_result_path=pipeline_path,
            validation_dir=validation_dir,
            gap_dir=gap_dir,
        )
        if (
            with_enzymes
            and profile == "full_a0"
            and metrics.get("formal_exact_solution_id") is not None
        ):
            config.solution = int(metrics["formal_exact_solution_id"])
            config.top_n = 10
            config.literature_search = False
            try:
                write_solution(config)
                enzyme_result = run_main_protein_selection(config)
                candidates_path = (
                    task_dir
                    / "main_protein_selection"
                    / "step_main_enzyme_candidates.csv"
                )
                if candidates_path.is_file():
                    result["artifacts"][candidates_path.name] = _artifact_record(
                        candidates_path, root=run_dir
                    )
                    metrics.update(
                        score_enzyme_artifacts(
                            case,
                            candidates_path=candidates_path,
                            source_unavailable=bool(
                                enzyme_result.get("source_unavailable")
                            ),
                        )
                    )
                else:
                    metrics["enzyme_evaluation_status"] = "no_artifact"
            except Exception as exc:  # Enzyme retrieval must not erase core results.
                metrics["enzyme_evaluation_status"] = "failed"
                result["enzyme_error"] = {
                    "type": type(exc).__name__,
                    "detail": str(exc),
                }
        result["metrics"] = metrics
        result["status"] = "completed"
        if validation_error is not None:
            result["outcome"] = "strict_validation_error"
            result["error"] = {
                "stage": "strict_validation",
                "code": type(validation_error).__name__,
                "detail": str(validation_error),
            }
        else:
            result["outcome"] = _outcome(metrics)
    except RetroPathPipelineError as exc:
        stage, code = _classify_pipeline_error(exc)
        if exc.result_path is not None and exc.result_path.is_file():
            result["artifacts"][PIPELINE_RESULT_FILE_NAME] = _artifact_record(
                exc.result_path, root=run_dir
            )
        result["error"] = {
            "stage": stage,
            "code": code,
            "detail": exc.detail,
        }
        result["outcome"] = stage
    except (FileNotFoundError, ValueError) as exc:
        result["error"] = {
            "stage": "validation_or_data",
            "code": type(exc).__name__,
            "detail": str(exc),
        }
        result["outcome"] = "data_or_validation_error"
    except Exception as exc:
        result["error"] = {
            "stage": "unexpected",
            "code": type(exc).__name__,
            "detail": str(exc),
        }
        result["outcome"] = "unexpected_error"
    result["runtime_seconds"] = round(time.perf_counter() - started, 6)
    result["completed_at"] = _utc_now()
    _atomic_write_json(task_dir / TASK_RESULT_FILE_NAME, result)
    return result


def _skip_paired_profile_after_infrastructure_failure(
    dataset: BenchmarkDataset,
    case: BenchmarkCase,
    profile: str,
    *,
    run_dir: Path,
    service_health: Mapping[str, Any],
    with_enzymes: bool,
    controlled_task: Mapping[str, Any],
) -> dict[str, Any]:
    """Record a terminal paired-profile skip without submitting another job.

    RetroPath expands the same target under both P11.1 sink profiles.  When the
    controlled run has already caused an infrastructure-level failure, an
    immediate full-A0 submission can repeat the same resource incident before
    the service has recovered.  The skip remains non-evaluable and is visible
    in the failure funnel; ``--retry-failed`` can retry both profiles later.
    """

    if profile != "full_a0":
        raise ValueError("paired infrastructure circuit breaker only applies to full_a0")
    task_dir = run_dir / "cases" / case.case_id / profile
    task_dir.mkdir(parents=True, exist_ok=True)
    a0_path = task_dir / "chassis_result" / "producible_kegg_compounds.csv"
    a0_sha = _prepare_a0(dataset, case, profile, a0_path)
    fingerprint = _task_fingerprint(
        dataset,
        case,
        profile,
        service_health=service_health,
        with_enzymes=with_enzymes,
    )
    controlled_error = controlled_task.get("error", {})
    controlled_code = (
        str(controlled_error.get("code") or "unknown")
        if isinstance(controlled_error, Mapping)
        else "unknown"
    )
    now = _utc_now()
    result: dict[str, Any] = {
        "schema_version": BENCHMARK_TASK_SCHEMA,
        "case_id": case.case_id,
        "ec_class": case.ec_class,
        "profile": profile,
        "target_kegg_id": case.target_kegg_id,
        "task_fingerprint": fingerprint,
        "started_at": now,
        "completed_at": now,
        "status": "skipped",
        "outcome": "paired_infrastructure_circuit_breaker",
        "runtime_seconds": 0.0,
        "a0_sha256": a0_sha,
        "search_input": case.search_dict(),
        "gold_standard": {
            "reaction_ids": list(case.gold_reaction_ids),
            "mnxr_ids": list(case.gold_mnxr_ids),
            "ec_numbers": list(case.gold_ec_numbers),
            "uniprot_accessions": list(case.gold_uniprot_accessions),
        },
        "parameters": dataset.defaults.to_dict(),
        "artifacts": {},
        "metrics": {},
        "error": {
            "stage": "infrastructure_circuit_breaker",
            "code": "controlled_profile_infrastructure_failure",
            "detail": (
                "full_a0 was not submitted because the paired controlled "
                f"profile ended with infrastructure error {controlled_code}"
            ),
        },
    }
    _atomic_write_json(task_dir / TASK_RESULT_FILE_NAME, result)
    return result


def _task_record(task: Mapping[str, Any], *, run_dir: Path) -> dict[str, Any]:
    path = (
        run_dir
        / "cases"
        / str(task["case_id"])
        / str(task["profile"])
        / TASK_RESULT_FILE_NAME
    )
    return {
        "case_id": task["case_id"],
        "profile": task["profile"],
        "status": task["status"],
        "result_path": str(path.relative_to(run_dir)),
        "sha256": sha256_file(path),
        "task_fingerprint": task["task_fingerprint"],
    }


def _new_run_manifest(
    dataset: BenchmarkDataset,
    *,
    run_id: str,
    profiles: Sequence[str],
    with_enzymes: bool,
    service_url: str,
    service_health: Mapping[str, Any],
) -> dict[str, Any]:
    root = repository_root()
    source_tree = _source_tree_record()
    runner_source = _runner_source_record()
    return {
        "schema_version": BENCHMARK_RUN_SCHEMA,
        "benchmark_id": dataset.benchmark_id,
        "run_id": run_id,
        "created_at": _utc_now(),
        "completed_at": None,
        "state": "running",
        "dataset_path": str(dataset.path),
        "dataset_sha256": dataset.dataset_sha256,
        "profiles": list(profiles),
        "with_enzymes": with_enzymes,
        "service_url": service_url,
        "service_health": dict(service_health),
        "service_identity": _runtime_identity(service_health),
        "defaults": dataset.defaults.to_dict(),
        "resources": {
            name: resource.to_dict()
            for name, resource in sorted(dataset.resources.items())
        },
        "git": _git_provenance(root),
        "source_tree_sha256": source_tree["sha256"],
        "source_files": source_tree["files"],
        "runner_source_sha256": runner_source["sha256"],
        "runner_source_files": runner_source["files"],
        "tasks": [],
        "reports": {},
    }


def run_benchmark(
    cases_path: str | Path = DEFAULT_CASES_PATH,
    *,
    output_root: str | Path = DEFAULT_OUTPUT_ROOT,
    profiles: Sequence[str] = BENCHMARK_PROFILES,
    with_enzymes: bool = False,
    service_url: str = "http://127.0.0.1:8765",
    resume: str | Path | None = None,
    run_id: str | None = None,
    retry_failed: bool = False,
) -> dict[str, Any]:
    dataset = load_benchmark_dataset(cases_path)
    validate_gold_against_mnxref(dataset)
    normalized_profiles = tuple(dict.fromkeys(str(item).strip() for item in profiles))
    invalid_profiles = sorted(set(normalized_profiles) - set(BENCHMARK_PROFILES))
    if not normalized_profiles or invalid_profiles:
        raise ValueError(f"invalid benchmark profiles: {invalid_profiles}")
    health = _service_health(
        service_url,
        request_timeout_seconds=dataset.defaults.request_timeout_seconds,
        get_attempts=dataset.defaults.get_attempts,
    )
    if resume is not None:
        run_dir = Path(resume).expanduser().resolve()
        manifest_path = run_dir / RUN_MANIFEST_FILE_NAME
        manifest = _read_json(manifest_path)
        if manifest.get("schema_version") != BENCHMARK_RUN_SCHEMA:
            raise ValueError("unsupported benchmark run manifest schema")
        if manifest.get("dataset_sha256") != dataset.dataset_sha256:
            if (
                manifest.get("benchmark_id") != dataset.benchmark_id
                or not _resources_compatible(manifest.get("resources"), dataset)
                or _biological_defaults(manifest.get("defaults", {}))
                != _biological_defaults(dataset.defaults.to_dict())
            ):
                raise ValueError("cannot resume: benchmark dataset changed")
            manifest.setdefault("dataset_migrations", []).append(
                {
                    "at": _utc_now(),
                    "from_sha256": manifest.get("dataset_sha256"),
                    "to_sha256": dataset.dataset_sha256,
                    "reason": "transport_observation_policy_only",
                }
            )
            manifest["dataset_sha256"] = dataset.dataset_sha256
            manifest["dataset_path"] = str(dataset.path)
            manifest["defaults"] = dataset.defaults.to_dict()
        if tuple(manifest.get("profiles") or ()) != normalized_profiles:
            raise ValueError("cannot resume with different profiles")
        if bool(manifest.get("with_enzymes")) != bool(with_enzymes):
            raise ValueError("cannot resume with a different enzyme mode")
        previous_identity = manifest.get("service_identity")
        if not isinstance(previous_identity, Mapping):
            previous_identity = _runtime_identity(manifest.get("service_health", {}))
        if dict(previous_identity) != _runtime_identity(health):
            raise ValueError("cannot resume: RetroPath runtime identity changed")
        records = {
            (str(item.get("case_id")), str(item.get("profile"))): item
            for item in manifest.get("tasks", [])
            if isinstance(item, Mapping)
        }
        manifest["state"] = "running"
        manifest["completed_at"] = None
        manifest["reports"] = {}
        manifest["service_health"] = dict(health)
        manifest["service_identity"] = _runtime_identity(health)
        source_tree = _source_tree_record()
        runner_source = _runner_source_record()
        manifest["git"] = _git_provenance(repository_root())
        manifest["source_tree_sha256"] = source_tree["sha256"]
        manifest["source_files"] = source_tree["files"]
        manifest["runner_source_sha256"] = runner_source["sha256"]
        manifest["runner_source_files"] = runner_source["files"]
    else:
        root = Path(output_root).expanduser().resolve()
        identifier = run_id or (
            datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
            + "_"
            + dataset.dataset_sha256[:8]
        )
        if not identifier or any(char in identifier for char in "\\/:*?\"<>|"):
            raise ValueError("run_id contains invalid path characters")
        run_dir = root / identifier
        if run_dir.exists():
            raise FileExistsError(f"benchmark run directory already exists: {run_dir}")
        run_dir.mkdir(parents=True)
        manifest = _new_run_manifest(
            dataset,
            run_id=identifier,
            profiles=normalized_profiles,
            with_enzymes=with_enzymes,
            service_url=service_url,
            service_health=health,
        )
        records = {}
    manifest_path = run_dir / RUN_MANIFEST_FILE_NAME
    manifest.setdefault("execution_attempts", []).append(
        {
            "started_at": _utc_now(),
            "resume": resume is not None,
            "retry_failed": retry_failed,
        }
    )
    _atomic_write_json(manifest_path, manifest)

    completed_records: dict[tuple[str, str], Mapping[str, Any]] = dict(records)
    terminal_tasks: dict[tuple[str, str], Mapping[str, Any]] = {}
    total = len(dataset.cases) * len(normalized_profiles)
    position = 0
    for case in dataset.cases:
        for profile in normalized_profiles:
            position += 1
            fingerprint = _task_fingerprint(
                dataset,
                case,
                profile,
                service_health=health,
                with_enzymes=with_enzymes,
            )
            key = (case.case_id, profile)
            existing = completed_records.get(key)
            if existing is not None and _existing_task_valid(
                existing,
                run_dir=run_dir,
                fingerprint=fingerprint,
                retry_failed=retry_failed,
            ):
                print(f"[{position}/{total}] resume {case.case_id}/{profile}")
                terminal_tasks[key] = _task_from_record(existing, run_dir=run_dir)
                continue
            controlled_task = terminal_tasks.get((case.case_id, "controlled"))
            if profile == "full_a0" and _is_infrastructure_failure(controlled_task):
                print(
                    f"[{position}/{total}] circuit-break "
                    f"{case.case_id}/{profile}"
                )
                task = _skip_paired_profile_after_infrastructure_failure(
                    dataset,
                    case,
                    profile,
                    run_dir=run_dir,
                    service_health=health,
                    with_enzymes=with_enzymes,
                    controlled_task=controlled_task,
                )
            else:
                print(f"[{position}/{total}] run {case.case_id}/{profile}")
                task = _run_task(
                    dataset,
                    case,
                    profile,
                    run_dir=run_dir,
                    service_url=service_url,
                    service_health=health,
                    with_enzymes=with_enzymes,
                    force_retropath=retry_failed,
                )
            terminal_tasks[key] = task
            completed_records[key] = _task_record(task, run_dir=run_dir)
            manifest["tasks"] = [
                completed_records[item]
                for item in sorted(completed_records)
            ]
            _atomic_write_json(manifest_path, manifest)

    task_payloads = [
        _read_json(run_dir / str(record["result_path"]))
        for record in completed_records.values()
    ]
    failed = sum(task.get("status") != "completed" for task in task_payloads)
    manifest["state"] = "core_complete" if failed == 0 else "core_complete_with_failures"
    manifest["completed_at"] = _utc_now()
    manifest["tasks"] = [completed_records[item] for item in sorted(completed_records)]
    _atomic_write_json(manifest_path, manifest)
    report = generate_benchmark_report(run_dir)
    manifest["reports"] = report["artifacts"]
    manifest["state"] = "completed" if failed == 0 else "completed_with_failures"
    _atomic_write_json(manifest_path, manifest)
    return {
        "ok": failed == 0,
        "run_id": manifest["run_id"],
        "run_dir": str(run_dir),
        "task_count": len(task_payloads),
        "failed_task_count": failed,
        "state": manifest["state"],
        "report": report,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="RetroPath P11.1 benchmark")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="validate dataset and gold identities")
    validate.add_argument("--cases", type=Path, default=DEFAULT_CASES_PATH)

    run = commands.add_parser("run", help="run controlled and full-A0 benchmark tasks")
    run.add_argument("--cases", type=Path, default=DEFAULT_CASES_PATH)
    run.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_ROOT)
    run.add_argument(
        "--profiles",
        nargs="+",
        choices=BENCHMARK_PROFILES,
        default=list(BENCHMARK_PROFILES),
    )
    run.add_argument("--with-enzymes", action="store_true")
    run.add_argument("--service-url", default="http://127.0.0.1:8765")
    run.add_argument("--resume", type=Path)
    run.add_argument("--run-id")
    run.add_argument(
        "--retry-failed",
        action="store_true",
        help="retry terminal failed/skipped tasks instead of resuming them",
    )

    report = commands.add_parser("report", help="regenerate deterministic reports")
    report.add_argument("--run", type=Path, required=True)
    report.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "validate":
        dataset = load_benchmark_dataset(args.cases)
        result = validate_gold_against_mnxref(dataset)
    elif args.command == "run":
        result = run_benchmark(
            args.cases,
            output_root=args.output,
            profiles=args.profiles,
            with_enzymes=args.with_enzymes,
            service_url=args.service_url,
            resume=args.resume,
            run_id=args.run_id,
            retry_failed=args.retry_failed,
        )
    else:
        result = generate_benchmark_report(args.run, output_path=args.output)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 0 if result.get("ok", False) else 1


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "BENCHMARK_PROFILES",
    "DEFAULT_CASES_PATH",
    "DEFAULT_OUTPUT_ROOT",
    "RUN_MANIFEST_FILE_NAME",
    "TASK_RESULT_FILE_NAME",
    "main",
    "run_benchmark",
    "validate_gold_against_mnxref",
]
