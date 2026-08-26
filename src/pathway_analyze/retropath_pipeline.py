"""User-facing orchestration for the isolated RetroPath candidate workflow.

The regular ``gap`` command remains the KEGG implementation.  This module is
only selected when the user explicitly passes ``--retropath`` and then runs
the already isolated P2 -> P3 -> P4 -> P5 stages without writing formal KEGG
solutions or design manifests.
"""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.pathway_analyze.expand_chassis_metabolites import load_expansion_bundle
from src.pathway_analyze.kegg_gap_analyze import (
    KeggRestClient,
    gap_depth_output_dir,
    run_gap,
)
from src.pathway_analyze.retropath_analyze import analyze_retropath_candidates
from src.pathway_analyze.retropath_client import (
    RetroPathClientError,
    RetroPathHttpClient,
    RetroPathJobParameters,
)
from src.pathway_analyze.retropath_input import (
    RetroPathInputBuildError,
    build_retropath_inputs,
)
from src.pathway_analyze.retropath_parser import RetroPathParseError
from src.pathway_analyze.retropath_routes import parse_and_enumerate_retropath
from src.pathway_analyze.retropath_structure import KeggMolStructureProvider
from src.pathway_analyze.target_id import validate_target_compound_id


RETROPATH_PIPELINE_SCHEMA = "retropath_pipeline_result.v2"
PIPELINE_RESULT_FILE_NAME = "pipeline_result.json"


class RetroPathPipelineError(RuntimeError):
    """Stable user-facing failure from one RetroPath pipeline stage."""

    def __init__(
        self,
        status: str,
        detail: str,
        *,
        stage: str,
        result_path: Path | None = None,
    ) -> None:
        self.status = str(status).strip()
        self.detail = str(detail).strip()
        self.stage = str(stage).strip()
        self.result_path = result_path
        suffix = f"; result: {result_path}" if result_path is not None else ""
        super().__init__(f"{self.status}: {self.detail}{suffix}")


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
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
            json.dump(
                payload,
                handle,
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
                allow_nan=False,
            )
            handle.write("\n")
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _pipeline_output_dir(config: Any, depth: int) -> Path:
    return (
        gap_depth_output_dir(Path(config.gap_output_path), depth) / "retropath"
    ).expanduser().resolve()


def _write_failure(
    *,
    output_dir: Path,
    target_compound: str,
    depth: int,
    status: str,
    detail: str,
    stage: str,
) -> RetroPathPipelineError:
    result_path = output_dir / PIPELINE_RESULT_FILE_NAME
    payload = {
        "schema_version": RETROPATH_PIPELINE_SCHEMA,
        "ok": False,
        "retropath_requested": True,
        "search_engine": "retropath",
        "status": status,
        "stage": stage,
        "detail": detail,
        "target_compound": target_compound,
        "expansion_depth": depth,
        "sink_source": (
            "chassis_A0" if depth == 0 else f"cumulative_expansion_A{depth}"
        ),
        "output_dir": str(output_dir),
        "pipeline_result_file": str(result_path),
    }
    _atomic_write_json(result_path, payload)
    return RetroPathPipelineError(
        status,
        detail,
        stage=stage,
        result_path=result_path,
    )


def _client_error_status(error: RetroPathClientError) -> str:
    if error.code == "client_poll_timeout":
        return "retropath_timeout"
    if error.code in {
        "runtime_not_ready",
        "service_unavailable",
        "queue_full",
        "submission_uncertain",
    }:
        return "retropath_service_unavailable"
    if error.code in {"input_invalid", "input_modified"}:
        return "retropath_input_invalid"
    return "retropath_execution_failed"


def _job_parameters(config: Any) -> RetroPathJobParameters:
    return RetroPathJobParameters(
        max_steps=getattr(config, "retropath_max_steps", 3),
        topx=getattr(config, "retropath_topx", 100),
        dmin=getattr(config, "retropath_dmin", 2),
        dmax=getattr(config, "retropath_dmax", 16),
        mwmax_source=getattr(config, "retropath_mwmax_source", 1000),
        msc_timeout=getattr(config, "retropath_msc_timeout", 10),
    )


def run_retropath_pipeline(config: Any) -> dict[str, Any]:
    """Run P2 -> P3 -> P4 -> P5 and return one auditable result envelope."""

    try:
        target_compound = validate_target_compound_id(config.target_name)
        depth = int(getattr(config, "depth", 0))
        if depth < 0:
            raise ValueError("depth must be greater than or equal to 0")
        output_dir = _pipeline_output_dir(config, depth)
        rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    except (AttributeError, TypeError, ValueError) as exc:
        raise RetroPathPipelineError(
            "retropath_configuration_invalid",
            str(exc),
            stage="configuration",
        ) from exc

    if not rules_path.is_file():
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_rules_missing",
            detail=f"RetroRules RR02 retro file does not exist: {rules_path}",
            stage="configuration",
        )

    try:
        expansion_bundle = load_expansion_bundle(
            base_path=Path(config.chassis_producible_csv),
            output_dir=Path(config.chassis_output_path),
            depth=depth,
        )
    except (OSError, ValueError) as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_expansion_missing",
            detail=str(exc),
            stage="expansion",
        ) from exc

    try:
        structure_provider = KeggMolStructureProvider(
            Path(config.cache_dir) / "kegg",
            timeout_seconds=getattr(config, "retropath_structure_timeout_seconds", 30.0),
            retries=getattr(config, "retropath_structure_retries", 3),
            request_sleep_seconds=getattr(
                config,
                "retropath_structure_request_sleep_seconds",
                0.2,
            ),
        )
        input_bundle = build_retropath_inputs(
            target_compound,
            expansion_bundle,
            structure_provider,
            output_dir / "input",
        )
    except RetroPathInputBuildError as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_input_invalid",
            detail=f"{exc.code}: {exc.detail}",
            stage="input",
        ) from exc
    except (OSError, ValueError) as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_input_invalid",
            detail=str(exc),
            stage="input",
        ) from exc

    if bool(getattr(input_bundle, "target_already_reachable", False)):
        result_path = output_dir / PIPELINE_RESULT_FILE_NAME
        payload = {
            "schema_version": RETROPATH_PIPELINE_SCHEMA,
            "ok": True,
            "retropath_requested": True,
            "search_engine": "retropath",
            "status": "retropath_target_already_reachable",
            "stage": "input",
            "detail": (
                "the exact target structure is already present in the "
                "cumulative chassis sink; RetroPath was not submitted"
            ),
            "target_compound": target_compound,
            "target_reachable_aliases": list(
                getattr(input_bundle, "target_reachable_aliases", tuple())
            ),
            "expansion_depth": depth,
            "sink_source": (
                "chassis_A0" if depth == 0 else f"cumulative_expansion_A{depth}"
            ),
            "output_dir": str(output_dir),
            "pipeline_result_file": str(result_path),
            "job_id": None,
            "service_status": "not_submitted",
            "return_code": None,
            "cache_hit": False,
            "complete_path_count": 0,
            "sink_match_count": 1,
            "scope_present": False,
            "candidate_count": 0,
            "rejection_count": input_bundle.rejected_compound_count,
            "input_summary": {
                "reachable_compound_count": input_bundle.reachable_compound_count,
                "sink_structure_count": input_bundle.sink_structure_count,
                "rejected_compound_count": input_bundle.rejected_compound_count,
            },
            "artifacts": {
                "expansion_source": {
                    "path": str(expansion_bundle.expanded_file.resolve()),
                    "sha256": _sha256_file(expansion_bundle.expanded_file),
                },
                "retropath_rules": {
                    "path": str(rules_path),
                    "sha256": _sha256_file(rules_path),
                },
                "target_source": {
                    "path": str(input_bundle.target_source_path.resolve()),
                    "sha256": input_bundle.target_source_sha256,
                },
                "chassis_sink": {
                    "path": str(input_bundle.chassis_sink_path.resolve()),
                    "sha256": input_bundle.chassis_sink_sha256,
                },
                "compound_mapping": {
                    "path": str(input_bundle.compound_mapping_path.resolve()),
                    "sha256": _sha256_file(input_bundle.compound_mapping_path),
                },
                "rejected_compounds": {
                    "path": str(input_bundle.rejected_compounds_path.resolve()),
                    "sha256": _sha256_file(input_bundle.rejected_compounds_path),
                },
            },
        }
        _atomic_write_json(result_path, payload)
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return payload

    try:
        job_parameters = _job_parameters(config)
        with RetroPathHttpClient(
            getattr(config, "retropath_service_url", "http://127.0.0.1:8765"),
            request_timeout_seconds=getattr(
                config,
                "retropath_request_timeout_seconds",
                30.0,
            ),
            get_attempts=getattr(config, "retropath_get_attempts", 3),
            retry_backoff_seconds=getattr(
                config,
                "retropath_retry_backoff_seconds",
                0.5,
            ),
            poll_interval_seconds=getattr(
                config,
                "retropath_poll_interval_seconds",
                1.0,
            ),
            wait_timeout_seconds=getattr(
                config,
                "retropath_wait_timeout_seconds",
                3900.0,
            ),
        ) as client:
            client_run = client.run(
                input_bundle,
                output_dir,
                parameters=job_parameters,
                force=getattr(config, "retropath_force", False),
            )
    except RetroPathClientError as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status=_client_error_status(exc),
            detail=f"{exc.code}: {exc.detail}",
            stage="client",
        ) from exc
    except (OSError, ValueError) as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_configuration_invalid",
            detail=str(exc),
            stage="client",
        ) from exc

    service_status = client_run.result.status
    if service_status in {"failed", "timed_out"}:
        status = (
            "retropath_timeout"
            if service_status == "timed_out"
            else "retropath_execution_failed"
        )
        detail = "; ".join(client_run.result.errors) or (
            f"RetroPath service ended with status {service_status}"
        )
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status=status,
            detail=detail,
            stage="client",
        )

    try:
        enumeration_result = parse_and_enumerate_retropath(
            client_run,
            input_bundle,
            rules_path,
            max_routes=getattr(config, "retropath_max_routes", 1000),
            max_search_states=getattr(
                config,
                "retropath_max_search_states",
                100000,
            ),
        )
    except RetroPathParseError as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_parse_failed",
            detail=f"{exc.code}: {exc.detail}",
            stage="enumeration",
        ) from exc
    except (OSError, ValueError) as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_parse_failed",
            detail=str(exc),
            stage="enumeration",
        ) from exc

    try:
        kegg_client = KeggRestClient(Path(config.cache_dir) / "kegg")
        candidate_artifacts = analyze_retropath_candidates(
            enumeration_result,
            expansion_bundle,
            kegg_client,
            output_dir,
            max_candidates=getattr(config, "retropath_max_candidates", 5),
            max_witness_plans=getattr(
                config,
                "retropath_max_witness_plans",
                3,
            ),
            max_total_steps=getattr(config, "max_total_steps", 20),
            max_new_enzymes=getattr(config, "max_new_enzymes", 20),
        )
    except (OSError, ValueError) as exc:
        raise _write_failure(
            output_dir=output_dir,
            target_compound=target_compound,
            depth=depth,
            status="retropath_merge_failed",
            detail=str(exc),
            stage="candidate_merge",
        ) from exc

    if candidate_artifacts.candidate_count:
        status = "retropath_candidates_found"
    elif service_status == "source_in_sink":
        status = "retropath_source_in_sink"
    else:
        status = "retropath_no_scope"

    result_path = output_dir / PIPELINE_RESULT_FILE_NAME
    provenance = client_run.result.provenance
    scope_artifact_names = {
        "target_scope.json",
        "scope.json",
        "target_scope.csv",
        "scope.csv",
        "results.csv",
    }
    scope_artifacts = (
        sorted(
            str(path.resolve())
            for path in client_run.raw_dir.iterdir()
            if path.is_file() and path.name in scope_artifact_names
        )
        if client_run.raw_dir.is_dir()
        else []
    )
    payload: dict[str, Any] = {
        "schema_version": RETROPATH_PIPELINE_SCHEMA,
        "ok": True,
        "retropath_requested": True,
        "search_engine": "retropath",
        "status": status,
        "target_compound": target_compound,
        "expansion_depth": depth,
        "sink_source": (
            "chassis_A0" if depth == 0 else f"cumulative_expansion_A{depth}"
        ),
        "output_dir": str(output_dir),
        "pipeline_result_file": str(result_path),
        "job_id": client_run.result.job_id,
        "service_status": service_status,
        "return_code": client_run.result.return_code,
        "failure_code": getattr(client_run.result, "failure_code", None),
        "cache_hit": client_run.cache_hit,
        "complete_path_count": enumeration_result.complete_path_count,
        "sink_match_count": len(enumeration_result.network.sink_matches),
        "scope_present": bool(scope_artifacts),
        "candidate_count": candidate_artifacts.candidate_count,
        "rejection_count": candidate_artifacts.rejection_count,
        "upstream_enumeration_truncated": enumeration_result.truncated,
        "candidate_top_k_truncated": candidate_artifacts.merge_result.truncated,
        "input_summary": {
            "reachable_compound_count": input_bundle.reachable_compound_count,
            "sink_structure_count": input_bundle.sink_structure_count,
            "rejected_compound_count": input_bundle.rejected_compound_count,
        },
        "parameters": {
            "retropath_job": job_parameters.to_dict(),
            "enumeration": {
                "max_routes": enumeration_result.max_routes,
                "max_search_states": enumeration_result.max_search_states,
            },
            "candidate_merge": {
                "max_candidates": candidate_artifacts.merge_result.max_candidates,
                "max_witness_plans": candidate_artifacts.merge_result.max_witness_plans,
                "max_total_steps": candidate_artifacts.merge_result.max_total_steps,
                "max_new_enzymes": candidate_artifacts.merge_result.max_new_enzymes,
            },
        },
        "provenance": None if provenance is None else provenance.to_dict(),
        "artifacts": {
            "expansion_source": {
                "path": str(expansion_bundle.expanded_file.resolve()),
                "sha256": _sha256_file(expansion_bundle.expanded_file),
            },
            "retropath_rules": {
                "path": str(rules_path),
                "sha256": None if provenance is None else provenance.rules_sha256,
            },
            "target_source": {
                "path": str(input_bundle.target_source_path.resolve()),
                "sha256": input_bundle.target_source_sha256,
            },
            "chassis_sink": {
                "path": str(input_bundle.chassis_sink_path.resolve()),
                "sha256": input_bundle.chassis_sink_sha256,
            },
            "compound_mapping": {
                "path": str(input_bundle.compound_mapping_path.resolve()),
                "sha256": _sha256_file(input_bundle.compound_mapping_path),
            },
            "rejected_compounds": {
                "path": str(input_bundle.rejected_compounds_path.resolve()),
                "sha256": _sha256_file(input_bundle.rejected_compounds_path),
            },
            "client_run_manifest": str(client_run.run_manifest_path.resolve()),
            "client_state": str(client_run.client_state_path.resolve()),
            "raw_directory": str(client_run.raw_dir.resolve()),
            "scope_artifacts": scope_artifacts,
            "candidate_routes": {
                "path": str(candidate_artifacts.candidate_routes_path.resolve()),
                "sha256": candidate_artifacts.candidate_routes_sha256,
            },
            "candidate_steps": {
                "path": str(candidate_artifacts.candidate_steps_path.resolve()),
                "sha256": candidate_artifacts.candidate_steps_sha256,
            },
            "rejected_routes": {
                "path": str(candidate_artifacts.rejected_routes_path.resolve()),
                "sha256": candidate_artifacts.rejected_routes_sha256,
            },
        },
    }
    _atomic_write_json(result_path, payload)
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return payload


def run_gap_command(config: Any) -> dict[str, Any]:
    """Dispatch ``gap`` while preserving the original KEGG default exactly."""

    if not getattr(config, "retropath", False):
        return run_gap(config)
    return run_retropath_pipeline(config)


__all__ = [
    "PIPELINE_RESULT_FILE_NAME",
    "RETROPATH_PIPELINE_SCHEMA",
    "RetroPathPipelineError",
    "run_gap_command",
    "run_retropath_pipeline",
]
