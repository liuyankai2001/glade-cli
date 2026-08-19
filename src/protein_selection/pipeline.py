"""Run auxiliary-protein research for a selected main-enzyme combination."""

from __future__ import annotations

from collections.abc import Awaitable, Callable, Mapping
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time
from typing import Any

import httpx
from langchain_core.language_models import BaseChatModel

from src.main_protein_selection.provenance import stable_json_hash
from src.protein_selection.agents.main_research_agent import (
    MainResearchResult,
    open_main_research_agent,
)
from src.protein_selection.graph import build_protein_supply_graph
from src.protein_selection.manifest_adapter import (
    build_main_enzyme_research_units,
)
from src.protein_selection.models import (
    AUXILIARY_PROTEIN_PIPELINE_VERSION,
    AuxiliaryProteinCombinationResult,
    MainEnzymeAuxiliaryResearchResult,
    MainEnzymeResearchStatus,
    auxiliary_protein_result_fingerprint,
)
from src.protein_selection.reaction_scope import normalize_reaction_ids
from src.protein_selection.services.cache import (
    PersistentTTLCache,
    RetrievalStats,
)
from src.protein_selection.services.kegg_client import KeggClient
from src.protein_selection.services.uniprot_client import UniProtClient
from src.protein_selection.state import MainEnzymeResearchUnit, ResearchMode
from src.protein_selection.tools.iml1515 import DEFAULT_IML1515_PATH
from src.write_manifest.store import read_design_manifest


AUXILIARY_PROTEIN_RESULT_FILENAME = "auxiliary_protein_research.json"
UnitResearchRunner = Callable[
    [MainEnzymeResearchUnit],
    Awaitable[Mapping[str, Any]],
]


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"design manifest omitted {field_name}")
    return value


def _string(value: Any, field_name: str) -> str:
    normalized = str(value or "").strip()
    if not normalized:
        raise ValueError(f"design manifest omitted {field_name}")
    return normalized


def _integer(value: Any, field_name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ValueError(f"design manifest {field_name} must be an integer")
    try:
        normalized = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"design manifest {field_name} must be an integer"
        ) from exc
    if normalized < minimum:
        raise ValueError(
            f"design manifest {field_name} must be >= {minimum}"
        )
    return normalized


def _unit_reaction_ids(unit: MainEnzymeResearchUnit) -> list[str]:
    values = list(
        dict.fromkeys(
            str(step["reaction_id"]).strip().upper()
            for step in unit["reaction_steps"]
        )
    )
    return normalize_reaction_ids(values)


def _unit_status(result: MainResearchResult) -> MainEnzymeResearchStatus:
    if result.outcome in {
        "independent",
        "host_supported",
        "supplement_required",
    }:
        return "complete"
    if result.outcome == "reaction_mismatch":
        return "blocked"
    return "review_required"


def _service_error_unit(
    unit: MainEnzymeResearchUnit,
    message: str,
) -> MainEnzymeAuxiliaryResearchResult:
    return MainEnzymeAuxiliaryResearchResult(
        accession=unit["accession"],
        sequence_sha256=unit["sequence_sha256"],
        reaction_scope=unit["reaction_scope"],
        assigned_step_indexes=unit["assigned_step_indexes"],
        reaction_ids=_unit_reaction_ids(unit),
        status="service_error",
        error=message.strip() or "unknown research service error",
    )


def _unit_result_from_state(
    unit: MainEnzymeResearchUnit,
    state: Mapping[str, Any],
) -> MainEnzymeAuxiliaryResearchResult:
    completion_status = state.get("completion_status")
    if completion_status == "service_error":
        return _service_error_unit(
            unit,
            str(state.get("research_error") or "research graph failed"),
        )
    raw_result = state.get("main_research_result")
    if raw_result is None:
        return _service_error_unit(
            unit,
            "research graph did not return main_research_result",
        )
    try:
        if isinstance(raw_result, MainResearchResult):
            result = raw_result
        elif isinstance(raw_result, Mapping):
            normalized_result = dict(raw_result)
            for field_name in MainResearchResult.model_computed_fields:
                normalized_result.pop(field_name, None)
            result = MainResearchResult.model_validate(normalized_result)
        else:
            raise TypeError("main_research_result must be an object")
        return MainEnzymeAuxiliaryResearchResult(
            accession=unit["accession"],
            sequence_sha256=unit["sequence_sha256"],
            reaction_scope=unit["reaction_scope"],
            assigned_step_indexes=unit["assigned_step_indexes"],
            reaction_ids=_unit_reaction_ids(unit),
            status=_unit_status(result),
            research_result=result,
        )
    except Exception as exc:
        return _service_error_unit(
            unit,
            f"invalid research result: {type(exc).__name__}: {exc}",
        )


def _input_fingerprint(
    manifest: Mapping[str, Any],
    units: list[MainEnzymeResearchUnit],
    research_mode: ResearchMode,
) -> str:
    solution = _mapping(manifest.get("solution"), "solution")
    selection = _mapping(
        manifest.get("main_enzyme_selection"),
        "main_enzyme_selection",
    )
    source = _mapping(selection.get("source"), "main_enzyme_selection.source")
    return stable_json_hash(
        {
            "pipeline_version": AUXILIARY_PROTEIN_PIPELINE_VERSION,
            "target_compound_id": _string(
                manifest.get("target_compound_id"),
                "target_compound_id",
            ),
            "selected_solution_id": _integer(
                solution.get("solution_id"),
                "solution.solution_id",
                minimum=1,
            ),
            "expansion_depth": _integer(
                solution.get("expansion_depth", 0),
                "solution.expansion_depth",
            ),
            "selected_set_id": _integer(
                selection.get("selected_set_id"),
                "main_enzyme_selection.selected_set_id",
                minimum=1,
            ),
            "selected_set_fingerprint": _string(
                selection.get("selected_set_fingerprint"),
                "main_enzyme_selection.selected_set_fingerprint",
            ).lower(),
            "solution_fingerprint": _string(
                source.get("solution_fingerprint"),
                "main_enzyme_selection.source.solution_fingerprint",
            ).lower(),
            "chassis_key": _string(
                selection.get("chassis_key"),
                "main_enzyme_selection.chassis_key",
            ),
            "research_mode": research_mode,
            "units": [
                {
                    "accession": unit["accession"],
                    "sequence_sha256": unit["sequence_sha256"],
                    "reaction_scope": unit["reaction_scope"],
                    "assigned_step_indexes": unit["assigned_step_indexes"],
                    "reaction_steps": unit["reaction_steps"],
                    "whole_reaction": unit.get("whole_reaction"),
                }
                for unit in units
            ],
        }
    )


def _stable_unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value.upper() for value in values))


def _aggregate_accessions(
    results: list[MainEnzymeAuxiliaryResearchResult],
) -> dict[str, list[str]]:
    required: list[str] = []
    recommended: list[str] = []
    host_available: list[str] = []
    to_introduce: list[str] = []
    candidates: list[str] = []
    for item in results:
        result = item.research_result
        if result is None:
            continue
        required.extend(
            protein.uniprot_id
            for protein in result.auxiliary_proteins
            if protein.necessity == "required"
        )
        recommended.extend(
            protein.uniprot_id
            for protein in result.auxiliary_proteins
            if protein.necessity == "enhancing"
        )
        host_available.extend(result.host_available_auxiliary_proteins)
        to_introduce.extend(result.auxiliary_proteins_to_introduce)
        candidates.extend(result.candidate_proteins)
    return {
        "required_auxiliary_protein_accessions": _stable_unique(required),
        "recommended_auxiliary_protein_accessions": _stable_unique(
            recommended
        ),
        "host_available_auxiliary_protein_accessions": _stable_unique(
            host_available
        ),
        "auxiliary_proteins_to_introduce": _stable_unique(to_introduce),
        "candidate_auxiliary_protein_accessions": _stable_unique(candidates),
    }


def _combination_status(
    results: list[MainEnzymeAuxiliaryResearchResult],
) -> MainEnzymeResearchStatus:
    if any(item.status == "service_error" for item in results):
        return "service_error"
    if any(item.status == "blocked" for item in results):
        return "blocked"
    if any(item.status == "review_required" for item in results):
        return "review_required"
    return "complete"


def _blocking_reasons(
    results: list[MainEnzymeAuxiliaryResearchResult],
) -> list[str]:
    reasons: list[str] = []
    for item in results:
        if item.status == "service_error":
            reasons.append(f"{item.accession}: {item.error}")
        elif item.status == "blocked":
            reasons.append(
                f"{item.accession}: selected main enzyme does not match its "
                "assigned reaction scope"
            )
        elif item.status == "review_required":
            reasons.append(
                f"{item.accession}: auxiliary-protein evidence remains unresolved"
            )
    return reasons


def _combination_result(
    manifest: Mapping[str, Any],
    units: list[MainEnzymeResearchUnit],
    unit_results: list[MainEnzymeAuxiliaryResearchResult],
    *,
    research_mode: ResearchMode,
    source_manifest: str,
) -> AuxiliaryProteinCombinationResult:
    solution = _mapping(manifest.get("solution"), "solution")
    selection = _mapping(
        manifest.get("main_enzyme_selection"),
        "main_enzyme_selection",
    )
    source = _mapping(selection.get("source"), "main_enzyme_selection.source")
    input_fingerprint = _input_fingerprint(manifest, units, research_mode)
    aggregate = _aggregate_accessions(unit_results)
    main_accessions = [item.accession for item in unit_results]
    status = _combination_status(unit_results)
    blocking_reasons = _blocking_reasons(unit_results)
    warnings = (
        []
        if status == "complete"
        else [
            "Successful unit results were retained, but this combination "
            "cannot advance until every blocking unit is resolved."
        ]
    )
    complete_protein_list = _stable_unique(
        [
            *main_accessions,
            *aggregate["required_auxiliary_protein_accessions"],
        ]
    )
    result_fingerprint = auxiliary_protein_result_fingerprint(
        input_fingerprint=input_fingerprint,
        status=status,
        main_enzyme_results=unit_results,
        main_enzyme_accessions=main_accessions,
        required_auxiliary_protein_accessions=(
            aggregate["required_auxiliary_protein_accessions"]
        ),
        recommended_auxiliary_protein_accessions=(
            aggregate["recommended_auxiliary_protein_accessions"]
        ),
        host_available_auxiliary_protein_accessions=(
            aggregate["host_available_auxiliary_protein_accessions"]
        ),
        auxiliary_proteins_to_introduce=(
            aggregate["auxiliary_proteins_to_introduce"]
        ),
        candidate_auxiliary_protein_accessions=(
            aggregate["candidate_auxiliary_protein_accessions"]
        ),
        complete_protein_list=complete_protein_list,
        blocking_reasons=blocking_reasons,
        warnings=warnings,
    )
    return AuxiliaryProteinCombinationResult(
        generated_at=datetime.now(timezone.utc).isoformat(),
        source_manifest=source_manifest,
        source_manifest_revision=_integer(
            manifest.get("revision", 0),
            "revision",
        ),
        target_compound_id=_string(
            manifest.get("target_compound_id"),
            "target_compound_id",
        ).upper(),
        selected_solution_id=_integer(
            solution.get("solution_id"),
            "solution.solution_id",
            minimum=1,
        ),
        expansion_depth=_integer(
            solution.get("expansion_depth", 0),
            "solution.expansion_depth",
        ),
        selected_set_id=_integer(
            selection.get("selected_set_id"),
            "main_enzyme_selection.selected_set_id",
            minimum=1,
        ),
        selected_set_fingerprint=_string(
            selection.get("selected_set_fingerprint"),
            "main_enzyme_selection.selected_set_fingerprint",
        ).lower(),
        solution_fingerprint=_string(
            source.get("solution_fingerprint"),
            "main_enzyme_selection.source.solution_fingerprint",
        ).lower(),
        chassis_key=_string(
            selection.get("chassis_key"),
            "main_enzyme_selection.chassis_key",
        ),
        research_mode=research_mode,
        input_fingerprint=input_fingerprint,
        result_fingerprint=result_fingerprint,
        status=status,
        can_advance=status == "complete",
        main_enzyme_results=unit_results,
        main_enzyme_accessions=main_accessions,
        complete_protein_list=complete_protein_list,
        blocking_reasons=blocking_reasons,
        warnings=warnings,
        **aggregate,
    )


async def _execute_validated_units(
    manifest: Mapping[str, Any],
    units: list[MainEnzymeResearchUnit],
    unit_runner: UnitResearchRunner,
    *,
    research_mode: ResearchMode,
    source_manifest: str,
) -> AuxiliaryProteinCombinationResult:
    unit_results: list[MainEnzymeAuxiliaryResearchResult] = []
    for unit in units:
        try:
            state = await unit_runner(unit)
            if not isinstance(state, Mapping):
                raise TypeError("unit runner returned a non-mapping state")
            unit_results.append(_unit_result_from_state(unit, state))
        except Exception as exc:
            unit_results.append(
                _service_error_unit(
                    unit,
                    f"{type(exc).__name__}: {exc}",
                )
            )
    return _combination_result(
        manifest,
        units,
        unit_results,
        research_mode=research_mode,
        source_manifest=source_manifest,
    )


async def execute_research_units(
    manifest: Mapping[str, Any],
    unit_runner: UnitResearchRunner,
    *,
    research_mode: ResearchMode = "balanced",
    source_manifest: str | Path | None = None,
) -> AuxiliaryProteinCombinationResult:
    """Run an injected unit researcher sequentially for one selected set."""

    units = build_main_enzyme_research_units(manifest)
    source = (
        str(Path(source_manifest).expanduser().resolve(strict=False))
        if source_manifest is not None
        else "<in-memory>"
    )
    return await _execute_validated_units(
        manifest,
        units,
        unit_runner,
        research_mode=research_mode,
        source_manifest=source,
    )


def auxiliary_protein_result_path(
    manifest_path: str | Path,
    output_dir: str | Path | None = None,
) -> Path:
    manifest = Path(manifest_path).expanduser().resolve(strict=False)
    directory = (
        Path(output_dir).expanduser().resolve(strict=False)
        if output_dir is not None
        else manifest.parent / "protein_selection"
    )
    return directory / AUXILIARY_PROTEIN_RESULT_FILENAME


def write_auxiliary_protein_result(
    result: AuxiliaryProteinCombinationResult,
    path: str | Path,
) -> Path:
    """Atomically write one canonical combination result."""

    target = Path(path).expanduser().resolve(strict=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(
        f".{target.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    try:
        temporary.write_text(
            json.dumps(
                result.model_dump(mode="json"),
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


async def run_auxiliary_protein_pipeline(
    manifest_path: str | Path,
    *,
    model: BaseChatModel,
    cache_dir: str | Path,
    output_dir: str | Path | None = None,
    research_mode: ResearchMode = "balanced",
    subagent_model: BaseChatModel | None = None,
    iml1515_path: str | Path = DEFAULT_IML1515_PATH,
    refresh_cache: bool = False,
    http_client: httpx.Client | None = None,
    http_timeout_seconds: float = 30.0,
) -> AuxiliaryProteinCombinationResult:
    """Run the live graph once per main enzyme and write a canonical result."""

    path = Path(manifest_path).expanduser().resolve(strict=False)
    manifest = read_design_manifest(path)
    units = build_main_enzyme_research_units(manifest)
    cache = PersistentTTLCache(
        Path(cache_dir).expanduser().resolve(strict=False),
        refresh=refresh_cache,
    )

    async def run_with_client(
        client: httpx.Client,
    ) -> AuxiliaryProteinCombinationResult:
        async def unit_runner(
            unit: MainEnzymeResearchUnit,
        ) -> Mapping[str, Any]:
            stats = RetrievalStats()
            uniprot_client = UniProtClient(
                client,
                cache=cache,
                stats=stats,
            )
            kegg_client = KeggClient(
                client,
                cache=cache,
                stats=stats,
            )

            def research_agent_context_factory():
                return open_main_research_agent(
                    model=model,
                    subagent_model=subagent_model,
                    iml1515_path=iml1515_path,
                    research_mode=research_mode,
                    response_cache=cache,
                    retrieval_stats=stats,
                )

            graph = build_protein_supply_graph(
                llm=model,
                uniprot_client=uniprot_client,
                kegg_client=kegg_client,
                research_agent_context_factory=(
                    research_agent_context_factory
                ),
            )
            return await graph.ainvoke(
                {
                    "research_unit": unit,
                    "research_mode": research_mode,
                }
            )

        return await _execute_validated_units(
            manifest,
            units,
            unit_runner,
            research_mode=research_mode,
            source_manifest=str(path),
        )

    if http_client is None:
        with httpx.Client(timeout=http_timeout_seconds) as owned_client:
            result = await run_with_client(owned_client)
    else:
        result = await run_with_client(http_client)

    write_auxiliary_protein_result(
        result,
        auxiliary_protein_result_path(path, output_dir),
    )
    return result


__all__ = [
    "AUXILIARY_PROTEIN_RESULT_FILENAME",
    "UnitResearchRunner",
    "auxiliary_protein_result_path",
    "execute_research_units",
    "run_auxiliary_protein_pipeline",
    "write_auxiliary_protein_result",
]
