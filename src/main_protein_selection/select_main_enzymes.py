from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import requests

from src.main_protein_selection.build_main_enzyme_sets import (
    shortlist_decision_fingerprint,
)
from src.main_protein_selection.common import (
    PROTEIN_CANDIDATE_COLUMNS,
    STEP_CANDIDATE_COLUMNS,
    candidate_rows_for_requirements,
    evidence_paths,
    get_solution_steps,
    heterologous_requirements,
    merge_step_candidates,
    read_manifest,
    recommend_uniprot_proteins,
    rel_or_abs,
    unique_main_ecs,
    write_csv,
    write_json_atomic,
)
from src.main_protein_selection.biochemical_realizability import (
    candidate_is_reaction_verified,
    ec_status,
    route_repair_requests_for_requirement,
)
from src.main_protein_selection.reaction_aware_retrieval import (
    RheaClient,
    reaction_evidence_for_requirements,
    retrieve_rhea_candidates_for_requirement,
)
from src.main_protein_selection.reaction_direction_verifier import (
    DIRECTION_CONTRADICTED,
    DIRECTION_UNKNOWN,
    DirectionEvidenceClient,
    direction_evidence_artifact,
    enrich_requirements_with_direction_context,
)
from src.main_protein_selection.kegg_ko_retrieval import (
    KeggKOClient,
    KeggKOSourceUnavailable,
    retrieve_ko_candidates,
)
from src.main_protein_selection.literature_activity import (
    run_literature_activity_search,
    write_source_unavailable_artifact,
)
from src.main_protein_selection.models import (
    MainEnzymeCandidate,
    MainEnzymeSelectionParameters,
    MainEnzymeSelectionResult,
)
from src.main_protein_selection.provenance import solution_fingerprint
from src.main_protein_selection.selenzyme_retrieval import (
    COMPLETE_EC_PATTERN,
    SelenzymeClient,
    SelenzymeSourceUnavailable,
    chassis_host_taxon_id,
    retrieve_selenzyme_candidates,
    selenzyme_target_count,
)


def _ko_requirements_for_retrieval(
    requirements: list[dict],
    verified_step_indexes: set[int],
) -> list[dict]:
    """Use KO only as a fallback source for steps without a verified main enzyme."""

    return [
        requirement
        for requirement in requirements
        if requirement.get("ko_ids")
        and int(requirement.get("step_index") or 0) not in verified_step_indexes
    ]


def _selenzyme_query_specs(requirement: dict[str, Any]) -> list[dict[str, str]]:
    """Build deterministic KEGG-reaction and complete-EC fallback queries."""

    reaction_id = str(requirement.get("reaction_id") or "").strip().upper()
    specs: list[dict[str, str]] = []
    if reaction_id:
        specs.append({
            "query_type": "selenzyme_by_kegg_reaction",
            "query_database": "kegg",
            "query_value": reaction_id,
            "reaction_id": reaction_id,
            "ec_number": "",
        })
    if ec_status(requirement) == "complete":
        seen: set[str] = set()
        for ec_number in (
            requirement.get("ec_numbers")
            or requirement.get("locked_ec_numbers")
            or []
        ):
            normalized = str(ec_number or "").strip()
            if (
                COMPLETE_EC_PATTERN.fullmatch(normalized)
                and normalized not in seen
            ):
                seen.add(normalized)
                specs.append({
                    "query_type": "selenzyme_by_ec_number",
                    "query_database": "ec",
                    "query_value": normalized,
                    "reaction_id": reaction_id,
                    "ec_number": normalized,
                })
    return specs


def _candidates_for_direction_analysis(
    requirements: list[dict],
    candidates_by_ec: dict[str, list],
    candidates_by_step: dict[int, list],
) -> dict[int, list]:
    result: dict[int, list] = {}
    for requirement in requirements:
        step_index = int(requirement.get("step_index") or 0)
        values = list(candidates_by_step.get(step_index, []))
        for ec_number in requirement.get("ec_numbers", []):
            values.extend(candidates_by_ec.get(ec_number, []))
        seen: set[tuple[str, str]] = set()
        result[step_index] = []
        for item in values:
            key = (
                str(getattr(item, "accession", "")),
                str(getattr(item, "retrieval_strategy", "")),
            )
            if key not in seen:
                seen.add(key)
                result[step_index].append(item)
    return result


def _route_repair_requests_v2(
    requirements: list[dict],
    candidate_rows: list[dict],
    solution_id: int,
    source_unavailable: bool = False,
) -> list[dict]:
    """生成可供通路自动回选消费的稳定修复请求模式。"""

    covered_steps = {
        int(row.get("step_index") or 0)
        for row in candidate_rows
        if candidate_is_reaction_verified(row)
    }
    requests: list[dict] = []
    for requirement in requirements:
        step_index = int(requirement.get("step_index") or 0)
        reaction_id = str(requirement.get("reaction_id") or "")
        for request in route_repair_requests_for_requirement(requirement):
            requests.append({
                **request,
                "solution_id": int(solution_id),
                "reason_code": str(
                    request.get("blocking_rule_id")
                    or "biochemical_direction_conflict"
                ),
                "retryable": False,
                "suggested_action": str(
                    request.get("recommended_action")
                    or "replace the reaction and rerun pathway validation"
                ),
            })
        if step_index in covered_steps:
            continue
        step_rows = [
            row for row in candidate_rows
            if int(row.get("step_index") or 0) == step_index
        ]
        direction_conflict = any(
            str(row.get("direction_verdict") or "") == DIRECTION_CONTRADICTED
            or "protein_direction_unsupported" in str(
                row.get("reaction_fit_rule_ids") or ""
            )
            for row in step_rows
        )
        reason_code = (
            "protein_evidence_source_unavailable"
            if source_unavailable
            else "protein_direction_unsupported"
            if direction_conflict
            else "no_exact_reaction_candidate"
        )
        requests.append({
            "solution_id": int(solution_id),
            "step_index": step_index,
            "reaction_id": reaction_id,
            "reason_code": reason_code,
            "retryable": bool(source_unavailable),
            "suggested_action": (
                "retry after the protein evidence source recovers"
                if source_unavailable
                else "reject this solution and evaluate the next flux-pass route"
            ),
            "status": "source_blocked" if source_unavailable else "route_reselection_required",
            "requires_gem_revalidation": False,
        })
    unique: list[dict] = []
    seen: set[tuple] = set()
    for request in requests:
        key = (
            request.get("solution_id"),
            request.get("step_index"),
            request.get("reaction_id"),
            request.get("reason_code"),
        )
        if key not in seen:
            seen.add(key)
            unique.append(request)
    return unique


def select_main_enzymes(
    *,
    manifest_path: str | Path,
    output_dir: str | Path,
    cache_dir: str | Path,
    chassis_key: str = "ecoli_mg1655",
    top_n: int = 5,
    max_results: int = 1000,
    allow_transmembrane: bool = False,
    fetch_proteins: bool = True,
    literature_search: bool = False,
) -> dict:
    """
    为已选 solution 的异源步骤选择主反应酶候选。

    限制：只写主酶候选与检索证据，不写 design_manifest.json。
    """

    try:
        if top_n < 1:
            raise ValueError("top_n must be at least 1")
        output_path = Path(output_dir).expanduser().resolve()
        cache_path = Path(cache_dir).expanduser().resolve()
        paths = evidence_paths(output_path)
        manifest = read_manifest(manifest_path)
        solution_id, steps = get_solution_steps(manifest)
        requirements = heterologous_requirements(steps)
        all_ecs = unique_main_ecs(requirements)
        ecs = list(dict.fromkeys(
            ec_number
            for requirement in requirements
            if ec_status(requirement) == "complete"
            for ec_number in requirement.get("ec_numbers", [])
        ))
        session = requests.Session()
        query_errors: dict[str, str] = {}
        ko_query_errors: dict[str, str] = {}
        candidates_by_ec = {}
        candidates_by_step: dict[int, list] = {
            int(requirement.get("step_index") or 0): []
            for requirement in requirements
        }

        for ec_number in ecs:
            if not fetch_proteins:
                candidates_by_ec[ec_number] = []
                continue
            try:
                candidates_by_ec[ec_number] = recommend_uniprot_proteins(
                    ec_number=ec_number,
                    chassis_key=chassis_key,
                    top_n=top_n,
                    max_results=max_results,
                    allow_transmembrane=allow_transmembrane,
                    session=session,
                )
            except Exception as exc:
                candidates_by_ec[ec_number] = []
                query_errors[ec_number] = str(exc)

        reaction_evidence = []
        reaction_query_ids: list[str] = []
        ko_evidence: list[dict] = []
        ko_query_ids: list[str] = []
        ko_results_by_id: dict[str, dict] = {}
        literature_status = "disabled"
        literature_query_errors: dict[str, str] = {}
        literature_evidence_count = 0
        literature_candidate_count = 0
        selenzyme_evidence: list[dict] = []
        selenzyme_query_ids: list[str] = []
        selenzyme_results_by_query: dict[tuple[str, str], dict] = {}
        uniprot_entry_cache: dict[str, dict | None] = {}
        selenzyme_source_unavailable = False
        selenzyme_circuit_error = ""
        selenzyme_client: SelenzymeClient | None = None
        fallback_requirements: list[dict] = []
        direction_context: list[dict] = []
        direction_client = DirectionEvidenceClient(
            session=session,
            cache_root=cache_path / "direction",
        )
        if fetch_proteins:
            reaction_evidence = reaction_evidence_for_requirements(
                requirements,
                client=RheaClient(
                    session=session,
                    cache_root=cache_path / "rhea",
                ),
            )
            if any(
                requirement.get("rhea_ids") or requirement.get("rhea_master_ids")
                for requirement in requirements
            ):
                direction_context = enrich_requirements_with_direction_context(
                    requirements,
                    direction_client,
                )
            else:
                for requirement in requirements:
                    requirement["direction_evidence_status"] = "unknown_no_rhea_mapping"
                direction_context = [{
                    "step_index": int(requirement.get("step_index") or 0),
                    "reaction_id": str(requirement.get("reaction_id") or ""),
                    "status": "unknown_no_rhea_mapping",
                } for requirement in requirements]
            host_taxon_id = chassis_host_taxon_id(chassis_key)
            targets = selenzyme_target_count(top_n)
            # First exhaust the normal complete-EC/Rhea evidence path.  Exact
            # KO annotations are then resolved before SelenzymeRF is used as
            # the final per-step fallback.
            for requirement in requirements:
                step_index = int(requirement.get("step_index") or 0)
                if ec_status(requirement) != "complete":
                    continue
                reaction_candidates, query_ids, errors = retrieve_rhea_candidates_for_requirement(
                    requirement,
                    chassis_key,
                    max_results=max_results,
                    top_n=top_n,
                    allow_transmembrane=allow_transmembrane,
                    session=session,
                )
                candidates_by_step[step_index].extend(reaction_candidates)
                reaction_query_ids.extend(query_ids)
                query_errors.update(errors)

            preliminary_rows = candidate_rows_for_requirements(
                requirements,
                candidates_by_ec,
                candidates_by_step,
            )
            verified_before_fallback = {
                int(row.get("step_index") or 0)
                for row in preliminary_rows
                if candidate_is_reaction_verified(row)
            }
            ko_requirements = _ko_requirements_for_retrieval(
                requirements, verified_before_fallback
            )

            ko_client = (
                KeggKOClient(
                    session=session,
                    cache_root=cache_path / "kegg_ko",
                )
                if ko_requirements
                else None
            )
            for requirement in ko_requirements:
                step_index = int(requirement.get("step_index") or 0)
                reaction_id = str(requirement.get("reaction_id") or "").strip().upper()
                annotation_status = ec_status(requirement)
                for raw_ko_id in requirement.get("ko_ids", []):
                    ko_id = str(raw_ko_id or "").strip().upper().removeprefix("KO:")
                    if not ko_id:
                        continue
                    try:
                        if ko_id not in ko_results_by_id:
                            ko_results_by_id[ko_id] = ko_client.proteins_for_ko(ko_id)
                        query_result = ko_results_by_id[ko_id]
                        reaction_candidates, audit_rows, query_ids, errors = (
                            retrieve_ko_candidates(
                                requirement,
                                query_result,
                                chassis_key,
                                top_n=top_n,
                                max_results=max_results,
                                allow_transmembrane=allow_transmembrane,
                                session=session,
                                entry_cache=uniprot_entry_cache,
                            )
                        )
                        candidates_by_step[step_index].extend(reaction_candidates)
                        ko_query_ids.append(str(query_result.get("query_id") or ""))
                        ko_query_ids.extend(query_ids)
                        ko_query_errors.update(errors)
                        ko_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "direction": str(requirement.get("direction") or ""),
                            "ec_status": annotation_status,
                            "ko_id": ko_id,
                            "source_status": str(query_result.get("status") or ""),
                            "query_id": str(query_result.get("query_id") or ""),
                            "cache_hit": bool(query_result.get("cache_hit")),
                            "response_sha256": str(
                                query_result.get("response_sha256") or ""
                            ),
                            "kegg_gene_count": len(query_result.get("gene_ids", [])),
                            "uniprot_mapping_count": len(
                                query_result.get("mappings", [])
                            ),
                            "candidate_count": len(reaction_candidates),
                            "rows": audit_rows,
                            "query_errors": errors,
                            "error": "",
                        })
                    except KeggKOSourceUnavailable as exc:
                        query_id = f"kegg_ko_{ko_id}"
                        ko_query_errors[query_id] = str(exc)
                        ko_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "direction": str(requirement.get("direction") or ""),
                            "ec_status": annotation_status,
                            "ko_id": ko_id,
                            "source_status": "source_unavailable",
                            "query_id": query_id,
                            "cache_hit": False,
                            "rows": [],
                            "error": str(exc),
                        })
                    except Exception as exc:
                        query_id = f"kegg_ko_{ko_id}"
                        ko_query_errors[query_id] = str(exc)
                        ko_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "direction": str(requirement.get("direction") or ""),
                            "ec_status": annotation_status,
                            "ko_id": ko_id,
                            "source_status": "retrieval_error",
                            "query_id": query_id,
                            "cache_hit": False,
                            "rows": [],
                            "error": str(exc),
                        })

            rows_after_ko = candidate_rows_for_requirements(
                requirements,
                candidates_by_ec,
                candidates_by_step,
            )
            verified_after_ko = {
                int(row.get("step_index") or 0)
                for row in rows_after_ko
                if candidate_is_reaction_verified(row)
            }
            literature_requirements = [
                requirement
                for requirement in requirements
                if int(requirement.get("step_index") or 0)
                not in verified_after_ko
            ]
            try:
                literature_result = run_literature_activity_search(
                    literature_requirements,
                    enabled=literature_search,
                    output_dir=output_path,
                    cache_dir=cache_path,
                    chassis_key=chassis_key,
                    top_n=top_n,
                    max_results=min(max_results, 25),
                    allow_transmembrane=allow_transmembrane,
                    session=session,
                )
                literature_status = str(literature_result.status)
                literature_query_errors.update(literature_result.query_errors)
                literature_evidence_count = int(
                    literature_result.artifact.summary.evidence_count
                )
                for step_index, literature_candidates in (
                    literature_result.candidates_by_step.items()
                ):
                    candidates_by_step.setdefault(int(step_index), []).extend(
                        literature_candidates
                    )
                    literature_candidate_count += len(literature_candidates)
            except Exception as exc:
                # Literature research is an optional evidence source.  Its
                # failure must never discard standard candidates or prevent
                # the deterministic Selenzyme fallback from running.
                failure_result = write_source_unavailable_artifact(
                    literature_requirements,
                    output_dir=output_path,
                    chassis_key=chassis_key,
                    message=f"{type(exc).__name__}: {exc}",
                    top_n=top_n,
                    max_results=min(max_results, 25),
                    allow_transmembrane=allow_transmembrane,
                )
                literature_status = str(failure_result.status)
                literature_query_errors.update(failure_result.query_errors)
                literature_evidence_count = int(
                    failure_result.artifact.summary.evidence_count
                )

            rows_after_literature = candidate_rows_for_requirements(
                requirements,
                candidates_by_ec,
                candidates_by_step,
            )
            verified_after_literature = {
                int(row.get("step_index") or 0)
                for row in rows_after_literature
                if candidate_is_reaction_verified(row)
            }
            fallback_requirements = [
                requirement
                for requirement in requirements
                if int(requirement.get("step_index") or 0)
                not in verified_after_literature
            ]

            if fallback_requirements:
                try:
                    selenzyme_client = SelenzymeClient(
                        session=session,
                        cache_root=cache_path / "selenzyme",
                    )
                except Exception as exc:
                    selenzyme_source_unavailable = True
                    selenzyme_circuit_error = str(exc)

            for requirement in fallback_requirements:
                step_index = int(requirement.get("step_index") or 0)
                annotation_status = ec_status(requirement)
                reaction_id = str(requirement.get("reaction_id") or "").strip().upper()
                query_specs = _selenzyme_query_specs(requirement)
                if not query_specs:
                    selenzyme_evidence.append({
                        "step_index": step_index,
                        "reaction_id": reaction_id,
                        "ec_number": "",
                        "ec_status": annotation_status,
                        "query_type": "",
                        "query_database": "",
                        "query_value": "",
                        "source_status": "no_query_input",
                        "query_id": "",
                        "rows": [],
                        "error": "Missing KEGG reaction ID and complete EC number",
                        "circuit_open": False,
                    })
                    continue
                if selenzyme_source_unavailable:
                    for spec in query_specs:
                        selenzyme_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "ec_number": spec["ec_number"],
                            "ec_status": annotation_status,
                            "query_type": spec["query_type"],
                            "query_database": spec["query_database"],
                            "query_value": spec["query_value"],
                            "source_status": "source_unavailable",
                            "query_id": (
                                f"selenzyme_{spec['query_database']}_"
                                f"{spec['query_value']}"
                            ),
                            "rows": [],
                            "error": selenzyme_circuit_error,
                            "circuit_open": True,
                        })
                    continue
                for spec in query_specs:
                    if selenzyme_source_unavailable:
                        selenzyme_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "ec_number": spec["ec_number"],
                            "ec_status": annotation_status,
                            "query_type": spec["query_type"],
                            "query_database": spec["query_database"],
                            "query_value": spec["query_value"],
                            "source_status": "source_unavailable",
                            "query_id": (
                                f"selenzyme_{spec['query_database']}_"
                                f"{spec['query_value']}"
                            ),
                            "rows": [],
                            "error": selenzyme_circuit_error,
                            "circuit_open": True,
                        })
                        continue
                    try:
                        if selenzyme_client is None:
                            raise SelenzymeSourceUnavailable(
                                selenzyme_circuit_error
                                or "Selenzyme client is unavailable"
                            )
                        query_key = (
                            spec["query_database"],
                            spec["query_value"],
                        )
                        if query_key not in selenzyme_results_by_query:
                            if spec["query_database"] == "ec":
                                query_result = selenzyme_client.query_ec_number(
                                    spec["query_value"],
                                    host_taxon_id=host_taxon_id,
                                    targets=targets,
                                )
                            else:
                                query_result = selenzyme_client.query_kegg_reaction(
                                    spec["query_value"],
                                    host_taxon_id=host_taxon_id,
                                    targets=targets,
                                )
                            selenzyme_results_by_query[query_key] = query_result
                        query_result = selenzyme_results_by_query[query_key]
                        reaction_candidates, audit_rows, query_ids, errors = (
                            retrieve_selenzyme_candidates(
                                requirement,
                                query_result,
                                chassis_key,
                                top_n=top_n,
                                allow_transmembrane=allow_transmembrane,
                                session=session,
                                entry_cache=uniprot_entry_cache,
                            )
                        )
                        candidates_by_step[step_index].extend(reaction_candidates)
                        selenzyme_query_ids.append(
                            str(query_result.get("query_id") or "")
                        )
                        selenzyme_query_ids.extend(query_ids)
                        query_errors.update(errors)
                        selenzyme_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "ec_number": spec["ec_number"],
                            "direction": str(requirement.get("direction") or ""),
                            "ec_status": annotation_status,
                            "query_type": spec["query_type"],
                            "query_database": spec["query_database"],
                            "query_value": spec["query_value"],
                            "source_status": str(query_result.get("status") or ""),
                            "query_id": str(query_result.get("query_id") or ""),
                            "app": str(query_result.get("app") or ""),
                            "version": str(query_result.get("version") or ""),
                            "cache_hit": bool(query_result.get("cache_hit")),
                            "response_sha256": str(
                                query_result.get("response_sha256") or ""
                            ),
                            "rows": audit_rows,
                            "error": "",
                            "circuit_open": False,
                        })
                    except SelenzymeSourceUnavailable as exc:
                        selenzyme_source_unavailable = True
                        selenzyme_circuit_error = str(exc)
                        query_id = (
                            f"selenzyme_{spec['query_database']}_"
                            f"{spec['query_value']}"
                        )
                        query_errors[query_id] = str(exc)
                        selenzyme_evidence.append({
                            "step_index": step_index,
                            "reaction_id": reaction_id,
                            "ec_number": spec["ec_number"],
                            "direction": str(requirement.get("direction") or ""),
                            "ec_status": annotation_status,
                            "query_type": spec["query_type"],
                            "query_database": spec["query_database"],
                            "query_value": spec["query_value"],
                            "source_status": "source_unavailable",
                            "query_id": query_id,
                            "rows": [],
                            "error": str(exc),
                            "circuit_open": True,
                        })
        else:
            literature_result = run_literature_activity_search(
                requirements,
                enabled=False,
                output_dir=output_path,
                cache_dir=cache_path,
                chassis_key=chassis_key,
                top_n=top_n,
                max_results=min(max_results, 25),
                allow_transmembrane=allow_transmembrane,
                session=session,
            )
            literature_status = str(literature_result.status)
            literature_evidence_count = int(
                literature_result.artifact.summary.evidence_count
            )
            reaction_evidence = [{
                "step_index": int(requirement.get("step_index") or 0),
                "reaction_id": str(requirement.get("reaction_id") or ""),
                "source_status": "fetch_disabled",
                "rhea_ids": [],
            } for requirement in requirements]
            for requirement in requirements:
                requirement["direction_evidence_status"] = "fetch_disabled"
            direction_context = [{
                "step_index": int(requirement.get("step_index") or 0),
                "reaction_id": str(requirement.get("reaction_id") or ""),
                "status": "fetch_disabled",
            } for requirement in requirements]
            ko_evidence = [{
                "step_index": int(requirement.get("step_index") or 0),
                "reaction_id": str(requirement.get("reaction_id") or ""),
                "ec_status": ec_status(requirement),
                "ko_ids": list(requirement.get("ko_ids") or []),
                "source_status": "fetch_disabled",
                "rows": [],
            } for requirement in requirements if requirement.get("ko_ids")]
            selenzyme_evidence = [{
                "step_index": int(requirement.get("step_index") or 0),
                "reaction_id": str(requirement.get("reaction_id") or ""),
                "ec_status": ec_status(requirement),
                "source_status": "fetch_disabled",
                "rows": [],
            } for requirement in requirements]

        step_rows = candidate_rows_for_requirements(
            requirements,
            candidates_by_ec,
            candidates_by_step,
            top_n=top_n,
        )
        verified_step_rows = [
            row for row in step_rows if candidate_is_reaction_verified(row)
        ]
        selected_step_rows = [
            row
            for row in step_rows
            if str(row.get("selection_status") or "") == "selected"
        ]
        merged_rows = merge_step_candidates(selected_step_rows)
        route_repair_requests = _route_repair_requests_v2(
            requirements,
            step_rows,
            solution_id,
            source_unavailable=selenzyme_source_unavailable,
        )
        write_csv(
            paths["step_main_enzyme_candidates_csv"],
            selected_step_rows,
            STEP_CANDIDATE_COLUMNS,
        )
        write_csv(
            paths["step_main_enzyme_candidate_audit_csv"],
            step_rows,
            STEP_CANDIDATE_COLUMNS,
        )
        write_csv(paths["main_enzyme_candidates_csv"], merged_rows, PROTEIN_CANDIDATE_COLUMNS)
        write_json_atomic(paths["reaction_evidence_json"], {
            "schema_version": "reaction_evidence.v1",
            "selected_solution_id": solution_id,
            "evidence": reaction_evidence,
            "query_ids": list(dict.fromkeys(reaction_query_ids)),
        })
        direction_artifact = direction_evidence_artifact(
            solution_id=solution_id,
            requirements=requirements,
            context_records=direction_context,
            candidates_by_step=_candidates_for_direction_analysis(
                requirements, candidates_by_ec, candidates_by_step
            ),
            client=direction_client,
        )
        write_json_atomic(paths["direction_evidence_json"], direction_artifact)
        write_json_atomic(paths["ko_evidence_json"], {
            "schema_version": "ko_evidence.v1",
            "selected_solution_id": solution_id,
            "policy": {
                "scope": "unverified_steps_with_kegg_ko",
                "acceptance": "exact_requirement_ko_intersection",
                "fallback": "literature_activity_then_selenzyme_rf",
            },
            "evidence": ko_evidence,
            "query_ids": list(dict.fromkeys(
                value for value in ko_query_ids if value
            )),
        })
        write_json_atomic(paths["selenzyme_evidence_json"], {
            "schema_version": "selenzyme_evidence.v4",
            "selected_solution_id": solution_id,
            "policy": {
                "ec_scope": "all_steps_without_verified_ec_rhea_ko_or_literature_candidate",
                "auto_accept": "kegg_query_valid_combined_reaction_similarity",
                "exact_match": "kegg_query_combined_reaction_similarity_equal_1",
                "risk_fallback": "kegg_similarity_or_complete_ec_association_requires_review",
                "exact_similarity_tolerance": 1e-6,
                "query_input": "kegg_reaction_id_and_complete_ec_numbers",
                "ec_query_policy": {
                    "candidate_status": "verified_with_risk",
                    "current_uniprot_ec_match_score": 69.0,
                    "shared_reaction_ec_overlap_score": 60.0,
                    "missing_current_uniprot_ec_score": 55.0,
                    "unrelated_current_uniprot_ec": "rejected",
                    "unit_similarity_is_not_locked_reaction_exactness": True,
                },
                "no_msa": True,
            },
            "evidence": selenzyme_evidence,
            "query_ids": list(dict.fromkeys(
                value for value in selenzyme_query_ids if value
            )),
        })
        write_json_atomic(paths["route_repair_requests_json"], {
            "schema_version": "route_repair_requests.v2",
            "selected_solution_id": solution_id,
            "requests": route_repair_requests,
        })

        covered_steps = {
            int(row.get("step_index") or 0) for row in verified_step_rows
        }
        expected_steps = {
            int(requirement.get("step_index") or 0) for requirement in requirements
        }
        rank_one_rows = [
            row for row in selected_step_rows
            if int(row.get("candidate_rank") or 0) == 1
        ]
        result = {
            "ok": not selenzyme_source_unavailable,
            "status": (
                "source_unavailable"
                if selenzyme_source_unavailable
                else "complete"
            ),
            "selected_solution_id": solution_id,
            "chassis_key": chassis_key,
            "ec_numbers": all_ecs,
            "complete_ec_query_numbers": ecs,
            "heterologous_step_count": len(requirements),
            "step_candidate_count": len(selected_step_rows),
            "evaluated_step_candidate_count": len(step_rows),
            "verified_step_candidate_count": len(verified_step_rows),
            "rejected_or_review_step_candidate_count": len(step_rows) - len(verified_step_rows),
            "protein_candidate_count": len(merged_rows),
            "reaction_evidence_count": len(reaction_evidence),
            "ko_evidence_count": len(ko_evidence),
            "ko_exact_candidate_count": sum(
                1
                for row in step_rows
                if str(row.get("reaction_confidence") or "") == "ko_exact"
            ),
            "literature_search_enabled": bool(literature_search and fetch_proteins),
            "literature_status": literature_status,
            "literature_evidence_count": literature_evidence_count,
            "literature_candidate_count": literature_candidate_count,
            "literature_query_errors": literature_query_errors,
            "selenzyme_evidence_count": len(selenzyme_evidence),
            "selenzyme_exact_candidate_count": sum(
                1
                for row in step_rows
                if str(row.get("reaction_confidence") or "") == "selenzyme_exact"
            ),
            "selenzyme_risk_candidate_count": sum(
                1
                for row in step_rows
                if str(row.get("reaction_confidence") or "")
                in {"selenzyme_risk", "selenzyme_ec_risk"}
            ),
            "selenzyme_ec_risk_candidate_count": sum(
                1
                for row in step_rows
                if str(row.get("reaction_confidence") or "")
                == "selenzyme_ec_risk"
            ),
            "route_repair_request_count": len(route_repair_requests),
            "uncovered_step_indexes": sorted(expected_steps - covered_steps),
            "direction_rejected_step_indexes": sorted({
                int(requirement.get("step_index") or 0)
                for requirement in requirements
                if int(requirement.get("step_index") or 0) not in covered_steps
                and any(
                    int(row.get("step_index") or 0) == int(requirement.get("step_index") or 0)
                    and str(row.get("direction_verdict") or "") == DIRECTION_CONTRADICTED
                    for row in step_rows
                )
            }),
            "direction_risk_step_indexes": sorted({
                int(row.get("step_index") or 0)
                for row in rank_one_rows
                if str(row.get("direction_verdict") or "") == DIRECTION_UNKNOWN
            }),
            "query_errors": query_errors,
            "ko_query_errors": ko_query_errors,
            "step_main_enzyme_candidates_csv": rel_or_abs(paths["step_main_enzyme_candidates_csv"]),
            "step_main_enzyme_candidate_audit_csv": rel_or_abs(
                paths["step_main_enzyme_candidate_audit_csv"]
            ),
            "main_enzyme_candidates_csv": rel_or_abs(paths["main_enzyme_candidates_csv"]),
            "reaction_evidence_json": rel_or_abs(paths["reaction_evidence_json"]),
            "direction_evidence_json": rel_or_abs(paths["direction_evidence_json"]),
            "ko_evidence_json": rel_or_abs(paths["ko_evidence_json"]),
            "literature_activity_evidence_json": rel_or_abs(
                paths["literature_activity_evidence_json"]
            ),
            "literature_activity_evidence_csv": rel_or_abs(
                paths["literature_activity_evidence_csv"]
            ),
            "selenzyme_evidence_json": rel_or_abs(paths["selenzyme_evidence_json"]),
            "route_repair_requests_json": rel_or_abs(paths["route_repair_requests_json"]),
        }
        candidates_by_step: dict[int, list[MainEnzymeCandidate]] = {}
        for row in selected_step_rows:
            candidate = MainEnzymeCandidate.from_candidate_row(row)
            candidates_by_step.setdefault(candidate.step_index, []).append(candidate)
        for candidates in candidates_by_step.values():
            candidates.sort(key=lambda item: item.candidate_rank)

        solution = manifest.get("solution")
        expansion_depth = int(
            solution.get("expansion_depth") or 0
            if isinstance(solution, dict)
            else 0
        )
        canonical_result = MainEnzymeSelectionResult(
            ok=result["ok"],
            status=result["status"],
            selected_solution_id=solution_id,
            expansion_depth=expansion_depth,
            solution_fingerprint=solution_fingerprint(solution_id, steps),
            chassis_key=chassis_key,
            parameters=MainEnzymeSelectionParameters(
                top_n=top_n,
                max_results=max_results,
                allow_transmembrane=allow_transmembrane,
                fetch_proteins=fetch_proteins,
                literature_search=literature_search,
            ),
            shortlist_decision_fingerprint=shortlist_decision_fingerprint(
                selected_step_rows
            ),
            candidates_by_step=candidates_by_step,
            uncovered_step_indexes=result["uncovered_step_indexes"],
            direction_rejected_step_indexes=result[
                "direction_rejected_step_indexes"
            ],
            direction_risk_step_indexes=result["direction_risk_step_indexes"],
            evidence_files={
                key: rel_or_abs(path)
                for key, path in paths.items()
                if key != "main_enzyme_selection_json"
            },
        )
        write_json_atomic(
            paths["main_enzyme_selection_json"],
            canonical_result.model_dump(mode="json"),
        )
        result["main_enzyme_selection_json"] = rel_or_abs(
            paths["main_enzyme_selection_json"]
        )
        result["schema_version"] = canonical_result.schema_version
        return result
    except Exception:
        raise


def run_main_protein_selection(config: Any, **selection_options: Any) -> dict:
    """Run main-enzyme selection using this project's ``RunConfig`` paths."""

    if hasattr(config, "top_n"):
        selection_options.setdefault("top_n", config.top_n)
    selection_options.setdefault(
        "literature_search", bool(getattr(config, "literature_search", False))
    )

    result = select_main_enzymes(
        manifest_path=Path(config.manifest_output_path),
        output_dir=Path(config.project_output_path) / "main_protein_selection",
        cache_dir=Path(config.cache_dir) / "main_protein_selection",
        **selection_options,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = ["run_main_protein_selection", "select_main_enzymes"]
