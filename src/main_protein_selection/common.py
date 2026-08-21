"""Shared I/O and candidate projection for main-enzyme selection."""

from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path
from typing import Any

from src.main_protein_selection.biochemical_realizability import (
    evaluate_candidate_reaction_fit,
)
from src.main_protein_selection.reaction_direction_verifier import (
    direction_decision_for_candidate,
)
from src.main_protein_selection.auxiliary_roles import (
    annotate_candidate_auxiliary_roles,
    merge_auxiliary_requirements,
)
from src.main_protein_selection.uniprot_protein_candidates import (
    PROTEIN_CANDIDATE_COLUMNS,
    STEP_CANDIDATE_COLUMNS,
    ProteinCandidate,
    _join,
    _json,
    _split_list_field,
    _step_candidate_row,
    _unique,
    recommend_uniprot_proteins,
)


STEP_MAIN_CANDIDATES_FILENAME = "step_main_enzyme_candidates.csv"
STEP_MAIN_CANDIDATE_AUDIT_FILENAME = "step_main_enzyme_candidate_audit.csv"
MAIN_CANDIDATES_FILENAME = "main_enzyme_candidates.csv"
MAIN_ENZYME_SELECTION_FILENAME = "main_enzyme_selection.json"
REACTION_EVIDENCE_FILENAME = "reaction_evidence.json"
DIRECTION_EVIDENCE_FILENAME = "direction_evidence.json"
KO_EVIDENCE_FILENAME = "ko_evidence.json"
SELENZYME_EVIDENCE_FILENAME = "selenzyme_evidence.json"
LITERATURE_ACTIVITY_EVIDENCE_FILENAME = "literature_activity_evidence.json"
LITERATURE_ACTIVITY_EVIDENCE_CSV_FILENAME = "literature_activity_evidence.csv"
ROUTE_REPAIR_REQUESTS_FILENAME = "route_repair_requests.json"


def read_manifest(path: str | Path) -> dict[str, Any]:
    manifest_path = Path(path)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")
    value = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Manifest root must be an object: {manifest_path}")
    return value


def write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=str(path.parent), delete=False, suffix=".tmp"
    ) as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: _json(row.get(key, ""))
                if isinstance(row.get(key), (dict, list))
                else row.get(key, "")
                for key in columns
            })


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def output_root(output_dir: str | Path, create: bool = False) -> Path:
    root = Path(output_dir)
    if create:
        root.mkdir(parents=True, exist_ok=True)
    return root


def evidence_paths(output_dir: str | Path) -> dict[str, Path]:
    root = output_root(output_dir)
    return {
        "step_main_enzyme_candidates_csv": root / STEP_MAIN_CANDIDATES_FILENAME,
        "step_main_enzyme_candidate_audit_csv": (
            root / STEP_MAIN_CANDIDATE_AUDIT_FILENAME
        ),
        "main_enzyme_candidates_csv": root / MAIN_CANDIDATES_FILENAME,
        "main_enzyme_selection_json": root / MAIN_ENZYME_SELECTION_FILENAME,
        "reaction_evidence_json": root / REACTION_EVIDENCE_FILENAME,
        "direction_evidence_json": root / DIRECTION_EVIDENCE_FILENAME,
        "ko_evidence_json": root / KO_EVIDENCE_FILENAME,
        "literature_activity_evidence_json": (
            root / LITERATURE_ACTIVITY_EVIDENCE_FILENAME
        ),
        "literature_activity_evidence_csv": (
            root / LITERATURE_ACTIVITY_EVIDENCE_CSV_FILENAME
        ),
        "selenzyme_evidence_json": root / SELENZYME_EVIDENCE_FILENAME,
        "route_repair_requests_json": root / ROUTE_REPAIR_REQUESTS_FILENAME,
    }


def rel_or_abs(path: Path) -> str:
    return str(path.resolve())


def get_solution_steps(manifest: dict[str, Any]) -> tuple[int, list[dict[str, Any]]]:
    solution = manifest.get("solution")
    if not isinstance(solution, dict):
        raise ValueError('manifest["solution"] is missing')
    solution_id = int(solution.get("solution_id") or 0)
    raw_steps = solution.get("steps")
    if solution_id < 1 or not isinstance(raw_steps, list) or not raw_steps:
        raise ValueError('manifest["solution"] has no selected route')
    return solution_id, [
        {**step, "solution_id": solution_id}
        for step in raw_steps
        if isinstance(step, dict)
    ]


def is_spontaneous_step(step: dict[str, Any]) -> bool:
    annotation = " ".join(
        str(step.get(key) or "")
        for key in ("reaction_comment", "resolution_action", "reaction_name")
    ).lower()
    return any(marker in annotation for marker in (
        "spontaneous", "non-enzymatic", "nonenzymatic", "non enzymatic"
    ))


def heterologous_requirements(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for step in steps:
        if str(step.get("status") or "").strip().lower() != "heterologous":
            continue
        if is_spontaneous_step(step):
            continue
        result.append({
            "solution_id": int(step.get("solution_id") or 0),
            "step_index": int(step.get("step_index") or 0),
            "reaction_id": str(step.get("reaction_id") or ""),
            "reaction_name": str(step.get("reaction_name") or ""),
            "reaction_comment": str(step.get("reaction_comment") or ""),
            "equation": str(step.get("equation") or ""),
            "direction": str(step.get("direction") or ""),
            "produced_compound_id": str(step.get("produced_compound_id") or ""),
            "produced_compound_name": str(step.get("produced_compound_name") or ""),
            "precursor_compound_ids": _split_list_field(step.get("precursor_compound_ids")),
            "precursor_compound_labels": _split_list_field(step.get("precursor_compound_labels")),
            "ko_ids": _split_list_field(step.get("ko_ids")),
            "module_ids": _split_list_field(step.get("module_ids")),
            "source_reaction_ids": _split_list_field(step.get("source_reaction_ids")),
            "resolution_action": str(step.get("resolution_action") or ""),
            "resolution_evidence": _split_list_field(step.get("resolution_evidence")),
            "rhea_ids": _split_list_field(step.get("rhea_ids")),
            "ec_numbers": _split_list_field(step.get("enzyme_ecs")),
            "locked_ec_numbers": _split_list_field(
                step.get("locked_enzyme_ecs") or step.get("enzyme_ecs")
            ),
            "ec_status": str(step.get("ec_status") or ""),
            "enzyme_search_eligible": str(step.get("enzyme_search_eligible") or ""),
            "auxiliary_requirements": (
                list(step.get("auxiliary_requirements") or [])
                if isinstance(step.get("auxiliary_requirements"), list)
                else []
            ),
        })
    return result


def unique_main_ecs(requirements: list[dict[str, Any]]) -> list[str]:
    return _unique([
        ec for requirement in requirements for ec in requirement.get("ec_numbers", [])
    ])


def candidate_rows_for_requirements(
    requirements: list[dict[str, Any]],
    candidates_by_ec: dict[str, list[ProteinCandidate]],
    candidates_by_step: dict[int, list[ProteinCandidate]] | None = None,
    *,
    top_n: int | None = None,
) -> list[dict[str, Any]]:
    """Project, deduplicate and rank catalytic-main candidates per step.

    A UniProt accession identifies a protein.  The same accession retrieved by
    multiple EC, Rhea, KO or Selenzyme queries therefore remains one candidate
    for a route step, with all retrieval evidence retained on that row.
    """

    if top_n is not None and top_n < 1:
        raise ValueError("top_n must be at least 1")

    by_step = candidates_by_step or {}
    rows: list[dict[str, Any]] = []
    for requirement in requirements:
        candidates: list[tuple[str, ProteinCandidate]] = []
        for ec_number in requirement.get("ec_numbers", []):
            candidates.extend((ec_number, item) for item in candidates_by_ec.get(ec_number, []))
        step_index = int(requirement.get("step_index") or 0)
        for item in by_step.get(step_index, []):
            ecs = item.ec_numbers or requirement.get("locked_ec_numbers", [])
            candidates.append((str(ecs[0]) if ecs else "", item))
        for ec_number, candidate in candidates:
            if str(candidate.candidate_role or "catalytic_main") != "catalytic_main":
                continue
            row = _step_candidate_row(requirement, ec_number, candidate)
            direction = direction_decision_for_candidate(requirement, row)
            row.update({
                "direction_support": direction["verdict"] or row.get("direction_support", ""),
                "direction_verdict": direction["verdict"],
                "direction_confidence": direction["confidence"],
                "direction_evidence_level": direction["evidence_level"],
                "direction_evidence_source_ids": _join(direction["source_ids"]),
                "direction_evidence": " | ".join(direction["evidence"]),
                "required_rhea_direction_ids": _join(
                    direction["required_rhea_direction_ids"]
                ),
            })
            fit = evaluate_candidate_reaction_fit(requirement, row)
            row.update({
                "reaction_fit_status": fit["status"],
                "reaction_fit_score": fit["score"],
                "reaction_fit_rule_ids": _join(fit["rule_ids"]),
                "reaction_fit_evidence": " | ".join(fit["evidence"]),
            })
            annotate_candidate_auxiliary_roles(requirement, row)
            rows.append(row)

        step_rows = [
            row for row in rows if _int(row.get("step_index")) == step_index
        ]
        rows = [
            row for row in rows if _int(row.get("step_index")) != step_index
        ]
        groups: dict[str, list[dict[str, Any]]] = {}
        for row in step_rows:
            accession = str(row.get("accession") or "").strip().upper()
            if not accession:
                raise ValueError(
                    f"Main-enzyme candidate for step {step_index} has no accession"
                )
            groups.setdefault(accession, []).append(row)
        rows.extend(
            _merge_step_accession_rows(requirement, matches)
            for matches in groups.values()
        )

    status_priority = {
        "verified": 0,
        "verified_with_risk": 1,
        "manual_review": 2,
        "rejected": 3,
    }

    def evidence_source_priority(row: dict[str, Any]) -> int:
        """Keep deterministic evidence precedence independent of raw scores.

        Literature-derived activities deliberately remain
        ``verified_with_risk``.  Selenzyme can nevertheless label a unit
        similarity hit ``verified``, so reaction-fit status alone cannot
        express the required source precedence.
        """

        strategies = {
            value.lower()
            for value in _split_list_field(row.get("retrieval_strategy"))
        }
        confidence = {
            value.lower()
            for value in _split_list_field(row.get("reaction_confidence"))
        }
        if strategies & {"ec_exact", "rhea_reaction", "kegg_ko_exact"}:
            return 0
        if "literature_experimental_activity" in strategies:
            if "literature_grade_a" in confidence:
                return 1
            if "literature_grade_b" in confidence:
                return 2
            return 2
        if "selenzyme_kegg_exact" in strategies:
            return 3
        if "selenzyme_ec_risk" in strategies:
            return 4
        if any(value.startswith("selenzyme_") for value in strategies):
            return 5
        return 6

    rows.sort(key=lambda row: (
        _int(row.get("step_index")),
        0
        if str(row.get("reaction_fit_status") or "")
        in {"verified", "verified_with_risk"}
        else 1,
        evidence_source_priority(row),
        status_priority.get(str(row.get("reaction_fit_status") or ""), 9),
        -_float(row.get("reaction_fit_score")),
        -_float(row.get("score")),
        0 if str(row.get("reviewed") or "").lower() in {"true", "1"} else 1,
        str(row.get("accession") or ""),
    ))
    evaluation_ranks: dict[int, int] = {}
    candidate_ranks: dict[int, int] = {}
    for row in rows:
        step_index = _int(row.get("step_index"))
        evaluation_ranks[step_index] = evaluation_ranks.get(step_index, 0) + 1
        row["evaluation_rank"] = evaluation_ranks[step_index]
        status = str(row.get("reaction_fit_status") or "")
        if status in {"verified", "verified_with_risk"}:
            next_rank = candidate_ranks.get(step_index, 0) + 1
            candidate_ranks[step_index] = next_rank
            if top_n is None or next_rank <= top_n:
                row["candidate_rank"] = next_rank
                row["selection_status"] = "selected"
            else:
                row["candidate_rank"] = ""
                row["selection_status"] = "eligible_not_selected"
        else:
            row["candidate_rank"] = ""
            row["selection_status"] = (
                "manual_review" if status == "manual_review" else "rejected"
            )
    return rows


_MERGED_LIST_FIELDS = (
    "ec_numbers",
    "retrieval_strategy",
    "retrieval_query_id",
    "matched_rhea_ids",
    "matched_ko_ids",
    "kegg_gene_ids",
    "direction_evidence_source_ids",
    "direction_evidence",
    "required_rhea_direction_ids",
    "reaction_confidence",
    "enzyme_system_type",
    "required_auxiliary_roles",
    "auxiliary_requirement_status",
    "selenzyme_risk_status",
    "selenzyme_matched_reaction_id",
    "selenzyme_direction_used",
    "selenzyme_direction_preferred",
    "gene_names",
    "catalytic_activities",
    "cofactors",
    "subunit",
    "function_comments",
    "ptm_comments",
    "feature_annotations",
    "domain_ids",
    "keywords",
    "aliases",
    "publication_ids",
    "cross_references",
    "subcellular_locations",
    "taxonomic_lineage",
    "rhea_ids",
    "reasons",
    "warnings",
)


def _values_from_rows(
    rows: list[dict[str, Any]],
    field: str,
) -> list[str]:
    return _unique([
        value
        for row in rows
        for value in _split_list_field(row.get(field))
    ])


def _normalized_sequence(value: Any) -> str:
    return "".join(str(value or "").split()).upper()


def _validate_duplicate_sequences(
    step_index: int,
    accession: str,
    rows: list[dict[str, Any]],
) -> None:
    sequences = {
        sequence
        for row in rows
        if (sequence := _normalized_sequence(row.get("sequence")))
    }
    hashes = {
        str(row.get("sequence_sha256") or "").strip().lower()
        for row in rows
        if str(row.get("sequence_sha256") or "").strip()
    }
    if len(sequences) > 1 or len(hashes) > 1:
        raise ValueError(
            "Conflicting sequence evidence for main-enzyme candidate "
            f"{accession} at step {step_index}"
        )


def _merge_step_accession_rows(
    requirement: dict[str, Any],
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    """Merge all evidence for one accession on one route step."""

    step_index = _int(requirement.get("step_index"))
    accession = str(rows[0].get("accession") or "").strip().upper()
    _validate_duplicate_sequences(step_index, accession, rows)

    best = max(
        rows,
        key=lambda row: (
            _float(row.get("score")),
            1 if str(row.get("reviewed") or "").lower() in {"true", "1"} else 0,
        ),
    )
    merged = dict(best)
    merged["accession"] = accession
    merged["score"] = max(_float(row.get("score")) for row in rows)
    merged["reviewed"] = any(
        str(row.get("reviewed") or "").lower() in {"true", "1"}
        for row in rows
    )
    merged["ec_number"] = _join(_unique([
        value
        for row in rows
        for value in _split_list_field(row.get("ec_number"))
    ]))
    for field in _MERGED_LIST_FIELDS:
        merged[field] = _join(_values_from_rows(rows, field))

    for field in ("sequence", "sequence_sha256"):
        merged[field] = next(
            (
                row.get(field)
                for row in rows
                if str(row.get(field) or "").strip()
            ),
            "",
        )
    sequence_versions = [
        _int(row.get("sequence_version"))
        for row in rows
        if _int(row.get("sequence_version"))
    ]
    merged["sequence_version"] = max(sequence_versions) if sequence_versions else ""

    direction_priority = {"contradicted": 0, "supported": 1, "unknown": 2, "": 3}
    direction_verdict = min(
        (str(row.get("direction_verdict") or "").strip().lower() for row in rows),
        key=lambda value: direction_priority.get(value, 4),
    )
    direction_rows = [
        row
        for row in rows
        if str(row.get("direction_verdict") or "").strip().lower()
        == direction_verdict
    ]
    confidence_priority = {"high": 0, "medium": 1, "low": 2, "": 3}
    merged["direction_support"] = direction_verdict
    merged["direction_verdict"] = direction_verdict
    merged["direction_confidence"] = min(
        (
            str(row.get("direction_confidence") or "").strip().lower()
            for row in direction_rows
        ),
        key=lambda value: confidence_priority.get(value, 4),
        default="",
    )
    merged["direction_evidence_level"] = _join(
        _values_from_rows(direction_rows, "direction_evidence_level")
    )

    ranks = [
        _int(row.get("selenzyme_rank"))
        for row in rows
        if _int(row.get("selenzyme_rank"))
    ]
    merged["selenzyme_rank"] = min(ranks) if ranks else ""
    for field in (
        "selenzyme_score",
        "selenzyme_sim_rf",
        "selenzyme_sim_2018",
        "selenzyme_reaction_similarity",
    ):
        values = [
            _float(row.get(field))
            for row in rows
            if str(row.get(field) or "").strip()
        ]
        merged[field] = max(values) if values else ""
    distances = [
        _float(row.get("selenzyme_taxonomic_distance"))
        for row in rows
        if str(row.get("selenzyme_taxonomic_distance") or "").strip()
    ]
    merged["selenzyme_taxonomic_distance"] = min(distances) if distances else ""

    fit = evaluate_candidate_reaction_fit(requirement, merged)
    fit_options = [
        {
            "status": str(row.get("reaction_fit_status") or ""),
            "score": _float(row.get("reaction_fit_score")),
        }
        for row in rows
    ]
    fit_options.append(fit)
    if direction_verdict == "contradicted":
        selected_fit = fit
    else:
        fit_priority = {
            "verified": 0,
            "verified_with_risk": 1,
            "manual_review": 2,
            "rejected": 3,
        }
        selected_fit = min(
            fit_options,
            key=lambda value: (
                fit_priority.get(str(value.get("status") or ""), 4),
                -_float(value.get("score")),
            ),
        )
    merged["reaction_fit_status"] = selected_fit["status"]
    merged["reaction_fit_score"] = selected_fit["score"]
    merged["reaction_fit_rule_ids"] = _join(_unique([
        *(
            value
            for row in rows
            for value in _split_list_field(row.get("reaction_fit_rule_ids"))
        ),
        *(str(value) for value in fit["rule_ids"]),
    ]))
    merged["reaction_fit_evidence"] = " | ".join(_unique([
        *(
            value
            for row in rows
            for value in _split_list_field(row.get("reaction_fit_evidence"))
        ),
        *(str(value) for value in fit["evidence"]),
    ]))
    return merged


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _int(value: Any) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return 0


def merge_step_candidates(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Merge repeated step hits without introducing component-system semantics."""

    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        accession = str(row.get("accession") or "").strip().upper()
        if accession:
            groups.setdefault(accession, []).append(row)
    merged: list[dict[str, Any]] = []
    for accession, matches in groups.items():
        first = max(
            matches,
            key=lambda row: (_float(row.get("reaction_fit_score")), _float(row.get("score"))),
        )
        scores = [_float(row.get("score")) for row in matches]
        step_indexes = sorted({_int(row.get("step_index")) for row in matches if _int(row.get("step_index"))})
        selection_evidence = []
        risk_evidence = []
        for row in matches:
            similarity = str(row.get("selenzyme_reaction_similarity") or "").strip()
            if similarity:
                selection_evidence.append({
                    "step_index": _int(row.get("step_index")),
                    "reaction_id": str(row.get("reaction_id") or ""),
                    "accession": accession,
                    "rank": _int(row.get("selenzyme_rank")),
                    "reaction_similarity": _float(similarity),
                    "match_type": str(row.get("reaction_confidence") or ""),
                })
            if str(row.get("selenzyme_risk_status") or "").strip():
                risk_evidence.append({
                    "step_index": _int(row.get("step_index")),
                    "reaction_id": str(row.get("reaction_id") or ""),
                    "accession": accession,
                    "risk_status": str(row.get("selenzyme_risk_status") or ""),
                })
        values = lambda field: _unique([
            value for row in matches for value in _split_list_field(row.get(field))
        ])
        ranks = sorted({_int(row.get("selenzyme_rank")) for row in matches if _int(row.get("selenzyme_rank"))})
        similarities = [
            _float(row.get("selenzyme_reaction_similarity"))
            for row in matches if str(row.get("selenzyme_reaction_similarity") or "").strip()
        ]
        auxiliary_requirement_rows: list[dict[str, Any]] = []
        for row in matches:
            try:
                raw_requirements = json.loads(
                    str(row.get("auxiliary_requirements_json") or "[]")
                )
            except json.JSONDecodeError as exc:
                raise ValueError(
                    "invalid auxiliary requirements in candidate rows"
                ) from exc
            for raw_requirement in raw_requirements:
                auxiliary_requirement_rows.append({
                    **dict(raw_requirement),
                    "step_indexes": [_int(row.get("step_index"))],
                    "main_enzyme_accessions": [accession],
                })
        auxiliary_requirements = merge_auxiliary_requirements(
            auxiliary_requirement_rows
        )
        merged.append({
            "accession": accession,
            "entry_name": first.get("entry_name", ""),
            "protein_name": first.get("protein_name", ""),
            "organism_name": first.get("organism_name", ""),
            "organism_id": first.get("organism_id", ""),
            "reviewed": first.get("reviewed", ""),
            "length": first.get("length", ""),
            "covered_step_indexes": _join([str(value) for value in step_indexes]),
            "covered_reaction_ids": _join(_unique([str(row.get("reaction_id") or "") for row in matches])),
            "covered_ec_numbers": _join(_unique([str(row.get("ec_number") or "") for row in matches])),
            "roles": "main",
            "enzyme_system_types": _join(values("enzyme_system_type")),
            "required_auxiliary_roles": _join(
                item["role"] for item in auxiliary_requirements
            ),
            "auxiliary_requirement_statuses": _join(
                values("auxiliary_requirement_status")
            ),
            "auxiliary_requirements_json": _json(auxiliary_requirements),
            "best_score": round(max(scores), 2),
            "mean_score": round(sum(scores) / len(scores), 2),
            "min_score": round(min(scores), 2),
            "match_count": len(matches),
            "candidate_roles": "catalytic_main",
            "retrieval_strategies": _join(values("retrieval_strategy")),
            "retrieval_query_ids": _join(values("retrieval_query_id")),
            "matched_rhea_ids": _join(values("matched_rhea_ids")),
            "matched_ko_ids": _join(values("matched_ko_ids")),
            "kegg_gene_ids": _join(values("kegg_gene_ids")),
            "direction_support": _join(values("direction_support")),
            "direction_verdict": _join(values("direction_verdict")),
            "direction_confidence": _join(values("direction_confidence")),
            "direction_evidence_level": _join(values("direction_evidence_level")),
            "direction_evidence_source_ids": _join(values("direction_evidence_source_ids")),
            "required_rhea_direction_ids": _join(values("required_rhea_direction_ids")),
            "reaction_confidence": _join(values("reaction_confidence")),
            "best_selenzyme_rank": min(ranks) if ranks else "",
            "selenzyme_ranks": _join([str(value) for value in ranks]),
            "best_selenzyme_score": max((_float(row.get("selenzyme_score")) for row in matches), default=0.0) or "",
            "best_selenzyme_reaction_similarity": max(similarities) if similarities else "",
            "selenzyme_match_types": _join([
                value for value in values("reaction_confidence") if value.startswith("selenzyme_")
            ]),
            "selenzyme_selection_evidence_json": _json(selection_evidence),
            "selenzyme_risk_statuses": _join(values("selenzyme_risk_status")),
            "selenzyme_risk_evidence_json": _json(risk_evidence),
            **{
                field: first.get(field, "")
                for field in (
                    "gene_names", "catalytic_activities",
                    "catalytic_activity_records_json", "cofactors", "subunit",
                    "function_comments", "ptm_comments", "feature_annotations",
                    "domain_ids", "keywords", "protein_existence", "sequence_version",
                    "sequence_sha256", "subcellular_locations", "taxonomic_lineage",
                    "rhea_ids", "sequence",
                )
            },
            "warnings": " | ".join(values("warnings")),
        })
    merged.sort(key=lambda row: (
        -len(_split_list_field(row.get("covered_step_indexes"))),
        -_float(row.get("mean_score")),
        str(row.get("accession") or ""),
    ))
    return merged


__all__ = [
    "KO_EVIDENCE_FILENAME",
    "LITERATURE_ACTIVITY_EVIDENCE_CSV_FILENAME",
    "LITERATURE_ACTIVITY_EVIDENCE_FILENAME",
    "MAIN_CANDIDATES_FILENAME",
    "MAIN_ENZYME_SELECTION_FILENAME",
    "PROTEIN_CANDIDATE_COLUMNS",
    "REACTION_EVIDENCE_FILENAME",
    "ROUTE_REPAIR_REQUESTS_FILENAME",
    "SELENZYME_EVIDENCE_FILENAME",
    "STEP_CANDIDATE_COLUMNS",
    "STEP_MAIN_CANDIDATE_AUDIT_FILENAME",
    "STEP_MAIN_CANDIDATES_FILENAME",
    "candidate_rows_for_requirements",
    "evidence_paths",
    "get_solution_steps",
    "heterologous_requirements",
    "is_spontaneous_step",
    "merge_step_candidates",
    "output_root",
    "read_csv",
    "read_manifest",
    "recommend_uniprot_proteins",
    "rel_or_abs",
    "unique_main_ecs",
    "write_csv",
    "write_json_atomic",
    "_split_list_field",
]
