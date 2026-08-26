"""Materialize P5 RetroPath candidates as ordinary, optionally validated solutions."""

from __future__ import annotations

import json
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Mapping, Sequence

from src.pathway_analyze.kegg_gap_analyze import KeggRestClient, gap_depth_output_dir
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEP_COLUMNS,
    CANDIDATE_STEPS_FILE_NAME,
)
from src.pathway_analyze.retropath_promotion import (
    _ELECTRON_STEP_COLUMNS,
    _ELECTRON_SUMMARY_COLUMNS,
    _STEP_BASE_COLUMNS,
    _STEP_RETROPATH_COLUMNS,
    _SUMMARY_BASE_COLUMNS,
    _SUMMARY_RETROPATH_COLUMNS,
    _as_int,
    _atomic_write_json,
    _atomic_write_text,
    _build_solution_rows,
    _candidate_steps_by_id,
    _canonical_hash,
    _complete_ecs,
    _electron_row,
    _electron_summary,
    _field_order,
    _format_equation,
    _labels,
    _parse_stoichiometry,
    _read_csv,
    _render_csv,
    _safe_compound_name,
    _sha256_file,
    _split,
)
from src.pathway_analyze.retropath_gem_validation import (
    HYPOTHESIS_COLUMNS,
    PASSING_ROUTE_VALIDATION_STATUSES,
    RETROPATH_GEM_VALIDATION_SCHEMA,
    STOICHIOMETRY_HYPOTHESES_FILE_NAME,
    STOICHIOMETRY_TERMS_FILE_NAME,
    SUMMARY_COLUMNS,
    TERM_COLUMNS,
    VALIDATION_MANIFEST_FILE_NAME,
    VALIDATION_SUMMARY_FILE_NAME,
)
from src.pathway_analyze.retropath_mnxref import MnxrefIndex
from src.pathway_analyze.target_id import validate_target_compound_id


RETROPATH_MATERIALIZATION_SCHEMA = "retropath_solution_materialization.v1"
MATERIALIZATION_MANIFEST_FILE_NAME = "solution_materialization.json"

_SUMMARY_STATE_COLUMNS = (
    "materialization_id",
    "validation_status",
    "stoichiometry_status",
    "gem_status",
    "cofactor_mode",
    "cofactor_relaxed",
    "opened_generic_compound_ids",
    "validation_issue",
)
_STEP_STATE_COLUMNS = (
    "validation_status",
    "stoichiometry_status",
    "gem_status",
)

_NOT_VALIDATED_WARNING = (
    "GEM validation has not been run; this predicted route requires "
    "manual review before experimental use."
)
_VALIDATION_FAILED_WARNING = (
    "GEM validation did not pass; this predicted route remains "
    "writable but requires manual review before experimental use."
)
_RELAXED_VALIDATION_WARNING = (
    "Relaxed cofactor validation opened generic carrier sinks; feasibility "
    "is diagnostic and requires manual review before experimental use."
)


def _core_terms(
    left: Sequence[tuple[str, float]],
    right: Sequence[tuple[str, float]],
) -> list[dict[str, Any]]:
    return [
        {
            "side": side,
            "coefficient": coefficient,
            "role": "predicted_core",
            "compound_id": compound_id,
            "source_mnxm_id": "",
            "name": compound_id,
            "formula": "",
            "charge": "",
            "inchi": "",
            "inchikey": "",
            "smiles": "",
            "xrefs": [],
        }
        for side, values in (("left", left), ("right", right))
        for compound_id, coefficient in values
    ]


def _raw_solution_rows(
    *,
    solution_id: int,
    route: Mapping[str, str],
    steps: Sequence[Mapping[str, str]],
    target: str,
    materialization_id: str,
    kegg: KeggRestClient,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    candidate_rank = _as_int(route.get("candidate_rank"), "candidate_rank", minimum=1)
    candidate_id = str(route.get("candidate_id") or "").strip()
    rows: list[dict[str, Any]] = []
    electron_rows: list[dict[str, Any]] = []
    heterologous_reactions: list[str] = []
    heterologous_ecs: list[str] = []
    heterologous_kos: list[str] = []

    for step in steps:
        step_index = _as_int(step.get("step_index"), "step_index", minimum=1)
        step_id = str(step.get("step_id") or "").strip()
        step_source = str(step.get("step_source") or "").strip()
        left = _parse_stoichiometry(step.get("substrate_stoichiometry_json"))
        right = _parse_stoichiometry(step.get("product_stoichiometry_json"))
        products = tuple(compound_id for compound_id, _ in right)
        if step_source == "kegg_expansion":
            anchors = set(_split(step.get("sink_anchor_kegg_ids")))
            anchored_products = tuple(
                compound_id for compound_id in products if compound_id in anchors
            )
            if len(anchored_products) != 1:
                raise ValueError(
                    f"{step_id} must have exactly one anchored KEGG product"
                )
            primary_product = anchored_products[0]
        elif step_source == "retropath":
            if len(products) != 1:
                raise ValueError(
                    f"{step_id} must have exactly one RetroPath core product"
                )
            primary_product = products[0]
        else:
            raise ValueError(f"unsupported hybrid step source: {step_source}")
        precursors = tuple(compound_id for compound_id, _ in left)
        common = {
            "solution_id": solution_id,
            "step_index": step_index,
            "produced_compound_id": primary_product,
            "produced_compound_name": _safe_compound_name(kegg, primary_product),
            "direction": "left_to_right",
            "oxygen_required": "false",
            "thermo_direction": "unknown",
            "screening_rule_hits": "",
            "precursor_compound_ids": ";".join(precursors),
            "precursor_compound_labels": _labels(kegg, precursors),
            "ko_ids": "",
            "module_ids": "",
            "expansion_depth": str(step.get("expansion_depth") or "0"),
            "expansion_anchor_compounds": str(step.get("sink_anchor_kegg_ids") or ""),
            "prediction_review_required": "true",
            "structure_match_quality": str(route.get("structure_match_quality") or "exact"),
            "stereo_review_required": str(route.get("stereo_review_required") or "false"),
            "stereo_resolution_status": (
                "unresolved"
                if str(route.get("stereo_review_required") or "false").lower() == "true"
                else "not_required"
            ),
            "stereo_resolution_source": "",
            "depends_on_step_ids": str(step.get("depends_on_step_ids") or ""),
            "validation_status": "not_run",
            "stoichiometry_status": "core_only",
            "gem_status": "not_run",
        }
        if step_source == "kegg_expansion":
            reaction_ids = _split(step.get("reaction_option_ids"))
            if len(reaction_ids) != 1:
                raise ValueError(f"{step_id} must have one KEGG reaction")
            reaction_id = reaction_ids[0]
            record = kegg.try_get_reaction(reaction_id)
            if record is None:
                name = reaction_id
                comment = "KEGG expansion reaction; full annotation unavailable"
                equation = _format_equation(left, right)
                ecs = _complete_ecs(_split(step.get("source_ec_numbers")))
                ko_ids: Sequence[str] = ()
                module_ids: Sequence[str] = ()
                rhea_ids: Sequence[str] = ()
            else:
                name = record.name
                comment = record.comment
                equation = record.equation
                ecs = record.enzyme_ecs
                ko_ids = record.ko_ids
                module_ids = record.module_ids
                rhea_ids = record.rhea_ids
            status = (
                "endogenous"
                if str(step.get("status") or "").lower() == "endogenous"
                else "heterologous"
            )
            row = {
                **common,
                "status": status,
                "reaction_id": reaction_id,
                "reaction_name": name,
                "reaction_comment": comment,
                "equation": equation,
                "ko_ids": ";".join(ko_ids),
                "module_ids": ";".join(module_ids),
                "enzyme_ecs": ";".join(ecs),
                "locked_enzyme_ecs": ";".join(ecs),
                "ec_status": "complete" if ecs else "missing",
                "enzyme_search_eligible": str(status == "heterologous").lower(),
                "source_reaction_ids": ";".join(
                    dict.fromkeys((*_split(step.get("source_reaction_ids")), reaction_id))
                ),
                "resolution_action": "none",
                "resolution_evidence": "KEGG expansion witness",
                "step_source": "kegg_expansion",
                "rhea_ids": ";".join(rhea_ids),
                "auxiliary_requirements_json": "[]",
            }
            electron_ecs = ecs
            if status == "heterologous":
                heterologous_reactions.append(reaction_id)
                heterologous_ecs.extend(ecs)
                heterologous_kos.extend(ko_ids)
        elif step_source == "retropath":
            rule_ids = _split(step.get("rule_ids"))
            source_ecs = _complete_ecs(_split(step.get("source_ec_numbers")))
            source_reactions = _split(step.get("source_reaction_ids"))
            source_uniprot = _split(step.get("source_uniprot_ids"))
            core_smiles = str(step.get("reaction_smiles") or "").strip()
            terms = _core_terms(left, right)
            provenance = {
                "materialization_id": materialization_id,
                "candidate_rank": candidate_rank,
                "candidate_id": candidate_id,
                "step_id": step_id,
                "rule_ids": list(rule_ids),
                "source_reaction_ids": list(source_reactions),
                "source_ec_numbers": list(source_ecs),
                "source_uniprot_ids": list(source_uniprot),
                "core_reaction_smiles": core_smiles,
                "validation_status": "not_run",
                "review_status": "pending",
            }
            row = {
                **common,
                "status": "heterologous",
                "reaction_id": step_id,
                "reaction_name": "RetroPath core prediction",
                "reaction_comment": (
                    "Predicted RR02 core transformation; stoichiometry/GEM validation not run"
                ),
                "equation": _format_equation(left, right),
                "enzyme_ecs": "",
                "locked_enzyme_ecs": "",
                "ec_status": "template_only" if source_ecs else "missing",
                "enzyme_search_eligible": "true",
                "source_reaction_ids": ";".join(source_reactions),
                "resolution_action": "retropath_prediction",
                "resolution_evidence": "P5 core Reaction SMILES; manual review required",
                "step_source": "retropath",
                "rhea_ids": "",
                "auxiliary_requirements_json": "[]",
                "retropath_step_id": step_id,
                "retropath_hypothesis_id": "",
                "retropath_rule_id": ";".join(rule_ids),
                "source_mnxr_id": "",
                "source_ec_numbers": ";".join(source_ecs),
                "source_uniprot_ids": ";".join(source_uniprot),
                "source_rhea_ids": "",
                "exact_kegg_reaction_ids": "",
                "exact_rhea_ids": "",
                "formal_mapping_exact": "false",
                "reaction_signature_sha256": "",
                "full_reaction_smiles": "",
                "core_reaction_smiles": core_smiles,
                "stoichiometry_terms_json": json.dumps(
                    terms, ensure_ascii=False, sort_keys=True, separators=(",", ":")
                ),
                "prediction_provenance_json": json.dumps(
                    provenance, ensure_ascii=False, sort_keys=True, separators=(",", ":")
                ),
            }
            electron_ecs = source_ecs
            heterologous_reactions.append(step_id)
        else:
            raise ValueError(f"unsupported hybrid step source: {step_source}")

        electron = _electron_row(
            solution_id=solution_id,
            step_index=step_index,
            reaction_id=str(row["reaction_id"]),
            left=left,
            right=right,
            reaction_name=str(row["reaction_name"]),
            reaction_comment=str(row["reaction_comment"]),
            equation=str(row["equation"]),
            enzyme_ecs=electron_ecs,
        )
        row["auxiliary_requirements_json"] = electron["auxiliary_requirements_json"]
        rows.append(row)
        electron_rows.append(electron)

    electron_summary = _electron_summary(solution_id, electron_rows)
    anchors = _split(route.get("sink_kegg_ids"))
    summary = {
        "solution_id": solution_id,
        "target_compound_id": target,
        "target_compound_name": _safe_compound_name(kegg, target),
        "total_steps": len(rows),
        "heterologous_steps": sum(row["status"] == "heterologous" for row in rows),
        "heterologous_reaction_ids": ";".join(dict.fromkeys(heterologous_reactions)),
        "heterologous_ko_ids": ";".join(dict.fromkeys(heterologous_kos)),
        "heterologous_enzyme_ecs": ";".join(dict.fromkeys(heterologous_ecs)),
        "reaction_resolution_status": "predicted_core_stoichiometry",
        "normalization_event_count": 0,
        "normalization_events": "",
        "blocking_reaction_count": 0,
        "blocking_reaction_ids": "",
        "eligible_for_recommendation": "true",
        "reachable_anchor_compounds": ";".join(anchors),
        "reachable_anchor_labels": _labels(kegg, anchors),
        "frontier_anchor_compounds": ";".join(anchors),
        "expansion_bridge_steps": str(route.get("kegg_prefix_steps") or "0"),
        "max_expansion_depth": str(route.get("maximum_sink_depth") or "0"),
        "route_total_nadph_burden": 0,
        "route_total_sam_burden": 0,
        "route_total_coa_burden": 0,
        "oxygen_required_steps": sum(
            "C00007" in _split(row.get("precursor_compound_ids")) for row in rows
        ),
        "thermo_disfavored_steps": 0,
        **{key: value for key, value in electron_summary.items() if key != "solution_id"},
        "solution_source": "retropath",
        "retropath_candidate_rank": candidate_rank,
        "retropath_candidate_id": candidate_id,
        "retropath_combination_id": "",
        "prediction_review_required": "true",
        "structure_match_quality": str(route.get("structure_match_quality") or "exact"),
        "stereo_review_required": str(route.get("stereo_review_required") or "false"),
        "stereo_resolution_status": (
            "unresolved"
            if str(route.get("stereo_review_required") or "false").lower() == "true"
            else "not_required"
        ),
        "stereo_resolution_source": "",
        "promotion_id": "",
        "materialization_id": materialization_id,
        "validation_status": "not_run",
        "stoichiometry_status": "core_only",
        "gem_status": "not_run",
        "cofactor_mode": "not_run",
        "cofactor_relaxed": "false",
        "opened_generic_compound_ids": "",
        "validation_issue": _NOT_VALIDATED_WARNING,
        "combination_truncated": "false",
        "upstream_enumeration_truncated": str(
            route.get("upstream_enumeration_truncated") or "false"
        ),
        "candidate_top_k_truncated": str(
            route.get("candidate_top_k_truncated") or "false"
        ),
    }
    return summary, rows, electron_summary, electron_rows


def materialize_retropath_candidate_solutions(config: Any) -> dict[str, Any]:
    """Replace the RetroPath solution slice with every current P5 Top-K candidate."""

    target = validate_target_compound_id(config.target_name)
    depth = _as_int(getattr(config, "depth", 0), "depth")
    gap_dir = gap_depth_output_dir(Path(config.gap_output_path), depth).resolve()
    retropath_dir = gap_dir / "retropath"
    routes_path = retropath_dir / CANDIDATE_ROUTES_FILE_NAME
    steps_path = retropath_dir / CANDIDATE_STEPS_FILE_NAME
    rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    if not rules_path.is_file():
        raise ValueError(f"RR02 rules file is missing: {rules_path}")
    routes = _read_csv(routes_path, CANDIDATE_ROUTE_COLUMNS)
    steps_by_candidate = _candidate_steps_by_id(
        _read_csv(steps_path, CANDIDATE_STEP_COLUMNS)
    )
    routes.sort(key=lambda row: _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1))
    seed = {
        "schema_version": RETROPATH_MATERIALIZATION_SCHEMA,
        "target_compound": target,
        "expansion_depth": depth,
        "candidate_routes_sha256": _sha256_file(routes_path),
        "candidate_steps_sha256": _sha256_file(steps_path),
        "rr02_sha256": _sha256_file(rules_path),
        "candidates": [
            {
                "candidate_rank": _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1),
                "candidate_id": str(row.get("candidate_id") or "").strip(),
            }
            for row in routes
        ],
    }
    materialization_id = "RP2MATERIALIZE:" + _canonical_hash(seed)

    solution_path = gap_dir / "solutions.csv"
    step_path = gap_dir / "all_solution_steps.csv"
    electron_summary_path = gap_dir / "solution_electron_summary.csv"
    electron_step_path = gap_dir / "route_electron_requirements.csv"
    existing_summaries = _read_csv(solution_path)
    existing_steps = _read_csv(step_path)
    existing_electron_summaries = _read_csv(electron_summary_path)
    existing_electron_steps = _read_csv(electron_step_path)
    kegg_summaries = [
        row for row in existing_summaries
        if str(row.get("solution_source") or "kegg").strip().lower() != "retropath"
    ]
    kegg_ids = {
        _as_int(row.get("solution_id"), "solution_id", minimum=1)
        for row in kegg_summaries
    }
    next_solution_id = max(kegg_ids, default=0) + 1
    kegg_steps = [
        row for row in existing_steps
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) in kegg_ids
    ]
    kegg_electron_summaries = [
        row for row in existing_electron_summaries
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) in kegg_ids
    ]
    kegg_electron_steps = [
        row for row in existing_electron_steps
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) in kegg_ids
    ]
    for row in kegg_summaries:
        row.setdefault("solution_source", "kegg")

    kegg = KeggRestClient(Path(config.cache_dir).resolve() / "kegg")
    rp_summaries: list[dict[str, Any]] = []
    rp_steps: list[dict[str, Any]] = []
    rp_electron_summaries: list[dict[str, Any]] = []
    rp_electron_steps: list[dict[str, Any]] = []
    mappings: list[dict[str, Any]] = []
    for offset, route in enumerate(routes):
        candidate_id = str(route.get("candidate_id") or "").strip()
        steps = steps_by_candidate.get(candidate_id)
        if not steps:
            raise ValueError(f"P5 candidate has no steps: {candidate_id}")
        solution_id = next_solution_id + offset
        summary, step_rows, electron_summary, electron_rows = _raw_solution_rows(
            solution_id=solution_id,
            route=route,
            steps=steps,
            target=target,
            materialization_id=materialization_id,
            kegg=kegg,
        )
        rp_summaries.append(summary)
        rp_steps.extend(step_rows)
        rp_electron_summaries.append(electron_summary)
        rp_electron_steps.extend(electron_rows)
        mappings.append(
            {
                "solution_id": solution_id,
                "candidate_rank": _as_int(route.get("candidate_rank"), "candidate_rank", minimum=1),
                "candidate_id": candidate_id,
                "validation": {
                    "status": "not_run",
                    "stoichiometry_status": "core_only",
                    "gem_status": "not_run",
                    "cofactor_mode": "not_run",
                    "cofactor_relaxed": False,
                    "opened_generic_compound_ids": [],
                    "combination_id": "",
                    "stoichiometry_hypothesis_ids": [],
                    "issues": [_NOT_VALIDATED_WARNING],
                },
            }
        )

    summary_rows = [*kegg_summaries, *rp_summaries]
    step_rows = [*kegg_steps, *rp_steps]
    electron_summary_rows = [*kegg_electron_summaries, *rp_electron_summaries]
    electron_step_rows = [*kegg_electron_steps, *rp_electron_steps]
    texts = {
        "solutions.csv": _render_csv(
            _field_order(
                summary_rows,
                _SUMMARY_BASE_COLUMNS,
                (*_SUMMARY_RETROPATH_COLUMNS, *_SUMMARY_STATE_COLUMNS),
            ),
            summary_rows,
        ),
        "all_solution_steps.csv": _render_csv(
            _field_order(
                step_rows,
                _STEP_BASE_COLUMNS,
                (*_STEP_RETROPATH_COLUMNS, *_STEP_STATE_COLUMNS),
            ),
            step_rows,
        ),
        "solution_electron_summary.csv": _render_csv(
            _field_order(electron_summary_rows, _ELECTRON_SUMMARY_COLUMNS, ()),
            electron_summary_rows,
        ),
        "route_electron_requirements.csv": _render_csv(
            _field_order(electron_step_rows, _ELECTRON_STEP_COLUMNS, ()),
            electron_step_rows,
        ),
    }
    output_hashes = {
        name: _atomic_write_text(gap_dir / name, text) for name, text in texts.items()
    }
    manifest_path = retropath_dir / MATERIALIZATION_MANIFEST_FILE_NAME
    manifest = {
        **seed,
        "materialization_id": materialization_id,
        "created_at": datetime.now(UTC).isoformat(),
        "solution_count": len(mappings),
        "solution_mappings": mappings,
        "validation_overlay": None,
        "inputs": {
            "candidate_routes": {"path": str(routes_path.resolve()), "sha256": seed["candidate_routes_sha256"]},
            "candidate_steps": {"path": str(steps_path.resolve()), "sha256": seed["candidate_steps_sha256"]},
            "rr02": {"path": str(rules_path), "sha256": seed["rr02_sha256"]},
        },
        "artifacts": {
            name: {"path": str((gap_dir / name).resolve()), "sha256": digest}
            for name, digest in output_hashes.items()
        },
    }
    manifest_sha256 = _atomic_write_json(manifest_path, manifest)
    return {
        "ok": True,
        "schema_version": RETROPATH_MATERIALIZATION_SCHEMA,
        "materialization_id": materialization_id,
        "materialization_manifest": str(manifest_path.resolve()),
        "materialization_manifest_sha256": manifest_sha256,
        "solution_count": len(mappings),
        "solution_ids": [item["solution_id"] for item in mappings],
        "solution_mappings": mappings,
        "output_hashes": output_hashes,
    }


def _load_materialization(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("RetroPath solution materialization manifest is invalid") from exc
    if not isinstance(payload, dict):
        raise ValueError("RetroPath solution materialization manifest is invalid")
    return payload


def verify_retropath_solution_materialization(
    *,
    gap_dir: str | Path,
    target_compound: str,
    expansion_depth: int,
    solution_id: int | None = None,
    replacement_validation_path: str | Path | None = None,
) -> dict[str, Any]:
    resolved_gap = Path(gap_dir).expanduser().resolve()
    path = resolved_gap / "retropath" / MATERIALIZATION_MANIFEST_FILE_NAME
    payload = _load_materialization(path)
    if (
        payload.get("schema_version") != RETROPATH_MATERIALIZATION_SCHEMA
        or payload.get("target_compound") != target_compound
        or _as_int(payload.get("expansion_depth"), "expansion_depth") != expansion_depth
    ):
        raise ValueError("RetroPath solution materialization identity mismatch")
    seed = {
        "schema_version": payload.get("schema_version"),
        "target_compound": payload.get("target_compound"),
        "expansion_depth": payload.get("expansion_depth"),
        "candidate_routes_sha256": payload.get("candidate_routes_sha256"),
        "candidate_steps_sha256": payload.get("candidate_steps_sha256"),
        "rr02_sha256": payload.get("rr02_sha256"),
        "candidates": payload.get("candidates"),
    }
    expected_id = "RP2MATERIALIZE:" + _canonical_hash(seed)
    if payload.get("materialization_id") != expected_id:
        raise ValueError("RetroPath solution materialization ID mismatch")
    for input_key, display_name in (
        ("candidate_routes", CANDIDATE_ROUTES_FILE_NAME),
        ("candidate_steps", CANDIDATE_STEPS_FILE_NAME),
        ("rr02", "RR02 rules"),
    ):
        record = payload.get("inputs", {}).get(input_key)
        if not isinstance(record, Mapping):
            raise ValueError(f"RetroPath materialization input missing: {display_name}")
        candidate = Path(str(record.get("path") or "")).expanduser().resolve()
        if not candidate.is_file() or record.get("sha256") != _sha256_file(candidate):
            raise ValueError(f"RetroPath materialization input stale: {display_name}")
        if record.get("sha256") != payload.get(f"{input_key}_sha256"):
            raise ValueError(f"RetroPath materialization input binding mismatch: {display_name}")
    for name in (
        "solutions.csv",
        "all_solution_steps.csv",
        "solution_electron_summary.csv",
        "route_electron_requirements.csv",
    ):
        record = payload.get("artifacts", {}).get(name)
        expected = (resolved_gap / name).resolve()
        if (
            not isinstance(record, Mapping)
            or Path(str(record.get("path") or "")).expanduser().resolve() != expected
            or not expected.is_file()
            or record.get("sha256") != _sha256_file(expected)
        ):
            raise ValueError(f"RetroPath materialization artifact mismatch: {name}")
    overlay = payload.get("validation_overlay")
    replacement_path = (
        Path(replacement_validation_path).expanduser().resolve()
        if replacement_validation_path is not None
        else None
    )
    if overlay is not None:
        if not isinstance(overlay, Mapping):
            raise ValueError("RetroPath validation overlay is invalid")
        validation_path = Path(str(overlay.get("path") or "")).expanduser().resolve()
        if (
            validation_path != replacement_path
            and (
                not validation_path.is_file()
                or overlay.get("sha256") != _sha256_file(validation_path)
            )
        ):
            raise ValueError("RetroPath validation overlay is stale")
    mappings = payload.get("solution_mappings", [])
    candidates = payload.get("candidates", [])
    if (
        not isinstance(mappings, list)
        or not isinstance(candidates, list)
        or _as_int(payload.get("solution_count"), "solution_count") != len(mappings)
        or len(mappings) != len(candidates)
    ):
        raise ValueError("RetroPath solution mapping count mismatch")
    mapping_ids: set[int] = set()
    for mapping, candidate_record in zip(mappings, candidates):
        if not isinstance(mapping, Mapping):
            raise ValueError("RetroPath solution mapping is invalid")
        if not isinstance(candidate_record, Mapping):
            raise ValueError("RetroPath candidate binding is invalid")
        solution_number = _as_int(
            mapping.get("solution_id"), "solution_id", minimum=1
        )
        if solution_number in mapping_ids:
            raise ValueError("RetroPath solution mapping contains duplicate IDs")
        mapping_ids.add(solution_number)
        if (
            mapping.get("candidate_id") != candidate_record.get("candidate_id")
            or _as_int(mapping.get("candidate_rank"), "candidate_rank", minimum=1)
            != _as_int(
                candidate_record.get("candidate_rank"),
                "candidate_rank",
                minimum=1,
            )
        ):
            raise ValueError("RetroPath solution candidate binding mismatch")
        state = mapping.get("validation")
        if not isinstance(state, Mapping):
            raise ValueError("RetroPath solution validation state is invalid")
        if state.get("status") not in {"passed", "failed"}:
            continue
        validation_path = Path(
            str(state.get("manifest_path") or "")
        ).expanduser().resolve()
        if validation_path != replacement_path and (
            not validation_path.is_file()
            or state.get("manifest_sha256") != _sha256_file(validation_path)
        ):
            raise ValueError("RetroPath stored validation binding is stale")
    summary_rows = _read_csv(resolved_gap / "solutions.csv")
    for mapping in mappings:
        matches = [
            row
            for row in summary_rows
            if _as_int(row.get("solution_id"), "solution_id", minimum=1)
            == _as_int(mapping.get("solution_id"), "solution_id", minimum=1)
        ]
        if (
            len(matches) != 1
            or str(matches[0].get("solution_source") or "").strip().lower()
            != "retropath"
            or matches[0].get("retropath_candidate_id")
            != mapping.get("candidate_id")
        ):
            raise ValueError("RetroPath materialized solution row binding mismatch")
    if solution_id is not None:
        matches = [
            item for item in payload.get("solution_mappings", [])
            if isinstance(item, Mapping)
            and _as_int(item.get("solution_id"), "solution_id", minimum=1) == solution_id
        ]
        if len(matches) != 1:
            raise ValueError(f"RetroPath solution {solution_id} is not materialized")
    return payload


def apply_retropath_validation_overlay(
    config: Any,
    *,
    validation_manifest_path: str | Path,
) -> dict[str, Any]:
    """Attach deterministic P8 results to existing solution IDs without fan-out."""

    target = validate_target_compound_id(config.target_name)
    depth = _as_int(getattr(config, "depth", 0), "depth")
    gap_dir = gap_depth_output_dir(Path(config.gap_output_path), depth).resolve()
    retropath_dir = gap_dir / "retropath"
    materialization_path = retropath_dir / MATERIALIZATION_MANIFEST_FILE_NAME
    p8_path = Path(validation_manifest_path).expanduser().resolve()
    if p8_path != (retropath_dir / "gem_validation" / VALIDATION_MANIFEST_FILE_NAME).resolve():
        raise ValueError("unexpected P8 validation manifest path")
    materialization = verify_retropath_solution_materialization(
        gap_dir=gap_dir,
        target_compound=target,
        expansion_depth=depth,
        replacement_validation_path=p8_path,
    )
    try:
        p8_manifest = json.loads(p8_path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("P8 validation manifest is invalid") from exc
    if (
        not isinstance(p8_manifest, Mapping)
        or p8_manifest.get("schema_version") != RETROPATH_GEM_VALIDATION_SCHEMA
        or p8_manifest.get("target_compound") != target
        or _as_int(p8_manifest.get("expansion_depth"), "expansion_depth") != depth
    ):
        raise ValueError("P8 validation manifest identity mismatch")
    for name, record in p8_manifest.get("artifacts", {}).items():
        if not isinstance(record, Mapping):
            raise ValueError(f"P8 artifact record invalid: {name}")
        artifact = Path(str(record.get("path") or "")).expanduser().resolve()
        if not artifact.is_file() or record.get("sha256") != _sha256_file(artifact):
            raise ValueError(f"P8 artifact mismatch: {name}")
    p8_inputs = p8_manifest.get("inputs", {})
    if (
        p8_inputs.get("candidate_routes_sha256")
        != materialization["candidate_routes_sha256"]
        or p8_inputs.get("candidate_steps_sha256")
        != materialization["candidate_steps_sha256"]
        or p8_inputs.get("rr02_sha256") != materialization["rr02_sha256"]
    ):
        raise ValueError(
            "P8 validation does not bind the materialized P5 candidates/RR02 rules"
        )

    p8_dir = p8_path.parent
    summary_rows = _read_csv(p8_dir / VALIDATION_SUMMARY_FILE_NAME, SUMMARY_COLUMNS)
    selected_ranks = {
        _as_int(value, "selected_candidate_rank", minimum=1)
        for value in p8_manifest.get("selected_candidate_ranks", [])
    }
    materialized_ranks = {
        _as_int(item.get("candidate_rank"), "candidate_rank", minimum=1)
        for item in materialization.get("solution_mappings", [])
        if isinstance(item, Mapping)
    }
    if not selected_ranks or not selected_ranks.issubset(materialized_ranks):
        raise ValueError("P8 selected candidate ranks do not match materialized solutions")
    hypotheses = {
        str(row.get("hypothesis_id") or ""): row
        for row in _read_csv(p8_dir / STOICHIOMETRY_HYPOTHESES_FILE_NAME, HYPOTHESIS_COLUMNS)
    }
    terms_by_hypothesis: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read_csv(p8_dir / STOICHIOMETRY_TERMS_FILE_NAME, TERM_COLUMNS):
        terms_by_hypothesis[str(row.get("hypothesis_id") or "")].append(row)
    routes = {
        str(row.get("candidate_id") or ""): row
        for row in _read_csv(retropath_dir / CANDIDATE_ROUTES_FILE_NAME, CANDIDATE_ROUTE_COLUMNS)
    }
    steps_by_candidate = _candidate_steps_by_id(
        _read_csv(retropath_dir / CANDIDATE_STEPS_FILE_NAME, CANDIDATE_STEP_COLUMNS)
    )
    rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    from src.pathway_analyze.retropath_promotion import _load_rules

    rules = _load_rules(rules_path)
    kegg = KeggRestClient(Path(config.cache_dir).resolve() / "kegg")
    mnxref_dir = Path(config.data_dir) / "retropath" / "mnxref" / "3.0"
    existing_solution_rows = _read_csv(gap_dir / "solutions.csv")
    existing_step_rows = _read_csv(gap_dir / "all_solution_steps.csv")
    existing_electron_summaries = _read_csv(
        gap_dir / "solution_electron_summary.csv"
    )
    existing_electron_steps = _read_csv(
        gap_dir / "route_electron_requirements.csv"
    )
    retropath_solution_ids = {
        _as_int(item.get("solution_id"), "solution_id", minimum=1)
        for item in materialization.get("solution_mappings", [])
        if isinstance(item, Mapping)
    }
    solution_rows = [
        row for row in existing_solution_rows
        if _as_int(row.get("solution_id"), "solution_id", minimum=1)
        not in retropath_solution_ids
    ]
    step_rows = [
        row for row in existing_step_rows
        if _as_int(row.get("solution_id"), "solution_id", minimum=1)
        not in retropath_solution_ids
    ]
    electron_summaries = [
        row for row in existing_electron_summaries
        if _as_int(row.get("solution_id"), "solution_id", minimum=1)
        not in retropath_solution_ids
    ]
    electron_steps = [
        row for row in existing_electron_steps
        if _as_int(row.get("solution_id"), "solution_id", minimum=1)
        not in retropath_solution_ids
    ]
    for original in materialization.get("solution_mappings", []):
        solution_id = _as_int(
            original.get("solution_id"), "solution_id", minimum=1
        )
        candidate_id = str(original.get("candidate_id") or "")
        route = routes.get(candidate_id)
        candidate_steps = steps_by_candidate.get(candidate_id)
        if route is None or not candidate_steps:
            raise ValueError(f"P5 candidate inputs missing for {candidate_id}")
        raw = _raw_solution_rows(
            solution_id=solution_id,
            route=route,
            steps=candidate_steps,
            target=target,
            materialization_id=str(materialization.get("materialization_id") or ""),
            kegg=kegg,
        )
        solution_rows.append(raw[0])
        step_rows.extend(raw[1])
        electron_summaries.append(raw[2])
        electron_steps.extend(raw[3])
    mapping_updates: list[dict[str, Any]] = []
    selected_by_solution: dict[int, tuple[dict[str, str], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]] = {}

    with MnxrefIndex(mnxref_dir, rules_path) as mnxref:
        for original in materialization.get("solution_mappings", []):
            solution_id = _as_int(
                original.get("solution_id"), "solution_id", minimum=1
            )
            candidate_id = str(original.get("candidate_id") or "")
            candidate_rank = _as_int(
                original.get("candidate_rank"), "candidate_rank", minimum=1
            )
            mapping = {
                "solution_id": solution_id,
                "candidate_rank": candidate_rank,
                "candidate_id": candidate_id,
                "validation": {
                    "status": "not_run",
                    "stoichiometry_status": "core_only",
                    "gem_status": "not_run",
                    "cofactor_mode": "not_run",
                    "cofactor_relaxed": False,
                    "opened_generic_compound_ids": [],
                    "combination_id": "",
                    "stoichiometry_hypothesis_ids": [],
                    "issues": [_NOT_VALIDATED_WARNING],
                },
            }
            if candidate_rank not in selected_ranks:
                mapping_updates.append(mapping)
                continue
            candidate_rows = [
                row for row in summary_rows
                if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1) == candidate_rank
                and str(row.get("candidate_id") or "") == candidate_id
            ]
            passing = [
                row for row in candidate_rows
                if row.get("validation_status")
                in PASSING_ROUTE_VALIDATION_STATUSES
                and str(row.get("combination_id") or "").strip()
            ]
            validation_mode = str(
                (candidate_rows[0].get("cofactor_mode") if candidate_rows else "")
                or p8_manifest.get("parameters", {}).get("cofactor_mode")
                or "strict_l1"
            ).strip().lower()
            cofactor_relaxed = validation_mode == "relaxed"
            opened_generic_compounds = sorted({
                compound_id
                for row in candidate_rows
                for compound_id in _split(
                    row.get("opened_generic_compound_ids")
                )
            })
            issues = tuple(dict.fromkeys(
                str(row.get("issues") or "").strip()
                for row in candidate_rows
                if str(row.get("issues") or "").strip()
            ))
            if cofactor_relaxed:
                issues = tuple(dict.fromkeys((*issues, _RELAXED_VALIDATION_WARNING)))
            if passing:
                selected = passing[0]
                status = "passed"
                stoichiometry_status = "completed"
                gem_status = "passed"
                hypothesis_ids = list(_split(selected.get("stoichiometry_hypothesis_ids")))
                route = routes.get(candidate_id)
                steps = steps_by_candidate.get(candidate_id)
                if route is None or not steps:
                    raise ValueError(f"P5 candidate inputs missing for {candidate_id}")
                validated = _build_solution_rows(
                    solution_id=solution_id,
                    passing=selected,
                    route=route,
                    steps=steps,
                    hypotheses=hypotheses,
                    terms_by_hypothesis=terms_by_hypothesis,
                    rules=rules,
                    mnxref=mnxref,
                    kegg=kegg,
                    target_compound=target,
                    promotion_id=str(materialization.get("materialization_id") or ""),
                )
                selected_by_solution[solution_id] = (selected, validated[1], validated[2], validated[3])
                combination_id = str(selected.get("combination_id") or "")
            else:
                status = "failed"
                gem_status = "failed"
                stoichiometry_status = "completed" if any(
                    str(row.get("combination_id") or "").strip() for row in candidate_rows
                ) else "incomplete"
                hypothesis_ids = []
                combination_id = ""
                issues = tuple(dict.fromkeys((*issues, _VALIDATION_FAILED_WARNING)))
            mapping["validation"] = {
                "status": status,
                "stoichiometry_status": stoichiometry_status,
                "gem_status": gem_status,
                "cofactor_mode": validation_mode,
                "cofactor_relaxed": cofactor_relaxed,
                "opened_generic_compound_ids": opened_generic_compounds,
                "combination_id": combination_id,
                "stoichiometry_hypothesis_ids": hypothesis_ids,
                "issues": list(issues),
                "manifest_path": str(p8_path),
                "manifest_sha256": _sha256_file(p8_path),
            }
            mapping_updates.append(mapping)

    mapping_by_solution = {
        _as_int(item["solution_id"], "solution_id", minimum=1): item
        for item in mapping_updates
    }
    for row in solution_rows:
        solution_id = _as_int(row.get("solution_id"), "solution_id", minimum=1)
        mapping = mapping_by_solution.get(solution_id)
        if mapping is None:
            continue
        validation = mapping["validation"]
        row["validation_status"] = validation["status"]
        row["stoichiometry_status"] = validation["stoichiometry_status"]
        row["gem_status"] = validation["gem_status"]
        row["cofactor_mode"] = validation["cofactor_mode"]
        row["cofactor_relaxed"] = str(
            validation["cofactor_relaxed"]
        ).lower()
        row["opened_generic_compound_ids"] = ";".join(
            validation["opened_generic_compound_ids"]
        )
        row["validation_issue"] = "; ".join(validation["issues"])
        row["retropath_combination_id"] = validation["combination_id"]
        if validation["status"] == "passed":
            row["reaction_resolution_status"] = (
                "predicted_relaxed_stoichiometry"
                if validation["cofactor_relaxed"]
                else "predicted_strict_stoichiometry"
            )
            validated_electron_summary = selected_by_solution[solution_id][2]
            for field, value in validated_electron_summary.items():
                if field not in {"solution_id", "solution_source"}:
                    row[field] = value

    validation_step_fields = {
        "retropath_hypothesis_id",
        "source_mnxr_id",
        "source_ec_numbers",
        "source_uniprot_ids",
        "source_rhea_ids",
        "exact_kegg_reaction_ids",
        "exact_rhea_ids",
        "formal_mapping_exact",
        "reaction_signature_sha256",
        "full_reaction_smiles",
        "stoichiometry_terms_json",
        "prediction_provenance_json",
        "enzyme_ecs",
        "locked_enzyme_ecs",
        "ec_status",
        "rhea_ids",
        "reaction_name",
        "reaction_comment",
        "equation",
        "resolution_evidence",
        "auxiliary_requirements_json",
    }
    for row in step_rows:
        solution_id = _as_int(row.get("solution_id"), "solution_id", minimum=1)
        mapping = mapping_by_solution.get(solution_id)
        if mapping is None:
            continue
        validation = mapping["validation"]
        row["validation_status"] = validation["status"]
        row["stoichiometry_status"] = validation["stoichiometry_status"]
        row["gem_status"] = validation["gem_status"]
    for solution_id, (_, validated_steps, _, _) in selected_by_solution.items():
        validated_by_step = {
            str(row.get("retropath_step_id") or ""): row
            for row in validated_steps
            if str(row.get("step_source") or "") == "retropath"
        }
        for row in step_rows:
            if _as_int(row.get("solution_id"), "solution_id", minimum=1) != solution_id:
                continue
            validated = validated_by_step.get(str(row.get("retropath_step_id") or ""))
            if validated is not None:
                for field in validation_step_fields:
                    row[field] = validated.get(field, "")

    validated_solution_ids = set(selected_by_solution)
    electron_summaries = [
        row for row in electron_summaries
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) not in validated_solution_ids
    ]
    electron_steps = [
        row for row in electron_steps
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) not in validated_solution_ids
    ]
    for _, (_, _, summary, rows) in selected_by_solution.items():
        electron_summaries.append(summary)
        electron_steps.extend(rows)

    texts = {
        "solutions.csv": _render_csv(
            _field_order(solution_rows, _SUMMARY_BASE_COLUMNS, (*_SUMMARY_RETROPATH_COLUMNS, *_SUMMARY_STATE_COLUMNS)),
            solution_rows,
        ),
        "all_solution_steps.csv": _render_csv(
            _field_order(step_rows, _STEP_BASE_COLUMNS, (*_STEP_RETROPATH_COLUMNS, *_STEP_STATE_COLUMNS)),
            step_rows,
        ),
        "solution_electron_summary.csv": _render_csv(
            _field_order(electron_summaries, _ELECTRON_SUMMARY_COLUMNS, ()), electron_summaries
        ),
        "route_electron_requirements.csv": _render_csv(
            _field_order(electron_steps, _ELECTRON_STEP_COLUMNS, ()), electron_steps
        ),
    }
    output_hashes = {
        name: _atomic_write_text(gap_dir / name, text) for name, text in texts.items()
    }
    materialization["updated_at"] = datetime.now(UTC).isoformat()
    materialization["solution_mappings"] = mapping_updates
    materialization["validation_overlay"] = {
        "path": str(p8_path),
        "sha256": _sha256_file(p8_path),
        "schema_version": RETROPATH_GEM_VALIDATION_SCHEMA,
    }
    materialization["artifacts"] = {
        name: {"path": str((gap_dir / name).resolve()), "sha256": digest}
        for name, digest in output_hashes.items()
    }
    manifest_sha256 = _atomic_write_json(materialization_path, materialization)
    passed_ids = [
        item["solution_id"] for item in mapping_updates
        if item["validation"]["status"] == "passed"
    ]
    return {
        "ok": True,
        "materialization_manifest": str(materialization_path),
        "materialization_manifest_sha256": manifest_sha256,
        "promotion_manifest": str(materialization_path),
        "promotion_manifest_sha256": manifest_sha256,
        "solution_ids": [item["solution_id"] for item in mapping_updates],
        "validated_solution_ids": passed_ids,
        "formal_solution_ids": passed_ids,
        "solution_mappings": mapping_updates,
        "output_hashes": output_hashes,
    }


__all__ = [
    "MATERIALIZATION_MANIFEST_FILE_NAME",
    "RETROPATH_MATERIALIZATION_SCHEMA",
    "apply_retropath_validation_overlay",
    "materialize_retropath_candidate_solutions",
    "verify_retropath_solution_materialization",
]
