"""Materialize strictly validated RetroPath combinations as formal solutions.

The raw RetroPath/P5/P8 artifacts remain the evidence source.  This module
projects every ``PASS_STRICT_ROUTE_FLUX`` combination into the same four CSV
tables consumed by the regular solution/info/manifest workflow.  Existing
KEGG rows keep their identifiers; the RetroPath slice is replaced atomically
enough to fail closed through a commit-marker manifest written last.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import tempfile
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from src.pathway_analyze.kegg_gap_analyze import (
    KeggRestClient,
    format_electron_carrier_net_changes,
    gap_depth_output_dir,
    infer_electron_requirement,
    summarize_solution_electron_requirements,
)
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEPS_FILE_NAME,
)
from src.pathway_analyze.retropath_gem_validation import (
    HYPOTHESIS_COLUMNS,
    STOICHIOMETRY_HYPOTHESES_FILE_NAME,
    STOICHIOMETRY_TERMS_FILE_NAME,
    SUMMARY_COLUMNS,
    TERM_COLUMNS,
    VALIDATION_MANIFEST_FILE_NAME,
    VALIDATION_SUMMARY_FILE_NAME,
)
from src.pathway_analyze.retropath_mnxref import MnxrefIndex
from src.pathway_analyze.target_id import validate_target_compound_id


RETROPATH_PROMOTION_SCHEMA = "retropath_solution_promotion.v1"
PROMOTION_MANIFEST_FILE_NAME = "formal_solution_promotion.json"

_KEGG_REACTION = re.compile(r"^R\d{5}$", re.IGNORECASE)
_KEGG_XREF = re.compile(r"^kegg:(R\d{5})$", re.IGNORECASE)
_RHEA_XREF = re.compile(r"^rhea:(\d+)$", re.IGNORECASE)
_COMPLETE_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

_SUMMARY_BASE_COLUMNS = (
    "solution_id",
    "target_compound_id",
    "target_compound_name",
    "total_steps",
    "heterologous_steps",
    "heterologous_reaction_ids",
    "heterologous_ko_ids",
    "heterologous_enzyme_ecs",
    "reaction_resolution_status",
    "normalization_event_count",
    "normalization_events",
    "blocking_reaction_count",
    "blocking_reaction_ids",
    "eligible_for_recommendation",
    "reachable_anchor_compounds",
    "reachable_anchor_labels",
    "frontier_anchor_compounds",
    "expansion_bridge_steps",
    "max_expansion_depth",
    "route_total_nadph_burden",
    "route_total_sam_burden",
    "route_total_coa_burden",
    "oxygen_required_steps",
    "thermo_disfavored_steps",
    "max_electron_risk_level",
    "max_electron_risk_score",
    "electron_system_status",
    "electron_balance_status",
    "requires_external_electron_regeneration",
    "requires_carrier_compatibility_check",
)
_SUMMARY_RETROPATH_COLUMNS = (
    "solution_source",
    "retropath_candidate_rank",
    "retropath_candidate_id",
    "retropath_combination_id",
    "prediction_review_required",
    "promotion_id",
    "combination_truncated",
    "upstream_enumeration_truncated",
    "candidate_top_k_truncated",
)
_STEP_BASE_COLUMNS = (
    "solution_id",
    "step_index",
    "status",
    "produced_compound_id",
    "produced_compound_name",
    "reaction_id",
    "reaction_name",
    "reaction_comment",
    "equation",
    "direction",
    "oxygen_required",
    "thermo_direction",
    "screening_rule_hits",
    "precursor_compound_ids",
    "precursor_compound_labels",
    "ko_ids",
    "module_ids",
    "enzyme_ecs",
    "locked_enzyme_ecs",
    "ec_status",
    "enzyme_search_eligible",
    "source_reaction_ids",
    "resolution_action",
    "resolution_evidence",
    "step_source",
    "expansion_depth",
    "expansion_anchor_compounds",
    "rhea_ids",
    "auxiliary_requirements_json",
)
_STEP_RETROPATH_COLUMNS = (
    "retropath_step_id",
    "retropath_hypothesis_id",
    "retropath_rule_id",
    "source_mnxr_id",
    "source_ec_numbers",
    "source_uniprot_ids",
    "source_rhea_ids",
    "exact_kegg_reaction_ids",
    "exact_rhea_ids",
    "formal_mapping_exact",
    "reaction_signature_sha256",
    "full_reaction_smiles",
    "core_reaction_smiles",
    "stoichiometry_terms_json",
    "prediction_provenance_json",
    "prediction_review_required",
    "depends_on_step_ids",
)
_ELECTRON_SUMMARY_COLUMNS = (
    "solution_id",
    "max_electron_risk_level",
    "max_electron_risk_score",
    "electron_system_status",
    "electron_balance_status",
    "requires_external_electron_regeneration",
    "requires_carrier_compatibility_check",
    "electron_carrier_ids",
    "electron_requirement_classes",
    "electron_carrier_net_changes",
    "balanced_electron_carrier_pairs",
    "unbalanced_electron_carrier_pairs",
    "unresolved_electron_carrier_ids",
    "annotation_only_electron_requirements",
    "required_auxiliary_roles",
    "auxiliary_requirement_status",
    "solution_source",
)
_ELECTRON_STEP_COLUMNS = (
    "solution_id",
    "step_index",
    "reaction_id",
    "electron_carrier_ids",
    "electron_requirement_classes",
    "electron_risk_level",
    "electron_risk_score",
    "electron_risk_evidence",
    "electron_carrier_net_changes",
    "required_auxiliary_roles",
    "auxiliary_requirement_status",
    "auxiliary_requirements_json",
    "solution_source",
)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _split(value: Any) -> tuple[str, ...]:
    if isinstance(value, (list, tuple, set)):
        raw = value
    else:
        raw = str(value or "").replace("|", ";").split(";")
    return tuple(dict.fromkeys(
        str(item).strip() for item in raw if str(item).strip()
    ))


def _as_int(value: Any, field: str, *, minimum: int = 0) -> int:
    try:
        result = int(float(str(value).strip()))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field}: {value!r}") from exc
    if result < minimum:
        raise ValueError(f"invalid {field}: {value!r}")
    return result


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def _read_csv(path: Path, required: Sequence[str] = ()) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        missing = sorted(set(required) - set(fields))
        if missing:
            raise ValueError(f"{path.name} missing columns: {missing}")
        return [dict(row) for row in reader]


def _field_order(
    existing: Sequence[Mapping[str, Any]],
    defaults: Sequence[str],
    additions: Sequence[str],
) -> tuple[str, ...]:
    result: list[str] = []
    for name in (
        *(existing[0].keys() if existing else ()),
        *defaults,
        *additions,
        *(key for row in existing for key in row),
    ):
        if name and name not in result:
            result.append(name)
    return tuple(result)


def _render_csv(columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=tuple(columns),
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({name: row.get(name, "") for name in columns})
    return buffer.getvalue()


def _atomic_write_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
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
    return _sha256_bytes(text.encode("utf-8"))


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> str:
    text = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        indent=2,
        allow_nan=False,
    ) + "\n"
    return _atomic_write_text(path, text)


def _parse_stoichiometry(value: Any) -> tuple[tuple[str, float], ...]:
    try:
        payload = json.loads(str(value or "[]"))
    except json.JSONDecodeError as exc:
        raise ValueError("invalid P5 stoichiometry JSON") from exc
    if not isinstance(payload, list):
        raise ValueError("invalid P5 stoichiometry JSON")
    result: list[tuple[str, float]] = []
    for item in payload:
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError("invalid P5 stoichiometry item")
        coefficient = _as_float(item[1], -1.0)
        compound_id = str(item[0] or "").strip()
        if not compound_id or coefficient <= 0:
            raise ValueError("invalid P5 stoichiometry item")
        result.append((compound_id, coefficient))
    return tuple(result)


def _format_amount(value: float) -> str:
    return "" if abs(value - 1.0) <= 1e-12 else f"{value:g} "


def _format_equation(
    left: Sequence[tuple[str, float]],
    right: Sequence[tuple[str, float]],
) -> str:
    def side(rows: Sequence[tuple[str, float]]) -> str:
        return " + ".join(
            f"{_format_amount(amount)}{compound_id}" for compound_id, amount in rows
        )

    return f"{side(left)} => {side(right)}"


def _safe_compound_name(client: KeggRestClient, compound_id: str) -> str:
    if not re.fullmatch(r"C\d{5}", compound_id, re.IGNORECASE):
        return compound_id
    try:
        return client.get_compound_name(compound_id)
    except Exception:
        return compound_id


def _labels(client: KeggRestClient, compound_ids: Sequence[str]) -> str:
    return ";".join(
        f"{compound_id} ({name})" if (name := _safe_compound_name(client, compound_id)) != compound_id else compound_id
        for compound_id in compound_ids
    )


def _load_rules(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {
            str(row.get("Rule ID") or "").strip(): dict(row)
            for row in csv.DictReader(handle)
            if str(row.get("Rule ID") or "").strip()
        }


def _complete_ecs(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({
        normalized
        for value in values
        if (normalized := str(value or "").strip().removeprefix("EC:"))
        and _COMPLETE_EC.fullmatch(normalized)
    }))


def _terms_signature(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[tuple[str, float, str, int], ...] | None:
    aggregated: dict[tuple[str, str, int], float] = defaultdict(float)
    for row in rows:
        side = str(row.get("side") or "").strip().lower()
        inchikey = str(row.get("inchikey") or "").strip().upper()
        compound_id = str(row.get("compound_id") or "").strip()
        identity = f"IK:{inchikey}" if inchikey else f"ID:{compound_id}"
        coefficient = _as_float(row.get("coefficient"), -1.0)
        try:
            charge = int(str(row.get("charge")).strip())
        except (TypeError, ValueError):
            return None
        if side not in {"left", "right"} or identity == "ID:" or coefficient <= 0:
            return None
        aggregated[(side, identity, charge)] += coefficient
    if not aggregated:
        return None
    scale = min(aggregated.values())
    return tuple(sorted(
        (side, round(value / scale, 12), identity, charge)
        for (side, identity, charge), value in aggregated.items()
    ))


def _template_signature(template: Any, mnxref: MnxrefIndex, orientation: str):
    chemicals = mnxref.chemicals(term.mnxm_id for term in template.terms)
    rows: list[dict[str, Any]] = []
    reverse = orientation == "right_to_left"
    for term in template.terms:
        chemical = chemicals.get(term.mnxm_id)
        if chemical is None or chemical.charge is None:
            return None
        side = term.side
        if reverse:
            side = "right" if side == "left" else "left"
        rows.append({
            "side": side,
            "coefficient": term.coefficient,
            "inchikey": chemical.inchikey,
            "compound_id": term.mnxm_id,
            "charge": chemical.charge,
        })
    return _terms_signature(rows)


def _xref_ids(values: Iterable[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    kegg: set[str] = set()
    rhea: set[str] = set()
    for value in values:
        normalized = str(value or "").strip()
        if match := _KEGG_XREF.fullmatch(normalized):
            kegg.add(match.group(1).upper())
        if match := _RHEA_XREF.fullmatch(normalized):
            rhea.add(match.group(1))
    return tuple(sorted(kegg)), tuple(sorted(rhea, key=int))


def _reaction_smiles(rows: Sequence[Mapping[str, Any]]) -> str:
    sides: dict[str, set[str]] = {"left": set(), "right": set()}
    for row in rows:
        side = str(row.get("side") or "").strip().lower()
        smiles = str(row.get("smiles") or "").strip()
        if side not in sides or not smiles:
            return ""
        sides[side].add(smiles)
    if not sides["left"] or not sides["right"]:
        return ""
    return f"{'.'.join(sorted(sides['left']))}>>{'.'.join(sorted(sides['right']))}"


def _term_pairs(
    rows: Sequence[Mapping[str, Any]], side: str
) -> tuple[tuple[str, float], ...]:
    return tuple(
        (str(row.get("compound_id") or "").strip(), _as_float(row.get("coefficient")))
        for row in rows
        if str(row.get("side") or "").strip().lower() == side
    )


def _auxiliary_json(requirement: Any) -> tuple[str, str, str]:
    payload = [
        {
            "role": item.role,
            "necessity": item.necessity,
            "confidence": item.confidence,
            "selection_status": item.selection_status,
            "carrier_ids": list(item.carrier_ids),
            "evidence": list(item.evidence),
        }
        for item in requirement.auxiliary_requirements
    ]
    roles = ";".join(item["role"] for item in payload)
    status = "pending_user_selection" if payload else "not_required"
    return (
        json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
        roles,
        status,
    )


def _electron_row(
    *,
    solution_id: int,
    step_index: int,
    reaction_id: str,
    left: Sequence[tuple[str, float]],
    right: Sequence[tuple[str, float]],
    reaction_name: str,
    reaction_comment: str,
    equation: str,
    enzyme_ecs: Sequence[str],
) -> dict[str, Any]:
    requirement = infer_electron_requirement(
        consumed_stoichiometry=left,
        produced_stoichiometry=right,
        reaction_name=reaction_name,
        reaction_comment=reaction_comment,
        equation=equation,
        enzyme_ecs=enzyme_ecs,
    )
    auxiliary_json, roles, status = _auxiliary_json(requirement)
    return {
        "solution_id": solution_id,
        "step_index": step_index,
        "reaction_id": reaction_id,
        "electron_carrier_ids": ";".join(requirement.carrier_ids),
        "electron_requirement_classes": ";".join(requirement.requirement_classes),
        "electron_risk_level": requirement.risk_level,
        "electron_risk_score": requirement.risk_score,
        "electron_risk_evidence": "; ".join(requirement.evidence),
        "electron_carrier_net_changes": format_electron_carrier_net_changes(
            requirement.carrier_net_stoichiometry
        ),
        "required_auxiliary_roles": roles,
        "auxiliary_requirement_status": status,
        "auxiliary_requirements_json": auxiliary_json,
        "solution_source": "retropath",
    }


def _electron_summary(solution_id: int, rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    summary = summarize_solution_electron_requirements(rows)
    return {"solution_id": solution_id, **summary, "solution_source": "retropath"}


def _candidate_steps_by_id(rows: Sequence[Mapping[str, str]]):
    grouped: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get("candidate_id") or "").strip()].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: _as_int(row.get("step_index"), "step_index", minimum=1))
    return grouped


def _build_solution_rows(
    *,
    solution_id: int,
    passing: Mapping[str, str],
    route: Mapping[str, str],
    steps: Sequence[Mapping[str, str]],
    hypotheses: Mapping[str, Mapping[str, str]],
    terms_by_hypothesis: Mapping[str, Sequence[Mapping[str, str]]],
    rules: Mapping[str, Mapping[str, str]],
    mnxref: MnxrefIndex,
    kegg: KeggRestClient,
    target_compound: str,
    promotion_id: str,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    candidate_rank = _as_int(passing.get("candidate_rank"), "candidate_rank", minimum=1)
    candidate_id = str(passing.get("candidate_id") or "").strip()
    combination_id = str(passing.get("combination_id") or "").strip()
    hypothesis_ids = _split(passing.get("stoichiometry_hypothesis_ids"))
    by_step = {
        str(hypotheses[hypothesis_id].get("step_id") or "").strip(): hypothesis_id
        for hypothesis_id in hypothesis_ids
        if hypothesis_id in hypotheses
    }
    if len(by_step) != len(hypothesis_ids):
        raise ValueError(f"combination {combination_id} references invalid hypotheses")

    step_rows: list[dict[str, Any]] = []
    electron_rows: list[dict[str, Any]] = []
    heterologous_reactions: list[str] = []
    heterologous_ecs: list[str] = []
    for step in steps:
        step_index = _as_int(step.get("step_index"), "step_index", minimum=1)
        step_id = str(step.get("step_id") or "").strip()
        step_source = str(step.get("step_source") or "").strip()
        raw_status = str(step.get("status") or "").strip().lower()
        precursors = _split(step.get("substrate_compound_ids"))
        products = _split(step.get("product_compound_ids"))
        if len(products) != 1:
            raise ValueError(f"{step_id} must have exactly one route product")
        produced = products[0]
        common: dict[str, Any] = {
            "solution_id": solution_id,
            "step_index": step_index,
            "produced_compound_id": produced,
            "produced_compound_name": _safe_compound_name(kegg, produced),
            "direction": "left_to_right",
            "oxygen_required": "false",
            "thermo_direction": "unknown",
            "screening_rule_hits": "",
            "precursor_compound_ids": ";".join(precursors),
            "precursor_compound_labels": _labels(kegg, precursors),
            "ko_ids": "",
            "module_ids": "",
            "resolution_evidence": "strict P8 route flux; exact selected stoichiometry",
            "expansion_depth": str(step.get("expansion_depth") or "0"),
            "expansion_anchor_compounds": str(step.get("sink_anchor_kegg_ids") or ""),
            "prediction_review_required": "true",
            "depends_on_step_ids": str(step.get("depends_on_step_ids") or ""),
        }
        if step_source == "kegg_expansion":
            reaction_ids = tuple(sorted({
                value.upper() for value in _split(step.get("reaction_option_ids"))
                if _KEGG_REACTION.fullmatch(value)
            }))
            if len(reaction_ids) != 1:
                raise ValueError(f"{step_id} has ambiguous KEGG reaction options")
            reaction_id = reaction_ids[0]
            record = kegg.try_get_reaction(reaction_id)
            left = _parse_stoichiometry(step.get("substrate_stoichiometry_json"))
            right = _parse_stoichiometry(step.get("product_stoichiometry_json"))
            if record is not None:
                direction = str(step.get("direction") or "left_to_right").lower()
                reverse = direction in {"right_to_left", "reverse", "backward", "rtl"}
                left = record.right_stoichiometry if reverse else record.left_stoichiometry
                right = record.left_stoichiometry if reverse else record.right_stoichiometry
                name = record.name
                comment = record.comment
                equation = record.equation
                ecs = record.enzyme_ecs
                ko_ids = record.ko_ids
                module_ids = record.module_ids
                rhea_ids = record.rhea_ids
            else:
                name = reaction_id
                comment = "KEGG expansion reaction; full annotation unavailable"
                equation = _format_equation(left, right)
                ecs = _complete_ecs(_split(step.get("source_ec_numbers")))
                ko_ids = tuple()
                module_ids = tuple()
                rhea_ids = tuple()
            status = "endogenous" if raw_status == "endogenous" else "heterologous"
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
                "step_source": "kegg_expansion",
                "rhea_ids": ";".join(rhea_ids),
                "auxiliary_requirements_json": "[]",
            }
        elif step_source == "retropath":
            hypothesis_id = by_step.get(step_id, "")
            if not hypothesis_id:
                raise ValueError(f"combination {combination_id} lacks {step_id}")
            hypothesis = hypotheses[hypothesis_id]
            term_rows = tuple(terms_by_hypothesis.get(hypothesis_id, ()))
            left = _term_pairs(term_rows, "left")
            right = _term_pairs(term_rows, "right")
            if not left or not right:
                raise ValueError(f"hypothesis {hypothesis_id} has incomplete terms")
            rule_id = str(hypothesis.get("rule_id") or "").strip()
            rule = rules.get(rule_id)
            if rule is None:
                raise ValueError(f"RR02 rule missing: {rule_id}")
            source_mnxr_id = str(hypothesis.get("source_mnxr_id") or "").strip()
            template = next(
                (
                    item for item in mnxref.templates_for_rules([rule_id])
                    if item.mnxr_id == source_mnxr_id
                ),
                None,
            )
            if template is None:
                raise ValueError(f"MNXref template missing: {rule_id}/{source_mnxr_id}")
            signature = _terms_signature(term_rows)
            template_signature = _template_signature(
                template,
                mnxref,
                str(hypothesis.get("source_orientation") or ""),
            )
            exact = bool(signature and template_signature and signature == template_signature)
            mapped_kegg, mapped_rhea = _xref_ids(template.reaction_xrefs)
            source_ecs = _complete_ecs((
                *_split(step.get("source_ec_numbers")),
                str(rule.get("EC number") or ""),
            ))
            reaction_signature = _canonical_hash({
                "direction": "biosynthetic",
                "terms": signature,
            })
            equation = _format_equation(left, right)
            terms_payload = [dict(item) for item in term_rows]
            provenance = {
                "promotion_id": promotion_id,
                "candidate_rank": candidate_rank,
                "candidate_id": candidate_id,
                "combination_id": combination_id,
                "step_id": step_id,
                "hypothesis_id": hypothesis_id,
                "rule_id": rule_id,
                "source_mnxr_id": source_mnxr_id,
                "source_reaction_ids": list(dict.fromkeys((
                    *_split(step.get("source_reaction_ids")), source_mnxr_id
                ))),
                "source_ec_numbers": list(source_ecs),
                "source_uniprot_ids": list(_split(step.get("source_uniprot_ids"))),
                "source_rhea_ids": list(mapped_rhea),
                "exact_kegg_reaction_ids": list(mapped_kegg if exact else ()),
                "exact_rhea_ids": list(mapped_rhea if exact else ()),
                "formal_mapping_exact": exact,
                "reaction_signature_sha256": reaction_signature,
                "full_reaction_smiles": _reaction_smiles(term_rows),
                "core_reaction_smiles": str(step.get("reaction_smiles") or "").strip(),
                "review_status": "pending",
            }
            formal_ecs = source_ecs if exact else tuple()
            row = {
                **common,
                "status": "heterologous",
                "reaction_id": hypothesis_id,
                "reaction_name": f"RetroPath prediction {rule_id}",
                "reaction_comment": "Predicted RR02 transformation; experimental review pending",
                "equation": equation,
                "enzyme_ecs": ";".join(formal_ecs),
                "locked_enzyme_ecs": ";".join(formal_ecs),
                "ec_status": "complete" if exact and formal_ecs else "template_only",
                "enzyme_search_eligible": "true",
                "source_reaction_ids": ";".join(provenance["source_reaction_ids"]),
                "resolution_action": "retropath_prediction",
                "step_source": "retropath",
                "rhea_ids": ";".join(provenance["exact_rhea_ids"]),
                "auxiliary_requirements_json": "[]",
                "retropath_step_id": step_id,
                "retropath_hypothesis_id": hypothesis_id,
                "retropath_rule_id": rule_id,
                "source_mnxr_id": source_mnxr_id,
                "source_ec_numbers": ";".join(source_ecs),
                "source_uniprot_ids": ";".join(provenance["source_uniprot_ids"]),
                "source_rhea_ids": ";".join(mapped_rhea),
                "exact_kegg_reaction_ids": ";".join(provenance["exact_kegg_reaction_ids"]),
                "exact_rhea_ids": ";".join(provenance["exact_rhea_ids"]),
                "formal_mapping_exact": str(exact).lower(),
                "reaction_signature_sha256": reaction_signature,
                "full_reaction_smiles": provenance["full_reaction_smiles"],
                "core_reaction_smiles": provenance["core_reaction_smiles"],
                "stoichiometry_terms_json": json.dumps(
                    terms_payload,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "prediction_provenance_json": json.dumps(
                    provenance,
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                ),
            }
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
            enzyme_ecs=_split(row.get("enzyme_ecs")),
        )
        row["auxiliary_requirements_json"] = electron["auxiliary_requirements_json"]
        step_rows.append(row)
        electron_rows.append(electron)
        if row["status"] == "heterologous":
            heterologous_reactions.append(str(row["reaction_id"]))
            heterologous_ecs.extend(_split(row.get("enzyme_ecs")))

    electron_summary = _electron_summary(solution_id, electron_rows)
    anchors = _split(route.get("sink_kegg_ids"))
    target_name = _safe_compound_name(kegg, target_compound)
    summary = {
        "solution_id": solution_id,
        "target_compound_id": target_compound,
        "target_compound_name": target_name,
        "total_steps": len(step_rows),
        "heterologous_steps": sum(row["status"] == "heterologous" for row in step_rows),
        "heterologous_reaction_ids": ";".join(dict.fromkeys(heterologous_reactions)),
        "heterologous_ko_ids": ";".join(dict.fromkeys(
            value for row in step_rows for value in _split(row.get("ko_ids"))
            if row["status"] == "heterologous"
        )),
        "heterologous_enzyme_ecs": ";".join(dict.fromkeys(heterologous_ecs)),
        "reaction_resolution_status": "predicted_strict_stoichiometry",
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
            "C00007" in _split(row.get("precursor_compound_ids")) for row in step_rows
        ),
        "thermo_disfavored_steps": 0,
        **{key: value for key, value in electron_summary.items() if key != "solution_id"},
        "solution_source": "retropath",
        "retropath_candidate_rank": candidate_rank,
        "retropath_candidate_id": candidate_id,
        "retropath_combination_id": combination_id,
        "prediction_review_required": "true",
        "promotion_id": promotion_id,
        "combination_truncated": str(passing.get("combination_truncated") or "false"),
        "upstream_enumeration_truncated": str(
            route.get("upstream_enumeration_truncated") or "false"
        ),
        "candidate_top_k_truncated": str(
            route.get("candidate_top_k_truncated") or "false"
        ),
    }
    return summary, step_rows, electron_summary, electron_rows


def _verified_artifact(
    p8_dir: Path,
    manifest: Mapping[str, Any],
    name: str,
) -> Path:
    artifacts = manifest.get("artifacts")
    record = artifacts.get(name) if isinstance(artifacts, Mapping) else None
    path = (p8_dir / name).resolve()
    if not isinstance(record, Mapping):
        raise ValueError(f"P8 artifact missing: {name}")
    if Path(str(record.get("path") or "")).expanduser().resolve() != path:
        raise ValueError(f"P8 artifact path stale: {name}")
    if not path.is_file() or record.get("sha256") != _sha256_file(path):
        raise ValueError(f"P8 artifact checksum mismatch: {name}")
    return path


def _load_validation_manifest(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("P8 validation manifest is invalid") from exc
    if not isinstance(payload, dict):
        raise ValueError("P8 validation manifest is invalid")
    return payload


def materialize_retropath_solutions(
    config: Any,
    *,
    validation_manifest_path: str | Path,
) -> dict[str, Any]:
    """Replace the formal RetroPath slice with the current strict P8 passes."""

    target = validate_target_compound_id(config.target_name)
    depth = _as_int(getattr(config, "depth", 0), "depth")
    gap_dir = gap_depth_output_dir(Path(config.gap_output_path), depth).resolve()
    retropath_dir = gap_dir / "retropath"
    p8_dir = retropath_dir / "gem_validation"
    manifest_path = Path(validation_manifest_path).expanduser().resolve()
    if manifest_path != (p8_dir / VALIDATION_MANIFEST_FILE_NAME).resolve():
        raise ValueError("unexpected P8 validation manifest path")
    p8_manifest = _load_validation_manifest(manifest_path)
    if (
        p8_manifest.get("target_compound") != target
        or _as_int(p8_manifest.get("expansion_depth"), "expansion_depth") != depth
    ):
        raise ValueError("P8 validation manifest identity mismatch")

    p8_summary_path = _verified_artifact(
        p8_dir, p8_manifest, VALIDATION_SUMMARY_FILE_NAME
    )
    hypothesis_path = _verified_artifact(
        p8_dir, p8_manifest, STOICHIOMETRY_HYPOTHESES_FILE_NAME
    )
    term_path = _verified_artifact(
        p8_dir, p8_manifest, STOICHIOMETRY_TERMS_FILE_NAME
    )
    passing = [
        row for row in _read_csv(p8_summary_path, SUMMARY_COLUMNS)
        if row.get("validation_status") == "PASS_STRICT_ROUTE_FLUX"
        and str(row.get("combination_id") or "").strip()
    ]
    passing.sort(key=lambda row: (
        _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1),
        str(row.get("combination_id") or ""),
    ))
    routes = _read_csv(retropath_dir / CANDIDATE_ROUTES_FILE_NAME, ("candidate_id",))
    route_by_id = {
        str(row.get("candidate_id") or "").strip(): row for row in routes
    }
    candidate_steps = _candidate_steps_by_id(
        _read_csv(retropath_dir / CANDIDATE_STEPS_FILE_NAME, ("candidate_id", "step_id"))
    )
    hypotheses = {
        str(row.get("hypothesis_id") or "").strip(): row
        for row in _read_csv(hypothesis_path, HYPOTHESIS_COLUMNS)
    }
    terms_by_hypothesis: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in _read_csv(term_path, TERM_COLUMNS):
        terms_by_hypothesis[str(row.get("hypothesis_id") or "").strip()].append(row)

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
    kegg_ids = [
        _as_int(row.get("solution_id"), "solution_id", minimum=1)
        for row in kegg_summaries
    ]
    next_solution_id = max(kegg_ids, default=0) + 1
    valid_kegg_ids = set(kegg_ids)
    kegg_steps = [
        row for row in existing_steps
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) in valid_kegg_ids
    ]
    kegg_electron_summaries = [
        row for row in existing_electron_summaries
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) in valid_kegg_ids
    ]
    kegg_electron_steps = [
        row for row in existing_electron_steps
        if _as_int(row.get("solution_id"), "solution_id", minimum=1) in valid_kegg_ids
    ]

    promotion_seed = {
        "schema_version": RETROPATH_PROMOTION_SCHEMA,
        "target_compound": target,
        "expansion_depth": depth,
        "p8_validation_sha256": _sha256_file(manifest_path),
        "passing_combinations": [
            {
                "candidate_rank": _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1),
                "candidate_id": row.get("candidate_id"),
                "combination_id": row.get("combination_id"),
                "hypothesis_ids": list(_split(row.get("stoichiometry_hypothesis_ids"))),
            }
            for row in passing
        ],
    }
    promotion_id = f"RP2PROMOTE:{_canonical_hash(promotion_seed)}"
    rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    rules = _load_rules(rules_path)
    kegg = KeggRestClient(Path(config.cache_dir).resolve() / "kegg")
    mnxref_dir = Path(config.data_dir) / "retropath" / "mnxref" / "3.0"
    rp_summaries: list[dict[str, Any]] = []
    rp_steps: list[dict[str, Any]] = []
    rp_electron_summaries: list[dict[str, Any]] = []
    rp_electron_steps: list[dict[str, Any]] = []
    mappings: list[dict[str, Any]] = []
    with MnxrefIndex(mnxref_dir, rules_path) as mnxref:
        for offset, row in enumerate(passing):
            candidate_id = str(row.get("candidate_id") or "").strip()
            route = route_by_id.get(candidate_id)
            steps = candidate_steps.get(candidate_id)
            if route is None or not steps:
                raise ValueError(f"P5 candidate inputs missing for {candidate_id}")
            solution_id = next_solution_id + offset
            summary, step_rows, electron_summary, electron_rows = _build_solution_rows(
                solution_id=solution_id,
                passing=row,
                route=route,
                steps=steps,
                hypotheses=hypotheses,
                terms_by_hypothesis=terms_by_hypothesis,
                rules=rules,
                mnxref=mnxref,
                kegg=kegg,
                target_compound=target,
                promotion_id=promotion_id,
            )
            rp_summaries.append(summary)
            rp_steps.extend(step_rows)
            rp_electron_summaries.append(electron_summary)
            rp_electron_steps.extend(electron_rows)
            mappings.append({
                "solution_id": solution_id,
                "candidate_rank": _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1),
                "candidate_id": candidate_id,
                "combination_id": str(row.get("combination_id") or ""),
                "stoichiometry_hypothesis_ids": list(
                    _split(row.get("stoichiometry_hypothesis_ids"))
                ),
            })

    for row in kegg_summaries:
        row.setdefault("solution_source", "kegg")
    summary_rows = [*kegg_summaries, *rp_summaries]
    step_rows = [*kegg_steps, *rp_steps]
    electron_summary_rows = [*kegg_electron_summaries, *rp_electron_summaries]
    electron_step_rows = [*kegg_electron_steps, *rp_electron_steps]
    texts = {
        "solutions.csv": _render_csv(
            _field_order(summary_rows, _SUMMARY_BASE_COLUMNS, _SUMMARY_RETROPATH_COLUMNS),
            summary_rows,
        ),
        "all_solution_steps.csv": _render_csv(
            _field_order(step_rows, _STEP_BASE_COLUMNS, _STEP_RETROPATH_COLUMNS),
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
    output_hashes: dict[str, str] = {}
    for name, text in texts.items():
        output_hashes[name] = _atomic_write_text(gap_dir / name, text)

    promotion_manifest_path = retropath_dir / PROMOTION_MANIFEST_FILE_NAME
    promotion_manifest = {
        **promotion_seed,
        "promotion_id": promotion_id,
        "created_at": datetime.now(UTC).isoformat(),
        "formal_promotion_allowed": bool(mappings),
        "formal_solution_count": len(mappings),
        "solution_mappings": mappings,
        "inputs": {
            "p8_validation_manifest": str(manifest_path),
            "p8_validation_sha256": _sha256_file(manifest_path),
            "rr02_path": str(rules_path),
            "rr02_sha256": _sha256_file(rules_path),
            "p8_inputs": p8_manifest.get("inputs", {}),
        },
        "artifacts": {
            name: {"path": str((gap_dir / name).resolve()), "sha256": digest}
            for name, digest in output_hashes.items()
        },
    }
    promotion_manifest_sha256 = _atomic_write_json(
        promotion_manifest_path,
        promotion_manifest,
    )
    return {
        "ok": True,
        "schema_version": RETROPATH_PROMOTION_SCHEMA,
        "promotion_id": promotion_id,
        "promotion_manifest": str(promotion_manifest_path.resolve()),
        "promotion_manifest_sha256": promotion_manifest_sha256,
        "formal_solution_count": len(mappings),
        "formal_solution_ids": [item["solution_id"] for item in mappings],
        "solution_mappings": mappings,
        "output_hashes": output_hashes,
    }


def verify_retropath_solution_promotion(
    *,
    gap_dir: str | Path,
    target_compound: str,
    expansion_depth: int,
    solution_id: int | None = None,
) -> dict[str, Any]:
    """Verify the commit marker and optionally require one mapped solution."""

    resolved_gap = Path(gap_dir).expanduser().resolve()
    path = resolved_gap / "retropath" / PROMOTION_MANIFEST_FILE_NAME
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("RetroPath formal solution promotion manifest is invalid") from exc
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != RETROPATH_PROMOTION_SCHEMA
        or payload.get("target_compound") != target_compound
        or _as_int(payload.get("expansion_depth"), "expansion_depth") != expansion_depth
    ):
        raise ValueError("RetroPath formal solution promotion identity mismatch")
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, Mapping):
        raise ValueError("RetroPath promotion artifacts are missing")
    for name in (
        "solutions.csv",
        "all_solution_steps.csv",
        "solution_electron_summary.csv",
        "route_electron_requirements.csv",
    ):
        record = artifacts.get(name)
        expected = (resolved_gap / name).resolve()
        if (
            not isinstance(record, Mapping)
            or Path(str(record.get("path") or "")).expanduser().resolve() != expected
            or not expected.is_file()
            or record.get("sha256") != _sha256_file(expected)
        ):
            raise ValueError(f"RetroPath promotion artifact mismatch: {name}")
    p8_record = payload.get("inputs", {}).get("p8_validation_manifest")
    p8_sha = payload.get("inputs", {}).get("p8_validation_sha256")
    p8_path = Path(str(p8_record or "")).expanduser().resolve()
    if not p8_path.is_file() or p8_sha != _sha256_file(p8_path):
        raise ValueError("RetroPath P8 validation binding is stale")
    p8_manifest = _load_validation_manifest(p8_path)
    p8_artifacts = p8_manifest.get("artifacts")
    if not isinstance(p8_artifacts, Mapping):
        raise ValueError("RetroPath P8 validation artifacts are missing")
    for name, record in p8_artifacts.items():
        if not isinstance(record, Mapping):
            raise ValueError(f"RetroPath P8 artifact record is invalid: {name}")
        artifact_path = Path(str(record.get("path") or "")).expanduser().resolve()
        if (
            not artifact_path.is_file()
            or record.get("sha256") != _sha256_file(artifact_path)
        ):
            raise ValueError(f"RetroPath P8 artifact mismatch: {name}")
    if solution_id is not None:
        matches = [
            row for row in payload.get("solution_mappings", [])
            if isinstance(row, Mapping)
            and _as_int(row.get("solution_id"), "solution_id", minimum=1) == solution_id
        ]
        if len(matches) != 1:
            raise ValueError(f"RetroPath solution {solution_id} is not committed")
    return payload


__all__ = [
    "PROMOTION_MANIFEST_FILE_NAME",
    "RETROPATH_PROMOTION_SCHEMA",
    "materialize_retropath_solutions",
    "verify_retropath_solution_promotion",
]
