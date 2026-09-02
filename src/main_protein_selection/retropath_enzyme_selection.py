"""Enzyme candidate retrieval for optionally P8-validated RetroPath routes."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import requests
from rdkit import Chem

from src.main_protein_selection.reaction_aware_retrieval import (
    RheaClient,
    retrieve_rhea_candidates_for_requirement,
)
from src.main_protein_selection.selenzyme_retrieval import (
    COMPLETE_EC_PATTERN,
    SelenzymeClient,
    SelenzymeSourceUnavailable,
    chassis_host_taxon_id,
    retrieve_selenzyme_candidates,
    selenzyme_target_count,
)
from src.main_protein_selection.sequence_quality import analyze_protein_sequence
from src.main_protein_selection.taxonomy_compatibility import (
    ChassisTaxonomyProfile,
    resolve_chassis_taxonomy,
)
from src.main_protein_selection.uniprot_protein_candidates import (
    ProteinCandidate,
    candidate_from_reaction_entry,
    recommend_uniprot_proteins,
    resolve_uniprot_accession_batches,
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
    _load_candidate_inputs,
)
from src.pathway_analyze.retropath_mnxref import (
    MnxrefIndex,
    MnxrefReactionTemplate,
)

RETROPATH_ENZYME_SELECTION_SCHEMA = "retropath_enzyme_selection.v1"
RETROPATH_ENZYME_REQUIREMENTS_SCHEMA = "retropath_enzyme_requirements.v1"
RETROPATH_SELENZYME_EVIDENCE_SCHEMA = "retropath_selenzyme_evidence.v1"

ENZYME_SELECTION_DIR_NAME = "enzyme_selection"
ENZYME_REQUIREMENTS_FILE_NAME = "retropath_enzyme_requirements.json"
STEP_ENZYME_CANDIDATES_FILE_NAME = "retropath_step_enzyme_candidates.csv"
STEP_ENZYME_AUDIT_FILE_NAME = "retropath_step_enzyme_candidate_audit.csv"
SELENZYME_EVIDENCE_FILE_NAME = "retropath_selenzyme_evidence.json"
ENZYME_SELECTION_FILE_NAME = "retropath_enzyme_selection.json"

EVIDENCE_TIER_ORDER = {
    "formal_reaction_exact": 0,
    "source_template_supported": 1,
    "full_reaction_similarity": 2,
    "core_reaction_similarity": 3,
    "rule_similarity": 4,
    "rejected": 9,
}

ENZYME_CANDIDATE_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "combination_rank",
    "combination_id",
    "step_index",
    "step_id",
    "step_source",
    "hypothesis_id",
    "protein_candidate_rank",
    "accession",
    "entry_name",
    "protein_name",
    "organism_name",
    "organism_id",
    "reviewed",
    "length",
    "evidence_tier",
    "evidence_tiers",
    "fit_status",
    "manual_review_required",
    "selection_status",
    "reaction_signature_sha256",
    "reaction_id",
    "source_mnxr_id",
    "rule_id",
    "ec_numbers",
    "matched_rhea_ids",
    "retrieval_strategies",
    "retrieval_query_ids",
    "reaction_similarity",
    "sim_rf",
    "sim_2018",
    "matched_reaction_id",
    "direction_verdict",
    "direction_confidence",
    "protein_score",
    "taxonomic_lineage",
    "taxonomic_lineage_ids",
    "taxonomic_shared_taxon_id",
    "taxonomic_shared_name",
    "taxonomic_shared_rank",
    "taxonomic_fit_status",
    "taxonomic_fit_score",
    "taxonomy_evidence_source",
    "taxonomic_distance",
    "sequence_sha256",
    "sequence",
    "sequence_version",
    "score_breakdown",
    "candidate_role",
    "publication_ids",
    "catalytic_activities",
    "catalytic_activity_records_json",
    "cofactors",
    "function_comments",
    "enzyme_system_type",
    "auxiliary_requirement_status",
    "auxiliary_requirements_json",
    "warnings",
    "reasons",
    "rejection_reasons",
)

_RHEA_XREF = re.compile(r"^rhea:(\d+)$", re.IGNORECASE)
_KEGG_XREF = re.compile(r"^kegg:(R\d{5})$", re.IGNORECASE)


@dataclass(frozen=True)
class RetropathEnzymeRequirement:
    """One enzyme-search task bound to an exact P8 route hypothesis."""

    candidate_rank: int
    candidate_id: str
    combination_id: str
    step_index: int
    step_id: str
    step_source: str
    step_status: str
    hypothesis_id: str
    reaction_signature_sha256: str
    full_reaction_smiles: str
    core_reaction_smiles: str
    rule_id: str
    rule_smarts: str
    source_mnxr_id: str
    source_reaction_ids: tuple[str, ...]
    source_ec_numbers: tuple[str, ...]
    source_uniprot_ids: tuple[str, ...]
    source_rhea_ids: tuple[str, ...]
    exact_kegg_reaction_ids: tuple[str, ...]
    exact_rhea_ids: tuple[str, ...]
    formal_mapping_exact: bool
    enzyme_required: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class _ValidatedP9Inputs:
    candidate_inputs: Any
    candidate_rank: int
    candidate_id: str
    route: Mapping[str, str]
    steps: tuple[Mapping[str, str], ...]
    p8_dir: Path
    p8_manifest_path: Path
    p8_manifest: Mapping[str, Any]
    passing_rows: tuple[Mapping[str, str], ...]
    hypotheses: tuple[Mapping[str, str], ...]
    terms: tuple[Mapping[str, str], ...]


@dataclass
class _SearchOutcome:
    candidates: list[dict[str, Any]]
    audit: list[dict[str, Any]]
    evidence: list[dict[str, Any]]
    source_unavailable: bool = False


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _atomic_write_text(path: Path, text: str) -> str:
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
            handle.write(text)
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
    return _sha256_file(path)


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> str:
    return _atomic_write_text(
        path,
        json.dumps(
            value,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n",
    )


def _render_csv(
    columns: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=tuple(columns),
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "") for key in columns})
    return buffer.getvalue()


def _read_csv(path: Path, columns: Sequence[str]) -> tuple[dict[str, str], ...]:
    if not path.is_file():
        raise FileNotFoundError(f"required P8 artifact not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise ValueError(
                f"P8 artifact schema mismatch; rerun RetroPath validation: {path}"
            )
        return tuple(dict(row) for row in reader)


def _as_int(value: Any, field_name: str, minimum: int = 0) -> int:
    try:
        result = int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer") from exc
    if result < minimum:
        raise ValueError(f"{field_name} must be >= {minimum}")
    return result


def _as_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _split(value: Any) -> tuple[str, ...]:
    if isinstance(value, (list, tuple, set)):
        values = value
    else:
        values = str(value or "").replace("|", ";").split(";")
    return tuple(dict.fromkeys(
        str(item).strip() for item in values if str(item).strip()
    ))


def _complete_ecs(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({
        str(value).strip().removeprefix("EC:").removeprefix("ec:")
        for value in values
        if COMPLETE_EC_PATTERN.fullmatch(
            str(value).strip().removeprefix("EC:").removeprefix("ec:")
        )
    }))


def _verified_p8_artifact(
    p8_dir: Path,
    manifest: Mapping[str, Any],
    name: str,
) -> Path:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, Mapping):
        raise ValueError("P8 manifest artifacts are missing")
    record = artifacts.get(name)
    if not isinstance(record, Mapping):
        raise ValueError(f"P8 artifact record is missing: {name}")
    expected = (p8_dir / name).resolve()
    if Path(str(record.get("path") or "")).expanduser().resolve() != expected:
        raise ValueError(f"P8 artifact path is stale: {name}")
    if not expected.is_file() or record.get("sha256") != _sha256_file(expected):
        raise ValueError(f"P8 artifact checksum mismatch: {name}")
    return expected


def _load_p9_inputs(config: Any) -> _ValidatedP9Inputs:
    candidate_inputs = _load_candidate_inputs(config)
    raw_rank = getattr(config, "retropath_candidate", None)
    if raw_rank is None:
        raise ValueError(
            "--retropath-candidate is required for RetroPath enzyme selection"
        )
    rank = _as_int(raw_rank, "retropath_candidate", minimum=1)
    route = next(
        (
            row
            for row in candidate_inputs.routes
            if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1)
            == rank
        ),
        None,
    )
    if route is None:
        available = [
            _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1)
            for row in candidate_inputs.routes
        ]
        raise ValueError(
            f"RetroPath candidate {rank} not found; available candidates: {available}"
        )
    candidate_id = str(route.get("candidate_id") or "").strip()
    steps = tuple(sorted(
        (
            row
            for row in candidate_inputs.steps
            if str(row.get("candidate_id") or "").strip() == candidate_id
        ),
        key=lambda row: _as_int(row.get("step_index"), "step_index", minimum=1),
    ))
    if not steps:
        raise ValueError("selected RetroPath candidate has no steps")

    p8_dir = candidate_inputs.retropath_dir / "gem_validation"
    manifest_path = p8_dir / VALIDATION_MANIFEST_FILE_NAME
    if not manifest_path.is_file():
        raise FileNotFoundError(
            "P8 validation is missing; run validate -s <solution_id> first"
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("P8 validation manifest is invalid") from exc
    if (
        not isinstance(manifest, Mapping)
        or manifest.get("schema_version") != RETROPATH_GEM_VALIDATION_SCHEMA
        or manifest.get("target_compound") != candidate_inputs.target_compound
        or _as_int(manifest.get("expansion_depth"), "P8 expansion_depth")
        != candidate_inputs.depth
    ):
        raise ValueError("P8 validation manifest identity mismatch")
    p8_inputs = manifest.get("inputs")
    if not isinstance(p8_inputs, Mapping):
        raise ValueError("P8 validation inputs are missing")
    pipeline_record = p8_inputs.get("pipeline_result")
    if (
        not isinstance(pipeline_record, Mapping)
        or pipeline_record.get("sha256")
        != _sha256_file(candidate_inputs.pipeline_path)
    ):
        raise ValueError("P8 pipeline input is stale")
    pipeline_artifacts = candidate_inputs.pipeline.get("artifacts")
    if not isinstance(pipeline_artifacts, Mapping):
        raise ValueError("RetroPath pipeline artifacts are missing")
    route_record = pipeline_artifacts.get("candidate_routes")
    step_record = pipeline_artifacts.get("candidate_steps")
    if (
        not isinstance(route_record, Mapping)
        or not isinstance(step_record, Mapping)
        or p8_inputs.get("candidate_routes_sha256") != route_record.get("sha256")
        or p8_inputs.get("candidate_steps_sha256") != step_record.get("sha256")
    ):
        raise ValueError("P8 candidate inputs are stale")
    statuses = manifest.get("candidate_statuses")
    if (
        not isinstance(statuses, Mapping)
        or str(statuses.get(str(rank)) or "") not in {
            "PASS_STRICT_HYPOTHESIS_EXISTS",
            "PASS_RELAXED_HYPOTHESIS_EXISTS",
        }
    ):
        raise ValueError(
            f"RetroPath candidate {rank} has no P8 GEM-passing hypothesis"
        )

    summary_path = _verified_p8_artifact(
        p8_dir,
        manifest,
        VALIDATION_SUMMARY_FILE_NAME,
    )
    hypotheses_path = _verified_p8_artifact(
        p8_dir,
        manifest,
        STOICHIOMETRY_HYPOTHESES_FILE_NAME,
    )
    terms_path = _verified_p8_artifact(
        p8_dir,
        manifest,
        STOICHIOMETRY_TERMS_FILE_NAME,
    )
    summary = _read_csv(summary_path, SUMMARY_COLUMNS)
    passing = tuple(
        row
        for row in summary
        if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1) == rank
        and row.get("candidate_id") == candidate_id
        and row.get("validation_status") in PASSING_ROUTE_VALIDATION_STATUSES
        and str(row.get("combination_id") or "").strip()
    )
    if not passing:
        raise ValueError(
            "P8 manifest claims success but contains no passing combination"
        )
    hypotheses = tuple(
        row
        for row in _read_csv(hypotheses_path, HYPOTHESIS_COLUMNS)
        if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1) == rank
        and row.get("candidate_id") == candidate_id
    )
    terms = tuple(
        row
        for row in _read_csv(terms_path, TERM_COLUMNS)
        if _as_int(row.get("candidate_rank"), "candidate_rank", minimum=1) == rank
        and row.get("candidate_id") == candidate_id
    )
    known_hypotheses = {row.get("hypothesis_id") for row in hypotheses}
    for row in passing:
        ids = _split(row.get("stoichiometry_hypothesis_ids"))
        if not ids or not set(ids).issubset(known_hypotheses):
            raise ValueError("P8 passing combination references unknown hypotheses")
    return _ValidatedP9Inputs(
        candidate_inputs=candidate_inputs,
        candidate_rank=rank,
        candidate_id=candidate_id,
        route=route,
        steps=steps,
        p8_dir=p8_dir,
        p8_manifest_path=manifest_path,
        p8_manifest=manifest,
        passing_rows=passing,
        hypotheses=hypotheses,
        terms=terms,
    )


def _canonical_smiles(smiles: str, inchi: str) -> str:
    molecule = Chem.MolFromSmiles(str(smiles or "").strip()) if smiles else None
    if molecule is None and inchi:
        molecule = Chem.MolFromInchi(str(inchi).strip())
    if molecule is None:
        return ""
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _structural_terms_signature(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[tuple[str, float, str, int], ...] | None:
    aggregated: dict[tuple[str, str, int], float] = defaultdict(float)
    for row in rows:
        side = str(row.get("side") or "").strip().lower()
        inchikey = str(row.get("inchikey") or "").strip().upper()
        coefficient = _as_float(row.get("coefficient"))
        charge_value = row.get("charge")
        try:
            charge = int(str(charge_value).strip())
        except (TypeError, ValueError):
            return None
        if side not in {"left", "right"} or not inchikey or not coefficient:
            return None
        if coefficient <= 0:
            return None
        aggregated[(side, inchikey, charge)] += coefficient
    if not aggregated:
        return None
    scale = min(aggregated.values())
    return tuple(sorted(
        (side, round(value / scale, 12), inchikey, charge)
        for (side, inchikey, charge), value in aggregated.items()
    ))


def _stable_terms_signature(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[tuple[str, float, str, int], ...] | None:
    """Return a stable signature even when a participant lacks structure."""

    aggregated: dict[tuple[str, str, int], float] = defaultdict(float)
    for row in rows:
        side = str(row.get("side") or "").strip().lower()
        inchikey = str(row.get("inchikey") or "").strip().upper()
        compound_id = str(row.get("compound_id") or "").strip()
        identity = f"IK:{inchikey}" if inchikey else f"ID:{compound_id}"
        coefficient = _as_float(row.get("coefficient"))
        try:
            charge = int(str(row.get("charge")).strip())
        except (TypeError, ValueError):
            return None
        if (
            side not in {"left", "right"}
            or identity == "ID:"
            or coefficient is None
            or coefficient <= 0
        ):
            return None
        aggregated[(side, identity, charge)] += coefficient
    if not aggregated:
        return None
    scale = min(aggregated.values())
    return tuple(sorted(
        (side, round(value / scale, 12), identity, charge)
        for (side, identity, charge), value in aggregated.items()
    ))


def _full_reaction_smiles(rows: Sequence[Mapping[str, Any]]) -> str:
    sides: dict[str, set[str]] = {"left": set(), "right": set()}
    for row in rows:
        side = str(row.get("side") or "").strip().lower()
        if side not in sides:
            return ""
        smiles = _canonical_smiles(
            str(row.get("smiles") or ""),
            str(row.get("inchi") or ""),
        )
        if not smiles:
            return ""
        sides[side].add(smiles)
    if not sides["left"] or not sides["right"]:
        return ""
    return f"{'.'.join(sorted(sides['left']))}>>{'.'.join(sorted(sides['right']))}"


def _template_signature(
    template: MnxrefReactionTemplate,
    mnxref: MnxrefIndex,
    orientation: str,
) -> tuple[tuple[str, float, str, int], ...] | None:
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
            "charge": chemical.charge,
        })
    return _structural_terms_signature(rows)


def _xref_ids(values: Iterable[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    kegg: set[str] = set()
    rhea: set[str] = set()
    for value in values:
        normalized = str(value or "").strip()
        kegg_match = _KEGG_XREF.fullmatch(normalized)
        if kegg_match:
            kegg.add(kegg_match.group(1).upper())
        rhea_match = _RHEA_XREF.fullmatch(normalized)
        if rhea_match:
            rhea.add(rhea_match.group(1))
    return tuple(sorted(kegg)), tuple(sorted(rhea, key=int))


def _load_rules(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        result: dict[str, dict[str, str]] = {}
        for row in rows:
            rule_id = str(row.get("Rule ID") or "").strip()
            if rule_id:
                result[rule_id] = dict(row)
    return result


def _build_requirements(
    inputs: _ValidatedP9Inputs,
    mnxref: MnxrefIndex,
    rules: Mapping[str, Mapping[str, str]],
) -> list[RetropathEnzymeRequirement]:
    hypotheses = {
        str(row.get("hypothesis_id") or ""): row for row in inputs.hypotheses
    }
    terms_by_hypothesis: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in inputs.terms:
        terms_by_hypothesis[str(row.get("hypothesis_id") or "")].append(row)

    requirements: list[RetropathEnzymeRequirement] = []
    for passing in inputs.passing_rows:
        combination_id = str(passing.get("combination_id") or "").strip()
        combination_hypotheses = {
            str(hypotheses[hypothesis_id].get("step_id") or ""): hypothesis_id
            for hypothesis_id in _split(
                passing.get("stoichiometry_hypothesis_ids")
            )
        }
        for step in inputs.steps:
            step_source = str(step.get("step_source") or "").strip()
            step_id = str(step.get("step_id") or "").strip()
            step_status = str(step.get("status") or "").strip()
            step_index = _as_int(step.get("step_index"), "step_index", minimum=1)
            if step_source == "kegg_expansion":
                reaction_ids = tuple(sorted({
                    value.upper()
                    for value in _split(step.get("reaction_option_ids"))
                    if re.fullmatch(r"R\d{5}", value.upper())
                }))
                signature = _canonical_sha256({
                    "step_id": step_id,
                    "reaction_ids": reaction_ids,
                    "direction": str(step.get("direction") or ""),
                })
                requirements.append(RetropathEnzymeRequirement(
                    candidate_rank=inputs.candidate_rank,
                    candidate_id=inputs.candidate_id,
                    combination_id=combination_id,
                    step_index=step_index,
                    step_id=step_id,
                    step_source=step_source,
                    step_status=step_status,
                    hypothesis_id="",
                    reaction_signature_sha256=signature,
                    full_reaction_smiles="",
                    core_reaction_smiles="",
                    rule_id="",
                    rule_smarts="",
                    source_mnxr_id="",
                    source_reaction_ids=_split(step.get("source_reaction_ids")),
                    source_ec_numbers=_complete_ecs(
                        _split(step.get("source_ec_numbers"))
                    ),
                    source_uniprot_ids=_split(step.get("source_uniprot_ids")),
                    source_rhea_ids=tuple(),
                    exact_kegg_reaction_ids=reaction_ids,
                    exact_rhea_ids=tuple(),
                    formal_mapping_exact=True,
                    enzyme_required=step_status != "endogenous",
                ))
                continue
            if step_source != "retropath":
                raise ValueError(f"unsupported hybrid step source: {step_source}")
            hypothesis_id = combination_hypotheses.get(step_id, "")
            if not hypothesis_id:
                raise ValueError(
                    f"P8 passing combination {combination_id} lacks hypothesis "
                    f"for {step_id}"
                )
            hypothesis = hypotheses[hypothesis_id]
            term_rows = terms_by_hypothesis.get(hypothesis_id, [])
            signature_terms = _structural_terms_signature(term_rows)
            stable_signature_terms = _stable_terms_signature(term_rows)
            if stable_signature_terms is None:
                raise ValueError(
                    f"P8 hypothesis identity is incomplete: {hypothesis_id}"
                )
            signature = _canonical_sha256({
                "direction": "biosynthetic",
                "terms": stable_signature_terms,
            })
            rule_id = str(hypothesis.get("rule_id") or "").strip()
            rule = rules.get(rule_id)
            if rule is None:
                raise ValueError(f"RR02 rule is missing for P8 hypothesis: {rule_id}")
            source_mnxr_id = str(hypothesis.get("source_mnxr_id") or "").strip()
            template = next(
                (
                    item
                    for item in mnxref.templates_for_rules([rule_id])
                    if item.mnxr_id == source_mnxr_id
                ),
                None,
            )
            if template is None:
                raise ValueError(
                    f"MNXref template is missing for {rule_id}/{source_mnxr_id}"
                )
            template_signature = _template_signature(
                template,
                mnxref,
                str(hypothesis.get("source_orientation") or ""),
            )
            exact = (
                signature_terms is not None
                and template_signature is not None
                and template_signature == signature_terms
            )
            mapped_kegg, mapped_rhea = _xref_ids(template.reaction_xrefs)
            source_ecs = _complete_ecs([
                *_split(step.get("source_ec_numbers")),
                str(rule.get("EC number") or ""),
            ])
            requirements.append(RetropathEnzymeRequirement(
                candidate_rank=inputs.candidate_rank,
                candidate_id=inputs.candidate_id,
                combination_id=combination_id,
                step_index=step_index,
                step_id=step_id,
                step_source=step_source,
                step_status=step_status,
                hypothesis_id=hypothesis_id,
                reaction_signature_sha256=signature,
                full_reaction_smiles=_full_reaction_smiles(term_rows),
                core_reaction_smiles=str(step.get("reaction_smiles") or "").strip(),
                rule_id=rule_id,
                rule_smarts=str(rule.get("Rule") or "").strip(),
                source_mnxr_id=source_mnxr_id,
                source_reaction_ids=tuple(dict.fromkeys([
                    *_split(step.get("source_reaction_ids")),
                    source_mnxr_id,
                ])),
                source_ec_numbers=source_ecs,
                source_uniprot_ids=_split(step.get("source_uniprot_ids")),
                source_rhea_ids=mapped_rhea,
                exact_kegg_reaction_ids=mapped_kegg if exact else tuple(),
                exact_rhea_ids=mapped_rhea if exact else tuple(),
                formal_mapping_exact=exact,
                enzyme_required=True,
            ))
    return requirements


def _candidate_record(
    candidate: ProteinCandidate,
    *,
    evidence_tier: str,
    fit_status: str,
    query_id: str = "",
) -> dict[str, Any]:
    return {
        "accession": candidate.accession.upper(),
        "entry_name": candidate.entry_name,
        "protein_name": candidate.protein_name,
        "organism_name": candidate.organism_name,
        "organism_id": candidate.organism_id or "",
        "reviewed": str(candidate.reviewed).lower(),
        "length": candidate.length or "",
        "evidence_tier": evidence_tier,
        "evidence_tiers": evidence_tier,
        "fit_status": fit_status,
        "manual_review_required": str(fit_status == "manual_review").lower(),
        "selection_status": "eligible",
        "ec_numbers": ";".join(sorted(set(candidate.ec_numbers))),
        "matched_rhea_ids": ";".join(sorted(set(candidate.matched_rhea_ids))),
        "retrieval_strategies": candidate.retrieval_strategy,
        "retrieval_query_ids": query_id or candidate.retrieval_query_id,
        "reaction_similarity": (
            ""
            if candidate.selenzyme_reaction_similarity is None
            else candidate.selenzyme_reaction_similarity
        ),
        "sim_rf": (
            "" if candidate.selenzyme_sim_rf is None else candidate.selenzyme_sim_rf
        ),
        "sim_2018": (
            "" if candidate.selenzyme_sim_2018 is None else candidate.selenzyme_sim_2018
        ),
        "matched_reaction_id": candidate.selenzyme_matched_reaction_id,
        "direction_verdict": candidate.direction_verdict or "unknown",
        "direction_confidence": candidate.direction_confidence,
        "protein_score": candidate.score,
        "taxonomic_fit_status": candidate.taxonomic_fit_status,
        "taxonomic_fit_score": candidate.taxonomic_fit_score,
        "taxonomic_shared_taxon_id": candidate.taxonomic_shared_taxon_id or "",
        "taxonomic_shared_name": candidate.taxonomic_shared_name,
        "taxonomic_shared_rank": candidate.taxonomic_shared_rank,
        "taxonomy_evidence_source": candidate.taxonomy_evidence_source,
        "taxonomic_lineage": ";".join(candidate.taxonomic_lineage),
        "taxonomic_lineage_ids": ";".join(
            str(value) for value in candidate.taxonomic_lineage_ids
        ),
        "taxonomic_distance": (
            ""
            if candidate.selenzyme_taxonomic_distance is None
            else candidate.selenzyme_taxonomic_distance
        ),
        "sequence_sha256": candidate.sequence_sha256,
        "sequence": candidate.sequence,
        "sequence_version": candidate.sequence_version or "",
        "score_breakdown": json.dumps(
            candidate.score_breakdown,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ),
        "candidate_role": candidate.candidate_role,
        "publication_ids": ";".join(candidate.publication_ids),
        "catalytic_activities": " | ".join(candidate.catalytic_activities),
        "catalytic_activity_records_json": json.dumps(
            candidate.catalytic_activity_records,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ),
        "cofactors": " | ".join(candidate.cofactors),
        "function_comments": " | ".join(candidate.function_comments),
        "enzyme_system_type": candidate.component_type,
        "auxiliary_requirement_status": "not_required",
        "auxiliary_requirements_json": "[]",
        "warnings": " | ".join(candidate.warnings),
        "reasons": " | ".join(candidate.reasons),
        "rejection_reasons": "",
    }


def _rejected_audit(
    reason: str,
    *,
    accession: str = "",
    evidence_tier: str = "rejected",
) -> dict[str, Any]:
    return {
        "accession": accession,
        "evidence_tier": evidence_tier,
        "evidence_tiers": evidence_tier,
        "fit_status": "rejected",
        "manual_review_required": "true",
        "selection_status": "rejected",
        "direction_verdict": "unknown",
        "rejection_reasons": reason,
    }


def _candidate_conflicts(
    requirement: RetropathEnzymeRequirement,
    candidate: Mapping[str, Any],
) -> str:
    sequence_quality = analyze_protein_sequence(candidate.get("sequence"))
    if not sequence_quality.normalized_sequence:
        return "missing_amino_acid_sequence"
    if sequence_quality.unsupported_positions:
        return sequence_quality.rejection_reason()
    if str(candidate.get("direction_verdict") or "").lower() in {
        "contradicted",
        "unsupported",
        "rejected",
    }:
        return "candidate direction contradicts the biosynthetic route"
    required_ecs = set(requirement.source_ec_numbers)
    candidate_ecs = set(_split(candidate.get("ec_numbers")))
    if required_ecs and candidate_ecs and not (required_ecs & candidate_ecs):
        return "candidate EC contradicts the RR02 source-template EC"
    similarity = _as_float(candidate.get("reaction_similarity"))
    if str(candidate.get("evidence_tier") or "") in {
        "full_reaction_similarity",
        "core_reaction_similarity",
        "rule_similarity",
    } and (similarity is None or not 0.0 < similarity <= 1.0):
        return "structural candidate lacks a positive valid reaction similarity"
    return ""


def _merge_candidates(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    by_accession: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in records:
        accession = str(row.get("accession") or "").strip().upper()
        if accession:
            by_accession[accession].append(row)
    merged: list[dict[str, Any]] = []
    for _accession, values in by_accession.items():
        ordered = sorted(
            values,
            key=lambda row: (
                EVIDENCE_TIER_ORDER.get(str(row.get("evidence_tier") or ""), 9),
                -(_as_float(row.get("reaction_similarity")) or -1.0),
                -(_as_float(row.get("protein_score")) or 0.0),
            ),
        )
        primary = dict(ordered[0])
        primary["evidence_tiers"] = ";".join(dict.fromkeys(
            str(row.get("evidence_tier") or "") for row in ordered
        ))
        primary["retrieval_strategies"] = ";".join(dict.fromkeys(
            value
            for row in ordered
            for value in _split(row.get("retrieval_strategies"))
        ))
        primary["retrieval_query_ids"] = ";".join(dict.fromkeys(
            value
            for row in ordered
            for value in _split(row.get("retrieval_query_ids"))
        ))
        primary["matched_rhea_ids"] = ";".join(sorted({
            value for row in ordered for value in _split(row.get("matched_rhea_ids"))
        }))
        primary["ec_numbers"] = ";".join(sorted({
            value for row in ordered for value in _split(row.get("ec_numbers"))
        }))
        primary["warnings"] = " | ".join(dict.fromkeys(
            value
            for row in ordered
            for value in str(row.get("warnings") or "").split(" | ")
            if value
        ))
        primary["reasons"] = " | ".join(dict.fromkeys(
            value
            for row in ordered
            for value in str(row.get("reasons") or "").split(" | ")
            if value
        ))
        similarities = [
            value
            for row in ordered
            if (value := _as_float(row.get("reaction_similarity"))) is not None
        ]
        if similarities:
            primary["reaction_similarity"] = max(similarities)
        merged.append(primary)
    return merged


def _candidate_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    direction = str(row.get("direction_verdict") or "").lower()
    return (
        EVIDENCE_TIER_ORDER.get(str(row.get("evidence_tier") or ""), 9),
        0 if direction == "supported" else 1,
        -(_as_float(row.get("reaction_similarity")) or -1.0),
        0 if str(row.get("reviewed") or "").lower() == "true" else 1,
        -(_as_float(row.get("protein_score")) or 0.0),
        -(_as_float(row.get("taxonomic_fit_score")) or 0.0),
        str(row.get("accession") or ""),
    )


def _prefetch_source_uniprot_entries(
    requirements: Sequence[RetropathEnzymeRequirement],
    *,
    session: requests.Session,
    entry_cache: dict[str, dict[str, Any] | None],
) -> dict[str, str]:
    accessions = [
        accession
        for requirement in requirements
        for accession in requirement.source_uniprot_ids
    ]
    resolution = resolve_uniprot_accession_batches(
        accessions,
        session=session,
        entry_cache=entry_cache,
        query_id_prefix="uniprot_rr02_source",
    )
    return resolution.accession_errors


def _source_uniprot_candidates(
    requirement: RetropathEnzymeRequirement,
    *,
    chassis_key: str,
    allow_transmembrane: bool,
    session: requests.Session,
    taxonomy_profile: ChassisTaxonomyProfile,
    entry_cache: dict[str, dict[str, Any] | None],
    accession_errors: Mapping[str, str] | None = None,
) -> tuple[list[ProteinCandidate], list[dict[str, Any]]]:
    candidates: list[ProteinCandidate] = []
    evidence: list[dict[str, Any]] = []
    lookup_errors = dict(accession_errors or {})
    if accession_errors is None:
        resolution = resolve_uniprot_accession_batches(
            requirement.source_uniprot_ids,
            session=session,
            entry_cache=entry_cache,
            query_id_prefix="uniprot_rr02_source",
        )
        lookup_errors.update(resolution.accession_errors)
    for accession in requirement.source_uniprot_ids:
        normalized = accession.upper()
        query_id = f"uniprot_source_accession_{normalized}"
        record: dict[str, Any] = {
            "query_id": query_id,
            "query_type": "rr02_source_uniprot",
            "query_value": normalized,
        }
        if normalized in lookup_errors and normalized not in entry_cache:
            record["status"] = "source_unavailable"
            record["error"] = lookup_errors[normalized]
            evidence.append(record)
            continue
        try:
            entry = entry_cache.get(normalized)
            candidate = (
                candidate_from_reaction_entry(
                    entry,
                    chassis_key,
                    retrieval_strategy="rr02_source_uniprot",
                    retrieval_query_id=query_id,
                    matched_rhea_ids=[],
                    allow_transmembrane=allow_transmembrane,
                    function_evidence_reason=(
                        "function: UniProt accession attached to the RR02 "
                        "source template"
                    ),
                    taxonomy_profile=taxonomy_profile,
                )
                if entry is not None
                else None
            )
            if candidate is not None:
                candidates.append(candidate)
                record["status"] = "ok"
            else:
                record["status"] = "no_hit"
        except Exception as exc:
            record["status"] = "source_unavailable"
            record["error"] = str(exc)
        evidence.append(record)
    return candidates, evidence


def _search_requirement(
    requirement: RetropathEnzymeRequirement,
    *,
    chassis_key: str,
    top_n: int,
    max_results: int,
    allow_transmembrane: bool,
    session: requests.Session,
    rhea_client: RheaClient,
    selenzyme_client: SelenzymeClient | None,
    selenzyme_circuit_error: str,
    entry_cache: dict[str, dict[str, Any] | None],
    taxonomy_profile: ChassisTaxonomyProfile | None = None,
    source_accession_errors: Mapping[str, str] | None = None,
) -> _SearchOutcome:
    taxonomy_profile = taxonomy_profile or resolve_chassis_taxonomy(
        chassis_key,
        session=session,
        allow_network=False,
    )
    if not requirement.enzyme_required:
        return _SearchOutcome([], [], [{"status": "not_required"}])
    records: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    source_unavailable = False
    tier = (
        "formal_reaction_exact"
        if requirement.formal_mapping_exact
        else "source_template_supported"
    )
    fit_status = "verified" if requirement.formal_mapping_exact else "manual_review"

    source_candidates, source_evidence = _source_uniprot_candidates(
        requirement,
        chassis_key=chassis_key,
        allow_transmembrane=allow_transmembrane,
        session=session,
        taxonomy_profile=taxonomy_profile,
        entry_cache=entry_cache,
        accession_errors=source_accession_errors,
    )
    evidence.extend(source_evidence)
    source_unavailable |= any(
        row.get("status") == "source_unavailable" for row in source_evidence
    )
    records.extend(
        _candidate_record(candidate, evidence_tier=tier, fit_status=fit_status)
        for candidate in source_candidates
    )

    for ec_number in requirement.source_ec_numbers:
        query_id = f"uniprot_ec_{ec_number}"
        try:
            candidates = recommend_uniprot_proteins(
                ec_number=ec_number,
                chassis_key=chassis_key,
                top_n=max(top_n * 3, 20),
                max_results=max_results,
                allow_transmembrane=allow_transmembrane,
                session=session,
                taxonomy_profile=taxonomy_profile,
            )
            records.extend(
                _candidate_record(
                    candidate,
                    evidence_tier=tier,
                    fit_status=fit_status,
                    query_id=query_id,
                )
                for candidate in candidates
            )
            evidence.append({
                "query_id": query_id,
                "query_type": "uniprot_by_source_ec",
                "query_value": ec_number,
                "status": "ok" if candidates else "no_hit",
                "candidate_count": len(candidates),
            })
        except Exception as exc:
            source_unavailable = True
            evidence.append({
                "query_id": query_id,
                "query_type": "uniprot_by_source_ec",
                "query_value": ec_number,
                "status": "source_unavailable",
                "error": str(exc),
            })

    rhea_ids: set[str] = set(requirement.source_rhea_ids)
    rhea_ids.update(requirement.exact_rhea_ids)
    if requirement.exact_kegg_reaction_ids:
        for reaction_id in requirement.exact_kegg_reaction_ids:
            source = rhea_client.reactions_for_kegg(reaction_id)
            evidence.append({
                "query_id": source.get("query_id", ""),
                "query_type": "rhea_by_exact_kegg_reaction",
                "query_value": reaction_id,
                "status": source.get("status", ""),
                "records": source.get("records", []),
                "error": source.get("error", ""),
            })
            if source.get("status") == "source_unavailable":
                source_unavailable = True
            rhea_ids.update(
                str(row.get("rhea_id") or "")
                for row in source.get("records", [])
                if str(row.get("rhea_id") or "")
            )
    if rhea_ids:
        rhea_requirement = {
            "rhea_retrieval_ids": sorted(rhea_ids),
        }
        candidates, query_ids, errors = retrieve_rhea_candidates_for_requirement(
            rhea_requirement,
            chassis_key,
            max_results=max_results,
            top_n=max(top_n * 3, 20),
            allow_transmembrane=allow_transmembrane,
            session=session,
            taxonomy_profile=taxonomy_profile,
        )
        records.extend(
            _candidate_record(
                candidate,
                evidence_tier=tier,
                fit_status=fit_status,
                query_id=";".join(query_ids),
            )
            for candidate in candidates
        )
        evidence.append({
            "query_ids": query_ids,
            "query_type": "uniprot_by_rhea",
            "query_value": sorted(rhea_ids),
            "status": "source_unavailable" if errors and not candidates else (
                "ok" if candidates else "no_hit"
            ),
            "candidate_count": len(candidates),
            "errors": errors,
        })
        source_unavailable |= bool(errors and not candidates)

    structural_query = ""
    query_kind = ""
    if requirement.step_source == "retropath" and not (
        requirement.formal_mapping_exact and records
    ):
        if requirement.full_reaction_smiles:
            structural_query = requirement.full_reaction_smiles
            query_kind = "full_reaction_smiles"
        elif requirement.core_reaction_smiles:
            structural_query = requirement.core_reaction_smiles
            query_kind = "core_reaction_smiles"
        elif requirement.rule_smarts:
            structural_query = requirement.rule_smarts
            query_kind = "rule_smarts"
    elif requirement.step_source == "kegg_expansion" and not records:
        query_kind = "kegg_reaction"

    if structural_query or query_kind == "kegg_reaction":
        if selenzyme_client is None:
            source_unavailable = True
            evidence.append({
                "query_type": f"selenzyme_by_{query_kind}",
                "query_value": structural_query,
                "status": "source_unavailable",
                "error": selenzyme_circuit_error or "Selenzyme client unavailable",
            })
        else:
            try:
                if query_kind == "kegg_reaction":
                    query_result = selenzyme_client.query_kegg_reaction(
                        requirement.exact_kegg_reaction_ids[0],
                        host_taxon_id=chassis_host_taxon_id(chassis_key),
                        targets=selenzyme_target_count(top_n),
                    )
                    structural_tier = "formal_reaction_exact"
                    structural_fit = "verified_with_risk"
                else:
                    query_result = selenzyme_client.query_reaction_smarts(
                        structural_query,
                        query_kind=query_kind,
                        host_taxon_id=chassis_host_taxon_id(chassis_key),
                        targets=selenzyme_target_count(top_n),
                    )
                    structural_tier = {
                        "full_reaction_smiles": "full_reaction_similarity",
                        "core_reaction_smiles": "core_reaction_similarity",
                        "rule_smarts": "rule_similarity",
                    }[query_kind]
                    structural_fit = "manual_review"
                candidates, raw_audit, query_ids, errors = (
                    retrieve_selenzyme_candidates(
                        {
                            "ec_numbers": list(requirement.source_ec_numbers),
                            "locked_ec_numbers": list(requirement.source_ec_numbers),
                        },
                        query_result,
                        chassis_key,
                        top_n=selenzyme_target_count(top_n),
                        allow_transmembrane=allow_transmembrane,
                        session=session,
                        entry_cache=entry_cache,
                        taxonomy_profile=taxonomy_profile,
                    )
                )
                records.extend(
                    _candidate_record(
                        candidate,
                        evidence_tier=structural_tier,
                        fit_status=structural_fit,
                        query_id=";".join(query_ids)
                        or str(query_result.get("query_id") or ""),
                    )
                    for candidate in candidates
                )
                for raw in raw_audit:
                    if raw.get("gate_status") != "passed":
                        audit.append(_rejected_audit(
                            ";".join(raw.get("rejection_reasons") or []),
                            accession=str(raw.get("accession") or ""),
                        ))
                evidence.append({
                    "query_id": query_result.get("query_id", ""),
                    "query_type": query_result.get("query_type", ""),
                    "query_kind": query_result.get("query_kind", query_kind),
                    "query_sha256": query_result.get("query_sha256", ""),
                    "query_value": query_result.get("query_value", structural_query),
                    "status": query_result.get("status", ""),
                    "app": query_result.get("app", ""),
                    "version": query_result.get("version", ""),
                    "cache_hit": bool(query_result.get("cache_hit")),
                    "response_sha256": query_result.get("response_sha256", ""),
                    "rows": raw_audit,
                    "errors": errors,
                })
            except SelenzymeSourceUnavailable as exc:
                source_unavailable = True
                evidence.append({
                    "query_type": f"selenzyme_by_{query_kind}",
                    "query_value": structural_query,
                    "status": "source_unavailable",
                    "error": str(exc),
                })

    eligible: list[dict[str, Any]] = []
    for row in records:
        conflict = _candidate_conflicts(requirement, row)
        if conflict:
            rejected = dict(row)
            rejected["selection_status"] = "rejected"
            rejected["fit_status"] = "rejected"
            rejected["rejection_reasons"] = conflict
            audit.append(rejected)
        else:
            eligible.append(row)
    merged = _merge_candidates(eligible)
    merged.sort(key=_candidate_sort_key)
    selected: list[dict[str, Any]] = []
    for index, row in enumerate(merged, start=1):
        item = dict(row)
        item["protein_candidate_rank"] = index
        item["selection_status"] = "selected" if index <= top_n else "not_selected"
        audit.append(item)
        if index <= top_n:
            selected.append(item)
    return _SearchOutcome(selected, audit, evidence, source_unavailable)


def _search_key(requirement: RetropathEnzymeRequirement) -> str:
    return _canonical_sha256({
        "step_source": requirement.step_source,
        "reaction_signature_sha256": requirement.reaction_signature_sha256,
        "source_ec_numbers": requirement.source_ec_numbers,
        "source_uniprot_ids": requirement.source_uniprot_ids,
        "source_rhea_ids": requirement.source_rhea_ids,
        "exact_kegg_reaction_ids": requirement.exact_kegg_reaction_ids,
        "exact_rhea_ids": requirement.exact_rhea_ids,
        "formal_mapping_exact": requirement.formal_mapping_exact,
        "full_reaction_smiles": requirement.full_reaction_smiles,
        "core_reaction_smiles": requirement.core_reaction_smiles,
        "rule_smarts": requirement.rule_smarts,
        "enzyme_required": requirement.enzyme_required,
    })


def _bind_candidate(
    requirement: RetropathEnzymeRequirement,
    row: Mapping[str, Any],
) -> dict[str, Any]:
    reaction_id = ";".join(requirement.exact_kegg_reaction_ids)
    return {
        "candidate_rank": requirement.candidate_rank,
        "candidate_id": requirement.candidate_id,
        "combination_rank": "",
        "combination_id": requirement.combination_id,
        "step_index": requirement.step_index,
        "step_id": requirement.step_id,
        "step_source": requirement.step_source,
        "hypothesis_id": requirement.hypothesis_id,
        "reaction_signature_sha256": requirement.reaction_signature_sha256,
        "reaction_id": reaction_id,
        "source_mnxr_id": requirement.source_mnxr_id,
        "rule_id": requirement.rule_id,
        **dict(row),
    }


def requirements_from_manifest(
    requirements: Sequence[Mapping[str, Any]],
    rules: Mapping[str, Mapping[str, str]],
) -> list[RetropathEnzymeRequirement]:
    """Build P9 search identities from already selected manifest steps."""

    result: list[RetropathEnzymeRequirement] = []
    for requirement in requirements:
        if str(requirement.get("step_source") or "").strip() != "retropath":
            continue
        provenance = requirement.get("prediction_provenance")
        if not isinstance(provenance, Mapping):
            raise ValueError("RetroPath manifest step lacks prediction provenance")
        rule_ids = _split(
            requirement.get("retropath_rule_id")
            or provenance.get("rule_id")
            or provenance.get("rule_ids")
        )
        if not rule_ids:
            raise ValueError("RetroPath manifest step has no RR02 rule evidence")
        missing_rules = [rule_id for rule_id in rule_ids if rule_id not in rules]
        if missing_rules:
            raise ValueError(
                "RR02 rule is missing for manifest step: "
                + ", ".join(missing_rules)
            )
        step_id = str(
            requirement.get("retropath_step_id")
            or provenance.get("step_id")
            or ""
        ).strip()
        hypothesis_id = str(
            requirement.get("retropath_hypothesis_id")
            or provenance.get("hypothesis_id")
            or step_id
            or ""
        ).strip()
        reaction_id = str(requirement.get("reaction_id") or "").strip()
        if not step_id or not hypothesis_id or reaction_id not in {
            hypothesis_id,
            step_id,
        }:
            raise ValueError("RetroPath manifest reaction identity is inconsistent")
        candidate_rank = _as_int(
            provenance.get("candidate_rank"),
            "candidate_rank",
            minimum=1,
        )
        candidate_id = str(provenance.get("candidate_id") or "").strip()
        combination_id = str(
            provenance.get("combination_id") or f"raw:{candidate_id}"
        ).strip()
        if not candidate_id:
            raise ValueError("RetroPath manifest route binding is incomplete")
        source_ecs = _complete_ecs(
            requirement.get("source_ec_numbers")
            or provenance.get("source_ec_numbers")
            or []
        )
        exact_kegg = tuple(sorted({
            value.upper()
            for value in _split(
                requirement.get("exact_kegg_reaction_ids")
                or provenance.get("exact_kegg_reaction_ids")
            )
            if re.fullmatch(r"R\d{5}", value.upper())
        }))
        exact_rhea = tuple(sorted(
            _split(
                requirement.get("exact_rhea_ids")
                or provenance.get("exact_rhea_ids")
            ),
            key=lambda value: int(value) if value.isdigit() else 10**12,
        ))
        formal_exact = bool(
            requirement.get("formal_mapping_exact")
            or provenance.get("formal_mapping_exact")
        )
        for rule_id in rule_ids:
            rule = rules[rule_id]
            result.append(RetropathEnzymeRequirement(
                candidate_rank=candidate_rank,
                candidate_id=candidate_id,
                combination_id=combination_id,
                step_index=_as_int(
                    requirement.get("step_index"),
                    "step_index",
                    minimum=1,
                ),
                step_id=step_id,
                step_source="retropath",
                step_status="heterologous",
                hypothesis_id=hypothesis_id,
                reaction_signature_sha256=str(
                    requirement.get("reaction_signature_sha256")
                    or provenance.get("reaction_signature_sha256")
                    or ""
                ),
                full_reaction_smiles=str(
                    requirement.get("full_reaction_smiles")
                    or provenance.get("full_reaction_smiles")
                    or ""
                ),
                core_reaction_smiles=str(
                    requirement.get("core_reaction_smiles")
                    or provenance.get("core_reaction_smiles")
                    or ""
                ),
                rule_id=rule_id,
                rule_smarts=str(rule.get("Rule") or "").strip(),
                source_mnxr_id=str(
                    requirement.get("source_mnxr_id")
                    or provenance.get("source_mnxr_id")
                    or ""
                ),
                source_reaction_ids=_split(
                    requirement.get("source_reaction_ids")
                    or provenance.get("source_reaction_ids")
                ),
                source_ec_numbers=source_ecs,
                source_uniprot_ids=_split(
                    requirement.get("source_uniprot_ids")
                    or provenance.get("source_uniprot_ids")
                ),
                source_rhea_ids=_split(
                    requirement.get("source_rhea_ids")
                    or provenance.get("source_rhea_ids")
                ),
                exact_kegg_reaction_ids=exact_kegg if formal_exact else tuple(),
                exact_rhea_ids=exact_rhea if formal_exact else tuple(),
                formal_mapping_exact=formal_exact,
                enzyme_required=True,
            ))
    return result


def _standard_candidate_row(
    requirement: Mapping[str, Any],
    search_requirement: RetropathEnzymeRequirement,
    row: Mapping[str, Any],
) -> dict[str, Any]:
    fit_status = str(row.get("fit_status") or "manual_review").strip()
    similarity = _as_float(row.get("reaction_similarity"))
    reaction_fit_score = (
        1.0
        if fit_status == "verified"
        else (similarity if similarity is not None and similarity > 0 else 0.0)
    )
    evidence_tier = str(row.get("evidence_tier") or "").strip()
    source_ids = ";".join(dict.fromkeys((
        *search_requirement.source_reaction_ids,
        search_requirement.source_mnxr_id,
    )))
    return {
        "solution_id": requirement.get("solution_id"),
        "step_index": requirement.get("step_index"),
        "reaction_id": requirement.get("reaction_id"),
        "reaction_name": requirement.get("reaction_name"),
        "produced_compound_id": requirement.get("produced_compound_id"),
        "produced_compound_name": requirement.get("produced_compound_name"),
        "role": "main",
        "candidate_role": row.get("candidate_role") or "catalytic_main",
        "ec_number": next(iter(_split(row.get("ec_numbers"))), ""),
        "ec_numbers": row.get("ec_numbers", ""),
        "accession": row.get("accession", ""),
        "entry_name": row.get("entry_name", ""),
        "protein_name": row.get("protein_name", ""),
        "organism_name": row.get("organism_name", ""),
        "organism_id": row.get("organism_id", ""),
        "taxonomic_lineage": row.get("taxonomic_lineage", ""),
        "taxonomic_lineage_ids": row.get("taxonomic_lineage_ids", ""),
        "taxonomic_shared_taxon_id": row.get(
            "taxonomic_shared_taxon_id", ""
        ),
        "taxonomic_shared_name": row.get("taxonomic_shared_name", ""),
        "taxonomic_shared_rank": row.get("taxonomic_shared_rank", ""),
        "taxonomic_fit_status": row.get("taxonomic_fit_status", "unknown"),
        "taxonomic_fit_score": row.get("taxonomic_fit_score", 50.0),
        "taxonomy_evidence_source": row.get("taxonomy_evidence_source", ""),
        "reviewed": row.get("reviewed", "false"),
        "length": row.get("length", ""),
        "score": row.get("protein_score", 0.0),
        "score_breakdown": row.get("score_breakdown", "{}"),
        "evaluation_rank": row.get("protein_candidate_rank", ""),
        "candidate_rank": row.get("protein_candidate_rank", ""),
        "selection_status": row.get("selection_status", ""),
        "retrieval_strategy": row.get("retrieval_strategies", ""),
        "retrieval_query_id": row.get("retrieval_query_ids", ""),
        "matched_rhea_ids": row.get("matched_rhea_ids", ""),
        "matched_ko_ids": "",
        "direction_support": row.get("direction_verdict", "unknown"),
        "direction_verdict": row.get("direction_verdict", "unknown") or "unknown",
        "direction_confidence": row.get("direction_confidence", "") or "low",
        "reaction_confidence": evidence_tier,
        "publication_ids": row.get("publication_ids", ""),
        "catalytic_activities": row.get("catalytic_activities", ""),
        "catalytic_activity_records_json": row.get(
            "catalytic_activity_records_json", "[]"
        ),
        "cofactors": row.get("cofactors", ""),
        "function_comments": row.get("function_comments", ""),
        "sequence_version": row.get("sequence_version", ""),
        "sequence_sha256": row.get("sequence_sha256", ""),
        "sequence": row.get("sequence", ""),
        "warnings": row.get("warnings", ""),
        "reasons": " | ".join(filter(None, (
            str(row.get("reasons") or ""),
            f"retropath_evidence_tier:{evidence_tier}",
            f"retropath_source_reactions:{source_ids}" if source_ids else "",
        ))),
        "rejection_reasons": row.get("rejection_reasons", ""),
        "reaction_fit_status": fit_status,
        "reaction_fit_score": reaction_fit_score,
        "reaction_fit_rule_ids": search_requirement.rule_id,
        "reaction_fit_evidence": (
            f"RetroPath evidence tier {evidence_tier}; "
            f"template {search_requirement.source_mnxr_id}"
        ),
        "specificity_status": (
            "supported" if search_requirement.formal_mapping_exact else "unknown"
        ),
        "enzyme_system_type": row.get("enzyme_system_type", ""),
        "auxiliary_requirement_status": row.get(
            "auxiliary_requirement_status", "not_required"
        ),
        "auxiliary_requirements_json": row.get(
            "auxiliary_requirements_json", "[]"
        ),
        "selenzyme_reaction_similarity": row.get("reaction_similarity", ""),
        "selenzyme_sim_rf": row.get("sim_rf", ""),
        "selenzyme_sim_2018": row.get("sim_2018", ""),
        "selenzyme_matched_reaction_id": row.get("matched_reaction_id", ""),
        "selenzyme_taxonomic_distance": row.get("taxonomic_distance", ""),
        "source_mnxr_id": search_requirement.source_mnxr_id,
        "retropath_rule_id": search_requirement.rule_id,
        "retropath_hypothesis_id": search_requirement.hypothesis_id,
        "retropath_evidence_tier": evidence_tier,
        "manual_review_required": row.get("manual_review_required", "true"),
    }


def retrieve_manifest_retropath_candidates(
    requirements: Sequence[Mapping[str, Any]],
    *,
    config: Any,
    top_n: int,
    max_results: int,
    allow_transmembrane: bool,
    session: requests.Session,
) -> dict[str, Any]:
    """Retrieve predicted-step candidates and project standard candidate rows."""

    rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    rules = _load_rules(rules_path)
    search_requirements = requirements_from_manifest(requirements, rules)
    if not search_requirements:
        return {
            "selected_rows": [],
            "audit_rows": [],
            "requirements": [],
            "evidence": [],
            "source_unavailable": False,
        }
    by_step = {
        _as_int(item.get("step_index"), "step_index", minimum=1): item
        for item in requirements
        if str(item.get("step_source") or "").strip() == "retropath"
    }
    chassis_key = str(getattr(config, "chassis_key", "ecoli_mg1655") or "ecoli_mg1655")
    cache_root = Path(config.cache_dir).resolve() / "main_protein_selection"
    taxonomy_profile = getattr(config, "taxonomy_profile", None)
    if not isinstance(taxonomy_profile, ChassisTaxonomyProfile):
        taxonomy_profile = resolve_chassis_taxonomy(
            chassis_key,
            session=session,
            cache_root=cache_root,
        )
    client_error = ""
    try:
        selenzyme_client: SelenzymeClient | None = SelenzymeClient(
            session=session,
            cache_root=cache_root / "selenzyme",
        )
    except Exception as exc:
        client_error = str(exc)
        selenzyme_client = None
    rhea_client = RheaClient(session=session, cache_root=cache_root / "rhea")
    entry_cache: dict[str, dict[str, Any] | None] = {}
    source_accession_errors = _prefetch_source_uniprot_entries(
        search_requirements,
        session=session,
        entry_cache=entry_cache,
    )
    search_cache: dict[str, _SearchOutcome] = {}
    selected: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    evidence: list[dict[str, Any]] = []
    source_unavailable = False
    for search_requirement in search_requirements:
        key = _search_key(search_requirement)
        if key not in search_cache:
            search_cache[key] = _search_requirement(
                search_requirement,
                chassis_key=chassis_key,
                top_n=top_n,
                max_results=max_results,
                allow_transmembrane=allow_transmembrane,
                session=session,
                rhea_client=rhea_client,
                selenzyme_client=selenzyme_client,
                selenzyme_circuit_error=client_error,
                entry_cache=entry_cache,
                taxonomy_profile=taxonomy_profile,
                source_accession_errors=source_accession_errors,
            )
        outcome = search_cache[key]
        source_unavailable |= outcome.source_unavailable
        if outcome.source_unavailable and any(
            row.get("status") == "source_unavailable"
            and str(row.get("query_type") or "").startswith("selenzyme")
            for row in outcome.evidence
        ):
            selenzyme_client = None
        base_requirement = by_step[search_requirement.step_index]
        selected.extend(
            _standard_candidate_row(base_requirement, search_requirement, row)
            for row in outcome.candidates
        )
        audit.extend(
            _standard_candidate_row(base_requirement, search_requirement, row)
            for row in outcome.audit
        )
        evidence.extend(
            {"search_key": key, "step_index": search_requirement.step_index, **row}
            for row in outcome.evidence
        )
    selected.sort(key=lambda row: (
        int(row.get("step_index") or 0),
        int(row.get("candidate_rank") or 0),
        str(row.get("accession") or ""),
    ))
    audit.sort(key=lambda row: (
        int(row.get("step_index") or 0),
        int(row.get("evaluation_rank") or 0),
        str(row.get("accession") or ""),
    ))
    return {
        "selected_rows": selected,
        "audit_rows": audit,
        "requirements": [item.to_dict() for item in search_requirements],
        "evidence": evidence,
        "source_unavailable": source_unavailable,
    }


def _combination_summaries(
    requirements: Sequence[RetropathEnzymeRequirement],
    selected: Sequence[Mapping[str, Any]],
    unavailable_keys: set[str],
) -> list[dict[str, Any]]:
    requirements_by_combination: dict[
        str,
        list[RetropathEnzymeRequirement],
    ] = defaultdict(list)
    rows_by_scope: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for requirement in requirements:
        requirements_by_combination[requirement.combination_id].append(requirement)
    for row in selected:
        rows_by_scope[(str(row["combination_id"]), str(row["step_id"]))].append(row)

    summaries: list[dict[str, Any]] = []
    for combination_id, combination_requirements in requirements_by_combination.items():
        # One raw RetroPath step can carry several RR02 rules.  They are
        # alternative search evidence for the same enzyme requirement, not
        # several independently required enzymes.
        required_by_step: dict[str, RetropathEnzymeRequirement] = {}
        for item in combination_requirements:
            if item.enzyme_required:
                required_by_step.setdefault(item.step_id, item)
        required = list(required_by_step.values())
        missing = [
            item.step_id
            for item in required
            if not rows_by_scope.get((combination_id, item.step_id))
        ]
        unavailable = [
            item.step_id
            for item in required
            if item.step_id in missing and _search_key(item) in unavailable_keys
        ]
        best_rows = [
            sorted(
                rows_by_scope[(combination_id, item.step_id)],
                key=_candidate_sort_key,
            )[0]
            for item in required
            if rows_by_scope.get((combination_id, item.step_id))
        ]
        if not missing:
            status = "ready_for_review"
        elif unavailable:
            status = "source_unavailable"
        else:
            status = "partial_no_candidate"
        worst_tier = max(
            (
                EVIDENCE_TIER_ORDER.get(str(row.get("evidence_tier") or ""), 9)
                for row in best_rows
            ),
            default=9,
        )
        similarities = [
            value
            for row in best_rows
            if (value := _as_float(row.get("reaction_similarity"))) is not None
        ]
        summaries.append({
            "combination_id": combination_id,
            "combination_rank": 0,
            "status": status,
            "required_step_count": len(required),
            "covered_step_count": len(required) - len(missing),
            "uncovered_step_ids": sorted(missing),
            "source_unavailable_step_ids": sorted(unavailable),
            "worst_evidence_tier_order": worst_tier,
            "formal_exact_step_count": sum(
                str(row.get("evidence_tier")) == "formal_reaction_exact"
                for row in best_rows
            ),
            "minimum_structural_similarity": (
                min(similarities) if similarities else None
            ),
        })
    status_order = {
        "ready_for_review": 0,
        "partial_no_candidate": 1,
        "source_unavailable": 2,
    }
    summaries.sort(key=lambda row: (
        status_order[row["status"]],
        row["worst_evidence_tier_order"],
        -row["formal_exact_step_count"],
        -(
            row["minimum_structural_similarity"]
            if row["minimum_structural_similarity"] is not None
            else 1.0
        ),
        row["combination_id"],
    ))
    for rank, row in enumerate(summaries, start=1):
        row["combination_rank"] = rank
    return summaries


def select_retropath_enzymes(
    config: Any,
    *,
    top_n: int = 5,
    max_results: int = 1000,
    allow_transmembrane: bool = False,
    session: requests.Session | None = None,
    selenzyme_client: SelenzymeClient | None = None,
) -> dict[str, Any]:
    """Generate P9 candidates without promoting a predicted route."""

    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    inputs = _load_p9_inputs(config)
    rules_path = Path(config.retropath_rules_path).expanduser().resolve()
    p8_rules_sha = str(inputs.p8_manifest.get("inputs", {}).get("rr02_sha256") or "")
    if not rules_path.is_file() or _sha256_file(rules_path) != p8_rules_sha:
        raise ValueError("RR02 rules differ from the P8-validated rules")
    rules = _load_rules(rules_path)
    mnxref_dir = Path(config.data_dir) / "retropath" / "mnxref" / "3.0"
    with MnxrefIndex(mnxref_dir, rules_path) as mnxref:
        requirements = _build_requirements(inputs, mnxref, rules)
        mnxref_manifest = mnxref.manifest

    http = session or requests.Session()
    chassis_key = str(getattr(config, "chassis_key", "ecoli_mg1655") or "ecoli_mg1655")
    cache_root = Path(config.cache_dir).resolve() / "main_protein_selection"
    taxonomy_profile = resolve_chassis_taxonomy(
        chassis_key,
        session=http,
        cache_root=cache_root,
    )
    client_error = ""
    if selenzyme_client is None:
        try:
            selenzyme_client = SelenzymeClient(
                session=http,
                cache_root=cache_root / "selenzyme",
            )
        except Exception as exc:
            client_error = str(exc)
            selenzyme_client = None
    rhea_client = RheaClient(session=http, cache_root=cache_root / "rhea")
    entry_cache: dict[str, dict[str, Any] | None] = {}
    source_accession_errors = _prefetch_source_uniprot_entries(
        requirements,
        session=http,
        entry_cache=entry_cache,
    )
    search_cache: dict[str, _SearchOutcome] = {}
    selected_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    query_evidence: list[dict[str, Any]] = []
    unavailable_keys: set[str] = set()

    for requirement in requirements:
        if not requirement.enzyme_required:
            continue
        key = _search_key(requirement)
        if key not in search_cache:
            outcome = _search_requirement(
                requirement,
                chassis_key=chassis_key,
                top_n=top_n,
                max_results=max_results,
                allow_transmembrane=allow_transmembrane,
                session=http,
                rhea_client=rhea_client,
                selenzyme_client=selenzyme_client,
                selenzyme_circuit_error=client_error,
                entry_cache=entry_cache,
                taxonomy_profile=taxonomy_profile,
                source_accession_errors=source_accession_errors,
            )
            search_cache[key] = outcome
            if outcome.source_unavailable:
                unavailable_keys.add(key)
                if any(
                    row.get("status") == "source_unavailable"
                    and str(row.get("query_type") or "").startswith("selenzyme")
                    for row in outcome.evidence
                ):
                    client_error = next(
                        (
                            str(row.get("error") or "")
                            for row in outcome.evidence
                            if row.get("status") == "source_unavailable"
                            and str(row.get("query_type") or "").startswith("selenzyme")
                        ),
                        client_error,
                    )
                    selenzyme_client = None
        outcome = search_cache[key]
        selected_rows.extend(
            _bind_candidate(requirement, row) for row in outcome.candidates
        )
        audit_rows.extend(_bind_candidate(requirement, row) for row in outcome.audit)
        if not any(item.get("search_key") == key for item in query_evidence):
            query_evidence.extend(
                {"search_key": key, **row} for row in outcome.evidence
            )

    combinations = _combination_summaries(
        requirements,
        selected_rows,
        unavailable_keys,
    )
    combination_ranks = {
        row["combination_id"]: row["combination_rank"] for row in combinations
    }
    for row in [*selected_rows, *audit_rows]:
        row["combination_rank"] = combination_ranks.get(row["combination_id"], "")
    selected_rows.sort(key=lambda row: (
        _as_int(row.get("combination_rank"), "combination_rank", minimum=1),
        _as_int(row.get("step_index"), "step_index", minimum=1),
        _as_int(row.get("protein_candidate_rank"), "protein_candidate_rank", minimum=1),
    ))
    audit_rows.sort(key=lambda row: (
        _as_int(row.get("combination_rank"), "combination_rank", minimum=1),
        _as_int(row.get("step_index"), "step_index", minimum=1),
        str(row.get("accession") or ""),
        str(row.get("selection_status") or ""),
    ))
    if any(row["status"] == "ready_for_review" for row in combinations):
        overall_status = "ready_for_review"
    elif any(row["status"] == "source_unavailable" for row in combinations):
        overall_status = "source_unavailable"
    else:
        overall_status = "partial_no_candidate"

    output_dir = (
        inputs.candidate_inputs.retropath_dir
        / ENZYME_SELECTION_DIR_NAME
        / f"candidate_{inputs.candidate_rank}"
    )
    requirements_path = output_dir / ENZYME_REQUIREMENTS_FILE_NAME
    candidate_path = output_dir / STEP_ENZYME_CANDIDATES_FILE_NAME
    audit_path = output_dir / STEP_ENZYME_AUDIT_FILE_NAME
    evidence_path = output_dir / SELENZYME_EVIDENCE_FILE_NAME
    selection_path = output_dir / ENZYME_SELECTION_FILE_NAME

    requirements_payload = {
        "schema_version": RETROPATH_ENZYME_REQUIREMENTS_SCHEMA,
        "candidate_rank": inputs.candidate_rank,
        "candidate_id": inputs.candidate_id,
        "requirements": [item.to_dict() for item in requirements],
    }
    evidence_payload = {
        "schema_version": RETROPATH_SELENZYME_EVIDENCE_SCHEMA,
        "policy": {
            "structural_matches_require_manual_review": True,
            "positive_similarity_required": True,
            "arbitrary_similarity_threshold": None,
            "query_precedence": [
                "full_reaction_smiles",
                "core_reaction_smiles",
                "rule_smarts",
            ],
        },
        "queries": query_evidence,
    }
    artifact_hashes = {
        ENZYME_REQUIREMENTS_FILE_NAME: _atomic_write_json(
            requirements_path,
            requirements_payload,
        ),
        STEP_ENZYME_CANDIDATES_FILE_NAME: _atomic_write_text(
            candidate_path,
            _render_csv(ENZYME_CANDIDATE_COLUMNS, selected_rows),
        ),
        STEP_ENZYME_AUDIT_FILE_NAME: _atomic_write_text(
            audit_path,
            _render_csv(ENZYME_CANDIDATE_COLUMNS, audit_rows),
        ),
        SELENZYME_EVIDENCE_FILE_NAME: _atomic_write_json(
            evidence_path,
            evidence_payload,
        ),
    }
    artifact_paths = {
        ENZYME_REQUIREMENTS_FILE_NAME: requirements_path,
        STEP_ENZYME_CANDIDATES_FILE_NAME: candidate_path,
        STEP_ENZYME_AUDIT_FILE_NAME: audit_path,
        SELENZYME_EVIDENCE_FILE_NAME: evidence_path,
    }
    selection: dict[str, Any] = {
        "schema_version": RETROPATH_ENZYME_SELECTION_SCHEMA,
        "created_at": datetime.now(UTC).isoformat(),
        "target_compound": inputs.candidate_inputs.target_compound,
        "expansion_depth": inputs.candidate_inputs.depth,
        "candidate_rank": inputs.candidate_rank,
        "candidate_id": inputs.candidate_id,
        "status": overall_status,
        "recommended_combination_id": (
            combinations[0]["combination_id"] if combinations else ""
        ),
        "combination_count": len(combinations),
        "combinations": combinations,
        "parameters": {
            "top_n": top_n,
            "max_results": max_results,
            "allow_transmembrane": allow_transmembrane,
            "chassis_key": chassis_key,
        },
        "inputs": {
            "pipeline_result": {
                "path": str(inputs.candidate_inputs.pipeline_path),
                "sha256": _sha256_file(inputs.candidate_inputs.pipeline_path),
            },
            "p8_validation_manifest": {
                "path": str(inputs.p8_manifest_path),
                "sha256": _sha256_file(inputs.p8_manifest_path),
            },
            "candidate_routes_sha256": inputs.candidate_inputs.pipeline["artifacts"][
                "candidate_routes"
            ]["sha256"],
            "candidate_steps_sha256": inputs.candidate_inputs.pipeline["artifacts"][
                "candidate_steps"
            ]["sha256"],
            "rr02_path": str(rules_path),
            "rr02_sha256": _sha256_file(rules_path),
            "mnxref_index_path": mnxref_manifest["index_path"],
            "mnxref_index_sha256": mnxref_manifest["index_sha256"],
        },
        "artifacts": {
            name: {
                "path": str(artifact_paths[name].resolve()),
                "sha256": sha256,
            }
            for name, sha256 in artifact_hashes.items()
        },
        "review_required": True,
        "formal_promotion_allowed": False,
    }
    selection_sha256 = _atomic_write_json(selection_path, selection)
    return {
        "ok": overall_status == "ready_for_review",
        "status": overall_status,
        "search_engine": "retropath",
        "target_compound": inputs.candidate_inputs.target_compound,
        "expansion_depth": inputs.candidate_inputs.depth,
        "candidate_rank": inputs.candidate_rank,
        "candidate_id": inputs.candidate_id,
        "recommended_combination_id": selection["recommended_combination_id"],
        "combination_count": len(combinations),
        "step_candidate_count": len(selected_rows),
        "output_dir": str(output_dir.resolve()),
        "selection_manifest": str(selection_path.resolve()),
        "selection_manifest_sha256": selection_sha256,
        "review_required": True,
        "formal_promotion_allowed": False,
    }


def run_retropath_enzyme_selection(config: Any) -> dict[str, Any]:
    if bool(getattr(config, "literature_search", False)):
        raise ValueError(
            "--literature-search is not supported for RetroPath P9; "
            "use the audited source-template and SelenzymeRF evidence"
        )
    result = select_retropath_enzymes(
        config,
        top_n=int(getattr(config, "top_n", 5)),
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "ENZYME_CANDIDATE_COLUMNS",
    "ENZYME_REQUIREMENTS_FILE_NAME",
    "ENZYME_SELECTION_DIR_NAME",
    "ENZYME_SELECTION_FILE_NAME",
    "RETROPATH_ENZYME_REQUIREMENTS_SCHEMA",
    "RETROPATH_ENZYME_SELECTION_SCHEMA",
    "RETROPATH_SELENZYME_EVIDENCE_SCHEMA",
    "RetropathEnzymeRequirement",
    "SELENZYME_EVIDENCE_FILE_NAME",
    "STEP_ENZYME_AUDIT_FILE_NAME",
    "STEP_ENZYME_CANDIDATES_FILE_NAME",
    "run_retropath_enzyme_selection",
    "requirements_from_manifest",
    "retrieve_manifest_retropath_candidates",
    "select_retropath_enzymes",
]
