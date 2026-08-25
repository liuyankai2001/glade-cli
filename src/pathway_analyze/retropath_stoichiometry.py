"""Conservative source-template stoichiometry reconstruction for RP2 steps."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdMolDescriptors

from src.pathway_analyze.retropath_mnxref import (
    MnxrefChemical,
    MnxrefIndex,
    MnxrefReactionTemplate,
)


STOICHIOMETRY_SCHEMA_VERSION = "retropath_stoichiometry.v1"
DEFAULT_MAX_STEP_HYPOTHESES = 8
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")
_ELEMENT_SYMBOLS = frozenset(
    Chem.GetPeriodicTable().GetElementSymbol(atomic_number)
    for atomic_number in range(1, 119)
)
_RP2_INCHIKEY_ID = re.compile(
    r"^RP2CPD:(?P<inchikey>[A-Z]{14}-[A-Z]{10}-[A-Z])$"
)


@dataclass(frozen=True)
class CompoundProperty:
    compound_id: str
    name: str
    formula: str
    charge: int
    inchi: str = ""
    inchikey: str = ""
    smiles: str = ""
    source: str = ""
    source_mnxm_id: str = ""
    xrefs: tuple[str, ...] = tuple()

    def to_dict(self) -> dict[str, Any]:
        return {
            "compound_id": self.compound_id,
            "name": self.name,
            "formula": self.formula,
            "charge": self.charge,
            "inchi": self.inchi,
            "inchikey": self.inchikey,
            "smiles": self.smiles,
            "source": self.source,
            "source_mnxm_id": self.source_mnxm_id,
            "xrefs": list(self.xrefs),
        }


@dataclass(frozen=True)
class CompletedTerm:
    side: str
    coefficient: float
    compound: CompoundProperty
    role: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "side": self.side,
            "coefficient": self.coefficient,
            "role": self.role,
            **self.compound.to_dict(),
        }


@dataclass(frozen=True)
class CompletedReactionHypothesis:
    hypothesis_id: str
    candidate_id: str
    step_id: str
    rule_id: str
    source_mnxr_id: str
    source_reference: str
    source_equation: str
    source_orientation: str
    evidence_grade: str
    terms: tuple[CompletedTerm, ...]
    balance_status: str
    cofactor_reconstruction_status: str

    @property
    def recovered_compound_ids(self) -> tuple[str, ...]:
        return tuple(
            sorted(
                {
                    term.compound.compound_id
                    for term in self.terms
                    if term.role == "recovered_template_participant"
                }
            )
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": STOICHIOMETRY_SCHEMA_VERSION,
            "hypothesis_id": self.hypothesis_id,
            "candidate_id": self.candidate_id,
            "step_id": self.step_id,
            "rule_id": self.rule_id,
            "source_mnxr_id": self.source_mnxr_id,
            "source_reference": self.source_reference,
            "source_equation": self.source_equation,
            "source_orientation": self.source_orientation,
            "evidence_grade": self.evidence_grade,
            "terms": [term.to_dict() for term in self.terms],
            "balance_status": self.balance_status,
            "cofactor_reconstruction_status": (
                self.cofactor_reconstruction_status
            ),
            "recovered_compound_ids": list(self.recovered_compound_ids),
        }


@dataclass(frozen=True)
class ReconstructionRejection:
    candidate_id: str
    step_id: str
    rule_id: str
    source_mnxr_id: str
    reason_code: str
    reason_detail: str

    def to_dict(self) -> dict[str, str]:
        return {
            "candidate_id": self.candidate_id,
            "step_id": self.step_id,
            "rule_id": self.rule_id,
            "source_mnxr_id": self.source_mnxr_id,
            "reason_code": self.reason_code,
            "reason_detail": self.reason_detail,
        }


@dataclass(frozen=True)
class StepReconstructionResult:
    candidate_id: str
    step_id: str
    status: str
    hypotheses: tuple[CompletedReactionHypothesis, ...]
    rejections: tuple[ReconstructionRejection, ...]
    truncated: bool


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def parse_formula(formula: str) -> dict[str, float]:
    normalized = str(formula or "").strip()
    if not normalized or "*" in normalized:
        raise ValueError("formula is missing or incomplete")
    position = 0
    elements: dict[str, float] = defaultdict(float)
    for match in _FORMULA_TOKEN.finditer(normalized):
        if match.start() != position:
            raise ValueError(f"unsupported formula syntax: {normalized}")
        amount = float(match.group(2) or "1")
        if match.group(1) not in _ELEMENT_SYMBOLS:
            raise ValueError(f"unsupported formula element: {match.group(1)}")
        if not math.isfinite(amount) or amount <= 0:
            raise ValueError(f"invalid formula coefficient: {normalized}")
        elements[match.group(1)] += amount
        position = match.end()
    if position != len(normalized) or not elements:
        raise ValueError(f"unsupported formula syntax: {normalized}")
    return dict(elements)


def _balance_vector(
    terms: Iterable[CompletedTerm],
) -> tuple[tuple[tuple[str, float], ...], float]:
    elements: dict[str, float] = defaultdict(float)
    charge = 0.0
    for term in terms:
        sign = -1.0 if term.side == "left" else 1.0
        for element, amount in parse_formula(term.compound.formula).items():
            elements[element] += sign * term.coefficient * amount
        charge += sign * term.coefficient * term.compound.charge
    normalized = tuple(
        sorted(
            (element, round(amount, 9))
            for element, amount in elements.items()
            if abs(amount) > 1e-8
        )
    )
    return normalized, round(charge, 9)


def is_balanced(terms: Iterable[CompletedTerm]) -> bool:
    elements, charge = _balance_vector(terms)
    return not elements and abs(charge) <= 1e-8


def _compound_from_molecule(
    compound_id: str,
    molecule: Chem.Mol,
    *,
    source: str,
) -> CompoundProperty:
    try:
        inchi = rd_inchi.MolToInchi(molecule).strip()
        inchikey = rd_inchi.InchiToInchiKey(inchi).strip().upper()
        formula = rdMolDescriptors.CalcMolFormula(molecule).strip()
        smiles = Chem.MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=True,
        ).strip()
    except (RuntimeError, ValueError) as exc:
        raise ValueError(f"cannot standardize structure for {compound_id}") from exc
    if not inchi.startswith("InChI=1S/") or not inchikey or not formula:
        raise ValueError(f"incomplete structure for {compound_id}")
    charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    return CompoundProperty(
        compound_id=compound_id,
        name=compound_id,
        formula=formula,
        charge=int(charge),
        inchi=inchi,
        inchikey=inchikey,
        smiles=smiles,
        source=source,
    )


def _compound_from_mnxref(chemical: MnxrefChemical) -> CompoundProperty:
    if chemical.charge is None:
        raise ValueError(f"MNXref chemical has no charge: {chemical.mnxm_id}")
    parse_formula(chemical.formula)
    return CompoundProperty(
        compound_id=chemical.mnxm_id,
        name=chemical.name,
        formula=chemical.formula,
        charge=chemical.charge,
        inchi=chemical.inchi,
        inchikey=chemical.inchikey,
        smiles=chemical.smiles,
        source="mnxref:3.0",
        source_mnxm_id=chemical.mnxm_id,
        xrefs=chemical.xrefs,
    )


def load_p2_compound_properties(path: str | Path) -> dict[str, CompoundProperty]:
    import csv

    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"P2 compound mapping not found: {resolved}")
    result: dict[str, CompoundProperty] = {}
    with resolved.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"kegg_id", "inchi", "inchikey", "formula", "charge"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"P2 compound mapping is missing columns: {missing}")
        for row in reader:
            compound_id = str(row.get("kegg_id") or "").strip().upper()
            formula = str(row.get("formula") or "").strip()
            charge_text = str(row.get("charge") or "").strip()
            if not compound_id or not formula or not charge_text:
                continue
            parse_formula(formula)
            result[compound_id] = CompoundProperty(
                compound_id=compound_id,
                name=compound_id,
                formula=formula,
                charge=int(charge_text),
                inchi=str(row.get("inchi") or "").strip(),
                inchikey=str(row.get("inchikey") or "").strip().upper(),
                smiles=str(row.get("isomeric_smiles") or "").strip(),
                source="p2_compound_mapping",
            )
    return result


def _parse_stoichiometry(value: Any, field_name: str) -> tuple[tuple[str, float], ...]:
    try:
        payload = json.loads(str(value or ""))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{field_name} is not valid JSON") from exc
    if not isinstance(payload, list) or not payload:
        raise ValueError(f"{field_name} must be a non-empty list")
    rows: list[tuple[str, float]] = []
    for item in payload:
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError(f"{field_name} contains an invalid item")
        compound_id = str(item[0] or "").strip()
        amount = float(item[1])
        if not compound_id or not math.isfinite(amount) or amount <= 0:
            raise ValueError(f"{field_name} contains an invalid item")
        rows.append((compound_id, amount))
    return tuple(rows)


def _reaction_smiles_sides(value: str) -> tuple[list[Chem.Mol], list[Chem.Mol]]:
    if str(value).count(">>") != 1:
        raise ValueError("Reaction SMILES must contain exactly one >>")
    raw_left, raw_right = str(value).split(">>", 1)

    def parse_side(side: str) -> list[Chem.Mol]:
        molecules: list[Chem.Mol] = []
        for component in (item for item in side.split(".") if item.strip()):
            molecule = Chem.MolFromSmiles(component.strip(), sanitize=True)
            if molecule is None:
                raise ValueError("Reaction SMILES contains an invalid component")
            molecules.append(molecule)
        if not molecules:
            raise ValueError("Reaction SMILES contains an empty side")
        return molecules

    return parse_side(raw_left), parse_side(raw_right)


def _assign_side_properties(
    stoichiometry: Sequence[tuple[str, float]],
    molecules: Sequence[Chem.Mol],
    known: Mapping[str, CompoundProperty],
) -> tuple[CompletedTerm, ...]:
    molecule_properties = [
        _compound_from_molecule(
            f"structure:{index}",
            molecule,
            source="retropath_reaction_smiles",
        )
        for index, molecule in enumerate(molecules)
    ]
    by_key: dict[str, list[CompoundProperty]] = defaultdict(list)
    for item in molecule_properties:
        by_key[item.inchikey].append(item)
    terms: list[CompletedTerm] = []
    for compound_id, coefficient in stoichiometry:
        property_value = known.get(compound_id)
        expected_key = property_value.inchikey if property_value is not None else ""
        match = _RP2_INCHIKEY_ID.fullmatch(compound_id)
        if match is not None:
            expected_key = match.group("inchikey")
        matching = by_key.get(expected_key, []) if expected_key else []
        if not matching:
            raise ValueError(
                f"cannot map {compound_id} to a Reaction SMILES component"
            )
        structure = matching.pop(0)
        if property_value is None:
            property_value = CompoundProperty(
                compound_id=compound_id,
                name=compound_id,
                formula=structure.formula,
                charge=structure.charge,
                inchi=structure.inchi,
                inchikey=structure.inchikey,
                smiles=structure.smiles,
                source="retropath_reaction_smiles",
            )
        elif (
            property_value.formula != structure.formula
            or property_value.charge != structure.charge
        ):
            raise ValueError(
                f"P2 and Reaction SMILES structures disagree for {compound_id}"
            )
        terms.append(
            CompletedTerm(
                side="",
                coefficient=coefficient,
                compound=property_value,
                role="predicted_core",
            )
        )
    if any(values for values in by_key.values()):
        raise ValueError("Reaction SMILES has structural components absent from P5")
    return tuple(terms)


def core_terms_from_step(
    step_row: Mapping[str, Any],
    known_properties: Mapping[str, CompoundProperty],
) -> tuple[tuple[CompletedTerm, ...], tuple[CompletedTerm, ...]]:
    left_stoichiometry = _parse_stoichiometry(
        step_row.get("substrate_stoichiometry_json"),
        "substrate_stoichiometry_json",
    )
    right_stoichiometry = _parse_stoichiometry(
        step_row.get("product_stoichiometry_json"),
        "product_stoichiometry_json",
    )
    left_molecules, right_molecules = _reaction_smiles_sides(
        str(step_row.get("reaction_smiles") or "")
    )
    left = tuple(
        CompletedTerm("left", item.coefficient, item.compound, item.role)
        for item in _assign_side_properties(
            left_stoichiometry,
            left_molecules,
            known_properties,
        )
    )
    right = tuple(
        CompletedTerm("right", item.coefficient, item.compound, item.role)
        for item in _assign_side_properties(
            right_stoichiometry,
            right_molecules,
            known_properties,
        )
    )
    return left, right


def _source_terms(
    template: MnxrefReactionTemplate,
    chemicals: Mapping[str, MnxrefChemical],
    *,
    reverse: bool,
) -> tuple[CompletedTerm, ...]:
    rows: list[CompletedTerm] = []
    for term in template.terms:
        chemical = chemicals.get(term.mnxm_id)
        if chemical is None:
            raise ValueError(f"MNXref chemical is missing: {term.mnxm_id}")
        side = term.side
        if reverse:
            side = "right" if side == "left" else "left"
        rows.append(
            CompletedTerm(
                side=side,
                coefficient=term.coefficient,
                compound=_compound_from_mnxref(chemical),
                role="source_core",
            )
        )
    return tuple(rows)


def _candidate_exclusions(
    source_terms: Sequence[CompletedTerm],
    main_mnxm_id: str,
) -> Iterable[set[int]]:
    left = [index for index, term in enumerate(source_terms) if term.side == "left"]
    right = [index for index, term in enumerate(source_terms) if term.side == "right"]
    main_indexes = [
        index
        for index, term in enumerate(source_terms)
        if term.compound.source_mnxm_id == main_mnxm_id
    ]
    for main_index in main_indexes:
        main_side = left if main_index in left else right
        opposite = right if main_index in left else left
        other_main = [index for index in main_side if index != main_index]
        for main_extra_count in range(len(other_main) + 1):
            for main_extra in itertools.combinations(other_main, main_extra_count):
                for opposite_count in range(1, len(opposite) + 1):
                    for opposite_values in itertools.combinations(
                        opposite,
                        opposite_count,
                    ):
                        yield {main_index, *main_extra, *opposite_values}


def _hypothesis(
    *,
    candidate_id: str,
    step_id: str,
    template: MnxrefReactionTemplate,
    reverse: bool,
    core_left: Sequence[CompletedTerm],
    core_right: Sequence[CompletedTerm],
    source_terms: Sequence[CompletedTerm],
    excluded: set[int],
) -> CompletedReactionHypothesis | None:
    recovered = tuple(
        CompletedTerm(
            side=term.side,
            coefficient=term.coefficient,
            compound=term.compound,
            role="recovered_template_participant",
        )
        for index, term in enumerate(source_terms)
        if index not in excluded
    )
    terms = tuple([*core_left, *core_right, *recovered])
    if not is_balanced(terms):
        return None
    canonical_terms = sorted(
        (
            term.side,
            round(term.coefficient, 12),
            term.compound.compound_id,
            term.role,
        )
        for term in terms
    )
    payload = {
        "candidate_id": candidate_id,
        "step_id": step_id,
        "rule_id": template.rule_id,
        "source_mnxr_id": template.mnxr_id,
        "source_orientation": "right_to_left" if reverse else "left_to_right",
        "terms": canonical_terms,
    }
    hypothesis_id = "RP2STOICH:" + _canonical_sha256(payload)
    return CompletedReactionHypothesis(
        hypothesis_id=hypothesis_id,
        candidate_id=candidate_id,
        step_id=step_id,
        rule_id=template.rule_id,
        source_mnxr_id=template.mnxr_id,
        source_reference=template.reference,
        source_equation=template.equation,
        source_orientation=("right_to_left" if reverse else "left_to_right"),
        evidence_grade="rr02_mnxref_v3_template_balanced",
        terms=tuple(
            sorted(
                terms,
                key=lambda item: (
                    0 if item.side == "left" else 1,
                    item.role,
                    item.compound.compound_id,
                ),
            )
        ),
        balance_status="balanced",
        cofactor_reconstruction_status=(
            "complete" if recovered else "not_applicable"
        ),
    )


def reconstruct_retropath_step(
    step_row: Mapping[str, Any],
    index: MnxrefIndex,
    known_properties: Mapping[str, CompoundProperty],
    *,
    max_hypotheses: int = DEFAULT_MAX_STEP_HYPOTHESES,
) -> StepReconstructionResult:
    candidate_id = str(step_row.get("candidate_id") or "").strip()
    step_id = str(step_row.get("step_id") or "").strip()
    if not candidate_id or not step_id:
        raise ValueError("candidate_id and step_id are required")
    if str(step_row.get("step_source") or "").strip() != "retropath":
        raise ValueError("only retropath steps require source-template reconstruction")
    if isinstance(max_hypotheses, bool) or max_hypotheses < 1:
        raise ValueError("max_hypotheses must be a positive integer")
    rule_ids = [
        value.strip()
        for value in str(step_row.get("rule_ids") or "").split(";")
        if value.strip()
    ]
    core_left, core_right = core_terms_from_step(step_row, known_properties)
    templates = index.templates_for_rules(rule_ids)
    hypotheses_by_terms: dict[
        tuple[tuple[str, float, str], ...], CompletedReactionHypothesis
    ] = {}
    rejections: list[ReconstructionRejection] = []
    if not templates:
        rejections.append(
            ReconstructionRejection(
                candidate_id,
                step_id,
                ";".join(rule_ids),
                "",
                "source_template_missing",
                "RR02 rules have no indexed MNXref v3 source template",
            )
        )
    for template in templates:
        if template.parse_status != "ok":
            rejections.append(
                ReconstructionRejection(
                    candidate_id,
                    step_id,
                    template.rule_id,
                    template.mnxr_id,
                    "source_equation_invalid",
                    f"MNXref equation parse status: {template.parse_status}",
                )
            )
            continue
        if not template.balanced:
            rejections.append(
                ReconstructionRejection(
                    candidate_id,
                    step_id,
                    template.rule_id,
                    template.mnxr_id,
                    "source_equation_unbalanced",
                    "MNXref v3 did not mark the source equation balanced",
                )
            )
            continue
        if template.transport:
            rejections.append(
                ReconstructionRejection(
                    candidate_id,
                    step_id,
                    template.rule_id,
                    template.mnxr_id,
                    "source_transport_unsupported",
                    "transport templates are not supported for cytosolic RP2 steps",
                )
            )
            continue
        chemicals = index.chemicals(term.mnxm_id for term in template.terms)
        template_found = False
        try:
            for reverse in (False, True):
                source_terms = _source_terms(template, chemicals, reverse=reverse)
                for excluded in _candidate_exclusions(
                    source_terms,
                    template.main_mnxm_id,
                ):
                    hypothesis = _hypothesis(
                        candidate_id=candidate_id,
                        step_id=step_id,
                        template=template,
                        reverse=reverse,
                        core_left=core_left,
                        core_right=core_right,
                        source_terms=source_terms,
                        excluded=excluded,
                    )
                    if hypothesis is None:
                        continue
                    template_found = True
                    key = tuple(
                        sorted(
                            (
                                term.side,
                                round(term.coefficient, 12),
                                term.compound.compound_id,
                            )
                            for term in hypothesis.terms
                        )
                    )
                    previous = hypotheses_by_terms.get(key)
                    if previous is None or (
                        hypothesis.rule_id,
                        hypothesis.source_mnxr_id,
                        hypothesis.hypothesis_id,
                    ) < (
                        previous.rule_id,
                        previous.source_mnxr_id,
                        previous.hypothesis_id,
                    ):
                        hypotheses_by_terms[key] = hypothesis
        except ValueError as exc:
            rejections.append(
                ReconstructionRejection(
                    candidate_id,
                    step_id,
                    template.rule_id,
                    template.mnxr_id,
                    "source_chemical_incomplete",
                    str(exc),
                )
            )
            continue
        if not template_found:
            rejections.append(
                ReconstructionRejection(
                    candidate_id,
                    step_id,
                    template.rule_id,
                    template.mnxr_id,
                    "no_balanced_template_transfer",
                    "no source-core partition produced a balanced predicted equation",
                )
            )

    ordered = sorted(
        hypotheses_by_terms.values(),
        key=lambda item: (
            len(item.recovered_compound_ids),
            item.rule_id,
            item.source_mnxr_id,
            item.hypothesis_id,
        ),
    )
    truncated = len(ordered) > max_hypotheses
    ordered = ordered[:max_hypotheses]
    status = "complete" if ordered else "incomplete"
    return StepReconstructionResult(
        candidate_id=candidate_id,
        step_id=step_id,
        status=status,
        hypotheses=tuple(ordered),
        rejections=tuple(
            sorted(
                rejections,
                key=lambda item: (
                    item.rule_id,
                    item.source_mnxr_id,
                    item.reason_code,
                    item.reason_detail,
                ),
            )
        ),
        truncated=truncated,
    )


def enumerate_candidate_hypotheses(
    step_results: Sequence[StepReconstructionResult],
    *,
    max_combinations: int = 32,
) -> tuple[tuple[tuple[CompletedReactionHypothesis, ...], ...], bool]:
    if isinstance(max_combinations, bool) or max_combinations < 1:
        raise ValueError("max_combinations must be a positive integer")
    if not step_results:
        return (tuple(),), False
    if any(not result.hypotheses for result in step_results):
        return tuple(), False
    combinations: list[tuple[CompletedReactionHypothesis, ...]] = []
    truncated = False
    for values in itertools.product(*(result.hypotheses for result in step_results)):
        if len(combinations) >= max_combinations:
            truncated = True
            break
        combinations.append(tuple(values))
    return tuple(combinations), truncated


__all__ = [
    "CompletedReactionHypothesis",
    "CompletedTerm",
    "CompoundProperty",
    "DEFAULT_MAX_STEP_HYPOTHESES",
    "ReconstructionRejection",
    "STOICHIOMETRY_SCHEMA_VERSION",
    "StepReconstructionResult",
    "core_terms_from_step",
    "enumerate_candidate_hypotheses",
    "is_balanced",
    "load_p2_compound_properties",
    "parse_formula",
    "reconstruct_retropath_step",
]
