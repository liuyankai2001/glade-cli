"""Conservative structural identity matching for RetroPath artifacts.

RetroPath/KNIME may drop stereochemical layers while preserving molecular
connectivity.  This module distinguishes that lossy representation from both
an exact match and an explicit stereochemical conflict.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdMolDescriptors


STRUCTURE_MATCH_TYPES = frozenset(
    {"exact", "stereo_missing", "stereo_conflict", "different"}
)
_FORMULA_CHARGE_SUFFIX = re.compile(r"[+-]\d*$")
_STEREO_LAYER = re.compile(r"/(?:t|b)[^/]+")


def _formula_without_charge(value: str) -> str:
    return _FORMULA_CHARGE_SUFFIX.sub("", str(value).strip())


@dataclass(frozen=True)
class StructureIdentity:
    inchi: str
    full_inchikey: str
    stereo_stripped_inchi: str
    stereo_stripped_inchikey: str
    formula: str
    charge: int
    stereo_specified: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "inchi": self.inchi,
            "full_inchikey": self.full_inchikey,
            "stereo_stripped_inchi": self.stereo_stripped_inchi,
            "stereo_stripped_inchikey": self.stereo_stripped_inchikey,
            "formula": self.formula,
            "charge": self.charge,
            "stereo_specified": self.stereo_specified,
        }


@dataclass(frozen=True)
class StructureMatch:
    match_type: str
    reference: StructureIdentity
    observed: StructureIdentity
    detail: str

    def __post_init__(self) -> None:
        if self.match_type not in STRUCTURE_MATCH_TYPES:
            raise ValueError(f"unsupported structure match type: {self.match_type}")

    @property
    def accepted(self) -> bool:
        return self.match_type in {"exact", "stereo_missing"}

    @property
    def review_required(self) -> bool:
        return self.match_type == "stereo_missing"

    def to_dict(self) -> dict[str, object]:
        return {
            "match_type": self.match_type,
            "detail": self.detail,
            "review_required": self.review_required,
            "reference": self.reference.to_dict(),
            "observed": self.observed.to_dict(),
        }


def structure_identity(inchi: str) -> StructureIdentity:
    normalized = str(inchi).strip()
    if not normalized.startswith("InChI=1S/"):
        raise ValueError("structure identity requires a standard InChI v1 value")
    try:
        molecule = Chem.MolFromInchi(normalized, sanitize=True, removeHs=True)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("RDKit could not parse structure identity InChI") from exc
    if molecule is None:
        raise ValueError("RDKit could not parse structure identity InChI")
    try:
        standard_inchi = rd_inchi.MolToInchi(molecule).strip()
        full_key = rd_inchi.InchiToInchiKey(standard_inchi).strip().upper()
        formula = _formula_without_charge(rdMolDescriptors.CalcMolFormula(molecule))
        charge = int(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()))
        achiral = Chem.Mol(molecule)
        Chem.RemoveStereochemistry(achiral)
        achiral_inchi = rd_inchi.MolToInchi(achiral).strip()
        achiral_key = rd_inchi.InchiToInchiKey(achiral_inchi).strip().upper()
    except (RuntimeError, ValueError) as exc:
        raise ValueError("RDKit could not standardize structure identity") from exc
    if (
        not standard_inchi.startswith("InChI=1S/")
        or not achiral_inchi.startswith("InChI=1S/")
        or not full_key
        or not achiral_key
    ):
        raise ValueError("RDKit returned an incomplete structure identity")
    return StructureIdentity(
        inchi=standard_inchi,
        full_inchikey=full_key,
        stereo_stripped_inchi=achiral_inchi,
        stereo_stripped_inchikey=achiral_key,
        formula=formula,
        charge=charge,
        stereo_specified=bool(_STEREO_LAYER.search(standard_inchi)),
    )


def compare_structure_identities(
    reference: StructureIdentity,
    observed: StructureIdentity,
) -> StructureMatch:
    same_formula_charge = (
        reference.formula == observed.formula
        and reference.charge == observed.charge
    )
    if reference.full_inchikey == observed.full_inchikey and same_formula_charge:
        return StructureMatch("exact", reference, observed, "full InChIKey match")
    if (
        reference.stereo_stripped_inchikey
        != observed.stereo_stripped_inchikey
        or not same_formula_charge
    ):
        return StructureMatch(
            "different",
            reference,
            observed,
            "connectivity, formula, or formal charge differs",
        )
    if reference.stereo_specified and observed.stereo_specified:
        return StructureMatch(
            "stereo_conflict",
            reference,
            observed,
            "both structures specify incompatible stereochemistry",
        )
    return StructureMatch(
        "stereo_missing",
        reference,
        observed,
        "connectivity/formula/charge match but stereochemistry is absent on one side",
    )


def compare_inchis(reference_inchi: str, observed_inchi: str) -> StructureMatch:
    return compare_structure_identities(
        structure_identity(reference_inchi),
        structure_identity(observed_inchi),
    )


__all__ = [
    "STRUCTURE_MATCH_TYPES",
    "StructureIdentity",
    "StructureMatch",
    "compare_inchis",
    "compare_structure_identities",
    "structure_identity",
]
