"""Parse RetroPath2 scope artifacts into auditable prediction objects.

The KNIME workflow emits a bipartite compound/transformation network rather
than ready-to-use pathways.  This module validates that network against the
P2 source/sink bundle and the exact RetroRules file used by the service.  It
does not reverse reactions or merge KEGG witnesses; those are later stages.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

import rdkit
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdChemReactions, rdMolDescriptors

from src.pathway_analyze.retropath_client import RetroPathClientRun
from src.pathway_analyze.retropath_input import RetroPathInputBundle
from src.pathway_analyze.retropath_identity import (
    StructureIdentity,
    compare_structure_identities,
    structure_identity,
)
from src.pathway_analyze.retropath_models import (
    PredictedCompound,
    PredictedReaction,
)

RETROPATH_RESULTS_COLUMNS = (
    "Initial source",
    "Transformation ID",
    "Reaction SMILES",
    "Substrate SMILES",
    "Substrate InChI",
    "Product SMILES",
    "Product InChI",
    "In Sink",
    "Sink name",
    "Diameter",
    "Rule ID",
    "EC number",
    "Score",
    "Iteration",
)
RETRO_RULES_REQUIRED_COLUMNS = (
    "Rule ID",
    "Rule",
    "EC number",
    "Reaction order",
    "Diameter",
    "Score",
)

SCOPE_JSON_FILE_NAMES = ("target_scope.json", "scope.json")
SCOPE_CSV_FILE_NAMES = ("target_scope.csv", "scope.csv")
RESULTS_CSV_FILE_NAMES = ("results.csv",)

_NON_STRUCTURAL_AUXILIARY_SMILES = frozenset(
    {"[H+]", "[H-]", "[H]", "[e-]"}
)
_RDKIT_UNSUPPORTED_REACTION_COMPONENTS = frozenset({"[e-]"})

_SOURCE_REACTION_PATTERN = re.compile(
    r"(?:MNXR\d+|RHEA(?::|_)\d+|(?<![A-Z0-9])R\d{5}(?!\d))",
    re.IGNORECASE,
)


class RetroPathParseError(ValueError):
    """A run-level artifact or provenance error that makes parsing unsafe."""

    def __init__(self, code: str, detail: str) -> None:
        self.code = str(code).strip()
        self.detail = str(detail).strip()
        super().__init__(f"{self.code}: {self.detail}")


@dataclass(frozen=True)
class RuleEvidence:
    """Rule-level evidence read from the exact RR02 file used for expansion."""

    rule_id: str
    diameter: int
    score_raw: float
    score_normalized: Optional[float]
    score_semantics: str
    ec_numbers: Tuple[str, ...]
    legacy_ids: Tuple[str, ...]
    source_reaction_ids: Tuple[str, ...]
    reaction_direction: Optional[str]
    rule_relative_direction: Optional[str]
    rule_usage: Optional[str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "rule_id": self.rule_id,
            "diameter": self.diameter,
            "score_raw": self.score_raw,
            "score_normalized": self.score_normalized,
            "score_semantics": self.score_semantics,
            "ec_numbers": list(self.ec_numbers),
            "legacy_ids": list(self.legacy_ids),
            "source_reaction_ids": list(self.source_reaction_ids),
            "reaction_direction": self.reaction_direction,
            "rule_relative_direction": self.rule_relative_direction,
            "rule_usage": self.rule_usage,
        }


@dataclass(frozen=True)
class SinkMatch:
    """A graded structural match to one or more P2 cumulative sink entries."""

    compound_id: str
    inchikey: str
    representative_kegg_id: str
    kegg_ids: Tuple[str, ...]
    minimum_depth: int
    wrapper_in_sink: bool
    wrapper_sink_names: Tuple[str, ...]
    match_type: str = "exact"
    stereo_review_required: bool = False
    stereo_stripped_inchikey: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "compound_id": self.compound_id,
            "inchikey": self.inchikey,
            "representative_kegg_id": self.representative_kegg_id,
            "kegg_ids": list(self.kegg_ids),
            "minimum_depth": self.minimum_depth,
            "wrapper_in_sink": self.wrapper_in_sink,
            "wrapper_sink_names": list(self.wrapper_sink_names),
            "match_type": self.match_type,
            "stereo_review_required": self.stereo_review_required,
            "stereo_stripped_inchikey": self.stereo_stripped_inchikey,
        }


@dataclass(frozen=True)
class AuxiliaryFragment:
    """A proton/electron-like fragment excluded from structural connectivity."""

    inchi: str
    smiles: str
    reason_code: str

    def to_dict(self) -> dict[str, str]:
        return {
            "inchi": self.inchi,
            "smiles": self.smiles,
            "reason_code": self.reason_code,
        }


@dataclass(frozen=True)
class ParsedCompoundNode:
    """One structure-addressable compound node in the RetroPath scope."""

    compound: PredictedCompound
    is_target: bool
    wrapper_in_sink: bool
    wrapper_sink_names: Tuple[str, ...]
    target_match_type: Optional[str] = None
    sink_match_type: Optional[str] = None
    stereo_review_required: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "compound": self.compound.to_dict(),
            "is_target": self.is_target,
            "wrapper_in_sink": self.wrapper_in_sink,
            "wrapper_sink_names": list(self.wrapper_sink_names),
            "target_match_type": self.target_match_type,
            "sink_match_type": self.sink_match_type,
            "stereo_review_required": self.stereo_review_required,
        }


@dataclass(frozen=True)
class NetworkRejection:
    """A transformation excluded while validating raw network artifacts."""

    reason_code: str
    reason_detail: str
    transformation_id: Optional[str] = None
    compound_id: Optional[str] = None

    def to_dict(self) -> dict[str, Optional[str]]:
        return {
            "reason_code": self.reason_code,
            "reason_detail": self.reason_detail,
            "transformation_id": self.transformation_id,
            "compound_id": self.compound_id,
        }


@dataclass(frozen=True)
class ParsedTransformation:
    """One concrete retrosynthetic transformation and its rule alternatives."""

    transformation_id: str
    substrate_compound_id: str
    product_compound_ids: Tuple[str, ...]
    reaction_smiles: str
    iteration: int
    diameter: int
    score_raw: float
    score_semantics: str
    rule_ids: Tuple[str, ...]
    reported_ec_numbers: Tuple[str, ...]
    reaction_variants: Tuple[PredictedReaction, ...]
    auxiliary_fragments: Tuple[AuxiliaryFragment, ...]

    @property
    def unique_product_compound_ids(self) -> Tuple[str, ...]:
        return tuple(sorted(set(self.product_compound_ids)))

    @property
    def minimum_rule_specificity(self) -> int:
        return min(
            reaction.rule_specificity
            for reaction in self.reaction_variants
            if reaction.rule_specificity is not None
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "transformation_id": self.transformation_id,
            "substrate_compound_id": self.substrate_compound_id,
            "product_compound_ids": list(self.product_compound_ids),
            "reaction_smiles": self.reaction_smiles,
            "iteration": self.iteration,
            "diameter": self.diameter,
            "score_raw": self.score_raw,
            "score_semantics": self.score_semantics,
            "rule_ids": list(self.rule_ids),
            "reported_ec_numbers": list(self.reported_ec_numbers),
            "reaction_variants": [item.to_dict() for item in self.reaction_variants],
            "auxiliary_fragments": [
                item.to_dict() for item in self.auxiliary_fragments
            ],
        }


@dataclass(frozen=True)
class ParsedRetroPathNetwork:
    """Validated retrosynthetic network ready for complete-path enumeration."""

    status: str
    target_compound_id: str
    max_steps: int
    compounds: Tuple[ParsedCompoundNode, ...]
    transformations: Tuple[ParsedTransformation, ...]
    sink_matches: Tuple[SinkMatch, ...]
    rule_evidence: Tuple[RuleEvidence, ...]
    rejections: Tuple[NetworkRejection, ...]
    warnings: Tuple[str, ...]
    source_in_sink: Optional[SinkMatch] = None

    @property
    def predicted_compounds(self) -> Tuple[PredictedCompound, ...]:
        return tuple(item.compound for item in self.compounds)

    @property
    def predicted_reactions(self) -> Tuple[PredictedReaction, ...]:
        return tuple(
            reaction
            for transformation in self.transformations
            for reaction in transformation.reaction_variants
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "target_compound_id": self.target_compound_id,
            "max_steps": self.max_steps,
            "compounds": [item.to_dict() for item in self.compounds],
            "transformations": [item.to_dict() for item in self.transformations],
            "sink_matches": [item.to_dict() for item in self.sink_matches],
            "rule_evidence": [item.to_dict() for item in self.rule_evidence],
            "rejections": [item.to_dict() for item in self.rejections],
            "warnings": list(self.warnings),
            "source_in_sink": (
                None if self.source_in_sink is None else self.source_in_sink.to_dict()
            ),
        }


@dataclass
class _NodeAccumulator:
    compound: PredictedCompound
    is_target: bool = False
    wrapper_flags: set[bool] | None = None
    sink_names: set[str] | None = None
    target_match_types: set[str] | None = None

    def __post_init__(self) -> None:
        if self.wrapper_flags is None:
            self.wrapper_flags = set()
        if self.sink_names is None:
            self.sink_names = set()
        if self.target_match_types is None:
            self.target_match_types = set()


@dataclass(frozen=True)
class _RawProduct:
    inchi: str
    smiles: str
    in_sink: bool
    sink_names: Tuple[str, ...]


@dataclass(frozen=True)
class _RawTransformation:
    transformation_id: str
    reaction_smiles: str
    substrate_inchi: str
    substrate_smiles: str
    products: Tuple[_RawProduct, ...]
    diameter: int
    rule_ids: Tuple[str, ...]
    ec_numbers: Tuple[str, ...]
    score_raw: float
    iteration: int


@dataclass(frozen=True)
class _ArtifactCompoundResolution:
    compound: Optional[PredictedCompound]
    auxiliary: Optional[AuxiliaryFragment]
    target_match_type: Optional[str] = None
    sink_candidates: Tuple[PredictedCompound, ...] = tuple()
    sink_match_type: Optional[str] = None
    observed_identity: Optional[StructureIdentity] = None


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _stable_unique(values: Iterable[str]) -> Tuple[str, ...]:
    return tuple(sorted({str(value).strip() for value in values if str(value).strip()}))


def _parse_bracket_list(
    value: Any, *, remove_empty_markers: bool = True
) -> Tuple[str, ...]:
    text = "" if value is None else str(value).strip()
    if not text:
        return tuple()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1].strip()
    if not text:
        return tuple()
    try:
        items = next(csv.reader([text], skipinitialspace=True))
    except csv.Error as exc:
        raise ValueError(f"invalid bracket-list value: {value!r}") from exc
    normalized = []
    for item in items:
        candidate = item.strip().strip("'\"").strip()
        if not candidate:
            continue
        if remove_empty_markers and candidate.upper() in {"NONE", "NULL", "NOEC"}:
            continue
        normalized.append(candidate)
    return _stable_unique(normalized)


def _parse_int(value: Any, field_name: str, *, minimum: int = 0) -> int:
    try:
        parsed = int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be an integer") from exc
    if parsed < minimum:
        raise ValueError(f"{field_name} must be at least {minimum}")
    return parsed


def _parse_float(value: Any, field_name: str) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{field_name} must be finite")
    return parsed


def _parse_sink_flag(value: Any) -> bool:
    text = str(value).strip().lower()
    if text in {"1", "true", "yes"}:
        return True
    if text in {"0", "false", "no"}:
        return False
    raise ValueError(f"In Sink must be 0/1 or false/true, got {value!r}")


def _validate_reaction_smiles(value: Any) -> str:
    reaction_smiles = str(value).strip()
    if reaction_smiles.count(">>") != 1:
        raise ValueError("Reaction SMILES must contain exactly one '>>'")
    left, right = reaction_smiles.split(">>", 1)
    if not left.strip() or not right.strip():
        raise ValueError("Reaction SMILES must contain non-empty sides")
    rdkit_left = ".".join(
        component
        for component in (item.strip() for item in left.split("."))
        if component not in _RDKIT_UNSUPPORTED_REACTION_COMPONENTS
    )
    rdkit_right = ".".join(
        component
        for component in (item.strip() for item in right.split("."))
        if component not in _RDKIT_UNSUPPORTED_REACTION_COMPONENTS
    )
    if not rdkit_left or not rdkit_right:
        raise ValueError(
            "Reaction SMILES must contain structural components on both sides"
        )
    try:
        reaction = rdChemReactions.ReactionFromSmarts(
            f"{rdkit_left}>>{rdkit_right}",
            useSmiles=True,
        )
    except (RuntimeError, ValueError) as exc:
        raise ValueError("Reaction SMILES is not parseable by RDKit") from exc
    if reaction is None:
        raise ValueError("Reaction SMILES is not parseable by RDKit")
    return reaction_smiles


def _artifact_path(raw_dir: Path, names: Sequence[str]) -> Optional[Path]:
    for name in names:
        direct = raw_dir / name
        if direct.is_file():
            return direct
    candidates = sorted(
        {
            path.resolve()
            for name in names
            for path in raw_dir.rglob(name)
            if path.is_file()
        }
    )
    if len(candidates) > 1:
        raise RetroPathParseError(
            "artifact_inconsistent",
            f"multiple candidate artifacts found for {', '.join(names)}",
        )
    return candidates[0] if candidates else None


def _read_csv(path: Path) -> tuple[Tuple[str, ...], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = tuple(reader.fieldnames or tuple())
            rows = [
                {
                    str(key): "" if value is None else str(value)
                    for key, value in row.items()
                }
                for row in reader
                if any(str(value or "").strip() for value in row.values())
            ]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise RetroPathParseError(
            "artifact_schema_invalid",
            f"cannot read CSV artifact {path.name}: {exc}",
        ) from exc
    return fieldnames, rows


def _select_scope_table(raw_dir: Path) -> tuple[Path, list[dict[str, str]]]:
    candidates = []
    for names in (SCOPE_CSV_FILE_NAMES, RESULTS_CSV_FILE_NAMES):
        path = _artifact_path(raw_dir, names)
        if path is None:
            continue
        fieldnames, rows = _read_csv(path)
        if not rows:
            continue
        missing = [
            column for column in RETROPATH_RESULTS_COLUMNS if column not in fieldnames
        ]
        if missing:
            raise RetroPathParseError(
                "artifact_schema_invalid",
                f"{path.name} is missing required columns: {', '.join(missing)}",
            )
        candidates.append((path, rows))
        if path.name in SCOPE_CSV_FILE_NAMES:
            break
    if not candidates:
        raise RetroPathParseError(
            "artifact_missing",
            "a non-empty target_scope.csv, scope.csv, or results.csv is required",
        )
    return candidates[0]


def _group_rows(
    rows: Sequence[Mapping[str, str]],
) -> tuple[
    Tuple[_RawTransformation, ...],
    Tuple[NetworkRejection, ...],
]:
    grouped: dict[str, list[Mapping[str, str]]] = {}
    rejections: list[NetworkRejection] = []
    for index, row in enumerate(rows, start=2):
        transformation_id = str(row.get("Transformation ID", "")).strip()
        if not transformation_id:
            rejections.append(
                NetworkRejection(
                    "artifact_schema_invalid",
                    f"row {index} has no Transformation ID",
                )
            )
            continue
        grouped.setdefault(transformation_id, []).append(row)

    transformations: list[_RawTransformation] = []
    common_fields = (
        "Reaction SMILES",
        "Substrate SMILES",
        "Substrate InChI",
        "Diameter",
        "Rule ID",
        "EC number",
        "Score",
        "Iteration",
    )
    for transformation_id in sorted(grouped):
        group = grouped[transformation_id]
        try:
            for field_name in common_fields:
                observed = {str(row.get(field_name, "")).strip() for row in group}
                if len(observed) != 1:
                    raise ValueError(f"rows disagree on {field_name}")
            first = group[0]
            reaction_smiles = _validate_reaction_smiles(first["Reaction SMILES"])
            rule_ids = _parse_bracket_list(first["Rule ID"])
            if not rule_ids:
                raise ValueError("Rule ID list is empty")
            products = tuple(
                _RawProduct(
                    inchi=str(row["Product InChI"]).strip(),
                    smiles=str(row["Product SMILES"]).strip(),
                    in_sink=_parse_sink_flag(row["In Sink"]),
                    sink_names=_parse_bracket_list(
                        row["Sink name"],
                        remove_empty_markers=False,
                    ),
                )
                for row in group
            )
            transformations.append(
                _RawTransformation(
                    transformation_id=transformation_id,
                    reaction_smiles=reaction_smiles,
                    substrate_inchi=str(first["Substrate InChI"]).strip(),
                    substrate_smiles=str(first["Substrate SMILES"]).strip(),
                    products=products,
                    diameter=_parse_int(first["Diameter"], "Diameter"),
                    rule_ids=rule_ids,
                    ec_numbers=_parse_bracket_list(first["EC number"]),
                    score_raw=_parse_float(first["Score"], "Score"),
                    iteration=_parse_int(first["Iteration"], "Iteration"),
                )
            )
        except ValueError as exc:
            rejections.append(
                NetworkRejection(
                    "artifact_schema_invalid",
                    str(exc),
                    transformation_id=transformation_id,
                )
            )
    return tuple(transformations), tuple(rejections)


def _extract_source_reactions(values: Iterable[str]) -> Tuple[str, ...]:
    matches = []
    for value in values:
        matches.extend(
            match.group(0).upper().replace("RHEA_", "RHEA:")
            for match in _SOURCE_REACTION_PATTERN.finditer(value)
        )
    return _stable_unique(matches)


def _read_rule_evidence(
    rules_path: Path,
    expected_sha256: str,
    requested_rule_ids: Iterable[str],
) -> tuple[dict[str, RuleEvidence], Tuple[str, ...], str]:
    if not rules_path.is_file():
        raise RetroPathParseError(
            "rules_file_missing",
            f"RetroRules file does not exist: {rules_path}",
        )
    observed_sha256 = _sha256_file(rules_path)
    if observed_sha256 != expected_sha256.lower():
        raise RetroPathParseError(
            "rules_checksum_mismatch",
            f"expected {expected_sha256.lower()}, observed {observed_sha256}",
        )

    requested = set(requested_rule_ids)
    raw_matches: dict[str, list[dict[str, str]]] = {
        rule_id: [] for rule_id in requested
    }
    any_score = False
    any_score_above_one = False
    try:
        with rules_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = tuple(reader.fieldnames or tuple())
            missing = [
                column
                for column in RETRO_RULES_REQUIRED_COLUMNS
                if column not in fieldnames
            ]
            if missing:
                raise RetroPathParseError(
                    "artifact_schema_invalid",
                    f"RetroRules file is missing columns: {', '.join(missing)}",
                )
            for row in reader:
                try:
                    score = float(str(row.get("Score", "")).strip())
                except ValueError:
                    score = math.nan
                if math.isfinite(score):
                    any_score = True
                    any_score_above_one = any_score_above_one or score > 1.0
                rule_id = str(row.get("Rule ID", "")).strip()
                if rule_id in requested:
                    raw_matches[rule_id].append(
                        {
                            str(key): "" if value is None else str(value)
                            for key, value in row.items()
                        }
                    )
    except RetroPathParseError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise RetroPathParseError(
            "artifact_schema_invalid",
            f"cannot read RetroRules file: {exc}",
        ) from exc
    if not any_score:
        raise RetroPathParseError(
            "artifact_schema_invalid",
            "RetroRules Score column contains no finite value",
        )
    score_semantics = "lower_is_better" if any_score_above_one else "higher_is_better"

    evidence: dict[str, RuleEvidence] = {}
    missing_ids = []
    for rule_id in sorted(requested):
        matches = raw_matches.get(rule_id, [])
        if not matches:
            missing_ids.append(rule_id)
            continue
        diameters = {_parse_int(row["Diameter"], "Diameter") for row in matches}
        scores = {_parse_float(row["Score"], "Score") for row in matches}
        if len(diameters) != 1 or len(scores) != 1:
            raise RetroPathParseError(
                "artifact_inconsistent",
                f"RetroRules rows disagree for {rule_id}",
            )
        normalized_values = {
            _parse_float(row["Score normalized"], "Score normalized")
            for row in matches
            if str(row.get("Score normalized", "")).strip()
        }
        if len(normalized_values) > 1:
            raise RetroPathParseError(
                "artifact_inconsistent",
                f"normalized scores disagree for {rule_id}",
            )
        ec_numbers = _stable_unique(
            ec
            for row in matches
            for ec in str(row.get("EC number", "")).split(";")
            if ec.strip().upper() not in {"", "NOEC", "NONE", "NULL"}
        )
        legacy_ids = _stable_unique(
            item
            for row in matches
            for item in str(row.get("Legacy ID", "")).split(";")
            if item.strip()
        )
        directions = _stable_unique(
            row.get("Reaction direction", "") for row in matches
        )
        relative_directions = _stable_unique(
            row.get("Rule relative direction", "") for row in matches
        )
        usages = _stable_unique(row.get("Rule usage", "") for row in matches)
        evidence[rule_id] = RuleEvidence(
            rule_id=rule_id,
            diameter=next(iter(diameters)),
            score_raw=next(iter(scores)),
            score_normalized=(
                None if not normalized_values else next(iter(normalized_values))
            ),
            score_semantics=score_semantics,
            ec_numbers=ec_numbers,
            legacy_ids=legacy_ids,
            source_reaction_ids=_extract_source_reactions(legacy_ids),
            reaction_direction=directions[0] if len(directions) == 1 else None,
            rule_relative_direction=(
                relative_directions[0] if len(relative_directions) == 1 else None
            ),
            rule_usage=usages[0] if len(usages) == 1 else None,
        )
    return evidence, tuple(missing_ids), score_semantics


def _pseudo_fragment(inchi: str, smiles: str) -> Optional[AuxiliaryFragment]:
    normalized_inchi = inchi.strip()
    normalized_smiles = smiles.strip()
    if normalized_inchi in {"InChI=1S", "InChI=1"} or normalized_inchi.startswith(
        ("InChI=1S/p", "InChI=1/p")
    ):
        return AuxiliaryFragment(
            normalized_inchi,
            normalized_smiles,
            "non_structural_auxiliary_fragment",
        )
    if normalized_smiles in _NON_STRUCTURAL_AUXILIARY_SMILES:
        return AuxiliaryFragment(
            normalized_inchi,
            normalized_smiles,
            "non_structural_auxiliary_fragment",
        )
    return None


def _compound_from_artifact(
    inchi: str,
    smiles: str,
    *,
    target: PredictedCompound,
    target_identity: StructureIdentity,
    sink_by_inchikey: Mapping[str, PredictedCompound],
    sink_by_stereo_stripped_inchikey: Mapping[str, Tuple[PredictedCompound, ...]],
    sink_identities: Mapping[str, StructureIdentity],
    provenance: str,
) -> _ArtifactCompoundResolution:
    pseudo = _pseudo_fragment(inchi, smiles)
    if pseudo is not None:
        return _ArtifactCompoundResolution(None, pseudo)
    normalized_inchi = str(inchi).strip()
    if not normalized_inchi.startswith("InChI=1S/"):
        raise ValueError("compound does not contain a standard InChI v1 structure")
    try:
        molecule = Chem.MolFromInchi(normalized_inchi, sanitize=True, removeHs=True)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("RDKit could not parse compound InChI") from exc
    if molecule is None:
        raise ValueError("RDKit could not parse compound InChI")
    if not any(atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms()):
        return _ArtifactCompoundResolution(
            None,
            AuxiliaryFragment(
                normalized_inchi,
                str(smiles).strip(),
                "no_heavy_atom_auxiliary_fragment",
            ),
        )
    try:
        standard_inchi = rd_inchi.MolToInchi(molecule).strip()
        inchikey = rd_inchi.InchiToInchiKey(standard_inchi).strip().upper()
        canonical_smiles = Chem.MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=True,
        ).strip()
        formula = rdMolDescriptors.CalcMolFormula(molecule).strip()
        charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    except (RuntimeError, ValueError) as exc:
        raise ValueError("RDKit could not standardize compound structure") from exc
    if not standard_inchi.startswith("InChI=1S/") or not inchikey:
        raise ValueError("RDKit did not produce a complete standard structure")

    observed_identity = structure_identity(standard_inchi)
    target_match = compare_structure_identities(target_identity, observed_identity)

    exact_sink = sink_by_inchikey.get(inchikey)
    sink_candidates: Tuple[PredictedCompound, ...] = tuple()
    sink_match_type: Optional[str] = None
    if exact_sink is not None:
        sink_candidates = (exact_sink,)
        sink_match_type = "exact"
    else:
        compatible = []
        for candidate in sink_by_stereo_stripped_inchikey.get(
            observed_identity.stereo_stripped_inchikey,
            tuple(),
        ):
            candidate_identity = sink_identities[candidate.compound_id]
            match = compare_structure_identities(candidate_identity, observed_identity)
            if match.match_type == "stereo_missing":
                compatible.append(candidate)
        if compatible:
            sink_candidates = tuple(
                sorted(
                    compatible,
                    key=lambda item: (
                        item.minimum_depth if item.minimum_depth is not None else 10**9,
                        item.compound_id,
                    ),
                )
            )
            sink_match_type = "stereo_missing"

    if target_match.match_type == "exact":
        return _ArtifactCompoundResolution(
            target,
            None,
            target_match_type="exact",
            sink_candidates=sink_candidates,
            sink_match_type=sink_match_type,
            observed_identity=observed_identity,
        )
    if target_match.match_type == "stereo_missing":
        return _ArtifactCompoundResolution(
            target,
            None,
            target_match_type="stereo_missing",
            sink_candidates=sink_candidates,
            sink_match_type=sink_match_type,
            observed_identity=observed_identity,
        )
    if exact_sink is not None:
        return _ArtifactCompoundResolution(
            exact_sink,
            None,
            sink_candidates=sink_candidates,
            sink_match_type="exact",
            observed_identity=observed_identity,
        )
    return _ArtifactCompoundResolution(
        PredictedCompound.create(
            inchi=standard_inchi,
            inchikey=inchikey,
            isomeric_smiles=canonical_smiles,
            formula=formula,
            charge=charge,
            structure_provenance=(
                "inchi:standard-v1",
                f"retropath_artifact:{provenance}",
                f"rdkit:{rdkit.__version__}",
            ),
        ),
        None,
        sink_candidates=sink_candidates,
        sink_match_type=sink_match_type,
        observed_identity=observed_identity,
    )


def _reaction_product_multiplicities(reaction_smiles: str) -> dict[str, int]:
    """Count structural right-hand components by stereo-stripped InChIKey.

    RetroPath's results table may contain one row per unique product structure,
    even when Reaction SMILES contains the same product more than once.  This
    helper restores only multiplicities for parseable heavy-atom components;
    components not represented by result rows remain P8 auxiliary evidence.
    """

    _, right = reaction_smiles.split(">>", 1)
    counts: dict[str, int] = {}
    for component in (item.strip() for item in right.split(".")):
        if not component or component in _NON_STRUCTURAL_AUXILIARY_SMILES:
            continue
        try:
            molecule = Chem.MolFromSmiles(component, sanitize=True)
        except (RuntimeError, ValueError):
            molecule = None
        if molecule is None or not any(
            atom.GetAtomicNum() > 1 for atom in molecule.GetAtoms()
        ):
            continue
        try:
            component_inchi = rd_inchi.MolToInchi(molecule).strip()
            key = structure_identity(component_inchi).stereo_stripped_inchikey
        except (RuntimeError, ValueError):
            continue
        if key:
            counts[key] = counts.get(key, 0) + 1
    return counts


def _json_topology(path: Path) -> dict[str, tuple[str, Tuple[str, ...]]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        elements = payload["elements"]
        raw_nodes = elements["nodes"]
        raw_edges = elements["edges"]
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise RetroPathParseError(
            "artifact_schema_invalid",
            f"invalid scope JSON: {exc}",
        ) from exc
    if not isinstance(raw_nodes, list) or not isinstance(raw_edges, list):
        raise RetroPathParseError(
            "artifact_schema_invalid",
            "scope JSON nodes and edges must be arrays",
        )

    node_types: dict[str, str] = {}
    compound_keys: dict[str, str] = {}
    for raw_node in raw_nodes:
        try:
            data = raw_node["data"]
            node_id = str(data["id"]).strip()
            node_type = str(data["type"]).strip().lower()
        except (KeyError, TypeError) as exc:
            raise RetroPathParseError(
                "artifact_schema_invalid",
                "scope JSON node is missing id or type",
            ) from exc
        if (
            not node_id
            or node_id in node_types
            or node_type not in {"compound", "reaction"}
        ):
            raise RetroPathParseError(
                "artifact_schema_invalid",
                f"invalid or duplicate scope JSON node: {node_id!r}",
            )
        node_types[node_id] = node_type
        if node_type == "compound":
            inchi = str(data.get("InChI", "")).strip()
            pseudo = _pseudo_fragment(inchi, str(data.get("SMILES", "")))
            if pseudo is not None:
                continue
            try:
                key = rd_inchi.InchiToInchiKey(inchi).strip().upper()
            except (RuntimeError, ValueError) as exc:
                raise RetroPathParseError(
                    "artifact_schema_invalid",
                    f"scope JSON compound {node_id} has invalid InChI",
                ) from exc
            if not key:
                raise RetroPathParseError(
                    "artifact_schema_invalid",
                    f"scope JSON compound {node_id} has no InChIKey",
                )
            compound_keys[node_id] = key

    incoming: dict[str, list[str]] = {}
    outgoing: dict[str, list[str]] = {}
    for raw_edge in raw_edges:
        try:
            data = raw_edge["data"]
            source = str(data["source"]).strip()
            target = str(data["target"]).strip()
        except (KeyError, TypeError) as exc:
            raise RetroPathParseError(
                "artifact_schema_invalid",
                "scope JSON edge is missing source or target",
            ) from exc
        if source not in node_types or target not in node_types:
            raise RetroPathParseError(
                "artifact_inconsistent",
                f"scope edge references an unknown node: {source}->{target}",
            )
        outgoing.setdefault(source, []).append(target)
        incoming.setdefault(target, []).append(source)

    topology: dict[str, tuple[str, Tuple[str, ...]]] = {}
    for node_id, node_type in node_types.items():
        if node_type != "reaction":
            continue
        substrate_nodes = [
            source
            for source in incoming.get(node_id, [])
            if node_types[source] == "compound" and source in compound_keys
        ]
        product_nodes = [
            target
            for target in outgoing.get(node_id, [])
            if node_types[target] == "compound" and target in compound_keys
        ]
        if len(substrate_nodes) != 1:
            raise RetroPathParseError(
                "artifact_inconsistent",
                f"reaction {node_id} must have exactly one structural substrate",
            )
        topology[node_id] = (
            compound_keys[substrate_nodes[0]],
            tuple(sorted({compound_keys[item] for item in product_nodes})),
        )
    return topology


def _max_steps(client_run: RetroPathClientRun) -> int:
    parameters = dict(client_run.result.parameters)
    value = parameters.get("max_steps", 3)
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 10:
        raise RetroPathParseError(
            "artifact_inconsistent",
            "P3 result contains an invalid max_steps value",
        )
    return value


def _sink_match(
    compound: PredictedCompound,
    *,
    wrapper_in_sink: bool,
    sink_names: Iterable[str],
    candidates: Optional[Sequence[PredictedCompound]] = None,
    match_type: str = "exact",
) -> SinkMatch:
    matched = tuple(candidates or (compound,))
    if compound.inchikey is None or not matched:
        raise ValueError("sink match lacks an observed InChIKey or candidate")
    if match_type not in {"exact", "stereo_missing"}:
        raise ValueError(f"unsupported accepted sink match type: {match_type}")
    if any(item.minimum_depth is None for item in matched):
        raise ValueError("sink candidate lacks minimum depth")
    representative = min(
        matched,
        key=lambda item: (
            item.minimum_depth if item.minimum_depth is not None else 10**9,
            item.compound_id,
        ),
    )
    candidate_kegg_ids = _stable_unique(
        kegg_id for item in matched for kegg_id in item.kegg_ids
    )
    observed_identity = structure_identity(compound.inchi or "")
    return SinkMatch(
        compound_id=compound.compound_id,
        inchikey=compound.inchikey,
        representative_kegg_id=representative.compound_id,
        kegg_ids=candidate_kegg_ids,
        minimum_depth=min(
            item.minimum_depth for item in matched if item.minimum_depth is not None
        ),
        wrapper_in_sink=wrapper_in_sink,
        wrapper_sink_names=_stable_unique(sink_names),
        match_type=match_type,
        stereo_review_required=match_type == "stereo_missing",
        stereo_stripped_inchikey=observed_identity.stereo_stripped_inchikey,
    )


def _empty_network(
    client_run: RetroPathClientRun,
    input_bundle: RetroPathInputBundle,
    *,
    source_in_sink: Optional[SinkMatch] = None,
) -> ParsedRetroPathNetwork:
    target = input_bundle.target_compound
    return ParsedRetroPathNetwork(
        status=client_run.result.status,
        target_compound_id=target.compound_id,
        max_steps=_max_steps(client_run),
        compounds=(
            ParsedCompoundNode(
                compound=target,
                is_target=True,
                wrapper_in_sink=source_in_sink is not None,
                wrapper_sink_names=(
                    tuple()
                    if source_in_sink is None
                    else source_in_sink.wrapper_sink_names
                ),
            ),
        ),
        transformations=tuple(),
        sink_matches=tuple() if source_in_sink is None else (source_in_sink,),
        rule_evidence=tuple(),
        rejections=tuple(),
        warnings=tuple(),
        source_in_sink=source_in_sink,
    )


def parse_retropath_network(
    client_run: RetroPathClientRun,
    input_bundle: RetroPathInputBundle,
    rules_path: str | Path,
) -> ParsedRetroPathNetwork:
    """Parse and validate one P3 run without claiming pathway completeness."""

    if not isinstance(client_run, RetroPathClientRun):
        raise ValueError("client_run must be a RetroPathClientRun")
    if not isinstance(input_bundle, RetroPathInputBundle):
        raise ValueError("input_bundle must be a RetroPathInputBundle")
    status = client_run.result.status
    if status in {"queued", "running"}:
        raise RetroPathParseError(
            "run_not_terminal",
            f"RetroPath run is still {status}",
        )
    if status in {"failed", "timed_out"}:
        raise RetroPathParseError(
            "retropath_execution_failed",
            f"RetroPath run ended with status {status}",
        )
    if status == "no_solution":
        return _empty_network(client_run, input_bundle)

    target = input_bundle.target_compound
    if target.inchikey is None:
        raise RetroPathParseError(
            "artifact_inconsistent",
            "P2 target has no full InChIKey",
        )
    target_identity = structure_identity(target.inchi or "")
    sink_by_inchikey = {
        compound.inchikey: compound
        for compound in input_bundle.sink_compounds
        if compound.inchikey is not None
    }
    sink_identities = {
        compound.compound_id: structure_identity(compound.inchi or "")
        for compound in input_bundle.sink_compounds
    }
    sink_by_stereo_stripped_inchikey_lists: dict[
        str, list[PredictedCompound]
    ] = {}
    for compound in input_bundle.sink_compounds:
        identity = sink_identities[compound.compound_id]
        sink_by_stereo_stripped_inchikey_lists.setdefault(
            identity.stereo_stripped_inchikey,
            [],
        ).append(compound)
    sink_by_stereo_stripped_inchikey = {
        key: tuple(sorted(values, key=lambda item: item.compound_id))
        for key, values in sink_by_stereo_stripped_inchikey_lists.items()
    }
    if status == "source_in_sink":
        sink = sink_by_inchikey.get(target.inchikey)
        if sink is None:
            raise RetroPathParseError(
                "sink_identity_mismatch",
                "service reported source_in_sink but target is not in the P2 sink",
            )
        match = _sink_match(
            sink,
            wrapper_in_sink=True,
            sink_names=(sink.compound_id,),
        )
        return _empty_network(client_run, input_bundle, source_in_sink=match)
    if status != "succeeded":
        raise RetroPathParseError(
            "artifact_inconsistent",
            f"unsupported RetroPath run status: {status}",
        )

    raw_dir = client_run.raw_dir.resolve()
    if not raw_dir.is_dir():
        raise RetroPathParseError(
            "artifact_missing",
            f"P3 raw artifact directory does not exist: {raw_dir}",
        )
    table_path, rows = _select_scope_table(raw_dir)
    raw_transformations, grouping_rejections = _group_rows(rows)
    json_path = _artifact_path(raw_dir, SCOPE_JSON_FILE_NAMES)
    json_topology = None if json_path is None else _json_topology(json_path)
    if json_topology is not None:
        table_by_id = {item.transformation_id: item for item in raw_transformations}
        missing_json_rows = sorted(set(json_topology) - set(table_by_id))
        if missing_json_rows:
            raise RetroPathParseError(
                "artifact_inconsistent",
                "scope JSON reactions are missing from CSV evidence: "
                + ", ".join(missing_json_rows[:10]),
            )
        raw_transformations = tuple(
            item
            for item in raw_transformations
            if item.transformation_id in json_topology
        )

    requested_rule_ids = _stable_unique(
        rule_id
        for transformation in raw_transformations
        for rule_id in transformation.rule_ids
    )
    evidence_by_id, missing_rule_ids, score_semantics = _read_rule_evidence(
        Path(rules_path).expanduser().resolve(),
        client_run.result.provenance.rules_sha256,
        requested_rule_ids,
    )
    missing_rule_set = set(missing_rule_ids)

    node_accumulators: dict[str, _NodeAccumulator] = {}
    accepted: list[ParsedTransformation] = []
    rejections = list(grouping_rejections)
    warnings: set[str] = set()

    def add_node(
        compound: PredictedCompound,
        *,
        is_target: bool,
        target_match_type: Optional[str],
        wrapper_in_sink: Optional[bool],
        sink_names: Iterable[str],
    ) -> None:
        existing = node_accumulators.get(compound.compound_id)
        if existing is None:
            existing = _NodeAccumulator(compound=compound)
            node_accumulators[compound.compound_id] = existing
        elif existing.compound.inchi != compound.inchi:
            raise RetroPathParseError(
                "artifact_inconsistent",
                f"compound identity collision for {compound.compound_id}",
            )
        existing.is_target = existing.is_target or is_target
        assert existing.wrapper_flags is not None
        assert existing.sink_names is not None
        assert existing.target_match_types is not None
        if wrapper_in_sink is not None:
            existing.wrapper_flags.add(wrapper_in_sink)
        existing.sink_names.update(sink_names)
        if target_match_type is not None:
            existing.target_match_types.add(target_match_type)

    for raw in raw_transformations:
        try:
            if any(rule_id in missing_rule_set for rule_id in raw.rule_ids):
                missing = sorted(set(raw.rule_ids) & missing_rule_set)
                raise ValueError(f"rules not found in RR02: {', '.join(missing)}")
            if any(
                evidence_by_id[item].diameter != raw.diameter for item in raw.rule_ids
            ):
                raise ValueError("reported Diameter disagrees with RR02 evidence")
            substrate_resolution = _compound_from_artifact(
                raw.substrate_inchi,
                raw.substrate_smiles,
                target=target,
                target_identity=target_identity,
                sink_by_inchikey=sink_by_inchikey,
                sink_by_stereo_stripped_inchikey=(
                    sink_by_stereo_stripped_inchikey
                ),
                sink_identities=sink_identities,
                provenance=table_path.name,
            )
            substrate = substrate_resolution.compound
            if substrate is None or substrate_resolution.auxiliary is not None:
                raise ValueError(
                    "transformation substrate is not a structural compound"
                )

            auxiliary_fragments: list[AuxiliaryFragment] = []
            product_keys: list[str] = []
            product_nodes: list[
                tuple[PredictedCompound, _RawProduct, _ArtifactCompoundResolution]
            ] = []
            product_groups: dict[
                str,
                list[
                    tuple[
                        PredictedCompound,
                        _RawProduct,
                        _ArtifactCompoundResolution,
                    ]
                ],
            ] = {}
            for raw_product in raw.products:
                product_resolution = _compound_from_artifact(
                    raw_product.inchi,
                    raw_product.smiles,
                    target=target,
                    target_identity=target_identity,
                    sink_by_inchikey=sink_by_inchikey,
                    sink_by_stereo_stripped_inchikey=(
                        sink_by_stereo_stripped_inchikey
                    ),
                    sink_identities=sink_identities,
                    provenance=table_path.name,
                )
                product = product_resolution.compound
                auxiliary = product_resolution.auxiliary
                if auxiliary is not None:
                    if raw_product.in_sink:
                        sink_names = ";".join(
                            sorted(set(raw_product.sink_names))
                        ) or "None"
                        warnings.add(
                            "auxiliary_sink_flag_ignored:"
                            f"{raw.transformation_id}:{sink_names}"
                        )
                    auxiliary_fragments.append(auxiliary)
                    continue
                if product is None or product.inchikey is None:
                    raise ValueError("product structure is incomplete")
                if raw_product.in_sink and not product_resolution.sink_candidates:
                    raise ValueError(
                        "wrapper sink flag has no accepted P2 structure match"
                    )
                observed_identity = product_resolution.observed_identity
                if observed_identity is None:
                    raise ValueError("product structure identity is incomplete")
                product_keys.append(observed_identity.full_inchikey)
                entry = (product, raw_product, product_resolution)
                product_nodes.append(entry)
                product_groups.setdefault(
                    observed_identity.stereo_stripped_inchikey,
                    [],
                ).append(entry)

            reaction_multiplicities = _reaction_product_multiplicities(
                raw.reaction_smiles
            )
            product_ids: list[str] = []
            for stereo_stripped_key, entries in sorted(product_groups.items()):
                reaction_count = reaction_multiplicities.get(stereo_stripped_key, 0)
                if reaction_count < 1:
                    raise ValueError(
                        "result product structure is absent from Reaction SMILES"
                    )
                unique_compound_ids = sorted(
                    {entry[0].compound_id for entry in entries}
                )
                if len(unique_compound_ids) == 1:
                    coefficient = max(len(entries), reaction_count)
                    product_ids.extend([unique_compound_ids[0]] * coefficient)
                else:
                    warnings.add(
                        "reaction_stereo_multiplicity_ambiguous:"
                        f"{raw.transformation_id}:{stereo_stripped_key}"
                    )
                    product_ids.extend(entry[0].compound_id for entry in entries)
            if not product_ids:
                raise ValueError("transformation has no structural product")
            if sorted(product_ids) == [substrate.compound_id]:
                raise ValueError("transformation is a structural no-op")
            normalized_product_ids = tuple(sorted(product_ids))
            normalized_auxiliary_fragments = tuple(
                sorted(
                    auxiliary_fragments,
                    key=lambda item: (item.reason_code, item.inchi, item.smiles),
                )
            )

            if json_topology is not None:
                expected_substrate_key, expected_product_keys = json_topology[
                    raw.transformation_id
                ]
                observed_substrate_key = (
                    substrate_resolution.observed_identity.full_inchikey
                    if substrate_resolution.observed_identity is not None
                    else substrate.inchikey
                )
                if (
                    observed_substrate_key != expected_substrate_key
                    or tuple(sorted(set(product_keys))) != expected_product_keys
                ):
                    raise ValueError(
                        "CSV transformation topology disagrees with scope JSON"
                    )

            variants = tuple(
                PredictedReaction.create(
                    rule_id=rule_id,
                    reaction_smiles=raw.reaction_smiles,
                    substrate_compounds=(substrate.compound_id,),
                    product_compounds=normalized_product_ids,
                    orientation="retrosynthetic",
                    source_reaction_ids=evidence_by_id[rule_id].source_reaction_ids,
                    source_ec_numbers=evidence_by_id[rule_id].ec_numbers,
                    rule_specificity=evidence_by_id[rule_id].diameter,
                    rule_specificity_semantics="diameter",
                    rule_score_raw=evidence_by_id[rule_id].score_raw,
                    score_semantics=evidence_by_id[rule_id].score_semantics,
                    cofactor_reconstruction_status=(
                        "incomplete"
                        if normalized_auxiliary_fragments
                        else "not_checked"
                    ),
                )
                for rule_id in raw.rule_ids
            )
            accepted.append(
                ParsedTransformation(
                    transformation_id=raw.transformation_id,
                    substrate_compound_id=substrate.compound_id,
                    product_compound_ids=normalized_product_ids,
                    reaction_smiles=raw.reaction_smiles,
                    iteration=raw.iteration,
                    diameter=raw.diameter,
                    score_raw=raw.score_raw,
                    score_semantics=score_semantics,
                    rule_ids=raw.rule_ids,
                    reported_ec_numbers=raw.ec_numbers,
                    reaction_variants=variants,
                    auxiliary_fragments=normalized_auxiliary_fragments,
                )
            )
            add_node(
                substrate,
                is_target=substrate_resolution.target_match_type
                in {"exact", "stereo_missing"},
                target_match_type=substrate_resolution.target_match_type,
                wrapper_in_sink=None,
                sink_names=tuple(),
            )
            for product, raw_product, product_resolution in product_nodes:
                add_node(
                    product,
                    is_target=product_resolution.target_match_type
                    in {"exact", "stereo_missing"},
                    target_match_type=product_resolution.target_match_type,
                    wrapper_in_sink=raw_product.in_sink,
                    sink_names=raw_product.sink_names,
                )
        except ValueError as exc:
            code = "artifact_inconsistent"
            detail = str(exc)
            if "no-op" in detail:
                code = "reaction_noop"
            elif "sink" in detail.lower():
                code = "sink_identity_mismatch"
            elif "rule" in detail.lower():
                code = "rule_evidence_missing"
            elif "structure" in detail.lower() or "compound" in detail.lower():
                code = "structure_invalid"
            rejections.append(
                NetworkRejection(
                    code,
                    detail,
                    transformation_id=raw.transformation_id,
                )
            )

    target_roots = [
        item
        for item in accepted
        if item.iteration == 0 and item.substrate_compound_id == target.compound_id
    ]
    if accepted and not target_roots:
        raise RetroPathParseError(
            "artifact_inconsistent",
            "scope has no iteration-0 transformation rooted at the P2 target",
        )

    nodes: list[ParsedCompoundNode] = []
    sink_matches: list[SinkMatch] = []
    for compound_id in sorted(node_accumulators):
        accumulator = node_accumulators[compound_id]
        assert accumulator.wrapper_flags is not None
        assert accumulator.sink_names is not None
        assert accumulator.target_match_types is not None
        wrapper_in_sink = True in accumulator.wrapper_flags
        if len(accumulator.wrapper_flags) > 1:
            warnings.add(f"sink_flag_disagreement:{compound_id}")
        compound = accumulator.compound
        compound_identity = structure_identity(compound.inchi or "")
        exact_sink = sink_by_inchikey.get(compound_identity.full_inchikey)
        accepted_sink_candidates: Tuple[PredictedCompound, ...] = tuple()
        sink_match_type: Optional[str] = None
        if exact_sink is not None:
            accepted_sink_candidates = (exact_sink,)
            sink_match_type = "exact"
        else:
            accepted_sink_candidates = tuple(
                item
                for item in sink_by_stereo_stripped_inchikey.get(
                    compound_identity.stereo_stripped_inchikey,
                    tuple(),
                )
                if compare_structure_identities(
                    sink_identities[item.compound_id],
                    compound_identity,
                ).match_type
                == "stereo_missing"
            )
            if accepted_sink_candidates:
                sink_match_type = "stereo_missing"
        accepted_sink = bool(accepted_sink_candidates)
        if accepted_sink and not wrapper_in_sink:
            warnings.add(f"sink_flag_disagreement:{compound_id}")
        if sink_match_type == "stereo_missing":
            warnings.add(f"stereo_missing:sink:{compound_id}")
        target_match_type = (
            "exact"
            if "exact" in accumulator.target_match_types
            else (
                "stereo_missing"
                if "stereo_missing" in accumulator.target_match_types
                else None
            )
        )
        if target_match_type == "stereo_missing":
            warnings.add(f"stereo_missing:target:{compound_id}")
        nodes.append(
            ParsedCompoundNode(
                compound=compound,
                is_target=accumulator.is_target,
                wrapper_in_sink=wrapper_in_sink,
                wrapper_sink_names=_stable_unique(accumulator.sink_names),
                target_match_type=target_match_type,
                sink_match_type=sink_match_type,
                stereo_review_required=(
                    target_match_type == "stereo_missing"
                    or sink_match_type == "stereo_missing"
                ),
            )
        )
        if accepted_sink:
            sink_matches.append(
                _sink_match(
                    compound,
                    wrapper_in_sink=wrapper_in_sink,
                    sink_names=accumulator.sink_names,
                    candidates=accepted_sink_candidates,
                    match_type=sink_match_type or "exact",
                )
            )

    if not any(node.is_target for node in nodes) and accepted:
        raise RetroPathParseError(
            "artifact_inconsistent",
            "P2 target structure is absent from the parsed scope",
        )
    return ParsedRetroPathNetwork(
        status=status,
        target_compound_id=target.compound_id,
        max_steps=_max_steps(client_run),
        compounds=tuple(nodes),
        transformations=tuple(
            sorted(accepted, key=lambda item: (item.iteration, item.transformation_id))
        ),
        sink_matches=tuple(
            sorted(
                sink_matches,
                key=lambda item: (
                    item.minimum_depth,
                    item.representative_kegg_id,
                    item.inchikey,
                ),
            )
        ),
        rule_evidence=tuple(evidence_by_id[key] for key in sorted(evidence_by_id)),
        rejections=tuple(
            sorted(
                rejections,
                key=lambda item: (
                    item.reason_code,
                    item.transformation_id or "",
                    item.compound_id or "",
                    item.reason_detail,
                ),
            )
        ),
        warnings=tuple(sorted(warnings)),
    )


__all__ = [
    "AuxiliaryFragment",
    "NetworkRejection",
    "ParsedCompoundNode",
    "ParsedRetroPathNetwork",
    "ParsedTransformation",
    "RETROPATH_RESULTS_COLUMNS",
    "RETRO_RULES_REQUIRED_COLUMNS",
    "RetroPathParseError",
    "RuleEvidence",
    "SinkMatch",
    "parse_retropath_network",
]
