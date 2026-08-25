"""RetroPath prediction models and deterministic namespaced identifiers.

The objects in this module deliberately keep predicted compounds and reactions
separate from KEGG identifiers.  They contain no RetroPath HTTP or CSV parsing
logic; later integration stages can therefore share one versioned, auditable
data contract.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Optional, Tuple, TypeAlias


RETROPATH_MODEL_SCHEMA_VERSION = 1

KEGG_COMPOUND_ID_PATTERN = re.compile(r"^C\d{5}$")
KEGG_REACTION_ID_PATTERN = re.compile(r"^R\d{5}$")
INCHIKEY_PATTERN = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
RP2_COMPOUND_ID_PATTERN = re.compile(
    r"^RP2CPD:(?:[A-Z]{14}-[A-Z]{10}-[A-Z]|[0-9a-f]{64})$"
)
RP2_REACTION_ID_PATTERN = re.compile(r"^RP2:[0-9a-f]{64}$")
RP2_ROUTE_ID_PATTERN = re.compile(r"^RP2ROUTE:[0-9a-f]{64}$")

REACTION_ORIENTATIONS = frozenset({"retrosynthetic", "biosynthetic"})
REACTION_EVIDENCE_TYPES = frozenset({"rule_predicted", "database_exact"})
RULE_SPECIFICITY_SEMANTICS = frozenset({"radius", "diameter"})
SCORE_SEMANTICS = frozenset({"higher_is_better", "lower_is_better"})
BALANCE_STATUSES = frozenset(
    {"not_checked", "balanced", "unbalanced", "incomplete"}
)
COFACTOR_RECONSTRUCTION_STATUSES = frozenset(
    {"not_checked", "complete", "incomplete", "not_applicable"}
)
ROUTE_VALIDATION_STATUSES = frozenset(
    {
        "raw",
        "structure",
        "stoichiometry",
        "gem",
        "enzyme",
        "promoted",
        "rejected",
    }
)
RETROPATH_RUN_STATUSES = frozenset(
    {
        "queued",
        "running",
        "succeeded",
        "no_solution",
        "source_in_sink",
        "failed",
        "timed_out",
    }
)
TERMINAL_RETROPATH_RUN_STATUSES = frozenset(
    {
        "succeeded",
        "no_solution",
        "source_in_sink",
        "failed",
        "timed_out",
    }
)
PROVENANCE_REQUIRED_RUN_STATUSES = frozenset(
    {"succeeded", "no_solution", "source_in_sink"}
)

JsonScalar: TypeAlias = str | int | float | bool | None


def _required(payload: Mapping[str, Any], key: str) -> Any:
    if key not in payload:
        raise ValueError(f"missing required field: {key}")
    return payload[key]


def _normalize_nonempty_string(value: Any, field_name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"{field_name} must not be empty")
    return normalized


def _normalize_optional_string(value: Any, field_name: str) -> Optional[str]:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a string or null")
    normalized = value.strip()
    return normalized or None


def _normalize_choice(value: Any, field_name: str, allowed: frozenset[str]) -> str:
    normalized = _normalize_nonempty_string(value, field_name).lower()
    if normalized not in allowed:
        choices = ", ".join(sorted(allowed))
        raise ValueError(f"{field_name} must be one of: {choices}")
    return normalized


def _normalize_optional_choice(
    value: Any,
    field_name: str,
    allowed: frozenset[str],
) -> Optional[str]:
    normalized = _normalize_optional_string(value, field_name)
    if normalized is None:
        return None
    normalized = normalized.lower()
    if normalized not in allowed:
        choices = ", ".join(sorted(allowed))
        raise ValueError(f"{field_name} must be null or one of: {choices}")
    return normalized


def _normalize_nonnegative_int(value: Any, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field_name} must be an integer")
    if value < 0:
        raise ValueError(f"{field_name} must be greater than or equal to 0")
    return value


def _normalize_optional_nonnegative_int(
    value: Any,
    field_name: str,
) -> Optional[int]:
    if value is None:
        return None
    return _normalize_nonnegative_int(value, field_name)


def _normalize_optional_finite_float(
    value: Any,
    field_name: str,
) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field_name} must be a finite number or null")
    normalized = float(value)
    if not math.isfinite(normalized):
        raise ValueError(f"{field_name} must be finite")
    return normalized


def _require_iterable(value: Any, field_name: str) -> Iterable[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise ValueError(f"{field_name} must be an iterable, not a string")
    return value


def _normalize_unique_strings(value: Any, field_name: str) -> Tuple[str, ...]:
    normalized = {
        _normalize_nonempty_string(item, field_name)
        for item in _require_iterable(value, field_name)
    }
    return tuple(sorted(normalized))


def _normalize_ordered_strings(value: Any, field_name: str) -> Tuple[str, ...]:
    return tuple(
        _normalize_nonempty_string(item, field_name)
        for item in _require_iterable(value, field_name)
    )


def _normalize_inchi(value: Any) -> str:
    inchi = _normalize_nonempty_string(value, "inchi")
    if not inchi.startswith("InChI=1S/"):
        raise ValueError("inchi must be a standard InChI beginning with 'InChI=1S/'")
    return inchi


def _normalize_inchikey(value: Any) -> Optional[str]:
    inchikey = _normalize_optional_string(value, "inchikey")
    if inchikey is None:
        return None
    inchikey = inchikey.upper()
    if not INCHIKEY_PATTERN.fullmatch(inchikey):
        raise ValueError("inchikey must use the standard 27-character InChIKey format")
    return inchikey


def _normalize_kegg_compound_id(value: Any, field_name: str) -> str:
    compound_id = _normalize_nonempty_string(value, field_name).upper()
    if not KEGG_COMPOUND_ID_PATTERN.fullmatch(compound_id):
        raise ValueError(f"{field_name} must be a KEGG Cxxxxx identifier")
    return compound_id


def _normalize_kegg_compound_ids(value: Any) -> Tuple[str, ...]:
    normalized = {
        _normalize_kegg_compound_id(item, "kegg_ids")
        for item in _require_iterable(value, "kegg_ids")
    }
    return tuple(sorted(normalized))


def _normalize_compound_id(value: Any, field_name: str = "compound_id") -> str:
    raw = _normalize_nonempty_string(value, field_name)
    upper = raw.upper()
    if KEGG_COMPOUND_ID_PATTERN.fullmatch(upper):
        return upper
    if upper.startswith("RP2CPD:"):
        suffix = raw.split(":", 1)[1]
        if INCHIKEY_PATTERN.fullmatch(suffix.upper()):
            normalized = f"RP2CPD:{suffix.upper()}"
        else:
            normalized = f"RP2CPD:{suffix.lower()}"
        if RP2_COMPOUND_ID_PATTERN.fullmatch(normalized):
            return normalized
    raise ValueError(
        f"{field_name} must be a KEGG Cxxxxx or RP2CPD:<InChIKey/SHA-256> identifier"
    )


def _normalize_reaction_id(value: Any, field_name: str = "reaction_id") -> str:
    raw = _normalize_nonempty_string(value, field_name)
    upper = raw.upper()
    if KEGG_REACTION_ID_PATTERN.fullmatch(upper):
        return upper
    if upper.startswith("RP2:"):
        normalized = f"RP2:{raw.split(':', 1)[1].lower()}"
        if RP2_REACTION_ID_PATTERN.fullmatch(normalized):
            return normalized
    raise ValueError(f"{field_name} must be a KEGG Rxxxxx or RP2:<SHA-256> identifier")


def _normalize_rp2_reaction_id(value: Any, field_name: str) -> str:
    reaction_id = _normalize_reaction_id(value, field_name)
    if not RP2_REACTION_ID_PATTERN.fullmatch(reaction_id):
        raise ValueError(f"{field_name} must be an RP2:<SHA-256> identifier")
    return reaction_id


def _normalize_candidate_id(value: Any) -> str:
    raw = _normalize_nonempty_string(value, "candidate_id")
    if raw.upper().startswith("RP2ROUTE:"):
        raw = f"RP2ROUTE:{raw.split(':', 1)[1].lower()}"
    if not RP2_ROUTE_ID_PATTERN.fullmatch(raw):
        raise ValueError("candidate_id must be an RP2ROUTE:<SHA-256> identifier")
    return raw


def _normalize_compound_multiset(value: Any, field_name: str) -> Tuple[str, ...]:
    normalized = [
        _normalize_compound_id(item, field_name)
        for item in _require_iterable(value, field_name)
    ]
    if not normalized:
        raise ValueError(f"{field_name} must contain at least one compound")
    return tuple(sorted(normalized))


def _stable_sha256(payload: Mapping[str, Any]) -> str:
    canonical_json = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()


def _predicted_compound_id(inchi: str, inchikey: Optional[str]) -> str:
    if inchikey is not None:
        return f"RP2CPD:{inchikey}"
    digest = _stable_sha256(
        {
            "entity_type": "compound",
            "inchi": inchi,
            "schema_version": RETROPATH_MODEL_SCHEMA_VERSION,
        }
    )
    return f"RP2CPD:{digest}"


def _predicted_reaction_id(
    *,
    rule_id: str,
    reaction_smiles: str,
    substrate_compounds: Tuple[str, ...],
    product_compounds: Tuple[str, ...],
    orientation: str,
) -> str:
    digest = _stable_sha256(
        {
            "entity_type": "reaction",
            "orientation": orientation,
            "product_compounds": list(product_compounds),
            "reaction_smiles": reaction_smiles,
            "rule_id": rule_id,
            "schema_version": RETROPATH_MODEL_SCHEMA_VERSION,
            "substrate_compounds": list(substrate_compounds),
        }
    )
    return f"RP2:{digest}"


def _candidate_route_id(
    *,
    target_compound_id: str,
    matched_sink_kegg_id: str,
    matched_sink_depth: int,
    kegg_prefix_reaction_ids: Tuple[str, ...],
    retropath_reaction_ids: Tuple[str, ...],
) -> str:
    digest = _stable_sha256(
        {
            "entity_type": "candidate_route",
            "kegg_prefix_reaction_ids": list(kegg_prefix_reaction_ids),
            "matched_sink_depth": matched_sink_depth,
            "matched_sink_kegg_id": matched_sink_kegg_id,
            "retropath_reaction_ids": list(retropath_reaction_ids),
            "schema_version": RETROPATH_MODEL_SCHEMA_VERSION,
            "target_compound_id": target_compound_id,
        }
    )
    return f"RP2ROUTE:{digest}"


@dataclass(frozen=True)
class PredictedCompound:
    """One structure-addressable KEGG or RetroPath compound."""

    compound_id: str
    inchi: str
    name: str = ""
    inchikey: Optional[str] = None
    isomeric_smiles: Optional[str] = None
    formula: Optional[str] = None
    charge: Optional[int] = None
    kegg_ids: Tuple[str, ...] = tuple()
    minimum_depth: Optional[int] = None
    structure_provenance: Tuple[str, ...] = tuple()

    def __post_init__(self) -> None:
        compound_id = _normalize_compound_id(self.compound_id)
        inchi = _normalize_inchi(self.inchi)
        inchikey = _normalize_inchikey(self.inchikey)
        kegg_ids = _normalize_kegg_compound_ids(self.kegg_ids)
        name = _normalize_optional_string(self.name, "name") or ""
        isomeric_smiles = _normalize_optional_string(
            self.isomeric_smiles,
            "isomeric_smiles",
        )
        formula = _normalize_optional_string(self.formula, "formula")
        charge = self.charge
        if charge is not None and (isinstance(charge, bool) or not isinstance(charge, int)):
            raise ValueError("charge must be an integer or null")
        minimum_depth = _normalize_optional_nonnegative_int(
            self.minimum_depth,
            "minimum_depth",
        )
        provenance = _normalize_unique_strings(
            self.structure_provenance,
            "structure_provenance",
        )

        if KEGG_COMPOUND_ID_PATTERN.fullmatch(compound_id):
            if compound_id not in kegg_ids:
                raise ValueError("a KEGG compound_id must also be present in kegg_ids")
        else:
            if kegg_ids:
                raise ValueError(
                    "a compound with KEGG mappings must use a KEGG Cxxxxx compound_id"
                )
            expected_id = _predicted_compound_id(inchi, inchikey)
            if compound_id != expected_id:
                raise ValueError(
                    "compound_id does not match the canonical InChI/InChIKey identity"
                )

        object.__setattr__(self, "compound_id", compound_id)
        object.__setattr__(self, "inchi", inchi)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "inchikey", inchikey)
        object.__setattr__(self, "isomeric_smiles", isomeric_smiles)
        object.__setattr__(self, "formula", formula)
        object.__setattr__(self, "charge", charge)
        object.__setattr__(self, "kegg_ids", kegg_ids)
        object.__setattr__(self, "minimum_depth", minimum_depth)
        object.__setattr__(self, "structure_provenance", provenance)

    @classmethod
    def create(
        cls,
        *,
        inchi: str,
        compound_id: Optional[str] = None,
        name: str = "",
        inchikey: Optional[str] = None,
        isomeric_smiles: Optional[str] = None,
        formula: Optional[str] = None,
        charge: Optional[int] = None,
        kegg_ids: Iterable[str] = tuple(),
        minimum_depth: Optional[int] = None,
        structure_provenance: Iterable[str] = tuple(),
    ) -> "PredictedCompound":
        normalized_inchi = _normalize_inchi(inchi)
        normalized_inchikey = _normalize_inchikey(inchikey)
        normalized_kegg_ids = _normalize_kegg_compound_ids(kegg_ids)
        if compound_id is None:
            if normalized_kegg_ids:
                compound_id = normalized_kegg_ids[0]
            else:
                compound_id = _predicted_compound_id(
                    normalized_inchi,
                    normalized_inchikey,
                )
        normalized_compound_id = _normalize_compound_id(compound_id)
        if KEGG_COMPOUND_ID_PATTERN.fullmatch(normalized_compound_id):
            normalized_kegg_ids = tuple(
                sorted({*normalized_kegg_ids, normalized_compound_id})
            )
        return cls(
            compound_id=normalized_compound_id,
            inchi=normalized_inchi,
            name=name,
            inchikey=normalized_inchikey,
            isomeric_smiles=isomeric_smiles,
            formula=formula,
            charge=charge,
            kegg_ids=normalized_kegg_ids,
            minimum_depth=minimum_depth,
            structure_provenance=tuple(structure_provenance),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "compound_id": self.compound_id,
            "name": self.name,
            "inchi": self.inchi,
            "inchikey": self.inchikey,
            "isomeric_smiles": self.isomeric_smiles,
            "formula": self.formula,
            "charge": self.charge,
            "kegg_ids": list(self.kegg_ids),
            "minimum_depth": self.minimum_depth,
            "structure_provenance": list(self.structure_provenance),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "PredictedCompound":
        if not isinstance(payload, Mapping):
            raise ValueError("PredictedCompound payload must be an object")
        return cls(
            compound_id=_required(payload, "compound_id"),
            inchi=_required(payload, "inchi"),
            name=payload.get("name", ""),
            inchikey=payload.get("inchikey"),
            isomeric_smiles=payload.get("isomeric_smiles"),
            formula=payload.get("formula"),
            charge=payload.get("charge"),
            kegg_ids=payload.get("kegg_ids", tuple()),
            minimum_depth=payload.get("minimum_depth"),
            structure_provenance=payload.get("structure_provenance", tuple()),
        )


@dataclass(frozen=True)
class PredictedReaction:
    """One explicitly directed reaction predicted from a RetroRules rule."""

    reaction_id: str
    rule_id: str
    reaction_smiles: str
    substrate_compounds: Tuple[str, ...]
    product_compounds: Tuple[str, ...]
    orientation: str
    reaction_source: str = "retropath"
    evidence_type: str = "rule_predicted"
    source_reaction_ids: Tuple[str, ...] = tuple()
    source_ec_numbers: Tuple[str, ...] = tuple()
    source_uniprot_ids: Tuple[str, ...] = tuple()
    rule_specificity: Optional[int] = None
    rule_specificity_semantics: Optional[str] = None
    rule_score_raw: Optional[float] = None
    score_semantics: Optional[str] = None
    balance_status: str = "not_checked"
    cofactor_reconstruction_status: str = "not_checked"

    def __post_init__(self) -> None:
        reaction_id = _normalize_rp2_reaction_id(self.reaction_id, "reaction_id")
        rule_id = _normalize_nonempty_string(self.rule_id, "rule_id")
        reaction_smiles = _normalize_nonempty_string(
            self.reaction_smiles,
            "reaction_smiles",
        )
        substrates = _normalize_compound_multiset(
            self.substrate_compounds,
            "substrate_compounds",
        )
        products = _normalize_compound_multiset(
            self.product_compounds,
            "product_compounds",
        )
        orientation = _normalize_choice(
            self.orientation,
            "orientation",
            REACTION_ORIENTATIONS,
        )
        reaction_source = _normalize_nonempty_string(
            self.reaction_source,
            "reaction_source",
        ).lower()
        if reaction_source != "retropath":
            raise ValueError("reaction_source must be 'retropath'")
        evidence_type = _normalize_choice(
            self.evidence_type,
            "evidence_type",
            REACTION_EVIDENCE_TYPES,
        )
        source_reaction_ids = _normalize_unique_strings(
            self.source_reaction_ids,
            "source_reaction_ids",
        )
        source_ec_numbers = _normalize_unique_strings(
            self.source_ec_numbers,
            "source_ec_numbers",
        )
        source_uniprot_ids = _normalize_unique_strings(
            self.source_uniprot_ids,
            "source_uniprot_ids",
        )
        specificity = _normalize_optional_nonnegative_int(
            self.rule_specificity,
            "rule_specificity",
        )
        specificity_semantics = _normalize_optional_choice(
            self.rule_specificity_semantics,
            "rule_specificity_semantics",
            RULE_SPECIFICITY_SEMANTICS,
        )
        if (specificity is None) != (specificity_semantics is None):
            raise ValueError(
                "rule_specificity and rule_specificity_semantics must be set together"
            )
        score = _normalize_optional_finite_float(
            self.rule_score_raw,
            "rule_score_raw",
        )
        score_semantics = _normalize_optional_choice(
            self.score_semantics,
            "score_semantics",
            SCORE_SEMANTICS,
        )
        if (score is None) != (score_semantics is None):
            raise ValueError("rule_score_raw and score_semantics must be set together")
        balance_status = _normalize_choice(
            self.balance_status,
            "balance_status",
            BALANCE_STATUSES,
        )
        cofactor_status = _normalize_choice(
            self.cofactor_reconstruction_status,
            "cofactor_reconstruction_status",
            COFACTOR_RECONSTRUCTION_STATUSES,
        )
        expected_id = _predicted_reaction_id(
            rule_id=rule_id,
            reaction_smiles=reaction_smiles,
            substrate_compounds=substrates,
            product_compounds=products,
            orientation=orientation,
        )
        if reaction_id != expected_id:
            raise ValueError("reaction_id does not match the canonical reaction identity")

        object.__setattr__(self, "reaction_id", reaction_id)
        object.__setattr__(self, "rule_id", rule_id)
        object.__setattr__(self, "reaction_smiles", reaction_smiles)
        object.__setattr__(self, "substrate_compounds", substrates)
        object.__setattr__(self, "product_compounds", products)
        object.__setattr__(self, "orientation", orientation)
        object.__setattr__(self, "reaction_source", reaction_source)
        object.__setattr__(self, "evidence_type", evidence_type)
        object.__setattr__(self, "source_reaction_ids", source_reaction_ids)
        object.__setattr__(self, "source_ec_numbers", source_ec_numbers)
        object.__setattr__(self, "source_uniprot_ids", source_uniprot_ids)
        object.__setattr__(self, "rule_specificity", specificity)
        object.__setattr__(
            self,
            "rule_specificity_semantics",
            specificity_semantics,
        )
        object.__setattr__(self, "rule_score_raw", score)
        object.__setattr__(self, "score_semantics", score_semantics)
        object.__setattr__(self, "balance_status", balance_status)
        object.__setattr__(self, "cofactor_reconstruction_status", cofactor_status)

    @classmethod
    def create(
        cls,
        *,
        rule_id: str,
        reaction_smiles: str,
        substrate_compounds: Iterable[str],
        product_compounds: Iterable[str],
        orientation: str,
        evidence_type: str = "rule_predicted",
        source_reaction_ids: Iterable[str] = tuple(),
        source_ec_numbers: Iterable[str] = tuple(),
        source_uniprot_ids: Iterable[str] = tuple(),
        rule_specificity: Optional[int] = None,
        rule_specificity_semantics: Optional[str] = None,
        rule_score_raw: Optional[float] = None,
        score_semantics: Optional[str] = None,
        balance_status: str = "not_checked",
        cofactor_reconstruction_status: str = "not_checked",
    ) -> "PredictedReaction":
        normalized_rule_id = _normalize_nonempty_string(rule_id, "rule_id")
        normalized_reaction_smiles = _normalize_nonempty_string(
            reaction_smiles,
            "reaction_smiles",
        )
        normalized_substrates = _normalize_compound_multiset(
            substrate_compounds,
            "substrate_compounds",
        )
        normalized_products = _normalize_compound_multiset(
            product_compounds,
            "product_compounds",
        )
        normalized_orientation = _normalize_choice(
            orientation,
            "orientation",
            REACTION_ORIENTATIONS,
        )
        reaction_id = _predicted_reaction_id(
            rule_id=normalized_rule_id,
            reaction_smiles=normalized_reaction_smiles,
            substrate_compounds=normalized_substrates,
            product_compounds=normalized_products,
            orientation=normalized_orientation,
        )
        return cls(
            reaction_id=reaction_id,
            rule_id=normalized_rule_id,
            reaction_smiles=normalized_reaction_smiles,
            substrate_compounds=normalized_substrates,
            product_compounds=normalized_products,
            orientation=normalized_orientation,
            evidence_type=evidence_type,
            source_reaction_ids=tuple(source_reaction_ids),
            source_ec_numbers=tuple(source_ec_numbers),
            source_uniprot_ids=tuple(source_uniprot_ids),
            rule_specificity=rule_specificity,
            rule_specificity_semantics=rule_specificity_semantics,
            rule_score_raw=rule_score_raw,
            score_semantics=score_semantics,
            balance_status=balance_status,
            cofactor_reconstruction_status=cofactor_reconstruction_status,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "reaction_id": self.reaction_id,
            "reaction_source": self.reaction_source,
            "evidence_type": self.evidence_type,
            "rule_id": self.rule_id,
            "reaction_smiles": self.reaction_smiles,
            "substrate_compounds": list(self.substrate_compounds),
            "product_compounds": list(self.product_compounds),
            "orientation": self.orientation,
            "source_reaction_ids": list(self.source_reaction_ids),
            "source_ec_numbers": list(self.source_ec_numbers),
            "source_uniprot_ids": list(self.source_uniprot_ids),
            "rule_specificity": self.rule_specificity,
            "rule_specificity_semantics": self.rule_specificity_semantics,
            "rule_score_raw": self.rule_score_raw,
            "score_semantics": self.score_semantics,
            "balance_status": self.balance_status,
            "cofactor_reconstruction_status": self.cofactor_reconstruction_status,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "PredictedReaction":
        if not isinstance(payload, Mapping):
            raise ValueError("PredictedReaction payload must be an object")
        return cls(
            reaction_id=_required(payload, "reaction_id"),
            reaction_source=payload.get("reaction_source", "retropath"),
            evidence_type=payload.get("evidence_type", "rule_predicted"),
            rule_id=_required(payload, "rule_id"),
            reaction_smiles=_required(payload, "reaction_smiles"),
            substrate_compounds=_required(payload, "substrate_compounds"),
            product_compounds=_required(payload, "product_compounds"),
            orientation=_required(payload, "orientation"),
            source_reaction_ids=payload.get("source_reaction_ids", tuple()),
            source_ec_numbers=payload.get("source_ec_numbers", tuple()),
            source_uniprot_ids=payload.get("source_uniprot_ids", tuple()),
            rule_specificity=payload.get("rule_specificity"),
            rule_specificity_semantics=payload.get("rule_specificity_semantics"),
            rule_score_raw=payload.get("rule_score_raw"),
            score_semantics=payload.get("score_semantics"),
            balance_status=payload.get("balance_status", "not_checked"),
            cofactor_reconstruction_status=payload.get(
                "cofactor_reconstruction_status",
                "not_checked",
            ),
        )


@dataclass(frozen=True)
class CandidateRoute:
    """A KEGG-prefix plus RetroPath-suffix candidate route."""

    candidate_id: str
    target_compound_id: str
    matched_sink_kegg_id: str
    matched_sink_depth: int
    kegg_prefix_reaction_ids: Tuple[str, ...]
    retropath_reaction_ids: Tuple[str, ...]
    minimum_rule_specificity: Optional[int] = None
    validation_status: str = "raw"
    review_required: bool = True
    rejection_reasons: Tuple[str, ...] = tuple()

    def __post_init__(self) -> None:
        candidate_id = _normalize_candidate_id(self.candidate_id)
        target_compound_id = _normalize_compound_id(
            self.target_compound_id,
            "target_compound_id",
        )
        matched_sink_kegg_id = _normalize_kegg_compound_id(
            self.matched_sink_kegg_id,
            "matched_sink_kegg_id",
        )
        matched_sink_depth = _normalize_nonnegative_int(
            self.matched_sink_depth,
            "matched_sink_depth",
        )
        kegg_prefix = tuple(
            _normalize_reaction_id(item, "kegg_prefix_reaction_ids")
            for item in _require_iterable(
                self.kegg_prefix_reaction_ids,
                "kegg_prefix_reaction_ids",
            )
        )
        if any(
            not KEGG_REACTION_ID_PATTERN.fullmatch(reaction_id)
            for reaction_id in kegg_prefix
        ):
            raise ValueError("kegg_prefix_reaction_ids may only contain KEGG Rxxxxx IDs")
        retropath_suffix = tuple(
            _normalize_rp2_reaction_id(item, "retropath_reaction_ids")
            for item in _require_iterable(
                self.retropath_reaction_ids,
                "retropath_reaction_ids",
            )
        )
        if not retropath_suffix:
            raise ValueError("retropath_reaction_ids must contain at least one reaction")
        minimum_specificity = _normalize_optional_nonnegative_int(
            self.minimum_rule_specificity,
            "minimum_rule_specificity",
        )
        validation_status = _normalize_choice(
            self.validation_status,
            "validation_status",
            ROUTE_VALIDATION_STATUSES,
        )
        if not isinstance(self.review_required, bool):
            raise ValueError("review_required must be a boolean")
        if not self.review_required:
            raise ValueError("a route containing RP2 steps must require manual review")
        rejection_reasons = _normalize_unique_strings(
            self.rejection_reasons,
            "rejection_reasons",
        )
        expected_id = _candidate_route_id(
            target_compound_id=target_compound_id,
            matched_sink_kegg_id=matched_sink_kegg_id,
            matched_sink_depth=matched_sink_depth,
            kegg_prefix_reaction_ids=kegg_prefix,
            retropath_reaction_ids=retropath_suffix,
        )
        if candidate_id != expected_id:
            raise ValueError("candidate_id does not match the canonical route identity")

        object.__setattr__(self, "candidate_id", candidate_id)
        object.__setattr__(self, "target_compound_id", target_compound_id)
        object.__setattr__(self, "matched_sink_kegg_id", matched_sink_kegg_id)
        object.__setattr__(self, "matched_sink_depth", matched_sink_depth)
        object.__setattr__(self, "kegg_prefix_reaction_ids", kegg_prefix)
        object.__setattr__(self, "retropath_reaction_ids", retropath_suffix)
        object.__setattr__(self, "minimum_rule_specificity", minimum_specificity)
        object.__setattr__(self, "validation_status", validation_status)
        object.__setattr__(self, "review_required", self.review_required)
        object.__setattr__(self, "rejection_reasons", rejection_reasons)

    @property
    def kegg_prefix_steps(self) -> int:
        return len(self.kegg_prefix_reaction_ids)

    @property
    def retropath_steps(self) -> int:
        return len(self.retropath_reaction_ids)

    @property
    def total_steps(self) -> int:
        return self.kegg_prefix_steps + self.retropath_steps

    @property
    def route_source(self) -> str:
        return "kegg_retropath"

    @property
    def contains_predicted_steps(self) -> bool:
        return True

    @classmethod
    def create(
        cls,
        *,
        target_compound_id: str,
        matched_sink_kegg_id: str,
        matched_sink_depth: int,
        kegg_prefix_reaction_ids: Iterable[str] = tuple(),
        retropath_reaction_ids: Iterable[str],
        minimum_rule_specificity: Optional[int] = None,
        validation_status: str = "raw",
        review_required: bool = True,
        rejection_reasons: Iterable[str] = tuple(),
    ) -> "CandidateRoute":
        normalized_target = _normalize_compound_id(
            target_compound_id,
            "target_compound_id",
        )
        normalized_sink = _normalize_kegg_compound_id(
            matched_sink_kegg_id,
            "matched_sink_kegg_id",
        )
        normalized_depth = _normalize_nonnegative_int(
            matched_sink_depth,
            "matched_sink_depth",
        )
        normalized_prefix = tuple(
            _normalize_reaction_id(item, "kegg_prefix_reaction_ids")
            for item in _require_iterable(
                kegg_prefix_reaction_ids,
                "kegg_prefix_reaction_ids",
            )
        )
        if any(
            not KEGG_REACTION_ID_PATTERN.fullmatch(reaction_id)
            for reaction_id in normalized_prefix
        ):
            raise ValueError("kegg_prefix_reaction_ids may only contain KEGG Rxxxxx IDs")
        normalized_suffix = tuple(
            _normalize_rp2_reaction_id(item, "retropath_reaction_ids")
            for item in _require_iterable(
                retropath_reaction_ids,
                "retropath_reaction_ids",
            )
        )
        if not normalized_suffix:
            raise ValueError("retropath_reaction_ids must contain at least one reaction")
        candidate_id = _candidate_route_id(
            target_compound_id=normalized_target,
            matched_sink_kegg_id=normalized_sink,
            matched_sink_depth=normalized_depth,
            kegg_prefix_reaction_ids=normalized_prefix,
            retropath_reaction_ids=normalized_suffix,
        )
        return cls(
            candidate_id=candidate_id,
            target_compound_id=normalized_target,
            matched_sink_kegg_id=normalized_sink,
            matched_sink_depth=normalized_depth,
            kegg_prefix_reaction_ids=normalized_prefix,
            retropath_reaction_ids=normalized_suffix,
            minimum_rule_specificity=minimum_rule_specificity,
            validation_status=validation_status,
            review_required=review_required,
            rejection_reasons=tuple(rejection_reasons),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "target_compound_id": self.target_compound_id,
            "matched_sink_kegg_id": self.matched_sink_kegg_id,
            "matched_sink_depth": self.matched_sink_depth,
            "kegg_prefix_reaction_ids": list(self.kegg_prefix_reaction_ids),
            "retropath_reaction_ids": list(self.retropath_reaction_ids),
            "kegg_prefix_steps": self.kegg_prefix_steps,
            "retropath_steps": self.retropath_steps,
            "total_steps": self.total_steps,
            "route_source": self.route_source,
            "contains_predicted_steps": self.contains_predicted_steps,
            "minimum_rule_specificity": self.minimum_rule_specificity,
            "validation_status": self.validation_status,
            "review_required": self.review_required,
            "rejection_reasons": list(self.rejection_reasons),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "CandidateRoute":
        if not isinstance(payload, Mapping):
            raise ValueError("CandidateRoute payload must be an object")
        route = cls(
            candidate_id=_required(payload, "candidate_id"),
            target_compound_id=_required(payload, "target_compound_id"),
            matched_sink_kegg_id=_required(payload, "matched_sink_kegg_id"),
            matched_sink_depth=_required(payload, "matched_sink_depth"),
            kegg_prefix_reaction_ids=payload.get(
                "kegg_prefix_reaction_ids",
                tuple(),
            ),
            retropath_reaction_ids=_required(payload, "retropath_reaction_ids"),
            minimum_rule_specificity=payload.get("minimum_rule_specificity"),
            validation_status=payload.get("validation_status", "raw"),
            review_required=payload.get("review_required", True),
            rejection_reasons=payload.get("rejection_reasons", tuple()),
        )
        derived_values = {
            "kegg_prefix_steps": route.kegg_prefix_steps,
            "retropath_steps": route.retropath_steps,
            "total_steps": route.total_steps,
            "route_source": route.route_source,
            "contains_predicted_steps": route.contains_predicted_steps,
        }
        for field_name, expected in derived_values.items():
            if field_name in payload and (
                type(payload[field_name]) is not type(expected)
                or payload[field_name] != expected
            ):
                raise ValueError(f"{field_name} does not match the route contents")
        return route


@dataclass(frozen=True)
class RetroPathRuntimeProvenance:
    """Pinned versions and rule checksum used by one RetroPath execution."""

    wrapper_version: str
    workflow_version: str
    knime_version: str
    rdkit_plugin_version: str
    rules_version: str
    rules_sha256: str
    wrapper_reported_version: Optional[str] = None

    def __post_init__(self) -> None:
        for field_name in (
            "wrapper_version",
            "workflow_version",
            "knime_version",
            "rdkit_plugin_version",
            "rules_version",
        ):
            object.__setattr__(
                self,
                field_name,
                _normalize_nonempty_string(getattr(self, field_name), field_name),
            )
        rules_sha256 = _normalize_nonempty_string(
            self.rules_sha256,
            "rules_sha256",
        ).lower()
        if not SHA256_PATTERN.fullmatch(rules_sha256):
            raise ValueError("rules_sha256 must be a full lowercase SHA-256 digest")
        object.__setattr__(self, "rules_sha256", rules_sha256)
        object.__setattr__(
            self,
            "wrapper_reported_version",
            _normalize_optional_string(
                self.wrapper_reported_version,
                "wrapper_reported_version",
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "wrapper_version": self.wrapper_version,
            "wrapper_reported_version": self.wrapper_reported_version,
            "workflow_version": self.workflow_version,
            "knime_version": self.knime_version,
            "rdkit_plugin_version": self.rdkit_plugin_version,
            "rules_version": self.rules_version,
            "rules_sha256": self.rules_sha256,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "RetroPathRuntimeProvenance":
        if not isinstance(payload, Mapping):
            raise ValueError("RetroPathRuntimeProvenance payload must be an object")
        return cls(
            wrapper_version=_required(payload, "wrapper_version"),
            wrapper_reported_version=payload.get("wrapper_reported_version"),
            workflow_version=_required(payload, "workflow_version"),
            knime_version=_required(payload, "knime_version"),
            rdkit_plugin_version=_required(payload, "rdkit_plugin_version"),
            rules_version=_required(payload, "rules_version"),
            rules_sha256=_required(payload, "rules_sha256"),
        )


def _normalize_parameters(value: Any) -> Tuple[Tuple[str, JsonScalar], ...]:
    if isinstance(value, Mapping):
        items = value.items()
    else:
        items = _require_iterable(value, "parameters")
    normalized: dict[str, JsonScalar] = {}
    for item in items:
        if isinstance(item, (str, bytes)) or not isinstance(item, Iterable):
            raise ValueError("each parameters entry must be a key/value pair")
        pair = tuple(item)
        if len(pair) != 2:
            raise ValueError("each parameters entry must be a key/value pair")
        key = _normalize_nonempty_string(pair[0], "parameters key")
        if key in normalized:
            raise ValueError(f"duplicate parameters key: {key}")
        parameter_value = pair[1]
        if not (
            parameter_value is None
            or isinstance(parameter_value, (str, int, float, bool))
        ):
            raise ValueError("parameters values must be JSON scalar values")
        if isinstance(parameter_value, float) and not math.isfinite(parameter_value):
            raise ValueError("parameters values must be finite")
        normalized[key] = parameter_value
    return tuple(sorted(normalized.items()))


@dataclass(frozen=True)
class RetroPathRunResult:
    """Versioned result envelope shared by the future client and parser."""

    job_id: str
    status: str
    return_code: Optional[int] = None
    provenance: Optional[RetroPathRuntimeProvenance] = None
    parameters: Tuple[Tuple[str, JsonScalar], ...] = tuple()
    artifacts: Tuple[str, ...] = tuple()
    compounds: Tuple[PredictedCompound, ...] = tuple()
    reactions: Tuple[PredictedReaction, ...] = tuple()
    candidate_routes: Tuple[CandidateRoute, ...] = tuple()
    errors: Tuple[str, ...] = tuple()
    schema_version: int = RETROPATH_MODEL_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != RETROPATH_MODEL_SCHEMA_VERSION
        ):
            raise ValueError(
                "unsupported RetroPath model schema_version: "
                f"{self.schema_version!r}"
            )
        job_id = _normalize_nonempty_string(self.job_id, "job_id")
        status = _normalize_choice(
            self.status,
            "status",
            RETROPATH_RUN_STATUSES,
        )
        return_code = self.return_code
        if return_code is not None and (
            isinstance(return_code, bool) or not isinstance(return_code, int)
        ):
            raise ValueError("return_code must be an integer or null")
        provenance = self.provenance
        if provenance is not None and not isinstance(
            provenance,
            RetroPathRuntimeProvenance,
        ):
            raise ValueError("provenance must be RetroPathRuntimeProvenance or null")
        if status in PROVENANCE_REQUIRED_RUN_STATUSES and provenance is None:
            raise ValueError(f"provenance is required when status is {status!r}")
        parameters = _normalize_parameters(self.parameters)
        artifacts = _normalize_unique_strings(self.artifacts, "artifacts")
        compounds = tuple(self.compounds)
        reactions = tuple(self.reactions)
        candidate_routes = tuple(self.candidate_routes)
        if not all(isinstance(item, PredictedCompound) for item in compounds):
            raise ValueError("compounds must contain only PredictedCompound objects")
        if not all(isinstance(item, PredictedReaction) for item in reactions):
            raise ValueError("reactions must contain only PredictedReaction objects")
        if not all(isinstance(item, CandidateRoute) for item in candidate_routes):
            raise ValueError("candidate_routes must contain only CandidateRoute objects")
        errors = _normalize_ordered_strings(self.errors, "errors")

        object.__setattr__(self, "job_id", job_id)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "return_code", return_code)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "artifacts", artifacts)
        object.__setattr__(self, "compounds", compounds)
        object.__setattr__(self, "reactions", reactions)
        object.__setattr__(self, "candidate_routes", candidate_routes)
        object.__setattr__(self, "errors", errors)

    @property
    def is_terminal(self) -> bool:
        return self.status in TERMINAL_RETROPATH_RUN_STATUSES

    @property
    def has_candidates(self) -> bool:
        return bool(self.candidate_routes)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "job_id": self.job_id,
            "status": self.status,
            "return_code": self.return_code,
            "provenance": (
                self.provenance.to_dict() if self.provenance is not None else None
            ),
            "parameters": dict(self.parameters),
            "artifacts": list(self.artifacts),
            "compounds": [item.to_dict() for item in self.compounds],
            "reactions": [item.to_dict() for item in self.reactions],
            "candidate_routes": [item.to_dict() for item in self.candidate_routes],
            "errors": list(self.errors),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "RetroPathRunResult":
        if not isinstance(payload, Mapping):
            raise ValueError("RetroPathRunResult payload must be an object")
        schema_version = _required(payload, "schema_version")
        provenance_payload = payload.get("provenance")
        provenance = (
            None
            if provenance_payload is None
            else RetroPathRuntimeProvenance.from_dict(provenance_payload)
        )
        compounds_payload = payload.get("compounds", tuple())
        reactions_payload = payload.get("reactions", tuple())
        routes_payload = payload.get("candidate_routes", tuple())
        return cls(
            schema_version=schema_version,
            job_id=_required(payload, "job_id"),
            status=_required(payload, "status"),
            return_code=payload.get("return_code"),
            provenance=provenance,
            parameters=payload.get("parameters", {}),
            artifacts=payload.get("artifacts", tuple()),
            compounds=tuple(
                PredictedCompound.from_dict(item)
                for item in _require_iterable(compounds_payload, "compounds")
            ),
            reactions=tuple(
                PredictedReaction.from_dict(item)
                for item in _require_iterable(reactions_payload, "reactions")
            ),
            candidate_routes=tuple(
                CandidateRoute.from_dict(item)
                for item in _require_iterable(routes_payload, "candidate_routes")
            ),
            errors=payload.get("errors", tuple()),
        )


def retropath_result_to_json(
    result: RetroPathRunResult,
    *,
    indent: Optional[int] = 2,
) -> str:
    """Serialize a run result to deterministic, UTF-8-friendly JSON text."""

    if not isinstance(result, RetroPathRunResult):
        raise ValueError("result must be a RetroPathRunResult")
    return json.dumps(
        result.to_dict(),
        ensure_ascii=False,
        sort_keys=True,
        indent=indent,
        allow_nan=False,
    )


def retropath_result_from_json(text: str) -> RetroPathRunResult:
    """Deserialize and validate a versioned RetroPath run result."""

    if not isinstance(text, str):
        raise ValueError("RetroPath result JSON must be text")
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid RetroPath result JSON: {exc.msg}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("RetroPath result JSON root must be an object")
    return RetroPathRunResult.from_dict(payload)


__all__ = [
    "BALANCE_STATUSES",
    "COFACTOR_RECONSTRUCTION_STATUSES",
    "CandidateRoute",
    "PredictedCompound",
    "PredictedReaction",
    "REACTION_EVIDENCE_TYPES",
    "REACTION_ORIENTATIONS",
    "RETROPATH_MODEL_SCHEMA_VERSION",
    "RETROPATH_RUN_STATUSES",
    "ROUTE_VALIDATION_STATUSES",
    "RetroPathRunResult",
    "RetroPathRuntimeProvenance",
    "RULE_SPECIFICITY_SEMANTICS",
    "SCORE_SEMANTICS",
    "retropath_result_from_json",
    "retropath_result_to_json",
]
