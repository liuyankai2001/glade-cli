from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Literal


MG1655_EXACT_HOST = "Escherichia coli K-12 MG1655"
MG1655_INFERRED_HOST = (
    "Escherichia coli K-12 MG1655 compatible (lineage-transfer inference)"
)
BW25113_HOST = "Escherichia coli BW25113 (K-12)"
XL1_BLUE_HOST = "Escherichia coli XL1-Blue (K-12 derivative)"

V3_HOST_LABELS_BY_KEY: dict[str, tuple[str, ...]] = {
    "ecoli_k12_mg1655": (MG1655_EXACT_HOST, MG1655_INFERRED_HOST),
    "ecoli_bw25113": (BW25113_HOST,),
    "ecoli_xl1_blue": (XL1_BLUE_HOST,),
    "ecoli_o157_sakai": (),
}
CODON_TRANSFORMER_ID_KEYS = {
    50: "ecoli_o157_sakai",
    51: "ecoli_general",
    52: "ecoli_k12_mg1655",
}
INFERRED_HOST_LABELS = {MG1655_INFERRED_HOST}


class ExpressionHostContextError(ValueError):
    """Base error for the read-only CDS-to-expression host adapter."""


class HostContextMissingError(ExpressionHostContextError):
    pass


class HostContextConflictError(ExpressionHostContextError):
    pass


class HostStrainRequiredError(ExpressionHostContextError):
    pass


class NoCompatibleExpressionPartsError(ExpressionHostContextError):
    pass


@dataclass(frozen=True)
class ExpressionHostResolution:
    codon_transformer_name: str
    codon_transformer_organism_id: int | None
    requested_name: str
    expression_host_key: str
    milvus_host_labels: tuple[str, ...]
    query_name: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": "cds_selection.host",
            "codon_transformer_name": self.codon_transformer_name,
            "codon_transformer_organism_id": self.codon_transformer_organism_id,
            "requested_name": self.requested_name or None,
            "expression_host_key": self.expression_host_key,
            "milvus_host_labels": list(self.milvus_host_labels),
            "mapping_mode": "read_only_expression_adapter",
        }


def _host_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    cds_selection = manifest.get("cds_selection")
    if not isinstance(cds_selection, dict):
        return {}
    host = cds_selection.get("host")
    return host if isinstance(host, dict) else {}


def _organism_id(value: Any) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _normalized_name(value: Any) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _compact_name(value: Any) -> str:
    return _normalized_name(value).replace(" ", "")


def _identity_from_name(value: Any) -> str | None:
    normalized = _normalized_name(value)
    compact = normalized.replace(" ", "")
    if not compact:
        return None
    if "mg1655" in compact:
        return "ecoli_k12_mg1655"
    if "bw25113" in compact:
        return "ecoli_bw25113"
    if "xl1blue" in compact:
        return "ecoli_xl1_blue"
    if "o157" in compact or "sakai" in compact:
        return "ecoli_o157_sakai"
    if "escherichiacoli" in compact or compact in {"ecoli", "ecoligeneral"}:
        return "ecoli_general"
    return f"exact::{normalized}"


def _identity_species(identity_key: str | None) -> str:
    if str(identity_key or "").startswith("ecoli_"):
        return "escherichia_coli"
    return str(identity_key or "")


def _is_general_identity(identity_key: str | None) -> bool:
    return identity_key == "ecoli_general"


def _codon_transformer_name_for_id(organism_id: int | None) -> str:
    if organism_id is None:
        return ""
    try:
        from CodonTransformer.CodonUtils import ORGANISM2ID
    except Exception:
        return ""
    for name, value in ORGANISM2ID.items():
        try:
            if int(value) == organism_id:
                return str(name)
        except (TypeError, ValueError):
            continue
    return ""


def _identity_from_id(organism_id: int | None, reference_name: str) -> str | None:
    if organism_id in CODON_TRANSFORMER_ID_KEYS:
        return CODON_TRANSFORMER_ID_KEYS[int(organism_id)]
    return _identity_from_name(reference_name)


def _choose_identity(
    *,
    name: str,
    organism_id: int | None,
    id_reference_name: str,
) -> str:
    name_identity = _identity_from_name(name)
    id_identity = _identity_from_id(organism_id, id_reference_name)
    if not name_identity and not id_identity:
        raise HostContextMissingError(
            "cds_selection.host must contain the original CodonTransformer host "
            "name or codon_transformer_organism_id"
        )
    if not name_identity:
        return str(id_identity)
    if not id_identity:
        return str(name_identity)
    if name_identity == id_identity:
        return name_identity

    same_species = _identity_species(name_identity) == _identity_species(id_identity)
    if same_species and _is_general_identity(name_identity):
        return id_identity
    if same_species and _is_general_identity(id_identity):
        return name_identity
    raise HostContextConflictError(
        "cds_selection.host name and codon_transformer_organism_id identify "
        f"different hosts: name={name!r}, organism_id={organism_id!r}"
    )


def _validate_requested_host(requested_name: str, expression_host_key: str) -> None:
    requested_identity = _identity_from_name(requested_name)
    if not requested_identity:
        return
    if requested_identity == expression_host_key:
        return
    same_species = (
        _identity_species(requested_identity)
        == _identity_species(expression_host_key)
    )
    if same_species and _is_general_identity(requested_identity):
        return
    raise HostContextConflictError(
        "requested expression host conflicts with cds_selection.host; "
        f"requested={requested_name!r}, upstream_identity={expression_host_key!r}"
    )


def resolve_expression_host(
    manifest: dict[str, Any],
    requested_host: str | None = None,
) -> ExpressionHostResolution:
    """Resolve V3 lookup labels without mutating the CodonTransformer host fields."""

    host = _host_payload(manifest)
    original_name = str(host.get("name") or "").strip()
    organism_id = _organism_id(host.get("codon_transformer_organism_id"))
    id_reference_name = _codon_transformer_name_for_id(organism_id)
    expression_host_key = _choose_identity(
        name=original_name,
        organism_id=organism_id,
        id_reference_name=id_reference_name,
    )
    requested_name = str(requested_host or "").strip()
    _validate_requested_host(requested_name, expression_host_key)

    if _is_general_identity(expression_host_key):
        raise HostStrainRequiredError(
            "cds_selection.host identifies Escherichia coli general but no "
            "specific strain. Confirm a specific CodonTransformer host in the "
            "protein_to_cds stage before formal expression-part retrieval."
        )

    labels = V3_HOST_LABELS_BY_KEY.get(expression_host_key)
    if labels is None:
        exact_name = original_name or id_reference_name
        labels = (exact_name,) if exact_name else ()
    if not labels:
        raise NoCompatibleExpressionPartsError(
            "expression_parts_v3 has no compatible host labels for "
            f"{original_name or id_reference_name or expression_host_key}"
        )

    return ExpressionHostResolution(
        codon_transformer_name=original_name,
        codon_transformer_organism_id=organism_id,
        requested_name=requested_name,
        expression_host_key=expression_host_key,
        milvus_host_labels=tuple(labels),
        query_name=original_name or id_reference_name or labels[0],
    )


def describe_manifest_host(manifest: dict[str, Any]) -> dict[str, Any]:
    """Describe upstream host completeness without enforcing retrieval readiness."""

    host = _host_payload(manifest)
    original_name = str(host.get("name") or "").strip()
    organism_id = _organism_id(host.get("codon_transformer_organism_id"))
    id_reference_name = _codon_transformer_name_for_id(organism_id)
    try:
        identity_key = _choose_identity(
            name=original_name,
            organism_id=organism_id,
            id_reference_name=id_reference_name,
        )
        conflict = False
    except HostContextMissingError:
        identity_key = ""
        conflict = False
    except HostContextConflictError:
        identity_key = ""
        conflict = True
    return {
        "available": bool(original_name or organism_id is not None),
        "name": original_name,
        "codon_transformer_organism_id": organism_id,
        "expression_host_key": identity_key,
        "requires_strain_confirmation": _is_general_identity(identity_key),
        "conflict": conflict,
        "mapping_mode": "read_only_expression_adapter",
    }


def candidate_host_match_kind(
    part_hosts: list[Any],
    resolution: ExpressionHostResolution,
) -> Literal["exact", "lineage_transfer_inference", "none"]:
    values = [str(value).strip() for value in part_hosts if str(value).strip()]
    for value in values:
        if (
            value in resolution.milvus_host_labels
            and value not in INFERRED_HOST_LABELS
        ):
            return "exact"
        if (
            value not in INFERRED_HOST_LABELS
            and _identity_from_name(value) == resolution.expression_host_key
        ):
            return "exact"
    for value in values:
        if value in INFERRED_HOST_LABELS and value in resolution.milvus_host_labels:
            return "lineage_transfer_inference"
    return "none"
