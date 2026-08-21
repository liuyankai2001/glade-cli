"""Deterministic auxiliary-role annotation for main-enzyme candidates."""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from typing import Any


P450_DOMAIN_IDS = frozenset({"PF00067", "IPR001128"})
P450_REDUCTASE_DOMAIN_IDS = frozenset({"PF00258", "PF00667", "PF00175"})
DOMAIN_NAMESPACES = frozenset({"pfam", "interpro"})


def _values(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        raw = [str(item).strip() for item in value]
    else:
        raw = [
            item.strip()
            for item in re.split(r"\s*[;|]\s*", str(value or ""))
        ]
    return list(dict.fromkeys(item for item in raw if item))


def _requirements(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        raw = value
    else:
        text = str(value or "").strip()
        if not text:
            return []
        try:
            raw = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError("auxiliary requirements are not valid JSON") from exc
    if not isinstance(raw, list) or any(not isinstance(item, Mapping) for item in raw):
        raise ValueError("auxiliary requirements must be a list of objects")
    normalized: list[dict[str, Any]] = []
    for item in raw:
        role = str(item.get("role") or "").strip()
        if not role:
            raise ValueError("auxiliary requirement role is missing")
        normalized.append({
            "role": role,
            "necessity": str(item.get("necessity") or "possibly_required"),
            "confidence": str(item.get("confidence") or "medium"),
            "selection_status": str(
                item.get("selection_status") or "pending_user_selection"
            ),
            "carrier_ids": _values(item.get("carrier_ids")),
            "evidence": _values(item.get("evidence")),
        })
    return normalized


def _domain_ids(value: Any) -> set[str]:
    """Return canonical domain IDs while preserving unknown namespaces.

    UniProt cross-references are stored as ``Pfam:PFxxxxx`` and
    ``InterPro:IPRxxxxx`` in candidate artifacts, while older callers and
    fixtures may provide the bare identifiers.  Normalize only the two known
    domain namespaces so arbitrary namespaced text cannot become domain
    evidence accidentally.
    """

    normalized: set[str] = set()
    for raw_value in _values(value):
        token = raw_value.strip()
        namespace, separator, identifier = token.partition(":")
        if separator and namespace.strip().casefold() in DOMAIN_NAMESPACES:
            token = identifier.strip()
        if token:
            normalized.add(token.upper())
    return normalized


def _candidate_text(candidate: Mapping[str, Any]) -> str:
    return " ".join(
        str(candidate.get(field) or "")
        for field in (
            "protein_name",
            "gene_names",
            "catalytic_activities",
            "function_comments",
            "feature_annotations",
            "subcellular_locations",
            "keywords",
        )
    ).lower()


def _protein_architecture(candidate: Mapping[str, Any]) -> dict[str, bool]:
    domains = _domain_ids(candidate.get("domain_ids"))
    lineage = {item.lower() for item in _values(candidate.get("taxonomic_lineage"))}
    text = _candidate_text(candidate)
    p450 = bool(domains & P450_DOMAIN_IDS) or "cytochrome p450" in text
    reductase_domain_count = len(domains & P450_REDUCTASE_DOMAIN_IDS)
    self_sufficient = p450 and (
        reductase_domain_count >= 2
        or "bifunctional cytochrome p450/nadph" in text
        or "self-sufficient p450" in text
    )
    membrane = "transmembrane" in text or "membrane" in text
    plant = bool(
        lineage
        & {
            "viridiplantae",
            "streptophyta",
            "embryophyta",
            "tracheophyta",
        }
    )
    eukaryote = "eukaryota" in lineage
    return {
        "p450": p450,
        "self_sufficient": self_sufficient,
        "membrane": membrane,
        "plant": plant,
        "eukaryote": eukaryote,
    }


def _merge_requirements(
    requirements: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, tuple[str, ...]], dict[str, Any]] = {}
    confidence_rank = {"low": 0, "medium": 1, "high": 2}
    necessity_rank = {"possibly_required": 0, "required": 1}
    status_rank = {
        "not_required": 0,
        "integrated_in_main_enzyme": 1,
        "pending_user_selection": 2,
    }
    for raw in requirements:
        item = dict(raw)
        role = str(item.get("role") or "").strip()
        if not role:
            continue
        accessions = tuple(sorted(set(
            value.upper()
            for value in _values(item.get("main_enzyme_accessions"))
        )))
        key = (role, accessions)
        current = grouped.get(key)
        if current is None:
            grouped[key] = {
                "role": role,
                "necessity": str(item.get("necessity") or "possibly_required"),
                "confidence": str(item.get("confidence") or "medium"),
                "selection_status": str(
                    item.get("selection_status") or "pending_user_selection"
                ),
                "carrier_ids": _values(item.get("carrier_ids")),
                "evidence": _values(item.get("evidence")),
                "step_indexes": sorted({
                    int(value)
                    for value in item.get("step_indexes", [])
                    if int(value) > 0
                }),
                "main_enzyme_accessions": list(accessions),
            }
            continue
        current["necessity"] = max(
            (current["necessity"], str(item.get("necessity") or "possibly_required")),
            key=lambda value: necessity_rank.get(value, -1),
        )
        current["confidence"] = max(
            (current["confidence"], str(item.get("confidence") or "medium")),
            key=lambda value: confidence_rank.get(value, -1),
        )
        current["selection_status"] = max(
            (
                current["selection_status"],
                str(item.get("selection_status") or "pending_user_selection"),
            ),
            key=lambda value: status_rank.get(value, -1),
        )
        current["carrier_ids"] = list(dict.fromkeys([
            *current["carrier_ids"],
            *_values(item.get("carrier_ids")),
        ]))
        current["evidence"] = list(dict.fromkeys([
            *current["evidence"],
            *_values(item.get("evidence")),
        ]))
        current["step_indexes"] = sorted(set([
            *current["step_indexes"],
            *[
                int(value)
                for value in item.get("step_indexes", [])
                if int(value) > 0
            ],
        ]))
        current["main_enzyme_accessions"] = sorted(set([
            *current["main_enzyme_accessions"],
            *_values(item.get("main_enzyme_accessions")),
        ]))
    result = [grouped[key] for key in sorted(grouped)]
    for item in result:
        item["carrier_ids"] = sorted(set(item["carrier_ids"]))
        item["step_indexes"] = sorted(set(item["step_indexes"]))
        item["main_enzyme_accessions"] = sorted(set(
            value.upper() for value in item["main_enzyme_accessions"]
        ))
    return result


def merge_auxiliary_requirements(
    requirements: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Public stable merge used by candidate and set aggregation."""

    return _merge_requirements(requirements)


def infer_candidate_auxiliary_requirements(
    requirement: Mapping[str, Any],
    candidate: Mapping[str, Any],
) -> tuple[str, list[dict[str, Any]]]:
    """Refine route roles using conservative UniProt architecture evidence."""

    route_requirements = _requirements(requirement.get("auxiliary_requirements"))
    architecture = _protein_architecture(candidate)
    candidate_only_inference_allowed = str(
        candidate.get("reaction_fit_status") or ""
    ).strip().lower() == "verified"
    candidate_evidence: list[str] = []
    if architecture["p450"]:
        candidate_evidence.append("protein_domain:cytochrome_p450")
    if architecture["membrane"]:
        candidate_evidence.append("protein_annotation:membrane")
    if architecture["plant"]:
        candidate_evidence.append("taxonomy:viridiplantae")

    if architecture["self_sufficient"] and (
        route_requirements or candidate_only_inference_allowed
    ):
        integrated = []
        found_p450_role = False
        for item in route_requirements:
            updated = dict(item)
            if updated["role"] in {
                "p450_reductase",
                "oxygenase_electron_partner",
            }:
                updated["role"] = "p450_reductase"
                updated["selection_status"] = "integrated_in_main_enzyme"
                updated["confidence"] = "high"
                updated["evidence"] = list(dict.fromkeys([
                    *updated["evidence"],
                    "protein_architecture:self_sufficient_p450_reductase",
                ]))
                found_p450_role = True
            integrated.append(updated)
        if not found_p450_role:
            integrated.append({
                "role": "p450_reductase",
                "necessity": "required",
                "confidence": "high",
                "selection_status": "integrated_in_main_enzyme",
                "carrier_ids": [],
                "evidence": [
                    "protein_architecture:self_sufficient_p450_reductase"
                ],
            })
        return "self_sufficient_p450", _merge_requirements(integrated)

    if architecture["p450"] and architecture["membrane"] and (
        architecture["plant"] or architecture["eukaryote"]
    ):
        if not route_requirements and not candidate_only_inference_allowed:
            return "p450_unresolved", []
        refined = []
        found_cpr = False
        for item in route_requirements:
            updated = dict(item)
            if updated["role"] == "oxygenase_electron_partner":
                updated["role"] = "p450_reductase"
                updated["necessity"] = "required"
                updated["confidence"] = "medium"
            if updated["role"] == "p450_reductase":
                found_cpr = True
                updated["evidence"] = list(dict.fromkeys([
                    *updated["evidence"],
                    *candidate_evidence,
                ]))
            refined.append(updated)
        if not found_cpr:
            refined.append({
                "role": "p450_reductase",
                "necessity": "required",
                "confidence": "medium",
                "selection_status": "pending_user_selection",
                "carrier_ids": [],
                "evidence": candidate_evidence,
            })
        return "eukaryotic_membrane_p450", _merge_requirements(refined)

    if (
        architecture["p450"]
        and not route_requirements
        and candidate_only_inference_allowed
    ):
        route_requirements.append({
            "role": "oxygenase_electron_partner",
            "necessity": "possibly_required",
            "confidence": "medium",
            "selection_status": "pending_user_selection",
            "carrier_ids": [],
            "evidence": candidate_evidence,
        })
    system_type = "p450_unresolved" if architecture["p450"] else "not_p450"
    return system_type, _merge_requirements(route_requirements)


def annotate_candidate_auxiliary_roles(
    requirement: Mapping[str, Any],
    candidate_row: dict[str, Any],
) -> dict[str, Any]:
    system_type, requirements = infer_candidate_auxiliary_requirements(
        requirement,
        candidate_row,
    )
    candidate_row["enzyme_system_type"] = system_type
    candidate_row["required_auxiliary_roles"] = ";".join(
        item["role"] for item in requirements
    )
    pending = any(
        item["selection_status"] == "pending_user_selection"
        for item in requirements
    )
    integrated = bool(requirements) and all(
        item["selection_status"] == "integrated_in_main_enzyme"
        for item in requirements
    )
    candidate_row["auxiliary_requirement_status"] = (
        "pending_user_selection"
        if pending
        else "integrated"
        if integrated
        else "not_required"
    )
    candidate_row["auxiliary_requirements_json"] = json.dumps(
        requirements,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return candidate_row


__all__ = [
    "P450_DOMAIN_IDS",
    "P450_REDUCTASE_DOMAIN_IDS",
    "annotate_candidate_auxiliary_roles",
    "infer_candidate_auxiliary_requirements",
    "merge_auxiliary_requirements",
]
