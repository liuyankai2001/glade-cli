"""Normalize route-step context and build bounded literature queries."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any

from src.main_protein_selection.literature_activity.models import (
    LITERATURE_ACTIVITY_ALGORITHM_VERSION,
    LiteratureActivityRequirement,
    LiteratureSearchQuery,
    ReactionCompound,
)


_COMPOUND_ID_PATTERN = re.compile(r"\bC\d{5}\b", re.IGNORECASE)
_EQUATION_SEPARATOR = re.compile(r"\s*(?:<=>|<->|=>|->|=)\s*")
_CHASSIS_SEARCH_NAMES = {
    "ecoli_mg1655": "Escherichia coli",
    "ecoli": "Escherichia coli",
}


def _unique(values: Sequence[str]) -> list[str]:
    answer: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        key = text.casefold()
        if text and key not in seen:
            seen.add(key)
            answer.append(text)
    return answer


def _values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return _unique(re.split(r"\s*[;|]\s*", value))
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        return _unique([str(item) for item in value])
    return _unique([str(value)])


def _compound_ids(value: Any) -> list[str]:
    return _unique([
        match.group(0).upper()
        for item in _values(value)
        for match in _COMPOUND_ID_PATTERN.finditer(item)
    ])


def _label_names(value: Any) -> dict[str, str]:
    names: dict[str, str] = {}
    for label in _values(value):
        match = _COMPOUND_ID_PATTERN.search(label)
        if match is None:
            continue
        compound_id = match.group(0).upper()
        suffix = label[match.end():].strip()
        suffix = re.sub(r"^\s*[:\-]\s*", "", suffix)
        # Route labels use ``Cxxxxx (name)``. Remove exactly that one outer
        # wrapper and preserve stereochemical prefixes which are themselves
        # meaningful parentheses, e.g. ``((+)-Borneol)`` -> ``(+)-Borneol``
        # and ``((R)-compound)`` -> ``(R)-compound``.
        if suffix.startswith("(") and suffix.endswith(")"):
            suffix = suffix[1:-1].strip()
        if suffix:
            names[compound_id] = suffix
    return names


def _normalize_direction(value: Any) -> str:
    text = str(value or "").strip().casefold().replace("-", "_")
    if text in {"right_to_left", "reverse", "rtl", "backward"}:
        return "right_to_left"
    if text in {"left_to_right", "forward", "ltr"}:
        return "left_to_right"
    if text in {"bidirectional", "reversible", "both"}:
        return "bidirectional"
    return "unknown"


def _equation_sides(equation: str) -> tuple[list[str], list[str]]:
    parts = _EQUATION_SEPARATOR.split(str(equation or ""), maxsplit=1)
    if len(parts) != 2:
        return [], []
    return _compound_ids(parts[0]), _compound_ids(parts[1])


def requirement_from_mapping(
    requirement: Mapping[str, Any],
) -> LiteratureActivityRequirement:
    """Project one existing selector requirement into a strict search context."""

    step_index = int(requirement.get("step_index") or 0)
    reaction_id = str(requirement.get("reaction_id") or "").strip().upper()
    equation = str(requirement.get("equation") or "").strip()
    direction = _normalize_direction(requirement.get("direction"))
    left_ids, right_ids = _equation_sides(equation)
    if direction == "right_to_left":
        substrate_ids, product_ids = right_ids, left_ids
    else:
        substrate_ids, product_ids = left_ids, right_ids

    core_substrates = _compound_ids(requirement.get("precursor_compound_ids"))
    core_products = _compound_ids(requirement.get("produced_compound_id"))
    if not core_substrates:
        core_substrates = substrate_ids[:1]
    if not core_products:
        core_products = product_ids[:1]
    substrate_ids = _unique([*substrate_ids, *core_substrates])
    product_ids = _unique([*product_ids, *core_products])

    names = _label_names(requirement.get("precursor_compound_labels"))
    produced_id = str(requirement.get("produced_compound_id") or "").strip().upper()
    produced_name = str(requirement.get("produced_compound_name") or "").strip()
    if produced_id and produced_name:
        names[produced_id] = produced_name

    return LiteratureActivityRequirement(
        step_index=step_index,
        reaction_id=reaction_id,
        reaction_name=str(requirement.get("reaction_name") or "").strip(),
        equation=equation,
        expected_direction=direction,
        ec_numbers=_values(
            requirement.get("locked_ec_numbers")
            or requirement.get("ec_numbers")
            or requirement.get("enzyme_ecs")
        ),
        substrates=[
            ReactionCompound(
                compound_id=compound_id,
                name=names.get(compound_id, ""),
                core=compound_id in core_substrates,
            )
            for compound_id in substrate_ids
        ],
        products=[
            ReactionCompound(
                compound_id=compound_id,
                name=names.get(compound_id, ""),
                core=compound_id in core_products,
            )
            for compound_id in product_ids
        ],
    )


def normalize_requirements(
    requirements: Sequence[Mapping[str, Any] | LiteratureActivityRequirement],
) -> list[LiteratureActivityRequirement]:
    """Normalize, de-duplicate and order unresolved step contexts."""

    by_step: dict[int, LiteratureActivityRequirement] = {}
    for item in requirements:
        normalized = (
            item
            if isinstance(item, LiteratureActivityRequirement)
            else requirement_from_mapping(item)
        )
        existing = by_step.get(normalized.step_index)
        if existing is not None and existing != normalized:
            raise ValueError(
                "conflicting literature requirements for step "
                f"{normalized.step_index}"
            )
        by_step[normalized.step_index] = normalized
    return [by_step[index] for index in sorted(by_step)]


def _quoted(value: str) -> str:
    return '"' + " ".join(str(value).replace('"', " ").split()) + '"'


def _preferred_name(compounds: Sequence[ReactionCompound]) -> str:
    core = next((item for item in compounds if item.core and item.name), None)
    if core is not None:
        return core.name
    named = next((item for item in compounds if item.name), None)
    if named is not None:
        return named.name
    core_id = next((item.compound_id for item in compounds if item.core), "")
    return core_id or (compounds[0].compound_id if compounds else "")


def build_literature_queries(
    requirement: LiteratureActivityRequirement,
    *,
    chassis_key: str,
    max_queries: int = 6,
) -> list[LiteratureSearchQuery]:
    """Build short complementary queries without an unconstrained LLM search."""

    if max_queries < 1:
        raise ValueError("max_queries must be positive")
    substrate = _preferred_name(requirement.substrates)
    product = _preferred_name(requirement.products)
    reaction_name = requirement.reaction_name
    ec = requirement.ec_numbers[0] if requirement.ec_numbers else ""
    chassis = _CHASSIS_SEARCH_NAMES.get(chassis_key, chassis_key.replace("_", " "))

    pairs: list[tuple[str, str, str]] = []
    if substrate and product:
        exact_pair = f"({_quoted(substrate)} AND {_quoted(product)})"
        pairs.extend([
            (
                "PubMed",
                exact_pair + " AND (enzyme OR phosphatase OR hydrolase OR synthase OR reductase OR dehydrogenase)",
                "Search the exact substrate/product pair with broad catalytic terms.",
            ),
            (
                "Europe PMC",
                exact_pair + " AND (activity OR catalysis OR bioconversion OR biosynthesis)",
                "Search full-text metadata for direct activity language.",
            ),
            (
                "PubMed",
                exact_pair + f" AND {_quoted(chassis)}",
                "Search for endogenous or heterologous activity demonstrated in the chassis.",
            ),
            (
                "PubMed",
                exact_pair + " AND (promiscuous OR noncanonical OR dephosphorylation OR conversion)",
                "Search for experimentally used non-standard or promiscuous activity.",
            ),
        ])
    if reaction_name:
        pairs.append((
            "Europe PMC",
            f"{_quoted(reaction_name)} AND (enzyme OR gene OR protein)",
            "Search the curated reaction name and protein identity terms.",
        ))
    reaction_terms = [_quoted(requirement.reaction_id)]
    if ec:
        reaction_terms.append(_quoted(ec))
    pairs.append((
        "PubMed",
        " OR ".join(reaction_terms),
        "Search database identifiers as a low-recall exact-reference rung.",
    ))

    unique: list[LiteratureSearchQuery] = []
    seen: set[tuple[str, str]] = set()
    for source, query, rationale in pairs:
        normalized_query = " ".join(query.split())
        key = (source, normalized_query.casefold())
        if key in seen:
            continue
        seen.add(key)
        digest = hashlib.sha256(
            f"{requirement.step_index}|{source}|{normalized_query}".encode("utf-8")
        ).hexdigest()[:12]
        unique.append(LiteratureSearchQuery(
            query_id=f"lit_s{requirement.step_index}_{digest}",
            step_index=requirement.step_index,
            reaction_id=requirement.reaction_id,
            source=source,
            query=normalized_query,
            rationale=rationale,
        ))
        if len(unique) >= max_queries:
            break
    return unique


def literature_component_versions(
    *,
    model_identity: str = "unconfigured",
) -> dict[str, str]:
    """Return every component version that can change cached evidence semantics."""

    # Lazy imports avoid a query-builder import cycle while keeping version
    # ownership next to the implementing component.
    from src.main_protein_selection.literature_activity.extractor import (
        EXTRACTOR_VERSION,
    )
    from src.main_protein_selection.literature_activity.identity import (
        IDENTITY_ADAPTER_VERSION,
    )
    from src.main_protein_selection.literature_activity.retriever import (
        RETRIEVER_VERSION,
    )
    from src.main_protein_selection.literature_activity.validator import (
        VALIDATOR_VERSION,
    )
    from src.main_protein_selection.taxonomy_compatibility import (
        TAXONOMY_SCORING_POLICY_VERSION,
    )

    return {
        "retriever": RETRIEVER_VERSION,
        "extractor": EXTRACTOR_VERSION,
        "validator": VALIDATOR_VERSION,
        "identity_adapter": IDENTITY_ADAPTER_VERSION,
        "taxonomy_scoring": TAXONOMY_SCORING_POLICY_VERSION,
        "model": str(model_identity or "unconfigured"),
    }


def request_fingerprint(
    requirements: Sequence[LiteratureActivityRequirement],
    *,
    chassis_key: str,
    top_n: int,
    max_results: int,
    allow_transmembrane: bool,
    model_identity: str = "unconfigured",
) -> str:
    """Bind cache entries to reaction direction, compounds and search controls."""

    payload = {
        "algorithm_version": LITERATURE_ACTIVITY_ALGORITHM_VERSION,
        "component_versions": literature_component_versions(
            model_identity=model_identity
        ),
        "chassis_key": chassis_key,
        "top_n": int(top_n),
        "max_results": int(max_results),
        "allow_transmembrane": bool(allow_transmembrane),
        "requirements": [
            item.model_dump(mode="json") for item in requirements
        ],
    }
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


__all__ = [
    "build_literature_queries",
    "literature_component_versions",
    "normalize_requirements",
    "request_fingerprint",
    "requirement_from_mapping",
]
