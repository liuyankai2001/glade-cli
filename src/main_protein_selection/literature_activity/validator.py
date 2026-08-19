"""Deterministic validation and evidence grading after model extraction."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Sequence
from typing import Any

from src.main_protein_selection.literature_activity.models import (
    ExtractedActivityClaim,
    LiteratureActivityEvidence,
    LiteratureActivityRequirement,
    RetrievedLiteraturePaper,
)


VALIDATOR_VERSION = "literature_activity_validator.v1"
_DIRECT_ASSAYS = {"purified_enzyme", "biochemical_reconstitution"}
_CELL_ASSAYS = {
    "whole_cell_overexpression",
    "engineered_whole_cell",
    "genetic_knockout",
    "genetic_complementation",
    "cell_free_extract",
}


def _unique(values: Sequence[str]) -> list[str]:
    answer: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            answer.append(text)
    return answer


def _normal_text(value: str) -> str:
    return " ".join(re.sub(r"[^\w]+", " ", str(value).casefold()).split())


def _excerpt_is_grounded(excerpt: str, source_text: str) -> bool:
    needle = " ".join(str(excerpt or "").split()).casefold()
    haystack = " ".join(str(source_text or "").split()).casefold()
    return bool(needle and needle in haystack)


_GENERIC_PROTEIN_TERMS = {
    "enzyme",
    "protein",
    "phosphatase",
    "hydrolase",
    "reductase",
    "dehydrogenase",
    "synthase",
    "synthetase",
    "transferase",
    "endogenous",
    "putative",
}
_TRANSFORMATION_VERB = (
    r"(?:convert\w*|transform\w*|cataly\w*|hydroly\w*|"
    r"dephosphorylat\w*|oxid\w*|reduc\w*|produc\w*|form\w*|"
    r"yield\w*|generat\w*|synthesi\w*)"
)
_STEREO_TRANSLATION = str.maketrans({
    "−": "-",  # U+2212 minus sign
    "–": "-",  # en dash
    "—": "-",  # em dash
    "＋": "+",  # full-width plus
})


def _normalize_stereo_symbols(value: str) -> str:
    return str(value or "").translate(_STEREO_TRANSLATION)


def _phrase_pattern(value: str) -> str:
    tokens = re.findall(r"\w+", str(value).casefold(), flags=re.UNICODE)
    return r"\W+".join(re.escape(token) for token in tokens)


def _compound_text_status(name: str, compound_id: str, source_text: str) -> str:
    """Return exact/generic/conflict/missing for immutable publication text."""

    name = _normalize_stereo_symbols(name)
    source_text = _normalize_stereo_symbols(source_text)
    if not name:
        return "exact" if re.search(
            rf"(?<!\w){re.escape(compound_id)}(?!\w)",
            source_text,
            re.IGNORECASE,
        ) else "missing"
    stereo = re.match(
        r"^\s*\(([+\-]|R|S|D|L)\)\s*-\s*(.+)$",
        name,
        re.IGNORECASE,
    )
    if stereo is not None:
        marker, base = stereo.groups()
        base_pattern = _phrase_pattern(base)
        if not base_pattern:
            return "missing"
        observed_markers = {
            match.group(1).upper()
            for match in re.finditer(
                rf"\(\s*([+\-]|R|S|D|L)\s*\)\s*-\s*{base_pattern}",
                source_text,
                re.IGNORECASE,
            )
        }
        expected_marker = marker.upper()
        if any(value != expected_marker for value in observed_markers):
            return "conflict"
        if expected_marker in observed_markers:
            return "exact"
        return "generic" if re.search(
            rf"(?<!\w){base_pattern}(?!\w)",
            source_text,
            re.IGNORECASE,
        ) else "missing"
    normalized = _normal_text(name)
    return (
        "exact"
        if normalized and normalized in _normal_text(source_text)
        else "missing"
    )


def _identifier_text_supported(
    claim: ExtractedActivityClaim,
    excerpt: str,
) -> bool:
    primary_identifiers = _unique([
        claim.protein_identifier,
        claim.gene_name,
        *claim.external_identifiers,
    ])
    candidates = (
        primary_identifiers
        if primary_identifiers
        else _unique([claim.protein_name])
    )
    normalized_excerpt = _normal_text(excerpt)
    for value in candidates:
        normalized = _normal_text(value)
        if len(normalized) < 3 or normalized in _GENERIC_PROTEIN_TERMS:
            continue
        if re.search(
            rf"(?<!\w){re.escape(normalized)}(?!\w)",
            normalized_excerpt,
        ):
            return True
    return False


def _compound_terms(compounds: Sequence[Any]) -> list[str]:
    terms: list[str] = []
    for item in compounds:
        name = _normalize_stereo_symbols(str(item.name or ""))
        stereo = re.match(
            r"^\s*\([+\-RSDL]\)\s*-\s*(.+)$",
            name,
            re.IGNORECASE,
        )
        normalized = _normal_text(stereo.group(1) if stereo else name)
        terms.append(normalized or item.compound_id.casefold())
    return _unique(terms)


def _explicit_conversion_support(
    substrates: Sequence[Any],
    products: Sequence[Any],
    excerpt: str,
) -> bool:
    """Require one bounded textual substrate-to-product transformation relation."""

    text = _normal_text(excerpt)
    for substrate in _compound_terms(substrates):
        for product in _compound_terms(products):
            substrate_pattern = re.escape(substrate)
            product_pattern = re.escape(product)
            patterns = (
                rf"{_TRANSFORMATION_VERB}.{{0,160}}{substrate_pattern}.{{0,120}}"
                rf"(?:to|into|yield\w*|form\w*|produc\w*|generat\w*).{{0,80}}{product_pattern}",
                rf"{substrate_pattern}.{{0,120}}{_TRANSFORMATION_VERB}.{{0,120}}{product_pattern}",
                rf"(?:production|formation|generation|synthesis).{{0,80}}{product_pattern}"
                rf".{{0,80}}(?:from|using).{{0,80}}{substrate_pattern}",
            )
            if any(re.search(pattern, text, re.IGNORECASE) for pattern in patterns):
                return True
    return False


def _explicit_activity_contradiction(excerpt: str) -> bool:
    text = _normal_text(excerpt)
    return bool(re.search(
        rf"(?:did not|does not|cannot|could not|failed to|unable to|no detectable|"
        rf"lacked)\s+(?:\w+\s+){{0,5}}{_TRANSFORMATION_VERB}",
        text,
        re.IGNORECASE,
    ))


def _assay_text_supported(assay_type: str, excerpt: str) -> bool:
    """Require source language consistent with the model-selected assay type."""

    text = _normal_text(excerpt)
    direct_preparation = bool(re.search(
        r"\b(?:purified|isolated|recombinant|reconstitut\w*)\b",
        text,
        re.IGNORECASE,
    ))
    direct_measurement = bool(re.search(
        r"\b(?:in vitro|enzyme assay|activity assay|specific activity|kinetic\w*|"
        r"turnover|incubat\w*|catalytic activity|reconstitut\w*)\b",
        text,
        re.IGNORECASE,
    ))
    if assay_type in _DIRECT_ASSAYS:
        return direct_preparation and direct_measurement
    patterns = {
        "whole_cell_overexpression": (
            r"\b(?:overexpress\w*|over expression|heterologous expression)\b"
        ),
        "engineered_whole_cell": (
            r"\b(?:engineered (?:strain|host|cell|escherichia coli)|whole cell|"
            r"resting cell|heterologous expression)\b"
        ),
        "genetic_knockout": (
            r"\b(?:knockout|knock out|gene deletion|deleted|disrupt\w*|null mutant)\b"
        ),
        "genetic_complementation": (
            r"\b(?:complement\w*|genetic rescue|restored|restoration)\b"
        ),
        "cell_free_extract": (
            r"\b(?:cell free|lysate|crude extract|cell extract)\b"
        ),
    }
    pattern = patterns.get(assay_type)
    return bool(pattern and re.search(pattern, text, re.IGNORECASE))


def _stable_evidence_id(payload: dict) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:20]
    return f"LIT-{payload['reaction_id']}-{digest}"


def validate_activity_claim(
    requirement: LiteratureActivityRequirement,
    paper: RetrievedLiteraturePaper,
    claim: ExtractedActivityClaim,
) -> LiteratureActivityEvidence:
    """Validate source grounding, reaction identity, direction and assay strength."""

    expected_substrates = {item.compound_id: item for item in requirement.substrates}
    expected_products = {item.compound_id: item for item in requirement.products}
    matched_substrates = _unique([
        value.upper() for value in claim.matched_substrate_ids
        if value.upper() in expected_substrates
    ])
    matched_products = _unique([
        value.upper() for value in claim.matched_product_ids
        if value.upper() in expected_products
    ])
    core_substrates = [item for item in requirement.substrates if item.core]
    core_products = [item for item in requirement.products if item.core]
    excerpt_grounded = _excerpt_is_grounded(
        claim.evidence_excerpt,
        paper.source_text,
    )
    substrate_text_statuses = {
        item.compound_id: _compound_text_status(
            item.name,
            item.compound_id,
            claim.evidence_excerpt,
        )
        for item in core_substrates
    }
    product_text_statuses = {
        item.compound_id: _compound_text_status(
            item.name,
            item.compound_id,
            claim.evidence_excerpt,
        )
        for item in core_products
    }
    substrate_id_alignment = all(
        compound_id in matched_substrates
        for compound_id in requirement.core_substrate_ids
    )
    product_id_alignment = all(
        compound_id in matched_products
        for compound_id in requirement.core_product_ids
    )
    substrate_text_grounded = all(
        status in {"exact", "generic"}
        for status in substrate_text_statuses.values()
    )
    product_text_grounded = all(
        status in {"exact", "generic"}
        for status in product_text_statuses.values()
    )
    stereochemistry_exact = all(
        status == "exact"
        for status in [
            *substrate_text_statuses.values(),
            *product_text_statuses.values(),
        ]
    )
    stereochemistry_generic = any(
        status == "generic"
        for status in [
            *substrate_text_statuses.values(),
            *product_text_statuses.values(),
        ]
    )
    stereochemistry_conflict = any(
        status == "conflict"
        for status in [
            *substrate_text_statuses.values(),
            *product_text_statuses.values(),
        ]
    )
    substrate_match = substrate_id_alignment and substrate_text_grounded
    product_match = product_id_alignment and product_text_grounded
    protein_identifier_grounded = (
        excerpt_grounded
        and _identifier_text_supported(claim, claim.evidence_excerpt)
    )
    assay_text_supported = (
        excerpt_grounded
        and _assay_text_supported(claim.assay_type, claim.evidence_excerpt)
    )
    forward_signal = (
        excerpt_grounded
        and substrate_text_grounded
        and product_text_grounded
        and _explicit_conversion_support(
            core_substrates,
            core_products,
            claim.evidence_excerpt,
        )
    )
    reverse_signal = (
        excerpt_grounded
        and substrate_text_grounded
        and product_text_grounded
        and _explicit_conversion_support(
            core_products,
            core_substrates,
            claim.evidence_excerpt,
        )
    )
    direction_supported = forward_signal and not reverse_signal
    direction_contradicted = reverse_signal and not forward_signal
    route_direction_known = requirement.expected_direction != "unknown"
    explicit_negative = (
        excerpt_grounded
        and substrate_text_grounded
        and product_text_grounded
        and _explicit_activity_contradiction(claim.evidence_excerpt)
    )
    reaction_supported = (
        substrate_match
        and product_match
        and protein_identifier_grounded
        and direction_supported
        and route_direction_known
    )

    rejection_reasons: list[str] = []
    limitations = list(claim.limitations)
    if paper.retraction_status == "retracted":
        rejection_reasons.append("source publication is retracted")
    if explicit_negative:
        rejection_reasons.append("publication text explicitly reports no proposed activity")
    if direction_contradicted:
        rejection_reasons.append("publication text explicitly supports the reverse reaction")
    if stereochemistry_conflict:
        rejection_reasons.append(
            "publication text explicitly names a conflicting stereoisomer"
        )
    if not substrate_match:
        limitations.append(
            "core substrate was not jointly supported by source text and aligned compound IDs"
        )
    if not product_match:
        limitations.append(
            "core product was not jointly supported by source text and aligned compound IDs"
        )
    if not direction_supported and not direction_contradicted:
        limitations.append(
            "publication excerpt did not explicitly establish substrate-to-product direction"
        )
    if not route_direction_known:
        limitations.append(
            "selected route direction is unknown, so literature activity cannot be promoted"
        )
    if not excerpt_grounded:
        limitations.append("supporting excerpt was absent or not verbatim-grounded")
    if not protein_identifier_grounded:
        limitations.append(
            "supporting excerpt did not directly name the extracted protein or gene"
        )
    if not assay_text_supported:
        limitations.append(
            "supporting excerpt did not establish the extracted experimental assay type"
        )
    if stereochemistry_generic:
        limitations.append(
            "publication excerpt did not explicitly state the route stereochemistry"
        )
    if claim.relationship == "contradicts" and not explicit_negative:
        limitations.append(
            "extracted contradiction label was not supported by explicit negative source text"
        )
    if (
        direction_supported
        and claim.direction not in {
            requirement.expected_direction,
            "bidirectional",
            "unknown",
        }
    ):
        limitations.append(
            "extracted direction field disagreed with the source-grounded route direction"
        )
    if paper.access_level == "abstract_only":
        limitations.append("evidence was extracted from an abstract only")
    if paper.retraction_status == "expression_of_concern":
        limitations.append("source publication has an expression of concern")
    if paper.retraction_status in {"not_checked", "uncertain"}:
        limitations.append("publication retraction status was not conclusively verified")

    if rejection_reasons:
        level = "Reject"
        fit_status = "rejected"
        review_status = "rejected"
    elif (
        claim.relationship == "supports"
        and reaction_supported
        and excerpt_grounded
        and protein_identifier_grounded
        and assay_text_supported
        and claim.assay_type in _DIRECT_ASSAYS
        and claim.direct_activity_measured
        and stereochemistry_exact
        and paper.access_level == "full_text_snippets"
        and paper.retraction_status == "not_retracted"
    ):
        level = "A"
        fit_status = "verified_with_risk"
        review_status = "pending"
    elif (
        claim.relationship == "supports"
        and reaction_supported
        and excerpt_grounded
        and protein_identifier_grounded
        and assay_text_supported
        and (
            claim.assay_type in _CELL_ASSAYS
            or (
                claim.assay_type in _DIRECT_ASSAYS
                and claim.direct_activity_measured
                and paper.access_level == "abstract_only"
            )
        )
        and paper.retraction_status not in {"retracted", "expression_of_concern"}
    ):
        level = "B"
        fit_status = "verified_with_risk"
        review_status = "pending"
    else:
        level = "C"
        fit_status = "audit_only"
        review_status = "pending"

    identifier = (
        claim.protein_identifier or claim.gene_name or claim.protein_name
    ).strip()
    id_payload = {
        "reaction_id": requirement.reaction_id,
        "step_index": requirement.step_index,
        "paper_id": paper.paper_id,
        "protein_identifier": identifier.casefold(),
        "tested_substrates": claim.tested_substrates,
        "tested_products": claim.tested_products,
        "matched_substrate_ids": matched_substrates,
        "matched_product_ids": matched_products,
        "assay_type": claim.assay_type,
        "direction": claim.direction,
        "relationship": claim.relationship,
        "excerpt_sha256": hashlib.sha256(
            claim.evidence_excerpt.encode("utf-8")
        ).hexdigest(),
        "source_text_sha256": paper.source_text_sha256,
        "level": level,
    }
    evidence_id = _stable_evidence_id(id_payload)
    publication_url = (
        f"https://doi.org/{paper.doi}" if paper.doi
        else f"https://pubmed.ncbi.nlm.nih.gov/{paper.pmid}/" if paper.pmid
        else ""
    )
    source_locator = claim.source_locator or (
        f"Europe PMC full-text snippets {paper.pmcid}"
        if paper.access_level == "full_text_snippets"
        else f"PubMed abstract {paper.pmid}"
        if paper.pmid else "publication metadata"
    )
    return LiteratureActivityEvidence(
        evidence_id=evidence_id,
        step_index=requirement.step_index,
        reaction_id=requirement.reaction_id,
        equation=requirement.equation,
        expected_direction=requirement.expected_direction,
        substrate_ids=[item.compound_id for item in requirement.substrates],
        substrate_names=[item.name for item in requirement.substrates if item.name],
        product_ids=[item.compound_id for item in requirement.products],
        product_names=[item.name for item in requirement.products if item.name],
        protein_identifier=identifier,
        protein_identifier_kind=claim.protein_identifier_kind,
        gene_name=claim.gene_name,
        protein_name=claim.protein_name,
        organism_name=claim.organism_name,
        taxon_id=claim.taxon_id,
        tested_substrates=claim.tested_substrates,
        tested_products=claim.tested_products,
        matched_substrate_ids=matched_substrates,
        matched_product_ids=matched_products,
        evidence_level=level,
        assay_type=claim.assay_type,
        direct_activity_measured=claim.direct_activity_measured,
        direction=claim.direction,
        relationship=claim.relationship,
        source_databases=paper.source_databases,
        title=paper.title,
        year=paper.year,
        doi=paper.doi,
        pmid=paper.pmid,
        pmcid=paper.pmcid,
        source_urls=_unique([*paper.source_urls, publication_url]),
        source_locator=source_locator,
        access_level=paper.access_level,
        retraction_status=paper.retraction_status,
        evidence_excerpt=claim.evidence_excerpt,
        evidence_summary=claim.evidence_summary,
        fit_status=fit_status,
        review_status=review_status,
        limitations=_unique(limitations),
        validation_checks={
            "core_substrate_match": substrate_match,
            "core_product_match": product_match,
            "core_substrate_id_alignment": substrate_id_alignment,
            "core_product_id_alignment": product_id_alignment,
            "core_substrate_text_grounded": substrate_text_grounded,
            "core_product_text_grounded": product_text_grounded,
            "substrate_text_statuses": substrate_text_statuses,
            "product_text_statuses": product_text_statuses,
            "stereochemistry_exact": stereochemistry_exact,
            "stereochemistry_conflict": stereochemistry_conflict,
            "direction_text_supported": direction_supported,
            "direction_text_contradicted": direction_contradicted,
            "route_direction_known": route_direction_known,
            "excerpt_grounded": excerpt_grounded,
            "protein_identifier_grounded": protein_identifier_grounded,
            "assay_text_supported": assay_text_supported,
            "retraction_acceptable": paper.retraction_status != "retracted",
        },
        rejection_reasons=_unique(rejection_reasons),
        raw_evidence_ids=paper.raw_evidence_ids,
        source_text_sha256=paper.source_text_sha256,
    )


__all__ = ["VALIDATOR_VERSION", "validate_activity_claim"]
