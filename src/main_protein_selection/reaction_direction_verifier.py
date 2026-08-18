"""Auditable reaction-direction evidence for main-enzyme selection.

The pathway record says which side of a KEGG reaction must be realised.  This
module maps that requirement onto the directional Rhea quartet, then evaluates
each UniProt candidate without treating an EC/KO/master-Rhea hit as proof that
the requested direction is catalysed.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import time
from pathlib import Path
from typing import Any

import requests

from src.main_protein_selection.settings import KEGG_HTTP_CONFIG, RHEA_HTTP_CONFIG


DIRECTION_SUPPORTED = "supported"
DIRECTION_CONTRADICTED = "contradicted"
DIRECTION_UNKNOWN = "unknown"
DIRECTION_EVIDENCE_SCHEMA_VERSION = "direction_evidence.v2"
DIRECTION_ANALYZER_VERSION = "official_reaction_direction.v2"
DIRECTION_PROMPT_VERSION = "disabled_deterministic"

RHEA_DIRECTIONS_URL = (
    "https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv"
)
RHEA_CHEBI_PH_URL = (
    "https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv"
)
RHEA_SPARQL_URL = "https://sparql.rhea-db.org/sparql"
KEGG_GET_URL = "https://rest.kegg.jp/get/{compound_id}"

_KEGG_ID_RE = re.compile(r"\bC\d{5}\b", re.IGNORECASE)
_RHEA_ID_RE = re.compile(r"\b(?:RHEA:)?(\d+)\b", re.IGNORECASE)
_CHEBI_ID_RE = re.compile(r"(?:CHEBI[:_])?(\d+)\b", re.IGNORECASE)
_IGNORED_CHEBI_IDS = {"15378"}  # proton; Rhea balances it at pH 7.3


def _unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _values(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return _unique([str(item).strip() for item in value if str(item).strip()])
    return _unique([
        item.strip()
        for item in str(value or "").replace("|", ";").split(";")
        if item.strip()
    ])


def normalize_rhea_ids(value: Any) -> list[str]:
    result: list[str] = []
    for item in _values(value):
        if "RHEA-COMP" in item.upper():
            continue
        match = _RHEA_ID_RE.search(item)
        if match:
            result.append(match.group(1))
    return _unique(result)


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _split_equation(equation: str, direction: str) -> tuple[list[str], list[str]]:
    parts = re.split(r"\s*(?:<=>|=>|=)\s*", str(equation or ""), maxsplit=1)
    if len(parts) != 2:
        return [], []
    left = [value.upper() for value in _KEGG_ID_RE.findall(parts[0])]
    right = [value.upper() for value in _KEGG_ID_RE.findall(parts[1])]
    if str(direction or "").strip().lower() == "right_to_left":
        return right, left
    return left, right


class DirectionEvidenceClient:
    """Cache official Rhea and KEGG evidence with content hashes."""

    def __init__(
        self,
        session: requests.Session | None = None,
        cache_root: Path | None = None,
    ) -> None:
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", "GLADE/2.3")
        if cache_root is None:
            raise ValueError("cache_root is required")
        self.cache_root = Path(cache_root)
        self.source_records: dict[str, dict[str, Any]] = {}

    def resolved_cache_root(self) -> Path:
        return self.cache_root

    def _cached_text(
        self,
        *,
        key: str,
        url: str,
        params: dict[str, Any] | None = None,
        timeout: float = 30,
        retries: int = 3,
    ) -> str:
        root = self.resolved_cache_root()
        cache_path = root / f"{key}.txt"
        if cache_path.exists():
            text = cache_path.read_text("utf-8")
            self.source_records[key] = {
                "url": url,
                "response_sha256": _sha256(text),
                "cache_hit": True,
            }
            return text
        last_error: Exception | None = None
        for attempt in range(max(1, retries)):
            try:
                response = self.session.get(url, params=params, timeout=timeout)
                response.raise_for_status()
                text = response.text
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                cache_path.write_text(text, encoding="utf-8")
                self.source_records[key] = {
                    "url": url,
                    "params": params or {},
                    "response_sha256": _sha256(text),
                    "cache_hit": False,
                }
                return text
            except requests.RequestException as exc:
                last_error = exc
                if attempt + 1 < retries:
                    time.sleep(min(2.0, 0.2 * (2**attempt)))
        self.source_records[key] = {
            "url": url,
            "params": params or {},
            "status": "source_unavailable",
            "error": str(last_error or "request failed"),
        }
        raise RuntimeError(f"official direction source unavailable: {url}: {last_error}")

    def direction_quartets(self) -> dict[str, dict[str, str]]:
        text = self._cached_text(
            key="rhea-directions",
            url=RHEA_DIRECTIONS_URL,
            timeout=RHEA_HTTP_CONFIG.timeout_seconds,
            retries=RHEA_HTTP_CONFIG.retries,
        )
        rows = csv.DictReader(io.StringIO(text), delimiter="\t")
        quartets: dict[str, dict[str, str]] = {}
        for row in rows:
            quartet = {
                "master": str(row.get("RHEA_ID_MASTER") or "").strip(),
                "left_to_right": str(row.get("RHEA_ID_LR") or "").strip(),
                "right_to_left": str(row.get("RHEA_ID_RL") or "").strip(),
                "bidirectional": str(row.get("RHEA_ID_BI") or "").strip(),
            }
            for rhea_id in quartet.values():
                if rhea_id:
                    quartets[rhea_id] = quartet
        return quartets

    def chebi_ph_mapping(self) -> dict[str, str]:
        text = self._cached_text(
            key="chebi-pH7_3-mapping",
            url=RHEA_CHEBI_PH_URL,
            timeout=RHEA_HTTP_CONFIG.timeout_seconds,
            retries=RHEA_HTTP_CONFIG.retries,
        )
        rows = csv.DictReader(io.StringIO(text), delimiter="\t")
        return {
            str(row.get("CHEBI") or "").strip():
            str(row.get("CHEBI_PH7_3") or "").strip()
            for row in rows
            if str(row.get("CHEBI") or "").strip()
        }

    def kegg_compound_chebi_ids(self, compound_id: str) -> list[str]:
        normalized = str(compound_id or "").strip().upper()
        text = self._cached_text(
            key=f"kegg-{normalized}",
            url=KEGG_GET_URL.format(compound_id=normalized),
            timeout=KEGG_HTTP_CONFIG.timeout_seconds,
            retries=KEGG_HTTP_CONFIG.retries,
        )
        line = re.search(
            r"^\s*(?:DBLINKS\s+)?ChEBI:\s*([0-9 ]+)\s*$",
            text,
            re.MULTILINE,
        )
        if not line:
            return []
        return _unique(re.findall(r"\d+", line.group(1)))

    def kegg_compound_names(self, compound_id: str) -> list[str]:
        normalized = str(compound_id or "").strip().upper()
        text = self._cached_text(
            key=f"kegg-{normalized}",
            url=KEGG_GET_URL.format(compound_id=normalized),
            timeout=KEGG_HTTP_CONFIG.timeout_seconds,
            retries=KEGG_HTTP_CONFIG.retries,
        )
        block = re.search(
            r"^NAME\s+(.+?)(?=^[A-Z][A-Z0-9_ -]*\s{2,}|\Z)",
            text,
            re.MULTILINE | re.DOTALL,
        )
        if not block:
            return []
        values: list[str] = []
        for line in block.group(1).splitlines():
            value = line.strip().rstrip(";").strip()
            if value:
                values.append(value)
        return _unique(values)

    def rhea_sides(self, master_id: str) -> dict[str, list[str]]:
        normalized = normalize_rhea_ids([master_id])[0]
        query = f"""
PREFIX rh: <http://rdf.rhea-db.org/>
SELECT ?side ?chebi WHERE {{
  rh:{normalized} rh:side ?side .
  ?side rh:contains ?participant .
  ?participant rh:compound ?compound .
  ?compound rh:chebi ?chebi .
}}
""".strip()
        text = self._cached_text(
            key=f"rhea-sides-{normalized}",
            url=RHEA_SPARQL_URL,
            params={"query": query, "format": "application/sparql-results+json"},
            timeout=RHEA_HTTP_CONFIG.timeout_seconds,
            retries=RHEA_HTTP_CONFIG.retries,
        )
        payload = json.loads(text)
        sides = {"left": [], "right": []}
        for binding in payload.get("results", {}).get("bindings", []):
            side_uri = str(binding.get("side", {}).get("value") or "")
            chebi_uri = str(binding.get("chebi", {}).get("value") or "")
            matches = _CHEBI_ID_RE.findall(chebi_uri)
            if not matches:
                continue
            side = "left" if side_uri.endswith("_L") else "right" if side_uri.endswith("_R") else ""
            if side:
                sides[side].append(matches[-1])
        return {key: _unique(values) for key, values in sides.items()}

def _mapped_kegg_side(
    compound_ids: list[str],
    client: DirectionEvidenceClient,
    ph_mapping: dict[str, str],
) -> tuple[list[str], list[str]]:
    mapped: list[str] = []
    missing: list[str] = []
    for compound_id in compound_ids:
        raw_ids = client.kegg_compound_chebi_ids(compound_id)
        if not raw_ids:
            missing.append(compound_id)
            continue
        mapped.extend(ph_mapping.get(value, value) for value in raw_ids)
    return sorted(set(mapped) - _IGNORED_CHEBI_IDS), missing


def _split_named_equation(equation: str) -> tuple[list[str], list[str]]:
    parts = re.split(r"\s*(?:<=>|=>|<=|=)\s*", str(equation or ""), maxsplit=1)
    if len(parts) != 2:
        return [], []

    def terms(side: str) -> list[str]:
        return [
            re.sub(r"^\s*\d+(?:\.\d+)?\s+", "", item).strip()
            for item in re.split(r"\s+\+\s+", side)
            if item.strip()
        ]

    return terms(parts[0]), terms(parts[1])


def _named_side_matches(
    compound_ids: list[str],
    reaction_terms: list[str],
    client: DirectionEvidenceClient,
) -> bool:
    if not compound_ids or len(compound_ids) != len(reaction_terms):
        return False
    unmatched = list(reaction_terms)
    for compound_id in compound_ids:
        names = client.kegg_compound_names(compound_id)
        matched_index: int | None = None
        for index, term in enumerate(unmatched):
            normalized_term = _normalized_text(term)
            if any(_name_span(normalized_term, name) for name in names):
                matched_index = index
                break
        if matched_index is None:
            return False
        unmatched.pop(matched_index)
    return not unmatched


def _named_equation_alignment(
    substrates: list[str],
    products: list[str],
    equations: Any,
    client: DirectionEvidenceClient,
) -> str:
    for equation in _values(equations):
        left, right = _split_named_equation(equation)
        if (
            _named_side_matches(substrates, left, client)
            and _named_side_matches(products, right, client)
        ):
            return "rhea_left_to_right"
        if (
            _named_side_matches(substrates, right, client)
            and _named_side_matches(products, left, client)
        ):
            return "rhea_right_to_left"
    return ""


def _named_compounds_match_subset(
    compound_ids: list[str],
    reaction_terms: list[str],
    client: DirectionEvidenceClient,
) -> bool:
    if not compound_ids or not reaction_terms:
        return False
    for compound_id in compound_ids:
        names = client.kegg_compound_names(compound_id)
        if not any(
            any(_name_span(_normalized_text(term), name) for name in names)
            for term in reaction_terms
        ):
            return False
    return True


def _named_backbone_alignment(
    requirement: dict[str, Any],
    equations: Any,
    client: DirectionEvidenceClient,
) -> str:
    precursors = [
        value.upper()
        for value in _values(requirement.get("precursor_compound_ids"))
        if _KEGG_ID_RE.fullmatch(value.upper())
    ]
    product = str(requirement.get("produced_compound_id") or "").upper()
    if not precursors or not _KEGG_ID_RE.fullmatch(product):
        return ""
    for equation in _values(equations):
        left, right = _split_named_equation(equation)
        if (
            _named_compounds_match_subset(precursors, left, client)
            and _named_compounds_match_subset([product], right, client)
        ):
            return "rhea_left_to_right_backbone"
        if (
            _named_compounds_match_subset(precursors, right, client)
            and _named_compounds_match_subset([product], left, client)
        ):
            return "rhea_right_to_left_backbone"
    return ""


def enrich_requirements_with_direction_context(
    requirements: list[dict[str, Any]],
    client: DirectionEvidenceClient | None = None,
) -> list[dict[str, Any]]:
    """Attach required/opposite Rhea directions; failures degrade to unknown."""

    evidence_client = client or DirectionEvidenceClient()
    try:
        quartets = evidence_client.direction_quartets()
        ph_mapping = evidence_client.chebi_ph_mapping()
    except Exception as exc:
        for requirement in requirements:
            requirement.update({
                "direction_evidence_status": "source_unavailable",
                "direction_context_error": str(exc),
                "required_rhea_direction_ids": [],
                "opposite_rhea_direction_ids": [],
                "rhea_bidirectional_ids": [],
            })
        return [{
            "step_index": int(item.get("step_index") or 0),
            "reaction_id": str(item.get("reaction_id") or ""),
            "status": "source_unavailable",
            "error": str(exc),
        } for item in requirements]

    records: list[dict[str, Any]] = []
    for requirement in requirements:
        original_ids = normalize_rhea_ids(
            requirement.get("kegg_rhea_ids") or requirement.get("rhea_ids")
        )
        master_ids = normalize_rhea_ids(requirement.get("rhea_master_ids"))
        quartet: dict[str, str] | None = None
        for rhea_id in [*master_ids, *original_ids]:
            if rhea_id in quartets:
                quartet = quartets[rhea_id]
                break
        record: dict[str, Any] = {
            "step_index": int(requirement.get("step_index") or 0),
            "reaction_id": str(requirement.get("reaction_id") or ""),
            "requested_direction": str(requirement.get("direction") or ""),
            "kegg_rhea_ids": original_ids,
            "status": "unknown",
        }
        if not quartet:
            requirement.update({
                "direction_evidence_status": "unknown_no_rhea_quartet",
                "required_rhea_direction_ids": [],
                "opposite_rhea_direction_ids": [],
                "rhea_bidirectional_ids": [],
                "rhea_retrieval_ids": original_ids,
            })
            record["reason"] = "No official Rhea quartet was mapped to the KEGG reaction"
            records.append(record)
            continue
        master = quartet["master"]
        try:
            substrates, products = _split_equation(
                str(requirement.get("equation") or ""),
                str(requirement.get("direction") or ""),
            )
            required_left, missing_left = _mapped_kegg_side(substrates, evidence_client, ph_mapping)
            required_right, missing_right = _mapped_kegg_side(products, evidence_client, ph_mapping)
            rhea_sides = evidence_client.rhea_sides(master)
            rhea_left = sorted(set(rhea_sides["left"]) - _IGNORED_CHEBI_IDS)
            rhea_right = sorted(set(rhea_sides["right"]) - _IGNORED_CHEBI_IDS)
            comparable = bool(required_left and required_right and rhea_left and rhea_right)
            forward = comparable and required_left == rhea_left and required_right == rhea_right
            reverse = comparable and required_left == rhea_right and required_right == rhea_left
            named_alignment = ""
            if not forward and not reverse:
                named_alignment = _named_equation_alignment(
                    substrates,
                    products,
                    requirement.get("rhea_master_equations"),
                    evidence_client,
                )
                if not named_alignment:
                    named_alignment = _named_backbone_alignment(
                        requirement,
                        requirement.get("rhea_master_equations"),
                        evidence_client,
                    )
                forward = named_alignment in {
                    "rhea_left_to_right",
                    "rhea_left_to_right_backbone",
                }
                reverse = named_alignment in {
                    "rhea_right_to_left",
                    "rhea_right_to_left_backbone",
                }
            if forward:
                required_id = quartet["left_to_right"]
                opposite_id = quartet["right_to_left"]
                alignment = named_alignment or "rhea_left_to_right"
            elif reverse:
                required_id = quartet["right_to_left"]
                opposite_id = quartet["left_to_right"]
                alignment = named_alignment or "rhea_right_to_left"
            else:
                required_id = ""
                opposite_id = ""
                alignment = "unresolved"
            status = (
                "resolved_backbone_names"
                if required_id and named_alignment.endswith("_backbone")
                else "resolved_carrier_aware_names"
                if required_id and named_alignment
                else "resolved"
                if required_id
                else "unknown_side_mismatch"
            )
            requirement.update({
                "direction_evidence_status": status,
                "rhea_master_ids": [master],
                "required_rhea_direction_ids": [required_id] if required_id else [],
                "opposite_rhea_direction_ids": [opposite_id] if opposite_id else [],
                "rhea_bidirectional_ids": [quartet["bidirectional"]],
                "rhea_retrieval_ids": _unique([
                    required_id, quartet["bidirectional"], master,
                ]),
            })
            record.update({
                "status": status,
                "rhea_quartet": quartet,
                "alignment": alignment,
                "required_kegg_substrates": substrates,
                "required_kegg_products": products,
                "required_chebi_left": required_left,
                "required_chebi_right": required_right,
                "rhea_chebi_left": rhea_left,
                "rhea_chebi_right": rhea_right,
                "missing_kegg_chebi_mappings": _unique(missing_left + missing_right),
                "required_rhea_direction_ids": requirement["required_rhea_direction_ids"],
                "opposite_rhea_direction_ids": requirement["opposite_rhea_direction_ids"],
                "rhea_bidirectional_ids": requirement["rhea_bidirectional_ids"],
            })
        except Exception as exc:
            requirement.update({
                "direction_evidence_status": "source_unavailable",
                "direction_context_error": str(exc),
                "required_rhea_direction_ids": [],
                "opposite_rhea_direction_ids": [],
                "rhea_bidirectional_ids": [quartet["bidirectional"]],
                "rhea_retrieval_ids": _unique([quartet["bidirectional"], master]),
            })
            record.update({"status": "source_unavailable", "error": str(exc)})
        records.append(record)
    return records


def _candidate_ecs(candidate: Any) -> list[str]:
    if isinstance(candidate, dict):
        return _values(candidate.get("ec_numbers") or candidate.get("ec_number"))
    return _values(getattr(candidate, "ec_numbers", []))


def _candidate_value(candidate: Any, name: str) -> Any:
    return candidate.get(name) if isinstance(candidate, dict) else getattr(candidate, name, None)


def _candidate_activity_records(candidate: Any) -> list[dict[str, Any]]:
    value = _candidate_value(candidate, "catalytic_activity_records_json")
    if not value:
        value = _candidate_value(candidate, "catalytic_activity_records")
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    text = str(value or "").strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return []
    if not isinstance(parsed, list):
        return []
    return [item for item in parsed if isinstance(item, dict)]


def _physiological_rhea_ids(candidate: Any) -> list[str]:
    values: list[str] = []
    for activity in _candidate_activity_records(candidate):
        physiological = activity.get("physiological_reactions")
        if not isinstance(physiological, list):
            continue
        for record in physiological:
            if not isinstance(record, dict):
                continue
            cross_references = record.get("reaction_cross_references")
            if not isinstance(cross_references, list):
                continue
            for cross_reference in cross_references:
                if not isinstance(cross_reference, dict):
                    continue
                if str(cross_reference.get("database") or "").lower() != "rhea":
                    continue
                values.append(str(cross_reference.get("id") or ""))
    return normalize_rhea_ids(values)


def _normalized_text(value: Any) -> str:
    text = str(value or "").lower()
    text = text.replace("α", "alpha").replace("β", "beta")
    text = re.sub(r"[-‐‑‒–—,;:/()\[\]{}]", " ", text)
    return " ".join(re.findall(r"[a-z0-9]+", text))


def _compound_names_from_labels(value: Any) -> list[str]:
    names: list[str] = []
    for item in _values(value):
        match = re.search(r"\((.+)\)\s*$", item)
        name = match.group(1) if match else item
        normalized = _normalized_text(name)
        if normalized and not re.fullmatch(r"c\d{5}", normalized):
            names.append(normalized)
    return _unique(names)


def _compound_name_variants(name: str) -> list[str]:
    normalized = _normalized_text(name)
    variants = [normalized]
    relaxed = re.sub(r"^(?:trans\s+){2}", "", normalized).strip()
    if relaxed and relaxed != normalized:
        variants.append(relaxed)
    return _unique(variants)


def _name_span(text: str, name: str) -> tuple[int, int] | None:
    for variant in sorted(_compound_name_variants(name), key=len, reverse=True):
        if not variant:
            continue
        match = re.search(
            rf"(?<![a-z0-9]){re.escape(variant)}(?![a-z0-9])",
            text,
        )
        if match:
            return match.span()
    return None


def _candidate_has_reviewed_literature(candidate: Any) -> bool:
    reviewed = _candidate_value(candidate, "reviewed")
    if isinstance(reviewed, bool):
        is_reviewed = reviewed
    else:
        is_reviewed = str(reviewed or "").strip().lower() in {
            "1",
            "true",
            "yes",
        }
    return is_reviewed and bool(_values(_candidate_value(candidate, "publication_ids")))


def _route_compound_names(
    requirement: dict[str, Any],
) -> tuple[list[str], str]:
    substrates = _compound_names_from_labels(
        requirement.get("precursor_compound_labels")
    )
    product = _normalized_text(requirement.get("produced_compound_name"))
    return substrates, product


def _via_chain_orientation(
    text: str,
    substrates: list[str],
    product: str,
) -> str:
    patterns = (
        re.compile(
            r"converts? (?P<source>.+?) into (?P<target>.+?) via "
            r"(?:the )?(?:intermediary|intermediates?) of (?P<via>.+?)"
            r"(?: by | through |$)"
        ),
        re.compile(
            r"converts? (?P<source>.+?) to (?P<target>.+?) via "
            r"(?:the )?(?:intermediary|intermediates?) of (?P<via>.+?)"
            r"(?: by | through |$)"
        ),
    )
    for pattern in patterns:
        match = pattern.search(text)
        if not match:
            continue
        source = match.group("source")
        target = match.group("target")
        via = match.group("via")

        def logical_position(name: str) -> int | None:
            if _name_span(source, name):
                return 0
            via_span = _name_span(via, name)
            if via_span:
                return 1 + via_span[0]
            if _name_span(target, name):
                return 1_000_000
            return None

        substrate_positions = [logical_position(name) for name in substrates]
        product_position = logical_position(product)
        if product_position is None or any(
            position is None for position in substrate_positions
        ):
            continue
        known_substrates = [
            int(position) for position in substrate_positions if position is not None
        ]
        if max(known_substrates) < product_position:
            return "forward"
        if product_position < min(known_substrates):
            return "reverse"
    return ""


def _direct_text_orientation(
    text: str,
    substrates: list[str],
    product: str,
) -> str:
    substrate_spans = [_name_span(text, name) for name in substrates]
    product_span = _name_span(text, product)
    if not product_span or any(span is None for span in substrate_spans):
        return ""
    known_substrates = [span for span in substrate_spans if span is not None]
    forward_markers = [
        match.span()
        for pattern in (r"\binto\b", r"\bto yield\b", r"\bto give\b")
        for match in re.finditer(pattern, text)
    ]
    for marker_start, marker_end in forward_markers:
        if (
            all(span[1] <= marker_start for span in known_substrates)
            and product_span[0] >= marker_end
        ):
            return "forward"
        if (
            product_span[1] <= marker_start
            and all(span[0] >= marker_end for span in known_substrates)
        ):
            return "reverse"

    for verb in ("produces", "produce", "forms", "form"):
        for marker in re.finditer(rf"\b{verb}\b", text):
            from_marker = re.search(r"\bfrom\b", text[marker.end():])
            if not from_marker:
                continue
            from_start = marker.end() + from_marker.start()
            if (
                product_span[0] >= marker.end()
                and product_span[1] <= from_start
                and all(span[0] >= from_start for span in known_substrates)
            ):
                return "forward"
    return ""


def _text_direction_decision(
    requirement: dict[str, Any],
    candidate: Any,
) -> dict[str, Any] | None:
    if not _candidate_has_reviewed_literature(candidate):
        return None
    substrates, product = _route_compound_names(requirement)
    if not substrates or not product:
        return None
    texts = _values(_candidate_value(candidate, "function_comments"))
    orientations: set[str] = set()
    multistep = False
    for raw_text in texts:
        text = _normalized_text(raw_text)
        if re.search(
            r"\b(?:does not|do not|cannot|unable to|may|might|possibly|"
            r"probably|probable)\b",
            text,
        ):
            continue
        orientation = _via_chain_orientation(text, substrates, product)
        if orientation:
            multistep = True
            orientations.add(orientation)
            continue
        orientation = _direct_text_orientation(text, substrates, product)
        if orientation:
            orientations.add(orientation)
    if not orientations or len(orientations) > 1:
        return None
    orientation = next(iter(orientations))
    verdict = (
        DIRECTION_SUPPORTED
        if orientation == "forward"
        else DIRECTION_CONTRADICTED
    )
    evidence_level = (
        "reviewed_uniprot_multistep_chain"
        if multistep
        else "reviewed_uniprot_direction_text"
    )
    source_ids = [
        f"UniProt:{str(_candidate_value(candidate, 'accession') or '').upper()}"
    ]
    source_ids.extend(_values(_candidate_value(candidate, "publication_ids")))
    return {
        "verdict": verdict,
        "confidence": "medium",
        "evidence_level": evidence_level,
        "source_ids": _unique(source_ids),
        "evidence": [
            "Reviewed UniProt function text explicitly supports the route "
            f"{orientation} direction"
        ],
        "required_rhea_direction_ids": normalize_rhea_ids(
            requirement.get("required_rhea_direction_ids")
        ),
    }


def direction_decision_for_candidate(
    requirement: dict[str, Any],
    candidate: Any,
) -> dict[str, Any]:
    """Return a deterministic tri-state direction verdict."""

    if not requirement.get("direction_evidence_status"):
        return {
            "verdict": "",
            "confidence": "",
            "evidence_level": "legacy_unassessed",
            "source_ids": [],
            "evidence": [],
            "required_rhea_direction_ids": [],
        }
    required = set(normalize_rhea_ids(requirement.get("required_rhea_direction_ids")))
    opposite = set(normalize_rhea_ids(requirement.get("opposite_rhea_direction_ids")))
    bidirectional = set(normalize_rhea_ids(requirement.get("rhea_bidirectional_ids")))
    annotated_rhea = normalize_rhea_ids(_candidate_value(candidate, "rhea_ids"))
    retrieval_rhea = normalize_rhea_ids(_candidate_value(candidate, "matched_rhea_ids"))
    # A retrieval query is not an annotation. Prefer UniProt's actual
    # catalytic-activity cross-reference whenever it is present.
    candidate_rhea = set(annotated_rhea or retrieval_rhea)
    candidate_rhea.update(_physiological_rhea_ids(candidate))
    sources = [f"RHEA:{value}" for value in sorted(candidate_rhea)]
    required_hits = candidate_rhea & (required | bidirectional)
    opposite_hits = candidate_rhea & opposite
    if required_hits and opposite_hits:
        return {
            "verdict": DIRECTION_UNKNOWN,
            "confidence": "low",
            "evidence_level": "conflicting_directional_rhea",
            "source_ids": [
                f"RHEA:{value}"
                for value in sorted(required_hits | opposite_hits)
            ],
            "evidence": [
                "Candidate contains both required and opposite directional "
                "Rhea annotations"
            ],
            "required_rhea_direction_ids": sorted(required),
        }
    if required_hits:
        matched = sorted(required_hits)
        return {
            "verdict": DIRECTION_SUPPORTED,
            "confidence": "high",
            "evidence_level": "candidate_directional_rhea",
            "source_ids": [f"RHEA:{value}" for value in matched],
            "evidence": ["Candidate is annotated to the required or bidirectional Rhea member"],
            "required_rhea_direction_ids": sorted(required),
        }
    if opposite_hits:
        matched = sorted(opposite_hits)
        return {
            "verdict": DIRECTION_CONTRADICTED,
            "confidence": "high",
            "evidence_level": "candidate_opposite_directional_rhea",
            "source_ids": [f"RHEA:{value}" for value in matched],
            "evidence": ["Candidate is annotated to the opposite member of the Rhea direction quartet"],
            "required_rhea_direction_ids": sorted(required),
        }

    text_decision = _text_direction_decision(requirement, candidate)
    if text_decision is not None:
        return text_decision

    substrates, products = _split_equation(
        str(requirement.get("equation") or ""),
        str(requirement.get("direction") or ""),
    )
    ecs = _candidate_ecs(candidate) or _values(requirement.get("ec_numbers"))
    if any(ec.startswith("3.") for ec in ecs):
        if "C00001" in products and "C00001" not in substrates:
            return {
                "verdict": DIRECTION_CONTRADICTED,
                "confidence": "high",
                "evidence_level": "ec_class_mechanism",
                "source_ids": [f"EC:{ec}" for ec in ecs if ec.startswith("3.")],
                "evidence": ["Hydrolase EC class requires water as a reactant, but the locked route produces water"],
                "required_rhea_direction_ids": sorted(required),
            }
        if "C00001" in substrates and "C00001" not in products:
            return {
                "verdict": DIRECTION_SUPPORTED,
                "confidence": "medium",
                "evidence_level": "ec_class_mechanism",
                "source_ids": [f"EC:{ec}" for ec in ecs if ec.startswith("3.")],
                "evidence": ["Hydrolase EC class and locked water-consuming direction are mechanistically aligned"],
                "required_rhea_direction_ids": sorted(required),
            }
    if any(ec.startswith("6.") for ec in ecs) and "C00002" in substrates:
        return {
            "verdict": DIRECTION_SUPPORTED,
            "confidence": "medium",
            "evidence_level": "ec_class_mechanism",
            "source_ids": [f"EC:{ec}" for ec in ecs if ec.startswith("6.")],
            "evidence": ["Ligase EC class is aligned with the locked ATP-consuming synthesis direction"],
            "required_rhea_direction_ids": sorted(required),
        }

    return {
        "verdict": DIRECTION_UNKNOWN,
        "confidence": "low",
        "evidence_level": "insufficient_direction_evidence",
        "source_ids": sources,
        "evidence": ["Official evidence does not establish that this candidate catalyses the locked direction"],
        "required_rhea_direction_ids": sorted(required),
    }


def direction_evidence_artifact(
    *,
    solution_id: int,
    requirements: list[dict[str, Any]],
    context_records: list[dict[str, Any]],
    candidates_by_step: dict[int, list[Any]],
    client: DirectionEvidenceClient,
) -> dict[str, Any]:
    steps: list[dict[str, Any]] = []
    for requirement in requirements:
        step_index = int(requirement.get("step_index") or 0)
        candidate_records = []
        for candidate in candidates_by_step.get(step_index, []):
            decision = direction_decision_for_candidate(requirement, candidate)
            candidate_records.append({
                "accession": str(_candidate_value(candidate, "accession") or ""),
                "protein_name": str(_candidate_value(candidate, "protein_name") or ""),
                "verdict": decision["verdict"] or DIRECTION_UNKNOWN,
                "confidence": decision["confidence"] or "low",
                "evidence_level": decision["evidence_level"],
                "source_ids": decision["source_ids"],
                "evidence": decision["evidence"],
                "required_rhea_direction_ids": decision["required_rhea_direction_ids"],
            })
        steps.append({
            "step_index": step_index,
            "reaction_id": str(requirement.get("reaction_id") or ""),
            "requested_direction": str(requirement.get("direction") or ""),
            "context_status": str(requirement.get("direction_evidence_status") or ""),
            "required_rhea_direction_ids": normalize_rhea_ids(
                requirement.get("required_rhea_direction_ids")
            ),
            "opposite_rhea_direction_ids": normalize_rhea_ids(
                requirement.get("opposite_rhea_direction_ids")
            ),
            "rhea_bidirectional_ids": normalize_rhea_ids(
                requirement.get("rhea_bidirectional_ids")
            ),
            "candidates": candidate_records,
        })
    return {
        "schema_version": DIRECTION_EVIDENCE_SCHEMA_VERSION,
        "analyzer_version": DIRECTION_ANALYZER_VERSION,
        "prompt_version": DIRECTION_PROMPT_VERSION,
        "selected_solution_id": int(solution_id),
        "policy": {
            "supported": "verified",
            "contradicted": "rejected",
            "unknown": "verified_with_risk",
            "clean_route_preference": True,
            "allowed_sources": [
                "KEGG compound names and mappings",
                "Rhea direction quartets and reaction sides",
                "UniProt physiological reactions",
                "reviewed UniProt function text with linked publications",
            ],
            "confidence_policy": {
                "directional_rhea": "high",
                "reviewed_explicit_function_text": "medium",
                "ec_ko_or_master_rhea_only": "low_unknown",
            },
            "free_text_policy": (
                "controlled directional templates only; ambiguous or "
                "negated text remains unknown"
            ),
        },
        "sources": client.source_records,
        "context": context_records,
        "agent": {
            "status": "disabled_deterministic",
            "assessments": [],
            "error": "",
        },
        "steps": steps,
    }


__all__ = [
    "DIRECTION_CONTRADICTED",
    "DIRECTION_EVIDENCE_SCHEMA_VERSION",
    "DIRECTION_SUPPORTED",
    "DIRECTION_UNKNOWN",
    "DirectionEvidenceClient",
    "direction_decision_for_candidate",
    "direction_evidence_artifact",
    "enrich_requirements_with_direction_context",
    "normalize_rhea_ids",
]
