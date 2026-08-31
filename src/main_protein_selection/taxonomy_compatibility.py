from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import requests


TAXONOMY_PROFILE_SCHEMA = "taxonomy_profile.v1"
TAXONOMY_SCORING_POLICY_VERSION = "taxonomy_lca_rank.v1"
UNIPROT_TAXONOMY_URL = "https://rest.uniprot.org/taxonomy"

SCORING_WEIGHTS = {
    "function": 0.40,
    "evidence": 0.15,
    "expression": 0.20,
    "host": 0.25,
}

UNKNOWN_TAXONOMY_SCORE = 50.0
RANK_SCORES = {
    "strain": 97.0,
    "subspecies": 97.0,
    "forma": 97.0,
    "varietas": 97.0,
    "species": 92.0,
    "genus": 82.0,
    "family": 70.0,
    "order": 58.0,
    "class": 46.0,
    "phylum": 34.0,
    "kingdom": 22.0,
    "domain": 12.0,
    "superkingdom": 12.0,
}
ROOT_ONLY_SCORE = 8.0


CHASSIS_TAXON_PRESETS: dict[str, dict[str, Any]] = {
    "ecoli": {
        "name": "Escherichia coli",
        "taxon_id": 562,
        "species_taxon_id": 562,
        "strain_taxon_id": None,
    },
    "ecoli_mg1655": {
        "name": "Escherichia coli K-12 MG1655",
        "taxon_id": 511145,
        "species_taxon_id": 562,
        "strain_taxon_id": 511145,
    },
}


_ECOLI_LINEAGE = (
    {"scientificName": "cellular organisms", "taxonId": 131567, "rank": "no rank"},
    {"scientificName": "Bacteria", "taxonId": 2, "rank": "domain"},
    {"scientificName": "Pseudomonadati", "taxonId": 3379134, "rank": "kingdom"},
    {"scientificName": "Pseudomonadota", "taxonId": 1224, "rank": "phylum"},
    {"scientificName": "Gammaproteobacteria", "taxonId": 1236, "rank": "class"},
    {"scientificName": "Enterobacterales", "taxonId": 91347, "rank": "order"},
    {"scientificName": "Enterobacteriaceae", "taxonId": 543, "rank": "family"},
    {"scientificName": "Escherichia", "taxonId": 561, "rank": "genus"},
)

BUILTIN_TAXONOMY_SNAPSHOTS: dict[int, dict[str, Any]] = {
    562: {
        "scientificName": "Escherichia coli",
        "taxonId": 562,
        "rank": "species",
        "lineage": list(_ECOLI_LINEAGE),
    },
    511145: {
        "scientificName": "Escherichia coli str. K-12 substr. MG1655",
        "taxonId": 511145,
        "rank": "no rank",
        "lineage": [
            *_ECOLI_LINEAGE,
            {"scientificName": "Escherichia coli", "taxonId": 562, "rank": "species"},
            {
                "scientificName": "Escherichia coli (strain K12)",
                "taxonId": 83333,
                "rank": "strain",
            },
        ],
    },
}


@dataclass(frozen=True, slots=True)
class TaxonNode:
    taxon_id: int | None
    scientific_name: str
    rank: str


@dataclass(frozen=True, slots=True)
class ChassisTaxonomyProfile:
    chassis_key: str
    taxon_id: int
    scientific_name: str
    rank: str
    lineage: tuple[TaxonNode, ...]
    source: str
    status: str = "resolved"
    message: str = ""

    @property
    def full_lineage(self) -> tuple[TaxonNode, ...]:
        current = TaxonNode(self.taxon_id, self.scientific_name, self.rank)
        if self.lineage and self.lineage[-1].taxon_id == self.taxon_id:
            return self.lineage
        return (*self.lineage, current)

    def to_evidence(self) -> dict[str, Any]:
        return {
            "schema_version": TAXONOMY_PROFILE_SCHEMA,
            "scoring_policy_version": TAXONOMY_SCORING_POLICY_VERSION,
            "chassis_key": self.chassis_key,
            "host_taxon_id": self.taxon_id,
            "host_scientific_name": self.scientific_name,
            "host_rank": self.rank,
            "status": self.status,
            "source": self.source,
            "message": self.message,
            "lineage": [asdict(node) for node in self.full_lineage],
            "rank_scores": dict(RANK_SCORES),
            "root_only_score": ROOT_ONLY_SCORE,
            "unknown_score": UNKNOWN_TAXONOMY_SCORE,
            "scoring_weights": dict(SCORING_WEIGHTS),
            "taxonomy_fingerprint": self.semantic_fingerprint(),
        }

    def semantic_fingerprint(self) -> str:
        payload = {
            "scoring_policy_version": TAXONOMY_SCORING_POLICY_VERSION,
            "chassis_key": self.chassis_key,
            "host_taxon_id": self.taxon_id,
            "host_scientific_name": self.scientific_name,
            "host_rank": self.rank,
            "status": self.status,
            "lineage": [asdict(node) for node in self.full_lineage],
            "rank_scores": dict(RANK_SCORES),
            "root_only_score": ROOT_ONLY_SCORE,
            "unknown_score": UNKNOWN_TAXONOMY_SCORE,
            "scoring_weights": dict(SCORING_WEIGHTS),
        }
        canonical = json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class TaxonomyFit:
    status: str
    score: float
    candidate_taxon_id: int | None = None
    shared_taxon_id: int | None = None
    shared_name: str = ""
    shared_rank: str = ""
    evidence_source: str = ""


def _text(value: Any) -> str:
    return str(value or "").strip()


def _integer(value: Any) -> int | None:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


def _normal_name(value: Any) -> str:
    return " ".join(_text(value).casefold().split())


def _node(value: Mapping[str, Any]) -> TaxonNode | None:
    name = _text(value.get("scientificName") or value.get("scientific_name"))
    if not name:
        return None
    return TaxonNode(
        taxon_id=_integer(value.get("taxonId") or value.get("taxon_id")),
        scientific_name=name,
        rank=_text(value.get("rank")).casefold() or "no rank",
    )


def _nodes(values: Any) -> tuple[TaxonNode, ...]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        return ()
    result: list[TaxonNode] = []
    for value in values:
        if isinstance(value, Mapping) and (item := _node(value)) is not None:
            result.append(item)
    return tuple(result)


def _profile_from_record(
    chassis_key: str,
    record: Mapping[str, Any],
    *,
    source: str,
) -> ChassisTaxonomyProfile:
    taxon_id = _integer(record.get("taxonId"))
    name = _text(record.get("scientificName"))
    if taxon_id is None or not name:
        raise ValueError("taxonomy record lacks taxonId or scientificName")
    return ChassisTaxonomyProfile(
        chassis_key=chassis_key,
        taxon_id=taxon_id,
        scientific_name=name,
        rank=_text(record.get("rank")).casefold() or "no rank",
        lineage=_nodes(record.get("lineage")),
        source=source,
    )


def _read_cached_record(path: Path, expected_taxon_id: int) -> Mapping[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, Mapping):
        return None
    if payload.get("schema_version") != TAXONOMY_PROFILE_SCHEMA:
        return None
    if _integer(payload.get("taxon_id")) != expected_taxon_id:
        return None
    record = payload.get("record")
    if not isinstance(record, Mapping):
        return None
    if _integer(record.get("taxonId")) != expected_taxon_id:
        return None
    return record


def _write_cached_record(path: Path, taxon_id: int, record: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            {
                "schema_version": TAXONOMY_PROFILE_SCHEMA,
                "taxon_id": taxon_id,
                "record": dict(record),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def resolve_chassis_taxonomy(
    chassis_key: str,
    *,
    session: requests.Session | None = None,
    cache_root: str | Path | None = None,
    timeout_seconds: float = 20.0,
    allow_network: bool = True,
) -> ChassisTaxonomyProfile:
    preset = CHASSIS_TAXON_PRESETS.get(chassis_key)
    if preset is None:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    taxon_id = int(preset["taxon_id"])
    cache_path = (
        Path(cache_root).expanduser().resolve() / "taxonomy" / f"{taxon_id}.json"
        if cache_root is not None
        else None
    )
    if cache_path is not None and cache_path.is_file():
        cached = _read_cached_record(cache_path, taxon_id)
        if cached is not None:
            return _profile_from_record(chassis_key, cached, source="cache")

    error = ""
    if allow_network:
        http = session or requests.Session()
        try:
            response = http.get(
                f"{UNIPROT_TAXONOMY_URL}/{taxon_id}",
                timeout=timeout_seconds,
            )
            response.raise_for_status()
            record = response.json()
            if not isinstance(record, Mapping):
                raise ValueError("taxonomy response is not an object")
            profile = _profile_from_record(chassis_key, record, source="uniprot_live")
            if cache_path is not None:
                _write_cached_record(cache_path, taxon_id, record)
            return profile
        except (
            requests.RequestException,
            ValueError,
            json.JSONDecodeError,
            AttributeError,
        ) as exc:
            error = f"{type(exc).__name__}: {exc}"
    else:
        error = "network disabled"

    snapshot = BUILTIN_TAXONOMY_SNAPSHOTS.get(taxon_id)
    if snapshot is not None:
        profile = _profile_from_record(chassis_key, snapshot, source="builtin_snapshot")
        return ChassisTaxonomyProfile(
            **{
                **asdict(profile),
                "lineage": profile.lineage,
                "message": error,
            }
        )

    return ChassisTaxonomyProfile(
        chassis_key=chassis_key,
        taxon_id=taxon_id,
        scientific_name=_text(preset.get("name")) or chassis_key,
        rank="unknown",
        lineage=(),
        source="unavailable",
        status="unknown",
        message=error or "taxonomy record unavailable",
    )


def _candidate_nodes(entry: Mapping[str, Any]) -> tuple[TaxonNode, ...]:
    ranked = _nodes(entry.get("lineages"))
    if ranked:
        return ranked
    organism = entry.get("organism")
    if not isinstance(organism, Mapping):
        return ()
    raw = organism.get("lineage")
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
        return ()
    return tuple(
        TaxonNode(None, _text(name), "unknown")
        for name in raw
        if _text(name)
    )


def _effective_rank(host_nodes: Sequence[TaxonNode], shared_index: int) -> str:
    rank = host_nodes[shared_index].rank
    if rank in RANK_SCORES:
        return rank
    for node in reversed(host_nodes[:shared_index]):
        if node.rank in RANK_SCORES:
            return node.rank
    return "root"


def score_taxonomic_fit(
    entry: Mapping[str, Any],
    profile: ChassisTaxonomyProfile,
) -> TaxonomyFit:
    organism = entry.get("organism")
    organism = organism if isinstance(organism, Mapping) else {}
    candidate_taxon_id = _integer(organism.get("taxonId"))
    if candidate_taxon_id == profile.taxon_id:
        return TaxonomyFit(
            status="exact_host_taxon",
            score=100.0,
            candidate_taxon_id=candidate_taxon_id,
            shared_taxon_id=profile.taxon_id,
            shared_name=profile.scientific_name,
            shared_rank=profile.rank,
            evidence_source=f"{profile.source}+taxon_id",
        )
    if profile.status != "resolved" or not profile.full_lineage:
        return TaxonomyFit(
            status="unknown",
            score=UNKNOWN_TAXONOMY_SCORE,
            candidate_taxon_id=candidate_taxon_id,
            evidence_source=profile.source,
        )

    candidate_nodes = _candidate_nodes(entry)
    if not candidate_nodes:
        return TaxonomyFit(
            status="unknown",
            score=UNKNOWN_TAXONOMY_SCORE,
            candidate_taxon_id=candidate_taxon_id,
            evidence_source=f"{profile.source}+candidate_lineage_missing",
        )
    candidate_ids = {
        node.taxon_id for node in candidate_nodes if node.taxon_id is not None
    }
    candidate_names = {
        _normal_name(node.scientific_name) for node in candidate_nodes
    }
    host_nodes = profile.full_lineage
    shared_index: int | None = None
    source = "lineage_name"
    for index in range(len(host_nodes) - 1, -1, -1):
        node = host_nodes[index]
        if node.taxon_id is not None and node.taxon_id in candidate_ids:
            shared_index = index
            source = "lineage_taxon_id"
            break
    if shared_index is None:
        for index in range(len(host_nodes) - 1, -1, -1):
            if _normal_name(host_nodes[index].scientific_name) in candidate_names:
                shared_index = index
                break
    if shared_index is None:
        return TaxonomyFit(
            status="cross_domain",
            score=ROOT_ONLY_SCORE,
            candidate_taxon_id=candidate_taxon_id,
            evidence_source=f"{profile.source}+no_shared_lineage",
        )

    shared = host_nodes[shared_index]
    rank = _effective_rank(host_nodes, shared_index)
    if rank == "root":
        status = "root_only"
        score = ROOT_ONLY_SCORE
    else:
        status = f"same_{rank}"
        score = RANK_SCORES[rank]
    return TaxonomyFit(
        status=status,
        score=score,
        candidate_taxon_id=candidate_taxon_id,
        shared_taxon_id=shared.taxon_id,
        shared_name=shared.scientific_name,
        shared_rank=rank,
        evidence_source=f"{profile.source}+{source}",
    )


def chassis_host_taxon_id(chassis_key: str) -> int:
    preset = CHASSIS_TAXON_PRESETS.get(chassis_key)
    if preset is None:
        raise ValueError(f"Unknown chassis_key: {chassis_key}")
    return int(preset["taxon_id"])


__all__ = [
    "BUILTIN_TAXONOMY_SNAPSHOTS",
    "CHASSIS_TAXON_PRESETS",
    "ChassisTaxonomyProfile",
    "RANK_SCORES",
    "ROOT_ONLY_SCORE",
    "SCORING_WEIGHTS",
    "TAXONOMY_PROFILE_SCHEMA",
    "TAXONOMY_SCORING_POLICY_VERSION",
    "TaxonNode",
    "TaxonomyFit",
    "UNKNOWN_TAXONOMY_SCORE",
    "chassis_host_taxon_id",
    "resolve_chassis_taxonomy",
    "score_taxonomic_fit",
]
