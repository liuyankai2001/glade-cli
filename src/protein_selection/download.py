"""Download one UniProt protein by accession or chassis-ranked name search."""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import requests

from src.config.config import INPUTS_DIR
from src.main_protein_selection.settings import UNIPROT_HTTP_CONFIG, UNIPROT_SEARCH_URL
from src.main_protein_selection.taxonomy_compatibility import score_taxonomic_fit
from src.main_protein_selection.uniprot_protein_candidates import UNIPROT_FIELDS, _score_evidence, extract_gene_names

from .host_context import resolve_download_host
from .sequence_store import accession_value, entry_sequence, publish_entry, read_cached_pair


def _json_get(session: requests.Session, url: str, *, params=None) -> tuple[dict[str, Any], dict]:
    try:
        response = session.get(url, params=params, timeout=UNIPROT_HTTP_CONFIG.timeout_seconds)
        response.raise_for_status()
        data = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise ValueError(f"UniProt 请求失败：{exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("UniProt 返回的 JSON 不是对象")
    return data, getattr(response, "links", {}) or {}


def _select_by_name(config: Any, name: str, session: requests.Session):
    profile = resolve_download_host(config, session)
    literal = json.dumps(name, ensure_ascii=False)
    params = {
        "query": f"(gene_exact:{literal} OR protein_name:{literal})",
        "format": "json", "fields": UNIPROT_FIELDS + ",sequence_version", "size": 500,
    }
    next_url, visited = UNIPROT_SEARCH_URL, set()
    best_gene, best_protein, found_gene = None, None, False
    while next_url:
        if next_url in visited:
            raise ValueError("UniProt 分页链接重复，已停止查询")
        parsed_url = urlparse(next_url)
        if parsed_url.scheme != "https" or parsed_url.netloc != "rest.uniprot.org":
            raise ValueError("UniProt 分页链接来源无效")
        visited.add(next_url)
        data, links = _json_get(session, next_url, params=params)
        entries = data.get("results")
        if not isinstance(entries, list):
            raise ValueError("UniProt 名称查询缺少 results 列表")
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            is_gene = any(gene.casefold() == name.casefold() for gene in extract_gene_names(entry))
            found_gene = found_gene or is_gene
            try:
                accession = accession_value(entry.get("primaryAccession", ""))
                entry_sequence(entry)
                fit = score_taxonomic_fit(entry, profile)
                evidence, _, _ = _score_evidence(entry)
            except (ValueError, TypeError, AttributeError):
                continue
            if fit.status == "unknown":
                continue
            order = (-fit.score, -evidence, accession)
            candidate = (order, entry, fit, evidence)
            if is_gene:
                if best_gene is None or order < best_gene[0]:
                    best_gene = candidate
            elif best_protein is None or order < best_protein[0]:
                best_protein = candidate
        next_url = links.get("next", {}).get("url")
        params = None
    best = best_gene if found_gene else best_protein
    if best is None:
        raise ValueError("未找到具有有效氨基酸序列且可比较底盘分类的候选，请检查名称或改用 --uniprot")
    _, entry, fit, evidence = best
    return entry, {
        "method": "name", "query": name, "host_taxon_id": profile.taxon_id,
        "match_type": "gene_exact" if found_gene else "protein_name",
        "host_scientific_name": profile.scientific_name,
        "taxonomy_fit": asdict(fit), "evidence_score": evidence,
    }


def _result(folder: Path, metadata: dict[str, Any], reused: bool, selection: dict) -> dict[str, Any]:
    accession = metadata["accession"]
    return {
        "ok": True,
        **{key: metadata[key] for key in ("accession", "protein_name", "organism_name", "organism_id", "length_aa")},
        "fasta_path": str(folder / f"{accession}.fasta"),
        "metadata_path": str(folder / f"{accession}.json"),
        "reused_existing": reused, "selection": selection,
    }


def download_protein_sequence(
    config: Any, *, uniprot: str | None = None, name: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Export a FASTA/JSON pair without selecting enzymes or touching a manifest.

    Name searches always resolve the current chassis and all result pages first;
    accession cache hits need neither network nor a selected chassis.
    """
    if (uniprot is None) == (name is None):
        raise ValueError("必须且只能提供一个参数：uniprot 或 name")
    if not isinstance(force, bool):
        raise ValueError("force 必须是布尔值")
    folder = Path(getattr(config, "inputs_dir", INPUTS_DIR)).expanduser().resolve() / "protein_sequences"
    if uniprot is not None:
        accession = accession_value(uniprot)
        selection = {"method": "accession", "requested_accession": accession}
        if not force and (metadata := read_cached_pair(folder, accession)) is not None:
            return _result(folder, metadata, True, selection)
        with requests.Session() as session:
            entry, _ = _json_get(session, f"https://rest.uniprot.org/uniprotkb/{accession}.json")
    else:
        if not isinstance(name, str) or not name.strip():
            raise ValueError("蛋白名称不能为空")
        with requests.Session() as session:
            entry, selection = _select_by_name(config, name.strip(), session)
    metadata, reused = publish_entry(folder, entry, force=force)
    return _result(folder, metadata, reused, selection)
