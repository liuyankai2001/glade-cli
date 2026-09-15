"""Case-insensitive enzyme selection and independent, double-stranded site audits."""

from __future__ import annotations

import re
from collections.abc import Iterable
from functools import lru_cache

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Seq import Seq
from dnachisel import AvoidPattern, EnzymeSitePattern


@lru_cache(maxsize=1)
def _enzyme_index() -> dict[str, str]:
    return {name.lower(): name for name in rest_dict}


def normalize_enzymes(values: Iterable[str] | None) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        values = [values]
    if not isinstance(values, (list, tuple, set)):
        raise ValueError("限制酶配置必须是酶名列表")
    names = set()
    for value in values:
        if not isinstance(value, str):
            raise ValueError("限制酶名称必须是字符串")
        name = _enzyme_index().get(value.strip().lower())
        if name is None:
            raise ValueError(f"未知限制酶：{value}；请检查酶名")
        site = rest_dict[name]["site"].upper()
        if not site or any(base not in ambiguous_dna_values for base in site):
            raise ValueError(f"限制酶 {name} 没有可用的 DNA 识别序列")
        names.add(name)
    return sorted(names, key=str.lower)


def enzyme_constraints(enzymes: Iterable[str]) -> list[AvoidPattern]:
    return [AvoidPattern(EnzymeSitePattern(name), strand=0) for name in normalize_enzymes(enzymes)]


def restriction_site_audit(sequence: str, enzymes: Iterable[str]) -> dict:
    """Scan IUPAC motifs with lookahead, independently of DNA Chisel's solver."""
    names = normalize_enzymes(enzymes)
    sites = []
    counts = {}
    positions = set()
    for name in names:
        motif = rest_dict[name]["site"].upper()
        found = {}
        for strand, pattern in ((1, motif), (-1, str(Seq(motif).reverse_complement()))):
            expression = "".join("[" + ambiguous_dna_values[base] + "]" for base in pattern)
            for match in re.finditer("(?=" + expression + ")", sequence.upper()):
                bounds = (match.start() + 1, match.start() + len(motif))
                found.setdefault(bounds, strand)
        counts[name] = len(found)
        positions.update(found)
        sites.extend(
            {"enzyme": name, "start_1based": start, "end_1based": end, "strand": strand}
            for (start, end), strand in sorted(found.items())
        )
    return {"enzymes": names, "counts": counts, "site_count": len(positions),
            "sites": sites, "passed": not sites}


def with_restriction_audit(metrics: dict, sequence: str, enzymes: Iterable[str]) -> dict:
    return {**metrics, "restriction_site_audit": restriction_site_audit(sequence, enzymes)}
