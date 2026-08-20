from __future__ import annotations

import hashlib
import math
import random
import re
import threading
from collections.abc import Iterable
from typing import Any

import numpy as np
from Bio.Data import CodonTable
from Bio.Seq import Seq
from dnachisel import (
    AvoidChanges,
    AvoidPattern,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
)

from src.protein_to_cds.config import CDS_CONSTRAINT_CONFIG

# 兼容原有模块字段；阈值的唯一默认值位于 config.py。
POLICY_VERSION = CDS_CONSTRAINT_CONFIG.policy_version
DNA_ALPHABET = frozenset("ACGT")
VALID_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

DEFAULT_FORBIDDEN_MOTIFS: dict[str, str] = {
    "EcoRI": "GAATTC",
    "XbaI": "TCTAGA",
    "SpeI": "ACTAGT",
    "PstI": "CTGCAG",
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "SapI": "GCTCTTC",
}

MG1655_CODON_WEIGHTS: dict[str, float] = {
    "GCA": 0.59275058,
    "GCC": 0.75588463,
    "GCG": 1.0,
    "GCT": 0.44754116,
    "TGC": 1.0,
    "TGT": 0.79396341,
    "GAC": 0.59793718,
    "GAT": 1.0,
    "GAA": 1.0,
    "GAG": 0.44822572,
    "TTC": 0.74341216,
    "TTT": 1.0,
    "GGA": 0.26143527,
    "GGC": 1.0,
    "GGG": 0.36819809,
    "GGT": 0.82961824,
    "CAC": 0.75514259,
    "CAT": 1.0,
    "ATA": 0.13511648,
    "ATC": 0.82842352,
    "ATT": 1.0,
    "AAA": 1.0,
    "AAG": 0.30317794,
    "CTA": 0.07272291,
    "CTC": 0.20946794,
    "CTG": 1.0,
    "CTT": 0.20605531,
    "TTA": 0.25849985,
    "TTG": 0.25635638,
    "ATG": 1.0,
    "AAC": 1.0,
    "AAT": 0.8079104,
    "CCA": 0.358749,
    "CCC": 0.23242983,
    "CCG": 1.0,
    "CCT": 0.29693665,
    "CAA": 0.53008113,
    "CAG": 1.0,
    "AGA": 0.08697126,
    "AGG": 0.04826464,
    "CGA": 0.15675841,
    "CGC": 1.0,
    "CGG": 0.24054366,
    "CGT": 0.95427739,
    "AGC": 1.0,
    "AGT": 0.53895921,
    "TCA": 0.43778715,
    "TCC": 0.53591186,
    "TCG": 0.55578997,
    "TCT": 0.52358181,
    "ACA": 0.29221653,
    "ACC": 1.0,
    "ACG": 0.61207559,
    "ACT": 0.37546445,
    "GTA": 0.41251389,
    "GTC": 0.58067641,
    "GTG": 1.0,
    "GTT": 0.69285124,
    "TGG": 1.0,
    "TAC": 0.76090636,
    "TAT": 1.0,
}

MG1655_PROFILE: dict[str, Any] = {
    "profile_id": "ecoli_k12_mg1655_nc_000913_3_v1",
    "host_name": "Escherichia coli K-12 MG1655",
    "organism_id": 52,
    "reference_accession": "NC_000913.3",
    "reference_cds_count": 4207,
    "global_gc_min": 0.37379921,
    "global_gc_max": 0.58017336,
    "global_gc_target": 0.52083333,
    "cai_min": 0.58377348,
    "codon_weights": MG1655_CODON_WEIGHTS,
}

GENERIC_PROFILE: dict[str, Any] = {
    "profile_id": "generic_coding_sequence_v1",
    "host_name": "generic",
    "reference_accession": "",
    "reference_cds_count": None,
    "global_gc_min": 0.30,
    "global_gc_max": 0.70,
    "global_gc_target": None,
    "cai_min": None,
    "codon_weights": None,
}

LOCAL_GC_WINDOW_NT = CDS_CONSTRAINT_CONFIG.local_gc_window_nt
LOCAL_GC_MIN = CDS_CONSTRAINT_CONFIG.local_gc_min
LOCAL_GC_MAX = CDS_CONSTRAINT_CONFIG.local_gc_max
HOMOPOLYMER_LIMIT = CDS_CONSTRAINT_CONFIG.homopolymer_limit
RARE_CODON_WEIGHT_LT = CDS_CONSTRAINT_CONFIG.rare_codon_weight_lt
RARE_CLUSTER_WINDOW_CODONS = CDS_CONSTRAINT_CONFIG.rare_cluster_window_codons
RARE_CLUSTER_MIN_CODONS = CDS_CONSTRAINT_CONFIG.rare_cluster_min_codons

_OPTIMIZER_LOCK = threading.RLock()


class CdsConstraintError(RuntimeError):
    """Raised when a CDS cannot be validated or constraint-repaired."""


def normalize_dna(sequence: str) -> str:
    return re.sub(r"\s+", "", str(sequence or "")).upper().replace("U", "T")


def normalize_protein(sequence: str) -> str:
    return re.sub(r"\s+", "", str(sequence or "")).upper().rstrip("*")


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return normalize_dna(sequence).translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _count_overlapping(sequence: str, motif: str) -> int:
    return len(re.findall(f"(?={re.escape(motif)})", sequence))


def motif_hits(sequence: str, motifs: dict[str, str]) -> dict[str, int]:
    dna = normalize_dna(sequence)
    hits: dict[str, int] = {}
    for name, motif in motifs.items():
        count = _count_overlapping(dna, motif)
        reverse = reverse_complement(motif)
        if reverse != motif:
            count += _count_overlapping(dna, reverse)
        if count:
            hits[name] = count
    return hits


def gc_fraction(sequence: str) -> float:
    dna = normalize_dna(sequence)
    if not dna:
        return float("nan")
    return (dna.count("G") + dna.count("C")) / len(dna)


def local_gc_values(sequence: str, window_nt: int = LOCAL_GC_WINDOW_NT) -> list[float]:
    dna = normalize_dna(sequence)
    if not dna:
        return []
    if len(dna) <= window_nt:
        return [gc_fraction(dna)]
    return [
        gc_fraction(dna[index : index + window_nt])
        for index in range(len(dna) - window_nt + 1)
    ]


def max_homopolymer_length(sequence: str) -> int:
    dna = normalize_dna(sequence)
    if not dna:
        return 0
    return max(len(match.group(0)) for match in re.finditer(r"(A+|C+|G+|T+)", dna))


def calculate_cai(sequence: str, weights: dict[str, float] | None) -> float | None:
    if not weights:
        return None
    dna = normalize_dna(sequence)
    table = CodonTable.unambiguous_dna_by_id[11]
    log_weights: list[float] = []
    for index in range(0, len(dna) - 2, 3):
        codon = dna[index : index + 3]
        amino_acid = table.forward_table.get(codon)
        if amino_acid is None or amino_acid in {"M", "W"}:
            continue
        log_weights.append(math.log(max(float(weights.get(codon, 0.0)), 1e-12)))
    return math.exp(sum(log_weights) / len(log_weights)) if log_weights else 1.0


def rare_codon_metrics(
    sequence: str, weights: dict[str, float] | None
) -> tuple[int | None, float | None, int | None]:
    if not weights:
        return None, None, None
    dna = normalize_dna(sequence)
    codons = [dna[index : index + 3] for index in range(0, len(dna) - 2, 3)]
    rare_flags = [
        codon not in VALID_STOP_CODONS
        and weights.get(codon, 0.0) < RARE_CODON_WEIGHT_LT
        for codon in codons
    ]
    rare_count = sum(rare_flags)
    coding_count = sum(codon not in VALID_STOP_CODONS for codon in codons)
    clusters = sum(
        sum(rare_flags[index : index + RARE_CLUSTER_WINDOW_CODONS])
        >= RARE_CLUSTER_MIN_CODONS
        for index in range(max(1, len(rare_flags) - RARE_CLUSTER_WINDOW_CODONS + 1))
    )
    return rare_count, rare_count / coding_count if coding_count else 0.0, clusters


def _profile_for_organism(organism_id: int) -> dict[str, Any]:
    if int(organism_id) == int(MG1655_PROFILE["organism_id"]):
        return dict(MG1655_PROFILE)
    profile = dict(GENERIC_PROFILE)
    profile["organism_id"] = int(organism_id)
    return profile


def _validate_additional_motifs(values: Iterable[str] | None) -> dict[str, str]:
    motifs: dict[str, str] = {}
    for index, raw in enumerate(values or (), start=1):
        motif = normalize_dna(raw)
        if len(motif) < 4 or not set(motif).issubset(DNA_ALPHABET):
            raise ValueError(
                f"additional_forbidden_motifs[{index}] must be an A/C/G/T sequence of at least 4 nt"
            )
        motifs[f"additional_{index}_{motif}"] = motif
    return motifs


def _translate_cds(sequence: str) -> tuple[str, int]:
    dna = normalize_dna(sequence)
    translated = str(Seq(dna).translate(table=11, to_stop=False))
    internal_stops = (
        translated[:-1].count("*")
        if translated.endswith("*")
        else translated.count("*")
    )
    return translated.rstrip("*"), internal_stops


def audit_cds(
    sequence: str,
    protein_sequence: str,
    profile: dict[str, Any],
    motifs: dict[str, str],
) -> dict[str, Any]:
    dna = normalize_dna(sequence)
    protein = normalize_protein(protein_sequence)
    translated, internal_stops = (
        _translate_cds(dna) if dna and len(dna) % 3 == 0 else ("", 0)
    )
    local_values = local_gc_values(dna)
    extreme_local_windows = sum(
        value < LOCAL_GC_MIN or value > LOCAL_GC_MAX for value in local_values
    )
    hits = motif_hits(dna, motifs)
    cai = calculate_cai(dna, profile.get("codon_weights"))
    rare_count, rare_fraction, rare_clusters = rare_codon_metrics(
        dna, profile.get("codon_weights")
    )
    gc_value = gc_fraction(dna)
    checks = {
        "valid_alphabet": bool(dna) and set(dna).issubset(DNA_ALPHABET),
        "length_multiple_of_three": bool(dna) and len(dna) % 3 == 0,
        "start_codon_valid": dna.startswith("ATG"),
        "stop_codon_valid": dna[-3:] in VALID_STOP_CODONS if len(dna) >= 3 else False,
        "no_internal_stop": internal_stops == 0,
        "amino_acid_identity_exact": translated == protein,
        "global_gc_pass": profile["global_gc_min"]
        <= gc_value
        <= profile["global_gc_max"],
        "local_gc_pass": extreme_local_windows == 0,
        "forbidden_motif_pass": not hits,
        "homopolymer_pass": max_homopolymer_length(dna) < HOMOPOLYMER_LIMIT,
        "cai_pass": cai is None
        or profile.get("cai_min") is None
        or cai >= profile["cai_min"],
        "rare_cluster_pass": rare_clusters is None or rare_clusters == 0,
    }
    return {
        "sequence_sha256": sha256_text(dna),
        "length_nt": len(dna),
        "gc_percent": round(100.0 * gc_value, 8),
        "local_gc_min_percent": round(100.0 * min(local_values), 8)
        if local_values
        else None,
        "local_gc_max_percent": round(100.0 * max(local_values), 8)
        if local_values
        else None,
        "extreme_local_gc_window_count": extreme_local_windows,
        "cai": round(cai, 8) if cai is not None else None,
        "rare_codon_count": rare_count,
        "rare_codon_fraction": round(rare_fraction, 8)
        if rare_fraction is not None
        else None,
        "rare_cluster_count": rare_clusters,
        "max_homopolymer": max_homopolymer_length(dna),
        "forbidden_site_count": sum(hits.values()),
        "forbidden_site_hits": hits,
        "translated_protein_sha256": sha256_text(translated),
        "checks": checks,
        "gate_status": "PASS" if all(checks.values()) else "FAIL",
        "failed_checks": [name for name, passed in checks.items() if not passed],
    }


def _rare_cluster_participating_indexes(
    sequence: str,
    weights: dict[str, float],
) -> list[int]:
    """Return rare-codon indexes that participate in at least one failing window."""

    dna = normalize_dna(sequence)
    codons = [dna[index : index + 3] for index in range(0, len(dna) - 2, 3)]
    rare_flags = [
        codon not in VALID_STOP_CODONS
        and weights.get(codon, 0.0) < RARE_CODON_WEIGHT_LT
        for codon in codons
    ]
    participating: set[int] = set()
    for start in range(max(1, len(rare_flags) - RARE_CLUSTER_WINDOW_CODONS + 1)):
        window = rare_flags[start : start + RARE_CLUSTER_WINDOW_CODONS]
        if sum(window) < RARE_CLUSTER_MIN_CODONS:
            continue
        participating.update(
            start + offset for offset, is_rare in enumerate(window) if is_rare
        )
    return sorted(participating)


def _repair_host_codon_usage(
    sequence: str,
    protein_sequence: str,
    profile: dict[str, Any],
    motifs: dict[str, str],
) -> tuple[str, dict[str, Any]]:
    """Repair CAI and rare clusters with deterministic synonymous edits.

    DNA Chisel enforces the sequence, GC, motif, and homopolymer constraints, but
    CAI and rare-cluster detection are independent host-specific gates.  A GC
    repair can therefore reduce CAI or introduce a rare cluster.  This pass only
    accepts higher-weight synonymous replacements that improve the host-usage
    gate while every sequence, translation, GC, motif, and homopolymer gate stays
    true.
    """

    current = normalize_dna(sequence)
    weights = profile.get("codon_weights")
    current_audit = audit_cds(current, protein_sequence, profile, motifs)
    initial_cluster_count = current_audit.get("rare_cluster_count")
    initial_cai = current_audit.get("cai")
    metadata: dict[str, Any] = {
        "applied": False,
        "initial_cai": initial_cai,
        "final_cai": initial_cai,
        "initial_rare_cluster_count": initial_cluster_count,
        "final_rare_cluster_count": initial_cluster_count,
        "replacement_count": 0,
        "replacements": [],
    }
    if not weights or (
        current_audit["checks"].get("cai_pass")
        and current_audit["checks"].get("rare_cluster_pass")
    ):
        return current, metadata

    usage_checks = {"cai_pass", "rare_cluster_pass"}
    non_usage_failures = [
        check
        for check in current_audit.get("failed_checks", [])
        if check not in usage_checks
    ]
    if non_usage_failures:
        return current, metadata

    table = CodonTable.unambiguous_dna_by_id[11]
    synonymous_codons: dict[str, list[str]] = {}
    for codon, amino_acid in table.forward_table.items():
        synonymous_codons.setdefault(amino_acid, []).append(codon)

    def usage_rank(audit: dict[str, Any]) -> tuple[int, int, float]:
        checks = audit["checks"]
        failed_count = sum(not checks.get(name, False) for name in usage_checks)
        cluster_count = int(audit.get("rare_cluster_count") or 0)
        cai = float(audit.get("cai") or 0.0)
        return failed_count, cluster_count, -cai

    replacements: list[dict[str, Any]] = []
    codon_count = len(current) // 3
    max_rounds = max(1, codon_count * 2)
    for _ in range(max_rounds):
        if current_audit["checks"].get("cai_pass") and current_audit["checks"].get(
            "rare_cluster_pass"
        ):
            break

        codons = [current[index : index + 3] for index in range(0, len(current), 3)]
        rare_indexes = set(_rare_cluster_participating_indexes(current, weights))
        ranked_indexes: list[tuple[tuple[Any, ...], int]] = []
        for codon_index, before in enumerate(codons):
            amino_acid = table.forward_table.get(before)
            if amino_acid is None:
                continue
            before_weight = float(weights.get(before, 0.0))
            best_weight = max(
                (
                    float(weights.get(codon, 0.0))
                    for codon in synonymous_codons.get(amino_acid, [])
                    if codon != before
                ),
                default=before_weight,
            )
            if best_weight <= before_weight:
                continue
            ranked_indexes.append(
                (
                    (
                        0 if codon_index in rare_indexes else 1,
                        -math.log(max(best_weight, 1e-12) / max(before_weight, 1e-12)),
                        codon_index,
                    ),
                    codon_index,
                )
            )

        accepted: tuple[str, dict[str, Any], dict[str, Any]] | None = None
        current_rank = usage_rank(current_audit)
        for _, codon_index in sorted(ranked_indexes):
            before = codons[codon_index]
            amino_acid = table.forward_table.get(before)
            if amino_acid is None:
                continue
            before_weight = float(weights.get(before, 0.0))
            alternatives = sorted(
                (
                    codon
                    for codon in synonymous_codons.get(amino_acid, [])
                    if codon != before
                    and float(weights.get(codon, 0.0)) > before_weight
                ),
                key=lambda codon: (-weights.get(codon, 0.0), codon),
            )
            for after in alternatives:
                candidate_codons = list(codons)
                candidate_codons[codon_index] = after
                candidate = "".join(candidate_codons)
                candidate_audit = audit_cds(
                    candidate, protein_sequence, profile, motifs
                )
                if any(
                    check not in usage_checks
                    for check in candidate_audit.get("failed_checks", [])
                ):
                    continue
                if usage_rank(candidate_audit) >= current_rank:
                    continue
                replacement = {
                    "codon_index_0": codon_index,
                    "amino_acid_position_1": codon_index + 1,
                    "before": before,
                    "after": after,
                    "amino_acid": amino_acid,
                    "before_weight": before_weight,
                    "after_weight": float(weights.get(after, 0.0)),
                    "cai_after": candidate_audit.get("cai"),
                    "rare_cluster_count_after": candidate_audit.get(
                        "rare_cluster_count"
                    ),
                }
                accepted = candidate, candidate_audit, replacement
                break
            if accepted is not None:
                break

        if accepted is None:
            break
        current, current_audit, replacement = accepted
        replacements.append(replacement)

    metadata.update(
        {
            "applied": bool(replacements),
            "final_cai": current_audit.get("cai"),
            "final_rare_cluster_count": current_audit.get("rare_cluster_count"),
            "replacement_count": len(replacements),
            "replacements": replacements,
        }
    )
    return current, metadata


def repair_cds(
    initial_sequence: str,
    protein_sequence: str,
    organism_id: int,
    accession: str = "",
    additional_forbidden_motifs: Iterable[str] | None = None,
) -> dict[str, Any]:
    initial = normalize_dna(initial_sequence)
    protein = normalize_protein(protein_sequence)
    profile = _profile_for_organism(organism_id)
    motifs = dict(DEFAULT_FORBIDDEN_MOTIFS)
    motifs.update(_validate_additional_motifs(additional_forbidden_motifs))
    initial_audit = audit_cds(initial, protein, profile, motifs)
    identity_checks = {
        key: initial_audit["checks"][key]
        for key in (
            "valid_alphabet",
            "length_multiple_of_three",
            "start_codon_valid",
            "stop_codon_valid",
            "no_internal_stop",
            "amino_acid_identity_exact",
        )
    }
    if not all(identity_checks.values()):
        raise CdsConstraintError(
            "CodonTransformer output failed encoding identity checks: "
            + ", ".join(name for name, passed in identity_checks.items() if not passed)
        )

    seed_material = "|".join(
        [
            accession.strip().upper(),
            sha256_text(protein),
            str(int(organism_id)),
            profile["profile_id"],
            POLICY_VERSION,
            ";".join(f"{key}:{value}" for key, value in sorted(motifs.items())),
        ]
    )
    seed = int(sha256_text(seed_material)[:8], 16)
    constraints: list[Any] = [
        EnforceTranslation(genetic_table="Bacterial", start_codon="keep"),
        EnforceGCContent(mini=profile["global_gc_min"], maxi=profile["global_gc_max"]),
        EnforceGCContent(
            mini=LOCAL_GC_MIN, maxi=LOCAL_GC_MAX, window=LOCAL_GC_WINDOW_NT
        ),
    ]
    constrained_motifs: set[str] = set()
    for motif in motifs.values():
        constrained_motifs.add(motif)
        constrained_motifs.add(reverse_complement(motif))
    constraints.extend(AvoidPattern(motif) for motif in sorted(constrained_motifs))
    constraints.extend(AvoidPattern(base * HOMOPOLYMER_LIMIT) for base in "ACGT")

    with _OPTIMIZER_LOCK:
        numpy_state = np.random.get_state()
        python_state = random.getstate()
        try:
            np.random.seed(seed)
            random.seed(seed)
            problem = DnaOptimizationProblem(
                initial,
                constraints=constraints,
                objectives=[AvoidChanges()],
                logger=None,
            )
            problem.resolve_constraints(final_check=True)
            problem.optimize()
            final = normalize_dna(problem.sequence)
        except Exception as exc:
            raise CdsConstraintError(
                f"constraint repair failed: {type(exc).__name__}: {exc}"
            ) from exc
        finally:
            np.random.set_state(numpy_state)
            random.setstate(python_state)

    final, host_codon_usage_repair = _repair_host_codon_usage(
        final,
        protein,
        profile,
        motifs,
    )
    final_audit = audit_cds(final, protein, profile, motifs)
    if final_audit["gate_status"] != "PASS":
        raise CdsConstraintError(
            "constraint repair produced a sequence that failed the independent gate: "
            + ", ".join(final_audit["failed_checks"])
        )

    codon_count = len(initial) // 3
    codon_changes = sum(
        initial[index : index + 3] != final[index : index + 3]
        for index in range(0, len(initial), 3)
    )
    nucleotide_changes = sum(before != after for before, after in zip(initial, final))
    return {
        "schema_version": "protein_to_cds.constraint_repair.v2",
        "policy_version": POLICY_VERSION,
        "status": "PASS",
        "host": {
            "organism_id": int(organism_id),
            "profile_id": profile["profile_id"],
            "profile_kind": "host_empirical" if int(organism_id) == 52 else "generic",
            "host_name": profile["host_name"],
            "reference_accession": profile["reference_accession"],
            "reference_cds_count": profile["reference_cds_count"],
            "global_gc_min_percent": 100.0 * profile["global_gc_min"],
            "global_gc_max_percent": 100.0 * profile["global_gc_max"],
            "cai_min": profile["cai_min"],
        },
        "optimizer": {
            "name": "DNA Chisel",
            "deterministic_seed": seed,
            "objective": "minimize nucleotide changes after satisfying hard constraints",
            "host_codon_usage_repair": host_codon_usage_repair,
        },
        "constraints": {
            "local_gc_window_nt": LOCAL_GC_WINDOW_NT,
            "local_gc_min_percent": 100.0 * LOCAL_GC_MIN,
            "local_gc_max_percent": 100.0 * LOCAL_GC_MAX,
            "homopolymer_limit": HOMOPOLYMER_LIMIT,
            "forbidden_motifs": motifs,
            "rare_codon_weight_lt": RARE_CODON_WEIGHT_LT,
            "rare_cluster_window_codons": RARE_CLUSTER_WINDOW_CODONS,
            "rare_cluster_min_codons": RARE_CLUSTER_MIN_CODONS,
        },
        "protein": {
            "accession": accession,
            "length_aa": len(protein),
            "sequence_sha256": sha256_text(protein),
        },
        "initial": initial_audit,
        "final": final_audit,
        "changes": {
            "nucleotide_change_count": nucleotide_changes,
            "nucleotide_change_fraction": round(nucleotide_changes / len(initial), 8),
            "codon_change_count": codon_changes,
            "codon_change_fraction": round(codon_changes / codon_count, 8),
        },
        "final_sequence": final,
    }
