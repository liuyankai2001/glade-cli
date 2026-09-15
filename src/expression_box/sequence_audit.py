"""Read-only expression DNA checks using the user's enzyme policy."""

from importlib.metadata import version

from dnachisel import AvoidPattern, DnaOptimizationProblem, EnforceGCContent

from src.protein_to_cds.restriction_sites import enzyme_constraints, restriction_site_audit
from src.protein_to_cds.sequence_constraints import (
    DNA_ALPHABET, HOMOPOLYMER_LIMIT, LOCAL_GC_MAX, LOCAL_GC_MIN,
    LOCAL_GC_WINDOW_NT, gc_fraction, local_gc_values, max_homopolymer_length, sha256_text,
)


def audit_expression_sequence(sequence: str, enzymes=()) -> dict:
    normalized = str(sequence or "").strip().upper()
    if not normalized or set(normalized) - DNA_ALPHABET:
        return {"engine": "DNA Chisel", "engine_version": version("dnachisel"),
                "gate_status": "FAIL", "checks": {"valid_alphabet": False}, "failed_checks": ["valid_alphabet"]}
    sites = restriction_site_audit(normalized, enzymes)
    constraints = [EnforceGCContent(mini=0.30, maxi=0.70),
                   EnforceGCContent(mini=LOCAL_GC_MIN, maxi=LOCAL_GC_MAX, window=LOCAL_GC_WINDOW_NT),
                   *enzyme_constraints(enzymes),
                   *(AvoidPattern(base * HOMOPOLYMER_LIMIT) for base in "ACGT")]
    problem = DnaOptimizationProblem(normalized, constraints=constraints, objectives=[], logger=None)
    local_values = local_gc_values(normalized)
    checks = {"valid_alphabet": True, "global_gc_pass": 0.30 <= gc_fraction(normalized) <= 0.70,
              "local_gc_pass": all(LOCAL_GC_MIN <= value <= LOCAL_GC_MAX for value in local_values),
              "forbidden_motif_pass": sites["passed"],
              "homopolymer_pass": max_homopolymer_length(normalized) < HOMOPOLYMER_LIMIT,
              "dnachisel_constraints_pass": problem.all_constraints_pass()}
    return {"engine": "DNA Chisel", "engine_version": version("dnachisel"),
            "sequence_sha256": sha256_text(normalized), "length_nt": len(normalized),
            "gc_percent": round(100 * gc_fraction(normalized), 8),
            "local_gc_min_percent": round(100 * min(local_values), 8),
            "local_gc_max_percent": round(100 * max(local_values), 8),
            "restriction_site_audit": sites,
            "forbidden_site_hits": {name: count for name, count in sites["counts"].items() if count},
            "max_homopolymer": max_homopolymer_length(normalized), "checks": checks,
            "failed_checks": [name for name, passed in checks.items() if not passed],
            "gate_status": "PASS" if all(checks.values()) else "FAIL"}
