"""Read-only expression DNA checks using the user's enzyme policy."""

from importlib.metadata import version

from dnachisel import DnaOptimizationProblem

from src.protein_to_cds.restriction_sites import enzyme_constraints, restriction_site_audit
from src.protein_to_cds.homopolymers import homopolymer_audit, homopolymer_constraints
from src.protein_to_cds.sequence_constraints import (
    DNA_ALPHABET, HOMOPOLYMER_LIMIT,
    LOCAL_GC_WINDOW_NT, gc_fraction, local_gc_values, max_homopolymer_length, sha256_text,
)


def audit_expression_sequence(sequence: str, enzymes=(), homopolymer_max: int | None = None) -> dict:
    normalized = str(sequence or "").strip().upper()
    if not normalized or set(normalized) - DNA_ALPHABET:
        return {"engine": "DNA Chisel", "engine_version": version("dnachisel"),
                "gate_status": "FAIL", "checks": {"valid_alphabet": False}, "failed_checks": ["valid_alphabet"]}
    sites = restriction_site_audit(normalized, enzymes)
    homopolymers = homopolymer_audit(normalized, HOMOPOLYMER_LIMIT - 1 if homopolymer_max is None else homopolymer_max)
    constraints = [*enzyme_constraints(enzymes),
                   *homopolymer_constraints(homopolymers["maximum_allowed"], len(normalized))]
    problem = DnaOptimizationProblem(normalized, constraints=constraints, objectives=[], logger=None)
    local_values = local_gc_values(normalized)
    checks = {"valid_alphabet": True,
              "forbidden_motif_pass": sites["passed"],
              "homopolymer_pass": homopolymers["passed"],
              "dnachisel_constraints_pass": problem.all_constraints_pass()}
    return {"engine": "DNA Chisel", "engine_version": version("dnachisel"),
            "sequence_sha256": sha256_text(normalized), "length_nt": len(normalized),
            "gc_percent": round(100 * gc_fraction(normalized), 8),
            "local_gc_window_nt": LOCAL_GC_WINDOW_NT,
            "local_gc_min_percent": round(100 * min(local_values), 8),
            "local_gc_max_percent": round(100 * max(local_values), 8),
            "restriction_site_audit": sites,
            "homopolymer_audit": homopolymers,
            "forbidden_site_hits": {name: count for name, count in sites["counts"].items() if count},
            "max_homopolymer": max_homopolymer_length(normalized), "checks": checks,
            "failed_checks": [name for name, passed in checks.items() if not passed],
            "gate_status": "PASS" if all(checks.values()) else "FAIL"}
