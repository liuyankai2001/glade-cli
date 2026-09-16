"""Inclusive homopolymer limits and independent maximal-run audits."""

import re

from dnachisel import AvoidPattern


def validate_homopolymer_max(value) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("同聚物最大长度必须是正整数，例如 --homopolymer-max 6")
    return value


def saved_homopolymer_max(selection) -> int | None:
    value = selection.get("homopolymer_max")
    return None if value is None else validate_homopolymer_max(value)


def homopolymer_audit(sequence: str, maximum: int) -> dict:
    maximum = validate_homopolymer_max(maximum)
    runs = list(re.finditer(r"A+|C+|G+|T+", sequence.upper()))
    violations = [{"base": match.group()[0], "length_nt": len(match.group()),
                   "start_1based": match.start() + 1, "end_1based": match.end()}
                  for match in runs if len(match.group()) > maximum]
    return {"maximum_allowed": maximum, "longest_run_length": max((len(run.group()) for run in runs), default=0),
            "violation_count": len(violations), "violations": violations, "passed": not violations}


def homopolymer_constraints(maximum: int | None, sequence_length: int | None = None) -> list:
    if maximum is None:
        return []
    maximum = validate_homopolymer_max(maximum)
    if sequence_length is not None and maximum >= sequence_length:
        return []
    return [AvoidPattern(base * (maximum + 1)) for base in "ACGT"]


def with_homopolymer_audit(metrics: dict, sequence: str, maximum: int | None) -> dict:
    if maximum is None:
        return metrics
    audit = homopolymer_audit(sequence, maximum)
    metrics = {**metrics, "homopolymer_audit": audit, "checks": {**metrics["checks"], "homopolymer_pass": audit["passed"]}}
    metrics["failed_checks"] = [name for name, passed in metrics["checks"].items() if not passed]
    metrics["gate_status"] = "FAIL" if metrics["failed_checks"] else "PASS"
    return metrics


def homopolymer_summary(before: dict, after: dict) -> dict:
    return {"maximum_allowed": after["maximum_allowed"], "input_max": before["longest_run_length"],
            "final_max": after["longest_run_length"], "input_count": before["violation_count"],
            "final_count": after["violation_count"]}
