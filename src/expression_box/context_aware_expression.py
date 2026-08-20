from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Literal

from pydantic import BaseModel, Field, model_validator

from src.tools.expression_cassette_assembly_tools.discrete_expression_optimizer import (
    candidate_space_size,
    cp_model,
    quantize_cost,
    solve_lexicographic_cp_sat,
)


StrengthTier = Literal["low", "medium", "high"]
EvidenceGrade = Literal["A", "B", "C", "D"]
AnalysisTier = Literal["primary", "exploratory"]
TerminatorDirection = Literal["forward", "bidirectional", "forward_or_bidirectional"]


class PromoterIntent(BaseModel):
    regulation: Literal["constitutive", "inducible", "repressible"]
    strength_tier: StrengthTier
    regulator_required: str = "not_required"
    evidence_id: str = ""
    evidence_grade: EvidenceGrade = "C"
    analysis_tier: AnalysisTier = "exploratory"


class TerminatorIntent(BaseModel):
    direction: TerminatorDirection = "forward_or_bidirectional"
    strength_tier: StrengthTier = "high"
    min_efficiency: float | None = Field(default=None, ge=0, le=1)
    dataset_id: str = ""
    evidence_id: str = ""
    evidence_grade: EvidenceGrade = "B"
    analysis_tier: AnalysisTier = "primary"


class CassetteExpressionTarget(BaseModel):
    target_id: str = Field(min_length=1)
    cassette_index: int = Field(ge=1)
    promoter: PromoterIntent
    terminator: TerminatorIntent
    avoid_reused_regulatory_parts: bool = True


class ExpressionTargetSet(BaseModel):
    target_set_id: str = ""
    cassette_targets: list[CassetteExpressionTarget] = Field(min_length=1)

    @model_validator(mode="after")
    def unique_cassettes(self) -> "ExpressionTargetSet":
        keys = [(row.target_id, row.cassette_index) for row in self.cassette_targets]
        if len(keys) != len(set(keys)):
            raise ValueError("cassette_targets contains duplicate target_id/cassette_index pairs")
        return self


def _tier_bounds(tier: StrengthTier) -> tuple[float, float]:
    return {
        "low": (0.0, 100.0 / 3.0),
        "medium": (100.0 / 3.0, 200.0 / 3.0),
        "high": (200.0 / 3.0, 100.0),
    }[tier]


def _tier_evaluation(percentile: float, tier: StrengthTier) -> dict[str, Any]:
    low, high = _tier_bounds(tier)
    if not math.isfinite(percentile) or percentile < 0:
        return {
            "status": "FAIL",
            "target_low": low,
            "target_high": high,
            "observed": percentile,
            "violation_loss": 2.0,
            "ranking_loss": 2.0,
        }
    status = "PASS" if low <= percentile <= high else "FAIL"
    boundary_loss = 0.0 if status == "PASS" else min(abs(percentile - low), abs(percentile - high)) / 100.0
    centre_loss = abs(percentile - (low + high) / 2.0) / 100.0
    return {
        "status": status,
        "target_low": low,
        "target_high": high,
        "observed": percentile,
        "violation_loss": boundary_loss,
        "ranking_loss": centre_loss,
    }


def _candidate_hard_status(part: dict[str, Any], role: str) -> tuple[int, list[str]]:
    reasons: list[str] = []
    if str(part.get("role") or "").lower() != role:
        reasons.append("role_mismatch")
    if not bool(part.get("sequence_available")):
        reasons.append("sequence_unavailable")
    if str(part.get("sequence_type") or role).lower() != role:
        reasons.append("sequence_type_mismatch")
    metadata = part.get("registry_metadata") if isinstance(part.get("registry_metadata"), dict) else {}
    if str(metadata.get("part_results") or "").lower() == "fails":
        reasons.append("registry_failure")
    return len(reasons), reasons


def evaluate_promoter_candidate(
    part: dict[str, Any],
    target: CassetteExpressionTarget,
) -> dict[str, Any]:
    hard, reasons = _candidate_hard_status(part, "promoter")
    observed_regulation = str(part.get("regulation") or "unknown").lower()
    regulation_status = "PASS" if observed_regulation == target.promoter.regulation else "FAIL"
    if regulation_status == "FAIL":
        hard += 1
        reasons.append("regulation_mismatch")
    regulator_required = str(part.get("regulator_required") or "unknown").lower()
    if target.promoter.regulator_required not in {"", "not_required"} and regulator_required in {"", "unknown"}:
        hard += 1
        reasons.append("required_regulator_unverified")
    percentile = float(part.get("activity_percentile", -1.0) or -1.0)
    strength = _tier_evaluation(percentile, target.promoter.strength_tier)
    evidence_grade = str(part.get("evidence_grade") or "D").upper()
    evidence_status = "PASS" if evidence_grade in {"A", "B"} else "REVIEW" if evidence_grade == "C" else "FAIL"
    if evidence_status == "FAIL":
        hard += 1
        reasons.append("evidence_ineligible")
    status = "PASS" if hard == 0 and strength["status"] == "PASS" else "FAIL"
    return {
        "role": "promoter",
        "part_id": str(part.get("part_id") or ""),
        "status": status,
        "hard_violations": hard,
        "intent_violations": int(strength["status"] != "PASS"),
        "violation_loss": float(strength["violation_loss"]),
        "ranking_loss": float(strength["ranking_loss"]),
        "regulation_status": regulation_status,
        "strength_status": strength["status"],
        "target_low": strength["target_low"],
        "target_high": strength["target_high"],
        "activity_percentile": percentile,
        "activity_value": float(part.get("activity_value", -1.0) or -1.0),
        "activity_unit": str(part.get("activity_unit") or ""),
        "activity_dataset": str(part.get("activity_dataset") or ""),
        "part_evidence_grade": evidence_grade,
        "target_evidence_grade": target.promoter.evidence_grade,
        "analysis_tier": target.promoter.analysis_tier,
        "reasons": reasons,
    }


def _direction_matches(observed: str, required: TerminatorDirection) -> bool:
    observed = observed.strip().lower()
    if required == "forward_or_bidirectional":
        return observed in {"forward", "bidirectional"}
    return observed == required


def evaluate_terminator_candidate(
    part: dict[str, Any],
    target: CassetteExpressionTarget,
) -> dict[str, Any]:
    hard, reasons = _candidate_hard_status(part, "terminator")
    observed_direction = str(part.get("direction") or "unknown").lower()
    direction_status = "PASS" if _direction_matches(observed_direction, target.terminator.direction) else "FAIL"
    if direction_status == "FAIL":
        hard += 1
        reasons.append("direction_mismatch")
    observed_dataset = str(part.get("activity_dataset") or "")
    observed_efficiency_raw = float(part.get("activity_value", -1.0) or -1.0)
    activity_unit = str(part.get("activity_unit") or "")
    unit_text = activity_unit.strip().lower()
    if "percent" in unit_text or "%" in unit_text:
        observed_efficiency = observed_efficiency_raw / 100.0
    elif 0.0 <= observed_efficiency_raw <= 1.0:
        observed_efficiency = observed_efficiency_raw
    else:
        observed_efficiency = -1.0
        hard += 1
        reasons.append("unsupported_efficiency_unit")
    dataset_status = "PASS" if not target.terminator.dataset_id or observed_dataset == target.terminator.dataset_id else "FAIL"
    if dataset_status == "FAIL":
        hard += 1
        reasons.append("incomparable_activity_dataset")
    if target.terminator.min_efficiency is not None:
        strength_status = "PASS" if observed_efficiency >= target.terminator.min_efficiency else "FAIL"
        violation_loss = max(0.0, float(target.terminator.min_efficiency) - observed_efficiency)
        ranking_loss = abs(1.0 - observed_efficiency) if observed_efficiency >= 0 else 2.0
        target_low = float(target.terminator.min_efficiency) * 100.0
        target_high = 100.0
        observed_percentile = observed_efficiency * 100.0
    else:
        tier = _tier_evaluation(float(part.get("activity_percentile", -1.0) or -1.0), target.terminator.strength_tier)
        strength_status = str(tier["status"])
        violation_loss = float(tier["violation_loss"])
        ranking_loss = float(tier["ranking_loss"])
        target_low = float(tier["target_low"])
        target_high = float(tier["target_high"])
        observed_percentile = float(tier["observed"])
    evidence_grade = str(part.get("evidence_grade") or "D").upper()
    evidence_status = "PASS" if evidence_grade in {"A", "B"} else "REVIEW" if evidence_grade == "C" else "FAIL"
    if evidence_status == "FAIL":
        hard += 1
        reasons.append("evidence_ineligible")
    status = "PASS" if hard == 0 and strength_status == "PASS" else "FAIL"
    return {
        "role": "terminator",
        "part_id": str(part.get("part_id") or ""),
        "status": status,
        "hard_violations": hard,
        "intent_violations": int(strength_status != "PASS"),
        "violation_loss": violation_loss,
        "ranking_loss": ranking_loss,
        "direction_status": direction_status,
        "dataset_status": dataset_status,
        "strength_status": strength_status,
        "target_low": target_low,
        "target_high": target_high,
        "activity_percentile": observed_percentile,
        "activity_value": observed_efficiency_raw,
        "activity_unit": activity_unit,
        "activity_dataset": observed_dataset,
        "part_evidence_grade": evidence_grade,
        "target_evidence_grade": target.terminator.evidence_grade,
        "analysis_tier": target.terminator.analysis_tier,
        "reasons": reasons,
    }


def _evidence_score(grade: str) -> int:
    return {"A": 4, "B": 3, "C": 1, "D": 0}.get(str(grade).upper(), 0)


def _optimize_role(
    *,
    role: Literal["promoter", "terminator"],
    candidates: list[dict[str, Any]],
    targets: list[CassetteExpressionTarget],
) -> dict[str, Any]:
    if not candidates:
        raise ValueError(f"no {role} candidates supplied")
    ordered_targets = sorted(targets, key=lambda row: row.cassette_index)
    ordered_candidates = sorted(candidates, key=lambda row: str(row.get("part_id") or ""))
    evaluator = evaluate_promoter_candidate if role == "promoter" else evaluate_terminator_candidate
    evaluations = {
        (target.cassette_index, str(part.get("part_id") or "")): evaluator(part, target)
        for target in ordered_targets
        for part in ordered_candidates
    }
    model = cp_model.CpModel()
    choice: dict[tuple[int, int], Any] = {}
    hard_terms: list[Any] = []
    intent_terms: list[Any] = []
    violation_loss_terms: list[Any] = []
    ranking_loss_terms: list[Any] = []
    evidence_terms: list[Any] = []
    score_terms: list[Any] = []

    for target_index, target in enumerate(ordered_targets):
        target_variables = []
        for candidate_index, part in enumerate(ordered_candidates):
            variable = model.new_bool_var(f"{role}_{target.cassette_index}_{candidate_index}")
            choice[(target_index, candidate_index)] = variable
            target_variables.append(variable)
            evaluation = evaluations[(target.cassette_index, str(part.get("part_id") or ""))]
            hard_terms.append(int(evaluation["hard_violations"]) * variable)
            intent_terms.append(int(evaluation["intent_violations"]) * variable)
            violation_loss_terms.append(quantize_cost(evaluation["violation_loss"]) * variable)
            ranking_loss_terms.append(quantize_cost(evaluation["ranking_loss"]) * variable)
            evidence_terms.append(_evidence_score(str(part.get("evidence_grade") or "D")) * variable)
            score_terms.append(int(part.get("score") or 0) * variable)
        model.add(sum(target_variables) == 1)

    duplicate_terms: list[Any] = []
    if any(target.avoid_reused_regulatory_parts for target in ordered_targets):
        for candidate_index, part in enumerate(ordered_candidates):
            use_count = sum(
                choice[(target_index, candidate_index)]
                for target_index in range(len(ordered_targets))
            )
            duplicate_count = model.new_int_var(
                0,
                max(0, len(ordered_targets) - 1),
                f"{role}_duplicate_{candidate_index}_{part.get('part_id')}",
            )
            model.add_max_equality(duplicate_count, [0, use_count - 1])
            duplicate_terms.append(duplicate_count)

    objectives = [
        sum(hard_terms),
        sum(intent_terms),
        sum(violation_loss_terms),
        sum(duplicate_terms),
        sum(ranking_loss_terms),
        -sum(evidence_terms),
        -sum(score_terms),
    ]
    tie_break_objectives = [
        sum(
            candidate_index * choice[(target_index, candidate_index)]
            for candidate_index in range(len(ordered_candidates))
        )
        for target_index in range(len(ordered_targets))
    ]
    solver, optimizer = solve_lexicographic_cp_sat(
        model=model,
        objectives=objectives,
        tie_break_objectives=tie_break_objectives,
        candidate_space=candidate_space_size(
            [len(ordered_candidates)] * len(ordered_targets)
        ),
        binary_variable_count=len(ordered_targets) * len(ordered_candidates),
    )

    selected: dict[int, dict[str, Any]] = {}
    selected_evaluations: dict[int, dict[str, Any]] = {}
    selected_rows: list[dict[str, Any]] = []
    selected_ids: list[str] = []
    for target_index, target in enumerate(ordered_targets):
        for candidate_index, part in enumerate(ordered_candidates):
            if solver.value(choice[(target_index, candidate_index)]):
                selected[target.cassette_index] = part
                evaluation = evaluations[
                    (target.cassette_index, str(part.get("part_id") or ""))
                ]
                selected_evaluations[target.cassette_index] = evaluation
                selected_rows.append(evaluation)
                selected_ids.append(str(part.get("part_id") or ""))
                break
        if target.cassette_index not in selected:
            raise RuntimeError(
                f"CP-SAT returned no selected {role} for cassette {target.cassette_index}"
            )

    objective = (
        sum(int(row["hard_violations"]) for row in selected_rows),
        sum(int(row["intent_violations"]) for row in selected_rows),
        round(sum(float(row["violation_loss"]) for row in selected_rows), 12),
        (
            len(selected_ids) - len(set(selected_ids))
            if any(target.avoid_reused_regulatory_parts for target in ordered_targets)
            else 0
        ),
        round(sum(float(row["ranking_loss"]) for row in selected_rows), 12),
        -sum(
            _evidence_score(str(part.get("evidence_grade") or "D"))
            for part in selected.values()
        ),
        -sum(int(part.get("score") or 0) for part in selected.values()),
        tuple(selected_ids),
    )
    return {
        "selected": selected,
        "selected_evaluations": selected_evaluations,
        "candidate_evaluations": [
            {"cassette_index": cassette_index, **evaluation}
            for (cassette_index, _), evaluation in sorted(evaluations.items())
        ],
        "objective": list(objective[:-1]) + [list(objective[-1])],
        "optimizer": optimizer.to_dict(),
    }


def optimize_expression_regulatory_parts(
    *,
    cassettes: list[dict[str, Any]],
    promoter_candidates: list[dict[str, Any]],
    terminator_candidates: list[dict[str, Any]],
    expression_targets: list[CassetteExpressionTarget],
) -> dict[str, Any]:
    cassette_ids = {int(row["cassette_index"]) for row in cassettes}
    target_ids = {row.cassette_index for row in expression_targets}
    if cassette_ids != target_ids:
        raise ValueError(
            f"expression target coverage mismatch; missing={sorted(cassette_ids-target_ids)}, "
            f"extra={sorted(target_ids-cassette_ids)}"
        )
    promoter = _optimize_role(role="promoter", candidates=promoter_candidates, targets=expression_targets)
    terminator = _optimize_role(role="terminator", candidates=terminator_candidates, targets=expression_targets)
    return {
        "promoter_by_cassette": promoter["selected"],
        "terminator_by_cassette": terminator["selected"],
        "promoter_evaluations": promoter["selected_evaluations"],
        "terminator_evaluations": terminator["selected_evaluations"],
        "candidate_evaluations": {
            "promoter": promoter["candidate_evaluations"],
            "terminator": terminator["candidate_evaluations"],
        },
        "objective": {"promoter": promoter["objective"], "terminator": terminator["objective"]},
        "optimizer": {
            "promoter": promoter["optimizer"],
            "terminator": terminator["optimizer"],
        },
    }


def group_expression_targets_by_target(
    targets: list[CassetteExpressionTarget],
) -> dict[str, list[CassetteExpressionTarget]]:
    grouped: dict[str, list[CassetteExpressionTarget]] = defaultdict(list)
    for row in targets:
        grouped[row.target_id].append(row)
    return {key: sorted(value, key=lambda row: row.cassette_index) for key, value in grouped.items()}
