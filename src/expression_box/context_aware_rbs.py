from __future__ import annotations

import hashlib
import json
import math
import time
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, model_validator

from src.tools.expression_cassette_assembly_tools.discrete_expression_optimizer import (
    OptimizerStats,
    candidate_space_size,
    combine_optimizer_stats,
    cp_model,
    quantize_cost,
    solve_lexicographic_cp_sat,
)


TargetKind = Literal["absolute_interval", "ordinal_band"]
TargetTier = Literal["low", "medium", "high"]
AnalysisTier = Literal["primary_gold", "exploratory"]
RelativeKind = Literal["balanced_group", "group_greater_than"]


class RegulatoryTarget(BaseModel):
    target_id: str = Field(default="", description="设计目标 ID；用于审计追踪")
    cassette_index: int = Field(ge=1, description="CDS 所在表达盒编号")
    cds_id: str = Field(min_length=1, description="目标 CDS ID")
    target_kind: TargetKind = Field(description="absolute_interval 或 ordinal_band")
    lower_tir: float | None = Field(default=None, gt=0)
    upper_tir: float | None = Field(default=None, gt=0)
    target_tier: TargetTier | None = None
    evidence_id: str = Field(default="", description="文献或作者锁定证据 ID")
    evidence_grade: str = Field(default="author_prior")
    analysis_tier: AnalysisTier = Field(default="exploratory")

    @model_validator(mode="after")
    def validate_target(self) -> "RegulatoryTarget":
        if self.target_kind == "absolute_interval":
            if self.lower_tir is None or self.upper_tir is None:
                raise ValueError("absolute_interval requires lower_tir and upper_tir")
            if self.lower_tir > self.upper_tir:
                raise ValueError("lower_tir cannot exceed upper_tir")
            if self.target_tier is not None:
                raise ValueError("absolute_interval cannot set target_tier")
        else:
            if self.target_tier is None:
                raise ValueError("ordinal_band requires target_tier")
            if self.lower_tir is not None or self.upper_tir is not None:
                raise ValueError("ordinal_band cannot set lower_tir or upper_tir")
        return self


class RelativeExpressionConstraint(BaseModel):
    constraint_id: str = Field(min_length=1)
    target_id: str = Field(default="")
    constraint_kind: RelativeKind
    left_cds_ids: list[str] = Field(min_length=1)
    right_cds_ids: list[str] = Field(default_factory=list)
    min_fold: float = Field(default=1.0, gt=0)
    max_fold: float = Field(default=2.0, gt=1)
    evidence_id: str = Field(default="")
    evidence_grade: str = Field(default="author_prior")
    analysis_tier: AnalysisTier = Field(default="exploratory")

    @model_validator(mode="after")
    def validate_constraint(self) -> "RelativeExpressionConstraint":
        self.left_cds_ids = [str(value).strip() for value in self.left_cds_ids if str(value).strip()]
        self.right_cds_ids = [str(value).strip() for value in self.right_cds_ids if str(value).strip()]
        if self.constraint_kind == "balanced_group":
            if len(self.left_cds_ids) < 2:
                raise ValueError("balanced_group requires at least two left_cds_ids")
            if self.right_cds_ids:
                raise ValueError("balanced_group does not use right_cds_ids")
        elif not self.right_cds_ids:
            raise ValueError("group_greater_than requires right_cds_ids")
        overlap = set(self.left_cds_ids) & set(self.right_cds_ids)
        if overlap:
            raise ValueError(f"relative constraint groups overlap: {sorted(overlap)}")
        return self


def _stable_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def normalize_dna(value: Any) -> str:
    sequence = "".join(
        line.strip()
        for line in str(value or "").splitlines()
        if line.strip() and not line.startswith(">")
    ).upper().replace("U", "T")
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("sequence must contain only A/C/G/T")
    return sequence


def predict_tir(sequence: str, start_bp: int, name: str) -> float:
    warnings.filterwarnings("ignore", message="RBS Calculator Vienna is missing dependency ViennaRNA")
    try:
        import RNA  # noqa: F401 - validates the ViennaRNA Python binding, not only package metadata
        import ostir
    except ImportError as exc:
        raise RuntimeError(
            "context-aware RBS selection requires ostir==1.1.3 and ViennaRNA==2.7.0"
        ) from exc

    results = ostir.run_ostir(
        normalize_dna(sequence),
        start=int(start_bp),
        end=int(start_bp),
        name=name,
        threads=1,
        decimal_places=6,
        circular=False,
        verbosity=0,
    )
    exact = [row for row in results if int(row.get("start_position") or 0) == int(start_bp)]
    if len(exact) != 1:
        raise ValueError(f"OSTIR returned {len(exact)} exact predictions at start {start_bp}")
    expression = float(exact[0]["expression"])
    if not math.isfinite(expression) or expression <= 0:
        raise ValueError(f"OSTIR returned invalid expression value: {expression}")
    return expression


def _part_selection(
    part: dict[str, Any],
    *,
    cassette_index: int,
    role: str,
    cds_id: str | None = None,
) -> dict[str, Any]:
    sequence_file = part.get("sequence_file")
    path = sequence_file.get("path") if isinstance(sequence_file, dict) else sequence_file
    if not path:
        raise ValueError(f"candidate lacks sequence_file: {part.get('part_id')}")
    payload: dict[str, Any] = {
        "cassette_index": cassette_index,
        "role": role,
        "part_id": str(part.get("part_id") or ""),
        "sequence_file": {"path": str(path)},
    }
    if cds_id:
        payload["cds_id"] = cds_id
    return payload


def _part_for_cassette(part_or_map: dict[Any, Any], cassette_index: int) -> dict[str, Any]:
    """Resolve either one shared part or a cassette-indexed part assignment."""
    if "part_id" in part_or_map:
        return part_or_map
    selected = part_or_map.get(cassette_index, part_or_map.get(str(cassette_index)))
    if not isinstance(selected, dict) or not selected.get("part_id"):
        raise ValueError(f"missing regulatory part for cassette {cassette_index}")
    return selected


def evaluate_candidate_contexts(
    *,
    manifest: dict[str, Any],
    output_dir: Path,
    cassettes: list[dict[str, Any]],
    promoter: dict[str, Any],
    terminator: dict[str, Any],
    baseline_rbs: dict[str, Any],
    rbs_candidates: list[dict[str, Any]],
    predictor: Any = predict_tir,
) -> list[dict[str, Any]]:
    """Predict every eligible RBS in the exact cassette assembly context."""

    from src.tools.expression_cassette_assembly_tools.assemble_expression_cassette import (
        _cds_by_id,
        assemble_cassette_sequence,
    )

    cds_items = _cds_by_id(manifest)
    rows: list[dict[str, Any]] = []
    for cassette in sorted(cassettes, key=lambda item: int(item["cassette_index"])):
        cassette_index = int(cassette["cassette_index"])
        cassette_promoter = _part_for_cassette(promoter, cassette_index)
        cassette_terminator = _part_for_cassette(terminator, cassette_index)
        cds_ids = [str(value) for value in cassette["cds_ids"]]
        for cds_order, cds_id in enumerate(cds_ids, start=1):
            placement_rows: list[dict[str, Any]] = []
            for candidate in sorted(rbs_candidates, key=lambda item: str(item.get("part_id") or "")):
                selected_parts = [
                    _part_selection(cassette_promoter, cassette_index=cassette_index, role="promoter"),
                    _part_selection(cassette_terminator, cassette_index=cassette_index, role="terminator"),
                ]
                for other_cds_id in cds_ids:
                    selected_parts.append(
                        _part_selection(
                            candidate if other_cds_id == cds_id else baseline_rbs,
                            cassette_index=cassette_index,
                            role="rbs",
                            cds_id=other_cds_id,
                        )
                    )
                try:
                    sequence, components = assemble_cassette_sequence(
                        output_dir=output_dir,
                        cassette=cassette,
                        cds_items=cds_items,
                        selected_parts=selected_parts,
                    )
                    cds_components = [
                        component
                        for component in components
                        if component.get("type") == "cds" and component.get("cds_id") == cds_id
                    ]
                    if len(cds_components) != 1:
                        raise ValueError(f"expected one assembled CDS component for {cds_id}")
                    start_bp = int(cds_components[0]["start_bp"])
                    tir = float(
                        predictor(
                            sequence,
                            start_bp,
                            f"cassette_{cassette_index}_{cds_id}_{candidate['part_id']}",
                        )
                    )
                    if not math.isfinite(tir) or tir <= 0:
                        raise ValueError(f"predictor returned invalid TIR: {tir}")
                    status = "COMPLETE"
                    error = ""
                except Exception as exc:
                    start_bp = 0
                    tir = float("nan")
                    status = "INCOMPLETE"
                    error = f"{type(exc).__name__}: {exc}"
                placement_rows.append(
                    {
                        "cassette_index": cassette_index,
                        "cds_order": cds_order,
                        "cds_id": cds_id,
                        "part_id": str(candidate.get("part_id") or ""),
                        "tir": tir,
                        "prediction_status": status,
                        "prediction_error": error,
                        "cds_start_bp": start_bp,
                        "metadata_score": int(candidate.get("score") or 0),
                    }
                )
            completed = sorted(
                float(row["tir"])
                for row in placement_rows
                if row["prediction_status"] == "COMPLETE"
            )
            if not completed:
                raise ValueError(f"all RBS predictions failed for {cds_id}")
            for row in placement_rows:
                if row["prediction_status"] == "COMPLETE":
                    row["percentile"] = 100.0 * sum(value <= float(row["tir"]) for value in completed) / len(completed)
                else:
                    row["percentile"] = float("nan")
                rows.append(row)
    return rows


def _ordinal_bounds(tier: TargetTier) -> tuple[float, float]:
    return {
        "low": (0.0, 100.0 / 3.0),
        "medium": (100.0 / 3.0, 200.0 / 3.0),
        "high": (200.0 / 3.0, 100.0),
    }[tier]


def evaluate_regulatory_target(target: RegulatoryTarget, row: dict[str, Any]) -> dict[str, Any]:
    tir = float(row["tir"])
    if target.target_kind == "absolute_interval":
        low = float(target.lower_tir)
        high = float(target.upper_tir)
        centre = math.sqrt(low * high)
        status = "PASS" if low <= tir <= high else "FAIL"
        boundary_loss = 0.0 if status == "PASS" else min(abs(math.log2(tir / low)), abs(math.log2(tir / high)))
        centre_loss = abs(math.log2(tir / centre))
        return {
            "status": status,
            "violation_loss": boundary_loss,
            "ranking_loss": centre_loss,
            "target_low": low,
            "target_high": high,
            "observed": tir,
        }

    percentile = float(row["percentile"])
    low, high = _ordinal_bounds(str(target.target_tier))
    status = "PASS" if low <= percentile <= high else "FAIL"
    boundary_loss = 0.0 if status == "PASS" else min(abs(percentile - low), abs(percentile - high)) / 100.0
    centre_loss = abs(percentile - ((low + high) / 2.0)) / 100.0
    return {
        "status": status,
        "violation_loss": boundary_loss,
        "ranking_loss": centre_loss,
        "target_low": low,
        "target_high": high,
        "observed": percentile,
    }


def _geometric_mean(values: list[float]) -> float:
    if not values or any(value <= 0 or not math.isfinite(value) for value in values):
        raise ValueError("geometric mean requires positive finite values")
    return math.exp(sum(math.log(value) for value in values) / len(values))


def evaluate_relative_constraint(
    constraint: RelativeExpressionConstraint,
    selected_by_cds: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    missing = sorted(
        set(constraint.left_cds_ids + constraint.right_cds_ids) - set(selected_by_cds)
    )
    if missing:
        raise ValueError(f"relative constraint {constraint.constraint_id} references missing CDS: {missing}")
    left = [float(selected_by_cds[cds_id]["tir"]) for cds_id in constraint.left_cds_ids]
    if constraint.constraint_kind == "balanced_group":
        observed_fold = max(left) / min(left)
        status = "PASS" if observed_fold <= constraint.max_fold else "FAIL"
        loss = max(0.0, math.log2(observed_fold / constraint.max_fold))
        ranking_loss = abs(math.log2(observed_fold))
        return {
            "constraint_id": constraint.constraint_id,
            "status": status,
            "observed_fold": observed_fold,
            "required_fold": constraint.max_fold,
            "violation_loss": loss,
            "ranking_loss": ranking_loss,
        }

    right = [float(selected_by_cds[cds_id]["tir"]) for cds_id in constraint.right_cds_ids]
    observed_fold = _geometric_mean(left) / _geometric_mean(right)
    status = "PASS" if observed_fold >= constraint.min_fold else "FAIL"
    loss = max(0.0, math.log2(constraint.min_fold / observed_fold))
    ranking_loss = max(0.0, math.log2(constraint.min_fold / observed_fold))
    return {
        "constraint_id": constraint.constraint_id,
        "status": status,
        "observed_fold": observed_fold,
        "required_fold": constraint.min_fold,
        "violation_loss": loss,
        "ranking_loss": ranking_loss,
    }


def _constraint_components(
    placement_order: list[str],
    constraints: list[RelativeExpressionConstraint],
) -> list[tuple[list[str], list[RelativeExpressionConstraint]]]:
    order_index = {cds_id: index for index, cds_id in enumerate(placement_order)}
    parent = {cds_id: cds_id for cds_id in placement_order}

    def find(cds_id: str) -> str:
        while parent[cds_id] != cds_id:
            parent[cds_id] = parent[parent[cds_id]]
            cds_id = parent[cds_id]
        return cds_id

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        if order_index[left_root] <= order_index[right_root]:
            parent[right_root] = left_root
        else:
            parent[left_root] = right_root

    placement_ids = set(placement_order)
    for constraint in constraints:
        involved = list(dict.fromkeys(constraint.left_cds_ids + constraint.right_cds_ids))
        missing = sorted(set(involved) - placement_ids)
        if missing:
            raise ValueError(
                f"relative constraint {constraint.constraint_id} references missing CDS: {missing}"
            )
        for cds_id in involved[1:]:
            union(involved[0], cds_id)

    grouped_ids: dict[str, list[str]] = defaultdict(list)
    for cds_id in placement_order:
        grouped_ids[find(cds_id)].append(cds_id)

    grouped_constraints: dict[str, list[RelativeExpressionConstraint]] = defaultdict(list)
    for constraint in constraints:
        grouped_constraints[find(constraint.left_cds_ids[0])].append(constraint)

    components = [
        (cds_ids, grouped_constraints.get(find(cds_ids[0]), []))
        for cds_ids in grouped_ids.values()
    ]
    return sorted(components, key=lambda item: min(order_index[cds_id] for cds_id in item[0]))


def _single_rbs_objective(
    *,
    cds_id: str,
    row: dict[str, Any],
    target: RegulatoryTarget | None,
) -> tuple[tuple[Any, ...], dict[str, Any] | None]:
    evaluation = evaluate_regulatory_target(target, row) if target is not None else None
    return (
        (
            int(evaluation is not None and evaluation["status"] != "PASS"),
            round(float(evaluation["violation_loss"]), 12) if evaluation is not None else 0.0,
            round(float(evaluation["ranking_loss"]), 12) if evaluation is not None else 0.0,
            -int(row.get("metadata_score") or 0),
            str(row["part_id"]),
        ),
        {"cds_id": cds_id, **evaluation} if evaluation is not None else None,
    )


def _solve_rbs_component(
    *,
    cds_ids: list[str],
    grouped: dict[str, list[dict[str, Any]]],
    targets_by_cds: dict[str, RegulatoryTarget],
    constraints: list[RelativeExpressionConstraint],
) -> tuple[dict[str, dict[str, Any]], OptimizerStats]:
    model = cp_model.CpModel()
    choice: dict[tuple[str, int], Any] = {}
    log_tir_values: dict[tuple[str, int], int] = {}
    y_by_cds: dict[str, Any] = {}
    y_bounds: dict[str, tuple[int, int]] = {}
    failure_terms: list[Any] = []
    violation_loss_terms: list[Any] = []
    ranking_loss_terms: list[Any] = []
    metadata_terms: list[Any] = []

    for cds_id in cds_ids:
        rows = grouped[cds_id]
        variables = []
        for candidate_index, row in enumerate(rows):
            variable = model.new_bool_var(f"rbs_{cds_id}_{candidate_index}")
            choice[(cds_id, candidate_index)] = variable
            variables.append(variable)
            log_tir = quantize_cost(math.log2(float(row["tir"])))
            log_tir_values[(cds_id, candidate_index)] = log_tir

            target = targets_by_cds.get(cds_id)
            if target is not None:
                evaluation = evaluate_regulatory_target(target, row)
                failure_terms.append(int(evaluation["status"] != "PASS") * variable)
                violation_loss_terms.append(quantize_cost(evaluation["violation_loss"]) * variable)
                ranking_loss_terms.append(quantize_cost(evaluation["ranking_loss"]) * variable)
            metadata_terms.append(int(row.get("metadata_score") or 0) * variable)

        model.add(sum(variables) == 1)
        lower = min(log_tir_values[(cds_id, index)] for index in range(len(rows)))
        upper = max(log_tir_values[(cds_id, index)] for index in range(len(rows)))
        y_value = model.new_int_var(lower, upper, f"log2_tir_{cds_id}")
        model.add(
            y_value
            == sum(
                log_tir_values[(cds_id, index)] * choice[(cds_id, index)]
                for index in range(len(rows))
            )
        )
        y_by_cds[cds_id] = y_value
        y_bounds[cds_id] = (lower, upper)

    for constraint_index, constraint in enumerate(constraints):
        if constraint.constraint_kind == "balanced_group":
            involved = constraint.left_cds_ids
            lower = min(y_bounds[cds_id][0] for cds_id in involved)
            upper = max(y_bounds[cds_id][1] for cds_id in involved)
            max_value = model.new_int_var(lower, upper, f"balanced_max_{constraint_index}")
            min_value = model.new_int_var(lower, upper, f"balanced_min_{constraint_index}")
            model.add_max_equality(max_value, [y_by_cds[cds_id] for cds_id in involved])
            model.add_min_equality(min_value, [y_by_cds[cds_id] for cds_id in involved])
            spread_upper = upper - lower
            spread = model.new_int_var(0, spread_upper, f"balanced_spread_{constraint_index}")
            model.add(spread == max_value - min_value)
            threshold = quantize_cost(math.log2(float(constraint.max_fold)))
            failed = model.new_bool_var(f"balanced_failed_{constraint_index}")
            model.add(spread <= threshold).only_enforce_if(failed.Not())
            model.add(spread >= threshold + 1).only_enforce_if(failed)
            loss_upper = max(0, spread_upper - threshold)
            loss = model.new_int_var(0, loss_upper, f"balanced_loss_{constraint_index}")
            model.add_max_equality(loss, [0, spread - threshold])
            failure_terms.append(failed)
            violation_loss_terms.append(loss)
            ranking_loss_terms.append(spread)
            continue

        left_ids = constraint.left_cds_ids
        right_ids = constraint.right_cds_ids
        left_count = len(left_ids)
        right_count = len(right_ids)
        denominator = left_count * right_count
        delta = (
            right_count * sum(y_by_cds[cds_id] for cds_id in left_ids)
            - left_count * sum(y_by_cds[cds_id] for cds_id in right_ids)
        )
        delta_min = (
            right_count * sum(y_bounds[cds_id][0] for cds_id in left_ids)
            - left_count * sum(y_bounds[cds_id][1] for cds_id in right_ids)
        )
        delta_max = (
            right_count * sum(y_bounds[cds_id][1] for cds_id in left_ids)
            - left_count * sum(y_bounds[cds_id][0] for cds_id in right_ids)
        )
        delta_value = model.new_int_var(delta_min, delta_max, f"group_delta_{constraint_index}")
        model.add(delta_value == delta)
        threshold = quantize_cost(math.log2(float(constraint.min_fold))) * denominator
        failed = model.new_bool_var(f"group_failed_{constraint_index}")
        model.add(delta_value >= threshold).only_enforce_if(failed.Not())
        model.add(delta_value <= threshold - 1).only_enforce_if(failed)
        shortfall_upper = max(0, threshold - delta_min)
        shortfall = model.new_int_var(0, shortfall_upper, f"group_shortfall_{constraint_index}")
        model.add_max_equality(shortfall, [0, threshold - delta_value])
        if denominator == 1:
            normalized_loss = shortfall
        else:
            rounded_upper = shortfall_upper + denominator // 2
            rounded_numerator = model.new_int_var(
                0,
                rounded_upper,
                f"group_rounded_numerator_{constraint_index}",
            )
            model.add(rounded_numerator == shortfall + denominator // 2)
            normalized_loss = model.new_int_var(
                0,
                math.ceil(rounded_upper / denominator),
                f"group_loss_{constraint_index}",
            )
            model.add_division_equality(normalized_loss, rounded_numerator, denominator)
        failure_terms.append(failed)
        violation_loss_terms.append(normalized_loss)
        ranking_loss_terms.append(normalized_loss)

    objectives = [
        sum(failure_terms),
        sum(violation_loss_terms),
        sum(ranking_loss_terms),
        -sum(metadata_terms),
    ]
    tie_break_objectives = [
        sum(index * choice[(cds_id, index)] for index in range(len(grouped[cds_id])))
        for cds_id in cds_ids
    ]
    binary_variable_count = sum(len(grouped[cds_id]) for cds_id in cds_ids)
    solver, stats = solve_lexicographic_cp_sat(
        model=model,
        objectives=objectives,
        tie_break_objectives=tie_break_objectives,
        candidate_space=candidate_space_size(len(grouped[cds_id]) for cds_id in cds_ids),
        binary_variable_count=binary_variable_count,
    )
    selected: dict[str, dict[str, Any]] = {}
    for cds_id in cds_ids:
        for index, row in enumerate(grouped[cds_id]):
            if solver.value(choice[(cds_id, index)]):
                selected[cds_id] = row
                break
        if cds_id not in selected:
            raise RuntimeError(f"CP-SAT returned no selected RBS for {cds_id}")
    return selected, stats


def optimize_rbs_combinations(
    *,
    candidate_rows: list[dict[str, Any]],
    regulatory_targets: list[RegulatoryTarget],
    relative_constraints: list[RelativeExpressionConstraint] | None = None,
    strict_context_coverage: bool = True,
) -> dict[str, Any]:
    constraints = list(relative_constraints or [])
    targets_by_cds = {target.cds_id: target for target in regulatory_targets}
    if len(targets_by_cds) != len(regulatory_targets):
        raise ValueError("regulatory_targets contains duplicated cds_id")

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in candidate_rows:
        if row.get("prediction_status") == "COMPLETE" and math.isfinite(float(row.get("tir") or float("nan"))):
            grouped[str(row["cds_id"])].append(row)
    placement_order = sorted(
        grouped,
        key=lambda cds_id: (
            int(grouped[cds_id][0].get("cassette_index") or 0),
            int(grouped[cds_id][0].get("cds_order") or 0),
            cds_id,
        ),
    )
    if not placement_order:
        raise ValueError("no completed RBS predictions were supplied")
    if strict_context_coverage:
        missing_targets = sorted(set(placement_order) - set(targets_by_cds))
        extra_targets = sorted(set(targets_by_cds) - set(placement_order))
        if missing_targets or extra_targets:
            raise ValueError(
                f"strict target coverage mismatch; missing={missing_targets}, extra={extra_targets}"
            )
    for cds_id in placement_order:
        grouped[cds_id].sort(key=lambda row: str(row["part_id"]))
        if not grouped[cds_id]:
            raise ValueError(f"no completed RBS prediction for {cds_id}")

    started = time.perf_counter()
    components = _constraint_components(placement_order, constraints)
    selected_by_cds: dict[str, dict[str, Any]] = {}
    component_stats: list[OptimizerStats] = []
    cp_sat_binary_variables = 0
    for cds_ids, component_constraints in components:
        if not component_constraints:
            cds_id = cds_ids[0]
            best_row = min(
                grouped[cds_id],
                key=lambda row: _single_rbs_objective(
                    cds_id=cds_id,
                    row=row,
                    target=targets_by_cds.get(cds_id),
                )[0],
            )
            selected_by_cds[cds_id] = best_row
            continue
        component_selection, stats = _solve_rbs_component(
            cds_ids=cds_ids,
            grouped=grouped,
            targets_by_cds=targets_by_cds,
            constraints=component_constraints,
        )
        selected_by_cds.update(component_selection)
        component_stats.append(stats)
        cp_sat_binary_variables += sum(len(grouped[cds_id]) for cds_id in cds_ids)

    target_evaluations = []
    violations = 0
    total_loss = 0.0
    ranking_loss = 0.0
    for cds_id in placement_order:
        target = targets_by_cds.get(cds_id)
        if target is None:
            continue
        evaluation = evaluate_regulatory_target(target, selected_by_cds[cds_id])
        target_evaluations.append({"cds_id": cds_id, **evaluation})
        violations += evaluation["status"] != "PASS"
        total_loss += float(evaluation["violation_loss"])
        ranking_loss += float(evaluation["ranking_loss"])
    relative_evaluations = []
    for constraint in constraints:
        evaluation = evaluate_relative_constraint(constraint, selected_by_cds)
        relative_evaluations.append(evaluation)
        violations += evaluation["status"] != "PASS"
        total_loss += float(evaluation["violation_loss"])
        ranking_loss += float(evaluation["ranking_loss"])
    metadata_score = sum(int(row.get("metadata_score") or 0) for row in selected_by_cds.values())
    part_ids = tuple(str(selected_by_cds[cds_id]["part_id"]) for cds_id in placement_order)
    objective = (
        int(violations),
        round(total_loss, 12),
        round(ranking_loss, 12),
        -metadata_score,
        part_ids,
    )
    optimizer = combine_optimizer_stats(
        stats=component_stats,
        candidate_counts=[len(grouped[cds_id]) for cds_id in placement_order],
        component_count=len(components),
        binary_variable_count=cp_sat_binary_variables,
        wall_time_seconds=time.perf_counter() - started,
    )

    target_payload = [target.model_dump(mode="json") for target in regulatory_targets]
    constraint_payload = [constraint.model_dump(mode="json") for constraint in constraints]
    return {
        "selected_by_cds": selected_by_cds,
        "target_evaluations": target_evaluations,
        "relative_evaluations": relative_evaluations,
        "objective": list(objective[:-1]) + [list(objective[-1])],
        "optimizer": optimizer.to_dict(),
        "target_set_sha256": _stable_sha256(
            {"regulatory_targets": target_payload, "relative_constraints": constraint_payload}
        ),
    }


def context_aware_reselection(
    *,
    manifest: dict[str, Any],
    output_dir: Path,
    cassettes: list[dict[str, Any]],
    promoter: dict[str, Any],
    terminator: dict[str, Any],
    baseline_rbs: dict[str, Any],
    rbs_candidates: list[dict[str, Any]],
    regulatory_targets: list[RegulatoryTarget],
    relative_constraints: list[RelativeExpressionConstraint] | None = None,
    strict_context_coverage: bool = True,
    predictor: Any = predict_tir,
) -> dict[str, Any]:
    candidate_rows = evaluate_candidate_contexts(
        manifest=manifest,
        output_dir=output_dir,
        cassettes=cassettes,
        promoter=promoter,
        terminator=terminator,
        baseline_rbs=baseline_rbs,
        rbs_candidates=rbs_candidates,
        predictor=predictor,
    )
    optimized = optimize_rbs_combinations(
        candidate_rows=candidate_rows,
        regulatory_targets=regulatory_targets,
        relative_constraints=relative_constraints,
        strict_context_coverage=strict_context_coverage,
    )
    return {"candidate_predictions": candidate_rows, **optimized}
