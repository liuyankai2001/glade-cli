from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass
from functools import reduce
from operator import mul
from typing import Any, Iterable, Sequence

from ortools.sat.python import cp_model

from src.config.expression_config import EXPRESSION_OPTIMIZER_CONFIG

# 兼容原有公共常量；数值含义和唯一默认值位于 expression_config.py。
COST_SCALE = EXPRESSION_OPTIMIZER_CONFIG.cost_scale
DEFAULT_COMPONENT_TIME_LIMIT_SECONDS = EXPRESSION_OPTIMIZER_CONFIG.component_time_limit_seconds
DEFAULT_RANDOM_SEED = EXPRESSION_OPTIMIZER_CONFIG.random_seed


class DiscreteExpressionOptimizerError(RuntimeError):
    """Base error for deterministic expression-part optimization."""


class DiscreteExpressionOptimizerInfeasibleError(DiscreteExpressionOptimizerError):
    """The discrete model has no feasible assignment."""


class DiscreteExpressionOptimizerTimeoutError(DiscreteExpressionOptimizerError):
    """The solver did not prove optimality within the configured deadline."""


@dataclass(frozen=True)
class OptimizerStats:
    backend: str
    status: str
    optimality_proven: bool
    candidate_space_size: str
    binary_variable_count: int
    constraint_count: int
    component_count: int
    branches_explored: int
    conflicts: int
    wall_time_seconds: float

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["wall_time_seconds"] = round(float(payload["wall_time_seconds"]), 6)
        return payload


def quantize_cost(value: float, *, scale: int = COST_SCALE) -> int:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"optimizer cost must be finite: {value}")
    return int(round(numeric * int(scale)))


def candidate_space_size(candidate_counts: Iterable[int]) -> int:
    counts = [int(value) for value in candidate_counts]
    if any(value <= 0 for value in counts):
        raise ValueError(f"candidate counts must be positive: {counts}")
    return reduce(mul, counts, 1)


def direct_optimizer_stats(
    *,
    candidate_counts: Iterable[int],
    component_count: int,
    wall_time_seconds: float,
) -> OptimizerStats:
    counts = [int(value) for value in candidate_counts]
    return OptimizerStats(
        backend="direct",
        status="OPTIMAL",
        optimality_proven=True,
        candidate_space_size=str(candidate_space_size(counts)),
        binary_variable_count=0,
        constraint_count=0,
        component_count=int(component_count),
        branches_explored=0,
        conflicts=0,
        wall_time_seconds=float(wall_time_seconds),
    )


def combine_optimizer_stats(
    *,
    stats: Sequence[OptimizerStats],
    candidate_counts: Iterable[int],
    component_count: int,
    binary_variable_count: int,
    wall_time_seconds: float,
) -> OptimizerStats:
    rows = list(stats)
    backend = "direct" if not rows else "cp_sat_componentized"
    return OptimizerStats(
        backend=backend,
        status="OPTIMAL",
        optimality_proven=True,
        candidate_space_size=str(candidate_space_size(candidate_counts)),
        binary_variable_count=int(binary_variable_count),
        constraint_count=sum(row.constraint_count for row in rows),
        component_count=int(component_count),
        branches_explored=sum(row.branches_explored for row in rows),
        conflicts=sum(row.conflicts for row in rows),
        wall_time_seconds=float(wall_time_seconds),
    )


def solve_lexicographic_cp_sat(
    *,
    model: cp_model.CpModel,
    objectives: Sequence[cp_model.LinearExpr | int],
    tie_break_objectives: Sequence[cp_model.LinearExpr | int],
    candidate_space: int,
    binary_variable_count: int,
    component_count: int = 1,
    time_limit_seconds: float = DEFAULT_COMPONENT_TIME_LIMIT_SECONDS,
    random_seed: int = DEFAULT_RANDOM_SEED,
) -> tuple[cp_model.CpSolver, OptimizerStats]:
    """
    Solve a CP-SAT model with exact staged lexicographic minimization.

    Each optimum is fixed before the next objective is introduced. A result is
    returned only when every stage is proven OPTIMAL.
    """

    phases = list(objectives) + list(tie_break_objectives)
    if not phases:
        phases = [0]

    started = time.perf_counter()
    deadline = started + float(time_limit_seconds)
    total_branches = 0
    total_conflicts = 0
    solver: cp_model.CpSolver | None = None

    for phase_index, objective in enumerate(phases, start=1):
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            raise DiscreteExpressionOptimizerTimeoutError(
                f"CP-SAT time limit exhausted before objective phase {phase_index}/{len(phases)}"
            )

        model.minimize(objective)
        phase_solver = cp_model.CpSolver()
        phase_solver.parameters.num_search_workers = 1
        phase_solver.parameters.random_seed = int(random_seed)
        phase_solver.parameters.max_time_in_seconds = max(0.001, remaining)
        status = phase_solver.solve(model)
        status_name = phase_solver.status_name(status)
        total_branches += int(phase_solver.num_branches)
        total_conflicts += int(phase_solver.num_conflicts)

        if status == cp_model.INFEASIBLE:
            raise DiscreteExpressionOptimizerInfeasibleError(
                f"CP-SAT model is infeasible at objective phase {phase_index}/{len(phases)}"
            )
        if status != cp_model.OPTIMAL:
            raise DiscreteExpressionOptimizerTimeoutError(
                "CP-SAT did not prove optimality "
                f"at objective phase {phase_index}/{len(phases)}; status={status_name}"
            )

        optimum = int(round(float(phase_solver.objective_value)))
        model.add(objective == optimum)
        solver = phase_solver

    if solver is None:
        raise DiscreteExpressionOptimizerError("CP-SAT returned no solver result")

    elapsed = time.perf_counter() - started
    stats = OptimizerStats(
        backend="cp_sat",
        status="OPTIMAL",
        optimality_proven=True,
        candidate_space_size=str(int(candidate_space)),
        binary_variable_count=int(binary_variable_count),
        constraint_count=len(model.proto.constraints),
        component_count=int(component_count),
        branches_explored=total_branches,
        conflicts=total_conflicts,
        wall_time_seconds=elapsed,
    )
    return solver, stats


__all__ = [
    "COST_SCALE",
    "DEFAULT_COMPONENT_TIME_LIMIT_SECONDS",
    "DiscreteExpressionOptimizerError",
    "DiscreteExpressionOptimizerInfeasibleError",
    "DiscreteExpressionOptimizerTimeoutError",
    "OptimizerStats",
    "candidate_space_size",
    "combine_optimizer_stats",
    "cp_model",
    "direct_optimizer_stats",
    "quantize_cost",
    "solve_lexicographic_cp_sat",
]
