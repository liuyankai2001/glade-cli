"""Scoring and deterministic reports for the P11.1 RetroPath benchmark."""

from __future__ import annotations

import csv
import io
import json
import math
import re
import statistics
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from src.pathway_analyze.retropath_benchmark_models import (
    BENCHMARK_RUN_SCHEMA,
    BENCHMARK_TASK_SCHEMA,
    BenchmarkCase,
    BenchmarkDataset,
    finite_float,
    load_benchmark_dataset,
    sha256_file,
    split_ids,
)


CASE_RESULTS_FILE_NAME = "case_results.jsonl"
CASE_METRICS_FILE_NAME = "case_metrics.csv"
SUMMARY_METRICS_FILE_NAME = "summary_metrics.csv"
FAILURE_FUNNEL_FILE_NAME = "failure_funnel.csv"
REPORT_SUMMARY_FILE_NAME = "report_summary.json"
REPORT_FILE_NAME = "report.md"


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid JSON artifact: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"JSON artifact must be an object: {path}")
    return payload


def _truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes"}


def _pipeline_artifact_path(
    pipeline: Mapping[str, Any],
    name: str,
) -> Path | None:
    artifacts = pipeline.get("artifacts")
    if not isinstance(artifacts, Mapping):
        return None
    record = artifacts.get(name)
    if isinstance(record, Mapping):
        value = record.get("path")
    else:
        value = record
    if not value:
        return None
    path = Path(str(value)).expanduser()
    return path.resolve() if path.is_absolute() else path


def _raw_gold_rule_rank(
    pipeline: Mapping[str, Any],
    gold_mnxr_ids: Sequence[str],
) -> int | None:
    raw_dir = _pipeline_artifact_path(pipeline, "raw_directory")
    rules_path = _pipeline_artifact_path(pipeline, "retropath_rules")
    if raw_dir is None or rules_path is None or not rules_path.is_file():
        return None
    scope_path = next(
        (
            raw_dir / name
            for name in ("target_scope.csv", "scope.csv", "results.csv")
            if (raw_dir / name).is_file()
        ),
        None,
    )
    if scope_path is None:
        return None
    observed_rule_ids = {
        rule_id.strip("[] '\"")
        for row in _read_csv(scope_path)
        for rule_id in split_ids(row.get("Rule ID"))
        if rule_id.strip("[] '\"")
    }
    observed_sources: set[str] = set()
    for row in _read_csv(rules_path):
        if str(row.get("Rule ID") or "").strip() not in observed_rule_ids:
            continue
        observed_sources.update(
            match.upper()
            for value in row.values()
            for match in re.findall(r"MNXR\d+", str(value or ""), re.IGNORECASE)
        )
    return 1 if set(gold_mnxr_ids) <= observed_sources else None


def _candidate_connectivity_ranks(
    pipeline: Mapping[str, Any],
    gold_mnxr_ids: Sequence[str],
) -> tuple[int | None, int | None]:
    routes_path = _pipeline_artifact_path(pipeline, "candidate_routes")
    steps_path = _pipeline_artifact_path(pipeline, "candidate_steps")
    if (
        routes_path is None
        or steps_path is None
        or not routes_path.is_file()
        or not steps_path.is_file()
    ):
        return None, None
    routes = _read_csv(routes_path)
    steps = _read_csv(steps_path)
    rank_by_candidate = {
        str(row.get("candidate_id") or "").strip(): _as_int(
            row.get("candidate_rank")
        )
        for row in routes
    }
    sources_by_candidate: dict[str, set[str]] = defaultdict(set)
    for row in steps:
        candidate_id = str(row.get("candidate_id") or "").strip()
        sources_by_candidate[candidate_id].update(
            item.upper() for item in split_ids(row.get("source_reaction_ids"))
        )
    gold = set(gold_mnxr_ids)
    connectivity = sorted(
        rank_by_candidate[candidate_id]
        for candidate_id, sources in sources_by_candidate.items()
        if rank_by_candidate.get(candidate_id, 0) > 0 and gold <= sources
    )
    exact_candidates = {
        str(row.get("candidate_id") or "").strip()
        for row in routes
        if not _truth(row.get("stereo_review_required"))
        and str(row.get("structure_match_quality") or "exact").strip()
        == "exact"
    }
    stereo_resolved = sorted(
        rank_by_candidate[candidate_id]
        for candidate_id, sources in sources_by_candidate.items()
        if candidate_id in exact_candidates
        and rank_by_candidate.get(candidate_id, 0) > 0
        and gold <= sources
    )
    return (
        connectivity[0] if connectivity else None,
        stereo_resolved[0] if stereo_resolved else None,
    )


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return default


def _candidate_rank_covering_gold(
    hypotheses: Sequence[Mapping[str, Any]],
    gold_mnxr_ids: Sequence[str],
    *,
    require_balanced: bool,
) -> int | None:
    gold = set(gold_mnxr_ids)
    sources_by_rank: dict[int, set[str]] = defaultdict(set)
    for row in hypotheses:
        if require_balanced and (
            str(row.get("balance_status") or "").strip() != "balanced"
            or str(row.get("cofactor_reconstruction_status") or "").strip()
            not in {"complete", "not_applicable"}
        ):
            continue
        rank = _as_int(row.get("candidate_rank"))
        source = str(row.get("source_mnxr_id") or "").strip().upper()
        if rank > 0 and source:
            sources_by_rank[rank].add(source)
    matches = [rank for rank, sources in sources_by_rank.items() if gold <= sources]
    return min(matches) if matches else None


def _strict_rank(
    summaries: Sequence[Mapping[str, Any]],
    hypotheses: Sequence[Mapping[str, Any]],
    gold_mnxr_ids: Sequence[str],
) -> int | None:
    by_id = {
        str(row.get("hypothesis_id") or "").strip(): str(
            row.get("source_mnxr_id") or ""
        ).strip().upper()
        for row in hypotheses
    }
    gold = set(gold_mnxr_ids)
    matches: list[int] = []
    for row in summaries:
        if str(row.get("validation_status") or "").strip() != "PASS_STRICT_ROUTE_FLUX":
            continue
        sources = {
            by_id[hypothesis_id]
            for hypothesis_id in split_ids(row.get("stoichiometry_hypothesis_ids"))
            if hypothesis_id in by_id
        }
        rank = _as_int(row.get("candidate_rank"))
        if rank > 0 and gold <= sources:
            matches.append(rank)
    return min(matches) if matches else None


def _formal_exact_rank(
    gap_dir: Path,
    case: BenchmarkCase,
) -> tuple[int | None, int | None]:
    solution_rows = _read_csv(gap_dir / "solutions.csv")
    step_rows = _read_csv(gap_dir / "all_solution_steps.csv")
    gold_reactions = set(case.gold_reaction_ids)
    gold_mnxr = set(case.gold_mnxr_ids)
    exact: list[tuple[int, int]] = []
    for solution in solution_rows:
        if str(solution.get("solution_source") or "").strip().lower() != "retropath":
            continue
        solution_id = _as_int(solution.get("solution_id"))
        rank = _as_int(solution.get("retropath_candidate_rank"))
        if solution_id < 1 or rank < 1:
            continue
        rows = [
            row
            for row in step_rows
            if _as_int(row.get("solution_id")) == solution_id
            and str(row.get("step_source") or "").strip() == "retropath"
        ]
        mapped_reactions = {
            item
            for row in rows
            if _truth(row.get("formal_mapping_exact"))
            for item in split_ids(row.get("exact_kegg_reaction_ids"))
        }
        mapped_mnxr = {
            str(row.get("source_mnxr_id") or "").strip().upper()
            for row in rows
            if _truth(row.get("formal_mapping_exact"))
        }
        if gold_reactions <= mapped_reactions and gold_mnxr <= mapped_mnxr:
            exact.append((rank, solution_id))
    if not exact:
        return None, None
    return min(exact)


def score_core_artifacts(
    case: BenchmarkCase,
    *,
    pipeline_result_path: str | Path,
    validation_dir: str | Path | None,
    gap_dir: str | Path,
) -> dict[str, Any]:
    """Score one case/profile without exposing gold data to the search code."""

    pipeline = _read_json(Path(pipeline_result_path))
    candidate_count = _as_int(pipeline.get("candidate_count"))
    scope_hit = bool(
        pipeline.get("ok") is True
        and _as_int(pipeline.get("sink_match_count")) > 0
        and bool(pipeline.get("scope_present"))
    )
    result: dict[str, Any] = {
        "pipeline_status": str(pipeline.get("status") or "unknown"),
        "scope_hit": scope_hit,
        "scope_present": bool(pipeline.get("scope_present")),
        "sink_match_count": _as_int(pipeline.get("sink_match_count")),
        "complete_path_count": _as_int(pipeline.get("complete_path_count")),
        "candidate_count": candidate_count,
        "rejection_count": _as_int(pipeline.get("rejection_count")),
        "upstream_enumeration_truncated": bool(
            pipeline.get("upstream_enumeration_truncated")
        ),
        "candidate_top_k_truncated": bool(
            pipeline.get("candidate_top_k_truncated")
        ),
        "raw_gold_rule_rank": _raw_gold_rule_rank(
            pipeline,
            case.gold_mnxr_ids,
        ),
        "connectivity_gold_rank": None,
        "stereo_resolved_gold_rank": None,
        "gold_template_rank": None,
        "balanced_gold_rank": None,
        "strict_gem_gold_rank": None,
        "formal_exact_gold_rank": None,
        "formal_exact_solution_id": None,
        "review_burden_before_gold": None,
        "strict_pass_combination_count": 0,
        "validation_combination_count": 0,
        "enzyme_ec_rank": None,
        "enzyme_accession_rank": None,
        "enzyme_evaluation_status": "not_run",
    }
    (
        result["connectivity_gold_rank"],
        result["stereo_resolved_gold_rank"],
    ) = _candidate_connectivity_ranks(pipeline, case.gold_mnxr_ids)
    if validation_dir is None:
        return result
    resolved_validation = Path(validation_dir)
    hypotheses = _read_csv(resolved_validation / "stoichiometry_hypotheses.csv")
    summaries = _read_csv(resolved_validation / "gem_validation_summary.csv")
    result["gold_template_rank"] = _candidate_rank_covering_gold(
        hypotheses,
        case.gold_mnxr_ids,
        require_balanced=False,
    )
    result["balanced_gold_rank"] = _candidate_rank_covering_gold(
        hypotheses,
        case.gold_mnxr_ids,
        require_balanced=True,
    )
    result["strict_gem_gold_rank"] = _strict_rank(
        summaries,
        hypotheses,
        case.gold_mnxr_ids,
    )
    result["validation_combination_count"] = len(
        [row for row in summaries if str(row.get("combination_id") or "").strip()]
    )
    result["strict_pass_combination_count"] = len(
        [
            row
            for row in summaries
            if str(row.get("validation_status") or "").strip()
            == "PASS_STRICT_ROUTE_FLUX"
        ]
    )
    exact_rank, solution_id = _formal_exact_rank(Path(gap_dir), case)
    result["formal_exact_gold_rank"] = exact_rank
    result["formal_exact_solution_id"] = solution_id
    if result["gold_template_rank"] is not None:
        result["review_burden_before_gold"] = max(
            0, int(result["gold_template_rank"]) - 1
        )
    return result


def score_enzyme_artifacts(
    case: BenchmarkCase,
    *,
    candidates_path: str | Path,
    source_unavailable: bool = False,
) -> dict[str, Any]:
    rows = _read_csv(Path(candidates_path))
    gold_ec = set(case.gold_ec_numbers)
    gold_accessions = {item.upper() for item in case.gold_uniprot_accessions}
    ec_ranks: list[int] = []
    accession_ranks: list[int] = []
    for row in rows:
        rank = _as_int(row.get("candidate_rank") or row.get("evaluation_rank"))
        if rank < 1:
            continue
        if gold_ec.intersection(split_ids(row.get("ec_numbers") or row.get("ec_number"))):
            ec_ranks.append(rank)
        accession = str(row.get("accession") or "").strip().upper()
        if accession and accession in gold_accessions:
            accession_ranks.append(rank)
    return {
        "enzyme_evaluation_status": (
            "source_unavailable" if source_unavailable else "completed"
        ),
        "enzyme_ec_rank": min(ec_ranks) if ec_ranks else None,
        "enzyme_accession_rank": min(accession_ranks) if accession_ranks else None,
        "gold_accession_available": bool(gold_accessions),
        "enzyme_candidate_count": len(rows),
    }


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    p = successes / total
    denominator = 1.0 + (z * z / total)
    centre = (p + z * z / (2.0 * total)) / denominator
    margin = (
        z
        * math.sqrt((p * (1.0 - p) / total) + z * z / (4.0 * total * total))
        / denominator
    )
    return max(0.0, centre - margin), min(1.0, centre + margin)


def _rank_hit(metrics: Mapping[str, Any], field: str, k: int) -> bool:
    rank = _as_int(metrics.get(field))
    return rank > 0 and rank <= k


def _percentile_summary(values: Sequence[float]) -> dict[str, float | None]:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return {"median": None, "q1": None, "q3": None}
    if len(finite) == 1:
        return {"median": finite[0], "q1": finite[0], "q3": finite[0]}
    q1, _, q3 = statistics.quantiles(finite, n=4, method="inclusive")
    return {"median": statistics.median(finite), "q1": q1, "q3": q3}


def _task_is_evaluable(task: Mapping[str, Any]) -> bool:
    return str(task.get("status") or "") == "completed"


def aggregate_results(
    tasks: Sequence[Mapping[str, Any]],
    *,
    top_k: Sequence[int],
) -> dict[str, Any]:
    profiles = sorted({str(task.get("profile") or "") for task in tasks})
    summary_rows: list[dict[str, Any]] = []
    profile_summaries: dict[str, Any] = {}
    rank_fields = {
        "raw_gold_rule_recall": "raw_gold_rule_rank",
        "connectivity_gold_recall": "connectivity_gold_rank",
        "stereo_resolved_gold_recall": "stereo_resolved_gold_rank",
        "gold_template_recall": "gold_template_rank",
        "balanced_gold_recall": "balanced_gold_rank",
        "strict_gem_gold_recall": "strict_gem_gold_rank",
        "formal_exact_gold_recall": "formal_exact_gold_rank",
        "enzyme_ec_recall": "enzyme_ec_rank",
        "enzyme_accession_recall": "enzyme_accession_rank",
    }
    for profile in profiles:
        selected = [task for task in tasks if task.get("profile") == profile]
        evaluable = [task for task in selected if _task_is_evaluable(task)]
        profile_summary: dict[str, Any] = {
            "selected_count": len(selected),
            "evaluable_count": len(evaluable),
            "operational_completion_rate": (
                len(evaluable) / len(selected) if selected else None
            ),
        }
        for denominator_name, values in (
            ("all_selected", selected),
            ("evaluable", evaluable),
        ):
            for metric_name, field in rank_fields.items():
                metric_values = values
                if metric_name == "enzyme_accession_recall":
                    metric_values = [
                        task
                        for task in values
                        if bool(task.get("metrics", {}).get("gold_accession_available"))
                    ]
                if metric_name.startswith("enzyme_"):
                    metric_values = [
                        task
                        for task in metric_values
                        if task.get("metrics", {}).get("enzyme_evaluation_status")
                        in {"completed", "source_unavailable"}
                    ]
                for k in top_k:
                    successes = sum(
                        _rank_hit(task.get("metrics", {}), field, k)
                        for task in metric_values
                    )
                    total = len(metric_values)
                    lower, upper = wilson_interval(successes, total)
                    summary_rows.append(
                        {
                            "profile": profile,
                            "denominator": denominator_name,
                            "metric": metric_name,
                            "k": k,
                            "numerator": successes,
                            "denominator_count": total,
                            "value": successes / total if total else None,
                            "ci95_lower": lower if total else None,
                            "ci95_upper": upper if total else None,
                        }
                    )
                reciprocal_sum = sum(
                    1.0 / rank
                    for task in metric_values
                    if (rank := _as_int(task.get("metrics", {}).get(field))) > 0
                )
                summary_rows.append(
                    {
                        "profile": profile,
                        "denominator": denominator_name,
                        "metric": f"{metric_name}_mrr",
                        "k": None,
                        "numerator": reciprocal_sum,
                        "denominator_count": len(metric_values),
                        "value": (
                            reciprocal_sum / len(metric_values)
                            if metric_values
                            else None
                        ),
                        "ci95_lower": None,
                        "ci95_upper": None,
                    }
                )
        completion = len(evaluable) / len(selected) if selected else None
        scope_hits = sum(bool(task.get("metrics", {}).get("scope_hit")) for task in evaluable)
        runtimes = [
            finite_float(task.get("runtime_seconds"))
            for task in evaluable
            if finite_float(task.get("runtime_seconds")) is not None
        ]
        burdens = [
            float(task["metrics"]["review_burden_before_gold"])
            for task in evaluable
            if task.get("metrics", {}).get("review_burden_before_gold") is not None
        ]
        profile_summary.update(
            {
                "operational_completion_rate": completion,
                "scope_hit_rate_evaluable": (
                    scope_hits / len(evaluable) if evaluable else None
                ),
                "runtime_seconds": _percentile_summary(runtimes),
                "review_burden": _percentile_summary(burdens),
                "truncation_rate_evaluable": (
                    sum(
                        bool(task.get("metrics", {}).get("upstream_enumeration_truncated"))
                        or bool(task.get("metrics", {}).get("candidate_top_k_truncated"))
                        for task in evaluable
                    )
                    / len(evaluable)
                    if evaluable
                    else None
                ),
            }
        )
        profile_summaries[profile] = profile_summary

    by_case: dict[str, dict[str, Mapping[str, Any]]] = defaultdict(dict)
    for task in tasks:
        by_case[str(task.get("case_id"))][str(task.get("profile"))] = task
    paired: dict[str, Any] = {
        "pair_count": 0,
        "controlled_only_exact_recovery": 0,
        "full_a0_only_exact_recovery": 0,
        "both_exact_recovery": 0,
        "neither_exact_recovery": 0,
        "rank_deltas_full_minus_controlled": [],
    }
    for profiles_by_case in by_case.values():
        controlled = profiles_by_case.get("controlled")
        full_a0 = profiles_by_case.get("full_a0")
        if not controlled or not full_a0:
            continue
        if not (_task_is_evaluable(controlled) and _task_is_evaluable(full_a0)):
            continue
        paired["pair_count"] += 1
        c_rank = _as_int(controlled.get("metrics", {}).get("formal_exact_gold_rank"))
        f_rank = _as_int(full_a0.get("metrics", {}).get("formal_exact_gold_rank"))
        if c_rank and f_rank:
            paired["both_exact_recovery"] += 1
            paired["rank_deltas_full_minus_controlled"].append(f_rank - c_rank)
        elif c_rank:
            paired["controlled_only_exact_recovery"] += 1
        elif f_rank:
            paired["full_a0_only_exact_recovery"] += 1
        else:
            paired["neither_exact_recovery"] += 1
    paired["rank_delta_summary"] = _percentile_summary(
        [float(value) for value in paired.pop("rank_deltas_full_minus_controlled")]
    )
    return {
        "profile_summaries": profile_summaries,
        "summary_rows": summary_rows,
        "paired_profile_comparison": paired,
    }


def _csv_text(rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(columns), lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                key: "" if row.get(key) is None else row.get(key)
                for key in columns
            }
        )
    return output.getvalue()


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(text)
            temporary = Path(handle.name)
        temporary.replace(path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _format_rate(value: Any) -> str:
    number = finite_float(value)
    return "N/A" if number is None else f"{100.0 * number:.1f}%"


def _summary_value(
    rows: Sequence[Mapping[str, Any]],
    *,
    profile: str,
    metric: str,
    k: int,
    denominator: str = "evaluable",
) -> Mapping[str, Any] | None:
    return next(
        (
            row
            for row in rows
            if row.get("profile") == profile
            and row.get("metric") == metric
            and _as_int(row.get("k")) == k
            and row.get("denominator") == denominator
        ),
        None,
    )


def render_markdown_report(
    run: Mapping[str, Any],
    tasks: Sequence[Mapping[str, Any]],
    aggregate: Mapping[str, Any],
) -> str:
    lines = [
        "# RetroPath P11.1 隐藏 KEGG 反应基线报告",
        "",
        f"- Benchmark：`{run.get('benchmark_id', '')}`",
        f"- Run ID：`{run.get('run_id', '')}`",
        f"- 创建时间：`{run.get('created_at', '')}`",
        f"- 数据集 SHA-256：`{run.get('dataset_sha256', '')}`",
        f"- Git commit：`{run.get('git', {}).get('commit', '')}`",
        f"- 工作树状态：`{'dirty' if run.get('git', {}).get('dirty') else 'clean'}`",
        "",
        "> 本报告是 P11.1 正向恢复能力测试。金标准反应的来源规则仍保留在 RR02 中，",
        "> 因此结果表示可恢复上限，不是未知反应假阳性率；来源规则排除属于 P11.2。",
        "",
        "## 执行概况",
        "",
        "| Profile | 计划任务 | 可评测任务 | 完成率 | Scope 命中率 | 截断率 | 运行时间中位数 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    profile_summaries = aggregate.get("profile_summaries", {})
    for profile in sorted(profile_summaries):
        item = profile_summaries[profile]
        runtime = item.get("runtime_seconds", {}).get("median")
        lines.append(
            "| "
            + " | ".join(
                (
                    profile,
                    str(item.get("selected_count", 0)),
                    str(item.get("evaluable_count", 0)),
                    _format_rate(item.get("operational_completion_rate")),
                    _format_rate(item.get("scope_hit_rate_evaluable")),
                    _format_rate(item.get("truncation_rate_evaluable")),
                    "N/A" if runtime is None else f"{float(runtime):.1f}s",
                )
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## 核心恢复率",
            "",
            "下表使用可评测案例作为分母；CSV 同时保留 `all_selected` 保守下限。",
            "",
            "| Profile | 指标 | @1 | @3 | @5 | @10 |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    summary_rows = aggregate.get("summary_rows", [])
    metric_labels = (
        ("raw_gold_rule_recall", "原始 RR02 规则出现"),
        ("connectivity_gold_recall", "结构连通路线恢复"),
        ("stereo_resolved_gold_recall", "立体身份精确恢复"),
        ("gold_template_recall", "来源模板恢复"),
        ("balanced_gold_recall", "平衡计量恢复"),
        ("strict_gem_gold_recall", "严格 GEM 恢复"),
        ("formal_exact_gold_recall", "精确正式 solution"),
    )
    for profile in sorted(profile_summaries):
        for metric, label in metric_labels:
            values = []
            for k in (1, 3, 5, 10):
                row = _summary_value(
                    summary_rows,
                    profile=profile,
                    metric=metric,
                    k=k,
                )
                if not row or not row.get("denominator_count"):
                    values.append("N/A")
                else:
                    values.append(
                        f"{_format_rate(row.get('value'))} "
                        f"({row.get('numerator')}/{row.get('denominator_count')})"
                    )
            lines.append(f"| {profile} | {label} | " + " | ".join(values) + " |")
    paired = aggregate.get("paired_profile_comparison", {})
    lines.extend(
        [
            "",
            "## 首次命中排名质量",
            "",
            "| Profile | 模板 MRR | 平衡 MRR | GEM MRR | 精确 solution MRR |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for profile in sorted(profile_summaries):
        mrr_values = []
        for metric in (
            "gold_template_recall_mrr",
            "balanced_gold_recall_mrr",
            "strict_gem_gold_recall_mrr",
            "formal_exact_gold_recall_mrr",
        ):
            row = next(
                (
                    item
                    for item in summary_rows
                    if item.get("profile") == profile
                    and item.get("denominator") == "evaluable"
                    and item.get("metric") == metric
                ),
                None,
            )
            value = finite_float(row.get("value")) if row else None
            mrr_values.append("N/A" if value is None else f"{value:.3f}")
        lines.append(f"| {profile} | " + " | ".join(mrr_values) + " |")
    lines.extend(
        [
            "",
            "## 双 Profile 配对比较",
            "",
            f"- 可比较案例：{paired.get('pair_count', 0)}",
            f"- 两种 profile 均精确恢复：{paired.get('both_exact_recovery', 0)}",
            f"- 仅 controlled 恢复：{paired.get('controlled_only_exact_recovery', 0)}",
            f"- 仅 full_a0 恢复：{paired.get('full_a0_only_exact_recovery', 0)}",
            f"- 两种 profile 均未恢复：{paired.get('neither_exact_recovery', 0)}",
            "",
            "## 单案例结果",
            "",
            "| Case | EC 类 | Profile | 状态 | Scope | 原始规则 | 连通排名 | 立体精确 | 模板排名 | 平衡排名 | GEM 排名 | 精确排名 | 耗时 |",
            "|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for task in sorted(tasks, key=lambda item: (str(item.get("case_id")), str(item.get("profile")))):
        metrics = task.get("metrics", {})
        values = []
        for field in (
            "raw_gold_rule_rank",
            "connectivity_gold_rank",
            "stereo_resolved_gold_rank",
            "gold_template_rank",
            "balanced_gold_rank",
            "strict_gem_gold_rank",
            "formal_exact_gold_rank",
        ):
            value = metrics.get(field)
            values.append("-" if value is None else str(value))
        runtime = finite_float(task.get("runtime_seconds"))
        lines.append(
            f"| {task.get('case_id', '')} | {task.get('ec_class', '')} | "
            f"{task.get('profile', '')} | {task.get('status', '')} | "
            f"{'是' if metrics.get('scope_hit') else '否'} | "
            + " | ".join(values)
            + f" | {'-' if runtime is None else f'{runtime:.1f}s'} |"
        )
    failures = [task for task in tasks if not _task_is_evaluable(task)]
    lines.extend(["", "## 失败与限制", ""])
    if not failures:
        lines.append("- 没有基础设施或数据层失败。")
    else:
        for task in sorted(failures, key=lambda item: (str(item.get("case_id")), str(item.get("profile")))):
            error = task.get("error", {})
            lines.append(
                f"- `{task.get('case_id')}/{task.get('profile')}`："
                f"{error.get('stage', 'unknown')} / {error.get('code', 'unknown')} — "
                f"{error.get('detail', '')}"
            )
    lines.extend(
        [
            "",
            "本试运行只有 12 个反应案例。比例使用 Wilson 95% 区间，运行时间使用中位数和 IQR；",
            "不进行显著性检验，也不基于本报告自动修改生产参数。",
            "",
        ]
    )
    return "\n".join(lines)


def _verified_tasks(run_dir: Path, run: Mapping[str, Any]) -> list[dict[str, Any]]:
    records = run.get("tasks")
    if not isinstance(records, list):
        raise ValueError("benchmark run manifest tasks must be an array")
    tasks: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ValueError(f"run task record {index} is invalid")
        relative = Path(str(record.get("result_path") or ""))
        path = relative if relative.is_absolute() else run_dir / relative
        if not path.is_file():
            raise FileNotFoundError(f"benchmark task result is missing: {path}")
        expected = str(record.get("sha256") or "")
        if sha256_file(path) != expected:
            raise ValueError(f"benchmark task result checksum mismatch: {path}")
        task = _read_json(path)
        if task.get("schema_version") != BENCHMARK_TASK_SCHEMA:
            raise ValueError(f"unsupported benchmark task schema: {path}")
        if task.get("case_id") != record.get("case_id") or task.get("profile") != record.get("profile"):
            raise ValueError(f"benchmark task identity mismatch: {path}")
        artifacts = task.get("artifacts")
        if not isinstance(artifacts, Mapping):
            raise ValueError(f"benchmark task artifacts are invalid: {path}")
        for name, artifact in artifacts.items():
            if not isinstance(artifact, Mapping):
                raise ValueError(f"benchmark artifact record is invalid: {name}")
            artifact_path = Path(str(artifact.get("path") or ""))
            artifact_path = (
                artifact_path
                if artifact_path.is_absolute()
                else run_dir / artifact_path
            )
            if not artifact_path.is_file():
                raise FileNotFoundError(
                    f"benchmark task artifact is missing: {artifact_path}"
                )
            if sha256_file(artifact_path) != artifact.get("sha256"):
                raise ValueError(
                    f"benchmark task artifact checksum mismatch: {artifact_path}"
                )
        tasks.append(task)
    return tasks


def generate_benchmark_report(
    run_dir: str | Path,
    *,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    resolved_run = Path(run_dir).expanduser().resolve()
    manifest_path = resolved_run / "benchmark_run_manifest.json"
    run = _read_json(manifest_path)
    if run.get("schema_version") != BENCHMARK_RUN_SCHEMA:
        raise ValueError("unsupported benchmark run schema")
    dataset_path = Path(str(run.get("dataset_path") or "")).expanduser().resolve()
    dataset = load_benchmark_dataset(dataset_path)
    if dataset.dataset_sha256 != run.get("dataset_sha256"):
        raise ValueError("benchmark dataset changed after the run")
    tasks = _verified_tasks(resolved_run, run)
    aggregate = aggregate_results(tasks, top_k=dataset.defaults.top_k)

    case_rows: list[dict[str, Any]] = []
    case_results_lines: list[str] = []
    for task in sorted(tasks, key=lambda item: (str(item.get("case_id")), str(item.get("profile")))):
        case_results_lines.append(
            json.dumps(task, ensure_ascii=False, sort_keys=True, allow_nan=False)
        )
        metrics = task.get("metrics", {})
        case_rows.append(
            {
                "case_id": task.get("case_id"),
                "ec_class": task.get("ec_class"),
                "profile": task.get("profile"),
                "status": task.get("status"),
                "outcome": task.get("outcome"),
                "runtime_seconds": task.get("runtime_seconds"),
                "scope_hit": metrics.get("scope_hit"),
                "candidate_count": metrics.get("candidate_count"),
                "raw_gold_rule_rank": metrics.get("raw_gold_rule_rank"),
                "connectivity_gold_rank": metrics.get("connectivity_gold_rank"),
                "stereo_resolved_gold_rank": metrics.get(
                    "stereo_resolved_gold_rank"
                ),
                "gold_template_rank": metrics.get("gold_template_rank"),
                "balanced_gold_rank": metrics.get("balanced_gold_rank"),
                "strict_gem_gold_rank": metrics.get("strict_gem_gold_rank"),
                "formal_exact_gold_rank": metrics.get("formal_exact_gold_rank"),
                "review_burden_before_gold": metrics.get("review_burden_before_gold"),
                "upstream_enumeration_truncated": metrics.get("upstream_enumeration_truncated"),
                "candidate_top_k_truncated": metrics.get("candidate_top_k_truncated"),
                "enzyme_evaluation_status": metrics.get("enzyme_evaluation_status"),
                "enzyme_ec_rank": metrics.get("enzyme_ec_rank"),
                "enzyme_accession_rank": metrics.get("enzyme_accession_rank"),
                "error_stage": task.get("error", {}).get("stage"),
                "error_code": task.get("error", {}).get("code"),
            }
        )

    failure_counts = Counter(
        (
            str(task.get("profile") or ""),
            str(task.get("error", {}).get("stage") or task.get("outcome") or "unknown"),
            str(task.get("error", {}).get("code") or ""),
        )
        for task in tasks
        if not _task_is_evaluable(task)
    )
    funnel_rows = [
        {"profile": profile, "stage": stage, "code": code, "count": count}
        for (profile, stage, code), count in sorted(failure_counts.items())
    ]
    funnel_metric_fields = (
        ("selected", None),
        ("operational_completed", "__completed__"),
        ("raw_gold_rule", "raw_gold_rule_rank"),
        ("connectivity", "connectivity_gold_rank"),
        ("stereo_resolved", "stereo_resolved_gold_rank"),
        ("balanced", "balanced_gold_rank"),
        ("strict_gem", "strict_gem_gold_rank"),
        ("formal_exact", "formal_exact_gold_rank"),
    )
    for profile in sorted({str(task.get("profile") or "") for task in tasks}):
        selected = [task for task in tasks if task.get("profile") == profile]
        for stage, field in funnel_metric_fields:
            if field is None:
                count = len(selected)
            elif field == "__completed__":
                count = sum(_task_is_evaluable(task) for task in selected)
            else:
                count = sum(
                    _as_int(task.get("metrics", {}).get(field)) > 0
                    for task in selected
                )
            funnel_rows.append(
                {
                    "profile": profile,
                    "stage": stage,
                    "code": "passed",
                    "count": count,
                }
            )
    funnel_rows.sort(
        key=lambda row: (
            str(row.get("profile") or ""),
            str(row.get("stage") or ""),
            str(row.get("code") or ""),
        )
    )
    markdown = render_markdown_report(run, tasks, aggregate)
    summary = {
        "schema_version": "retropath_benchmark_report.v1",
        "benchmark_id": dataset.benchmark_id,
        "run_id": run.get("run_id"),
        "dataset_sha256": dataset.dataset_sha256,
        "task_count": len(tasks),
        **aggregate,
    }
    artifact_texts = {
        CASE_RESULTS_FILE_NAME: "\n".join(case_results_lines) + ("\n" if case_results_lines else ""),
        CASE_METRICS_FILE_NAME: _csv_text(case_rows, tuple(case_rows[0]) if case_rows else ("case_id",)),
        SUMMARY_METRICS_FILE_NAME: _csv_text(
            aggregate["summary_rows"],
            (
                "profile",
                "denominator",
                "metric",
                "k",
                "numerator",
                "denominator_count",
                "value",
                "ci95_lower",
                "ci95_upper",
            ),
        ),
        FAILURE_FUNNEL_FILE_NAME: _csv_text(
            funnel_rows,
            ("profile", "stage", "code", "count"),
        ),
        REPORT_SUMMARY_FILE_NAME: json.dumps(
            summary,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        + "\n",
        REPORT_FILE_NAME: markdown,
    }
    artifacts: dict[str, dict[str, str]] = {}
    for name, text in artifact_texts.items():
        path = resolved_run / name
        _atomic_write(path, text)
        artifacts[name] = {"path": str(path), "sha256": sha256_file(path)}
    if output_path is not None:
        resolved_output = Path(output_path).expanduser().resolve()
        _atomic_write(resolved_output, markdown)
        artifacts["external_report"] = {
            "path": str(resolved_output),
            "sha256": sha256_file(resolved_output),
        }
    return {
        "ok": True,
        "run_id": run.get("run_id"),
        "task_count": len(tasks),
        "artifacts": artifacts,
        "profile_summaries": aggregate["profile_summaries"],
        "paired_profile_comparison": aggregate["paired_profile_comparison"],
    }


__all__ = [
    "CASE_METRICS_FILE_NAME",
    "CASE_RESULTS_FILE_NAME",
    "FAILURE_FUNNEL_FILE_NAME",
    "REPORT_FILE_NAME",
    "REPORT_SUMMARY_FILE_NAME",
    "SUMMARY_METRICS_FILE_NAME",
    "aggregate_results",
    "generate_benchmark_report",
    "render_markdown_report",
    "score_core_artifacts",
    "score_enzyme_artifacts",
    "wilson_interval",
]
