"""Orchestrate P5 merging and write isolated RetroPath candidate artifacts."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from src.pathway_analyze.expand_chassis_metabolites import ExpansionBundle
from src.pathway_analyze.kegg_gap_analyze import KeggRestClient
from src.pathway_analyze.retropath_merge import (
    DEFAULT_MAX_CANDIDATES,
    DEFAULT_MAX_NEW_ENZYMES,
    DEFAULT_MAX_TOTAL_STEPS,
    DEFAULT_MAX_WITNESS_PLANS,
    HybridCandidateRoute,
    RetroPathMergeRejection,
    RetroPathMergeResult,
    merge_retropath_candidates,
)
from src.pathway_analyze.retropath_routes import RetroPathEnumerationResult

CANDIDATE_ROUTES_FILE_NAME = "candidate_routes.csv"
CANDIDATE_STEPS_FILE_NAME = "candidate_steps.csv"
REJECTED_ROUTES_FILE_NAME = "rejected_routes.csv"

CANDIDATE_ROUTE_COLUMNS = (
    "candidate_rank",
    "candidate_id",
    "source_retrosynthetic_path_id",
    "target_compound_id",
    "sink_kegg_ids",
    "sink_depths",
    "sink_inchikeys",
    "kegg_prefix_reaction_ids",
    "retropath_step_ids",
    "retropath_reaction_option_ids",
    "kegg_prefix_steps",
    "retropath_steps",
    "total_steps",
    "maximum_sink_depth",
    "minimum_rule_specificity",
    "worst_rule_score",
    "score_semantics",
    "contains_auxiliary_fragments",
    "structure_match_quality",
    "stereo_review_required",
    "route_source",
    "contains_predicted_steps",
    "validation_status",
    "review_required",
    "upstream_enumeration_truncated",
    "candidate_top_k_truncated",
)

CANDIDATE_STEP_COLUMNS = (
    "candidate_id",
    "step_index",
    "step_id",
    "step_source",
    "status",
    "orientation",
    "direction",
    "reaction_option_ids",
    "reaction_smiles",
    "substrate_compound_ids",
    "product_compound_ids",
    "substrate_stoichiometry_json",
    "product_stoichiometry_json",
    "depends_on_step_ids",
    "source_transformation_ids",
    "sink_anchor_kegg_ids",
    "expansion_depth",
    "is_endogenous",
    "rule_ids",
    "source_reaction_ids",
    "source_ec_numbers",
    "source_uniprot_ids",
    "minimum_rule_specificity",
    "worst_rule_score",
    "score_semantics",
    "balance_status",
    "cofactor_reconstruction_status",
)

REJECTED_ROUTE_COLUMNS = (
    "source_stage",
    "source_path_id",
    "reason_code",
    "reason_detail",
    "sink_kegg_ids",
    "compound_id",
    "transformation_id",
)


@dataclass(frozen=True)
class RetroPathCandidateArtifacts:
    """Written P5 files, checksums, and the in-memory merge result."""

    merge_result: RetroPathMergeResult
    output_dir: Path
    candidate_routes_path: Path
    candidate_steps_path: Path
    rejected_routes_path: Path
    candidate_routes_sha256: str
    candidate_steps_sha256: str
    rejected_routes_sha256: str

    @property
    def candidate_count(self) -> int:
        return self.merge_result.candidate_count

    @property
    def rejection_count(self) -> int:
        return len(self.merge_result.rejections)


def _compact_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _render_csv(
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=tuple(fieldnames),
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def _atomic_write_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
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
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _join(values: Iterable[str]) -> str:
    return ";".join(str(value) for value in values)


def _candidate_route_row(
    rank: int,
    candidate: HybridCandidateRoute,
    *,
    upstream_truncated: bool,
    top_k_truncated: bool,
) -> dict[str, Any]:
    return {
        "candidate_rank": rank,
        "candidate_id": candidate.candidate_id,
        "source_retrosynthetic_path_id": (candidate.source_retrosynthetic_path_id),
        "target_compound_id": candidate.target_compound_id,
        "sink_kegg_ids": _join(
            item.representative_kegg_id for item in candidate.sink_matches
        ),
        "sink_depths": _join(
            f"{item.representative_kegg_id}:{item.minimum_depth}"
            for item in candidate.sink_matches
        ),
        "sink_inchikeys": _join(item.inchikey for item in candidate.sink_matches),
        "kegg_prefix_reaction_ids": _join(candidate.kegg_prefix_reaction_ids),
        "retropath_step_ids": _join(candidate.retropath_step_ids),
        "retropath_reaction_option_ids": _join(candidate.retropath_reaction_option_ids),
        "kegg_prefix_steps": candidate.kegg_prefix_steps,
        "retropath_steps": candidate.retropath_steps,
        "total_steps": candidate.total_steps,
        "maximum_sink_depth": candidate.maximum_sink_depth,
        "minimum_rule_specificity": candidate.minimum_rule_specificity,
        "worst_rule_score": candidate.worst_rule_score,
        "score_semantics": candidate.score_semantics,
        "contains_auxiliary_fragments": str(
            candidate.contains_auxiliary_fragments
        ).lower(),
        "structure_match_quality": candidate.structure_match_quality,
        "stereo_review_required": str(
            candidate.stereo_review_required
        ).lower(),
        "route_source": candidate.route_source,
        "contains_predicted_steps": "true",
        "validation_status": candidate.validation_status,
        "review_required": str(candidate.review_required).lower(),
        "upstream_enumeration_truncated": str(upstream_truncated).lower(),
        "candidate_top_k_truncated": str(top_k_truncated).lower(),
    }


def _candidate_step_rows(
    candidate: HybridCandidateRoute,
) -> Iterable[dict[str, Any]]:
    for index, step in enumerate(candidate.steps, start=1):
        yield {
            "candidate_id": candidate.candidate_id,
            "step_index": index,
            "step_id": step.step_id,
            "step_source": step.step_source,
            "status": step.status,
            "orientation": step.orientation,
            "direction": step.direction,
            "reaction_option_ids": _join(step.reaction_option_ids),
            "reaction_smiles": step.reaction_smiles,
            "substrate_compound_ids": _join(step.substrate_compound_ids),
            "product_compound_ids": _join(step.product_compound_ids),
            "substrate_stoichiometry_json": _compact_json(
                [list(item) for item in step.substrate_stoichiometry]
            ),
            "product_stoichiometry_json": _compact_json(
                [list(item) for item in step.product_stoichiometry]
            ),
            "depends_on_step_ids": _join(step.depends_on_step_ids),
            "source_transformation_ids": _join(step.source_transformation_ids),
            "sink_anchor_kegg_ids": _join(step.sink_anchor_kegg_ids),
            "expansion_depth": step.expansion_depth,
            "is_endogenous": (
                "" if step.is_endogenous is None else str(step.is_endogenous).lower()
            ),
            "rule_ids": _join(step.rule_ids),
            "source_reaction_ids": _join(step.source_reaction_ids),
            "source_ec_numbers": _join(step.source_ec_numbers),
            "source_uniprot_ids": _join(step.source_uniprot_ids),
            "minimum_rule_specificity": (
                ""
                if step.minimum_rule_specificity is None
                else step.minimum_rule_specificity
            ),
            "worst_rule_score": (
                "" if step.worst_rule_score is None else step.worst_rule_score
            ),
            "score_semantics": step.score_semantics or "",
            "balance_status": step.balance_status,
            "cofactor_reconstruction_status": (step.cofactor_reconstruction_status),
        }


def _rejected_route_row(
    rejection: RetroPathMergeRejection,
) -> dict[str, Any]:
    return {
        "source_stage": rejection.source_stage,
        "source_path_id": rejection.source_path_id or "",
        "reason_code": rejection.reason_code,
        "reason_detail": rejection.reason_detail,
        "sink_kegg_ids": _join(rejection.sink_kegg_ids),
        "compound_id": rejection.compound_id or "",
        "transformation_id": rejection.transformation_id or "",
    }


def write_retropath_candidate_artifacts(
    merge_result: RetroPathMergeResult,
    output_dir: str | Path,
) -> RetroPathCandidateArtifacts:
    """Write deterministic P5 CSVs without touching official KEGG outputs."""

    if not isinstance(merge_result, RetroPathMergeResult):
        raise ValueError("merge_result must be a RetroPathMergeResult")
    resolved_output_dir = Path(output_dir).expanduser().resolve()
    candidate_routes_path = resolved_output_dir / CANDIDATE_ROUTES_FILE_NAME
    candidate_steps_path = resolved_output_dir / CANDIDATE_STEPS_FILE_NAME
    rejected_routes_path = resolved_output_dir / REJECTED_ROUTES_FILE_NAME

    route_text = _render_csv(
        CANDIDATE_ROUTE_COLUMNS,
        (
            _candidate_route_row(
                rank,
                candidate,
                upstream_truncated=merge_result.upstream_truncated,
                top_k_truncated=merge_result.truncated,
            )
            for rank, candidate in enumerate(merge_result.candidates, start=1)
        ),
    )
    step_text = _render_csv(
        CANDIDATE_STEP_COLUMNS,
        (
            row
            for candidate in merge_result.candidates
            for row in _candidate_step_rows(candidate)
        ),
    )
    rejected_text = _render_csv(
        REJECTED_ROUTE_COLUMNS,
        (_rejected_route_row(item) for item in merge_result.rejections),
    )
    return RetroPathCandidateArtifacts(
        merge_result=merge_result,
        output_dir=resolved_output_dir,
        candidate_routes_path=candidate_routes_path,
        candidate_steps_path=candidate_steps_path,
        rejected_routes_path=rejected_routes_path,
        candidate_routes_sha256=_atomic_write_text(
            candidate_routes_path,
            route_text,
        ),
        candidate_steps_sha256=_atomic_write_text(
            candidate_steps_path,
            step_text,
        ),
        rejected_routes_sha256=_atomic_write_text(
            rejected_routes_path,
            rejected_text,
        ),
    )


def analyze_retropath_candidates(
    enumeration_result: RetroPathEnumerationResult,
    expansion_bundle: ExpansionBundle,
    kegg_client: KeggRestClient,
    output_dir: str | Path,
    *,
    max_candidates: int = DEFAULT_MAX_CANDIDATES,
    max_witness_plans: int = DEFAULT_MAX_WITNESS_PLANS,
    max_total_steps: int = DEFAULT_MAX_TOTAL_STEPS,
    max_new_enzymes: int = DEFAULT_MAX_NEW_ENZYMES,
    ignored_common_compounds: set[str] | None = None,
) -> RetroPathCandidateArtifacts:
    """Merge P4 routes and immediately write the isolated candidate files."""

    result = merge_retropath_candidates(
        enumeration_result,
        expansion_bundle,
        kegg_client,
        max_candidates=max_candidates,
        max_witness_plans=max_witness_plans,
        max_total_steps=max_total_steps,
        max_new_enzymes=max_new_enzymes,
        ignored_common_compounds=ignored_common_compounds,
    )
    return write_retropath_candidate_artifacts(result, output_dir)


__all__ = [
    "CANDIDATE_ROUTES_FILE_NAME",
    "CANDIDATE_ROUTE_COLUMNS",
    "CANDIDATE_STEPS_FILE_NAME",
    "CANDIDATE_STEP_COLUMNS",
    "REJECTED_ROUTES_FILE_NAME",
    "REJECTED_ROUTE_COLUMNS",
    "RetroPathCandidateArtifacts",
    "analyze_retropath_candidates",
    "write_retropath_candidate_artifacts",
]
