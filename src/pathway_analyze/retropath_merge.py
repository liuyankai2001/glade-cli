"""Flip RetroPath predictions and merge them with KEGG expansion witnesses.

P4 deliberately returns a retrosynthetic reaction graph.  This module turns
that graph into a biosynthetic, dependency-aware candidate without flattening
multi-sink branches or expanding alternative RetroRules annotations into
duplicate chemical routes.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

from src.pathway_analyze.expand_chassis_metabolites import ExpansionBundle
from src.pathway_analyze.kegg_gap_analyze import (
    IGNORED_COMMON_COMPOUNDS,
    KeggRestClient,
    PlanStep,
    Solution,
    build_frontier_bridge_plans,
    directional_stoichiometry,
    materialize_frontier_solution,
)
from src.pathway_analyze.retropath_models import PredictedReaction
from src.pathway_analyze.retropath_parser import (
    ParsedTransformation,
    SinkMatch,
)
from src.pathway_analyze.retropath_routes import (
    RetroPathEnumerationResult,
    RetrosyntheticPath,
)

HYBRID_CANDIDATE_SCHEMA_VERSION = 1
DEFAULT_MAX_CANDIDATES = 5
DEFAULT_MAX_WITNESS_PLANS = 3
DEFAULT_MAX_TOTAL_STEPS = 10
DEFAULT_MAX_NEW_ENZYMES = 10

HYBRID_STEP_SOURCES = frozenset({"kegg_expansion", "retropath"})
HYBRID_STEP_STATUSES = frozenset({"endogenous", "heterologous", "predicted"})


@dataclass(frozen=True)
class HybridCandidateStep:
    """One biosynthetic KEGG or RetroPath step in a candidate reaction DAG."""

    step_id: str
    step_source: str
    status: str
    orientation: str
    direction: str
    reaction_option_ids: Tuple[str, ...]
    reaction_smiles: str
    substrate_stoichiometry: Tuple[Tuple[str, float], ...]
    product_stoichiometry: Tuple[Tuple[str, float], ...]
    depends_on_step_ids: Tuple[str, ...] = tuple()
    source_transformation_ids: Tuple[str, ...] = tuple()
    sink_anchor_kegg_ids: Tuple[str, ...] = tuple()
    expansion_depth: int = 0
    is_endogenous: Optional[bool] = None
    rule_ids: Tuple[str, ...] = tuple()
    source_reaction_ids: Tuple[str, ...] = tuple()
    source_ec_numbers: Tuple[str, ...] = tuple()
    minimum_rule_specificity: Optional[int] = None
    worst_rule_score: Optional[float] = None
    score_semantics: Optional[str] = None
    balance_status: str = "not_checked"
    cofactor_reconstruction_status: str = "not_checked"

    @property
    def substrate_compound_ids(self) -> Tuple[str, ...]:
        return tuple(item[0] for item in self.substrate_stoichiometry)

    @property
    def product_compound_ids(self) -> Tuple[str, ...]:
        return tuple(item[0] for item in self.product_stoichiometry)

    def to_dict(self) -> dict[str, Any]:
        return {
            "step_id": self.step_id,
            "step_source": self.step_source,
            "status": self.status,
            "orientation": self.orientation,
            "direction": self.direction,
            "reaction_option_ids": list(self.reaction_option_ids),
            "reaction_smiles": self.reaction_smiles,
            "substrate_stoichiometry": [
                [compound_id, amount]
                for compound_id, amount in self.substrate_stoichiometry
            ],
            "product_stoichiometry": [
                [compound_id, amount]
                for compound_id, amount in self.product_stoichiometry
            ],
            "depends_on_step_ids": list(self.depends_on_step_ids),
            "source_transformation_ids": list(self.source_transformation_ids),
            "sink_anchor_kegg_ids": list(self.sink_anchor_kegg_ids),
            "expansion_depth": self.expansion_depth,
            "is_endogenous": self.is_endogenous,
            "rule_ids": list(self.rule_ids),
            "source_reaction_ids": list(self.source_reaction_ids),
            "source_ec_numbers": list(self.source_ec_numbers),
            "minimum_rule_specificity": self.minimum_rule_specificity,
            "worst_rule_score": self.worst_rule_score,
            "score_semantics": self.score_semantics,
            "balance_status": self.balance_status,
            "cofactor_reconstruction_status": self.cofactor_reconstruction_status,
        }


@dataclass(frozen=True)
class HybridCandidateRoute:
    """One complete biosynthetic DAG from A0 to the requested target."""

    candidate_id: str
    source_retrosynthetic_path_id: str
    target_compound_id: str
    sink_matches: Tuple[SinkMatch, ...]
    steps: Tuple[HybridCandidateStep, ...]
    minimum_rule_specificity: int
    worst_rule_score: float
    score_semantics: str
    contains_auxiliary_fragments: bool
    validation_status: str = "raw"
    review_required: bool = True

    @property
    def kegg_prefix_steps(self) -> int:
        return sum(item.step_source == "kegg_expansion" for item in self.steps)

    @property
    def retropath_steps(self) -> int:
        return sum(item.step_source == "retropath" for item in self.steps)

    @property
    def total_steps(self) -> int:
        return len(self.steps)

    @property
    def kegg_prefix_reaction_ids(self) -> Tuple[str, ...]:
        return tuple(
            item.reaction_option_ids[0]
            for item in self.steps
            if item.step_source == "kegg_expansion"
        )

    @property
    def retropath_step_ids(self) -> Tuple[str, ...]:
        return tuple(
            item.step_id for item in self.steps if item.step_source == "retropath"
        )

    @property
    def retropath_reaction_option_ids(self) -> Tuple[str, ...]:
        return tuple(
            reaction_id
            for item in self.steps
            if item.step_source == "retropath"
            for reaction_id in item.reaction_option_ids
        )

    @property
    def maximum_sink_depth(self) -> int:
        return max(item.minimum_depth for item in self.sink_matches)

    @property
    def route_source(self) -> str:
        return "kegg_retropath"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": HYBRID_CANDIDATE_SCHEMA_VERSION,
            "candidate_id": self.candidate_id,
            "source_retrosynthetic_path_id": self.source_retrosynthetic_path_id,
            "target_compound_id": self.target_compound_id,
            "sink_matches": [item.to_dict() for item in self.sink_matches],
            "steps": [item.to_dict() for item in self.steps],
            "kegg_prefix_reaction_ids": list(self.kegg_prefix_reaction_ids),
            "retropath_step_ids": list(self.retropath_step_ids),
            "retropath_reaction_option_ids": list(self.retropath_reaction_option_ids),
            "kegg_prefix_steps": self.kegg_prefix_steps,
            "retropath_steps": self.retropath_steps,
            "total_steps": self.total_steps,
            "maximum_sink_depth": self.maximum_sink_depth,
            "minimum_rule_specificity": self.minimum_rule_specificity,
            "worst_rule_score": self.worst_rule_score,
            "score_semantics": self.score_semantics,
            "contains_auxiliary_fragments": self.contains_auxiliary_fragments,
            "route_source": self.route_source,
            "contains_predicted_steps": True,
            "validation_status": self.validation_status,
            "review_required": self.review_required,
        }


@dataclass(frozen=True)
class RetroPathMergeRejection:
    """One P4 or P5 reason for excluding a candidate or branch."""

    source_stage: str
    reason_code: str
    reason_detail: str
    source_path_id: Optional[str] = None
    sink_kegg_ids: Tuple[str, ...] = tuple()
    compound_id: Optional[str] = None
    transformation_id: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_stage": self.source_stage,
            "source_path_id": self.source_path_id,
            "reason_code": self.reason_code,
            "reason_detail": self.reason_detail,
            "sink_kegg_ids": list(self.sink_kegg_ids),
            "compound_id": self.compound_id,
            "transformation_id": self.transformation_id,
        }


@dataclass(frozen=True)
class RetroPathMergeResult:
    """Accepted hybrid candidates, flipped reactions, and rejection audit."""

    candidates: Tuple[HybridCandidateRoute, ...]
    biosynthetic_reactions: Tuple[PredictedReaction, ...]
    rejections: Tuple[RetroPathMergeRejection, ...]
    upstream_truncated: bool
    truncated: bool
    max_candidates: int
    max_witness_plans: int
    max_total_steps: int
    max_new_enzymes: int

    @property
    def candidate_count(self) -> int:
        return len(self.candidates)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": HYBRID_CANDIDATE_SCHEMA_VERSION,
            "candidates": [item.to_dict() for item in self.candidates],
            "biosynthetic_reactions": [
                item.to_dict() for item in self.biosynthetic_reactions
            ],
            "rejections": [item.to_dict() for item in self.rejections],
            "upstream_truncated": self.upstream_truncated,
            "truncated": self.truncated,
            "max_candidates": self.max_candidates,
            "max_witness_plans": self.max_witness_plans,
            "max_total_steps": self.max_total_steps,
            "max_new_enzymes": self.max_new_enzymes,
            "candidate_count": self.candidate_count,
        }


def _positive_int(value: Any, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{field_name} must be a positive integer")
    return value


def _nonnegative_int(value: Any, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field_name} must be a non-negative integer")
    return value


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    text = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _stable_unique(values: Iterable[str]) -> Tuple[str, ...]:
    return tuple(sorted({str(item).strip() for item in values if str(item).strip()}))


def _normalize_stoichiometry(
    values: Iterable[Tuple[str, float]],
) -> Tuple[Tuple[str, float], ...]:
    totals: dict[str, float] = defaultdict(float)
    for compound_id, amount in values:
        normalized_id = str(compound_id).strip()
        normalized_amount = float(amount)
        if not normalized_id or not math.isfinite(normalized_amount):
            raise ValueError("stoichiometry must contain IDs and finite amounts")
        if normalized_amount <= 0:
            raise ValueError("stoichiometric amounts must be positive")
        totals[normalized_id] += normalized_amount
    if not totals:
        raise ValueError("reaction side must contain at least one compound")
    return tuple(sorted((key, round(value, 12)) for key, value in totals.items()))


def _multiset_stoichiometry(values: Iterable[str]) -> Tuple[Tuple[str, float], ...]:
    counts = Counter(str(item).strip() for item in values if str(item).strip())
    return tuple(
        sorted((compound_id, float(count)) for compound_id, count in counts.items())
    )


def _reverse_reaction_smiles(reaction_smiles: str) -> str:
    if reaction_smiles.count(">>") != 1:
        raise ValueError("Reaction SMILES must contain exactly one '>>'")
    left, right = reaction_smiles.split(">>", 1)
    if not left.strip() or not right.strip():
        raise ValueError("Reaction SMILES must contain two non-empty sides")
    return f"{right.strip()}>>{left.strip()}"


def flip_predicted_reaction(reaction: PredictedReaction) -> PredictedReaction:
    """Create the biosynthetic identity of one retrosynthetic P1 reaction."""

    if not isinstance(reaction, PredictedReaction):
        raise ValueError("reaction must be a PredictedReaction")
    if reaction.orientation != "retrosynthetic":
        raise ValueError("only retrosynthetic reactions can be flipped")
    return PredictedReaction.create(
        rule_id=reaction.rule_id,
        reaction_smiles=_reverse_reaction_smiles(reaction.reaction_smiles),
        substrate_compounds=reaction.product_compounds,
        product_compounds=reaction.substrate_compounds,
        orientation="biosynthetic",
        evidence_type=reaction.evidence_type,
        source_reaction_ids=reaction.source_reaction_ids,
        source_ec_numbers=reaction.source_ec_numbers,
        source_uniprot_ids=reaction.source_uniprot_ids,
        rule_specificity=reaction.rule_specificity,
        rule_specificity_semantics=reaction.rule_specificity_semantics,
        rule_score_raw=reaction.rule_score_raw,
        score_semantics=reaction.score_semantics,
        balance_status=reaction.balance_status,
        cofactor_reconstruction_status=reaction.cofactor_reconstruction_status,
    )


def _worst_score(
    values: Sequence[float],
    score_semantics: str,
) -> float:
    if not values:
        raise ValueError("at least one rule score is required")
    if score_semantics == "lower_is_better":
        return max(values)
    if score_semantics == "higher_is_better":
        return min(values)
    raise ValueError(f"unsupported score semantics: {score_semantics}")


def _status_union(values: Iterable[str], *, incomplete_value: str) -> str:
    normalized = set(values)
    if incomplete_value in normalized:
        return incomplete_value
    if len(normalized) == 1:
        return next(iter(normalized))
    return "not_checked"


def _step_id(payload: Mapping[str, Any]) -> str:
    return "RP2STEP:" + _canonical_sha256(
        {
            "entity_type": "hybrid_candidate_step",
            "schema_version": HYBRID_CANDIDATE_SCHEMA_VERSION,
            **payload,
        }
    )


def _retropath_step(
    transformation: ParsedTransformation,
) -> tuple[HybridCandidateStep, Tuple[PredictedReaction, ...]]:
    flipped = tuple(
        sorted(
            (
                flip_predicted_reaction(item)
                for item in transformation.reaction_variants
            ),
            key=lambda item: item.reaction_id,
        )
    )
    if not flipped:
        raise ValueError("transformation has no reaction variants")
    first = flipped[0]
    expected_sides = (first.substrate_compounds, first.product_compounds)
    if any(
        (item.substrate_compounds, item.product_compounds) != expected_sides
        for item in flipped[1:]
    ):
        raise ValueError("rule variants disagree on concrete reaction sides")
    semantics = _stable_unique(
        item.score_semantics or ""
        for item in flipped
        if item.rule_score_raw is not None
    )
    if len(semantics) != 1:
        raise ValueError("rule variants disagree on score semantics")
    scores = [
        float(item.rule_score_raw)
        for item in flipped
        if item.rule_score_raw is not None
    ]
    specificities = [
        int(item.rule_specificity)
        for item in flipped
        if item.rule_specificity is not None
    ]
    if len(scores) != len(flipped) or len(specificities) != len(flipped):
        raise ValueError("rule variants lack specificity or score evidence")
    reaction_ids = tuple(item.reaction_id for item in flipped)
    substrate_stoichiometry = _multiset_stoichiometry(first.substrate_compounds)
    product_stoichiometry = _multiset_stoichiometry(first.product_compounds)
    payload = {
        "step_source": "retropath",
        "reaction_option_ids": list(reaction_ids),
        "substrate_stoichiometry": [list(item) for item in substrate_stoichiometry],
        "product_stoichiometry": [list(item) for item in product_stoichiometry],
    }
    step = HybridCandidateStep(
        step_id=_step_id(payload),
        step_source="retropath",
        status="predicted",
        orientation="biosynthetic",
        direction="biosynthetic",
        reaction_option_ids=reaction_ids,
        reaction_smiles=first.reaction_smiles,
        substrate_stoichiometry=substrate_stoichiometry,
        product_stoichiometry=product_stoichiometry,
        source_transformation_ids=(transformation.transformation_id,),
        rule_ids=tuple(item.rule_id for item in flipped),
        source_reaction_ids=_stable_unique(
            value for item in flipped for value in item.source_reaction_ids
        ),
        source_ec_numbers=_stable_unique(
            value for item in flipped for value in item.source_ec_numbers
        ),
        minimum_rule_specificity=min(specificities),
        worst_rule_score=_worst_score(scores, semantics[0]),
        score_semantics=semantics[0],
        balance_status=_status_union(
            (item.balance_status for item in flipped),
            incomplete_value="incomplete",
        ),
        cofactor_reconstruction_status=_status_union(
            (item.cofactor_reconstruction_status for item in flipped),
            incomplete_value="incomplete",
        ),
    )
    return step, flipped


def _kegg_step(step: PlanStep) -> HybridCandidateStep:
    consumed, produced = directional_stoichiometry(
        step.option.reaction,
        step.option.direction,
    )
    substrates = _normalize_stoichiometry(consumed)
    products = _normalize_stoichiometry(produced)
    reaction_id = step.option.reaction.reaction_id
    payload = {
        "step_source": "kegg_expansion",
        "reaction_option_ids": [reaction_id],
        "direction": step.option.direction,
        "substrate_stoichiometry": [list(item) for item in substrates],
        "product_stoichiometry": [list(item) for item in products],
    }
    return HybridCandidateStep(
        step_id=_step_id(payload),
        step_source="kegg_expansion",
        status="endogenous" if step.is_endogenous else "heterologous",
        orientation="biosynthetic",
        direction=step.option.direction,
        reaction_option_ids=(reaction_id,),
        reaction_smiles="",
        substrate_stoichiometry=substrates,
        product_stoichiometry=products,
        sink_anchor_kegg_ids=_stable_unique(step.expansion_anchor_compounds),
        expansion_depth=int(step.expansion_depth),
        is_endogenous=bool(step.is_endogenous),
        source_reaction_ids=_stable_unique(step.source_reaction_ids or (reaction_id,)),
        source_ec_numbers=_stable_unique(step.option.reaction.enzyme_ecs),
    )


def _attach_dependencies(
    kegg_pairs: Sequence[tuple[PlanStep, HybridCandidateStep]],
    retropath_steps: Sequence[HybridCandidateStep],
    sink_matches: Sequence[SinkMatch],
) -> Tuple[HybridCandidateStep, ...]:
    kegg_producer_ids: dict[str, set[str]] = defaultdict(set)
    retropath_producer_ids: dict[str, set[str]] = defaultdict(set)

    for plan_step, candidate_step in kegg_pairs:
        kegg_producer_ids[plan_step.option.produced_compound].add(
            candidate_step.step_id
        )
    for candidate_step in retropath_steps:
        for compound_id in candidate_step.product_compound_ids:
            retropath_producer_ids[compound_id].add(candidate_step.step_id)

    sink_ids = {item.representative_kegg_id for item in sink_matches}
    updated: dict[str, HybridCandidateStep] = {}
    for plan_step, candidate_step in kegg_pairs:
        dependencies = set()
        for precursor_id in plan_step.option.precursor_compounds:
            dependencies.update(kegg_producer_ids.get(precursor_id, set()))
        dependencies.discard(candidate_step.step_id)
        updated[candidate_step.step_id] = replace(
            candidate_step,
            depends_on_step_ids=tuple(sorted(dependencies)),
        )

    for candidate_step in retropath_steps:
        dependencies = set()
        for substrate_id in candidate_step.substrate_compound_ids:
            producers = set(retropath_producer_ids.get(substrate_id, set()))
            if substrate_id in sink_ids:
                producers.update(kegg_producer_ids.get(substrate_id, set()))
            dependencies.update(producers)
            if not producers and substrate_id not in sink_ids:
                raise ValueError(
                    "biosynthetic RP2 substrate has no producer or sink: "
                    f"{substrate_id}"
                )
        dependencies.discard(candidate_step.step_id)
        updated[candidate_step.step_id] = replace(
            candidate_step,
            depends_on_step_ids=tuple(sorted(dependencies)),
        )

    indegree = {step_id: 0 for step_id in updated}
    consumers: dict[str, set[str]] = defaultdict(set)
    for step_id, candidate_step in updated.items():
        for dependency in candidate_step.depends_on_step_ids:
            if dependency not in updated:
                raise ValueError(f"step dependency is not in candidate: {dependency}")
            indegree[step_id] += 1
            consumers[dependency].add(step_id)

    def ready_key(step_id: str) -> tuple[Any, ...]:
        item = updated[step_id]
        return (
            0 if item.step_source == "kegg_expansion" else 1,
            item.expansion_depth,
            item.reaction_option_ids,
            item.step_id,
        )

    ready = sorted(
        (step_id for step_id, value in indegree.items() if value == 0),
        key=ready_key,
    )
    ordered: list[HybridCandidateStep] = []
    while ready:
        step_id = ready.pop(0)
        ordered.append(updated[step_id])
        for consumer in sorted(consumers.get(step_id, set()), key=ready_key):
            indegree[consumer] -= 1
            if indegree[consumer] == 0:
                ready.append(consumer)
                ready.sort(key=ready_key)
    if len(ordered) != len(updated):
        raise ValueError(
            "combined KEGG/RetroPath candidate contains a dependency cycle"
        )
    return tuple(ordered)


def _candidate_id(
    path: RetrosyntheticPath,
    sink_matches: Sequence[SinkMatch],
    steps: Sequence[HybridCandidateStep],
) -> str:
    ordered_sink_matches = sorted(
        sink_matches,
        key=lambda item: (
            item.minimum_depth,
            item.representative_kegg_id,
            item.inchikey,
        ),
    )
    payload = {
        "entity_type": "hybrid_candidate_route",
        "schema_version": HYBRID_CANDIDATE_SCHEMA_VERSION,
        "target_compound_id": path.target_compound_id,
        "sink_matches": [
            {
                "representative_kegg_id": item.representative_kegg_id,
                "inchikey": item.inchikey,
                "minimum_depth": item.minimum_depth,
            }
            for item in ordered_sink_matches
        ],
        "steps": [
            {
                "step_id": item.step_id,
                "depends_on_step_ids": list(item.depends_on_step_ids),
            }
            for item in steps
        ],
    }
    return "RP2ROUTE:" + _canonical_sha256(payload)


def _candidate_rank(candidate: HybridCandidateRoute) -> tuple[Any, ...]:
    score_rank = (
        candidate.worst_rule_score
        if candidate.score_semantics == "lower_is_better"
        else -candidate.worst_rule_score
    )
    return (
        candidate.retropath_steps,
        candidate.maximum_sink_depth,
        -candidate.minimum_rule_specificity,
        score_rank,
        candidate.total_steps,
        candidate.kegg_prefix_steps,
        candidate.candidate_id,
    )


def _stable_rejections(
    values: Iterable[RetroPathMergeRejection],
) -> Tuple[RetroPathMergeRejection, ...]:
    unique = {
        (
            item.source_stage,
            item.reason_code,
            item.reason_detail,
            item.source_path_id,
            item.sink_kegg_ids,
            item.compound_id,
            item.transformation_id,
        ): item
        for item in values
    }
    return tuple(
        unique[key]
        for key in sorted(
            unique,
            key=lambda item: (
                item[0],
                item[1],
                item[3] or "",
                item[5] or "",
                item[6] or "",
                item[2],
            ),
        )
    )


def _validate_sink_matches(
    path: RetrosyntheticPath,
    expansion_bundle: ExpansionBundle,
) -> Optional[RetroPathMergeRejection]:
    for match in path.sink_matches:
        compound_id = match.representative_kegg_id
        if compound_id not in expansion_bundle.reachable_compounds:
            return RetroPathMergeRejection(
                source_stage="p5",
                source_path_id=path.path_id,
                reason_code="sink_not_in_expansion_bundle",
                reason_detail=f"sink {compound_id} is not cumulatively reachable",
                sink_kegg_ids=tuple(
                    item.representative_kegg_id for item in path.sink_matches
                ),
                compound_id=compound_id,
            )
        observed_depth = expansion_bundle.depth_by_compound.get(compound_id)
        if observed_depth != match.minimum_depth:
            return RetroPathMergeRejection(
                source_stage="p5",
                source_path_id=path.path_id,
                reason_code="sink_depth_mismatch",
                reason_detail=(
                    f"sink {compound_id} P4 depth={match.minimum_depth}, "
                    f"expansion depth={observed_depth!r}"
                ),
                sink_kegg_ids=tuple(
                    item.representative_kegg_id for item in path.sink_matches
                ),
                compound_id=compound_id,
            )
        if (
            match.minimum_depth == 0
            and compound_id not in expansion_bundle.base_compounds
        ):
            return RetroPathMergeRejection(
                source_stage="p5",
                source_path_id=path.path_id,
                reason_code="sink_depth_mismatch",
                reason_detail=f"depth-0 sink {compound_id} is absent from A0",
                sink_kegg_ids=tuple(
                    item.representative_kegg_id for item in path.sink_matches
                ),
                compound_id=compound_id,
            )
    return None


def _path_transformations(
    path: RetrosyntheticPath,
    transformations_by_id: Mapping[str, ParsedTransformation],
) -> Tuple[ParsedTransformation, ...]:
    selected = []
    for transformation_id in path.transformation_ids:
        transformation = transformations_by_id.get(transformation_id)
        if transformation is None:
            raise ValueError(
                f"P4 path references unknown transformation: {transformation_id}"
            )
        selected.append(transformation)
    return tuple(selected)


def _probe_witnesses(
    sink_matches: Sequence[SinkMatch],
    expansion_bundle: ExpansionBundle,
    client: KeggRestClient,
    ignored_common_compounds: set[str],
    max_witness_plans: int,
) -> Optional[str]:
    memo: dict[str, Tuple[Tuple[PlanStep, ...], ...]] = {}
    for match in sink_matches:
        compound_id = match.representative_kegg_id
        if compound_id in expansion_bundle.base_compounds:
            continue
        plans = build_frontier_bridge_plans(
            compound_id=compound_id,
            expansion_bundle=expansion_bundle,
            client=client,
            ignored_common_compounds=ignored_common_compounds,
            max_plans=max_witness_plans,
            memo=memo,
        )
        if not plans:
            return compound_id
    return None


def merge_retropath_candidates(
    enumeration_result: RetroPathEnumerationResult,
    expansion_bundle: ExpansionBundle,
    kegg_client: KeggRestClient,
    *,
    max_candidates: int = DEFAULT_MAX_CANDIDATES,
    max_witness_plans: int = DEFAULT_MAX_WITNESS_PLANS,
    max_total_steps: int = DEFAULT_MAX_TOTAL_STEPS,
    max_new_enzymes: int = DEFAULT_MAX_NEW_ENZYMES,
    ignored_common_compounds: Optional[set[str]] = None,
) -> RetroPathMergeResult:
    """Create Top-K biosynthetic DAG candidates from one P4 enumeration."""

    if not isinstance(enumeration_result, RetroPathEnumerationResult):
        raise ValueError("enumeration_result must be RetroPathEnumerationResult")
    if not isinstance(expansion_bundle, ExpansionBundle):
        raise ValueError("expansion_bundle must be an ExpansionBundle")
    candidate_limit = _positive_int(max_candidates, "max_candidates")
    witness_limit = _positive_int(max_witness_plans, "max_witness_plans")
    total_step_limit = _positive_int(max_total_steps, "max_total_steps")
    enzyme_limit = _nonnegative_int(max_new_enzymes, "max_new_enzymes")
    ignored = (
        set(IGNORED_COMMON_COMPOUNDS)
        if ignored_common_compounds is None
        else set(ignored_common_compounds)
    )

    rejections = [
        RetroPathMergeRejection(
            source_stage="p4",
            reason_code=item.reason_code,
            reason_detail=item.reason_detail,
            compound_id=item.compound_id,
            transformation_id=item.transformation_id,
        )
        for item in enumeration_result.rejections
    ]
    if enumeration_result.network.status in {"no_solution", "source_in_sink"}:
        return RetroPathMergeResult(
            candidates=tuple(),
            biosynthetic_reactions=tuple(),
            rejections=_stable_rejections(rejections),
            upstream_truncated=enumeration_result.truncated,
            truncated=False,
            max_candidates=candidate_limit,
            max_witness_plans=witness_limit,
            max_total_steps=total_step_limit,
            max_new_enzymes=enzyme_limit,
        )
    if enumeration_result.network.status != "succeeded":
        raise ValueError("P5 requires a terminal P4 network")

    transformations_by_id = {
        item.transformation_id: item
        for item in enumeration_result.network.transformations
    }
    generated: list[HybridCandidateRoute] = []
    flipped_by_id: dict[str, PredictedReaction] = {}

    for path in enumeration_result.paths:
        normalized_sink_matches = tuple(
            sorted(
                path.sink_matches,
                key=lambda item: (
                    item.minimum_depth,
                    item.representative_kegg_id,
                    item.inchikey,
                ),
            )
        )
        sink_ids = tuple(
            item.representative_kegg_id for item in normalized_sink_matches
        )
        sink_rejection = _validate_sink_matches(path, expansion_bundle)
        if sink_rejection is not None:
            rejections.append(sink_rejection)
            continue
        if path.reaction_count > total_step_limit:
            rejections.append(
                RetroPathMergeRejection(
                    source_stage="p5",
                    source_path_id=path.path_id,
                    reason_code="candidate_step_limit_exceeded",
                    reason_detail=(
                        f"RP2 steps {path.reaction_count} exceed "
                        f"max_total_steps={total_step_limit}"
                    ),
                    sink_kegg_ids=sink_ids,
                )
            )
            continue
        if path.reaction_count > enzyme_limit:
            rejections.append(
                RetroPathMergeRejection(
                    source_stage="p5",
                    source_path_id=path.path_id,
                    reason_code="candidate_enzyme_limit_exceeded",
                    reason_detail=(
                        f"RP2 steps {path.reaction_count} exceed "
                        f"max_new_enzymes={enzyme_limit}"
                    ),
                    sink_kegg_ids=sink_ids,
                )
            )
            continue

        try:
            selected_transformations = _path_transformations(
                path,
                transformations_by_id,
            )
            retropath_steps = []
            path_flipped: list[PredictedReaction] = []
            for transformation in selected_transformations:
                candidate_step, flipped = _retropath_step(transformation)
                retropath_steps.append(candidate_step)
                path_flipped.extend(flipped)
            retropath_substrates = {
                compound_id
                for item in retropath_steps
                for compound_id in item.substrate_compound_ids
            }
            missing_boundaries = set(sink_ids) - retropath_substrates
            if missing_boundaries:
                raise ValueError(
                    "P4 sink boundaries are not consumed by the biosynthetic "
                    "RP2 graph: " + ", ".join(sorted(missing_boundaries))
                )
            retropath_products = {
                compound_id
                for item in retropath_steps
                for compound_id in item.product_compound_ids
            }
            if path.target_compound_id not in retropath_products:
                raise ValueError(
                    "biosynthetic RP2 graph does not produce the requested target"
                )
            missing_witness = _probe_witnesses(
                normalized_sink_matches,
                expansion_bundle,
                kegg_client,
                ignored,
                witness_limit,
            )
            if missing_witness is not None:
                rejections.append(
                    RetroPathMergeRejection(
                        source_stage="p5",
                        source_path_id=path.path_id,
                        reason_code="expansion_witness_missing",
                        reason_detail=(
                            f"no complete A0 witness for sink {missing_witness}"
                        ),
                        sink_kegg_ids=sink_ids,
                        compound_id=missing_witness,
                    )
                )
                continue

            prefix_step_budget = total_step_limit - path.reaction_count
            prefix_enzyme_budget = enzyme_limit - path.reaction_count
            prefix_solutions = materialize_frontier_solution(
                solution=Solution(steps=tuple()),
                explicit_frontier_anchors=sink_ids,
                expansion_bundle=expansion_bundle,
                client=kegg_client,
                ignored_common_compounds=ignored,
                max_plans=witness_limit,
                max_total_steps=prefix_step_budget,
                max_new_enzymes=prefix_enzyme_budget,
            )
            if not prefix_solutions:
                rejections.append(
                    RetroPathMergeRejection(
                        source_stage="p5",
                        source_path_id=path.path_id,
                        reason_code="candidate_limit_exceeded",
                        reason_detail=(
                            "complete witness exists but no merged route satisfies "
                            "the remaining step/enzyme limits"
                        ),
                        sink_kegg_ids=sink_ids,
                    )
                )
                continue

            for prefix_solution in prefix_solutions:
                kegg_pairs = tuple(
                    (plan_step, _kegg_step(plan_step))
                    for plan_step in prefix_solution.steps
                )
                ordered_steps = _attach_dependencies(
                    kegg_pairs,
                    tuple(retropath_steps),
                    normalized_sink_matches,
                )
                if len(ordered_steps) > total_step_limit:
                    continue
                new_enzymes = sum(
                    item.status in {"heterologous", "predicted"}
                    for item in ordered_steps
                )
                if new_enzymes > enzyme_limit:
                    continue
                candidate = HybridCandidateRoute(
                    candidate_id=_candidate_id(
                        path,
                        normalized_sink_matches,
                        ordered_steps,
                    ),
                    source_retrosynthetic_path_id=path.path_id,
                    target_compound_id=path.target_compound_id,
                    sink_matches=normalized_sink_matches,
                    steps=ordered_steps,
                    minimum_rule_specificity=path.minimum_rule_specificity,
                    worst_rule_score=path.worst_rule_score,
                    score_semantics=path.score_semantics,
                    contains_auxiliary_fragments=path.contains_auxiliary_fragments,
                )
                generated.append(candidate)
                for reaction in path_flipped:
                    flipped_by_id[reaction.reaction_id] = reaction
        except ValueError as exc:
            rejections.append(
                RetroPathMergeRejection(
                    source_stage="p5",
                    source_path_id=path.path_id,
                    reason_code=(
                        "candidate_dependency_cycle"
                        if "cycle" in str(exc).lower()
                        else "candidate_merge_invalid"
                    ),
                    reason_detail=str(exc),
                    sink_kegg_ids=sink_ids,
                )
            )

    unique_candidates = {item.candidate_id: item for item in generated}
    ordered_candidates = sorted(unique_candidates.values(), key=_candidate_rank)
    truncated = len(ordered_candidates) > candidate_limit
    pruned_count = max(0, len(ordered_candidates) - candidate_limit)
    kept = tuple(ordered_candidates[:candidate_limit])
    if pruned_count:
        rejections.append(
            RetroPathMergeRejection(
                source_stage="p5",
                reason_code="top_k_pruned",
                reason_detail=(
                    f"{pruned_count} valid candidates omitted by "
                    f"max_candidates={candidate_limit}"
                ),
            )
        )
    retained_reaction_ids = {
        reaction_id
        for candidate in kept
        for reaction_id in candidate.retropath_reaction_option_ids
    }
    return RetroPathMergeResult(
        candidates=kept,
        biosynthetic_reactions=tuple(
            flipped_by_id[key]
            for key in sorted(flipped_by_id)
            if key in retained_reaction_ids
        ),
        rejections=_stable_rejections(rejections),
        upstream_truncated=enumeration_result.truncated,
        truncated=truncated,
        max_candidates=candidate_limit,
        max_witness_plans=witness_limit,
        max_total_steps=total_step_limit,
        max_new_enzymes=enzyme_limit,
    )


__all__ = [
    "DEFAULT_MAX_CANDIDATES",
    "DEFAULT_MAX_NEW_ENZYMES",
    "DEFAULT_MAX_TOTAL_STEPS",
    "DEFAULT_MAX_WITNESS_PLANS",
    "HYBRID_CANDIDATE_SCHEMA_VERSION",
    "HybridCandidateRoute",
    "HybridCandidateStep",
    "RetroPathMergeRejection",
    "RetroPathMergeResult",
    "flip_predicted_reaction",
    "merge_retropath_candidates",
]
