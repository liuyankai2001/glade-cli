"""Enumerate complete sink-closed routes from a parsed RetroPath network."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Optional, Tuple

from src.pathway_analyze.retropath_client import RetroPathClientRun
from src.pathway_analyze.retropath_input import RetroPathInputBundle
from src.pathway_analyze.retropath_parser import (
    ParsedRetroPathNetwork,
    ParsedTransformation,
    SinkMatch,
    parse_retropath_network,
)

RETROPATH_PATH_SCHEMA_VERSION = 1
DEFAULT_MAX_ROUTES = 1000
DEFAULT_MAX_SEARCH_STATES = 100_000


@dataclass(frozen=True)
class RouteRejection:
    """One stable reason why a network branch could not form a complete route."""

    reason_code: str
    reason_detail: str
    compound_id: Optional[str] = None
    transformation_id: Optional[str] = None

    def to_dict(self) -> dict[str, Optional[str]]:
        return {
            "reason_code": self.reason_code,
            "reason_detail": self.reason_detail,
            "compound_id": self.compound_id,
            "transformation_id": self.transformation_id,
        }


@dataclass(frozen=True)
class RetrosyntheticPath:
    """A complete Target-to-sink retrosynthetic reaction subgraph."""

    path_id: str
    target_compound_id: str
    transformation_ids: Tuple[str, ...]
    reaction_ids_by_transformation: Tuple[Tuple[str, Tuple[str, ...]], ...]
    sink_matches: Tuple[SinkMatch, ...]
    reaction_count: int
    maximum_branch_depth: int
    minimum_rule_specificity: int
    worst_rule_score: float
    score_semantics: str
    contains_auxiliary_fragments: bool
    review_required: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "path_id": self.path_id,
            "target_compound_id": self.target_compound_id,
            "orientation": "retrosynthetic",
            "transformation_ids": list(self.transformation_ids),
            "reaction_ids_by_transformation": [
                {
                    "transformation_id": transformation_id,
                    "reaction_ids": list(reaction_ids),
                }
                for transformation_id, reaction_ids in (
                    self.reaction_ids_by_transformation
                )
            ],
            "sink_matches": [item.to_dict() for item in self.sink_matches],
            "reaction_count": self.reaction_count,
            "maximum_branch_depth": self.maximum_branch_depth,
            "minimum_rule_specificity": self.minimum_rule_specificity,
            "worst_rule_score": self.worst_rule_score,
            "score_semantics": self.score_semantics,
            "contains_auxiliary_fragments": self.contains_auxiliary_fragments,
            "review_required": self.review_required,
        }


@dataclass(frozen=True)
class RetroPathEnumerationResult:
    """P4 network plus accepted routes and complete rejection audit."""

    network: ParsedRetroPathNetwork
    paths: Tuple[RetrosyntheticPath, ...]
    rejections: Tuple[RouteRejection, ...]
    explored_states: int
    max_routes: int
    max_search_states: int
    truncated: bool

    @property
    def complete_path_count(self) -> int:
        return len(self.paths)

    def to_dict(self) -> dict[str, Any]:
        return {
            "network": self.network.to_dict(),
            "paths": [item.to_dict() for item in self.paths],
            "rejections": [item.to_dict() for item in self.rejections],
            "explored_states": self.explored_states,
            "max_routes": self.max_routes,
            "max_search_states": self.max_search_states,
            "truncated": self.truncated,
            "complete_path_count": self.complete_path_count,
        }


@dataclass(frozen=True)
class _PartialRoute:
    choices: Tuple[Tuple[str, str], ...]
    sink_compound_ids: Tuple[str, ...]


def _positive_limit(value: Any, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{field_name} must be a positive integer")
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


def _stable_rejections(values: Iterable[RouteRejection]) -> Tuple[RouteRejection, ...]:
    unique = {
        (
            item.reason_code,
            item.reason_detail,
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
                item[2] or "",
                item[3] or "",
                item[1],
            ),
        )
    )


def _merge_partial(
    left: _PartialRoute, right: _PartialRoute
) -> Optional[_PartialRoute]:
    choices = dict(left.choices)
    for compound_id, transformation_id in right.choices:
        existing = choices.get(compound_id)
        if existing is not None and existing != transformation_id:
            return None
        choices[compound_id] = transformation_id
    return _PartialRoute(
        choices=tuple(sorted(choices.items())),
        sink_compound_ids=tuple(
            sorted(set(left.sink_compound_ids) | set(right.sink_compound_ids))
        ),
    )


def _partial_signature(partial: _PartialRoute) -> tuple[Any, ...]:
    return (len(partial.choices), partial.choices, partial.sink_compound_ids)


def _deduplicate_partials(
    values: Iterable[_PartialRoute],
    limit: int,
    rank: Callable[[_PartialRoute], tuple[Any, ...]] = _partial_signature,
) -> tuple[list[_PartialRoute], bool]:
    unique = {(item.choices, item.sink_compound_ids): item for item in values}
    ordered = sorted(unique.values(), key=rank)
    return ordered[:limit], len(ordered) > limit


def _path_from_partial(
    network: ParsedRetroPathNetwork,
    partial: _PartialRoute,
    transformations_by_id: Mapping[str, ParsedTransformation],
    sink_by_compound_id: Mapping[str, SinkMatch],
) -> RetrosyntheticPath:
    selected = [
        transformations_by_id[transformation_id]
        for _, transformation_id in partial.choices
    ]
    selected.sort(
        key=lambda item: (
            item.iteration,
            item.substrate_compound_id,
            item.transformation_id,
        )
    )
    sink_matches = tuple(
        sorted(
            (sink_by_compound_id[item] for item in partial.sink_compound_ids),
            key=lambda item: (
                item.minimum_depth,
                item.representative_kegg_id,
                item.inchikey,
            ),
        )
    )
    steps = [
        {
            "iteration": item.iteration,
            "substrate_compound_id": item.substrate_compound_id,
            "product_compound_ids": list(item.product_compound_ids),
            "reaction_options": sorted(
                reaction.reaction_id for reaction in item.reaction_variants
            ),
        }
        for item in selected
    ]
    path_id = "RP2PATH:" + _canonical_sha256(
        {
            "entity_type": "retrosynthetic_path",
            "schema_version": RETROPATH_PATH_SCHEMA_VERSION,
            "target_compound_id": network.target_compound_id,
            "steps": steps,
            "sink_inchikeys": [item.inchikey for item in sink_matches],
        }
    )
    score_semantics = selected[0].score_semantics
    scores = [item.score_raw for item in selected]
    worst_score = max(scores) if score_semantics == "lower_is_better" else min(scores)
    return RetrosyntheticPath(
        path_id=path_id,
        target_compound_id=network.target_compound_id,
        transformation_ids=tuple(item.transformation_id for item in selected),
        reaction_ids_by_transformation=tuple(
            (
                item.transformation_id,
                tuple(
                    sorted(reaction.reaction_id for reaction in item.reaction_variants)
                ),
            )
            for item in selected
        ),
        sink_matches=sink_matches,
        reaction_count=len(selected),
        maximum_branch_depth=max(item.iteration for item in selected) + 1,
        minimum_rule_specificity=min(
            item.minimum_rule_specificity for item in selected
        ),
        worst_rule_score=worst_score,
        score_semantics=score_semantics,
        contains_auxiliary_fragments=any(item.auxiliary_fragments for item in selected),
    )


def _path_rank(path: RetrosyntheticPath) -> tuple[Any, ...]:
    score_rank = (
        path.worst_rule_score
        if path.score_semantics == "lower_is_better"
        else -path.worst_rule_score
    )
    return (
        path.reaction_count,
        path.maximum_branch_depth,
        -path.minimum_rule_specificity,
        score_rank,
        path.path_id,
    )


def enumerate_sink_routes(
    network: ParsedRetroPathNetwork,
    *,
    max_routes: int = DEFAULT_MAX_ROUTES,
    max_search_states: int = DEFAULT_MAX_SEARCH_STATES,
) -> RetroPathEnumerationResult:
    """Enumerate AND/OR routes whose every structural leaf is a P2 sink."""

    if not isinstance(network, ParsedRetroPathNetwork):
        raise ValueError("network must be a ParsedRetroPathNetwork")
    route_limit = _positive_limit(max_routes, "max_routes")
    state_limit = _positive_limit(max_search_states, "max_search_states")
    inherited_rejections = [
        RouteRejection(
            item.reason_code,
            item.reason_detail,
            compound_id=item.compound_id,
            transformation_id=item.transformation_id,
        )
        for item in network.rejections
    ]
    if network.status in {"no_solution", "source_in_sink"}:
        return RetroPathEnumerationResult(
            network=network,
            paths=tuple(),
            rejections=_stable_rejections(inherited_rejections),
            explored_states=0,
            max_routes=route_limit,
            max_search_states=state_limit,
            truncated=False,
        )
    if network.status != "succeeded":
        raise ValueError(
            "only succeeded, no_solution, or source_in_sink networks can be enumerated"
        )

    transformations_by_id = {
        item.transformation_id: item for item in network.transformations
    }
    outgoing: dict[str, list[ParsedTransformation]] = {}
    for transformation in network.transformations:
        outgoing.setdefault(transformation.substrate_compound_id, []).append(
            transformation
        )
    for values in outgoing.values():
        values.sort(
            key=lambda item: (
                item.iteration,
                -item.minimum_rule_specificity,
                (
                    item.score_raw
                    if item.score_semantics == "lower_is_better"
                    else -item.score_raw
                ),
                item.transformation_id,
            )
        )
    sink_by_compound_id = {item.compound_id: item for item in network.sink_matches}

    def partial_rank(partial: _PartialRoute) -> tuple[Any, ...]:
        selected = [
            transformations_by_id[transformation_id]
            for _, transformation_id in partial.choices
        ]
        if not selected:
            return (0, 0, 0, 0.0, partial.choices, partial.sink_compound_ids)
        score_semantics = selected[0].score_semantics
        scores = [item.score_raw for item in selected]
        worst_score = (
            max(scores) if score_semantics == "lower_is_better" else -min(scores)
        )
        return (
            len(selected),
            max(item.iteration for item in selected) + 1,
            -min(item.minimum_rule_specificity for item in selected),
            worst_score,
            partial.choices,
            partial.sink_compound_ids,
        )

    rejections = inherited_rejections
    explored_states = 0
    truncated = False
    state_limit_reached = False

    def tick() -> bool:
        nonlocal explored_states, state_limit_reached, truncated
        if explored_states >= state_limit:
            state_limit_reached = True
            truncated = True
            return False
        explored_states += 1
        return True

    def reject(
        code: str,
        detail: str,
        *,
        compound_id: Optional[str] = None,
        transformation_id: Optional[str] = None,
    ) -> None:
        rejections.append(
            RouteRejection(
                code,
                detail,
                compound_id=compound_id,
                transformation_id=transformation_id,
            )
        )

    def solve(
        compound_id: str,
        expected_iteration: int,
        ancestors: Tuple[str, ...],
    ) -> list[_PartialRoute]:
        nonlocal truncated
        if not tick():
            return []
        if compound_id in sink_by_compound_id:
            return [
                _PartialRoute(
                    choices=tuple(),
                    sink_compound_ids=(compound_id,),
                )
            ]
        if compound_id in ancestors:
            reject(
                "cycle_detected",
                "compound re-enters its own retrosynthetic ancestry",
                compound_id=compound_id,
            )
            return []
        if expected_iteration >= network.max_steps:
            reject(
                "depth_exceeded",
                f"branch exceeds max_steps={network.max_steps}",
                compound_id=compound_id,
            )
            return []

        candidates = outgoing.get(compound_id, [])
        matching = [item for item in candidates if item.iteration == expected_iteration]
        if not matching:
            if candidates:
                observed = ",".join(str(item.iteration) for item in candidates)
                reject(
                    "reaction_direction_invalid",
                    (
                        f"expected iteration {expected_iteration}, observed "
                        f"transformations at {observed}"
                    ),
                    compound_id=compound_id,
                )
            else:
                reject(
                    "unresolved_non_sink_leaf",
                    "structural leaf has no outgoing transformation and is not "
                    "a P2 sink",
                    compound_id=compound_id,
                )
            return []

        solutions: list[_PartialRoute] = []
        next_ancestors = (*ancestors, compound_id)
        for transformation in matching:
            products = transformation.unique_product_compound_ids
            if compound_id in products:
                reject(
                    "cycle_detected",
                    "transformation regenerates its own substrate",
                    compound_id=compound_id,
                    transformation_id=transformation.transformation_id,
                )
                continue
            combinations = [
                _PartialRoute(
                    choices=((compound_id, transformation.transformation_id),),
                    sink_compound_ids=tuple(),
                )
            ]
            branch_failed = False
            for product_id in products:
                child_solutions = solve(
                    product_id,
                    expected_iteration + 1,
                    next_ancestors,
                )
                if not child_solutions:
                    branch_failed = True
                    break
                merged_values = []
                merge_limit_reached = False
                for existing in combinations:
                    for child in child_solutions:
                        if not tick():
                            merge_limit_reached = True
                            break
                        merged = _merge_partial(existing, child)
                        if merged is not None:
                            merged_values.append(merged)
                    if merge_limit_reached:
                        break
                combinations, was_trimmed = _deduplicate_partials(
                    merged_values,
                    route_limit,
                    partial_rank,
                )
                truncated = truncated or was_trimmed
                if not combinations:
                    branch_failed = True
                    break
            if not branch_failed:
                solutions.extend(combinations)
                solutions, was_trimmed = _deduplicate_partials(
                    solutions,
                    route_limit,
                    partial_rank,
                )
                truncated = truncated or was_trimmed
        return solutions

    partials = solve(network.target_compound_id, 0, tuple())

    paths_by_id: dict[str, RetrosyntheticPath] = {}
    for partial in partials:
        if not partial.choices:
            continue
        path = _path_from_partial(
            network,
            partial,
            transformations_by_id,
            sink_by_compound_id,
        )
        paths_by_id[path.path_id] = path
    ordered_paths = sorted(paths_by_id.values(), key=_path_rank)
    if len(ordered_paths) > route_limit:
        ordered_paths = ordered_paths[:route_limit]
        truncated = True
    if truncated and not any(
        item.reason_code == "enumeration_limit_reached" for item in rejections
    ):
        detail = (
            f"search stopped after max_search_states={state_limit}"
            if state_limit_reached
            else f"route combinations exceeded max_routes={route_limit}"
        )
        reject(
            "enumeration_limit_reached",
            detail,
            compound_id=network.target_compound_id,
        )
    return RetroPathEnumerationResult(
        network=network,
        paths=tuple(ordered_paths),
        rejections=_stable_rejections(rejections),
        explored_states=explored_states,
        max_routes=route_limit,
        max_search_states=state_limit,
        truncated=truncated,
    )


def parse_and_enumerate_retropath(
    client_run: RetroPathClientRun,
    input_bundle: RetroPathInputBundle,
    rules_path: str | Path,
    *,
    max_routes: int = DEFAULT_MAX_ROUTES,
    max_search_states: int = DEFAULT_MAX_SEARCH_STATES,
) -> RetroPathEnumerationResult:
    """Run the complete P4 parse-and-enumerate pipeline."""

    network = parse_retropath_network(client_run, input_bundle, rules_path)
    return enumerate_sink_routes(
        network,
        max_routes=max_routes,
        max_search_states=max_search_states,
    )


__all__ = [
    "DEFAULT_MAX_ROUTES",
    "DEFAULT_MAX_SEARCH_STATES",
    "RETROPATH_PATH_SCHEMA_VERSION",
    "RetroPathEnumerationResult",
    "RetrosyntheticPath",
    "RouteRejection",
    "enumerate_sink_routes",
    "parse_and_enumerate_retropath",
]
