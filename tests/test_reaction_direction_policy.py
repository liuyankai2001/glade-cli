from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from src.config.run_config import RunConfig
from src.pathway_analyze import expand_chassis_metabolites as expansion
from src.pathway_analyze import kegg_gap_analyze as gap
from src.pathway_analyze.reaction_direction_policy import (
    reaction_direction_allowed,
    validate_reaction_direction_mode,
)


def reaction(identifier: str, equation: str, *, name: str = "", comment: str = ""):
    left, right = gap.split_equation(equation)
    return gap.ReactionRecord(
        reaction_id=identifier,
        name=name or identifier,
        comment=comment,
        equation=equation,
        annotation_text="",
        left_stoichiometry=gap.parse_equation_side(left),
        right_stoichiometry=gap.parse_equation_side(right),
        enzyme_ecs=("1.1.1.1",),
        ko_ids=(),
        pathway_ids=(),
        module_ids=(),
    )


class LocalKeggClient:
    """Small offline KEGG network; unexpected retrieval never reaches the web."""

    def __init__(self, records, *, module_ids=()):
        self.records = {record.reaction_id: record for record in records}
        self.module_ids = module_ids

    def get_compound_record(self, identifier):
        return gap.CompoundRecord(identifier, (), self.module_ids)

    def get_compound_name(self, identifier):
        return identifier

    def get_reaction_ids_for_compound(self, identifier):
        return tuple(
            record.reaction_id
            for record in self.records.values()
            if identifier in gap.reaction_compounds(record)
        )

    def get_compound_reaction_index(self):
        return {
            identifier: self.get_reaction_ids_for_compound(identifier)
            for record in self.records.values()
            for identifier in gap.reaction_compounds(record)
        }

    def try_get_reaction(self, identifier):
        return self.records.get(identifier)

    def prefetch_reactions(self, identifiers):
        pass

    def prefetch_modules(self, identifiers):
        pass

    def get_module(self, identifier):
        forward_ids = tuple(
            record.reaction_id
            for record in self.records.values()
            if "C00246" in dict(record.right_stoichiometry)
        )
        return gap.ModuleRecord("test module", forward_ids, (), (), ())


def search(client, *, target="C00246", seeds=None, mode=3, **overrides):
    arguments = dict(
        target_compound=target,
        reachable_compounds=seeds or {"C00136", "C00332", "C00042", "C00033"},
        endogenous_reactions=set(),
        client=client,
        max_total_steps=3,
        max_new_enzymes=3,
        max_solutions=10,
        max_reactions_per_compound=15,
        max_major_precursors=4,
        max_routes_per_state=10,
        ignored_common_compounds=set(gap.IGNORED_COMMON_COMPOUNDS),
        module_filter_mode="off",
        unpre_re_reaction=mode,
    )
    arguments.update(overrides)
    return gap.run_search_gap_analysis(**arguments)


@pytest.mark.parametrize(
    "mode,expected",
    [
        (0, {("left_to_right", "favored"), ("left_to_right", "neutral"),
             ("left_to_right", "disfavored"), ("right_to_left", "favored"),
             ("right_to_left", "neutral"), ("right_to_left", "disfavored")}),
        (1, {("left_to_right", "favored"), ("left_to_right", "neutral"),
             ("right_to_left", "favored"), ("right_to_left", "neutral")}),
        (2, {("left_to_right", "favored"), ("right_to_left", "favored")}),
        (3, {("left_to_right", "favored")}),
    ],
)
def test_four_progressively_stricter_policies(mode, expected):
    actual = {
        (direction, label)
        for direction in ("left_to_right", "right_to_left")
        for label in ("favored", "neutral", "disfavored")
        if reaction_direction_allowed(direction, label, mode)
    }
    assert actual == expected


def test_allowed_sets_shrink_with_each_level():
    combinations = {
        (direction, label)
        for direction in ("left_to_right", "right_to_left")
        for label in ("favored", "neutral", "disfavored")
    }
    allowed = [
        {item for item in combinations if reaction_direction_allowed(*item, mode)}
        for mode in range(4)
    ]
    assert allowed[3] < allowed[2] < allowed[1] < allowed[0]


@pytest.mark.parametrize("value", [True, False, None, "1", 1.0, -1, 4, [], {}])
def test_invalid_policy_fails_before_file_access(value):
    with pytest.raises(ValueError, match="unpre_re_reaction"):
        validate_reaction_direction_mode(value)
    config = SimpleNamespace(target_name="C00246", unpre_re_reaction=value)
    with pytest.raises(ValueError, match="unpre_re_reaction"):
        gap.kegg_gap_analyze(config)
    with pytest.raises(ValueError, match="unpre_re_reaction"):
        expansion.expand_chassis_metabolites(config)


def test_default_is_three():
    assert RunConfig(target_name="C00246").unpre_re_reaction == 3


@pytest.fixture
def butyrate_client():
    # Equations from the archived butyrate KEGG case, with no network dependency.
    return LocalKeggClient([
        reaction("R01179", "C00136 + C00033 <=> C00246 + C00024"),
        reaction("R01176", "C00002 + C00246 + C00010 <=> C00020 + C00013 + C00136",
                 name="butanoate:CoA ligase (AMP-forming)"),
        reaction("R01365", "C00136 + C00164 <=> C00246 + C00332"),
        reaction("R01358", "C00332 + C00001 <=> C00164 + C00010"),
        reaction("R00410", "C00091 + C00164 <=> C00042 + C00332"),
        reaction("R01359", "C00332 + C00033 <=> C00164 + C00024"),
    ])


@pytest.mark.parametrize("mode,expected_count", [(0, 5), (1, 4), (2, 0), (3, 0)])
def test_butyrate_routes_and_persisted_policy(butyrate_client, tmp_path, mode, expected_count):
    result = search(butyrate_client, mode=mode)
    assert len(result.solutions) == expected_count
    reaction_ids = {
        step.option.reaction.reaction_id
        for solution in result.solutions
        for step in solution.steps
    }
    assert ("R01176" in reaction_ids) == (mode == 0)
    assert ("R00410" in reaction_ids) == (mode in {0, 1})
    gap.write_outputs("C00246", result, butyrate_client, tmp_path, {})
    saved = json.loads((tmp_path / "run_config.json").read_text(encoding="utf-8"))
    assert saved["unpre_re_reaction"] == mode
    assert saved["reaction_direction_policy_version"] == "kegg_direction_policy.v2"


def test_policy_filters_before_candidate_limit():
    client = LocalKeggClient([
        reaction("R10001", "C10000 <=> C10001"),
        reaction("R10002", "C10002 + C00002 <=> C10000 + C00008"),
    ])
    result = search(client, target="C10000", seeds={"C10001", "C10002"},
                    max_reactions_per_compound=1)
    assert len(result.solutions) == 1
    assert result.solutions[0].steps[0].option.reaction.reaction_id == "R10002"


def test_module_and_electron_fallback_do_not_relax_direction_policy():
    client = LocalKeggClient([
        reaction("R01176", "C00002 + C00246 + C00010 <=> C00020 + C00013 + C00136"),
        reaction("R10001", "C10000 <=> C00246"),
    ], module_ids=("M00001",))
    result = search(client, mode=1, module_filter_mode="prefer",
                    electron_avoidance_mode="strict_with_fallback")
    assert result.solutions == ()
    assert result.did_fallback
    assert result.electron_avoidance_fallback
    assert search(client, mode=0, module_filter_mode="prefer").solutions


@pytest.mark.parametrize("mode,expected", [(0, True), (1, True), (2, False), (3, False)])
def test_multistep_summary_cannot_hide_forbidden_reverse_step(mode, expected):
    summary = reaction("R19999", "C10002 <=> C10000",
                       comment="two-step reaction R10001 R10002")
    client = LocalKeggClient([
        summary,
        reaction("R10001", "C10001 <=> C10000"),
        reaction("R10002", "C10001 <=> C10002"),
    ])
    result = search(client, target="C10000", seeds={"C10002"}, mode=mode)
    assert bool(result.solutions) == expected
    if expected:
        assert any(
            step.option.reaction.reaction_id == "R10002"
            and step.option.direction == "right_to_left"
            for solution in result.solutions for step in solution.steps
        )
    option = gap.make_reaction_option("C10000", summary, "left_to_right", set())
    components = gap.decompose_multistep_option(option, client, set(), mode)
    assert bool(components) == expected
    assert all(item.reaction.reaction_id != "R19999" for item in components)


def test_normalization_cannot_replace_pair_with_forbidden_direction():
    first = reaction("R10001", "C10001 <=> C10000", comment="incomplete reaction")
    second = reaction("R10002", "C10002 <=> C10001", comment="incomplete reaction")
    replacement = reaction("R10003", "C10000 <=> C10002")
    client = LocalKeggClient([first, second, replacement])
    solution = gap.Solution(tuple(
        gap.make_plan_step(gap.make_reaction_option(target, record, "left_to_right", set()),
                           (), False)
        for target, record in [("C10000", first), ("C10001", second)]
    ))
    allowed = gap.normalize_solution_reactions(solution, client, {"C10002"}, set(), set(), None, 0)
    blocked = gap.normalize_solution_reactions(solution, client, {"C10002"}, set(), set(), None, 3)
    assert len(allowed.steps) == 1
    assert allowed.steps[0].option.direction == "right_to_left"
    assert len(blocked.steps) == 2


def build_expansion(tmp_path, *, mode=0):
    base = tmp_path / "base.csv"
    if not base.exists():
        base.write_text("kegg_id\nC10001\n", encoding="utf-8")
    client = LocalKeggClient([
        reaction("R10001", "C10000 <=> C10001"),
        reaction("R10002", "C10001 <=> C10002"),
    ])
    arguments = dict(base_path=base, output_dir=tmp_path / "expanded",
                     model_path=tmp_path / "model.json", cache_dir=tmp_path / "cache",
                     requested_depth=1, client=client, unpre_re_reaction=mode)
    with patch.object(expansion, "load_endogenous_direction_index", return_value={}):
        bundle = expansion.ensure_expansion_depth(**arguments)
    return bundle, arguments


def test_expansion_cache_requires_matching_policy_and_rebuilds(tmp_path):
    bundle, arguments = build_expansion(tmp_path, mode=0)
    assert bundle.reachable_compounds == {"C10000", "C10001", "C10002"}
    before = bundle.expanded_file.read_bytes()
    with pytest.raises(ValueError, match="unpre_re_reaction"):
        expansion.load_expansion_bundle(base_path=arguments["base_path"],
                                         output_dir=arguments["output_dir"], depth=1)
    assert bundle.expanded_file.read_bytes() == before
    bundle, arguments = build_expansion(tmp_path, mode=3)
    assert bundle.reachable_compounds == {"C10001"}
    assert bundle.manifest["unpre_re_reaction"] == 3
    arguments["client"] = Mock(side_effect=AssertionError("cache must be reused"))
    with patch.object(expansion, "load_endogenous_direction_index",
                      side_effect=AssertionError("cache must be reused")):
        assert expansion.ensure_expansion_depth(**arguments).manifest == bundle.manifest


def test_old_expansion_manifest_is_rejected(tmp_path):
    bundle, arguments = build_expansion(tmp_path)
    manifest = dict(bundle.manifest)
    manifest.pop("unpre_re_reaction")
    manifest_path = arguments["output_dir"] / expansion.MANIFEST_FILE_NAME
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="expand command"):
        expansion.load_expansion_bundle(base_path=arguments["base_path"],
                                         output_dir=arguments["output_dir"], depth=1,
                                         unpre_re_reaction=0)


@pytest.mark.parametrize("mode", [0, 3])
def test_previous_policy_version_is_rejected_and_rebuilt(tmp_path, mode):
    bundle, arguments = build_expansion(tmp_path, mode=mode)
    manifest_path = arguments["output_dir"] / expansion.MANIFEST_FILE_NAME
    legacy_manifest = {**bundle.manifest, "direction_policy": "kegg_direction_policy.v1"}
    manifest_path.write_text(json.dumps(legacy_manifest), encoding="utf-8")
    before = bundle.expanded_file.read_bytes()
    with pytest.raises(ValueError, match="expand command"):
        expansion.load_expansion_bundle(
            base_path=arguments["base_path"], output_dir=arguments["output_dir"],
            depth=1, unpre_re_reaction=mode,
        )
    assert bundle.expanded_file.read_bytes() == before
    with patch.object(expansion, "load_endogenous_direction_index", return_value={}):
        rebuilt = expansion.ensure_expansion_depth(**arguments)
    assert rebuilt.manifest["direction_policy"] == "kegg_direction_policy.v2"
    assert rebuilt.reachable_compounds == bundle.reachable_compounds


@pytest.mark.parametrize("mode,expected", [(0, True), (1, True), (2, False), (3, False)])
def test_restoring_expansion_routes_cannot_bypass_policy(tmp_path, mode, expected):
    bundle, arguments = build_expansion(tmp_path)
    restored = gap.materialize_frontier_solution(
        solution=gap.Solution(()), explicit_frontier_anchors=("C10000",),
        expansion_bundle=bundle, client=arguments["client"],
        ignored_common_compounds=set(), max_plans=2, max_total_steps=2,
        max_new_enzymes=2, unpre_re_reaction=mode,
    )
    assert bool(restored) == expected
    # Also cover the target-already-in-expanded-set shortcut in gap.
    result = search(arguments["client"], target="C10000", seeds=set(bundle.reachable_compounds),
                    base_reachable_compounds=set(bundle.base_compounds),
                    expansion_bundle=bundle, mode=mode)
    assert bool(result.solutions) == expected


@pytest.mark.parametrize("mode,expected", [(0, True), (1, False), (2, False), (3, False)])
def test_expansion_mode_zero_includes_previously_filtered_disfavored(butyrate_client, mode, expected):
    witnesses = expansion._reaction_witnesses(
        reaction=butyrate_client.records["R01176"], depth=1,
        available_compounds={"C00136", "C00020", "C00013"},
        previous_frontier={"C00136"}, endogenous_reactions=set(),
        endogenous_direction_index={}, unpre_re_reaction=mode,
    )
    assert any(item.product_compound == "C00246" for item in witnesses) == expected


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_favored_forward_expansion_is_retained(mode):
    record = reaction("R10001", "C10000 + C00002 <=> C10001 + C00008")
    witnesses = expansion._reaction_witnesses(
        reaction=record, depth=1, available_compounds={"C10000", "C00002"},
        previous_frontier={"C10000"}, endogenous_reactions=set(),
        endogenous_direction_index={}, unpre_re_reaction=mode,
    )
    assert any(item.product_compound == "C10001" for item in witnesses)


def test_output_guard_runs_before_removing_previous_results(butyrate_client, tmp_path):
    result = search(butyrate_client, mode=0)
    previous = tmp_path / "solutions.csv"
    previous.write_text("keep previous results", encoding="utf-8")
    with pytest.raises(ValueError, match="violates unpre_re_reaction"):
        gap.write_outputs("C00246", replace(result, unpre_re_reaction=3),
                          butyrate_client, tmp_path, {})
    assert previous.read_text(encoding="utf-8") == "keep previous results"


@pytest.mark.parametrize("mode,expected_count", [(0, 5), (3, 0)])
def test_gap_entrypoint_passes_and_records_configuration(butyrate_client, tmp_path, mode, expected_count):
    config = RunConfig(target_name="C00246")
    config.unpre_re_reaction = mode
    config.module_filter_mode = "off"
    config.electron_avoidance_mode = "off"
    config.model_path = tmp_path / "model.json"
    config.model_path.write_text("{}", encoding="utf-8")
    config.chassis_output_path = tmp_path / "chassis"
    config.chassis_output_path.mkdir()
    config.chassis_producible_csv = config.chassis_output_path / "base.csv"
    pd.DataFrame({"kegg_id": ["C00136", "C00332", "C00042", "C00033"]}).to_csv(
        config.chassis_producible_csv, index=False)
    config.gap_output_path = tmp_path / "gap"
    config.cache_dir = tmp_path / "cache"
    with patch.object(gap, "load_endogenous_direction_index", return_value={}), \
         patch.object(gap, "KeggRestClient", return_value=butyrate_client):
        result = gap.kegg_gap_analyze(config)
    assert result["solution_count"] == expected_count
    assert result["unpre_re_reaction"] == mode
    saved = json.loads((config.gap_output_path / "depth0" / "run_config.json").read_text(encoding="utf-8"))
    assert saved["unpre_re_reaction"] == mode


def test_forbidden_summary_does_not_consume_candidate_limit():
    client = LocalKeggClient([
        reaction("R09999", "C10002 + C00002 <=> C10000 + C00008",
                 comment="two-step reaction R10001 R10002"),
        reaction("R10001", "C10001 + C00002 <=> C10000 + C00008"),
        reaction("R10002", "C10001 <=> C10002"),
        reaction("R10003", "C10002 + C00002 <=> C10000 + C00008"),
    ])
    result = search(client, target="C10000", seeds={"C10002"},
                    max_reactions_per_compound=1)
    assert len(result.solutions) == 1
    assert result.solutions[0].steps[0].option.reaction.reaction_id == "R10003"


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_favored_forward_route_is_retained(mode):
    client = LocalKeggClient([
        reaction("R10001", "C10000 + C00002 <=> C10001 + C00008"),
    ])
    result = search(client, target="C10001", seeds={"C10000"}, mode=mode)
    assert len(result.solutions) == 1
    assert result.solutions[0].steps[0].option.screening.thermo_direction == "favored"


def test_multistep_summary_cannot_hide_disfavored_component():
    summary = reaction("R19999", "C10002 <=> C10000",
                       comment="two-step reaction R10001 R10002")
    client = LocalKeggClient([
        summary,
        reaction("R10001", "C10001 <=> C10000"),
        reaction("R10002", "C10001 + C00002 <=> C10002 + C00020 + C00013"),
    ])
    option = gap.make_reaction_option("C10000", summary, "left_to_right", set())
    ignored = set(gap.IGNORED_COMMON_COMPOUNDS)
    assert gap.decompose_multistep_option(option, client, ignored, 0)
    assert not gap.decompose_multistep_option(option, client, ignored, 1)


def test_expansion_rejects_witnesses_that_violate_manifest_policy(tmp_path):
    bundle, arguments = build_expansion(tmp_path, mode=0)
    # A forged policy value must not make a reverse witness reusable in mode 3.
    manifest = {**bundle.manifest, "unpre_re_reaction": 3}
    (arguments["output_dir"] / expansion.MANIFEST_FILE_NAME).write_text(
        json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="witnesses violate"):
        expansion.load_expansion_bundle(base_path=arguments["base_path"],
                                         output_dir=arguments["output_dir"], depth=1,
                                         unpre_re_reaction=3)


@pytest.mark.parametrize("mode,expected_count", [(0, 2), (2, 0), (3, 0)])
def test_expand_entrypoint_passes_configuration(tmp_path, mode, expected_count):
    base = tmp_path / "base.csv"
    base.write_text("kegg_id\nC10001\n", encoding="utf-8")
    model = tmp_path / "model.json"
    model.write_text("{}", encoding="utf-8")
    client = LocalKeggClient([
        reaction("R10001", "C10000 <=> C10001"),
        reaction("R10002", "C10001 <=> C10002"),
    ])
    config = SimpleNamespace(depth=1, chassis_producible_csv=base,
                             chassis_output_path=tmp_path / "expanded",
                             model_path=model, cache_dir=tmp_path / "cache",
                             unpre_re_reaction=mode)
    with patch.object(expansion, "load_endogenous_direction_index", return_value={}), \
         patch.object(expansion, "KeggRestClient", return_value=client):
        result = expansion.expand_chassis_metabolites(config)
    assert result["unpre_re_reaction"] == mode
    assert result["expanded_compound_count"] == expected_count
