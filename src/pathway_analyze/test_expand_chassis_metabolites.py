from __future__ import annotations

import shutil
import unittest
import uuid
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

from src.pathway_analyze.expand_chassis_metabolites import (
    ensure_expansion_depth,
    load_expansion_bundle,
)
from src.pathway_analyze.gem_validation import validation_depth_output_dir
from src.pathway_analyze.kegg_gap_analyze import (
    KeggRestClient,
    ReactionRecord,
    build_frontier_bridge_plans,
    gap_depth_output_dir,
    materialize_frontier_solution,
    search_gap_solutions_once,
    Solution,
)


def reaction(
    reaction_id: str,
    left: tuple[str, ...],
    right: tuple[str, ...],
    *,
    comment: str = "",
) -> ReactionRecord:
    return ReactionRecord(
        reaction_id=reaction_id,
        name=reaction_id,
        comment=comment,
        equation=(" + ".join(left) + " <=> " + " + ".join(right)),
        annotation_text="",
        left_stoichiometry=tuple((item, 1.0) for item in left),
        right_stoichiometry=tuple((item, 1.0) for item in right),
        enzyme_ecs=("1.1.1.1",),
        ko_ids=tuple(),
        pathway_ids=tuple(),
        module_ids=tuple(),
    )


class FakeKeggClient:
    def __init__(self, reactions: dict[str, ReactionRecord]) -> None:
        self.reactions = reactions
        mutable: dict[str, list[str]] = {}
        for reaction_id, record in reactions.items():
            compounds = {
                item for item, _ in (*record.left_stoichiometry, *record.right_stoichiometry)
            }
            for compound_id in compounds:
                mutable.setdefault(compound_id, []).append(reaction_id)
        self.index = {
            compound_id: tuple(sorted(reaction_ids))
            for compound_id, reaction_ids in mutable.items()
        }

    def get_compound_reaction_index(self):
        return self.index

    def prefetch_reactions(self, reaction_ids):
        return None

    def try_get_reaction(self, reaction_id):
        return self.reactions.get(reaction_id)

    def get_compound_record(self, compound_id):
        return SimpleNamespace(pathway_ids=tuple())


class ExpansionAlgorithmTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = (
            Path.cwd()
            / "outputs"
            / ".expand_test_runtime"
            / uuid.uuid4().hex
        )
        self.root.mkdir(parents=True, exist_ok=False)
        self.base_path = self.root / "producible_kegg_compounds.csv"
        pd.DataFrame(
            [
                {
                    "source": "producible",
                    "met_id": "a_c",
                    "met_name": "A",
                    "compartment": "c",
                    "kegg_id": "C10000",
                }
            ]
        ).to_csv(self.base_path, index=False)
        self.reactions = {
            "R10000": reaction("R10000", ("C10000",), ("C10001",)),
            "R10001": reaction("R10001", ("C10001",), ("C10002",)),
            # 同一 F1 产物的第二条证据，用于验证多 witness 保留。
            "R10003": reaction("R10003", ("C10000",), ("C10001",)),
            # 缺少 C19999，因此不能在任何已有层触发。
            "R10004": reaction(
                "R10004",
                ("C10000", "C19999"),
                ("C10004",),
            ),
            # 泛化载体参与的反应必须被严格策略拒绝。
            "R10005": reaction("R10005", ("C10000", "C00028"), ("C10005",)),
            # 多步汇总条目不能冒充一步扩展。
            "R10006": reaction(
                "R10006",
                ("C10000",),
                ("C10006",),
                comment="two-step reaction",
            ),
        }
        self.client = FakeKeggClient(self.reactions)

    def tearDown(self) -> None:
        runtime_parent = self.root.parent
        shutil.rmtree(self.root, ignore_errors=True)
        try:
            runtime_parent.rmdir()
        except OSError:
            pass

    def _expand(self, depth: int):
        with patch(
            "src.pathway_analyze.expand_chassis_metabolites.load_endogenous_direction_index",
            return_value={},
        ):
            return ensure_expansion_depth(
                base_path=self.base_path,
                output_dir=self.root,
                model_path=self.root / "unused_model.json",
                cache_dir=self.root / "cache",
                requested_depth=depth,
                client=self.client,
            )

    def test_synchronous_layers_and_strict_substrate_policy(self) -> None:
        bundle = self._expand(1)
        self.assertIn("C10001", bundle.reachable_compounds)
        self.assertNotIn("C10002", bundle.reachable_compounds)
        self.assertNotIn("C10004", bundle.reachable_compounds)
        self.assertNotIn("C10005", bundle.reachable_compounds)
        self.assertNotIn("C10006", bundle.reachable_compounds)
        self.assertEqual(2, len(bundle.witnesses_by_product["C10001"]))

        bundle = self._expand(2)
        self.assertIn("C10002", bundle.reachable_compounds)
        self.assertEqual(2, bundle.depth_by_compound["C10002"])
        self.assertTrue((self.root / "chassis_frontier_depth_1.csv").is_file())
        self.assertTrue((self.root / "chassis_frontier_depth_2.csv").is_file())
        self.assertTrue(
            (self.root / "chassis_expanded_reachable_depth_2.csv").is_file()
        )

    def test_manifest_detects_changed_a0(self) -> None:
        self._expand(1)
        with self.base_path.open("a", encoding="utf-8") as handle:
            handle.write("producible,b_c,B,c,C10009\n")
        with self.assertRaisesRegex(ValueError, "missing or stale"):
            load_expansion_bundle(
                base_path=self.base_path,
                output_dir=self.root,
                depth=1,
            )

    def test_recursive_bridge_materialization_counts_real_steps(self) -> None:
        bundle = self._expand(2)
        plans = build_frontier_bridge_plans(
            compound_id="C10002",
            expansion_bundle=bundle,
            client=self.client,
            ignored_common_compounds=set(),
            max_plans=3,
        )
        self.assertEqual(2, len(plans))
        self.assertTrue(all(len(plan) == 2 for plan in plans))
        self.assertTrue(
            all(
                {step.expansion_depth for step in plan} == {1, 2}
                for plan in plans
            )
        )

        solutions = materialize_frontier_solution(
            solution=Solution(steps=tuple()),
            explicit_frontier_anchors=("C10002",),
            expansion_bundle=bundle,
            client=self.client,
            ignored_common_compounds=set(),
            max_plans=3,
            max_total_steps=10,
            max_new_enzymes=10,
        )
        self.assertEqual(2, len(solutions))
        self.assertTrue(all(item.total_steps == 2 for item in solutions))
        self.assertTrue(all(item.heterologous_steps == 2 for item in solutions))

    def test_global_compound_reaction_index_parser(self) -> None:
        client = KeggRestClient(cache_dir=None, use_shared_cache=False)
        client._fetch_text = lambda url, cache_key: (
            "cpd:C10000\trn:R10000\n"
            "cpd:C10000\trn:R10003\n"
            "cpd:C10001\trn:R10001\n"
        )
        index = client.get_compound_reaction_index()
        self.assertEqual(("R10000", "R10003"), index["C10000"])
        self.assertEqual(("R10001",), index["C10001"])

    def test_gap_output_directory_is_partitioned_by_depth(self) -> None:
        gap_root = self.root / "kegg_gap_C10000"
        self.assertEqual(gap_root / "depth0", gap_depth_output_dir(gap_root, 0))
        self.assertEqual(gap_root / "depth2", gap_depth_output_dir(gap_root, 2))
        self.assertEqual(
            gap_root / "depth2" / "gem_validation",
            validation_depth_output_dir(gap_root, 2),
        )
        with self.assertRaisesRegex(ValueError, "greater than or equal to 0"):
            gap_depth_output_dir(gap_root, -1)

    def test_prefetch_reuses_individual_reaction_files_across_batches(self) -> None:
        cache_dir = self.root / "kegg_cache"

        def reaction_entry(reaction_id: str) -> str:
            return (
                f"ENTRY       {reaction_id}                      Reaction\n"
                f"NAME        {reaction_id}\n"
                "EQUATION     C10000 <=> C10001\n"
                "///\n"
            )

        first_client = KeggRestClient(
            cache_dir=cache_dir,
            use_shared_cache=False,
        )
        fetch_disk_cache_flags: list[bool] = []

        def fake_batch_fetch(
            url: str,
            cache_key: str,
            *,
            use_disk_cache: bool = True,
        ) -> str:
            fetch_disk_cache_flags.append(use_disk_cache)
            return reaction_entry("R10000") + reaction_entry("R10001")

        first_client._fetch_text = fake_batch_fetch
        first_client.prefetch_reactions(
            ["R10000", "R10001"],
            batch_size=2,
        )

        self.assertEqual([False], fetch_disk_cache_flags)
        self.assertTrue((cache_dir / "reaction" / "R10000.txt").is_file())
        self.assertTrue((cache_dir / "reaction" / "R10001.txt").is_file())
        self.assertFalse((cache_dir / "reaction_batch").exists())

        second_client = KeggRestClient(
            cache_dir=cache_dir,
            use_shared_cache=False,
        )

        def unexpected_fetch(url: str, cache_key: str) -> str:
            raise AssertionError(
                f"disk-cached reactions must not be fetched again: {cache_key}"
            )

        second_client._fetch_text = unexpected_fetch
        second_client.prefetch_reactions(
            ["R10001", "R10000"],
            batch_size=1,
        )

        self.assertEqual(
            {"R10000", "R10001"},
            set(second_client._reaction_cache),
        )

    def test_prefetch_downloads_only_r1_minus_r2(self) -> None:
        cache_dir = self.root / "set_difference_kegg_cache"
        reaction_dir = cache_dir / "reaction"
        reaction_dir.mkdir(parents=True)

        def reaction_entry(reaction_id: str) -> str:
            return (
                f"ENTRY       {reaction_id}                      Reaction\n"
                "EQUATION     C10000 <=> C10001\n"
                "///\n"
            )

        for reaction_id in ("R10000", "R10001"):
            (reaction_dir / f"{reaction_id}.txt").write_text(
                reaction_entry(reaction_id),
                encoding="utf-8",
            )

        client = KeggRestClient(cache_dir=cache_dir, use_shared_cache=False)
        requested_batches: list[str] = []

        def fetch_missing(
            url: str,
            cache_key: str,
            *,
            use_disk_cache: bool = True,
        ) -> str:
            requested_batches.append(url)
            self.assertFalse(use_disk_cache)
            return reaction_entry("R10002") + reaction_entry("R10003")

        client._fetch_text = fetch_missing
        client.prefetch_reactions(
            ["R10000", "R10001", "R10002", "R10003"],
            batch_size=10,
        )

        self.assertEqual(1, len(requested_batches))
        self.assertNotIn("R10000", requested_batches[0])
        self.assertNotIn("R10001", requested_batches[0])
        self.assertIn("R10002", requested_batches[0])
        self.assertIn("R10003", requested_batches[0])
        self.assertEqual(
            {"R10000", "R10001", "R10002", "R10003"},
            set(client._reaction_cache),
        )
        self.assertFalse(
            any(key.startswith("reaction_batch:") for key in client._text_cache)
        )
        self.assertFalse((cache_dir / "reaction_batch").exists())

    def test_gap_depth_zero_and_frontier_target_completion(self) -> None:
        depth_zero = search_gap_solutions_once(
            target_compound="C10000",
            reachable_compounds={"C10000"},
            endogenous_reactions=set(),
            client=self.client,
            max_total_steps=10,
            max_new_enzymes=10,
            max_solutions=5,
            max_reactions_per_compound=5,
            max_major_precursors=4,
            max_routes_per_state=3,
            ignored_common_compounds=set(),
            allowed_reaction_ids=None,
            electron_avoidance_mode="off",
        )
        self.assertEqual(1, len(depth_zero.solutions))
        self.assertEqual(0, depth_zero.solutions[0].total_steps)

        bundle = self._expand(2)
        frontier_target = search_gap_solutions_once(
            target_compound="C10002",
            reachable_compounds=set(bundle.reachable_compounds),
            base_reachable_compounds=set(bundle.base_compounds),
            endogenous_reactions=set(),
            client=self.client,
            max_total_steps=10,
            max_new_enzymes=10,
            max_solutions=5,
            max_reactions_per_compound=5,
            max_major_precursors=4,
            max_routes_per_state=3,
            ignored_common_compounds=set(),
            allowed_reaction_ids=None,
            electron_avoidance_mode="off",
            expansion_bundle=bundle,
        )
        self.assertEqual(2, len(frontier_target.solutions))
        self.assertTrue(
            all(solution.total_steps == 2 for solution in frontier_target.solutions)
        )


if __name__ == "__main__":
    unittest.main()
