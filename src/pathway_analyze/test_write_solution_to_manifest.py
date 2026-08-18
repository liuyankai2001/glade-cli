from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from src.pathway_analyze.kegg_gap_analyze import gap_depth_output_dir
from src.pathway_analyze.write_solution_to_manifest import select_solution


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


class SelectSolutionDepthTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.gap_root = self.root / "kegg_gap_C12345"
        self.manifest_path = self.root / "design_manifest.json"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def config(self, *, depth: int | None) -> SimpleNamespace:
        values = {
            "target_name": "C12345",
            "solution_id": 1,
            "gap_output_path": self.gap_root,
            "validation_output_path": self.gap_root / "legacy_validation",
            "manifest_output_path": self.manifest_path,
        }
        if depth is not None:
            values["depth"] = depth
        return SimpleNamespace(**values)

    def write_depth(self, depth: int, reaction_id: str) -> Path:
        depth_dir = gap_depth_output_dir(self.gap_root, depth)
        write_csv(
            depth_dir / "solutions.csv",
            [{
                "solution_id": 1,
                "target_compound_id": "C12345",
                "target_compound_name": "Depth test target",
                "total_steps": 1,
                "heterologous_steps": 1,
                "heterologous_reaction_ids": reaction_id,
                "heterologous_ko_ids": "K00001",
                "heterologous_enzyme_ecs": "1.1.1.1",
                "reaction_resolution_status": "resolved",
                "normalization_event_count": 0,
                "blocking_reaction_count": 0,
                "eligible_for_recommendation": True,
                "reachable_anchor_compounds": "C00001",
                "reachable_anchor_labels": "C00001 (Water)",
            }],
        )
        write_csv(
            depth_dir / "all_solution_steps.csv",
            [{
                "solution_id": 1,
                "step_index": 1,
                "status": "heterologous",
                "reaction_id": reaction_id,
                "reaction_name": f"Reaction {reaction_id}",
                "equation": "C00001 => C12345",
                "direction": "left_to_right",
                "produced_compound_id": "C12345",
                "produced_compound_name": "Depth test target",
                "precursor_compound_ids": "C00001",
                "precursor_compound_labels": "C00001 (Water)",
                "ko_ids": "K00001",
                "enzyme_ecs": "1.1.1.1",
                "source_reaction_ids": reaction_id,
                "resolution_action": "chassis_forward_expansion",
                "resolution_evidence": f"expansion_depth:{depth}",
                "step_source": "chassis_forward_expansion",
                "expansion_depth": depth,
                "expansion_anchor_compounds": "C12345",
            }],
        )
        write_csv(
            depth_dir / "gem_validation" / "gem_validation_summary.csv",
            [{
                "validation_mode": "per-solution",
                "solution_ids": "1",
                "validation_status": "PASS_ROUTE_REQUIRED_AT_TARGET_FLUX",
                "fba_status": "optimal",
                "fba_product_flux": 1.0,
                "cofactor_mode": "strict_l1",
                "cofactor_relaxed": False,
            }],
        )
        return depth_dir

    def read_manifest(self) -> dict[str, object]:
        return json.loads(self.manifest_path.read_text(encoding="utf-8"))

    def test_missing_depth_defaults_to_depth_zero(self) -> None:
        depth_dir = self.write_depth(0, "R00000")

        result = select_solution(self.config(depth=None))
        solution = self.read_manifest()["solution"]

        self.assertEqual(0, result["扩展深度"])
        self.assertEqual(0, solution["expansion_depth"])
        self.assertEqual(str(depth_dir.resolve()), solution["gap_dir"])
        self.assertEqual("R00000", solution["steps"][0]["reaction_id"])

    def test_selected_depth_isolated_and_written_to_manifest(self) -> None:
        self.write_depth(0, "R00000")
        depth_dir = self.write_depth(2, "R00002")

        result = select_solution(self.config(depth=2))
        solution = self.read_manifest()["solution"]
        selected_step = solution["steps"][0]

        self.assertEqual(2, result["扩展深度"])
        self.assertEqual(2, solution["expansion_depth"])
        self.assertEqual(str(depth_dir.resolve()), solution["gap_dir"])
        self.assertEqual("R00002", selected_step["reaction_id"])
        self.assertNotIn("step_source", selected_step)
        self.assertNotIn("expansion_depth", selected_step)
        self.assertNotIn("expansion_anchor_compounds", selected_step)

    def test_negative_depth_is_rejected_before_file_access(self) -> None:
        with self.assertRaisesRegex(ValueError, "greater than or equal to 0"):
            select_solution(self.config(depth=-1))

    def test_depth_zero_does_not_fall_back_to_legacy_root(self) -> None:
        write_csv(
            self.gap_root / "solutions.csv",
            [{"solution_id": 1, "target_compound_id": "C12345"}],
        )

        with self.assertRaises(FileNotFoundError) as raised:
            select_solution(self.config(depth=0))

        self.assertIn(str(self.gap_root / "depth0" / "solutions.csv"), str(raised.exception))


if __name__ == "__main__":
    unittest.main()
