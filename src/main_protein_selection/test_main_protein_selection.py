from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from src.main_protein_selection.reaction_direction_verifier import (
    DIRECTION_UNKNOWN,
    direction_decision_for_candidate,
)
from src.main_protein_selection.select_main_enzymes import (
    run_main_protein_selection,
    select_main_enzymes,
)
from src.main_protein_selection.uniprot_protein_candidates import ProteinCandidate


class MainProteinSelectionRuntimeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.manifest_path = self.root / "design_manifest.json"
        self.output_dir = self.root / "main_protein_selection"
        self.cache_dir = self.root / "cache"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_manifest(self, steps: list[dict[str, object]]) -> None:
        self.manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": "design_manifest.v1",
                    "revision": 1,
                    "solution": {
                        "solution_id": 1,
                        "expansion_depth": 0,
                        "steps": steps,
                    },
                }
            ),
            encoding="utf-8",
        )

    @staticmethod
    def step(
        step_index: int,
        reaction_id: str,
        enzyme_ecs: str = "",
    ) -> dict[str, object]:
        return {
            "step_index": step_index,
            "status": "heterologous",
            "reaction_id": reaction_id,
            "reaction_name": f"Reaction {reaction_id}",
            "equation": "C00001 => C00002",
            "direction": "left_to_right",
            "produced_compound_id": "C00002",
            "produced_compound_name": "Product",
            "precursor_compound_ids": "C00001",
            "precursor_compound_labels": "C00001 (Water)",
            "enzyme_ecs": enzyme_ecs,
            "locked_enzyme_ecs": enzyme_ecs,
            "ec_status": "complete" if enzyme_ecs else "missing",
            "ko_ids": "",
            "rhea_ids": "",
        }

    def test_fetch_disabled_uses_explicit_paths_and_writes_artifacts(self) -> None:
        self.write_manifest([self.step(1, "R00001", "1.1.1.1")])

        result = select_main_enzymes(
            manifest_path=self.manifest_path,
            output_dir=self.output_dir,
            cache_dir=self.cache_dir,
            fetch_proteins=False,
        )

        self.assertTrue(result["ok"])
        self.assertEqual("complete", result["status"])
        self.assertEqual(0, result["step_candidate_count"])
        expected_files = {
            "step_main_enzyme_candidates.csv",
            "main_enzyme_candidates.csv",
            "reaction_evidence.json",
            "direction_evidence.json",
            "ko_evidence.json",
            "selenzyme_evidence.json",
            "route_repair_requests.json",
        }
        self.assertEqual(
            expected_files,
            {path.name for path in self.output_dir.iterdir()},
        )
        direction = json.loads(
            (self.output_dir / "direction_evidence.json").read_text("utf-8")
        )
        self.assertEqual("disabled_deterministic", direction["agent"]["status"])

    def test_selenzyme_configuration_failure_preserves_existing_candidates(self) -> None:
        self.write_manifest(
            [
                self.step(1, "R00001", "1.1.1.1"),
                self.step(2, "R00002"),
            ]
        )
        candidate = ProteinCandidate(
            accession="P12345",
            entry_name="TEST_ECOLI",
            protein_name="Test enzyme",
            organism_name="Escherichia coli",
            organism_id=562,
            reviewed=True,
            length=300,
            ec_numbers=["1.1.1.1"],
            score=90.0,
            reasons=["test candidate"],
            sequence="M" * 300,
        )

        with (
            patch(
                "src.main_protein_selection.select_main_enzymes."
                "recommend_uniprot_proteins",
                return_value=[candidate],
            ),
            patch(
                "src.main_protein_selection.select_main_enzymes."
                "reaction_evidence_for_requirements",
                return_value=[],
            ),
            patch(
                "src.main_protein_selection.select_main_enzymes.SelenzymeClient",
                side_effect=RuntimeError("SELENZYME_REST_URL is required"),
            ),
        ):
            result = select_main_enzymes(
                manifest_path=self.manifest_path,
                output_dir=self.output_dir,
                cache_dir=self.cache_dir,
            )

        self.assertFalse(result["ok"])
        self.assertEqual("source_unavailable", result["status"])
        self.assertEqual(1, result["step_candidate_count"])
        with (
            self.output_dir / "step_main_enzyme_candidates.csv"
        ).open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(["P12345"], [row["accession"] for row in rows])
        self.assertEqual([2], result["uncovered_step_indexes"])

    def test_direction_agent_payload_cannot_override_unknown(self) -> None:
        requirement = {
            "direction_evidence_status": "resolved",
            "equation": "C00001 => C00002",
            "direction": "left_to_right",
            "ec_numbers": ["1.1.1.1"],
            "required_rhea_direction_ids": [],
            "opposite_rhea_direction_ids": [],
            "rhea_bidirectional_ids": [],
            "direction_agent_assessment": {
                "verdict": "supported",
                "confidence": "high",
            },
        }
        candidate = {"ec_numbers": ["1.1.1.1"], "rhea_ids": []}

        decision = direction_decision_for_candidate(requirement, candidate)

        self.assertEqual(DIRECTION_UNKNOWN, decision["verdict"])
        self.assertEqual("insufficient_direction_evidence", decision["evidence_level"])

    def test_run_config_wrapper_maps_project_paths(self) -> None:
        config = SimpleNamespace(
            manifest_output_path=self.manifest_path,
            project_output_path=self.root / "project",
            cache_dir=self.root / "cache_root",
        )
        expected = {"ok": True}
        with patch(
            "src.main_protein_selection.select_main_enzymes.select_main_enzymes",
            return_value=expected,
        ) as mocked:
            result = run_main_protein_selection(config, fetch_proteins=False)

        self.assertEqual(expected, result)
        mocked.assert_called_once_with(
            manifest_path=self.manifest_path,
            output_dir=config.project_output_path / "main_protein_selection",
            cache_dir=config.cache_dir / "main_protein_selection",
            fetch_proteins=False,
        )


if __name__ == "__main__":
    unittest.main()
