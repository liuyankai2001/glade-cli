from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from pydantic import ValidationError

from src.main_protein_selection.common import candidate_rows_for_requirements
from src.main_protein_selection.models import (
    MainEnzymeCandidate,
    MainEnzymeSelectionResult,
)
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
            "step_main_enzyme_candidate_audit.csv",
            "main_enzyme_candidates.csv",
            "reaction_evidence.json",
            "direction_evidence.json",
            "ko_evidence.json",
            "selenzyme_evidence.json",
            "route_repair_requests.json",
            "main_enzyme_selection.json",
        }
        self.assertEqual(
            expected_files,
            {path.name for path in self.output_dir.iterdir()},
        )
        direction = json.loads(
            (self.output_dir / "direction_evidence.json").read_text("utf-8")
        )
        self.assertEqual("disabled_deterministic", direction["agent"]["status"])
        canonical = MainEnzymeSelectionResult.model_validate_json(
            (self.output_dir / "main_enzyme_selection.json").read_text("utf-8")
        )
        self.assertEqual("main_enzyme_selection.v1", canonical.schema_version)
        self.assertEqual(0, canonical.expansion_depth)
        self.assertEqual({}, canonical.candidates_by_step)
        self.assertEqual(64, len(canonical.solution_fingerprint))

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
        canonical = MainEnzymeSelectionResult.model_validate_json(
            (self.output_dir / "main_enzyme_selection.json").read_text("utf-8")
        )
        self.assertEqual("source_unavailable", canonical.status)
        self.assertEqual(
            ["P12345"],
            [item.accession for item in canonical.candidates_by_step[1]],
        )

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
            top_n=7,
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
            top_n=7,
        )

    def test_candidate_model_rejects_invalid_or_nonusable_rows(self) -> None:
        valid = {
            "step_index": 1,
            "reaction_id": "R00001",
            "accession": "P12345",
            "reviewed": True,
            "candidate_rank": 1,
            "protein_score": 90.0,
            "reaction_fit_status": "verified",
            "reaction_fit_score": 100.0,
        }
        with self.assertRaises(ValidationError):
            MainEnzymeCandidate.model_validate({**valid, "accession": ""})
        with self.assertRaises(ValidationError):
            MainEnzymeCandidate.model_validate({**valid, "candidate_rank": 0})
        with self.assertRaises(ValidationError):
            MainEnzymeCandidate.model_validate(
                {**valid, "reaction_fit_status": "rejected"}
            )
        with self.assertRaises(ValidationError):
            MainEnzymeCandidate.model_validate({**valid, "unexpected": True})

    @staticmethod
    def requirement(
        step_index: int,
        ec_numbers: list[str],
        **overrides: object,
    ) -> dict[str, object]:
        value: dict[str, object] = {
            "solution_id": 1,
            "step_index": step_index,
            "reaction_id": f"R{step_index:05d}",
            "reaction_name": f"Reaction {step_index}",
            "produced_compound_id": "C00002",
            "produced_compound_name": "Product",
            "equation": "C00001 => C00002",
            "direction": "left_to_right",
            "ec_numbers": ec_numbers,
            "locked_ec_numbers": ec_numbers,
            "ec_status": "complete",
            "ko_ids": [],
            "rhea_ids": [],
        }
        value.update(overrides)
        return value

    @staticmethod
    def candidate(
        accession: str,
        ec_number: str,
        score: float,
        **overrides: object,
    ) -> ProteinCandidate:
        values: dict[str, object] = {
            "accession": accession,
            "entry_name": f"{accession}_ECOLI",
            "protein_name": f"Enzyme {accession}",
            "organism_name": "Escherichia coli",
            "organism_id": 562,
            "reviewed": True,
            "length": 100,
            "ec_numbers": [ec_number],
            "score": score,
            "reasons": [f"evidence for {ec_number}"],
            "sequence": "M" * 100,
        }
        values.update(overrides)
        return ProteinCandidate(**values)

    def test_same_accession_across_ec_and_sources_is_merged(self) -> None:
        requirement = self.requirement(
            1,
            ["1.1.1.1", "2.2.2.2"],
            ko_ids=["K00001"],
        )
        ec_candidate = self.candidate(
            "P12345",
            "1.1.1.1",
            80.0,
            retrieval_strategy="ec_exact",
            retrieval_query_id="uniprot_ec_1",
            matched_rhea_ids=["10000"],
            sequence_sha256="same-hash",
        )
        ko_candidate = self.candidate(
            "p12345",
            "2.2.2.2",
            90.0,
            retrieval_strategy="kegg_ko_exact",
            retrieval_query_id="kegg_ko_1",
            matched_ko_ids=["K00001"],
            kegg_gene_ids=["eco:b0001"],
            sequence_sha256="same-hash",
            warnings=["KO evidence warning"],
        )

        rows = candidate_rows_for_requirements(
            [requirement],
            {
                "1.1.1.1": [ec_candidate],
                "2.2.2.2": [ko_candidate],
            },
            top_n=5,
        )

        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual("P12345", row["accession"])
        self.assertEqual("1.1.1.1;2.2.2.2", row["ec_number"])
        self.assertEqual("ec_exact;kegg_ko_exact", row["retrieval_strategy"])
        self.assertEqual("uniprot_ec_1;kegg_ko_1", row["retrieval_query_id"])
        self.assertEqual("K00001", row["matched_ko_ids"])
        self.assertEqual("eco:b0001", row["kegg_gene_ids"])
        self.assertIn("KO evidence warning", row["warnings"])
        self.assertEqual(1, row["evaluation_rank"])
        self.assertEqual(1, row["candidate_rank"])
        self.assertEqual("selected", row["selection_status"])

    def test_top_n_is_applied_after_cross_ec_deduplication(self) -> None:
        self.write_manifest(
            [self.step(1, "R00737", "4.3.1.23;4.3.1.25")]
        )

        def candidates_for_ec(*, ec_number: str, **_: object) -> list[ProteinCandidate]:
            indexes = range(1, 6) if ec_number == "4.3.1.23" else range(4, 9)
            return [
                self.candidate(
                    f"P{index:05d}",
                    ec_number,
                    100.0 - index,
                )
                for index in indexes
            ]

        with (
            patch(
                "src.main_protein_selection.select_main_enzymes."
                "recommend_uniprot_proteins",
                side_effect=candidates_for_ec,
            ),
            patch(
                "src.main_protein_selection.select_main_enzymes."
                "reaction_evidence_for_requirements",
                return_value=[],
            ),
            patch(
                "src.main_protein_selection.select_main_enzymes."
                "retrieve_rhea_candidates_for_requirement",
                return_value=([], [], {}),
            ),
        ):
            result = select_main_enzymes(
                manifest_path=self.manifest_path,
                output_dir=self.output_dir,
                cache_dir=self.cache_dir,
                top_n=5,
            )

        self.assertEqual(5, result["step_candidate_count"])
        self.assertEqual(8, result["evaluated_step_candidate_count"])
        with (
            self.output_dir / "step_main_enzyme_candidates.csv"
        ).open("r", encoding="utf-8-sig", newline="") as handle:
            selected = list(csv.DictReader(handle))
        with (
            self.output_dir / "step_main_enzyme_candidate_audit.csv"
        ).open("r", encoding="utf-8-sig", newline="") as handle:
            audit = list(csv.DictReader(handle))
        self.assertEqual(5, len(selected))
        self.assertEqual(["1", "2", "3", "4", "5"], [
            row["candidate_rank"] for row in selected
        ])
        self.assertEqual(8, len(audit))
        self.assertEqual(8, len({row["accession"] for row in audit}))
        self.assertEqual(
            3,
            sum(
                row["selection_status"] == "eligible_not_selected"
                for row in audit
            ),
        )
        canonical = MainEnzymeSelectionResult.model_validate_json(
            (self.output_dir / "main_enzyme_selection.json").read_text("utf-8")
        )
        self.assertEqual(5, len(canonical.candidates_by_step[1]))

    def test_verified_candidates_rank_before_higher_scoring_risk_candidates(self) -> None:
        requirement = self.requirement(1, ["1.1.1.1"])
        verified = self.candidate("P00001", "1.1.1.1", 10.0)
        risk = self.candidate(
            "P00002",
            "1.1.1.1",
            99.0,
            retrieval_strategy="selenzyme_kegg_risk",
            reaction_confidence="selenzyme_risk",
            selenzyme_reaction_similarity=0.99,
        )

        rows = candidate_rows_for_requirements(
            [requirement],
            {"1.1.1.1": [risk, verified]},
            top_n=2,
        )

        self.assertEqual(["P00001", "P00002"], [row["accession"] for row in rows])
        self.assertEqual(
            ["verified", "verified_with_risk"],
            [row["reaction_fit_status"] for row in rows],
        )

    def test_direction_conflict_rejects_merged_accession(self) -> None:
        requirement = self.requirement(
            1,
            ["1.1.1.1", "2.2.2.2"],
            direction_evidence_status="resolved",
            required_rhea_direction_ids=["123"],
            opposite_rhea_direction_ids=["456"],
            rhea_bidirectional_ids=[],
        )
        supported = self.candidate(
            "P12345",
            "1.1.1.1",
            80.0,
            rhea_ids=["123"],
        )
        contradicted = self.candidate(
            "P12345",
            "2.2.2.2",
            90.0,
            rhea_ids=["456"],
        )

        rows = candidate_rows_for_requirements(
            [requirement],
            {"1.1.1.1": [supported], "2.2.2.2": [contradicted]},
            top_n=5,
        )

        self.assertEqual(1, len(rows))
        self.assertEqual("contradicted", rows[0]["direction_verdict"])
        self.assertEqual("rejected", rows[0]["reaction_fit_status"])
        self.assertEqual("", rows[0]["candidate_rank"])

    def test_ineligible_candidates_only_receive_evaluation_rank(self) -> None:
        manual_requirement = self.requirement(
            1,
            ["1.1.-.-"],
            ec_status="partial",
        )
        rejected_requirement = self.requirement(
            2,
            ["1.1.1.1"],
            direction_evidence_status="resolved",
            required_rhea_direction_ids=["123"],
            opposite_rhea_direction_ids=["456"],
            rhea_bidirectional_ids=[],
        )
        manual = self.candidate("P00001", "1.1.-.-", 50.0)
        rejected = self.candidate(
            "P00002",
            "1.1.1.1",
            50.0,
            rhea_ids=["456"],
        )

        rows = candidate_rows_for_requirements(
            [manual_requirement, rejected_requirement],
            {
                "1.1.-.-": [manual],
                "1.1.1.1": [rejected],
            },
            top_n=5,
        )

        self.assertEqual([1, 1], [row["evaluation_rank"] for row in rows])
        self.assertEqual(["", ""], [row["candidate_rank"] for row in rows])
        self.assertEqual(
            ["manual_review", "rejected"],
            [row["selection_status"] for row in rows],
        )

    def test_conflicting_sequences_for_same_accession_raise(self) -> None:
        requirement = self.requirement(1, ["1.1.1.1", "2.2.2.2"])
        first = self.candidate(
            "P12345",
            "1.1.1.1",
            80.0,
            sequence="M" * 100,
        )
        second = self.candidate(
            "P12345",
            "2.2.2.2",
            90.0,
            sequence="A" * 100,
        )

        with self.assertRaisesRegex(ValueError, "Conflicting sequence evidence"):
            candidate_rows_for_requirements(
                [requirement],
                {"1.1.1.1": [first], "2.2.2.2": [second]},
                top_n=5,
            )


if __name__ == "__main__":
    unittest.main()
