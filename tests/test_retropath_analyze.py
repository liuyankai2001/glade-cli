from __future__ import annotations

import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_STEP_COLUMNS,
    REJECTED_ROUTE_COLUMNS,
    analyze_retropath_candidates,
    write_retropath_candidate_artifacts,
)
from src.pathway_analyze.retropath_merge import (
    HybridCandidateRoute,
    HybridCandidateStep,
    RetroPathMergeRejection,
    RetroPathMergeResult,
)
from src.pathway_analyze.retropath_parser import SinkMatch


class RetroPathAnalyzeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temporary.name) / "retropath"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def sink(
        compound_id: str,
        inchikey: str,
        depth: int,
    ) -> SinkMatch:
        return SinkMatch(
            compound_id=compound_id,
            inchikey=inchikey,
            representative_kegg_id=compound_id,
            kegg_ids=(compound_id,),
            minimum_depth=depth,
            wrapper_in_sink=True,
            wrapper_sink_names=(compound_id,),
        )

    def merge_result(self) -> RetroPathMergeResult:
        sink_one = self.sink(
            "C10001",
            "OKKJLVBELUTLKV-UHFFFAOYSA-N",
            0,
        )
        sink_two = self.sink(
            "C10002",
            "VNWKTOKETHGBQD-UHFFFAOYSA-N",
            2,
        )
        kegg_step = HybridCandidateStep(
            step_id=f"RP2STEP:{'1' * 64}",
            step_source="kegg_expansion",
            status="endogenous",
            orientation="biosynthetic",
            direction="left_to_right",
            reaction_option_ids=("R10001",),
            reaction_smiles="",
            substrate_stoichiometry=(("C10000", 1.0),),
            product_stoichiometry=(("C10002", 1.0),),
            sink_anchor_kegg_ids=("C10002",),
            expansion_depth=2,
            is_endogenous=True,
            source_reaction_ids=("R10001",),
            source_ec_numbers=("1.1.1.1",),
        )
        rp2_step = HybridCandidateStep(
            step_id=f"RP2STEP:{'2' * 64}",
            step_source="retropath",
            status="predicted",
            orientation="biosynthetic",
            direction="biosynthetic",
            reaction_option_ids=(f"RP2:{'3' * 64}", f"RP2:{'4' * 64}"),
            reaction_smiles="C.CO>>CCO",
            substrate_stoichiometry=(("C10001", 1.0), ("C10002", 1.0)),
            product_stoichiometry=(("C90000", 1.0),),
            depends_on_step_ids=(kegg_step.step_id,),
            source_transformation_ids=("TRS_ROOT",),
            rule_ids=("RULE-A", "RULE-B"),
            source_reaction_ids=("MNXR1", "MNXR2"),
            source_ec_numbers=("1.1.1.1", "1.1.1.2"),
            minimum_rule_specificity=8,
            worst_rule_score=0.5,
            score_semantics="lower_is_better",
            cofactor_reconstruction_status="incomplete",
        )
        candidate = HybridCandidateRoute(
            candidate_id=f"RP2ROUTE:{'5' * 64}",
            source_retrosynthetic_path_id=f"RP2PATH:{'6' * 64}",
            target_compound_id="C90000",
            sink_matches=(sink_one, sink_two),
            steps=(kegg_step, rp2_step),
            minimum_rule_specificity=8,
            worst_rule_score=0.5,
            score_semantics="lower_is_better",
            contains_auxiliary_fragments=True,
        )
        rejection = RetroPathMergeRejection(
            source_stage="p4",
            source_path_id=f"RP2PATH:{'7' * 64}",
            reason_code="cycle_detected",
            reason_detail="test cycle",
            sink_kegg_ids=("C10001",),
            compound_id="RP2CPD:test",
            transformation_id="TRS_BAD",
        )
        return RetroPathMergeResult(
            candidates=(candidate,),
            biosynthetic_reactions=tuple(),
            rejections=(rejection,),
            upstream_truncated=False,
            truncated=False,
            max_candidates=5,
            max_witness_plans=3,
            max_total_steps=10,
            max_new_enzymes=10,
        )

    @staticmethod
    def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            return list(reader.fieldnames or []), list(reader)

    def test_writes_all_candidate_artifacts_with_fixed_schemas(self) -> None:
        artifacts = write_retropath_candidate_artifacts(
            self.merge_result(),
            self.output_dir,
        )

        route_columns, route_rows = self.read_rows(artifacts.candidate_routes_path)
        step_columns, step_rows = self.read_rows(artifacts.candidate_steps_path)
        rejection_columns, rejection_rows = self.read_rows(
            artifacts.rejected_routes_path
        )
        self.assertEqual(route_columns, list(CANDIDATE_ROUTE_COLUMNS))
        self.assertEqual(step_columns, list(CANDIDATE_STEP_COLUMNS))
        self.assertEqual(rejection_columns, list(REJECTED_ROUTE_COLUMNS))
        self.assertEqual(len(route_rows), 1)
        self.assertEqual(len(step_rows), 2)
        self.assertEqual(len(rejection_rows), 1)
        self.assertEqual(route_rows[0]["sink_kegg_ids"], "C10001;C10002")
        self.assertEqual(
            route_rows[0]["upstream_enumeration_truncated"],
            "false",
        )
        self.assertEqual(route_rows[0]["candidate_top_k_truncated"], "false")
        self.assertEqual(
            step_rows[1]["reaction_option_ids"],
            f"RP2:{'3' * 64};RP2:{'4' * 64}",
        )
        self.assertEqual(step_rows[1]["depends_on_step_ids"], step_rows[0]["step_id"])
        self.assertEqual(rejection_rows[0]["source_stage"], "p4")

    def test_output_bytes_and_hashes_are_deterministic_lf_utf8(self) -> None:
        first = write_retropath_candidate_artifacts(
            self.merge_result(),
            self.output_dir,
        )
        first_bytes = {
            path.name: path.read_bytes()
            for path in (
                first.candidate_routes_path,
                first.candidate_steps_path,
                first.rejected_routes_path,
            )
        }
        second = write_retropath_candidate_artifacts(
            self.merge_result(),
            self.output_dir,
        )

        for path, expected_hash in (
            (second.candidate_routes_path, second.candidate_routes_sha256),
            (second.candidate_steps_path, second.candidate_steps_sha256),
            (second.rejected_routes_path, second.rejected_routes_sha256),
        ):
            content = path.read_bytes()
            self.assertEqual(content, first_bytes[path.name])
            self.assertNotIn(b"\r\n", content)
            self.assertEqual(hashlib.sha256(content).hexdigest(), expected_hash)

    def test_empty_result_still_writes_parseable_headers(self) -> None:
        populated = self.merge_result()
        empty = RetroPathMergeResult(
            candidates=tuple(),
            biosynthetic_reactions=tuple(),
            rejections=tuple(),
            upstream_truncated=False,
            truncated=False,
            max_candidates=populated.max_candidates,
            max_witness_plans=populated.max_witness_plans,
            max_total_steps=populated.max_total_steps,
            max_new_enzymes=populated.max_new_enzymes,
        )

        artifacts = write_retropath_candidate_artifacts(empty, self.output_dir)

        for path, columns in (
            (artifacts.candidate_routes_path, CANDIDATE_ROUTE_COLUMNS),
            (artifacts.candidate_steps_path, CANDIDATE_STEP_COLUMNS),
            (artifacts.rejected_routes_path, REJECTED_ROUTE_COLUMNS),
        ):
            fieldnames, rows = self.read_rows(path)
            self.assertEqual(fieldnames, list(columns))
            self.assertFalse(rows)

    def test_analyze_orchestrates_merge_then_write(self) -> None:
        merge_result = self.merge_result()
        with patch(
            "src.pathway_analyze.retropath_analyze.merge_retropath_candidates",
            return_value=merge_result,
        ) as merge:
            artifacts = analyze_retropath_candidates(
                object(),  # type: ignore[arg-type]
                object(),  # type: ignore[arg-type]
                object(),  # type: ignore[arg-type]
                self.output_dir,
                max_candidates=7,
            )

        self.assertEqual(artifacts.merge_result, merge_result)
        self.assertTrue(artifacts.candidate_routes_path.is_file())
        self.assertEqual(merge.call_args.kwargs["max_candidates"], 7)


if __name__ == "__main__":
    unittest.main()
