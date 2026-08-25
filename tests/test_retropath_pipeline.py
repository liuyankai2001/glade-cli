from __future__ import annotations

import argparse
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from src.cli.commands import gap
from src.cli.common import apply_args_to_config
from src.config.run_config import RunConfig
from src.pathway_analyze.retropath_client import RetroPathClientError
from src.pathway_analyze.retropath_pipeline import (
    RetroPathPipelineError,
    run_gap_command,
    run_retropath_pipeline,
)


class RetroPathPipelineTests(unittest.TestCase):
    def make_parser(self) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(dest="command", required=True)
        gap.register(subparsers)
        return parser

    def test_cli_is_opt_in_and_does_not_override_config_when_absent(self) -> None:
        parser = self.make_parser()

        regular_args = parser.parse_args(["gap", "--input", "case.json"])
        self.assertIsNone(regular_args.retropath)
        self.assertIs(regular_args.func, run_gap_command)
        config = SimpleNamespace(retropath=False)
        apply_args_to_config(config, regular_args)
        self.assertFalse(config.retropath)
        self.assertEqual(config.depth, 0)

        retropath_args = parser.parse_args(
            ["gap", "--input", "case.json", "--depth", "0", "--retropath"]
        )
        self.assertTrue(retropath_args.retropath)

    def test_run_config_contains_reproducible_retropath_defaults(self) -> None:
        config = RunConfig("C00001")

        self.assertFalse(config.retropath)
        self.assertEqual(config.retropath_service_url, "http://127.0.0.1:8765")
        self.assertEqual(
            config.retropath_rules_path.name,
            "retrorules_rr02_rp2_flat_retro.csv",
        )
        self.assertEqual(config.retropath_max_steps, 3)
        self.assertEqual(config.retropath_max_candidates, 5)
        self.assertEqual(config.retropath_wait_timeout_seconds, 3900.0)

    @patch("src.pathway_analyze.retropath_pipeline.run_gap")
    def test_default_dispatch_is_the_original_kegg_entrypoint(
        self,
        run_gap: Mock,
    ) -> None:
        expected = {"ok": True, "search_engine": "kegg"}
        run_gap.return_value = expected
        config = SimpleNamespace(retropath=False)

        self.assertIs(run_gap_command(config), expected)
        run_gap.assert_called_once_with(config)

    def make_config(self, root: Path, *, depth: int = 0) -> SimpleNamespace:
        rules_path = root / "rules.csv"
        rules_path.write_text("Rule ID,Rule\n", encoding="utf-8")
        base_path = root / "producible.csv"
        base_path.write_text("kegg_id\nC00001\n", encoding="utf-8")
        return SimpleNamespace(
            retropath=True,
            target_name="C12345",
            depth=depth,
            gap_output_path=root / "gap",
            retropath_rules_path=rules_path,
            chassis_producible_csv=base_path,
            chassis_output_path=root / "chassis",
            cache_dir=root / "cache",
            retropath_service_url="http://127.0.0.1:8765",
            retropath_structure_timeout_seconds=30.0,
            retropath_structure_retries=3,
            retropath_structure_request_sleep_seconds=0.0,
            retropath_max_steps=3,
            retropath_topx=100,
            retropath_dmin=2,
            retropath_dmax=16,
            retropath_mwmax_source=1000,
            retropath_msc_timeout=10,
            retropath_request_timeout_seconds=30.0,
            retropath_get_attempts=3,
            retropath_retry_backoff_seconds=0.5,
            retropath_poll_interval_seconds=1.0,
            retropath_wait_timeout_seconds=3900.0,
            retropath_force=False,
            retropath_max_routes=1000,
            retropath_max_search_states=100000,
            retropath_max_candidates=5,
            retropath_max_witness_plans=3,
            max_total_steps=20,
            max_new_enzymes=20,
        )

    def fake_stage_objects(self, root: Path, *, status: str = "succeeded"):
        output_dir = root / "gap" / "depth0" / "retropath"
        input_dir = output_dir / "input"
        input_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        expanded_file = root / "producible.csv"

        input_paths = {
            name: input_dir / name
            for name in (
                "target_source.csv",
                "chassis_sink.csv",
                "compound_mapping.csv",
                "rejected_compounds.csv",
            )
        }
        for path in input_paths.values():
            path.write_text("header\n", encoding="utf-8")
        candidate_paths = {
            name: output_dir / name
            for name in (
                "candidate_routes.csv",
                "candidate_steps.csv",
                "rejected_routes.csv",
            )
        }
        for path in candidate_paths.values():
            path.write_text("header\n", encoding="utf-8")

        expansion = SimpleNamespace(
            depth=0,
            expanded_file=expanded_file,
        )
        input_bundle = SimpleNamespace(
            expansion_depth=0,
            reachable_compound_count=20,
            sink_structure_count=18,
            rejected_compound_count=2,
            target_source_path=input_paths["target_source.csv"],
            chassis_sink_path=input_paths["chassis_sink.csv"],
            compound_mapping_path=input_paths["compound_mapping.csv"],
            rejected_compounds_path=input_paths["rejected_compounds.csv"],
            target_source_sha256="a" * 64,
            chassis_sink_sha256="b" * 64,
        )
        provenance = Mock()
        provenance.rules_sha256 = "c" * 64
        provenance.to_dict.return_value = {
            "wrapper_version": "3.9.1",
            "rules_sha256": "c" * 64,
        }
        client_run = SimpleNamespace(
            result=SimpleNamespace(
                job_id="rp2-" + "1" * 32,
                status=status,
                return_code=0,
                errors=tuple(),
                provenance=provenance,
            ),
            cache_hit=False,
            raw_dir=output_dir / "raw",
            run_manifest_path=output_dir / "run_manifest.json",
            client_state_path=output_dir / "client_state.json",
        )
        enumeration = SimpleNamespace(
            complete_path_count=2,
            network=SimpleNamespace(sink_matches=("C00001",)),
            max_routes=1000,
            max_search_states=100000,
            truncated=True,
        )
        merge_result = SimpleNamespace(
            truncated=True,
            max_candidates=5,
            max_witness_plans=3,
            max_total_steps=20,
            max_new_enzymes=20,
        )
        candidates = SimpleNamespace(
            candidate_count=2,
            rejection_count=1,
            merge_result=merge_result,
            candidate_routes_path=candidate_paths["candidate_routes.csv"],
            candidate_steps_path=candidate_paths["candidate_steps.csv"],
            rejected_routes_path=candidate_paths["rejected_routes.csv"],
            candidate_routes_sha256="d" * 64,
            candidate_steps_sha256="e" * 64,
            rejected_routes_sha256="f" * 64,
        )
        return expansion, input_bundle, client_run, enumeration, candidates

    def test_retropath_orchestrates_p2_to_p5_and_writes_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config = self.make_config(root)
            expansion, inputs, client_run, enumeration, candidates = (
                self.fake_stage_objects(root)
            )
            http_client = Mock()
            http_client.__enter__ = Mock(return_value=http_client)
            http_client.__exit__ = Mock(return_value=None)
            http_client.run.return_value = client_run

            with (
                patch(
                    "src.pathway_analyze.retropath_pipeline.load_expansion_bundle",
                    return_value=expansion,
                ) as load_expansion,
                patch("src.pathway_analyze.retropath_pipeline.KeggMolStructureProvider"),
                patch(
                    "src.pathway_analyze.retropath_pipeline.build_retropath_inputs",
                    return_value=inputs,
                ) as build_inputs,
                patch(
                    "src.pathway_analyze.retropath_pipeline.RetroPathHttpClient",
                    return_value=http_client,
                ),
                patch(
                    "src.pathway_analyze.retropath_pipeline.parse_and_enumerate_retropath",
                    return_value=enumeration,
                ) as enumerate_routes,
                patch("src.pathway_analyze.retropath_pipeline.KeggRestClient"),
                patch(
                    "src.pathway_analyze.retropath_pipeline.analyze_retropath_candidates",
                    return_value=candidates,
                ) as analyze_candidates,
                patch("builtins.print"),
            ):
                result = run_retropath_pipeline(config)

            self.assertEqual(result["status"], "retropath_candidates_found")
            self.assertEqual(result["search_engine"], "retropath")
            self.assertEqual(result["sink_source"], "chassis_A0")
            self.assertEqual(result["candidate_count"], 2)
            self.assertTrue(result["upstream_enumeration_truncated"])
            self.assertTrue(result["candidate_top_k_truncated"])
            self.assertEqual(
                Path(result["output_dir"]),
                (root / "gap" / "depth0" / "retropath").resolve(),
            )
            result_path = Path(result["pipeline_result_file"])
            self.assertEqual(json.loads(result_path.read_text(encoding="utf-8")), result)

            load_expansion.assert_called_once_with(
                base_path=config.chassis_producible_csv,
                output_dir=config.chassis_output_path,
                depth=0,
            )
            self.assertEqual(build_inputs.call_args.args[3], result_path.parent / "input")
            self.assertEqual(http_client.run.call_args.args[1], result_path.parent)
            self.assertEqual(enumerate_routes.call_args.args[2], config.retropath_rules_path)
            self.assertEqual(
                analyze_candidates.call_args.kwargs["max_total_steps"],
                20,
            )

    def test_depth_greater_than_zero_uses_cumulative_sink_label(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config = self.make_config(root, depth=3)
            expansion, inputs, client_run, enumeration, candidates = (
                self.fake_stage_objects(root)
            )
            expansion.depth = 3
            expansion.expanded_file = root / "expanded_A3.csv"
            expansion.expanded_file.write_text("kegg_id\nC00001\n", encoding="utf-8")
            inputs.expansion_depth = 3
            candidates.candidate_count = 0
            candidates.merge_result.truncated = False
            client_run.result.status = "no_solution"
            http_client = Mock()
            http_client.__enter__ = Mock(return_value=http_client)
            http_client.__exit__ = Mock(return_value=None)
            http_client.run.return_value = client_run

            with (
                patch(
                    "src.pathway_analyze.retropath_pipeline.load_expansion_bundle",
                    return_value=expansion,
                ),
                patch("src.pathway_analyze.retropath_pipeline.KeggMolStructureProvider"),
                patch(
                    "src.pathway_analyze.retropath_pipeline.build_retropath_inputs",
                    return_value=inputs,
                ),
                patch(
                    "src.pathway_analyze.retropath_pipeline.RetroPathHttpClient",
                    return_value=http_client,
                ),
                patch(
                    "src.pathway_analyze.retropath_pipeline.parse_and_enumerate_retropath",
                    return_value=enumeration,
                ),
                patch("src.pathway_analyze.retropath_pipeline.KeggRestClient"),
                patch(
                    "src.pathway_analyze.retropath_pipeline.analyze_retropath_candidates",
                    return_value=candidates,
                ),
                patch("builtins.print"),
            ):
                result = run_retropath_pipeline(config)

            self.assertEqual(result["status"], "retropath_no_scope")
            self.assertEqual(result["sink_source"], "cumulative_expansion_A3")
            self.assertIn("depth3", result["output_dir"])

    def test_missing_expansion_has_stable_failure_and_audit_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config = self.make_config(root, depth=2)
            with patch(
                "src.pathway_analyze.retropath_pipeline.load_expansion_bundle",
                side_effect=ValueError("Run the expand command for depth 2 first."),
            ):
                with self.assertRaises(RetroPathPipelineError) as caught:
                    run_retropath_pipeline(config)

            self.assertEqual(caught.exception.status, "retropath_expansion_missing")
            payload = json.loads(
                caught.exception.result_path.read_text(encoding="utf-8")
            )
            self.assertFalse(payload["ok"])
            self.assertEqual(payload["stage"], "expansion")
            self.assertEqual(payload["sink_source"], "cumulative_expansion_A2")

    def test_client_errors_map_to_stable_pipeline_statuses(self) -> None:
        cases = (
            ("service_unavailable", "retropath_service_unavailable"),
            ("client_poll_timeout", "retropath_timeout"),
            ("input_invalid", "retropath_input_invalid"),
        )
        for client_code, expected_status in cases:
            with self.subTest(client_code=client_code):
                with tempfile.TemporaryDirectory() as temporary_directory:
                    root = Path(temporary_directory)
                    config = self.make_config(root)
                    expansion, inputs, _, _, _ = self.fake_stage_objects(root)
                    http_client = Mock()
                    http_client.__enter__ = Mock(return_value=http_client)
                    http_client.__exit__ = Mock(return_value=None)
                    http_client.run.side_effect = RetroPathClientError(
                        client_code,
                        "client failure",
                    )

                    with (
                        patch(
                            "src.pathway_analyze.retropath_pipeline.load_expansion_bundle",
                            return_value=expansion,
                        ),
                        patch(
                            "src.pathway_analyze.retropath_pipeline.KeggMolStructureProvider"
                        ),
                        patch(
                            "src.pathway_analyze.retropath_pipeline.build_retropath_inputs",
                            return_value=inputs,
                        ),
                        patch(
                            "src.pathway_analyze.retropath_pipeline.RetroPathHttpClient",
                            return_value=http_client,
                        ),
                    ):
                        with self.assertRaises(RetroPathPipelineError) as caught:
                            run_retropath_pipeline(config)

                    self.assertEqual(caught.exception.status, expected_status)
                    self.assertEqual(caught.exception.stage, "client")


if __name__ == "__main__":
    unittest.main()
