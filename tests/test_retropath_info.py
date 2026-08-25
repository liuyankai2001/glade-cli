from __future__ import annotations

import csv
import hashlib
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from src.cli.app import build_parser
from src.info_show import (
    get_retropath_candidate_info,
    get_retropath_info,
    run_info,
)
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEP_COLUMNS,
    CANDIDATE_STEPS_FILE_NAME,
    REJECTED_ROUTE_COLUMNS,
    REJECTED_ROUTES_FILE_NAME,
)
from src.pathway_analyze.retropath_pipeline import (
    PIPELINE_RESULT_FILE_NAME,
    RETROPATH_PIPELINE_SCHEMA,
)


class RetroPathInfoTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.retropath_dir = self.root / "depth2" / "retropath"
        self.retropath_dir.mkdir(parents=True)
        self.config = SimpleNamespace(
            target_name="C12345",
            gap_output_path=self.root,
            depth=2,
            chassis=False,
            gap=False,
            retropath=True,
            retropath_candidate=None,
            solution=None,
            step=None,
        )
        self.route_rows = [self.route_row()]
        self.step_rows = self.candidate_step_rows()
        self.rejection_rows = [
            {
                "source_stage": "p4",
                "source_path_id": "",
                "reason_code": "unresolved_non_sink_leaf",
                "reason_detail": "leaf is not a trusted sink",
                "sink_kegg_ids": "",
                "compound_id": "RP2CPD:leaf",
                "transformation_id": "TRS-REJECTED",
            },
            {
                "source_stage": "p4",
                "source_path_id": "",
                "reason_code": "unresolved_non_sink_leaf",
                "reason_detail": "another leaf is not a trusted sink",
                "sink_kegg_ids": "",
                "compound_id": "RP2CPD:leaf-2",
                "transformation_id": "TRS-REJECTED-2",
            },
        ]
        self.write_success_result()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def route_row() -> dict[str, object]:
        return {
            "candidate_rank": 1,
            "candidate_id": "RP2ROUTE:" + "1" * 64,
            "source_retrosynthetic_path_id": "RP2PATH:" + "2" * 64,
            "target_compound_id": "C12345",
            "sink_kegg_ids": "C00001;C00002",
            "sink_depths": "C00001:0;C00002:2",
            "sink_inchikeys": "AAAA;BBBB",
            "kegg_prefix_reaction_ids": "R00001;R00002",
            "retropath_step_ids": "RP2STEP:3",
            "retropath_reaction_option_ids": "RP2:" + "3" * 64,
            "kegg_prefix_steps": 2,
            "retropath_steps": 1,
            "total_steps": 3,
            "maximum_sink_depth": 2,
            "minimum_rule_specificity": 8,
            "worst_rule_score": 0.25,
            "score_semantics": "lower_is_better",
            "contains_auxiliary_fragments": "true",
            "route_source": "kegg_retropath",
            "contains_predicted_steps": "true",
            "validation_status": "raw",
            "review_required": "true",
            "upstream_enumeration_truncated": "true",
            "candidate_top_k_truncated": "false",
        }

    @staticmethod
    def step_row(**overrides: object) -> dict[str, object]:
        row: dict[str, object] = {
            "candidate_id": "RP2ROUTE:" + "1" * 64,
            "step_index": 1,
            "step_id": "KEGG:R00001:left_to_right",
            "step_source": "kegg_expansion",
            "status": "heterologous",
            "orientation": "biosynthetic",
            "direction": "left_to_right",
            "reaction_option_ids": "R00001",
            "reaction_smiles": "",
            "substrate_compound_ids": "C00010",
            "product_compound_ids": "C00001",
            "substrate_stoichiometry_json": '[["C00010",1.0]]',
            "product_stoichiometry_json": '[["C00001",1.0]]',
            "depends_on_step_ids": "",
            "source_transformation_ids": "",
            "sink_anchor_kegg_ids": "C00001",
            "expansion_depth": 1,
            "is_endogenous": "false",
            "rule_ids": "",
            "source_reaction_ids": "R00001",
            "source_ec_numbers": "1.1.1.1",
            "minimum_rule_specificity": "",
            "worst_rule_score": "",
            "score_semantics": "",
            "balance_status": "not_checked",
            "cofactor_reconstruction_status": "not_checked",
        }
        row.update(overrides)
        return row

    def candidate_step_rows(self) -> list[dict[str, object]]:
        return [
            self.step_row(),
            self.step_row(
                step_index=2,
                step_id="KEGG:R00002:left_to_right",
                reaction_option_ids="R00002",
                substrate_compound_ids="C00020",
                product_compound_ids="C00002",
                substrate_stoichiometry_json='[["C00020",1.0]]',
                product_stoichiometry_json='[["C00002",1.0]]',
                sink_anchor_kegg_ids="C00002",
                expansion_depth=2,
                source_reaction_ids="R00002",
                source_ec_numbers="2.2.2.2",
            ),
            self.step_row(
                step_index=3,
                step_id="RP2STEP:3",
                step_source="retropath",
                status="predicted",
                reaction_option_ids="RP2:" + "3" * 64,
                reaction_smiles="CC.O>>CCO",
                substrate_compound_ids="C00001;C00002",
                product_compound_ids="C12345",
                substrate_stoichiometry_json=(
                    '[["C00001",1.0],["C00002",1.0]]'
                ),
                product_stoichiometry_json='[["C12345",1.0]]',
                depends_on_step_ids=(
                    "KEGG:R00001:left_to_right;KEGG:R00002:left_to_right"
                ),
                source_transformation_ids="TRS-1",
                sink_anchor_kegg_ids="C00001;C00002",
                expansion_depth=0,
                is_endogenous="",
                rule_ids="RR-01",
                source_reaction_ids="MNXR1",
                source_ec_numbers="3.3.3.3",
                minimum_rule_specificity=8,
                worst_rule_score=0.25,
                score_semantics="lower_is_better",
                balance_status="not_checked",
                cofactor_reconstruction_status="incomplete",
            ),
        ]

    @staticmethod
    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def write_csv(
        self,
        name: str,
        columns: tuple[str, ...],
        rows: list[dict[str, object]],
    ) -> Path:
        path = self.retropath_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def write_success_result(
        self,
        *,
        status: str = "retropath_candidates_found",
    ) -> dict[str, object]:
        route_path = self.write_csv(
            CANDIDATE_ROUTES_FILE_NAME,
            CANDIDATE_ROUTE_COLUMNS,
            self.route_rows,
        )
        step_path = self.write_csv(
            CANDIDATE_STEPS_FILE_NAME,
            CANDIDATE_STEP_COLUMNS,
            self.step_rows,
        )
        rejected_path = self.write_csv(
            REJECTED_ROUTES_FILE_NAME,
            REJECTED_ROUTE_COLUMNS,
            self.rejection_rows,
        )
        payload: dict[str, object] = {
            "schema_version": RETROPATH_PIPELINE_SCHEMA,
            "ok": True,
            "retropath_requested": True,
            "search_engine": "retropath",
            "status": status,
            "target_compound": "C12345",
            "expansion_depth": 2,
            "sink_source": "cumulative_expansion_A2",
            "output_dir": str(self.retropath_dir.resolve()),
            "pipeline_result_file": str(
                (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).resolve()
            ),
            "job_id": "rp2-" + "4" * 32,
            "service_status": "succeeded",
            "return_code": 0,
            "cache_hit": False,
            "scope_present": True,
            "sink_match_count": 2,
            "complete_path_count": 3,
            "candidate_count": len(self.route_rows),
            "rejection_count": len(self.rejection_rows),
            "upstream_enumeration_truncated": True,
            "candidate_top_k_truncated": False,
            "input_summary": {
                "reachable_compound_count": 1200,
                "sink_structure_count": 1100,
                "rejected_compound_count": 100,
            },
            "artifacts": {
                "candidate_routes": {
                    "path": str(route_path.resolve()),
                    "sha256": self.sha256(route_path),
                },
                "candidate_steps": {
                    "path": str(step_path.resolve()),
                    "sha256": self.sha256(step_path),
                },
                "rejected_routes": {
                    "path": str(rejected_path.resolve()),
                    "sha256": self.sha256(rejected_path),
                },
            },
        }
        (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        return payload

    def test_summary_is_audited_compact_and_risk_aware(self) -> None:
        result = get_retropath_info(self.config)

        self.assertTrue(result["运行成功"])
        self.assertEqual("cumulative_expansion_A2", result["Sink来源"])
        self.assertEqual(2, result["命中Sink数"])
        self.assertEqual(1, result["候选路线数"])
        self.assertEqual(2, result["拒绝路线数"])
        summary = result["候选路线摘要"][0]
        self.assertEqual(2, len(summary["命中边界化合物"]))
        self.assertEqual(["R00001", "R00002"], summary["KEGG前缀反应"])
        self.assertEqual(1, summary["RetroPath预测步骤数"])
        self.assertEqual("分数越低越好", summary["分数含义"])
        self.assertEqual(2, result["拒绝原因统计"][0]["数量"])
        self.assertIn("未完成计量", " ".join(result["警告"]))
        self.assertIn("并非穷尽", " ".join(result["警告"]))

    def test_candidate_and_single_step_show_dag_evidence(self) -> None:
        self.config.retropath = False
        self.config.retropath_candidate = 1

        candidate = get_retropath_candidate_info(self.config)
        self.assertEqual([1, 2], candidate["KEGG前缀步骤编号"])
        self.assertEqual([3], candidate["RetroPath预测步骤编号"])
        predicted = candidate["反应DAG步骤"][2]
        self.assertEqual("RetroPath 预测步骤", predicted["步骤来源"])
        self.assertEqual(["RR-01"], predicted["RetroRules规则ID"])
        self.assertEqual(["3.3.3.3"], predicted["来源EC编号"])
        self.assertEqual(2, len(predicted["依赖步骤ID"]))
        self.assertEqual(2, len(predicted["底物计量"]))
        self.assertIn("辅因子恢复状态", " ".join(predicted["风险提示"]))

        self.config.step = 3
        step = get_retropath_candidate_info(self.config)
        self.assertEqual(3, step["步骤编号"])
        self.assertEqual("RP2STEP:3", step["步骤详情"]["步骤ID"])

    def test_cli_and_run_info_use_independent_retropath_views(self) -> None:
        parser = build_parser()
        summary_args = parser.parse_args(
            ["info", "--input", "example.json", "--retropath", "--depth", "2"]
        )
        self.assertTrue(summary_args.retropath)
        candidate_args = parser.parse_args(
            [
                "info",
                "--input",
                "example.json",
                "--retropath-candidate",
                "1",
                "--step",
                "3",
                "--depth",
                "2",
            ]
        )
        self.assertEqual(1, candidate_args.retropath_candidate)
        self.assertEqual(3, candidate_args.step)
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(
                    [
                        "info",
                        "--input",
                        "example.json",
                        "--gap",
                        "--retropath",
                    ]
                )

        with patch("builtins.print"):
            result = run_info(self.config)
        self.assertEqual(1, result["候选路线数"])
        self.config.retropath = False
        self.config.retropath_candidate = 1
        with patch("builtins.print"):
            detail = run_info(self.config)
        self.assertEqual(1, detail["候选排名"])

        self.config.retropath_candidate = None
        self.config.retropath = True
        self.config.step = 1
        with self.assertRaisesRegex(ValueError, "--retropath-candidate"):
            run_info(self.config)

    def test_failure_manifest_is_viewable_without_candidate_files(self) -> None:
        for name in (
            CANDIDATE_ROUTES_FILE_NAME,
            CANDIDATE_STEPS_FILE_NAME,
            REJECTED_ROUTES_FILE_NAME,
        ):
            (self.retropath_dir / name).unlink()
        failure = {
            "schema_version": RETROPATH_PIPELINE_SCHEMA,
            "ok": False,
            "retropath_requested": True,
            "search_engine": "retropath",
            "status": "retropath_service_unavailable",
            "stage": "client",
            "detail": "connection refused",
            "target_compound": "C12345",
            "expansion_depth": 2,
            "sink_source": "cumulative_expansion_A2",
            "output_dir": str(self.retropath_dir.resolve()),
        }
        (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).write_text(
            json.dumps(failure),
            encoding="utf-8",
        )

        result = get_retropath_info(self.config)
        self.assertFalse(result["运行成功"])
        self.assertEqual("client", result["失败阶段"])
        self.assertEqual("connection refused", result["失败详情"])

        self.config.retropath_candidate = 1
        with self.assertRaisesRegex(ValueError, "运行未成功"):
            get_retropath_candidate_info(self.config)

    def test_source_in_sink_and_empty_candidates_are_explained(self) -> None:
        self.route_rows = []
        self.step_rows = []
        self.rejection_rows = []
        self.write_success_result(status="retropath_source_in_sink")

        result = get_retropath_info(self.config)
        self.assertEqual(0, result["候选路线数"])
        self.assertIn("目标化合物已属于", " ".join(result["警告"]))

        self.write_success_result(status="retropath_no_scope")
        no_scope = get_retropath_info(self.config)
        self.assertIn("没有得到命中可信", " ".join(no_scope["警告"]))

        self.config.retropath_candidate = 1
        with self.assertRaisesRegex(ValueError, r"可用候选排名：\[\]"):
            get_retropath_candidate_info(self.config)

    def test_hash_schema_target_and_relationship_mismatches_are_rejected(self) -> None:
        route_path = self.retropath_dir / CANDIDATE_ROUTES_FILE_NAME
        route_path.write_text(
            route_path.read_text(encoding="utf-8") + "tampered\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "校验失败"):
            get_retropath_info(self.config)

        payload = self.write_success_result()
        payload["schema_version"] = "retropath_pipeline_result.v0"
        (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).write_text(
            json.dumps(payload),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "版本不兼容"):
            get_retropath_info(self.config)

        payload = self.write_success_result()
        payload["target_compound"] = "C54321"
        (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).write_text(
            json.dumps(payload),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "运行目标"):
            get_retropath_info(self.config)

        payload = self.write_success_result()
        payload["candidate_count"] = 2
        (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).write_text(
            json.dumps(payload),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "候选数量"):
            get_retropath_info(self.config)

        payload = self.write_success_result()
        payload["expansion_depth"] = 1
        (self.retropath_dir / PIPELINE_RESULT_FILE_NAME).write_text(
            json.dumps(payload),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "运行深度"):
            get_retropath_info(self.config)

        self.step_rows[2]["depends_on_step_ids"] = "UNKNOWN-STEP"
        self.write_success_result()
        with self.assertRaisesRegex(ValueError, "有效拓扑顺序"):
            get_retropath_info(self.config)

        self.step_rows = self.candidate_step_rows()
        self.write_success_result()
        (self.retropath_dir / REJECTED_ROUTES_FILE_NAME).unlink()
        with self.assertRaisesRegex(FileNotFoundError, "缺少 RetroPath 候选文件"):
            get_retropath_info(self.config)

    def test_candidate_and_step_indexes_must_be_positive_and_available(self) -> None:
        self.config.retropath_candidate = 0
        with self.assertRaisesRegex(ValueError, "大于等于 1"):
            get_retropath_candidate_info(self.config)

        self.config.retropath_candidate = 2
        with self.assertRaisesRegex(ValueError, "可用候选排名"):
            get_retropath_candidate_info(self.config)

        self.config.retropath_candidate = 1
        self.config.step = 4
        with self.assertRaisesRegex(ValueError, "可用步骤"):
            get_retropath_candidate_info(self.config)


if __name__ == "__main__":
    unittest.main()
