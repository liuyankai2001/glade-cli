from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from src.info_show import get_chassis_info, run_info


class ChassisInfoTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.summary_path = self.root / "analyze_chassis_metabolites_summary.csv"
        self.compounds_path = self.root / "producible_kegg_compounds.csv"
        self.config = SimpleNamespace(
            target_name="C00811",
            model_path=self.root / "model.json",
            medium_path=self.root / "medium.json",
            chassis_metabolites_summary_csv=self.summary_path,
            chassis_producible_csv=self.compounds_path,
            show_chassis=True,
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_summary(self, rows: list[tuple[str, object]]) -> None:
        with self.summary_path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["item", "value"])
            writer.writerows(rows)

    def write_compounds(self) -> None:
        with self.compounds_path.open(
            "w", encoding="utf-8-sig", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["source", "met_id", "met_name", "compartment", "kegg_id"],
            )
            writer.writeheader()
            writer.writerows([
                {
                    "source": "producible",
                    "met_id": "h2o_c",
                    "met_name": "Water",
                    "compartment": "c",
                    "kegg_id": "C00001",
                },
                {
                    "source": "producible",
                    "met_id": "target_c",
                    "met_name": "Target",
                    "compartment": "c",
                    "kegg_id": "C00811",
                },
            ])

    def test_reads_summary_and_reports_target_availability(self) -> None:
        self.write_summary([
            ("baseline_growth", 0.87),
            ("required_growth", 0.087),
            ("growth_fraction", 0.1),
            ("tested_metabolites", 10),
            ("producible_metabolites", 5),
            ("producible_with_kegg", 2),
            ("producible_without_kegg", 3),
            ("optimization_failed", 0),
            ("kegg_mapping_rows", 2),
            ("unique_producible_kegg_compounds", 2),
        ])
        self.write_compounds()

        result = get_chassis_info(self.config)

        self.assertTrue(result["ok"])
        self.assertEqual("chassis_info.v1", result["schema_version"])
        self.assertTrue(result["target_producible_by_chassis"])
        self.assertEqual("target_c", result["target_matches"][0]["met_id"])
        self.assertEqual(10, result["screening"]["tested_metabolites"])
        self.assertEqual(2, result["kegg_compounds"]["mapping_rows"])
        self.assertEqual(2, result["kegg_compounds"]["unique_compounds"])
        self.assertEqual(1, len(result["warnings"]))

    def test_legacy_summary_remains_readable(self) -> None:
        self.write_summary([
            ("baseline_growth", 0.87),
            ("growth_fraction", 0.1),
            ("flux_threshold", "1e-08"),
            ("producible_metabolites", 2),
            ("producible_kegg_compounds", 2),
        ])
        self.write_compounds()

        result = run_info(self.config)

        self.assertEqual(2, result["可生成KEGG化合物数"])
        self.assertIsNone(result["检测代谢物数"])
        self.assertEqual("model.json", result["底盘模型"])
        self.assertNotIn("数据文件", result)

    def test_missing_results_instructs_user_to_run_chassis(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, "请先运行 chassis 命令"):
            get_chassis_info(self.config)

    def test_info_requires_an_explicit_view(self) -> None:
        self.config.show_chassis = False
        with self.assertRaisesRegex(ValueError, "请使用 --chassis"):
            run_info(self.config)


if __name__ == "__main__":
    unittest.main()
