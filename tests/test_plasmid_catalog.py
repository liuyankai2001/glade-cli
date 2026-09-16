import csv
import hashlib
import shutil
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"


class ModuleCatalogContractTests(unittest.TestCase):
    def test_catalog_exposes_all_curated_modules_with_web_payloads(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(DATA_DIR)

        self.assertEqual(
            {module.id: module.length_bp for module in catalog.resistance.values()},
            {
                "basic_seva_ap": 1168,
                "basic_seva_km": 1056,
                "basic_seva_cm": 912,
                "basic_seva_sm_sp": 1126,
                "basic_seva_tet_5a": 1404,
                "basic_seva_gm": 934,
                "basic_seva_gm_11": 953,
            },
        )
        self.assertEqual(
            {module.id: module.length_bp for module in catalog.replication.values()},
            {
                "basic_seva_rsf1010": 3671,
                "basic_seva_p15a": 732,
                "basic_seva_psc101": 1460,
                "basic_seva_pbr322_rop": 1385,
                "basic_seva_psc101_pkd46_ts": 1551,
            },
        )
        module = catalog.get_resistance("basic_seva_ap")
        self.assertEqual(module.length_bp, 1168)
        self.assertEqual(
            hashlib.sha256(module.sequence.encode()).hexdigest(),
            "e589d17e939c40824023f92436804adf95f353c91ccf30eb9780155b470bc57f",
        )
        self.assertEqual(module.to_payload()["resistance_gene"], "bla")
        self.assertEqual(module.to_payload()["type"], "resistance")
        self.assertEqual(
            set(module.to_payload()),
            {
                "id",
                "name",
                "type",
                "sequence",
                "antibiotic",
                "resistance_gene",
                "notes",
                "length_bp",
                "gc_percent",
            },
        )
        self.assertNotIn("features", module.to_payload())
        self.assertIsInstance(module.features, tuple)
        self.assertEqual(module.notes, ())
        gm = catalog.get_resistance("basic_seva_gm")
        self.assertIsInstance(gm.notes, tuple)
        self.assertTrue(any("p15A/pBR322" in note for note in gm.notes))
        self.assertEqual(gm.to_payload()["notes"], list(gm.notes))
        replication_payload = catalog.get_replication("basic_seva_p15a").to_payload()
        self.assertEqual(
            set(replication_payload),
            {
                "id",
                "name",
                "type",
                "sequence",
                "host_range",
                "copy_number",
                "notes",
                "length_bp",
                "gc_percent",
            },
        )

    def test_catalog_collections_are_read_only(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(DATA_DIR)
        with self.assertRaises(TypeError):
            catalog.resistance["replacement"] = catalog.get_resistance("basic_seva_ap")
        with self.assertRaises(AttributeError):
            catalog.replication = {}

    def test_scaffold_preserves_source_gap_and_t1_boundaries(self):
        from src.plasmid_design.catalog import ModuleCatalog

        scaffold = ModuleCatalog(DATA_DIR).scaffold

        self.assertEqual(len(scaffold["resistance_to_replication"]), 76)
        self.assertEqual(len(scaffold["replication_to_t1"]), 14)
        self.assertEqual(len(scaffold["t1"]), 105)
        self.assertEqual(len(scaffold["landing_pad_spacer"]), 24)
        self.assertEqual(
            scaffold["landing_pad_spacer"], scaffold["resistance_to_replication"][12:36]
        )
        self.assertEqual(
            scaffold["provenance"]["resistance_to_replication"]["start_1based"], 1232
        )
        self.assertEqual(scaffold["provenance"]["t1"]["end_1based"], 5097)

    def test_module_annotations_are_bounded_by_the_module_sequence(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(DATA_DIR)
        modules = tuple(catalog.resistance.values()) + tuple(
            catalog.replication.values()
        )
        self.assertTrue(any(module.features for module in modules))
        for module in modules:
            for feature in module.features:
                for part in feature["parts"]:
                    self.assertGreaterEqual(part["start"], 0)
                    self.assertGreater(part["end"], part["start"])
                    self.assertLessEqual(part["end"], module.length_bp)
        ap_features = catalog.get_resistance("basic_seva_ap").features
        self.assertEqual(ap_features[0]["parts"][0]["end"], 103)
        self.assertEqual(ap_features[1]["parts"][0]["start"], 129)
        pkd46_features = catalog.get_replication("basic_seva_psc101_pkd46_ts").features
        self.assertEqual(pkd46_features[0]["type"], "CDS")
        self.assertEqual(pkd46_features[0]["parts"][0]["start"], 61)
        self.assertEqual(pkd46_features[1]["parts"][0]["end"], 1282)

    def test_invalid_csv_and_unknown_ids_are_rejected(self):
        from src.plasmid_design.catalog import ModuleCatalog

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in ("resistance_modules.csv", "replication_modules.csv"):
                shutil.copy2(DATA_DIR / name, staged / name)
            with (staged / "resistance_modules.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["sequence"] = "NOT_DNA"
            with (staged / "resistance_modules.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(ValueError):
                ModuleCatalog(staged)

        catalog = ModuleCatalog(DATA_DIR)
        with self.assertRaises(KeyError):
            catalog.get_replication("does-not-exist")

    def test_required_catalog_metadata_and_extra_csv_cells_are_rejected(self):
        from src.plasmid_design.catalog import ModuleCatalog

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in ("resistance_modules.csv", "replication_modules.csv"):
                shutil.copy2(DATA_DIR / name, staged / name)
            with (staged / "resistance_modules.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["antibiotic"] = ""
            with (staged / "resistance_modules.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "invalid resistance module row"):
                ModuleCatalog(staged)

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in ("resistance_modules.csv", "replication_modules.csv"):
                shutil.copy2(DATA_DIR / name, staged / name)
            path = staged / "replication_modules.csv"
            lines = path.read_text(encoding="utf-8").splitlines()
            path.write_text(
                "\n".join([lines[0], lines[1] + ",extra", *lines[2:]]) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "malformed CSV row"):
                ModuleCatalog(staged)

    def test_real_sequence_change_is_rejected_by_its_feature_snapshot(self):
        from src.plasmid_design.catalog import ModuleCatalog

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in ("resistance_modules.csv", "replication_modules.csv"):
                shutil.copy2(DATA_DIR / name, staged / name)
            with (staged / "replication_modules.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["sequence"] = "A" + rows[0]["sequence"][1:]
            with (staged / "replication_modules.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "source feature snapshot mismatch"):
                ModuleCatalog(staged)

    def test_fingerprint_changes_for_a_valid_catalog_name_change(self):
        from src.plasmid_design.catalog import ModuleCatalog

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in ("resistance_modules.csv", "replication_modules.csv"):
                shutil.copy2(DATA_DIR / name, staged / name)
            before = ModuleCatalog(staged).fingerprint
            with (staged / "replication_modules.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["name"] = "Renamed RSF1010"
            with (staged / "replication_modules.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            catalog = ModuleCatalog(staged)
            self.assertEqual(
                catalog.get_replication("basic_seva_rsf1010").name, "Renamed RSF1010"
            )
            self.assertNotEqual(catalog.fingerprint, before)


if __name__ == "__main__":
    unittest.main()
