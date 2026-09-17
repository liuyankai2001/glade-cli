import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"


class ModuleCatalogContractTests(unittest.TestCase):
    def test_catalog_exposes_all_curated_modules_with_web_payloads(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(DATA_DIR)

        self.assertEqual(
            {module.id: module.length_bp for module in catalog.resistance.values()},
            {
                "basic_seva_ap": 1039,
                "basic_seva_km": 927,
                "basic_seva_cm": 783,
                "basic_seva_sm_sp": 989,
                "basic_seva_tet_5a": 1267,
                "basic_seva_gm": 805,
                "basic_seva_gm_11": 824,
            },
        )
        self.assertEqual(
            {module.id: module.length_bp for module in catalog.replication.values()},
            {
                "basic_seva_rsf1010": 3671,
                "basic_seva_p15a": 732,
                "basic_seva_psc101": 1460,
                "basic_seva_pbr322_rop": 1385,
                "basic_seva_psc101_pkd46_ts": 1537,
            },
        )
        module = catalog.get_resistance("basic_seva_ap")
        self.assertEqual(module.length_bp, 1039)
        self.assertEqual(
            hashlib.sha256(module.sequence.encode()).hexdigest(),
            "ab14bcddcec1db800438d7ff11882f3a970a286db2ef0e2e9d5b09284b5cefa2",
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

    def test_terminator_payloads_and_annotations_are_immutable(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(DATA_DIR)
        self.assertEqual(set(catalog.terminator), {"basic_seva_t0", "basic_seva_t1"})
        for module_id, role, length in (
            ("basic_seva_t0", "t0", 103),
            ("basic_seva_t1", "t1", 105),
        ):
            module = catalog.get_terminator(module_id)
            self.assertEqual(
                (module.type, module.role, module.length_bp),
                ("terminator", role, length),
            )
            self.assertEqual(module.name, role.upper())
            self.assertEqual(
                set(module.to_payload()),
                {
                    "id",
                    "name",
                    "type",
                    "sequence",
                    "notes",
                    "length_bp",
                    "gc_percent",
                    "role",
                },
            )
            self.assertEqual(module.to_payload()["role"], role)
            self.assertEqual(module.features[0]["parts"][0]["start"], 0)
            self.assertEqual(module.features[0]["parts"][0]["end"], length)
            self.assertEqual(
                module.features[0]["qualifiers"]["label"], ("SEVA_" + role.upper(),)
            )
            with self.assertRaises(FrozenInstanceError):
                module.sequence = "A"
            with self.assertRaises(TypeError):
                module.features[0]["parts"][0]["end"] = 1
        with self.assertRaises(TypeError):
            catalog.terminator["x"] = catalog.get_terminator("basic_seva_t0")
        with self.assertRaises(AttributeError):
            catalog.terminator = {}
        with self.assertRaises(KeyError):
            catalog.get_terminator("missing")

    def test_split_resistance_retains_exact_parent_dna_and_source_features(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(DATA_DIR)
        snapshot = json.loads(
            (ROOT / "src/plasmid_design/source_features.json").read_text(
                encoding="utf-8"
            )
        )
        expected_parents = {
            "basic_seva_ap": "e589d17e939c40824023f92436804adf95f353c91ccf30eb9780155b470bc57f",
            "basic_seva_km": "ba94804c214aff91d0ba37d9e986569d58fe514e4d2abe9732766ec66f92a026",
            "basic_seva_cm": "d828a3cb5cf7eee9e64292b7ca6ce0331ec8edcc9059bae0030eaf9247fae588",
            "basic_seva_sm_sp": "80ca647ebd5f4bb318d56c116f647bb80a63c08b9577af864bdfbc792ef34f9a",
            "basic_seva_tet_5a": "84c764ccf85a46acfd548abf4ff76ffc6f191bcf09e61b37f0523510faa37a4f",
            "basic_seva_gm": "1c23b4979add17b384977980684c801431a9763cab93fb9d7424d654715cdb33",
            "basic_seva_gm_11": "431fc598ac83caf7cf4c7942e331e93025cebaef797fd9833dc2ed2c13c9f2b8",
        }
        expected_extracts = {
            "basic_seva_ap": "ab14bcddcec1db800438d7ff11882f3a970a286db2ef0e2e9d5b09284b5cefa2",
            "basic_seva_km": "8eac5e46ff21210298cf843254f0ff98422ecdf598b55ee789a2d26e6f11b65e",
            "basic_seva_cm": "34afaa3076169e6b4d717fcfa9c697f1efd7a66525002b9aff0044adb87ebc2a",
            "basic_seva_sm_sp": "da24b066e4ccfb64865951ef3656879a7bbf3628291c1e7f9c01e9bccac26795",
            "basic_seva_tet_5a": "8879e9f7ef95caf830df263edce06760e5c250689d0167ad6b67bc1280c29eb3",
            "basic_seva_gm": "78c6afbf2e0fdf54f6aadbdb56c7c9629cac0a147eb530186cd3ffbbeabd64dc",
            "basic_seva_gm_11": "e08cf2b57a0e62727c88560376c2d5d063402cc42e64200631cc59804aca2436",
        }
        t0 = catalog.get_terminator("basic_seva_t0").sequence
        self.assertEqual(
            t0,
            "CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG",
        )
        for module_id, parent_hash in expected_parents.items():
            module = catalog.get_resistance(module_id)
            parent = t0 + snapshot[module_id]["gap_extraction"]["original_sequence"]
            self.assertEqual(hashlib.sha256(parent.encode()).hexdigest(), parent_hash)
            metadata = snapshot[module_id]
            self.assertEqual(metadata["derivation"]["parent_sha256"], parent_hash)
            self.assertEqual(metadata["derivation"]["removed_prefix_bp"], 103)
            self.assertEqual(
                metadata["provenance"]["end_1based"]
                - metadata["provenance"]["start_1based"]
                + 1,
                len(metadata["gap_extraction"]["original_sequence"]),
            )
            for feature in module.features:
                self.assertNotIn("SEVA_T0", feature["qualifiers"].get("label", ()))
                for part in feature["parts"]:
                    self.assertEqual(
                        hashlib.sha256(
                            module.sequence[part["start"] : part["end"]].encode()
                        ).hexdigest(),
                        expected_extracts[module_id],
                    )
        t1_source = snapshot["basic_seva_t1"]["provenance"]
        self.assertEqual(
            (t1_source["start_1based"], t1_source["end_1based"]), (4993, 5097)
        )
        self.assertEqual(
            t1_source["source_sha256"],
            "a4cdb1eb8917f65f3727ec54108b604c97caa9768f546b9cf2ed3d2d9fc29917",
        )

    def test_terminator_csv_is_three_columns_and_affects_fingerprint(self):
        from src.plasmid_design.catalog import ModuleCatalog

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
                shutil.copy2(DATA_DIR / name, staged / name)
            path = staged / "terminator_modules.csv"
            with path.open(encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(reader.fieldnames, ["id", "name", "sequence"])
                rows = list(reader)
            before = ModuleCatalog(staged).fingerprint
            rows[0]["name"] = "Renamed T0"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            catalog = ModuleCatalog(staged)
            self.assertEqual(catalog.get_terminator("basic_seva_t0").name, "Renamed T0")
            self.assertNotEqual(catalog.fingerprint, before)
            rows[0]["sequence"] = "A" + rows[0]["sequence"][1:]
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "source feature snapshot mismatch"):
                ModuleCatalog(staged)

    def test_terminator_role_and_source_provenance_are_validated(self):
        from src.plasmid_design.catalog import ModuleCatalog

        read_json = ModuleCatalog._read_json
        for field, invalid in (
            ("role", "other"),
            ("source_sha256", "broken"),
            ("end_1based", 5096),
            ("source_record", 123),
            ("source_key", []),
        ):
            metadata = read_json(ROOT / "src/plasmid_design/source_features.json")
            target = (
                metadata["basic_seva_t1"]
                if field == "role"
                else metadata["basic_seva_t1"]["provenance"]
            )
            target[field] = invalid

            def staged_read(path, metadata=metadata):
                return (
                    metadata if path.name == "source_features.json" else read_json(path)
                )

            with (
                self.subTest(field=field),
                patch.object(ModuleCatalog, "_read_json", side_effect=staged_read),
                self.assertRaises(ValueError),
            ):
                ModuleCatalog(DATA_DIR)

    def test_scaffold_source_hash_and_coordinates_are_validated(self):
        from src.plasmid_design.catalog import ModuleCatalog

        read_json = ModuleCatalog._read_json
        for field in ("source_sha256", "end_1based"):
            metadata = read_json(ROOT / "src/plasmid_design/scaffold.json")
            if field == "source_sha256":
                metadata["provenance"][field] = "0" * 64
            else:
                metadata["provenance"]["replication_to_t1"][field] -= 1

            def staged_read(path, metadata=metadata):
                return metadata if path.name == "scaffold.json" else read_json(path)

            with (
                self.subTest(field=field),
                patch.object(ModuleCatalog, "_read_json", side_effect=staged_read),
                self.assertRaises(ValueError),
            ):
                ModuleCatalog(DATA_DIR)

    def test_identifiers_are_unique_across_all_three_module_kinds(self):
        from src.plasmid_design.catalog import ModuleCatalog

        for collision in ("basic_seva_ap", "basic_seva_p15a"):
            with (
                self.subTest(collision=collision),
                tempfile.TemporaryDirectory() as tmp,
            ):
                staged = Path(tmp)
                for name in (
                    "resistance_modules.csv",
                    "replication_modules.csv",
                    "terminator_modules.csv",
                    "gap_modules.csv",
                ):
                    shutil.copy2(DATA_DIR / name, staged / name)
                path = staged / "terminator_modules.csv"
                with path.open(encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle))
                rows[0]["id"] = collision
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                    writer.writeheader()
                    writer.writerows(rows)
                with self.assertRaisesRegex(ValueError, "unique across types"):
                    ModuleCatalog(staged)

    def test_loading_verified_csv_dna_never_removes_an_arbitrary_prefix(self):
        from src.plasmid_design.catalog import ModuleCatalog

        read_json = ModuleCatalog._read_json
        metadata = read_json(ROOT / "src/plasmid_design/source_features.json")
        original = ModuleCatalog(DATA_DIR)
        prefix = original.get_terminator("basic_seva_t0").sequence
        sequence = prefix + original.get_resistance("basic_seva_ap").sequence
        item = metadata["basic_seva_ap"]
        item.pop("derivation")
        item["sha256"] = hashlib.sha256(sequence.encode()).hexdigest()
        item["provenance"]["start_1based"] -= 103
        for feature in item["features"]:
            for part in feature["parts"]:
                part["start"] += 103
                part["end"] += 103

        item["provenance"]["end_1based"] = (
            item["provenance"]["start_1based"] + len(sequence) - 1
        )
        item["gap_extraction"] = {
            "status": "no_confirmed_internal_gap",
            "original_sequence": sequence,
            "original_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "original_features": item["features"],
            "removed_gaps": [],
            "retained_segments": [
                {
                    "original_start_1based": 1,
                    "original_end_1based": len(sequence),
                    "module_start_1based": 1,
                    "module_end_1based": len(sequence),
                    "source_start_1based": item["provenance"]["start_1based"],
                    "source_end_1based": item["provenance"]["end_1based"],
                }
            ],
        }

        def staged_read(path):
            return metadata if path.name == "source_features.json" else read_json(path)

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
                shutil.copy2(DATA_DIR / name, staged / name)
            path = staged / "resistance_modules.csv"
            with path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["sequence"] = sequence
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            with patch.object(ModuleCatalog, "_read_json", side_effect=staged_read):
                self.assertEqual(
                    ModuleCatalog(staged).get_resistance("basic_seva_ap").sequence,
                    sequence,
                )

    def test_recorded_split_derivation_rejects_changed_parent_or_prefix_hash(self):
        from src.plasmid_design.catalog import ModuleCatalog

        read_json = ModuleCatalog._read_json
        for field in ("parent_sha256", "removed_prefix_sha256"):
            metadata = read_json(ROOT / "src/plasmid_design/source_features.json")
            metadata["basic_seva_ap"]["derivation"][field] = "0" * 64

            def staged_read(path, metadata=metadata):
                return (
                    metadata if path.name == "source_features.json" else read_json(path)
                )

            with (
                self.subTest(field=field),
                patch.object(ModuleCatalog, "_read_json", side_effect=staged_read),
                self.assertRaisesRegex(ValueError, "source derivation mismatch"),
            ):
                ModuleCatalog(DATA_DIR)

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
        self.assertNotIn("t1", scaffold)
        self.assertEqual(len(scaffold["landing_pad_spacer"]), 24)
        self.assertEqual(
            scaffold["landing_pad_spacer"], scaffold["resistance_to_replication"][12:36]
        )
        self.assertEqual(
            scaffold["provenance"]["resistance_to_replication"]["start_1based"], 1232
        )
        self.assertNotIn("t1", scaffold["provenance"])

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
        self.assertEqual(len(ap_features), 1)
        self.assertEqual(ap_features[0]["parts"][0]["start"], 0)
        self.assertEqual(ap_features[0]["parts"][0]["end"], 1039)
        pkd46_features = catalog.get_replication("basic_seva_psc101_pkd46_ts").features
        self.assertEqual(pkd46_features[0]["type"], "CDS")
        self.assertEqual(pkd46_features[0]["parts"][0]["start"], 61)
        self.assertEqual(pkd46_features[1]["parts"][0]["end"], 1282)

    def test_invalid_csv_and_unknown_ids_are_rejected(self):
        from src.plasmid_design.catalog import ModuleCatalog

        with tempfile.TemporaryDirectory() as tmp:
            staged = Path(tmp)
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
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
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
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
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
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
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
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
            for name in (
                "resistance_modules.csv",
                "replication_modules.csv",
                "terminator_modules.csv",
                "gap_modules.csv",
            ):
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
