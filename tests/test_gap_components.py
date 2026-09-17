import csv
import hashlib
import importlib
import json
import tempfile
import time
import unittest
from pathlib import Path

from Bio import SeqIO
from fastapi.testclient import TestClient
from plasmid_fixtures import make_project

ROOT = Path(__file__).resolve().parents[1]


class GapAssemblyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.config = make_project(Path(self.tmp.name))
        create_app = importlib.import_module("src.web_api.app").create_app
        self.client = TestClient(create_app(self.config), base_url="http://127.0.0.1")
        self.client.__enter__()
        self.addCleanup(lambda: self.client.__exit__(None, None, None))

    def selection(self, *, version=2):
        context = self.client.get("/api/context").json()
        return {
            "assembly_schema_version": version,
            "expected_revision": context["project"]["manifest_revision"],
            "source_fingerprint": context["project"]["source_fingerprint"],
            "components": [
                {
                    "instance_id": "amp",
                    "component_type": "resistance",
                    "module_id": "basic_seva_ap",
                },
                {
                    "instance_id": "ori",
                    "component_type": "replication",
                    "module_id": "basic_seva_p15a",
                },
                {"instance_id": "source", "component_type": "expression"},
            ],
        }

    def test_gapless_v2_and_legacy_intervals_have_distinct_exact_lengths(self):
        new = self.client.post("/api/preview", json=self.selection())
        self.assertEqual(new.status_code, 200, new.text)
        body = new.json()
        self.assertTrue(body["valid"], body["issues"])
        self.assertEqual(body["assembly_schema_version"], 2)
        self.assertEqual(body["length_bp"], 1039 + 732 + 31 + 12)
        self.assertFalse(any(s["kind"] in ("gap", "linker") for s in body["segments"]))
        old_request = self.selection(version=1)
        old = self.client.post("/api/preview", json=old_request)
        self.assertEqual(old.status_code, 200, old.text)
        self.assertEqual(old.json()["length_bp"], body["length_bp"] + 76 + 14)
        self.assertEqual(
            [c["component_type"] for c in old.json()["components"]],
            ["resistance", "gap", "replication", "gap", "expression"],
        )
        self.assertEqual(
            [
                c["instance_id"]
                for c in old.json()["components"]
                if c["component_type"] == "gap"
            ],
            ["legacy-gap-0-resistance", "legacy-gap-1-replication"],
        )
        restored = {**self.selection(), "components": old.json()["components"]}
        self.assertEqual(
            self.client.post("/api/preview", json=restored).json()["sequence"],
            old.json()["sequence"],
        )

    def test_repeated_gap_dna_survives_generate_annotations_context_and_report(self):
        modules = self.client.get("/api/modules").json()
        self.assertIn("gap", modules)
        gap = next(g for g in modules["gap"] if "replication_to_t1" in g["aliases"])
        request = self.selection()
        copies = [
            {"instance_id": f"gap-{i}", "component_type": "gap", "module_id": gap["id"]}
            for i in range(3)
        ]
        request["components"][1:1] = copies
        response = self.client.post("/api/preview", json=request)
        self.assertEqual(response.status_code, 200, response.text)
        preview = response.json()
        self.assertTrue(preview["valid"], preview["issues"])
        self.assertEqual(preview["length_bp"], 1039 + 732 + 31 + 12 + 3 * 14)
        for segment in (s for s in preview["segments"] if s["kind"] == "gap"):
            self.assertEqual(
                preview["sequence"][segment["start_bp"] - 1 : segment["end_bp"]],
                "GGCGCGCCCAGCTG",
            )
        result = self.client.post("/api/generate", json=request)
        self.assertEqual(result.status_code, 202, result.text)
        for _ in range(200):
            job = self.client.get(f"/api/jobs/{result.json()['id']}").json()
            if job["status"] in ("failed", "succeeded"):
                break
            time.sleep(0.01)
        self.assertEqual(job["status"], "succeeded", job)
        context = self.client.get("/api/context").json()
        self.assertEqual(context["selection"]["assembly_schema_version"], 2)
        self.assertEqual(context["result"]["components"], request["components"])
        service = self.client.app.state.design_service
        final = SeqIO.read(
            service.download(job["result"]["id"], "final_genbank"), "genbank"
        )
        self.assertEqual(str(final.seq), preview["sequence"])
        self.assertEqual(
            {
                "".join(f.qualifiers["web_instance_id"])
                for f in final.features
                if f.qualifiers.get("web_kind") == ["gap"]
            },
            {"gap-0", "gap-1", "gap-2"},
        )
        report = service.download(job["result"]["id"], "report").read_text(
            encoding="utf-8"
        )
        self.assertIn(gap["id"], report)

    def test_unknown_gap_and_layout_version_rejected_without_writing(self):
        before = self.config.manifest_output_path.read_bytes()
        request = self.selection()
        request["components"].append(
            {"instance_id": "gap", "component_type": "gap", "module_id": "missing"}
        )
        self.assertEqual(
            self.client.post("/api/preview", json=request).status_code, 422
        )
        for version in (0, 3, True, "2", None):
            self.assertEqual(
                self.client.post(
                    "/api/preview",
                    json={**self.selection(), "assembly_schema_version": version},
                ).status_code,
                422,
            )
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)

    def test_legacy_migration_handles_id_collisions_and_explicit_unversioned_gaps(self):
        request = self.selection(version=1)
        request["components"][2]["instance_id"] = "legacy-gap-0-resistance"
        response = self.client.post("/api/preview", json=request)
        self.assertEqual(response.status_code, 200, response.text)
        migrated = response.json()["components"]
        self.assertEqual(migrated[1]["instance_id"], "legacy-gap-0-resistance-1")
        explicit = {**request, "components": migrated}
        explicit.pop("assembly_schema_version")
        restored = self.client.post("/api/preview", json=explicit)
        self.assertEqual(restored.json()["components"], migrated)
        self.assertEqual(restored.json()["sequence"], response.json()["sequence"])

    def test_context_migrates_old_selection_without_writes_and_preserves_v2_gaplessness(
        self,
    ):
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        components = self.selection()["components"]
        manifest["plasmid_selection"] = {"component_design": {"components": components}}
        self.config.manifest_output_path.write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        before = self.config.manifest_output_path.read_bytes()
        selection = self.client.get("/api/context").json()["selection"]
        self.assertEqual(selection["assembly_schema_version"], 2)
        self.assertEqual(len(selection["components"]), 5)
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)
        manifest["plasmid_selection"]["component_design"]["assembly_schema_version"] = 2
        self.config.manifest_output_path.write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        selection = self.client.get("/api/context").json()["selection"]
        self.assertEqual(selection["components"], components)


class GapLibraryTests(unittest.TestCase):
    def test_archived_pinned_sources_match_every_gap_and_original_module(self):
        import ast
        import io

        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(ROOT / "data")
        metadata = json.loads(
            (ROOT / "src/plasmid_design/source_features.json").read_text(
                encoding="utf-8"
            )
        )
        records = {}
        for key, source in metadata["_sources"].items():
            archive = (
                ROOT
                / "data/gap_sources"
                / (key + (".py.txt" if key == "standard_linkers" else ".gb"))
            )
            raw = archive.read_bytes()
            self.assertEqual(hashlib.sha256(raw).hexdigest(), source["source_sha256"])
            if key == "standard_linkers":
                tree = ast.parse(raw.decode())
                standard = next(
                    ast.literal_eval(node.value)
                    for node in tree.body
                    if isinstance(node, ast.Assign)
                    and any(
                        isinstance(target, ast.Name) and target.id == "STANDARD_LINKERS"
                        for target in node.targets
                    )
                )
            else:
                records[key] = {
                    record.name: record
                    for record in SeqIO.parse(io.StringIO(raw.decode()), "genbank")
                }
        for gap in catalog.gap.values():
            for source in metadata[gap.id]["sources"]:
                sequence = (
                    standard[source["source_record"].split(".")[-1]]
                    if source["source_key"] == "standard_linkers"
                    else str(
                        records[source["source_key"]][source["source_record"]].seq[
                            source["start_1based"] - 1 : source["end_1based"]
                        ]
                    ).upper()
                )
                self.assertEqual(sequence, gap.sequence)
        for module in [
            *catalog.resistance.values(),
            *catalog.replication.values(),
            *catalog.terminator.values(),
        ]:
            item = metadata[module.id]
            source = item["provenance"]
            original = item["gap_extraction"]["original_sequence"]
            self.assertEqual(
                str(
                    records[source["source_key"]][source["source_record"]].seq[
                        source["start_1based"] - 1 : source["end_1based"]
                    ]
                ).upper(),
                original,
            )
            for feature, remapped in zip(
                item["gap_extraction"]["original_features"],
                item["features"],
                strict=True,
            ):
                for part, mapped in zip(
                    feature["parts"], remapped["parts"], strict=True
                ):
                    self.assertEqual(
                        original[part["start"] : part["end"]],
                        module.sequence[mapped["start"] : mapped["end"]],
                    )

    def test_gap_library_is_unique_source_traced_and_covers_all_modules(self):
        from src.plasmid_design.catalog import ModuleCatalog

        catalog = ModuleCatalog(ROOT / "data")
        self.assertEqual(len(catalog.gap), 12)
        aliases = {alias for gap in catalog.gap.values() for alias in gap.aliases}
        self.assertTrue(
            {
                "L1",
                "L2",
                "L3",
                "L4",
                "L5",
                "L6",
                "BSEVA_L1",
                "resistance_to_replication",
                "replication_to_t1",
                "landing_pad_spacer",
                "resistance_prefix_26",
                "resistance_prefix_34",
            }.issubset(aliases)
        )
        sequences = [g.sequence for g in catalog.gap.values()]
        self.assertEqual(len(sequences), len(set(sequences)))
        metadata = json.loads(
            (ROOT / "src/plasmid_design/source_features.json").read_text(
                encoding="utf-8"
            )
        )
        for gap in catalog.gap.values():
            self.assertEqual(
                gap.id, "gap_" + hashlib.sha256(gap.sequence.encode()).hexdigest()
            )
            self.assertTrue(gap.to_payload()["aliases"])
            self.assertTrue(metadata[gap.id]["sources"])
        for module in [
            *catalog.resistance.values(),
            *catalog.replication.values(),
            *catalog.terminator.values(),
        ]:
            extraction = metadata[module.id]["gap_extraction"]
            self.assertIn(
                extraction["status"],
                ("confirmed_removed", "no_confirmed_internal_gap", "pending_review"),
            )
            self.assertEqual(
                hashlib.sha256(extraction["original_sequence"].encode()).hexdigest(),
                extraction["original_sha256"],
            )
