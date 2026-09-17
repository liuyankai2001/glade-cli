import importlib
import importlib.util
import json
import tempfile
import unittest
from itertools import permutations
from pathlib import Path

from Bio import SeqIO
from plasmid_fixtures import make_project


class DesignServiceTests(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.find_spec("src.plasmid_design.service")
        self.assertIsNotNone(spec, "project design service must exist")
        self.Service = importlib.import_module(
            "src.plasmid_design.service"
        ).DesignService
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.config = make_project(Path(self.tmp.name))
        self.service = self.Service(self.config)

    def request(self):
        context = self.service.context()
        return {
            "resistance_id": "basic_seva_ap",
            "replication_id": "basic_seva_p15a",
            "t0_id": "basic_seva_t0",
            "t1_id": "basic_seva_t1",
            "expected_revision": context["project"]["manifest_revision"],
            "source_fingerprint": context["project"]["source_fingerprint"],
        }

    def test_preview_is_real_and_has_no_manifest_side_effect(self):
        before = self.config.manifest_output_path.read_bytes()
        context = self.service.context()
        self.assertTrue(context["ready"], context["issues"])
        self.assertEqual(context["construct"]["length_bp"], 31)
        preview = self.service.preview(self.request())
        self.assertTrue(preview["valid"], preview["issues"])
        self.assertEqual(preview["length_bp"], len(preview["sequence"]))
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)
        self.assertFalse((self.config.project_output_path / "plasmid_designs").exists())

    def test_missing_terminators_can_generate_and_record_absence_and_warnings(self):
        request = {**self.request(), "t0_id": None, "t1_id": None}
        preview = self.service.preview(request)
        self.assertTrue(preview["valid"], preview["issues"])
        self.assertEqual(len(preview["terminator_warnings"]), 2)
        result = self.service.generate(request)
        self.assertIsNone(result["t0_id"])
        self.assertIsNone(result["t1_id"])
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        component = manifest["plasmid_selection"]["component_design"]
        self.assertIsNone(component["t0_id"])
        self.assertIsNone(component["t1_id"])
        self.assertTrue(
            all(w in result["warnings"] for w in preview["terminator_warnings"])
        )
        final = SeqIO.read(
            self.service.download(result["id"], "final_genbank"), "genbank"
        )
        self.assertEqual(str(final.seq), preview["sequence"])
        self.assertFalse(
            any(
                f.qualifiers.get("label", [""])[0] in ("SEVA_T0", "SEVA_T1", "T0", "T1")
                for f in final.features
            )
        )

    def test_reordered_generation_matches_manifest_context_and_legacy_execution(self):
        from src.final_assemble_execute.common import (
            assemble_sequence,
            load_execution_context,
        )

        for order in permutations(["resistance", "replication", "expression"]):
            with self.subTest(order=order):
                expanded = [
                    piece
                    for kind in order
                    for piece in (
                        ["t1", "expression", "t0"] if kind == "expression" else [kind]
                    )
                ]
                request = {**self.request(), "component_order": expanded}
                preview = self.service.preview(request)
                expanded = [
                    piece
                    for i, kind in enumerate(expanded)
                    for piece in (
                        [kind, f"legacy-gap-{i}-{kind}"]
                        if kind in ("resistance", "replication")
                        else [kind]
                    )
                ]
                self.assertEqual(preview["component_order"], expanded)
                result = self.service.generate(request)
                self.assertEqual(result["component_order"], expanded)
                context = self.service.context()
                self.assertEqual(context["selection"]["component_order"], expanded)
                self.assertEqual(context["result"]["component_order"], expanded)
                manifest = json.loads(
                    self.config.manifest_output_path.read_text(encoding="utf-8")
                )
                self.assertEqual(
                    manifest["plasmid_selection"]["component_design"][
                        "component_order"
                    ],
                    expanded,
                )
                final = SeqIO.read(
                    self.service.download(result["id"], "final_genbank"), "genbank"
                )
                self.assertEqual(str(final.seq), preview["sequence"])
                execution = load_execution_context(self.config)
                plan = execution.plans[0]
                backbone = SeqIO.read(
                    self.service.download(result["id"], "backbone_genbank"), "genbank"
                )
                insert = SeqIO.read(
                    self.config.project_output_path / "expression_1.gb", "genbank"
                )
                self.assertEqual(
                    assemble_sequence(
                        str(backbone.seq), str(insert.seq), plan
                    ).sequence,
                    str(final.seq),
                )

    def test_generation_commits_existing_schemas_and_real_export_files(self):
        before_source = (
            self.config.project_output_path / "expression_1.gb"
        ).read_bytes()
        result = self.service.generate(self.request())
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["revision"], 2)
        for section, schema in [
            ("plasmid_selection", "plasmid_selection.v2"),
            ("final_assembly_plan", "final_assembly_plan.v2"),
            ("final_assembly", "final_assembly.v2"),
        ]:
            self.assertEqual(manifest[section]["schema_version"], schema)
        for file in result["files"]:
            path = self.service.download(result["id"], file["id"])
            self.assertTrue(path.is_file())
        final_file = next(f for f in result["files"] if f["id"] == "final_genbank")
        record = SeqIO.read(
            self.service.download(result["id"], final_file["id"]), "genbank"
        )
        self.assertEqual(len(record), result["length_bp"])
        self.assertEqual(record.annotations["topology"], "circular")
        self.assertEqual(
            (self.config.project_output_path / "expression_1.gb").read_bytes(),
            before_source,
        )
        # The existing CLI executor can authenticate a module-derived web plan.
        from src.final_assemble_execute.common import load_execution_context

        context = load_execution_context(self.config)
        self.assertEqual(len(context.plans), 1)
        self.assertEqual(context.plans[0]["restriction"]["left_enzyme"], "EcoRI")
        self.assertEqual(self.service.context()["result"]["id"], result["id"])

    def test_changed_revision_rejects_generation_and_preserves_old_results(self):
        result = self.service.generate(self.request())
        request = self.request()
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        manifest["revision"] += 1
        self.config.manifest_output_path.write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        before = self.config.manifest_output_path.read_bytes()
        with self.assertRaises(ValueError):
            self.service.generate(request)
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)
        self.assertTrue(self.service.download(result["id"], "final_genbank").is_file())

    def test_changed_source_or_tampered_download_is_rejected(self):
        result = self.service.generate(self.request())
        path = self.service.download(result["id"], "final_genbank")
        path.write_text("broken", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.service.download(result["id"], "final_genbank")
        source = self.config.project_output_path / "expression_1.gb"
        source.write_text(source.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        self.assertFalse(self.service.context()["ready"])
        with self.assertRaises(ValueError):
            self.service.preview(self.request())

    def test_multiple_confirmed_constructs_are_not_silently_selected(self):
        extra = Path(self.tmp.name) / "multi"
        config = make_project(extra, count=2)
        context = self.Service(config).context()
        self.assertFalse(context["ready"])
        self.assertTrue(any(i["code"] == "construct_count" for i in context["issues"]))
        self.assertIsNone(context["construct"])

    def test_real_export_failure_rolls_back_only_the_new_version(self):
        result = self.service.generate(self.request())
        old_file = self.service.download(result["id"], "final_genbank")
        before = self.config.manifest_output_path.read_bytes()
        parent = self.config.project_output_path / "plasmid_designs"
        initial_versions = set(parent.iterdir())

        def obstruct(stage, message):
            if stage == "exporting":
                newly_created = (set(parent.iterdir()) - initial_versions).pop()
                (newly_created / "backbone.gb").mkdir()

        with self.assertRaises(OSError):
            self.service.generate(self.request(), progress=obstruct)
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)
        self.assertEqual(set(parent.iterdir()), initial_versions)
        self.assertTrue(old_file.is_file())

    def test_upstream_change_while_exporting_never_commits_stale_result(self):
        result = self.service.generate(self.request())
        old_file = self.service.download(result["id"], "final_genbank")
        parent = self.config.project_output_path / "plasmid_designs"
        initial_versions = set(parent.iterdir())
        request = self.request()

        def update_source(stage, message):
            if stage == "registering":
                data = json.loads(
                    self.config.manifest_output_path.read_text(encoding="utf-8")
                )
                data["revision"] += 1
                data["cds_selection"]["restriction_enzymes"] = ["EcoRI", "BamHI"]
                self.config.manifest_output_path.write_text(
                    json.dumps(data), encoding="utf-8"
                )

        with self.assertRaises(ValueError):
            self.service.generate(request, progress=update_source)
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        self.assertEqual(
            manifest["plasmid_selection"]["component_design"]["generation_id"],
            result["id"],
        )
        self.assertEqual(set(parent.iterdir()), initial_versions)
        self.assertTrue(old_file.is_file())

    def test_missing_enzyme_configuration_keeps_construct_visible(self):
        data = json.loads(self.config.manifest_output_path.read_text(encoding="utf-8"))
        data["cds_selection"]["restriction_enzymes"] = []
        self.config.manifest_output_path.write_text(json.dumps(data), encoding="utf-8")
        context = self.service.context()
        self.assertFalse(context["ready"])
        self.assertIsNotNone(context["construct"])
        self.assertEqual(context["issues"][0]["code"], "enzyme_count")

    def test_malformed_manifest_sections_return_a_recoverable_context(self):
        original = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        for value in [None, [], "invalid"]:
            data = dict(original)
            data["cds_selection"] = value
            self.config.manifest_output_path.write_text(
                json.dumps(data), encoding="utf-8"
            )
            context = self.service.context()
            self.assertFalse(context["ready"])
            self.assertTrue(context["issues"])
        self.config.manifest_output_path.write_text(
            json.dumps(original), encoding="utf-8"
        )
        self.service.generate(self.request())
        data = json.loads(self.config.manifest_output_path.read_text(encoding="utf-8"))
        data["plasmid_selection"] = None
        self.config.manifest_output_path.write_text(json.dumps(data), encoding="utf-8")
        context = self.service.context()
        self.assertIsNone(context["result"])
        self.assertTrue(context["issues"])

    def test_bad_registered_path_is_reported_without_losing_source_context(self):
        self.service.generate(self.request())
        data = json.loads(self.config.manifest_output_path.read_text(encoding="utf-8"))
        selection = data["plasmid_selection"]
        selection["vector"]["selected_sequence_file"]["path"] = "../outside.gb"
        from src.plasmid_selection.get_plasmid_context import stable_json_hash

        selection.pop("selection_fingerprint")
        selection["selection_fingerprint"] = stable_json_hash(selection)
        final = data["final_assembly"]
        final["source_fingerprints"]["plasmid_selection_fingerprint"] = selection[
            "selection_fingerprint"
        ]
        final.pop("selection_fingerprint")
        final["selection_fingerprint"] = stable_json_hash(final)
        self.config.manifest_output_path.write_text(json.dumps(data), encoding="utf-8")
        context = self.service.context()
        self.assertTrue(context["ready"])
        self.assertIsNotNone(context["construct"])
        self.assertIsNone(context["result"])
        self.assertTrue(any(i["code"] == "artifact_invalid" for i in context["issues"]))


if __name__ == "__main__":
    unittest.main()
