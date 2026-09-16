import hashlib
import importlib
import importlib.util
import json
import tempfile
import time
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from fastapi.testclient import TestClient
from plasmid_fixtures import make_project


class LocalApiTests(unittest.TestCase):
    def setUp(self):
        try:
            spec = importlib.util.find_spec("src.web_api.app")
        except ModuleNotFoundError:
            spec = None
        self.assertIsNotNone(spec, "local project API must exist")
        create_app = importlib.import_module("src.web_api.app").create_app
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.config = make_project(Path(self.tmp.name))
        self.client = TestClient(create_app(self.config), base_url="http://127.0.0.1")
        self.client.__enter__()
        self.addCleanup(lambda: self.client.__exit__(None, None, None))

    def request(self):
        project = self.client.get("/api/context").json()["project"]
        return {
            "resistance_id": "basic_seva_ap",
            "replication_id": "basic_seva_p15a",
            "expected_revision": project["manifest_revision"],
            "source_fingerprint": project["source_fingerprint"],
        }

    def test_library_and_preview_are_real_and_generate_download_is_a_genbank_file(self):
        modules = self.client.get("/api/modules").json()
        self.assertEqual(len(modules["resistance"]), 7)
        self.assertEqual(len(modules["replication"]), 5)
        self.assertEqual(
            {m["role"]: m["length_bp"] for m in modules["terminator"]},
            {"t0": 103, "t1": 105},
        )
        self.assertIsInstance(modules["resistance"][0]["notes"], list)
        request = self.request()
        preview = self.client.post("/api/preview", json=request)
        self.assertEqual(preview.status_code, 200)
        self.assertTrue(preview.json()["valid"])
        self.assertEqual(len(preview.json()["terminator_warnings"]), 2)
        response = self.client.post("/api/generate", json=request)
        self.assertEqual(response.status_code, 202)
        job_id = response.json()["id"]
        for _ in range(200):
            job = self.client.get(f"/api/jobs/{job_id}").json()
            if job["status"] in ("succeeded", "failed"):
                break
            time.sleep(0.01)
        self.assertEqual(job["status"], "succeeded", job)
        result = job["result"]
        self.assertIsNone(result["t0_id"])
        self.assertIsNone(result["t1_id"])
        self.assertTrue(
            all(w in result["warnings"] for w in preview.json()["terminator_warnings"])
        )
        file = next(f for f in result["files"] if f["id"] == "final_genbank")
        downloaded = self.client.get(file["url"])
        self.assertEqual(downloaded.status_code, 200)
        self.assertTrue(downloaded.text.startswith("LOCUS"))
        self.assertIn("attachment", downloaded.headers["content-disposition"])
        self.assertEqual(
            self.client.post("/api/generate", json=request).status_code, 409
        )
        self.assertEqual(
            self.client.get("/api/context").json()["result"]["id"], result["id"]
        )

    def test_strict_requests_and_foreign_browser_origin_do_not_modify_project(self):
        before = self.config.manifest_output_path.read_bytes()
        request = self.request()
        request["resistance_id"] = "unknown"
        response = self.client.post("/api/preview", json=request)
        self.assertEqual(response.status_code, 422)
        self.assertEqual(response.json()["error"]["code"], "unknown_module")
        request = self.request()
        request["expected_revision"] = True
        self.assertEqual(
            self.client.post("/api/generate", json=request).status_code, 422
        )
        response = self.client.post(
            "/api/generate",
            json=self.request(),
            headers={"Origin": "https://foreign.example"},
        )
        self.assertEqual(response.status_code, 403)
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)

    def test_order_is_validated_and_used_by_preview(self):
        request = {
            **self.request(),
            "component_order": ["expression", "replication", "resistance"],
        }
        preview = self.client.post("/api/preview", json=request)
        self.assertEqual(preview.status_code, 200)
        self.assertEqual(
            preview.json()["component_order"],
            ["t1", "expression", "t0", "replication", "resistance"],
        )
        self.assertEqual(preview.json()["segments"][0]["kind"], "restriction")
        before = self.config.manifest_output_path.read_bytes()
        for order in (
            [],
            ["expression", "expression", "resistance"],
            ["resistance", "replication", "unknown"],
            None,
            "expression",
        ):
            for endpoint in ("preview", "generate"):
                with self.subTest(order=order, endpoint=endpoint):
                    response = self.client.post(
                        "/api/" + endpoint, json={**request, "component_order": order}
                    )
                    self.assertEqual(response.status_code, 422)
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)

    def component_request(self):
        project = self.client.get("/api/context").json()["project"]
        return {
            "expected_revision": project["manifest_revision"],
            "source_fingerprint": project["source_fingerprint"],
            "components": [
                {
                    "instance_id": "t1-a",
                    "component_type": "t1",
                    "module_id": "basic_seva_t1",
                },
                {"instance_id": "expression", "component_type": "expression"},
                {
                    "instance_id": "t0-a",
                    "component_type": "t0",
                    "module_id": "basic_seva_t0",
                },
                {
                    "instance_id": "t0-b",
                    "component_type": "t0",
                    "module_id": "basic_seva_t0",
                },
                {
                    "instance_id": "amp-a",
                    "component_type": "resistance",
                    "module_id": "basic_seva_ap",
                },
                {
                    "instance_id": "amp-b",
                    "component_type": "resistance",
                    "module_id": "basic_seva_ap",
                },
                {
                    "instance_id": "replication",
                    "component_type": "replication",
                    "module_id": "basic_seva_p15a",
                },
            ],
        }

    def test_repeated_components_survive_generation_download_and_context(self):
        request = self.component_request()
        response = self.client.post("/api/preview", json=request)
        self.assertEqual(response.status_code, 200, response.text)
        preview = response.json()
        self.assertTrue(preview["valid"], preview["issues"])
        self.assertEqual(preview["components"], request["components"])
        self.assertEqual(preview["terminator_warnings"], [])
        self.assertEqual(
            [
                s["instance_id"]
                for s in preview["segments"]
                if s["kind"] == "resistance"
            ],
            ["amp-a", "amp-b"],
        )
        modules = self.client.get("/api/modules").json()
        amp = next(m for m in modules["resistance"] if m["id"] == "basic_seva_ap")
        for segment in preview["segments"]:
            if segment["kind"] == "resistance":
                self.assertEqual(
                    preview["sequence"][segment["start_bp"] - 1 : segment["end_bp"]],
                    amp["sequence"],
                )
        generation = self.client.post("/api/generate", json=request)
        self.assertEqual(generation.status_code, 202, generation.text)
        for _ in range(200):
            job = self.client.get(f"/api/jobs/{generation.json()['id']}").json()
            if job["status"] in ("succeeded", "failed"):
                break
            time.sleep(0.01)
        self.assertEqual(job["status"], "succeeded", job)
        self.assertEqual(job["result"]["components"], request["components"])
        result = job["result"]
        final_path = self.client.app.state.design_service.download(
            result["id"], "final_genbank"
        )
        final = SeqIO.read(final_path, "genbank")
        self.assertEqual(str(final.seq), preview["sequence"])
        self.assertEqual(
            {
                f.qualifiers["web_instance_id"][0]
                for f in final.features
                if f.qualifiers.get("web_kind") == ["resistance"]
                and f.qualifiers.get("web_instance_id")
            },
            {"amp-a", "amp-b"},
        )
        context = self.client.get("/api/context").json()
        self.assertEqual(context["selection"]["components"], request["components"])
        self.assertEqual(context["result"]["components"], request["components"])
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        self.assertEqual(
            manifest["plasmid_selection"]["component_design"]["components"],
            request["components"],
        )
        self.assertEqual(
            len(manifest["plasmid_selection"]["vector"]["selection_markers"]), 2
        )
        report = self.client.app.state.design_service.download(
            result["id"], "report"
        ).read_text(encoding="utf-8")
        for instance in ("amp-a", "amp-b", "t0-a", "t0-b"):
            self.assertIn(instance, report)

    def test_long_instance_ids_roundtrip_in_all_genbank_exports(self):
        request = self.component_request()
        request["components"][2]["instance_id"] = "11223344-5566-7788-9900-aabbccddeeff"
        request["components"][4][
            "instance_id"
        ] = "resistance:11223344-5566-7788-9900-aabbccddeeff"
        request["components"][5]["instance_id"] = "long-instance-" + "a" * 106
        generation = self.client.post("/api/generate", json=request)
        self.assertEqual(generation.status_code, 202, generation.text)
        for _ in range(200):
            job = self.client.get(f"/api/jobs/{generation.json()['id']}").json()
            if job["status"] in ("succeeded", "failed"):
                break
            time.sleep(0.01)
        self.assertEqual(job["status"], "succeeded", job)
        expected = {
            c["instance_id"]
            for c in request["components"]
            if c["component_type"] != "expression"
        }
        service = self.client.app.state.design_service
        for file_id in (
            "backbone_genbank",
            "final_genbank",
            "backbone_preparation_genbank",
        ):
            with self.subTest(file_id=file_id):
                record = SeqIO.read(
                    service.download(job["result"]["id"], file_id), "genbank"
                )
                restored = {
                    "".join(f.qualifiers["web_instance_id"])
                    for f in record.features
                    if f.qualifiers.get("web_instance_id")
                }
                self.assertEqual(restored, expected)
                from src.plasmid_design.sequence import feature_payloads

                self.assertEqual(
                    {
                        f["instance_id"]
                        for f in feature_payloads(record)
                        if "instance_id" in f
                    },
                    expected,
                )

    def test_component_requests_reject_ambiguous_or_invalid_instances_without_writing(
        self,
    ):
        before = self.config.manifest_output_path.read_bytes()
        base = self.component_request()
        variants = [
            {**base, "resistance_id": "basic_seva_ap"},
            {**base, "component_order": ["resistance", "replication", "expression"]},
            {**base, "components": []},
            {**base, "components": base["components"] + [base["components"][0]]},
            {
                **base,
                "components": [
                    c for c in base["components"] if c["component_type"] != "resistance"
                ],
            },
            {
                **base,
                "components": [
                    c for c in base["components"] if c["component_type"] != "expression"
                ],
            },
            {
                **base,
                "components": base["components"]
                + [
                    {
                        "instance_id": "ori-2",
                        "component_type": "replication",
                        "module_id": "basic_seva_p15a",
                    }
                ],
            },
            {
                **base,
                "components": [
                    {**c, "module_id": "missing"} if c["instance_id"] == "amp-a" else c
                    for c in base["components"]
                ],
            },
            {
                **base,
                "components": [
                    (
                        {**c, "module_id": "basic_seva_t1"}
                        if c["instance_id"] == "t0-a"
                        else c
                    )
                    for c in base["components"]
                ],
            },
            {
                **base,
                "components": [
                    (
                        {**c, "module_id": "basic_seva_ap"}
                        if c["component_type"] == "expression"
                        else c
                    )
                    for c in base["components"]
                ],
            },
        ]
        for request in variants:
            for endpoint in ("preview", "generate"):
                with self.subTest(request=request, endpoint=endpoint):
                    self.assertEqual(
                        self.client.post("/api/" + endpoint, json=request).status_code,
                        422,
                    )
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)

    def test_terminator_warnings_do_not_allow_restriction_conflicts_or_invalid_roles(
        self,
    ):
        request = {**self.request(), "t0_id": "basic_seva_t1"}
        self.assertEqual(
            self.client.post("/api/preview", json=request).status_code, 422
        )
        request = {
            **self.request(),
            "t0_id": "basic_seva_t0",
            "t1_id": "basic_seva_t1",
            "component_order": ["resistance", "replication", "t0", "expression", "t1"],
        }
        preview = self.client.post("/api/preview", json=request)
        self.assertEqual(preview.status_code, 200)
        self.assertTrue(preview.json()["valid"])
        self.assertEqual(len(preview.json()["terminator_warnings"]), 2)
        request["t1_id"] = "unknown"
        self.assertEqual(
            self.client.post("/api/generate", json=request).status_code, 422
        )

    def test_unknown_artifact_and_api_path_never_fall_back_to_frontend_html(self):
        response = self.client.get("/api/files/" + "a" * 32 + "/manifest")
        self.assertNotEqual(response.status_code, 200)
        self.assertIn("error", response.json())
        response = self.client.get("/api/not-a-route")
        self.assertEqual(response.status_code, 404)
        self.assertIn("error", response.json())
        index = self.client.get("/")
        self.assertEqual(index.status_code, 200)
        self.assertIn('<div id="root"></div>', index.text)
        self.assertNotEqual(self.client.get("/design_manifest.json").status_code, 200)

    def test_internal_site_still_blocks_with_missing_terminators(self):
        path = self.config.project_output_path / "expression_1.gb"
        record = SeqIO.read(path, "genbank")
        record.seq = Seq("GAATTC" + str(record.seq))
        for feature in record.features:
            feature.location += 6
        SeqIO.write(record, path, "genbank")
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        construct = manifest["assembled_expression_constructs"]["constructs"][0]
        construct.update(
            length_bp=len(record),
            sequence_sha256=hashlib.sha256(str(record.seq).encode()).hexdigest(),
            file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )
        self.config.manifest_output_path.write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        before = self.config.manifest_output_path.read_bytes()
        request = self.request()
        preview = self.client.post("/api/preview", json=request)
        self.assertEqual(preview.status_code, 200)
        self.assertFalse(preview.json()["valid"])
        self.assertEqual(len(preview.json()["terminator_warnings"]), 2)
        self.assertTrue(
            any(
                issue["code"] == "restriction_conflict"
                for issue in preview.json()["issues"]
            )
        )
        self.assertEqual(
            self.client.post("/api/generate", json=request).status_code, 422
        )
        self.assertEqual(self.config.manifest_output_path.read_bytes(), before)
        self.assertFalse((self.config.project_output_path / "plasmid_designs").exists())


if __name__ == "__main__":
    unittest.main()
