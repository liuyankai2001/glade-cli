import importlib
import importlib.util
import tempfile
import time
import unittest
from pathlib import Path

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
        self.assertIsInstance(modules["resistance"][0]["notes"], list)
        request = self.request()
        preview = self.client.post("/api/preview", json=request)
        self.assertEqual(preview.status_code, 200)
        self.assertTrue(preview.json()["valid"])
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
        self.assertEqual(preview.json()["component_order"], request["component_order"])
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


if __name__ == "__main__":
    unittest.main()
