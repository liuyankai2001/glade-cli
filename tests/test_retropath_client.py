from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable, Optional

import httpx

from src.pathway_analyze.retropath_client import (
    RetroPathClientError,
    RetroPathHttpClient,
    RetroPathJobParameters,
)
from src.pathway_analyze.retropath_input import RetroPathInputBundle
from src.pathway_analyze.retropath_models import PredictedCompound


RULES_SHA256 = "a" * 64
DEFAULT_JOB_ID = f"rp2-{'1' * 32}"


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def make_input_bundle(root: Path) -> RetroPathInputBundle:
    input_dir = root / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    source = (
        "Name,InChI\n"
        'C00999,"InChI=1S/CH4/h1H4"\n'
    ).encode()
    sink = (
        "Name,InChI\n"
        'C00001,"InChI=1S/H2O/h1H2"\n'
    ).encode()
    source_path = input_dir / "target_source.csv"
    sink_path = input_dir / "chassis_sink.csv"
    source_path.write_bytes(source)
    sink_path.write_bytes(sink)
    target = PredictedCompound.create(
        compound_id="C00999",
        inchi="InChI=1S/CH4/h1H4",
        inchikey="VNWKTOKETHGBQD-UHFFFAOYSA-N",
        isomeric_smiles="C",
        formula="CH4",
        charge=0,
        kegg_ids=("C00999",),
        structure_provenance=("test",),
    )
    water = PredictedCompound.create(
        compound_id="C00001",
        inchi="InChI=1S/H2O/h1H2",
        inchikey="XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        isomeric_smiles="O",
        formula="H2O",
        charge=0,
        kegg_ids=("C00001",),
        minimum_depth=0,
        structure_provenance=("test",),
    )
    return RetroPathInputBundle(
        expansion_depth=0,
        reachable_compound_count=1,
        target_compound=target,
        sink_compounds=(water,),
        mappings=tuple(),
        rejected_compounds=tuple(),
        target_source_path=source_path,
        chassis_sink_path=sink_path,
        compound_mapping_path=input_dir / "compound_mapping.csv",
        rejected_compounds_path=input_dir / "rejected_compounds.csv",
        target_source_sha256=sha256_bytes(source),
        chassis_sink_sha256=sha256_bytes(sink),
    )


class MockRetroPathService:
    def __init__(self, bundle: RetroPathInputBundle) -> None:
        self.bundle = bundle
        self.parameters = RetroPathJobParameters()
        self.ready = True
        self.rules_version = "rr02-rp2-hs"
        self.rules_sha256 = RULES_SHA256
        self.health_statuses: list[int] = []
        self.status_sequence = ["queued", "running", "source_in_sink"]
        self.pending_statuses: list[str] = []
        self.final_status = "source_in_sink"
        self.post_status: Optional[int] = None
        self.post_detail = "invalid input"
        self.raise_post_once = False
        self.results_status: Optional[int] = None
        self.bad_json_operation: Optional[str] = None
        self.job_id_mismatch = False
        self.manifest_mutator: Optional[Callable[[dict[str, Any]], None]] = None
        self.artifact_paths = ["stdout.log", "stderr.log", "raw/results.csv"]
        self.artifact_bytes = {
            "stdout.log": b"KNIME-Worker completed\n",
            "stderr.log": b"",
            "raw/results.csv": b"transformation_id,rule_id\n",
        }
        self.artifact_statuses: dict[str, list[int]] = {}
        self.health_calls = 0
        self.submissions = 0
        self.job_get_calls = 0
        self.results_calls = 0
        self.artifact_calls: dict[str, int] = {}
        self.last_post_body = b""
        self.current_job_id = DEFAULT_JOB_ID

    def health_payload(self) -> dict[str, Any]:
        return {
            "ready": self.ready,
            "wrapper_version": "3.9.1",
            "wrapper_reported_version": "3.9.0",
            "workflow_version": "r20260212",
            "knime_version": "4.7.0",
            "rdkit_plugin_version": "4.9.1",
            "rules_version": self.rules_version,
            "rules_sha256": self.rules_sha256,
            "worker_concurrency": 1,
            "queue_active": 0,
            "errors": [] if self.ready else ["runtime missing"],
            "service_version": "1.0.0",
        }

    @staticmethod
    def return_code(status: str) -> Optional[int]:
        return {
            "queued": None,
            "running": None,
            "succeeded": 0,
            "source_in_sink": 10,
            "no_solution": 11,
            "failed": 1,
            "timed_out": -15,
        }[status]

    @staticmethod
    def error(status: str) -> Optional[str]:
        if status == "failed":
            return "retropath2_wrapper exited with code 1"
        if status == "timed_out":
            return "RetroPath execution exceeded 3600 seconds"
        return None

    def job_payload(
        self,
        status: str,
        *,
        job_id: Optional[str] = None,
    ) -> dict[str, Any]:
        return {
            "job_id": job_id or self.current_job_id,
            "status": status,
            "created_at": "2026-08-25T01:00:00+00:00",
            "started_at": (
                None if status == "queued" else "2026-08-25T01:00:01+00:00"
            ),
            "finished_at": (
                "2026-08-25T01:00:02+00:00"
                if status
                in {"succeeded", "source_in_sink", "no_solution", "failed", "timed_out"}
                else None
            ),
            "return_code": self.return_code(status),
            "error": self.error(status),
            "parameters": self.parameters.to_dict(),
        }

    def manifest(self) -> dict[str, Any]:
        manifest = {
            "schema_version": 1,
            "job_id": self.current_job_id,
            "status": self.final_status,
            "return_code": self.return_code(self.final_status),
            "parameters": self.parameters.to_dict(),
            "versions": {
                "retropath2_wrapper": "3.9.1",
                "retropath2_wrapper_reported": "3.9.0",
                "workflow": "r20260212",
                "knime": "4.7.0",
                "knime_rdkit_nodes": "4.9.1",
                "rules": self.rules_version,
            },
            "rules_sha256": self.rules_sha256,
            "input_sha256": {
                "source.csv": self.bundle.target_source_sha256,
                "sink.csv": self.bundle.chassis_sink_sha256,
            },
            "artifacts": list(self.artifact_paths),
        }
        if self.manifest_mutator is not None:
            self.manifest_mutator(manifest)
        return manifest

    def prepare_next_job(
        self,
        *,
        status: Optional[str] = None,
        sequence: Optional[list[str]] = None,
    ) -> None:
        if status is not None:
            self.final_status = status
        if sequence is not None:
            self.status_sequence = list(sequence)
        else:
            self.status_sequence = ["queued", self.final_status]
        self.pending_statuses = []

    def __call__(self, request: httpx.Request) -> httpx.Response:
        path = request.url.path
        if path == "/health" and request.method == "GET":
            self.health_calls += 1
            if self.health_statuses:
                status = self.health_statuses.pop(0)
                if status != 200:
                    return httpx.Response(status, json={"detail": "temporarily unavailable"})
            if self.bad_json_operation == "health":
                return httpx.Response(200, content=b"{")
            return httpx.Response(200, json=self.health_payload())

        if path == "/v1/jobs" and request.method == "POST":
            self.submissions += 1
            self.last_post_body = request.read()
            if self.raise_post_once:
                self.raise_post_once = False
                raise httpx.ReadError("response lost", request=request)
            if self.post_status is not None:
                return httpx.Response(
                    self.post_status,
                    json={"detail": self.post_detail},
                )
            self.current_job_id = f"rp2-{self.submissions:032x}"
            self.pending_statuses = list(self.status_sequence)
            return httpx.Response(202, json=self.job_payload("queued"))

        prefix = f"/v1/jobs/{self.current_job_id}"
        if path == prefix and request.method == "GET":
            self.job_get_calls += 1
            status = (
                self.pending_statuses.pop(0)
                if self.pending_statuses
                else self.final_status
            )
            if self.bad_json_operation == "job":
                return httpx.Response(200, content=b"not-json")
            if status == "unknown":
                payload = self.job_payload("queued")
                payload["status"] = "mystery"
                return httpx.Response(200, json=payload)
            job_id = f"rp2-{'f' * 32}" if self.job_id_mismatch else None
            return httpx.Response(200, json=self.job_payload(status, job_id=job_id))

        if path == f"{prefix}/results" and request.method == "GET":
            self.results_calls += 1
            if self.results_status is not None:
                return httpx.Response(
                    self.results_status,
                    json={"detail": "job has not reached a terminal state"},
                )
            payload = {
                **self.job_payload(self.final_status),
                "manifest": self.manifest(),
                "artifacts": ["run_manifest.json", *self.artifact_paths],
            }
            if self.bad_json_operation == "results":
                return httpx.Response(200, content=b"[")
            return httpx.Response(200, json=payload)

        artifact_prefix = f"{prefix}/artifacts/"
        if path.startswith(artifact_prefix) and request.method == "GET":
            remote_path = path[len(artifact_prefix) :]
            self.artifact_calls[remote_path] = (
                self.artifact_calls.get(remote_path, 0) + 1
            )
            statuses = self.artifact_statuses.get(remote_path, [])
            if statuses:
                status = statuses.pop(0)
                if status != 200:
                    return httpx.Response(status, json={"detail": "artifact unavailable"})
            if remote_path == "run_manifest.json":
                return httpx.Response(
                    200,
                    content=(
                        json.dumps(self.manifest(), sort_keys=True) + "\n"
                    ).encode(),
                )
            if remote_path not in self.artifact_bytes:
                return httpx.Response(404, json={"detail": "artifact not found"})
            return httpx.Response(200, content=self.artifact_bytes[remote_path])

        return httpx.Response(404, json={"detail": "job not found"})


class FakeClock:
    def __init__(self) -> None:
        self.value = 0.0

    def monotonic(self) -> float:
        return self.value

    def sleep(self, seconds: float) -> None:
        self.value += seconds


class RetroPathClientTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.bundle = make_input_bundle(self.root / "retropath")
        self.http_clients: list[httpx.Client] = []

    def tearDown(self) -> None:
        for client in self.http_clients:
            client.close()
        self.temporary_directory.cleanup()

    def make_client(
        self,
        service: MockRetroPathService,
        **kwargs: Any,
    ) -> RetroPathHttpClient:
        http_client = httpx.Client(transport=httpx.MockTransport(service))
        self.http_clients.append(http_client)
        return RetroPathHttpClient(
            client=http_client,
            sleep=kwargs.pop("sleep", lambda _: None),
            **kwargs,
        )

    def test_end_to_end_download_and_health_consistent_cache(self) -> None:
        service = MockRetroPathService(self.bundle)
        client = self.make_client(service)
        statuses: list[str] = []
        output_dir = self.root / "retropath"

        first = client.run(
            self.bundle,
            output_dir,
            on_status=lambda state: statuses.append(state.status),
        )
        submissions_after_first = service.submissions
        second = client.run(self.bundle, output_dir)

        self.assertEqual(first.result.status, "source_in_sink")
        self.assertEqual(first.result.return_code, 10)
        self.assertFalse(first.cache_hit)
        self.assertTrue(second.cache_hit)
        self.assertEqual(service.submissions, submissions_after_first)
        self.assertGreaterEqual(service.health_calls, 2)
        self.assertEqual(statuses, ["queued", "running", "source_in_sink"])
        self.assertIn(b'name="source_file"', service.last_post_body)
        self.assertIn(b"InChI=1S/CH4/h1H4", service.last_post_body)
        self.assertTrue((first.raw_dir / "service_run_manifest.json").is_file())
        self.assertTrue((first.raw_dir / "results.csv").is_file())
        manifest = json.loads(first.run_manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["run_result"]["status"], "source_in_sink")
        self.assertEqual(
            manifest["health"]["rules_sha256"],
            RULES_SHA256,
        )

    def test_all_service_terminal_statuses_are_returned(self) -> None:
        cases = {
            "succeeded": 0,
            "no_solution": 11,
            "failed": 1,
            "timed_out": -15,
        }
        for index, (status, return_code) in enumerate(cases.items()):
            with self.subTest(status=status):
                service = MockRetroPathService(self.bundle)
                service.prepare_next_job(status=status)
                client = self.make_client(service)
                result = client.run(
                    self.bundle,
                    self.root / f"terminal-{index}",
                )
                self.assertEqual(result.result.status, status)
                self.assertEqual(result.result.return_code, return_code)
                if status in {"failed", "timed_out"}:
                    self.assertTrue(result.result.errors)

    def test_failed_and_timed_out_jobs_are_not_reused(self) -> None:
        for index, status in enumerate(("failed", "timed_out")):
            with self.subTest(status=status):
                service = MockRetroPathService(self.bundle)
                service.prepare_next_job(status=status)
                client = self.make_client(service)
                output_dir = self.root / f"non-cacheable-{index}"

                first = client.run(self.bundle, output_dir)
                service.prepare_next_job(status=status)
                second = client.run(self.bundle, output_dir)

                self.assertFalse(first.cache_hit)
                self.assertFalse(second.cache_hit)
                self.assertEqual(service.submissions, 2)

    def test_modified_p2_input_is_rejected_before_health(self) -> None:
        service = MockRetroPathService(self.bundle)
        client = self.make_client(service)
        self.bundle.target_source_path.write_text(
            "Name,InChI\nchanged,InChI=1S/H2O/h1H2\n",
            encoding="utf-8",
        )

        with self.assertRaises(RetroPathClientError) as caught:
            client.run(self.bundle, self.root / "modified")

        self.assertEqual(caught.exception.code, "input_modified")
        self.assertEqual(service.health_calls, 0)

    def test_health_get_retries_and_runtime_not_ready_are_distinct(self) -> None:
        delays: list[float] = []
        service = MockRetroPathService(self.bundle)
        service.health_statuses = [503, 503, 200]
        client = self.make_client(service, sleep=delays.append)

        health = client.health()

        self.assertEqual(health.provenance.wrapper_version, "3.9.1")
        self.assertEqual(service.health_calls, 3)
        self.assertEqual(delays, [0.5, 1.0])

        unavailable = MockRetroPathService(self.bundle)
        unavailable.ready = False
        client = self.make_client(unavailable)
        with self.assertRaises(RetroPathClientError) as caught:
            client.health()
        self.assertEqual(caught.exception.code, "runtime_not_ready")

        down = MockRetroPathService(self.bundle)
        down.health_statuses = [502, 502, 502]
        client = self.make_client(down)
        with self.assertRaises(RetroPathClientError) as caught:
            client.health()
        self.assertEqual(caught.exception.code, "service_unavailable")

    def test_post_is_not_retried_and_uncertain_state_requires_force(self) -> None:
        service = MockRetroPathService(self.bundle)
        service.raise_post_once = True
        client = self.make_client(service)
        output_dir = self.root / "uncertain"

        with self.assertRaises(RetroPathClientError) as first_error:
            client.run(self.bundle, output_dir)
        with self.assertRaises(RetroPathClientError) as second_error:
            client.run(self.bundle, output_dir)
        forced = client.run(self.bundle, output_dir, force=True)

        self.assertEqual(first_error.exception.code, "submission_uncertain")
        self.assertEqual(second_error.exception.code, "submission_uncertain")
        self.assertEqual(service.submissions, 2)
        self.assertEqual(forced.result.status, "source_in_sink")

    def test_poll_timeout_preserves_job_for_resume(self) -> None:
        clock = FakeClock()
        service = MockRetroPathService(self.bundle)
        service.status_sequence = ["running"] * 20
        service.final_status = "running"
        client = self.make_client(
            service,
            sleep=clock.sleep,
            monotonic=clock.monotonic,
            wait_timeout_seconds=2,
            poll_interval_seconds=1,
        )
        output_dir = self.root / "resume"

        with self.assertRaises(RetroPathClientError) as caught:
            client.run(self.bundle, output_dir)
        submissions = service.submissions
        state = json.loads(
            (output_dir / "client_state.json").read_text(encoding="utf-8")
        )
        service.final_status = "source_in_sink"
        service.pending_statuses = ["source_in_sink"]
        resumed = client.run(self.bundle, output_dir)

        self.assertEqual(caught.exception.code, "client_poll_timeout")
        self.assertEqual(state["phase"], "client_poll_timeout")
        self.assertIsNotNone(state["job"])
        self.assertEqual(state["job"]["status"], "running")
        self.assertEqual(service.submissions, submissions)
        self.assertEqual(resumed.result.status, "source_in_sink")

    def test_cache_invalidates_for_artifact_parameters_and_runtime(self) -> None:
        service = MockRetroPathService(self.bundle)
        client = self.make_client(service)
        output_dir = self.root / "cache-invalidations"
        first = client.run(self.bundle, output_dir)
        first_artifact = output_dir / first.result.artifacts[0]
        first_artifact.write_bytes(b"corrupt")

        service.prepare_next_job()
        repaired = client.run(self.bundle, output_dir)
        after_corruption = service.submissions

        new_parameters = RetroPathJobParameters(topx=2)
        service.parameters = new_parameters
        service.prepare_next_job()
        client.run(self.bundle, output_dir, parameters=new_parameters)
        after_parameters = service.submissions

        service.rules_version = "rr02-rp2-hs-updated"
        service.rules_sha256 = "b" * 64
        service.prepare_next_job()
        client.run(self.bundle, output_dir, parameters=new_parameters)

        self.assertFalse(repaired.cache_hit)
        self.assertNotEqual(first_artifact.read_bytes(), b"corrupt")
        self.assertEqual(after_corruption, 1)
        self.assertEqual(after_parameters, 2)
        self.assertEqual(service.submissions, 3)

    def test_manifest_input_or_version_mismatch_is_protocol_error(self) -> None:
        mutators = (
            lambda manifest: manifest["input_sha256"].update(
                {"source.csv": "0" * 64}
            ),
            lambda manifest: manifest["versions"].update(
                {"workflow": "unexpected"}
            ),
            lambda manifest: manifest["parameters"].update({"topx": 999}),
            lambda manifest: manifest.update({"schema_version": 99}),
        )
        for index, mutator in enumerate(mutators):
            with self.subTest(index=index):
                service = MockRetroPathService(self.bundle)
                service.manifest_mutator = mutator
                client = self.make_client(service)
                with self.assertRaises(RetroPathClientError) as caught:
                    client.run(self.bundle, self.root / f"manifest-{index}")
                self.assertEqual(caught.exception.code, "protocol_error")

    def test_unsafe_artifact_path_and_download_failure_are_rejected(self) -> None:
        unsafe = MockRetroPathService(self.bundle)
        unsafe.artifact_paths = ["../escape.csv"]
        client = self.make_client(unsafe)
        with self.assertRaises(RetroPathClientError) as caught:
            client.run(self.bundle, self.root / "unsafe")
        self.assertEqual(caught.exception.code, "artifact_path_unsafe")

        missing = MockRetroPathService(self.bundle)
        missing.artifact_paths = ["raw/missing.csv"]
        client = self.make_client(missing)
        with self.assertRaises(RetroPathClientError) as caught:
            client.run(self.bundle, self.root / "missing")
        self.assertEqual(caught.exception.code, "artifact_download_failed")

    def test_artifact_transient_failure_is_retried(self) -> None:
        service = MockRetroPathService(self.bundle)
        service.artifact_statuses["raw/results.csv"] = [503, 200]
        client = self.make_client(service)

        result = client.run(self.bundle, self.root / "artifact-retry")

        self.assertEqual(result.result.status, "source_in_sink")
        self.assertEqual(service.artifact_calls["raw/results.csv"], 2)

    def test_bad_json_unknown_status_and_job_id_mismatch_are_protocol_errors(self) -> None:
        cases = ("health", "job", "job_id", "results")
        for index, case in enumerate(cases):
            with self.subTest(case=case):
                service = MockRetroPathService(self.bundle)
                if case == "health":
                    service.bad_json_operation = "health"
                elif case == "job":
                    service.status_sequence = ["unknown"]
                elif case == "results":
                    service.bad_json_operation = "results"
                else:
                    service.job_id_mismatch = True
                client = self.make_client(service)
                with self.assertRaises(RetroPathClientError) as caught:
                    client.run(self.bundle, self.root / f"protocol-{index}")
                self.assertEqual(caught.exception.code, "protocol_error")

    def test_http_error_mapping_for_submit_and_results(self) -> None:
        cases = (
            (400, "bad csv", "input_invalid"),
            (503, "RetroPath queue already contains 8 active jobs", "queue_full"),
        )
        for index, (status, detail, code) in enumerate(cases):
            with self.subTest(status=status):
                service = MockRetroPathService(self.bundle)
                service.post_status = status
                service.post_detail = detail
                client = self.make_client(service)
                with self.assertRaises(RetroPathClientError) as caught:
                    client.run(self.bundle, self.root / f"submit-error-{index}")
                self.assertEqual(caught.exception.code, code)

        service = MockRetroPathService(self.bundle)
        service.results_status = 409
        client = self.make_client(service)
        with self.assertRaises(RetroPathClientError) as caught:
            client.run(self.bundle, self.root / "results-not-ready")
        self.assertEqual(caught.exception.code, "protocol_error")

        service = MockRetroPathService(self.bundle)
        client = self.make_client(service)
        with self.assertRaises(RetroPathClientError) as caught:
            client.get_job(f"rp2-{'e' * 32}")
        self.assertEqual(caught.exception.code, "job_not_found")

    def test_parameter_and_loopback_validation(self) -> None:
        invalid_parameter_factories = (
            lambda: RetroPathJobParameters(max_steps=0),
            lambda: RetroPathJobParameters(topx=1001),
            lambda: RetroPathJobParameters(dmin=10, dmax=2),
            lambda: RetroPathJobParameters(msc_timeout=0),
        )
        for factory in invalid_parameter_factories:
            with self.subTest(factory=factory), self.assertRaises(ValueError):
                factory()
        with self.assertRaisesRegex(ValueError, "loopback"):
            RetroPathHttpClient("http://example.com:8765")


if __name__ == "__main__":
    unittest.main()
