import hashlib
import os
import time
from pathlib import Path

from fastapi.testclient import TestClient

from app.main import app


SOURCE = b'Name,InChI\n"target","InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"\n'
RULES = b"Rule ID,Rule,Score\nrule-1,[#6:1]>>[#6:1],1\n"


def test_health_and_source_in_sink_job(tmp_path: Path, monkeypatch):
    rules_path = tmp_path / "rules.csv"
    rules_path.write_bytes(RULES)
    knime_dir = tmp_path / "knime"
    knime_dir.mkdir()
    knime_exec = knime_dir / "knime"
    knime_exec.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    os.chmod(knime_exec, 0o755)
    (knime_dir / "p2" / "plugins").mkdir(parents=True)
    (
        knime_dir / "p2" / "plugins" / "org.rdkit.knime.nodes_4.9.1.jar"
    ).touch()
    (knime_dir / "p2" / "features").mkdir()
    (
        knime_dir
        / "p2"
        / "features"
        / "org.knime.features.chem.types_4.7.1"
    ).mkdir()
    openssl10_dir = tmp_path / "openssl10"
    openssl10_dir.mkdir()
    (openssl10_dir / "libssl.so.10").touch()
    (openssl10_dir / "libcrypto.so.10").touch()

    monkeypatch.setenv("RETROPATH_RULES_PATH", str(rules_path))
    monkeypatch.setenv(
        "RETROPATH_RULES_SHA256", hashlib.sha256(RULES).hexdigest()
    )
    monkeypatch.setenv("RETROPATH_KNIME_DIR", str(knime_dir))
    monkeypatch.setenv("RETROPATH_OPENSSL10_DIR", str(openssl10_dir))
    monkeypatch.setenv("RETROPATH_JOBS_DIR", str(tmp_path / "jobs"))

    with TestClient(app) as client:
        response = client.get("/health")
        assert response.status_code == 200
        assert response.json()["ready"] is True, response.json()

        response = client.post(
            "/v1/jobs",
            files={
                "source_file": ("source.csv", SOURCE, "text/csv"),
                "sink_file": ("sink.csv", SOURCE, "text/csv"),
            },
            data={"max_steps": "1", "topx": "1", "dmin": "2", "dmax": "2"},
        )
        assert response.status_code == 202
        job_id = response.json()["job_id"]

        payload = None
        for _ in range(100):
            payload = client.get(f"/v1/jobs/{job_id}").json()
            if payload["status"] not in {"queued", "running"}:
                break
            time.sleep(0.05)
        assert payload is not None
        assert payload["status"] == "source_in_sink"

        results = client.get(f"/v1/jobs/{job_id}/results")
        assert results.status_code == 200
        assert "run_manifest.json" in results.json()["artifacts"]

        blocked = client.get(
            f"/v1/jobs/{job_id}/artifacts/../../jobs.sqlite3"
        )
        assert blocked.status_code == 404
