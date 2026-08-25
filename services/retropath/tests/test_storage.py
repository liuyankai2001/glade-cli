from pathlib import Path

from app.storage import JobStorage


def test_storage_persists_and_recovers_running_jobs(tmp_path: Path):
    storage = JobStorage(tmp_path)
    storage.initialize()
    storage.create_job("rp2-test", {"max_steps": 1}, tmp_path / "rp2-test")
    claimed = storage.claim_job("rp2-test")
    assert claimed is not None
    assert claimed["status"] == "running"

    restarted = JobStorage(tmp_path)
    restarted.initialize()
    job = restarted.get_job("rp2-test")
    assert job is not None
    assert job["status"] == "failed"
    assert "restarted" in str(job["error"])

