from __future__ import annotations

import json
import sqlite3
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


class JobStorage:
    def __init__(self, jobs_dir: Path):
        self.jobs_dir = jobs_dir
        self.database_path = jobs_dir / "jobs.sqlite3"
        self._write_lock = threading.Lock()

    def initialize(self) -> None:
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        with self._connect() as connection:
            connection.execute("PRAGMA journal_mode=WAL")
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    job_id TEXT PRIMARY KEY,
                    status TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    started_at TEXT,
                    finished_at TEXT,
                    return_code INTEGER,
                    failure_code TEXT,
                    error TEXT,
                    parameters_json TEXT NOT NULL,
                    job_dir TEXT NOT NULL
                )
                """
            )
            columns = {
                str(row[1])
                for row in connection.execute("PRAGMA table_info(jobs)").fetchall()
            }
            if "failure_code" not in columns:
                connection.execute(
                    "ALTER TABLE jobs ADD COLUMN failure_code TEXT"
                )
            connection.execute(
                """
                UPDATE jobs
                   SET status = 'failed',
                       finished_at = ?,
                       failure_code = 'service_restarted',
                       error = 'service restarted while the job was running'
                 WHERE status = 'running'
                """,
                (utcnow(),),
            )

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.database_path, timeout=30)
        connection.row_factory = sqlite3.Row
        return connection

    def active_count(self) -> int:
        with self._connect() as connection:
            row = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs WHERE status IN ('queued', 'running')"
            ).fetchone()
        return int(row["count"])

    def create_job(
        self, job_id: str, parameters: dict[str, int], job_dir: Path
    ) -> dict[str, object]:
        created_at = utcnow()
        with self._write_lock, self._connect() as connection:
            connection.execute(
                """
                INSERT INTO jobs (
                    job_id, status, created_at, parameters_json, job_dir
                ) VALUES (?, 'queued', ?, ?, ?)
                """,
                (job_id, created_at, json.dumps(parameters), str(job_dir)),
            )
        job = self.get_job(job_id)
        assert job is not None
        return job

    def queued_job_ids(self) -> Iterable[str]:
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT job_id FROM jobs WHERE status = 'queued' ORDER BY created_at"
            ).fetchall()
        return [str(row["job_id"]) for row in rows]

    def claim_job(self, job_id: str) -> dict[str, object] | None:
        started_at = utcnow()
        with self._write_lock, self._connect() as connection:
            cursor = connection.execute(
                """
                UPDATE jobs
                   SET status = 'running', started_at = ?,
                       failure_code = NULL, error = NULL
                 WHERE job_id = ? AND status = 'queued'
                """,
                (started_at, job_id),
            )
            if cursor.rowcount != 1:
                return None
        return self.get_job(job_id)

    def finish_job(
        self,
        job_id: str,
        *,
        status: str,
        return_code: int | None,
        error: str | None,
        failure_code: str | None = None,
    ) -> None:
        with self._write_lock, self._connect() as connection:
            connection.execute(
                """
                UPDATE jobs
                   SET status = ?, finished_at = ?, return_code = ?,
                       failure_code = ?, error = ?
                 WHERE job_id = ?
                """,
                (
                    status,
                    utcnow(),
                    return_code,
                    failure_code,
                    error,
                    job_id,
                ),
            )

    def get_job(self, job_id: str) -> dict[str, object] | None:
        with self._connect() as connection:
            row = connection.execute(
                "SELECT * FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
        if row is None:
            return None
        job = dict(row)
        job["parameters"] = json.loads(str(job.pop("parameters_json")))
        return job
