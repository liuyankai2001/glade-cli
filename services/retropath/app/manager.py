from __future__ import annotations

import queue
import threading
import uuid
from pathlib import Path

from .models import JobParameters
from .runner import RetroPathRunner
from .storage import JobStorage


class QueueFullError(RuntimeError):
    pass


class JobManager:
    def __init__(
        self,
        storage: JobStorage,
        runner: RetroPathRunner,
        *,
        max_queue: int,
    ):
        self.storage = storage
        self.runner = runner
        self.max_queue = max_queue
        self._queue: queue.Queue[str | None] = queue.Queue()
        self._stop = threading.Event()
        self._thread = threading.Thread(
            target=self._worker_loop,
            name="retropath-worker",
            daemon=True,
        )

    def start(self) -> None:
        self.storage.initialize()
        for job_id in self.storage.queued_job_ids():
            self._queue.put(job_id)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        self._queue.put(None)
        self._thread.join(timeout=5)

    def submit(
        self,
        source_bytes: bytes,
        sink_bytes: bytes,
        parameters: JobParameters,
    ) -> dict[str, object]:
        if self.storage.active_count() >= self.max_queue:
            raise QueueFullError(
                f"RetroPath queue already contains {self.max_queue} active jobs"
            )
        job_id = f"rp2-{uuid.uuid4().hex}"
        job_dir = self.storage.jobs_dir / job_id
        input_dir = job_dir / "input"
        input_dir.mkdir(parents=True, exist_ok=False)
        (input_dir / "source.csv").write_bytes(source_bytes)
        (input_dir / "sink.csv").write_bytes(sink_bytes)
        job = self.storage.create_job(
            job_id,
            parameters.model_dump(),
            job_dir,
        )
        self._queue.put(job_id)
        return job

    def _worker_loop(self) -> None:
        while not self._stop.is_set():
            try:
                job_id = self._queue.get(timeout=1)
            except queue.Empty:
                continue
            if job_id is None:
                self._queue.task_done()
                break
            job = self.storage.claim_job(job_id)
            if job is None:
                self._queue.task_done()
                continue
            try:
                result = self.runner.run(job)
                self.storage.finish_job(
                    job_id,
                    status=result.status,
                    return_code=result.return_code,
                    error=result.error,
                )
            except Exception as exc:  # pragma: no cover - defensive boundary
                self.storage.finish_job(
                    job_id,
                    status="failed",
                    return_code=None,
                    error=f"internal worker error: {type(exc).__name__}: {exc}",
                )
            finally:
                self._queue.task_done()

