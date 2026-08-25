from __future__ import annotations

import json
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Annotated

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile, status
from fastapi.responses import FileResponse
from pydantic import ValidationError

from .config import Settings, inspect_runtime
from .manager import JobManager, QueueFullError
from .models import JobParameters, TERMINAL_STATUSES, public_job
from .runner import RetroPathRunner
from .storage import JobStorage
from .validation import InputValidationError, validate_compound_csv


class ServiceState:
    def __init__(self) -> None:
        self.settings = Settings.from_env()
        self.runtime = inspect_runtime(self.settings)
        self.storage = JobStorage(self.settings.jobs_dir)
        self.manager = JobManager(
            self.storage,
            RetroPathRunner(self.settings, self.runtime),
            max_queue=self.settings.max_queue,
        )

    def start(self) -> None:
        self.manager.start()

    def stop(self) -> None:
        self.manager.stop()


@asynccontextmanager
async def lifespan(app: FastAPI):
    service = ServiceState()
    service.start()
    app.state.retropath = service
    try:
        yield
    finally:
        service.stop()


app = FastAPI(
    title="GLADE RetroPath local service",
    version="1.0.0",
    lifespan=lifespan,
)


def get_service(request: Request) -> ServiceState:
    return request.app.state.retropath


@app.get("/health")
def health(request: Request) -> dict[str, object]:
    service = get_service(request)
    result = service.runtime.public_dict(service.settings)
    result["service_version"] = app.version
    result["queue_active"] = service.storage.active_count()
    return result


@app.post("/v1/jobs", status_code=status.HTTP_202_ACCEPTED)
async def create_job(
    request: Request,
    source_file: Annotated[UploadFile, File(...)],
    sink_file: Annotated[UploadFile, File(...)],
    max_steps: Annotated[int, Form()] = 3,
    topx: Annotated[int, Form()] = 100,
    dmin: Annotated[int, Form()] = 2,
    dmax: Annotated[int, Form()] = 16,
    mwmax_source: Annotated[int, Form()] = 1000,
    msc_timeout: Annotated[int, Form()] = 10,
) -> dict[str, object]:
    service = get_service(request)
    if not service.runtime.ready:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail={"message": "RetroPath runtime is not ready", "errors": service.runtime.errors},
        )

    source_bytes = await source_file.read()
    sink_bytes = await sink_file.read()
    try:
        validate_compound_csv(source_bytes, kind="source")
        validate_compound_csv(sink_bytes, kind="sink")
        parameters = JobParameters(
            max_steps=max_steps,
            topx=topx,
            dmin=dmin,
            dmax=dmax,
            mwmax_source=mwmax_source,
            msc_timeout=msc_timeout,
        )
    except (InputValidationError, ValidationError) as exc:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(exc),
        ) from exc

    try:
        job = service.manager.submit(source_bytes, sink_bytes, parameters)
    except QueueFullError as exc:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail=str(exc),
        ) from exc
    return public_job(job)


@app.get("/v1/jobs/{job_id}")
def get_job(job_id: str, request: Request) -> dict[str, object]:
    service = get_service(request)
    job = service.storage.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    return public_job(job)


@app.get("/v1/jobs/{job_id}/results")
def get_results(job_id: str, request: Request) -> dict[str, object]:
    service = get_service(request)
    job = service.storage.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    if job["status"] not in TERMINAL_STATUSES:
        raise HTTPException(status_code=409, detail="job has not reached a terminal state")
    manifest_path = Path(str(job["job_dir"])) / "run_manifest.json"
    manifest = None
    artifacts: list[str] = []
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        artifacts = ["run_manifest.json", *manifest.get("artifacts", [])]
    return {
        **public_job(job),
        "manifest": manifest,
        "artifacts": artifacts,
    }


@app.get("/v1/jobs/{job_id}/artifacts/{artifact_path:path}")
def get_artifact(job_id: str, artifact_path: str, request: Request) -> FileResponse:
    service = get_service(request)
    job = service.storage.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    job_dir = Path(str(job["job_dir"])).resolve()
    candidate = (job_dir / artifact_path).resolve()
    if candidate != job_dir and job_dir not in candidate.parents:
        raise HTTPException(status_code=404, detail="artifact not found")

    results = get_results(job_id, request)
    allowed = set(results["artifacts"])
    normalized = str(candidate.relative_to(job_dir)).replace("\\", "/")
    if normalized not in allowed or not candidate.is_file():
        raise HTTPException(status_code=404, detail="artifact not found")
    return FileResponse(candidate, filename=candidate.name)

