import json
import shutil
import uuid
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse

from .config import Settings, inspect_runtime, sha256_file
from .manager import JobManager
from .models import INPUTS, TERMINAL, public_job, validate_parameters


def create_app(settings=None, runner=None, runtime=None):
    @asynccontextmanager
    async def lifespan(app):
        config = settings or Settings.from_env()
        manager = JobManager(config,runtime or inspect_runtime(config),runner)
        app.state.bactabolize = manager
        manager.start()
        try:
            yield
        finally:
            manager.stop()

    app = FastAPI(title='GLADE Bactabolize service',version='1.0.0',lifespan=lifespan)

    def service(request):
        return request.app.state.bactabolize

    def lookup(request, job_id):
        job = service(request).storage.get_job(job_id)
        if job is None:
            raise HTTPException(404,'job not found')
        return job

    @app.get('/health')
    def health(request: Request):
        manager = service(request)
        errors = manager.health_errors()
        return dict(ready=not errors,service_version='1.0.0',versions=manager.runtime.versions,worker_concurrency=1,
                    job_timeout_seconds=manager.settings.job_timeout_seconds,queue_active=manager.storage.active_count(),errors=errors)

    @app.post('/v1/jobs',status_code=202)
    async def submit(request: Request, reference_model_file: UploadFile = File(...), reference_genome_file: UploadFile = File(...), target_genome_file: UploadFile = File(...), parameters: str = Form('{}')):
        manager = service(request)
        if manager.health_errors():
            raise HTTPException(503,{'message':'Bactabolize service is not ready','errors':manager.health_errors()})
        try:
            params = validate_parameters(json.loads(parameters))
        except (ValueError,TypeError) as exc:
            raise HTTPException(400,str(exc))
        uploads = dict(reference_model_file=reference_model_file,reference_genome_file=reference_genome_file,target_genome_file=target_genome_file)
        for field in ('reference_genome_file','target_genome_file'):
            if not (uploads[field].filename or '').lower().endswith('.gbff'):
                raise HTTPException(400,'genomes must be annotated GBFF files')
        if not manager.reserve():
            raise HTTPException(503,'Bactabolize queue is full')
        job_id = 'bact-' + uuid.uuid4().hex
        job_dir = manager.storage.jobs_dir / job_id
        created = False
        try:
            (job_dir/'input').mkdir(parents=True)
            for field, name in INPUTS.items():
                size = 0
                with (job_dir/'input'/name).open('wb') as stream:
                    while True:
                        chunk = await uploads[field].read(1024 * 1024)
                        if not chunk:
                            break
                        size += len(chunk)
                        if size > manager.settings.upload_limit_bytes:
                            raise HTTPException(413,'file exceeds upload limit')
                        stream.write(chunk)
                if size == 0:
                    raise HTTPException(400,'empty input file')
            with manager.lock:
                job = manager.storage.create_job(job_id,params,job_dir)
                created = True
                manager.queue.put(job_id)
            return public_job(job)
        finally:
            with manager.lock:
                manager.reservations -= 1
            if not created:
                shutil.rmtree(job_dir,ignore_errors=True)
            for upload in uploads.values():
                await upload.close()

    @app.get('/v1/jobs/{job_id}')
    def get_job(job_id: str, request: Request):
        return public_job(lookup(request,job_id))

    @app.post('/v1/jobs/{job_id}/cancel')
    def cancel(job_id: str, request: Request):
        lookup(request,job_id)
        return public_job(service(request).cancel(job_id))

    @app.get('/v1/jobs/{job_id}/results')
    def results(job_id: str, request: Request):
        job = lookup(request,job_id)
        if job['status'] not in TERMINAL:
            raise HTTPException(409,'job is not terminal')
        path = Path(job['job_dir'])/'run_manifest.json'
        manifest = json.loads(path.read_text(encoding='utf-8')) if path.is_file() else None
        return dict(**public_job(job),manifest=manifest,artifacts=['run_manifest.json'] + manifest['artifacts'] if manifest else [])

    @app.get('/v1/jobs/{job_id}/artifacts/{artifact_path:path}')
    def artifact(job_id: str, artifact_path: str, request: Request):
        result = results(job_id,request)
        if artifact_path not in result['artifacts'] or '\\' in artifact_path:
            raise HTTPException(404,'artifact not found')
        root = Path(lookup(request,job_id)['job_dir']).resolve()
        path = (root/artifact_path).resolve()
        if root not in path.parents or not path.is_file():
            raise HTTPException(404,'artifact not found')
        if artifact_path != 'run_manifest.json' and sha256_file(path) != result['manifest']['artifact_sha256'].get(artifact_path):
            raise HTTPException(409,'artifact checksum mismatch')
        return FileResponse(path,filename=path.name)

    return app


app = create_app()
