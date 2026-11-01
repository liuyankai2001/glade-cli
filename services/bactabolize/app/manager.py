import json
import os
import queue
import subprocess
import threading
import time
from pathlib import Path

from .config import sha256_file
from .models import INPUTS, TERMINAL, validate_parameters
from .runner import worker_command, worker_environment, terminate_group
from .storage import JobStorage


class JobManager:
    def __init__(self, settings, runtime, runner=None):
        self.settings, self.runtime = settings, runtime
        self.storage = JobStorage(settings.jobs_dir)
        self.command = runner or worker_command
        self.queue = queue.Queue()
        self.lock = threading.RLock()
        self.reservations = 0
        self.cancelled = set()
        self.stop_event = threading.Event()
        self.lifecycle_error = None
        self.thread = threading.Thread(target=self.loop,daemon=True,name='bactabolize-worker')

    def start(self):
        self.storage.initialize()
        # Interrupted runs acquire a fresh terminal manifest after restart.
        with self.storage.connect() as db:
            ids = [r[0] for r in db.execute("SELECT job_id FROM jobs WHERE failure_code='service_restarted'")]
        for job_id in ids:
            job = self.storage.get_job(job_id)
            self.finish(job,job['status'],job['return_code'],job['failure_code'],job['error'])
        for job_id in self.storage.queued_job_ids():
            self.queue.put(job_id)
        self.thread.start()

    def stop(self):
        self.stop_event.set()
        self.queue.put(None)
        self.thread.join(timeout=self.settings.shutdown_grace_seconds + 12)

    def reserve(self):
        with self.lock:
            if self.health_errors():
                return False
            if self.storage.active_count() + self.reservations >= self.settings.max_queue:
                return False
            self.reservations += 1
            return True

    def health_errors(self):
        errors = list(self.runtime.errors)
        if not self.thread.is_alive():
            errors.append('Bactabolize worker thread is not alive')
        if self.stop_event.is_set():
            errors.append('Bactabolize worker is shutting down')
        if self.lifecycle_error:
            errors.append(self.lifecycle_error)
        return errors

    def cancel(self, job_id):
        with self.lock:
            job = self.storage.get_job(job_id)
            if job and job['status'] not in TERMINAL:
                self.cancelled.add(job_id)
                if job['status'] == 'queued':
                    self.finish(job,'cancelled',failure_code='cancelled',error='cancelled by request')
            return self.storage.get_job(job_id)

    def finish(self, job, status, return_code=None, failure_code=None, error=None):
        # Write terminal manifest first; public status becomes terminal only when
        # every downloadable artifact and its hash are finalized.
        final = dict(job,status=status,return_code=return_code,failure_code=failure_code,error=error)
        try:
            self.manifest(final)
        except Exception as exc:
            detail = 'manifest finalization failed: ' + str(exc)
            error = (error + '; ' if error else '') + detail
            if status == 'succeeded':
                status, failure_code = 'failed','manifest_finalization_failed'
            final.update(status=status,failure_code=failure_code,error=error)
            # A failed job may have partial scientific metadata. Try a tolerant
            # terminal manifest; persistence must still finish if disk writing
            # fails too. Never recurse through the same failing error handler.
            try:
                self.manifest(final)
            except Exception as fallback_exc:
                error += '; fallback manifest failed: ' + str(fallback_exc)
        for attempt in range(2):
            try:
                self.storage.finish_job(job['job_id'],status,return_code,failure_code,error)
                break
            except Exception as exc:
                if attempt == 1:
                    self.lifecycle_error = 'terminal job persistence failed: ' + str(exc)
        self.cancelled.discard(job['job_id'])

    def manifest(self, job):
        path = Path(job['job_dir'])
        strict = job['status'] == 'succeeded'
        warnings = []
        allowed = ['bactabolize.log','input_manifest.json','service_validation.json']
        if (path / 'raw').exists():
            allowed.extend(p.relative_to(path).as_posix() for p in sorted((path/'raw').rglob('*')) if p.is_file() and not p.is_symlink())
        artifacts, artifact_hashes = [], {}
        for artifact in allowed:
            try:
                if (path/artifact).is_file():
                    artifact_hashes[artifact] = sha256_file(path/artifact)
                    artifacts.append(artifact)
            except OSError as exc:
                if strict:
                    raise
                warnings.append('artifact unavailable: ' + artifact + ': ' + str(exc))
        resolved = path / 'resolved_parameters.json'
        def reject_constant(value):
            raise ValueError('nonfinite JSON constant: ' + value)
        try:
            resolved_parameters = json.loads(resolved.read_text(encoding='utf-8'),parse_constant=reject_constant)
            validated = validate_parameters(resolved_parameters)
            if not validated['biomass_reaction_id'] or resolved_parameters.get('no_reannotation') is not True:
                raise ValueError('resolved parameters require biomass_reaction_id and no_reannotation=true')
        except (OSError,ValueError,TypeError) as exc:
            if strict:
                raise ValueError('invalid resolved_parameters.json: ' + str(exc)) from exc
            resolved_parameters = None
            warnings.append('resolved parameters unavailable: ' + str(exc))
        input_hashes = {}
        for name in INPUTS.values():
            try:
                input_hashes[name] = sha256_file(path/'input'/name)
            except OSError as exc:
                if strict:
                    raise
                warnings.append('input unavailable: ' + name + ': ' + str(exc))
        manifest = dict(schema_version='bactabolize_run.v1',job_id=job['job_id'],status=job['status'],parameters=job['parameters'],
                        resolved_parameters=resolved_parameters,
                        versions=self.runtime.versions,input_sha256=input_hashes,
                        artifacts=artifacts,artifact_sha256=artifact_hashes,warnings=warnings)
        temporary = path / 'run_manifest.tmp'
        temporary.write_text(json.dumps(manifest,indent=2,allow_nan=False),encoding='utf-8')
        temporary.replace(path/'run_manifest.json')

    def loop(self):
        while not self.stop_event.is_set():
            job_id = self.queue.get()
            if job_id is None:
                break
            with self.lock:
                try:
                    job = self.storage.claim_job(job_id)
                except Exception as exc:
                    self.lifecycle_error = 'worker could not claim queued job: ' + str(exc)
                    continue
            if job is None:
                continue
            try:
                self.run(job)
            except Exception as exc:
                self.finish(job,'failed',failure_code='worker_error',error=str(exc))

    def run(self, job):
        path = Path(job['job_dir'])
        (path/'parameters.json').write_text(json.dumps(job['parameters']),encoding='utf-8')
        start = time.monotonic()
        status, code, error, failure = 'failed', None, None, None
        with (path/'bactabolize.log').open('wb') as log:
            kwargs = dict(cwd=str(path),env=worker_environment(),stdout=log,stderr=subprocess.STDOUT)
            if os.name == 'nt':
                kwargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP
            else:
                kwargs['start_new_session'] = True
            process = subprocess.Popen(self.command(path),**kwargs)
            try:
                while process.poll() is None:
                    if job['job_id'] in self.cancelled or self.stop_event.is_set():
                        status, failure, error = 'cancelled','cancelled','cancelled by shutdown' if self.stop_event.is_set() else 'cancelled by request'
                        terminate_group(process,min(1,self.settings.shutdown_grace_seconds))
                        break
                    if time.monotonic() - start >= self.settings.job_timeout_seconds:
                        status, failure, error = 'timed_out','wall_timeout','job exceeded wall timeout'
                        terminate_group(process)
                        break
                    self.stop_event.wait(.05)
                code = process.wait()
                if failure is None:
                    if code:
                        failure, error = 'bactabolize_failed','worker exited with code ' + str(code)
                    else:
                        validation = json.loads((path/'service_validation.json').read_text(encoding='utf-8'))
                        if validation.get('model_readable') is not True:
                            raise ValueError('worker model validation failed')
                        status = 'succeeded'
            finally:
                if process.poll() is None:
                    terminate_group(process)
        with self.lock:
            if job['job_id'] in self.cancelled and status == 'succeeded':
                status, failure, error = 'cancelled','cancelled','cancelled by request'
            self.finish(job,status,code,failure,error)
