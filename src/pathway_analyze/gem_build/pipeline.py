"""Rebuild through Bactabolize and publish only the original model JSON."""
from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path, PurePosixPath
from urllib.parse import quote

import httpx
from filelock import FileLock, Timeout as LockTimeout

from .client import BactabolizeClient, GemBuildError
from .config import validate_managed_paths
from .validation import validate_model

REQUIRED = {'raw/draft_model.json', 'service_validation.json'}
SHA = re.compile(r'^[0-9a-f]{64}$')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def atomic_bytes(path, data):
    """Replace the sole public output only after its complete bytes are durable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.' + path.name + '-', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'wb') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def safe_artifact(name):
    if (not isinstance(name, str) or not name or '\\' in name or ':' in name or '%' in name
            or name.startswith('/') or any(part in {'', '.', '..'} for part in name.split('/'))
            or PurePosixPath(name).is_absolute()):
        raise GemBuildError(f'unsafe_artifact: {name!r}')
    for part in name.split('/'):
        if (part.endswith(('.', ' ')) or any(ord(char) < 32 for char in part)
                or part.split('.')[0].upper() in {'CON', 'PRN', 'AUX', 'NUL',
                    *(f'COM{i}' for i in range(1, 10)), *(f'LPT{i}' for i in range(1, 10))}):
            raise GemBuildError(f'unsafe_artifact: {name!r}')
    return name


def _verified_manifest(results, job_id, config, inputs, health):
    manifest = results.get('manifest')
    if not isinstance(manifest, dict):
        raise GemBuildError('protocol_error: missing run manifest')
    for key, value in [('schema_version', 'bactabolize_run.v1'), ('job_id', job_id),
                       ('status', 'succeeded'), ('parameters', config.parameters),
                       ('input_sha256', inputs), ('versions', health['versions'])]:
        if manifest.get(key) != value:
            raise GemBuildError(f'provenance_mismatch: manifest {key}')
    resolved = manifest.get('resolved_parameters')
    if (not isinstance(resolved, dict) or any(resolved.get(key) != value
            for key, value in config.parameters.items() if key != 'biomass_reaction_id')
            or not isinstance(resolved.get('biomass_reaction_id'), str)
            or not resolved['biomass_reaction_id']):
        raise GemBuildError('protocol_error: invalid resolved_parameters')
    if (config.parameters['biomass_reaction_id'] is not None
            and resolved['biomass_reaction_id'] != config.parameters['biomass_reaction_id']):
        raise GemBuildError('provenance_mismatch: resolved biomass reaction')
    artifacts, hashes = manifest.get('artifacts'), manifest.get('artifact_sha256')
    if not isinstance(artifacts, list) or not isinstance(hashes, dict):
        raise GemBuildError('protocol_error: invalid artifact manifest')
    for name in artifacts:
        safe_artifact(name)
    if (len(artifacts) != len({name.casefold() for name in artifacts})
            or 'run_manifest.json' in artifacts or not REQUIRED.issubset(artifacts)
            or set(hashes) != set(artifacts)
            or any(not isinstance(value, str) or not SHA.fullmatch(value) for value in hashes.values())):
        raise GemBuildError('protocol_error: missing/duplicate/unhashed artifacts')
    offered = results.get('artifacts')
    if (not isinstance(offered, list) or any(not isinstance(name, str) for name in offered)
            or len(offered) != len(set(offered)) or set(offered) != {'run_manifest.json', *artifacts}):
        raise GemBuildError('protocol_error: results artifact whitelist differs from manifest')
    return manifest


def _lock_path(output):
    identity = os.path.normcase(str(output.resolve())).encode('utf-8')
    lock = Path(tempfile.gettempdir()) / 'glade-bactabolize-locks' / (sha(identity) + '.lock')
    if lock.parent.resolve().is_relative_to(output.parent.resolve()):
        raise GemBuildError('output_temp_conflict: use a dedicated model output directory '
                            'that does not contain the system temporary lock directory')
    return lock


def _validate_sources(config, lock):
    sources = [config.reference_model_file, config.reference_genome_file,
               config.target_genome_file, config.medium_path]
    validate_managed_paths(config.model_file, sources)
    for source in sources:
        if source.resolve() == lock.resolve() or (source.is_file() and lock.is_file() and source.samefile(lock)):
            raise ValueError(f'lock must not alias an input or medium file: {source}')


def _download(transport, job_id, name, manifest):
    data = transport.request('GET', f'/v1/jobs/{job_id}/artifacts/' + quote(name, safe='/')).content
    if sha(data) != manifest['artifact_sha256'][name]:
        raise GemBuildError(f'artifact_hash_mismatch: {name}')
    return data


def _build(config, client):
    transport = BactabolizeClient(config, client)
    health = transport.health()
    inputs = {name: sha(path.read_bytes()) for name, path in [
        ('reference_model.json', config.reference_model_file),
        ('reference.gbff', config.reference_genome_file), ('target.gbff', config.target_genome_file)]}
    job_id, terminal = None, False
    try:
        job = transport.submit()
        job_id = job['job_id']
        print(f'[INFO] Bactabolize job submitted: {job_id}')
        job = transport.wait(job_id, lambda _: None)
        terminal = True
        if job['status'] != 'succeeded':
            raise GemBuildError(f'build_{job["status"]}: {job.get("failure_code")}: {job.get("error")} (job_id={job_id})')
        results = transport.job(transport.json('GET', f'/v1/jobs/{job_id}/results'), job_id)
        if results.get('status') != 'succeeded':
            raise GemBuildError('protocol_error: results do not match successful terminal state')
        manifest = _verified_manifest(results, job_id, config, inputs, health)
        # Temporary validation data never live beside the public JSON. Cleanup
        # finishes before publication, so a cleanup error cannot replace an old model.
        with tempfile.TemporaryDirectory(prefix='glade-bactabolize-run-') as temporary:
            staging = Path(temporary)
            json_data = _download(transport, job_id, 'raw/draft_model.json', manifest)
            service_data = _download(transport, job_id, 'service_validation.json', manifest)
            service_validation = json.loads(service_data)
            if not isinstance(service_validation, dict) or service_validation.get('model_readable') is not True:
                raise GemBuildError('invalid_model: service did not validate a readable model')
            if (not isinstance(service_validation.get('growth'), dict)
                    or not isinstance(service_validation['growth'].get('passed'), bool)):
                raise GemBuildError('protocol_error: invalid service growth validation report')
            model = staging / 'model.json'
            model.write_bytes(json_data)
            validation = validate_model(model, config.medium_path)
            if manifest['resolved_parameters']['biomass_reaction_id'] not in validation['objective_reaction_ids']:
                raise GemBuildError('invalid_model: resolved biomass is not a model objective')
            warnings = []
            if not validation['growth']['passed']:
                warnings.append('Model is a readable draft; growth on the project default medium was not validated.')
            if service_validation['growth']['passed'] is not True:
                warnings.append('Bactabolize media/atmosphere growth validation did not pass; retain as an unvalidated draft.')
        atomic_bytes(config.model_file, json_data)
        return {'status': 'succeeded', 'job_id': job_id, 'model_file': str(config.model_file),
                'validation': validation, 'warnings': warnings}
    except BaseException:
        if job_id is not None and not terminal:
            try:
                transport.request('POST', f'/v1/jobs/{job_id}/cancel')
            except (GemBuildError, OSError):
                pass
        raise


def build_model(config, *, client=None):
    """Always submit a new job. ``force`` is accepted only for CLI compatibility."""
    try:
        lock = _lock_path(config.model_file)
        _validate_sources(config, lock)
        lock.parent.mkdir(parents=True, exist_ok=True)
        with FileLock(str(lock), timeout=config.wait_timeout_seconds):
            _validate_sources(config, lock)
            if client is not None:
                return _build(config, client)
            with httpx.Client(follow_redirects=False) as owned:
                return _build(config, owned)
    except LockTimeout as exc:
        raise GemBuildError(f'output_busy: another gem-build owns {config.model_file}') from exc
    except (OSError, ValueError, TypeError, KeyError) as exc:
        raise GemBuildError(f'build_io_or_protocol_error: {exc}') from exc


def run_gem_build(config):
    try:
        result = build_model(config)
    except GemBuildError as exc:
        print(f'[ERROR] gem-build: {exc}')
        raise SystemExit(1) from exc
    except KeyboardInterrupt:
        print('[ERROR] gem-build interrupted; cancellation requested for any running accepted job.')
        raise SystemExit(130)
    print(f'[INFO] model: {result["model_file"]} (job {result["job_id"]})')
    for warning in result['warnings']:
        print(f'[WARNING] {warning}')
    return result


__all__ = ['build_model', 'run_gem_build', 'GemBuildError']
