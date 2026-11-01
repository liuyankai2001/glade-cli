"""Bounded synchronous HTTP transport; never launches a local container."""
from __future__ import annotations

import json
import re
import time
from contextlib import ExitStack

import httpx

JOB_ID = re.compile(r'^bact-[0-9a-f]{32}$')
TERMINAL = {'succeeded', 'failed', 'timed_out', 'cancelled'}
STATUSES = TERMINAL | {'queued', 'running'}
VERSION_KEYS = {'python', 'bactabolize', 'cobra', 'biopython', 'blast', 'prodigal'}


class GemBuildError(RuntimeError):
    """Operational or integrity failure, leaving the formal model unchanged."""


class BactabolizeClient:
    def __init__(self, config, client):
        self.config = config
        self.client = client

    def request(self, method, path, **kwargs):
        attempts = self.config.get_attempts if method == 'GET' else 1
        for attempt in range(attempts):
            try:
                response = self.client.request(method, self.config.service_url + path,
                    timeout=self.config.request_timeout_seconds, **kwargs)
            except httpx.TransportError as exc:
                if method == 'POST' and path == '/v1/jobs':
                    raise GemBuildError('submission_uncertain: request was not retried; check service jobs before resubmitting') from exc
                if attempt + 1 == attempts:
                    raise GemBuildError(f'transport_error: {method} {path}: {exc}') from exc
            else:
                if response.status_code < 400:
                    return response
                if response.status_code not in {408, 429, 500, 502, 503, 504} or attempt + 1 == attempts:
                    raise GemBuildError(f'http_error: {method} {path}: HTTP {response.status_code}: {response.text[:500]}')
            time.sleep(min(0.25 * 2 ** attempt, 2))
        raise GemBuildError('transport_error: no GET attempts configured')

    def json(self, method, path, **kwargs):
        response = self.request(method, path, **kwargs)
        try:
            value = response.json()
            if not isinstance(value, dict):
                raise ValueError('expected object')
            return value
        except (ValueError, TypeError) as exc:
            raise GemBuildError(f'protocol_error: invalid JSON at {path}') from exc

    def health(self):
        result = self.json('GET', '/health')
        if result.get('ready') is not True:
            raise GemBuildError(f'service_not_ready: {result.get("errors", [])}')
        versions = result.get('versions')
        if (result.get('service_version') != '1.0.0' or not isinstance(versions, dict)
                or not VERSION_KEYS.issubset(versions)
                or any(not isinstance(versions[k], str) or not versions[k] for k in VERSION_KEYS)):
            raise GemBuildError('protocol_error: missing or unsupported service/runtime versions')
        return result

    def job(self, value, expected=None):
        job_id = value.get('job_id')
        if not isinstance(job_id, str) or not JOB_ID.fullmatch(job_id) or (expected and job_id != expected):
            raise GemBuildError('protocol_error: invalid job_id')
        if value.get('status') not in STATUSES or value.get('parameters') != self.config.parameters:
            raise GemBuildError('protocol_error: unexpected job status/parameters')
        return value

    def submit(self):
        with ExitStack() as stack:
            files = {}
            for field, filename in [('reference_model_file', 'reference_model.json'),
                                    ('reference_genome_file', 'reference.gbff'), ('target_genome_file', 'target.gbff')]:
                files[field] = (filename, stack.enter_context(getattr(self.config, field).open('rb')), 'application/octet-stream')
            return self.job(self.json('POST', '/v1/jobs', files=files,
                                     data={'parameters': json.dumps(self.config.parameters, allow_nan=False)}))

    def wait(self, job_id, record):
        deadline = time.monotonic() + self.config.wait_timeout_seconds
        while True:
            job = self.job(self.json('GET', f'/v1/jobs/{job_id}'), job_id)
            record(job)
            if job['status'] in TERMINAL:
                return job
            if time.monotonic() >= deadline:
                raise GemBuildError(f'client_poll_timeout: job {job_id} exceeded client wait limit')
            time.sleep(min(self.config.poll_interval_seconds, max(0, deadline - time.monotonic())))
