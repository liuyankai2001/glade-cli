from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Optional
import hashlib
import os
import platform
import subprocess


@dataclass(frozen=True)
class Settings:
    jobs_dir: Path = Path('/var/lib/bactabolize/jobs')
    max_queue: int = 8
    job_timeout_seconds: float = 1800
    upload_limit_bytes: int = 67108864
    shutdown_grace_seconds: float = 30

    @classmethod
    def from_env(cls):
        return cls(Path(os.getenv('BACTABOLIZE_JOBS_DIR', '/var/lib/bactabolize/jobs')),
                   int(os.getenv('BACTABOLIZE_MAX_QUEUE', '8')),
                   float(os.getenv('BACTABOLIZE_JOB_TIMEOUT_SECONDS', '1800')),
                   int(os.getenv('BACTABOLIZE_UPLOAD_LIMIT_BYTES', '67108864')),
                   float(os.getenv('BACTABOLIZE_SHUTDOWN_GRACE_SECONDS', '30')))


@dataclass(frozen=True)
class RuntimeInfo:
    versions: dict
    errors: tuple

    @property
    def ready(self):
        return not self.errors


def inspect_runtime(settings):
    versions = {'python': platform.python_version()}
    errors = []
    for key, package in [('bactabolize','bactabolize'), ('cobra','cobra'), ('biopython','biopython')]:
        try:
            versions[key] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[key] = 'missing'
            errors.append(package + ' is not installed')
    for key, command in [('blast',['blastp','-version']), ('prodigal',['prodigal','-v'])]:
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=10)
            versions[key] = (result.stdout + result.stderr).strip().splitlines()[0]
            if result.returncode:
                errors.append(key + ' version probe failed')
        except (OSError, subprocess.SubprocessError, IndexError) as exc:
            versions[key] = 'missing'
            errors.append(key + ': ' + str(exc))
    return RuntimeInfo(versions, tuple(errors))


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()
