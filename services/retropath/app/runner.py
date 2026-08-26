from __future__ import annotations

import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from .config import RuntimeInfo, Settings
from .storage import utcnow


@dataclass(frozen=True)
class RunResult:
    status: str
    return_code: int | None
    error: str | None
    failure_code: str | None = None


def status_for_return_code(return_code: int) -> str:
    if return_code == 0:
        return "succeeded"
    if return_code == 10:
        return "source_in_sink"
    if return_code == 11:
        return "no_solution"
    return "failed"


def hash_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def wrapper_environment() -> dict[str, str]:
    environment = os.environ.copy()
    # retropath2_wrapper appends an old CentOS 7 sysroot whenever
    # CONDA_PREFIX is present. That makes KNIME load an ATK 2.28 library ahead
    # of Ubuntu's matching ATK/GTK stack. The Python executable and libraries
    # are already absolute/in LD_LIBRARY_PATH, so the wrapper does not need
    # CONDA_PREFIX at runtime.
    environment.pop("CONDA_PREFIX", None)
    return environment


class RetroPathRunner:
    def __init__(
        self,
        settings: Settings,
        runtime: RuntimeInfo,
    ):
        self.settings = settings
        self.runtime = runtime

    def run(self, job: dict[str, object]) -> RunResult:
        job_dir = Path(str(job["job_dir"]))
        input_dir = job_dir / "input"
        source_path = input_dir / "source.csv"
        sink_path = input_dir / "sink.csv"
        raw_dir = job_dir / "raw"
        stdout_path = job_dir / "stdout.log"
        stderr_path = job_dir / "stderr.log"
        parameters = dict(job["parameters"])
        started_at = utcnow()

        command = [
            sys.executable,
            "-m",
            "retropath2_wrapper",
            str(sink_path),
            str(self.settings.rules_path),
            str(raw_dir),
            "--source_file",
            str(source_path),
            "--kinstall",
            str(self.settings.knime_dir),
            "--rp2_version",
            self.settings.workflow_version,
            "--max_steps",
            str(parameters["max_steps"]),
            "--topx",
            str(parameters["topx"]),
            "--dmin",
            str(parameters["dmin"]),
            "--dmax",
            str(parameters["dmax"]),
            "--mwmax_source",
            str(parameters["mwmax_source"]),
            "--msc_timeout",
            str(parameters["msc_timeout"]),
            "--std_hydrogen",
            "auto",
            "--score_mode",
            "auto",
            "--log",
            "info",
        ]

        return_code: int | None = None
        status = "failed"
        error: str | None = None
        failure_code: str | None = None
        timed_out = False
        resource_exhausted = False
        peak_memory_bytes: int | None = None
        peak_working_set_bytes: int | None = None
        memory_samples = 0
        consecutive_high_memory_samples = 0
        maximum_consecutive_high_memory_samples = 0
        started_monotonic = time.monotonic()
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            process = subprocess.Popen(
                command,
                cwd=job_dir,
                env=wrapper_environment(),
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
            while process.poll() is None:
                elapsed = time.monotonic() - started_monotonic
                if elapsed >= self.settings.job_timeout_seconds:
                    timed_out = True
                    status = "timed_out"
                    failure_code = "wall_timeout"
                    error = (
                        "RetroPath execution exceeded "
                        f"{self.settings.job_timeout_seconds} seconds"
                    )
                    self._terminate_process_group(process)
                    break

                memory_bytes, working_set_bytes = self._read_memory_sample(
                    self.settings.cgroup_memory_current_path,
                    self.settings.cgroup_memory_stat_path,
                )
                if memory_bytes is not None:
                    memory_samples += 1
                    peak_memory_bytes = max(peak_memory_bytes or 0, memory_bytes)
                    monitored_bytes = (
                        working_set_bytes
                        if working_set_bytes is not None
                        else memory_bytes
                    )
                    peak_working_set_bytes = max(
                        peak_working_set_bytes or 0,
                        monitored_bytes,
                    )
                    if monitored_bytes > self.settings.memory_limit_bytes:
                        consecutive_high_memory_samples += 1
                        maximum_consecutive_high_memory_samples = max(
                            maximum_consecutive_high_memory_samples,
                            consecutive_high_memory_samples,
                        )
                    else:
                        consecutive_high_memory_samples = 0
                    if (
                        consecutive_high_memory_samples
                        >= self.settings.memory_limit_consecutive_samples
                    ):
                        resource_exhausted = True
                        status = "failed"
                        failure_code = "resource_exhausted"
                        error = (
                            "RetroPath cgroup working set exceeded "
                            f"{self.settings.memory_limit_bytes} bytes for "
                            f"{consecutive_high_memory_samples} consecutive samples"
                        )
                        self._terminate_process_group(process)
                        break
                time.sleep(self.settings.resource_poll_seconds)

            return_code = process.returncode
            if not timed_out and not resource_exhausted:
                if return_code is None:
                    return_code = process.wait(timeout=10)
                status = status_for_return_code(return_code)
                if status == "failed":
                    failure_code = "knime_execution_failed"
                    error = f"retropath2_wrapper exited with code {return_code}"

        wall_seconds = round(time.monotonic() - started_monotonic, 6)

        artifacts = self._discover_artifacts(job_dir)
        manifest = {
            "schema_version": 2,
            "job_id": job["job_id"],
            "status": status,
            "created_at": job["created_at"],
            "started_at": started_at,
            "finished_at": utcnow(),
            "return_code": return_code,
            "timed_out": timed_out,
            "failure_code": failure_code,
            "error": error,
            "parameters": parameters,
            "versions": {
                "retropath2_wrapper": self.runtime.wrapper_version,
                "retropath2_wrapper_reported": self.runtime.wrapper_reported_version,
                "workflow": self.settings.workflow_version,
                "knime": self.settings.knime_version,
                "knime_rdkit_nodes": self.settings.rdkit_plugin_version,
                "rules": self.settings.rules_version,
            },
            "rules_sha256": self.runtime.rules_sha256,
            "input_sha256": {
                "source.csv": hash_bytes(source_path.read_bytes()),
                "sink.csv": hash_bytes(sink_path.read_bytes()),
            },
            "command": command,
            "resource_telemetry": {
                "wall_seconds": wall_seconds,
                "memory_current_path": str(
                    self.settings.cgroup_memory_current_path
                ),
                "memory_stat_path": str(self.settings.cgroup_memory_stat_path),
                "memory_samples": memory_samples,
                "peak_memory_bytes": peak_memory_bytes,
                "peak_working_set_bytes": peak_working_set_bytes,
                "memory_limit_bytes": self.settings.memory_limit_bytes,
                "memory_limit_consecutive_samples": (
                    self.settings.memory_limit_consecutive_samples
                ),
                "maximum_consecutive_high_memory_samples": (
                    maximum_consecutive_high_memory_samples
                ),
                "resource_poll_seconds": self.settings.resource_poll_seconds,
                "resource_exhausted": resource_exhausted,
            },
            "artifacts": artifacts,
        }
        (job_dir / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        return RunResult(
            status=status,
            return_code=return_code,
            error=error,
            failure_code=failure_code,
        )

    @staticmethod
    def _read_memory_bytes(path: Path) -> int | None:
        try:
            value = int(path.read_text(encoding="ascii").strip())
        except (OSError, UnicodeError, ValueError):
            return None
        return value if value >= 0 else None

    @classmethod
    def _read_memory_sample(
        cls,
        current_path: Path,
        stat_path: Path,
    ) -> tuple[int | None, int | None]:
        current = cls._read_memory_bytes(current_path)
        if current is None:
            return None, None
        inactive_file: int | None = None
        try:
            for line in stat_path.read_text(encoding="ascii").splitlines():
                key, value = line.split(maxsplit=1)
                if key == "inactive_file":
                    inactive_file = int(value)
                    break
        except (OSError, UnicodeError, ValueError):
            inactive_file = None
        if inactive_file is None or inactive_file < 0:
            return current, None
        return current, max(0, current - inactive_file)

    @staticmethod
    def _terminate_process_group(process: subprocess.Popen[bytes]) -> None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait(timeout=10)
        except (ProcessLookupError, subprocess.TimeoutExpired):
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait(timeout=10)

    @staticmethod
    def _discover_artifacts(job_dir: Path) -> list[str]:
        artifacts: list[str] = []
        for filename in ("stdout.log", "stderr.log"):
            path = job_dir / filename
            if path.is_file():
                artifacts.append(filename)
        raw_dir = job_dir / "raw"
        if raw_dir.is_dir():
            artifacts.extend(
                str(path.relative_to(job_dir)).replace("\\", "/")
                for path in sorted(raw_dir.rglob("*"))
                if path.is_file()
            )
        return artifacts
