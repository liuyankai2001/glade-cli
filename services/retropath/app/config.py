from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    rules_path: Path
    rules_version: str
    expected_rules_sha256: str
    knime_dir: Path
    knime_version: str
    rdkit_plugin_version: str
    openssl10_dir: Path
    workflow_version: str
    jobs_dir: Path
    max_queue: int
    job_timeout_seconds: int

    @classmethod
    def from_env(cls) -> "Settings":
        return cls(
            rules_path=Path(
                os.environ.get(
                    "RETROPATH_RULES_PATH",
                    "/opt/retropath/rules/rr02.csv",
                )
            ),
            rules_version=os.environ.get(
                "RETROPATH_RULES_VERSION", "rr02-rp2-hs"
            ),
            expected_rules_sha256=os.environ.get(
                "RETROPATH_RULES_SHA256",
                "e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f",
            ).lower(),
            knime_dir=Path(os.environ.get("RETROPATH_KNIME_DIR", "/opt/knime")),
            knime_version=os.environ.get("RETROPATH_KNIME_VERSION", "4.7.0"),
            rdkit_plugin_version=os.environ.get(
                "RETROPATH_RDKIT_PLUGIN_VERSION", "4.9.1"
            ),
            openssl10_dir=Path(
                os.environ.get(
                    "RETROPATH_OPENSSL10_DIR", "/home/mambauser/.openssl10/lib"
                )
            ),
            workflow_version=os.environ.get(
                "RETROPATH_WORKFLOW_VERSION", "r20260212"
            ),
            jobs_dir=Path(
                os.environ.get("RETROPATH_JOBS_DIR", "/var/lib/retropath/jobs")
            ),
            max_queue=int(os.environ.get("RETROPATH_MAX_QUEUE", "8")),
            job_timeout_seconds=int(
                os.environ.get("RETROPATH_JOB_TIMEOUT_SECONDS", "3600")
            ),
        )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_versions() -> tuple[str, str]:
    reported_version = metadata.version("retropath2_wrapper")
    conda_prefix = Path(os.environ.get("CONDA_PREFIX", "/opt/conda"))
    conda_meta = conda_prefix / "conda-meta"
    if conda_meta.is_dir():
        for path in sorted(conda_meta.glob("retropath2_wrapper-*.json")):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if payload.get("name") == "retropath2_wrapper" and payload.get("version"):
                return str(payload["version"]), reported_version
    return reported_version, reported_version


def workflow_path(workflow_version: str) -> Path:
    import retropath2_wrapper

    package_dir = Path(retropath2_wrapper.__file__).resolve().parent
    return package_dir / "workflows" / f"RetroPath2.0_{workflow_version}.knwf"


def find_knime_executable(knime_dir: Path) -> Path | None:
    if not knime_dir.is_dir():
        return None
    candidates = sorted(
        path
        for path in knime_dir.rglob("knime")
        if path.is_file() and os.access(path, os.X_OK)
    )
    return candidates[0] if candidates else None


@dataclass(frozen=True)
class RuntimeInfo:
    wrapper_version: str
    wrapper_reported_version: str
    workflow_path: Path
    knime_executable: Path | None
    rules_sha256: str | None
    errors: tuple[str, ...]

    @property
    def ready(self) -> bool:
        return not self.errors

    def public_dict(self, settings: Settings) -> dict[str, object]:
        return {
            "ready": self.ready,
            "wrapper_version": self.wrapper_version,
            "wrapper_reported_version": self.wrapper_reported_version,
            "workflow_version": settings.workflow_version,
            "knime_version": settings.knime_version,
            "rdkit_plugin_version": settings.rdkit_plugin_version,
            "rules_version": settings.rules_version,
            "rules_sha256": self.rules_sha256,
            "worker_concurrency": 1,
            "errors": list(self.errors),
        }


def inspect_runtime(settings: Settings) -> RuntimeInfo:
    errors: list[str] = []
    try:
        wrapper_version, wrapper_reported_version = package_versions()
    except metadata.PackageNotFoundError:
        wrapper_version = "missing"
        wrapper_reported_version = "missing"
        errors.append("retropath2_wrapper is not installed")

    expected_wrapper = "3.9.1"
    if wrapper_version != expected_wrapper:
        errors.append(
            f"retropath2_wrapper version mismatch: expected {expected_wrapper}, got {wrapper_version}"
        )

    workflow = workflow_path(settings.workflow_version)
    if not workflow.is_file():
        errors.append(f"workflow is missing: {workflow}")

    knime_executable = find_knime_executable(settings.knime_dir)
    if knime_executable is None:
        errors.append(f"KNIME executable is missing under {settings.knime_dir}")
    else:
        knime_root = knime_executable.parent
        if not list(
            (knime_root / "p2" / "plugins").glob("org.rdkit.knime.nodes_*.jar")
        ):
            errors.append("KNIME RDKit nodes plugin is missing")
        if not list(
            (knime_root / "p2" / "features").glob(
                "org.knime.features.chem.types_*"
            )
        ):
            errors.append("KNIME chemistry feature is missing")

    for library in ("libssl.so.10", "libcrypto.so.10"):
        if not (settings.openssl10_dir / library).is_file():
            errors.append(f"KNIME OpenSSL compatibility library is missing: {library}")

    rules_sha256: str | None = None
    if not settings.rules_path.is_file():
        errors.append(f"RR02 rules file is missing: {settings.rules_path}")
    else:
        rules_sha256 = sha256_file(settings.rules_path)
        if rules_sha256.lower() != settings.expected_rules_sha256:
            errors.append(
                "RR02 SHA-256 mismatch: "
                f"expected {settings.expected_rules_sha256}, got {rules_sha256}"
            )

    return RuntimeInfo(
        wrapper_version=wrapper_version,
        wrapper_reported_version=wrapper_reported_version,
        workflow_path=workflow,
        knime_executable=knime_executable,
        rules_sha256=rules_sha256,
        errors=tuple(errors),
    )
