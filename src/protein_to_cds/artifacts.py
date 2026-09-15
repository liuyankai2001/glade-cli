"""Rollback CDS artifacts until their corresponding manifest commit succeeds."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


def write_bytes_atomic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False, suffix=".tmp") as handle:
            temporary = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class ArtifactTransaction:
    def __init__(self, paths: Iterable[Path], *, lock_path: Path | None = None):
        self.paths = tuple(paths)
        self.backups = {}
        self.lock_path = lock_path
        self.lock_fd = None
        self.committed = False

    def __enter__(self):
        if self.lock_path is not None:
            self.lock_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                self.lock_fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError as exc:
                raise ValueError("当前项目正在生成或提交 CDS，请稍后重试") from exc
        try:
            self.backups = {path: path.read_bytes() if path.exists() else None for path in self.paths}
        except BaseException:
            self._unlock()
            raise
        return self

    def commit(self) -> None:
        self.committed = True

    def __exit__(self, exc_type, exc, traceback):
        try:
            if not self.committed:
                for path, content in self.backups.items():
                    current = path.read_bytes() if path.exists() else None
                    if current == content:
                        continue
                    if content is None:
                        path.unlink(missing_ok=True)
                    else:
                        write_bytes_atomic(path, content)
        finally:
            self._unlock()
        return False

    def _unlock(self):
        if self.lock_fd is not None:
            os.close(self.lock_fd)
            self.lock_fd = None
            self.lock_path.unlink(missing_ok=True)


def raw_cds_metadata(item: Mapping[str, Any], report: Mapping[str, Any]) -> dict[str, Any] | None:
    """Read the immutable baseline from a report, accepting legacy manifests."""
    metadata = report.get("raw_cds", item.get("raw_cds"))
    if isinstance(metadata, Mapping) and metadata.get("path"):
        return dict(metadata)
    metrics = report.get("raw")
    if not isinstance(metrics, Mapping) or not metrics.get("sequence_sha256"):
        return None
    return {
        "path": f"protein_to_cds/raw_cds/{item['accession']}.raw.fasta",
        "sequence_sha256": metrics["sequence_sha256"],
        "length_nt": metrics.get("length_nt"),
    }
