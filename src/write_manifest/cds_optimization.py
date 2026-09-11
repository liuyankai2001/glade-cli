"""Commit CDS files and manifest together, rolling files back on write failure."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any

from src.write_manifest.store import read_design_manifest, update_design_manifest


def _stage(path: Path, content: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".gc-", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        return temporary
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def commit_cds_optimization(
    *, manifest_path: Path, project_root: Path, target: str, revision: int,
    selection: dict[str, Any], discard_sections: tuple[str, ...],
    files: dict[Path, bytes], guards: dict[Path, bytes | None],
) -> dict[str, Any]:
    """Stage outputs, check input snapshots, install files, then update manifest."""
    root = project_root.resolve()
    if any(not path.resolve().is_relative_to(root) for path in (*files, *guards, manifest_path)):
        raise ValueError("CDS 文件路径超出当前项目输出目录")
    lock_path = root / "protein_to_cds" / ".gc_optimization.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise ValueError("当前项目正在提交 CDS 优化结果，请稍后重试") from exc
    staged: dict[Path, Path] = {}
    backups: dict[Path, bytes | None] = {}
    installed: list[Path] = []
    try:
        for path, content in files.items():
            backups[path] = path.read_bytes() if path.exists() else None
            staged[path] = _stage(path, content)
        current = read_design_manifest(manifest_path)
        if int(current.get("revision", 0)) != revision:
            raise ValueError("优化期间 manifest revision 已变化，请重试")
        for path, expected in guards.items():
            observed = path.read_bytes() if path.exists() else None
            if observed != expected:
                raise ValueError(f"优化期间来源或目标文件已变化：{path.name}")
        for path, temporary in staged.items():
            os.replace(temporary, path)
            installed.append(path)
        return update_design_manifest(
            manifest_path, target_compound_id=target,
            sections={"cds_selection": selection}, discard_sections=discard_sections,
            expected_revision=revision,
        )
    except BaseException:
        for path in reversed(installed):
            previous = backups[path]
            if previous is None:
                path.unlink(missing_ok=True)
            else:
                restore = _stage(path, previous)
                try:
                    os.replace(restore, path)
                finally:
                    restore.unlink(missing_ok=True)
        raise
    finally:
        try:
            for temporary in staged.values():
                temporary.unlink(missing_ok=True)
        finally:
            os.close(lock_fd)
            lock_path.unlink(missing_ok=True)


__all__ = ["commit_cds_optimization"]
