from __future__ import annotations

import json
import tempfile
from collections.abc import Mapping
from functools import cache
from pathlib import Path
from typing import Any

from filelock import FileLock

SCHEMA_VERSION = "design_manifest.v1"

# Compatibility only: discard obsolete write counters without interpreting them.
_LEGACY_WRITE_METADATA = frozenset({
    "revision", "expected_revision", "manifest_revision",
    "source_manifest_revision", "artifact_manifest_revision",
})


def _without_legacy_write_metadata(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            key: _without_legacy_write_metadata(child)
            for key, child in value.items()
            if key not in _LEGACY_WRITE_METADATA
        }
    if isinstance(value, (list, tuple)):
        return [_without_legacy_write_metadata(child) for child in value]
    return value


@cache
def _manifest_lock(path: Path) -> FileLock:
    return FileLock(str(path) + ".lock", timeout=30)


def manifest_update_lock(path: str | Path) -> FileLock:
    """Serialize manifest writes across CLI/API processes."""
    resolved = Path(path).expanduser().resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return _manifest_lock(resolved)


def read_design_manifest(path: str | Path) -> dict[str, Any]:
    """读取 design manifest；文件不存在时返回尚未写入的初始结构。"""

    manifest_path = Path(path).expanduser()
    if not manifest_path.exists():
        return {
            "schema_version": SCHEMA_VERSION,
        }

    # Windows readers must close the file before an atomic replacement begins.
    with _manifest_lock(manifest_path.resolve()):
        return _read_existing_manifest(manifest_path)


def _read_existing_manifest(manifest_path: Path) -> dict[str, Any]:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"manifest 不是有效 JSON：{manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise ValueError(f"manifest 根节点必须是 JSON 对象：{manifest_path}")  # noqa: TRY004 - invalid persisted JSON value

    schema_version = manifest.get("schema_version", SCHEMA_VERSION)
    if schema_version != SCHEMA_VERSION:
        raise ValueError(
            f"不支持的 manifest schema_version：{schema_version!r}；"
            f"当前仅支持 {SCHEMA_VERSION}"
        )
    manifest["schema_version"] = SCHEMA_VERSION
    return _without_legacy_write_metadata(manifest)


def _write_json_atomic(path: Path, data: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            delete=False,
            suffix=".tmp",
        ) as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def update_design_manifest(
    path: str | Path,
    *,
    target_compound_id: str,
    sections: Mapping[str, Any],
    discard_sections: tuple[str, ...] = (),
) -> dict[str, Any]:
    """原子更新 manifest 的指定区段。"""

    with manifest_update_lock(path):
        return _update_locked(
            path,
            target_compound_id=target_compound_id,
            sections=sections,
            discard_sections=discard_sections,
        )


def _update_locked(
    path: str | Path,
    *,
    target_compound_id: str,
    sections: Mapping[str, Any],
    discard_sections: tuple[str, ...],
) -> dict[str, Any]:

    manifest_path = Path(path).expanduser()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if recorded_target and recorded_target != target_compound_id:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target}，"
            f"不能写入目标 {target_compound_id} 的结果"
        )

    for section_name in discard_sections:
        if section_name not in sections:
            manifest.pop(section_name, None)
    manifest["schema_version"] = SCHEMA_VERSION
    manifest["target_compound_id"] = target_compound_id
    manifest.update(sections)
    manifest = _without_legacy_write_metadata(manifest)
    _write_json_atomic(manifest_path, manifest)
    return manifest


__all__ = [
    "SCHEMA_VERSION",
    "manifest_update_lock",
    "read_design_manifest",
    "update_design_manifest",
]
