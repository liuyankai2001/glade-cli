from __future__ import annotations

import json
import tempfile
from collections.abc import Mapping
from functools import cache
from pathlib import Path
from typing import Any

from filelock import FileLock

SCHEMA_VERSION = "design_manifest.v1"


@cache
def _manifest_lock(path: Path) -> FileLock:
    return FileLock(str(path) + ".lock", timeout=30)


def manifest_update_lock(path: str | Path) -> FileLock:
    """Serialize revision checks and writes across CLI/API processes."""
    resolved = Path(path).expanduser().resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return _manifest_lock(resolved)


def read_design_manifest(path: str | Path) -> dict[str, Any]:
    """读取 design manifest；文件不存在时返回尚未写入的初始结构。"""

    manifest_path = Path(path).expanduser()
    if not manifest_path.exists():
        return {
            "schema_version": SCHEMA_VERSION,
            "revision": 0,
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
    manifest.setdefault("revision", 0)
    return manifest


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
    expected_revision: int | None = None,
) -> dict[str, Any]:
    """原子更新 manifest 的指定区段并自动递增 revision。"""

    with manifest_update_lock(path):
        return _update_locked(
            path,
            target_compound_id=target_compound_id,
            sections=sections,
            discard_sections=discard_sections,
            expected_revision=expected_revision,
        )


def _update_locked(
    path: str | Path,
    *,
    target_compound_id: str,
    sections: Mapping[str, Any],
    discard_sections: tuple[str, ...],
    expected_revision: int | None,
) -> dict[str, Any]:

    manifest_path = Path(path).expanduser()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if recorded_target and recorded_target != target_compound_id:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target}，"
            f"不能写入目标 {target_compound_id} 的结果"
        )
    try:
        current_revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"manifest revision 无效：{manifest.get('revision')!r}"
        ) from exc
    if current_revision < 0:
        raise ValueError(f"manifest revision 不能为负数：{current_revision}")
    if expected_revision is not None and current_revision != expected_revision:
        raise ValueError(
            "manifest revision 已变化："
            f"预期 {expected_revision}，当前 {current_revision}"
        )

    for section_name in discard_sections:
        if section_name not in sections:
            manifest.pop(section_name, None)
    manifest["schema_version"] = SCHEMA_VERSION
    manifest["target_compound_id"] = target_compound_id
    manifest.update(sections)
    manifest["revision"] = current_revision + 1
    _write_json_atomic(manifest_path, manifest)
    return manifest


__all__ = [
    "SCHEMA_VERSION",
    "manifest_update_lock",
    "read_design_manifest",
    "update_design_manifest",
]
