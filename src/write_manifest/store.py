from __future__ import annotations

import json
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "design_manifest.v1"


def read_design_manifest(path: str | Path) -> dict[str, Any]:
    """读取 design manifest；文件不存在时返回尚未写入的初始结构。"""

    manifest_path = Path(path).expanduser()
    if not manifest_path.exists():
        return {
            "schema_version": SCHEMA_VERSION,
            "revision": 0,
        }
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"manifest 不是有效 JSON：{manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise ValueError(f"manifest 根节点必须是 JSON 对象：{manifest_path}")

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
) -> dict[str, Any]:
    """原子更新 manifest 的指定区段并自动递增 revision。"""

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
    "read_design_manifest",
    "update_design_manifest",
]
