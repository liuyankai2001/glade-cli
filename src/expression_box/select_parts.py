from __future__ import annotations

import hashlib
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field

from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import design_manifest_file, outputs_dir as resolve_outputs_dir, parts_dir as resolve_parts_dir
from src.tools.expression_cassette_assembly_tools.part_identifiers import (
    normalize_part_id,
    part_sequence_path,
)


ALLOWED_ROLES = {"promoter", "rbs", "terminator"}


class SelectedPartArgs(BaseModel):
    cassette_index: int = Field(description="表达盒编号，从 1 开始")
    role: str = Field(description="元件角色：promoter、rbs 或 terminator")
    part_id: str = Field(description="BioBrick 或来源限定元件 ID，例如 BBa_R0085、KOSURI2013:R001")
    cds_id: str | None = Field(
        default=None,
        description="RBS 对应的优化 CDS cds_id；promoter 和 terminator 可不填",
    )
    name: str = Field(default="", description="可选；get_parts_facts 返回的元件名称")
    part_short_desc: str = Field(default="", description="可选；get_parts_facts 返回的元件简短描述")
    sequence_type: str = Field(default="", description="可选；get_parts_facts 返回的 sequence_type")


class SelectPartsArgs(BaseModel):
    selections: list[SelectedPartArgs] = Field(description="用户选择的表达元件列表")
    expected_revision: int | None = Field(default=None, description="可选，用于防止覆盖旧版本")


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"design_manifest.json root must be a JSON object: {path}")
    return data


def _write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
        suffix=".tmp",
    ) as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        tmp_path = Path(handle.name)
    tmp_path.replace(path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    return "".join(lines).upper()


def _relpath_or_abs(base_dir: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _normalize_part_id(part_id: str) -> str:
    return normalize_part_id(part_id)


def _normalize_role(role: str) -> str:
    normalized = str(role or "").strip().lower()
    if normalized not in ALLOWED_ROLES:
        raise ValueError(f"role must be one of {sorted(ALLOWED_ROLES)}: {role}")
    return normalized


def _cassette_map(manifest: dict[str, Any]) -> dict[int, set[str]]:
    assembly = manifest.get("expression_cassette_assembly")
    if not isinstance(assembly, dict):
        raise ValueError('design_manifest.json missing "expression_cassette_assembly"')

    cassettes = assembly.get("cassettes")
    if not isinstance(cassettes, list) or not cassettes:
        raise ValueError('design_manifest.json field "expression_cassette_assembly.cassettes" must be a non-empty list')

    result: dict[int, set[str]] = {}
    for cassette in cassettes:
        if not isinstance(cassette, dict):
            continue
        cassette_index = int(cassette.get("cassette_index") or 0)
        cds_ids = cassette.get("cds_ids")
        if cassette_index <= 0 or not isinstance(cds_ids, list):
            continue
        result[cassette_index] = {str(cds_id).strip() for cds_id in cds_ids if str(cds_id).strip()}

    if not result:
        raise ValueError("expression_cassette_assembly has no valid cassettes")
    return result


def _part_sequence_path(parts_dir: Path, part_id: str) -> Path:
    return part_sequence_path(parts_dir, part_id)


def _part_payload(
    *,
    output_dir: Path,
    parts_dir: Path,
    selection: SelectedPartArgs,
    cassette_cds_ids: set[str],
) -> dict[str, Any]:
    part_id = _normalize_part_id(selection.part_id)
    role = _normalize_role(selection.role)
    cds_id = str(selection.cds_id or "").strip() or None

    if role == "rbs":
        if not cds_id:
            raise ValueError(f"RBS part {part_id} must provide cds_id")
        if cds_id not in cassette_cds_ids:
            raise ValueError(f"cds_id {cds_id} is not in cassette {selection.cassette_index}")
    elif cds_id and cds_id not in cassette_cds_ids:
        raise ValueError(f"cds_id {cds_id} is not in cassette {selection.cassette_index}")

    sequence_path = _part_sequence_path(parts_dir, part_id)
    if not sequence_path.exists():
        raise FileNotFoundError(
            f"part sequence file not found: {sequence_path}. Call get_parts_facts first."
        )

    sequence = _read_sequence(sequence_path)
    if not sequence:
        raise ValueError(f"part sequence file is empty: {sequence_path}")

    payload = {
        "cassette_index": selection.cassette_index,
        "cds_id": cds_id,
        "role": role,
        "part_id": part_id,
        "sequence_file": {
            "path": _relpath_or_abs(output_dir, sequence_path),
            "sha256": _sha256_file(sequence_path),
            "length_bp": len(sequence),
        },
    }
    if selection.name:
        payload["name"] = selection.name
    if selection.part_short_desc:
        payload["part_short_desc"] = selection.part_short_desc
    if selection.sequence_type:
        payload["sequence_type"] = selection.sequence_type
    return payload


def _validate_required_roles(selected_parts: list[dict[str, Any]], cassette_map: dict[int, set[str]]) -> None:
    for cassette_index, cds_ids in cassette_map.items():
        cassette_parts = [
            part
            for part in selected_parts
            if int(part.get("cassette_index") or 0) == cassette_index
        ]
        roles = {str(part.get("role") or "") for part in cassette_parts}
        missing_roles = {"promoter", "terminator"} - roles
        if missing_roles:
            raise ValueError(f"cassette {cassette_index} missing parts: {sorted(missing_roles)}")

        rbs_cds_ids = {
            str(part.get("cds_id") or "")
            for part in cassette_parts
            if part.get("role") == "rbs"
        }
        missing_rbs = sorted(cds_ids - rbs_cds_ids)
        if missing_rbs:
            raise ValueError(f"cassette {cassette_index} missing RBS for cds_id: {missing_rbs}")


@tool(args_schema=SelectPartsArgs)
def select_parts(
    selections: list[SelectedPartArgs],
    expected_revision: int | None = None,
) -> str:
    """
    把用户确认的 promoter/RBS/terminator 选择写入 design_manifest.json。
    
    调用时机：recommend_expression_parts 返回候选后，用户确认具体元件组合。
    输入：selections，按 cassette_index、cds_id、role、part_id 指定元件。
    返回：ok、parts_selection 摘要、缺失序列提示和 manifest revision。
    限制：不下载 parts 序列；缺序列时先调用 get_parts_facts。
    """

    tool_name = "select_parts"
    monitor.report_start(tool_name, {"selection_count": len(selections) if selections else 0})
    try:
        output_path = resolve_outputs_dir()
        manifest_path = design_manifest_file()
        manifest = _read_manifest(manifest_path)

        current_revision = int(manifest.get("revision", 0))
        if expected_revision is not None and int(expected_revision) != current_revision:
            raise ValueError(
                f"manifest revision mismatch: expected {expected_revision}, current {current_revision}"
            )

        if not selections:
            raise ValueError("selections must be a non-empty list")

        cassette_map = _cassette_map(manifest)
        resolved_parts_dir = resolve_parts_dir()

        selected_parts = []
        unique_keys = set()
        for selection in selections:
            cassette_index = int(selection.cassette_index)
            if cassette_index not in cassette_map:
                raise ValueError(f"cassette_index not found in expression_cassette_assembly: {cassette_index}")

            payload = _part_payload(
                output_dir=output_path,
                parts_dir=resolved_parts_dir,
                selection=selection,
                cassette_cds_ids=cassette_map[cassette_index],
            )
            unique_key = (
                payload["cassette_index"],
                payload["role"],
                payload.get("cds_id") or "",
            )
            if unique_key in unique_keys:
                raise ValueError(
                    "duplicate part selection for cassette/role/cds_id: "
                    f"{unique_key}"
                )
            unique_keys.add(unique_key)
            selected_parts.append(payload)

        _validate_required_roles(selected_parts, cassette_map)

        manifest["parts_selection"] = {
            "status": "selected",
            "source": "select_parts",
            "selected_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "parts_dir": _relpath_or_abs(output_path, resolved_parts_dir),
            "selected_parts": selected_parts,
        }
        clear_downstream_fields(manifest, "parts_selection")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_path, manifest)

        monitor.report_end(tool_name, {"selected_part_count": len(selected_parts), "revision": manifest["revision"]})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": str(manifest_path),
                "revision": manifest["revision"],
                "selected_part_count": len(selected_parts),
                "selected_parts": selected_parts,
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(select_parts.invoke({
        "selections": [
            {"cassette_index": 1, "role": "promoter", "part_id": "BBa_R0085"},
            {"cassette_index": 1, "role": "rbs", "part_id": "BBa_B0034", "cds_id": "P21685_optimized"},
            {"cassette_index": 1, "role": "terminator", "part_id": "BBa_B0015"},
        ],
    }))
