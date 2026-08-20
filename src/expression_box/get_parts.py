from __future__ import annotations

import re
import shutil
import subprocess
from typing import Any
from xml.etree import ElementTree
from langchain.tools import tool
import requests

from src.config.service_config import (
    IGEM_FALLBACK_PROCESS_TIMEOUT_SECONDS,
    IGEM_PART_URL,
    IGEM_REQUEST_TIMEOUT,
)
from src.runtime.monitor import monitor
from src.tools.common.session_paths import parts_dir as resolve_parts_dir


PART_ID_RE = re.compile(r"^BBa_[A-Za-z0-9_]+$", re.IGNORECASE)
# 保留原字段名，允许测试替换；唯一默认值和中文说明位于 service_config.py。
REQUEST_TIMEOUT = IGEM_REQUEST_TIMEOUT


def _normalize_part_id(part_id: str) -> str:
    normalized = str(part_id or "").strip()
    if normalized.lower().startswith("igem:"):
        normalized = normalized.split(":", 1)[1]
    return normalized


def _fetch_part_payload(part_id: str) -> bytes:
    url = f"{IGEM_PART_URL}{requests.utils.quote(part_id)}"
    try:
        response = requests.get(url, headers={"Accept": "application/xml"}, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        return response.content
    except Exception:
        return _fetch_part_payload_via_powershell(part_id)


def _fetch_part_payload_via_powershell(part_id: str) -> bytes:
    safe_part_id = part_id.replace("'", "''")
    https_url = f"{IGEM_PART_URL}{safe_part_id}"
    command = (
        "[Console]::OutputEncoding=[System.Text.Encoding]::UTF8; "
        f"(Invoke-WebRequest -UseBasicParsing -AllowInsecureRedirect '{https_url}').Content"
    )
    shell = shutil.which("pwsh") or shutil.which("powershell") or "powershell"
    completed = subprocess.run(
        [shell, "-NoProfile", "-Command", command],
        capture_output=True,
        text=True,
        timeout=IGEM_FALLBACK_PROCESS_TIMEOUT_SECONDS,
        check=False,
    )
    if completed.returncode != 0 or not completed.stdout.strip():
        stderr = completed.stderr.strip() if completed.stderr else "empty response"
        raise RuntimeError(f"PowerShell iGEM fetch failed: {stderr}")
    return completed.stdout.encode("utf-8")


def _text(element: ElementTree.Element | None, tag: str) -> str:
    if element is None:
        return ""
    child = element.find(tag)
    if child is None or child.text is None:
        return ""
    return child.text.strip()


def _parse_part_xml(payload: bytes) -> dict[str, Any]:
    root = ElementTree.fromstring(payload)
    part = root.find(".//part")
    if part is None:
        raise ValueError("iGEM XML response does not contain a <part> element.")

    parameters: dict[str, str] = {}
    for parameter in part.findall(".//parameter"):
        name = _text(parameter, "name")
        value = _text(parameter, "value")
        if name:
            parameters[name] = value

    sequence = _text(part.find("sequences"), "seq_data")
    normalized_sequence = "".join(sequence.split()).upper()
    return {
        "id": _text(part, "part_name"),
        "name": _text(part, "part_short_name") or _text(part, "part_name"),
        "part_short_desc": _text(part, "part_short_desc"),
        "sequence": normalized_sequence,
        "source": "iGEM Registry",
        "registry_metadata": {
            "registry_part_type": _text(part, "part_type"),
            "release_status": _text(part, "release_status"),
            "sample_status": _text(part, "sample_status"),
            "part_results": _text(part, "part_results"),
            "author": _text(part, "part_author"),
            "parameters": parameters,
        },
    }


def _infer_part_sequence_type(part: dict[str, Any]) -> str:
    registry_type = str((part.get("registry_metadata") or {}).get("registry_part_type") or "").strip().lower()
    if registry_type in {"promoter", "rbs"}:
        return registry_type
    if "terminator" in registry_type:
        return "terminator"
    return "part"

@tool
def get_parts_facts(*, part_id: str) -> dict[str, Any]:
    """
    下载并保存指定 iGEM part 的事实信息和序列。
    
    调用时机：推荐或选择 parts 时发现某个 part 缺少本地序列。
    输入：part_id。
    返回：part 元数据、sequence_file、length、hash 和保存路径。
    限制：只获取单个 part；不写 parts_selection，不组装表达盒。
    """

    tool_name = "get_parts_facts"
    monitor.report_start(tool_name, {"part_id": part_id})
    try:
        normalized_part_id = _normalize_part_id(part_id)
        if not normalized_part_id:
            raise ValueError("part_id 不能为空。")
        if not PART_ID_RE.fullmatch(normalized_part_id):
            raise ValueError("part_id 必须是合法的 BBa_* 标识。")

        monitor.report_running(tool_name, f"正在获取 iGEM part: {normalized_part_id}", progress=0.35)
        part = _parse_part_xml(_fetch_part_payload(normalized_part_id))
    except Exception:
        monitor.report_error(tool_name, f"未能从 iGEM 获取 {part_id} 的序列。")
        return {
            "success": False,
            "part_id": _normalize_part_id(part_id),
            "error": f"未能从 iGEM 获取 {_normalize_part_id(part_id)} 的序列。",
        }

    sequence = str(part.get("sequence") or "").strip().upper()
    if not sequence:
        return {
            "success": False,
            "part_id": normalized_part_id,
            "error": f"{normalized_part_id} 未返回可用序列。",
        }

    sequence_type = _infer_part_sequence_type(part)

    parts_path = resolve_parts_dir()
    parts_path.mkdir(parents=True, exist_ok=True)
    write_path = parts_path / f"{normalized_part_id}.txt"
    write_path.write_text(sequence, encoding="utf-8")
    monitor.report_end(tool_name, {"part_id": normalized_part_id, "sequence_file": str(write_path.resolve())})

    return {
        "success": True,
        "part_id": normalized_part_id,
        "name": str(part.get("name") or normalized_part_id).strip(),
        "part_short_desc": str(part.get("part_short_desc") or "").strip(),
        "sequence": sequence,
        "sequence_type": sequence_type,
    }


if __name__ == '__main__':
    id = 'BBa_R0085'
    print(get_parts_facts(part_id=id))
