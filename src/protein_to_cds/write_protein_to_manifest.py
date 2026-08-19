from __future__ import annotations

import hashlib
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field
from Bio.Seq import Seq
from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import (
    design_manifest_file,
    protein_optimized_sequence_dir,
    protein_sequence_dir,
    session_dir as resolve_session_dir,
)


class WriteSelectedOptimizedCdsArgs(BaseModel):
    cds_id: str = Field(description="CDS 选择 ID，例如 step_1_P12345_ecoli")
    protein_accession: str = Field(default="", description="来源蛋白 accession")
    protein_name: str = Field(default="", description="来源蛋白名称")
    organism_name: str = Field(default="", description="来源蛋白物种名称")
    step_index: int | None = Field(default=None, description="对应 pathway step_index")
    reaction_id: str = Field(default="", description="对应 KEGG reaction ID")
    ec_number: str = Field(default="", description="对应 EC 编号")
    host_name: str = Field(default="", description="密码子优化宿主名称")
    host_organism_id: int | None = Field(default=None, description="CodonTransformer 宿主 organism id")
    expected_revision: int | None = Field(default=None, description="可选，用于防止覆盖旧版本")


def _session_root() -> Path:
    return resolve_session_dir()


def _safe_path(session_root: Path, path_value: str | Path) -> Path:
    path = Path(path_value)
    if not path.is_absolute():
        path = session_root / path
    path = path.resolve()

    if path != session_root and session_root not in path.parents:
        raise ValueError(f"path escapes session_dir: {path}")
    return path


def _relpath(session_root: Path, path: Path) -> str:
    return path.resolve().relative_to(session_root.resolve()).as_posix()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _read_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    return "".join(lines).upper()


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return value


def _translated_protein(sequence: str) -> str:
    if len(sequence) % 3:
        return ""
    return str(Seq(sequence).translate(table=11, to_stop=False)).rstrip("*")


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


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False)


@tool(args_schema=WriteSelectedOptimizedCdsArgs)
def write_selected_optimized_cds_to_manifest(
    cds_id: str,
    protein_accession: str = "",
    protein_name: str = "",
    organism_name: str = "",
    step_index: int | None = None,
    reaction_id: str = "",
    ec_number: str = "",
    host_name: str = "",
    host_organism_id: int | None = None,
    expected_revision: int | None = None,
) -> str:
    """
    把用户确认的优化 CDS 写入 design_manifest.json 的 cds_selection。
    
    调用时机：search_protein_sequence 和 codon_optimization 产物经用户确认后。
    输入：accession、optimized_cds_path 或序列、宿主和步骤信息。
    返回：ok、selected_cds 摘要、文件 hash、manifest_path 和 revision。
    限制：不搜索蛋白、不运行优化模型；只登记已确认 CDS。
    """

    tool_name = "write_selected_optimized_cds_to_manifest"
    monitor.report_start(tool_name, {"cds_id": cds_id, "protein_accession": protein_accession})
    try:
        session_root = _session_root()
        optimized_cds_path = protein_optimized_sequence_dir() / f"{cds_id}_optimized.txt"
        if not optimized_cds_path.exists() and protein_accession:
            optimized_cds_path = protein_optimized_sequence_dir() / f"{protein_accession}_optimized.txt"
        cds_path = _safe_path(session_root, optimized_cds_path)
        if not cds_path.exists():
            raise FileNotFoundError(cds_path)

        sequence = _read_sequence(cds_path)
        if not sequence:
            raise ValueError("optimized CDS sequence is empty")

        report_path = protein_optimized_sequence_dir() / f"{cds_id}_optimization_report.json"
        if not report_path.exists() and protein_accession:
            report_path = protein_optimized_sequence_dir() / f"{protein_accession}_optimization_report.json"
        report_path = _safe_path(session_root, report_path)
        if not report_path.exists():
            raise ValueError(
                "constraint optimization report is missing; rerun codon_optimization before writing the CDS"
            )
        report = _load_json(report_path)
        final_report = report.get("final") if isinstance(report.get("final"), dict) else {}
        if report.get("status") != "PASS" or final_report.get("gate_status") != "PASS":
            raise ValueError("constraint optimization report did not pass the CDS gate")
        sequence_sha256 = _sha256_text(sequence)
        if final_report.get("sequence_sha256") != sequence_sha256:
            raise ValueError("optimized CDS hash does not match the constraint optimization report")

        raw_report = report.get("initial") if isinstance(report.get("initial"), dict) else {}
        raw_path_value = raw_report.get("sequence_path")
        if not raw_path_value:
            raise ValueError("constraint optimization report is missing the raw CDS path")
        raw_path = _safe_path(session_root, raw_path_value)
        if not raw_path.exists():
            raise FileNotFoundError(raw_path)
        raw_sequence = _read_sequence(raw_path)
        if raw_report.get("sequence_sha256") != _sha256_text(raw_sequence):
            raise ValueError("raw CDS hash does not match the constraint optimization report")

        if not protein_accession:
            raise ValueError("protein_accession is required for independent CDS translation validation")
        protein_path = _safe_path(session_root, protein_sequence_dir() / f"{protein_accession}.fasta")
        if not protein_path.exists():
            raise FileNotFoundError(protein_path)
        protein_sequence = _read_sequence(protein_path).rstrip("*")
        if _translated_protein(sequence) != protein_sequence:
            raise ValueError("optimized CDS does not exactly recover the selected protein")
        report_host = report.get("host") if isinstance(report.get("host"), dict) else {}
        if host_organism_id is not None and report_host.get("organism_id") != host_organism_id:
            raise ValueError("host organism id does not match the constraint optimization report")

        manifest_path = design_manifest_file()

        if manifest_path.exists():
            with manifest_path.open("r", encoding="utf-8") as file:
                manifest = json.load(file)
        else:
            manifest = {
                "schema_version": "design_manifest.v1",
                "revision": 0,
            }

        if not isinstance(manifest, dict):
            raise ValueError("design_manifest.json root must be a JSON object")
        manifest.setdefault("schema_version", "design_manifest.v1")

        current_revision = int(manifest.get("revision", 0))
        if expected_revision is not None and expected_revision != current_revision:
            raise ValueError(
                f"manifest revision mismatch: expected {expected_revision}, current {current_revision}"
            )

        selected_item = {
            "cds_id": cds_id,
            "step_index": step_index,
            "reaction_id": reaction_id,
            "ec_number": ec_number,
            "protein": {
                "accession": protein_accession,
                "protein_name": protein_name,
                "organism_name": organism_name,
            },
            "optimized_cds": {
                "sequence_file": {
                    "path": _relpath(session_root, cds_path),
                    "sha256": _sha256_file(cds_path),
                },
                "sequence_sha256": sequence_sha256,
                "length_nt": len(sequence),
                "gc_percent": round(
                    (sequence.count("G") + sequence.count("C")) * 100 / len(sequence),
                    2,
                ),
                "start_codon": sequence[:3],
                "stop_codon": sequence[-3:],
                "length_multiple_of_3": len(sequence) % 3 == 0,
                "constraint_gate_status": "PASS",
                "constraint_policy_version": report.get("policy_version"),
                "optimization_report": {
                    "path": _relpath(session_root, report_path),
                    "sha256": _sha256_file(report_path),
                },
                "initial_cds": {
                    "sequence_file": {
                        "path": _relpath(session_root, raw_path),
                        "sha256": _sha256_file(raw_path),
                    },
                    "sequence_sha256": _sha256_text(raw_sequence),
                },
                "metrics": {
                    "initial": raw_report,
                    "final": final_report,
                    "changes": report.get("changes", {}),
                },
            },
        }

        cds_selection = manifest.setdefault("cds_selection", {})
        if not isinstance(cds_selection, dict):
            raise ValueError('design_manifest.json field "cds_selection" must be an object')

        cds_selection.update({
            "status": "selected",
            "source": "codon_optimization_v2_constraint_repair",
            "selected_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "host": {
                "name": host_name,
                "codon_transformer_organism_id": host_organism_id,
            },
        })

        existing_items = cds_selection.get("selected_cds", [])
        if not isinstance(existing_items, list):
            raise ValueError('design_manifest.json field "cds_selection.selected_cds" must be a list')

        next_items = [
            item
            for item in existing_items
            if not isinstance(item, dict) or item.get("cds_id") != cds_id
        ]
        next_items.append(selected_item)
        cds_selection["selected_cds"] = next_items

        clear_downstream_fields(manifest, "cds_selection")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_path, manifest)

        result = {
            "ok": True,
            "manifest_path": _relpath(session_root, manifest_path),
            "revision": manifest["revision"],
            "cds_id": cds_id,
            "sequence_file": selected_item["optimized_cds"]["sequence_file"]["path"],
        }
        monitor.report_end(tool_name, {"cds_id": cds_id, "revision": manifest["revision"]})
        return json.dumps(
            result,
            ensure_ascii=False,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    from src.runtime.context import set_session_context

    test_session_dir = (
        Path(__file__).resolve().parents[3]
        / "agent_workspace"
        / "users"
        / "admin"
        / "sessions"
        / "0e09e4680bb24cb4a364f3d4c6c316a5"
    )
    set_session_context(str(test_session_dir))

    result = write_selected_optimized_cds_to_manifest.invoke({
        "cds_id": "P21685_optimized",
        "protein_accession": "P21685",
        "protein_name": "",
        "organism_name": "",
        "step_index": None,
        "reaction_id": "",
        "ec_number": "",
        "host_name": "",
        "host_organism_id": None,
    })
    print(result)
