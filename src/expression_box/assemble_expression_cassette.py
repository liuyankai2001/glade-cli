from __future__ import annotations

import hashlib
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel

from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import (
    design_manifest_file,
    expression_cassettes_dir as resolve_expression_cassettes_dir,
    outputs_dir as resolve_outputs_dir,
)


class AssembleExpressionCassetteArgs(BaseModel):
    pass


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _default_expression_cassette_dir() -> Path:
    return resolve_expression_cassettes_dir()


def _resolve_expression_cassette_dir(expression_cassette_dir: str | None) -> Path:
    if not expression_cassette_dir:
        return _default_expression_cassette_dir()

    path = Path(expression_cassette_dir).expanduser().resolve()
    if path.name == "outputs":
        return resolve_expression_cassettes_dir(path)
    return path


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


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _read_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    sequence = "".join(lines).upper()
    if not sequence:
        raise ValueError(f"sequence file is empty: {path}")
    return sequence


def _escape_qualifier(value: Any) -> str:
    return str(value or "").replace('"', "'")


def _genbank_feature_key(component_type: str) -> str:
    if component_type == "rbs":
        return "RBS"
    if component_type == "cds":
        return "CDS"
    return component_type


def _biopython_feature(component: dict[str, Any]) -> Any:
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    feature_type = _genbank_feature_key(str(component.get("type") or "misc_feature"))
    if len(feature_type) > 15:
        feature_type = "misc_feature"

    start = int(component.get("start_bp") or 0)
    end = int(component.get("end_bp") or 0)
    if start <= 0 or end < start:
        raise ValueError(f"invalid GenBank feature coordinates: {component}")

    source_id = _escape_qualifier(component.get("source_id"))
    label = _escape_qualifier(component.get("label") or source_id or feature_type)
    note = _escape_qualifier(component.get("note") or source_id or label)
    qualifiers: dict[str, list[str]] = {
        "label": [label],
        "note": [note],
    }
    if component.get("type") == "cds":
        product = _escape_qualifier(component.get("protein_name") or source_id or label)
        qualifiers["product"] = [product]
        qualifiers["codon_start"] = ["1"]

    return SeqFeature(
        SimpleLocation(start - 1, end),
        type=feature_type,
        qualifiers=qualifiers,
    )


def _write_genbank(path: Path, record_name: str, sequence: str, components: list[dict[str, Any]]) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(
        Seq(sequence),
        id=record_name[:16],
        name=record_name[:16],
        description=f"{record_name}.",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["data_file_division"] = "SYN"
    record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
    record.annotations["source"] = "synthetic DNA construct"
    record.annotations["organism"] = "synthetic DNA construct"
    record.features = [_biopython_feature(component) for component in components]
    SeqIO.write(record, str(path), "genbank")


def _relpath_or_abs(base_dir: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _resolve_sequence_path(output_dir: Path, path_value: Any) -> Path:
    if not path_value:
        raise ValueError("sequence file path is empty")

    path = Path(str(path_value))
    if path.is_absolute():
        return path.resolve()

    output_relative = (output_dir / path).resolve()
    if output_relative.exists():
        return output_relative

    session_root = output_dir.parent
    session_relative = (session_root / path).resolve()
    if session_relative.exists():
        return session_relative

    return output_relative


def _cds_by_id(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    cds_selection = manifest.get("cds_selection")
    if not isinstance(cds_selection, dict):
        raise ValueError('design_manifest.json missing "cds_selection"')

    selected_cds = cds_selection.get("selected_cds")
    if not isinstance(selected_cds, list) or not selected_cds:
        raise ValueError('design_manifest.json field "cds_selection.selected_cds" must be a non-empty list')

    result = {}
    for item in selected_cds:
        if not isinstance(item, dict):
            continue
        cds_id = str(item.get("cds_id") or "").strip()
        if cds_id:
            result[cds_id] = item
    if not result:
        raise ValueError("selected_cds has no valid cds_id")
    return result


def _assembly_cassettes(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    assembly = manifest.get("expression_cassette_assembly")
    if not isinstance(assembly, dict):
        raise ValueError('design_manifest.json missing "expression_cassette_assembly"')

    cassettes = assembly.get("cassettes")
    if not isinstance(cassettes, list) or not cassettes:
        raise ValueError('design_manifest.json field "expression_cassette_assembly.cassettes" must be a non-empty list')

    valid_cassettes = [cassette for cassette in cassettes if isinstance(cassette, dict)]
    valid_cassettes.sort(key=lambda item: int(item.get("cassette_index") or 0))
    return valid_cassettes


def _selected_parts(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    parts_selection = manifest.get("parts_selection")
    if not isinstance(parts_selection, dict):
        raise ValueError('design_manifest.json missing "parts_selection"')

    selected_parts = parts_selection.get("selected_parts")
    if not isinstance(selected_parts, list) or not selected_parts:
        raise ValueError('design_manifest.json field "parts_selection.selected_parts" must be a non-empty list')
    return [part for part in selected_parts if isinstance(part, dict)]


def _part_for(
    parts: list[dict[str, Any]],
    *,
    cassette_index: int,
    role: str,
    cds_id: str | None = None,
) -> dict[str, Any]:
    matches = [
        part
        for part in parts
        if int(part.get("cassette_index") or 0) == cassette_index
        and str(part.get("role") or "").strip().lower() == role
        and (cds_id is None or str(part.get("cds_id") or "").strip() == cds_id)
    ]
    if len(matches) != 1:
        target = f"cassette {cassette_index} {role}"
        if cds_id is not None:
            target += f" for {cds_id}"
        raise ValueError(f"expected exactly one part for {target}, found {len(matches)}")
    return matches[0]


def _sequence_file_path(item: dict[str, Any]) -> Any:
    sequence_file = item.get("sequence_file")
    if not isinstance(sequence_file, dict):
        raise ValueError(f"missing sequence_file in item: {item}")
    return sequence_file.get("path")


def _cds_sequence_file_path(item: dict[str, Any]) -> Any:
    optimized_cds = item.get("optimized_cds")
    if not isinstance(optimized_cds, dict):
        raise ValueError(f"missing optimized_cds for cds_id: {item.get('cds_id')}")
    sequence_file = optimized_cds.get("sequence_file")
    if not isinstance(sequence_file, dict):
        raise ValueError(f"missing optimized_cds.sequence_file for cds_id: {item.get('cds_id')}")
    return sequence_file.get("path")


def _component_payload(
    *,
    output_dir: Path,
    component_type: str,
    source_id: str,
    path: Path,
    sequence: str,
    start_bp: int,
    end_bp: int,
    cds_id: str | None = None,
    label: str = "",
    note: str = "",
    protein_name: str = "",
) -> dict[str, Any]:
    payload = {
        "type": component_type,
        "source_id": source_id,
        "start_bp": start_bp,
        "end_bp": end_bp,
        "label": label or source_id,
        "note": note or source_id,
        "sequence_file": {
            "path": _relpath_or_abs(output_dir, path),
            "sha256": _sha256_file(path),
            "length_bp": len(sequence),
        },
    }
    if cds_id:
        payload["cds_id"] = cds_id
    if protein_name:
        payload["protein_name"] = protein_name
    return payload


def _protein_name(cds_item: dict[str, Any], cds_id: str) -> str:
    protein = cds_item.get("protein")
    if not isinstance(protein, dict):
        return cds_id
    return (
        str(protein.get("protein_name") or "").strip()
        or str(protein.get("accession") or "").strip()
        or cds_id
    )


def _assemble_one_cassette(
    *,
    output_dir: Path,
    cassette: dict[str, Any],
    cds_items: dict[str, dict[str, Any]],
    selected_parts: list[dict[str, Any]],
) -> tuple[str, list[dict[str, Any]]]:
    cassette_index = int(cassette.get("cassette_index") or 0)
    cds_ids = cassette.get("cds_ids")
    if cassette_index <= 0 or not isinstance(cds_ids, list) or not cds_ids:
        raise ValueError(f"invalid cassette entry: {cassette}")

    components = []
    sequence_chunks = []
    cursor = 1

    promoter = _part_for(selected_parts, cassette_index=cassette_index, role="promoter")
    promoter_path = _resolve_sequence_path(output_dir, _sequence_file_path(promoter))
    promoter_sequence = _read_sequence(promoter_path)
    sequence_chunks.append(promoter_sequence)
    promoter_start = cursor
    promoter_end = cursor + len(promoter_sequence) - 1
    cursor = promoter_end + 1
    promoter_id = str(promoter.get("part_id") or "")
    components.append(_component_payload(
        output_dir=output_dir,
        component_type="promoter",
        source_id=promoter_id,
        path=promoter_path,
        sequence=promoter_sequence,
        start_bp=promoter_start,
        end_bp=promoter_end,
        label=promoter_id,
        note=f"promoter {promoter_id}",
    ))

    for cds_id_value in cds_ids:
        cds_id = str(cds_id_value or "").strip()
        if cds_id not in cds_items:
            raise ValueError(f"cds_id not found in selected_cds: {cds_id}")

        rbs = _part_for(selected_parts, cassette_index=cassette_index, role="rbs", cds_id=cds_id)
        rbs_path = _resolve_sequence_path(output_dir, _sequence_file_path(rbs))
        rbs_sequence = _read_sequence(rbs_path)
        sequence_chunks.append(rbs_sequence)
        rbs_start = cursor
        rbs_end = cursor + len(rbs_sequence) - 1
        cursor = rbs_end + 1
        rbs_id = str(rbs.get("part_id") or "")
        components.append(_component_payload(
            output_dir=output_dir,
            component_type="rbs",
            source_id=rbs_id,
            path=rbs_path,
            sequence=rbs_sequence,
            start_bp=rbs_start,
            end_bp=rbs_end,
            cds_id=cds_id,
            label=rbs_id,
            note=f"RBS {rbs_id} for {cds_id}",
        ))

        cds_item = cds_items[cds_id]
        cds_path = _resolve_sequence_path(output_dir, _cds_sequence_file_path(cds_item))
        cds_sequence = _read_sequence(cds_path)
        sequence_chunks.append(cds_sequence)
        cds_start = cursor
        cds_end = cursor + len(cds_sequence) - 1
        cursor = cds_end + 1
        protein_name = _protein_name(cds_item, cds_id)
        components.append(_component_payload(
            output_dir=output_dir,
            component_type="cds",
            source_id=cds_id,
            path=cds_path,
            sequence=cds_sequence,
            start_bp=cds_start,
            end_bp=cds_end,
            cds_id=cds_id,
            label=cds_id,
            note=f"CDS {cds_id}; protein {protein_name}",
            protein_name=protein_name,
        ))

    terminator = _part_for(selected_parts, cassette_index=cassette_index, role="terminator")
    terminator_path = _resolve_sequence_path(output_dir, _sequence_file_path(terminator))
    terminator_sequence = _read_sequence(terminator_path)
    sequence_chunks.append(terminator_sequence)
    terminator_start = cursor
    terminator_end = cursor + len(terminator_sequence) - 1
    terminator_id = str(terminator.get("part_id") or "")
    components.append(_component_payload(
        output_dir=output_dir,
        component_type="terminator",
        source_id=terminator_id,
        path=terminator_path,
        sequence=terminator_sequence,
        start_bp=terminator_start,
        end_bp=terminator_end,
        label=terminator_id,
        note=f"terminator {terminator_id}",
    ))

    return "".join(sequence_chunks), components


def assemble_cassette_sequence(
    *,
    output_dir: Path,
    cassette: dict[str, Any],
    cds_items: dict[str, dict[str, Any]],
    selected_parts: list[dict[str, Any]],
) -> tuple[str, list[dict[str, Any]]]:
    """Build one cassette with the exact sequence rules used by the file exporter."""

    return _assemble_one_cassette(
        output_dir=output_dir,
        cassette=cassette,
        cds_items=cds_items,
        selected_parts=selected_parts,
    )


@tool(args_schema=AssembleExpressionCassetteArgs)
def assemble_expression_cassette() -> str:
    """
    根据已选择的 CDS 和表达元件组装表达盒 GenBank 文件。
    
    调用时机：expression_cassette_assembly 和 parts_selection 均已准备好。
    返回：ok、表达盒文件列表、长度、hash、manifest revision 和 warnings。
    限制：只写 assembled_expression_cassettes 相关输出；不选择 parts、不做最终质粒组装。
    """

    tool_name = "assemble_expression_cassette"
    monitor.report_start(tool_name)
    try:
        cassette_dir = _resolve_expression_cassette_dir(None)
        output_dir = resolve_outputs_dir()
        manifest_path = design_manifest_file()
        monitor.report_running(tool_name, "正在读取表达盒组装所需 manifest 数据...", progress=0.2)
        manifest = _read_manifest(manifest_path)

        cds_items = _cds_by_id(manifest)
        cassettes = _assembly_cassettes(manifest)
        selected_parts = _selected_parts(manifest)

        cassette_records = []
        for cassette in cassettes:
            cassette_index = int(cassette.get("cassette_index") or 0)
            monitor.report_running(
                tool_name,
                f"正在组装表达盒 {cassette_index}...",
                data={"cassette_index": cassette_index},
            )
            cassette_sequence, components = _assemble_one_cassette(
                output_dir=output_dir,
                cassette=cassette,
                cds_items=cds_items,
                selected_parts=selected_parts,
            )
            genbank_path = cassette_dir / f"cassette_{cassette_index}.gb"
            _write_genbank(
                genbank_path,
                record_name=f"expression_cassette_{cassette_index}",
                sequence=cassette_sequence,
                components=components,
            )
            cassette_records.append({
                "cassette_index": cassette_index,
                "format": "genbank",
                "path": _relpath_or_abs(output_dir, genbank_path),
                "sha256": _sha256_file(genbank_path),
                "sequence_sha256": _sha256_text(cassette_sequence),
                "length_bp": len(cassette_sequence),
                "component_count": len(components),
                "components": components,
            })

        current_revision = int(manifest.get("revision", 0))
        manifest["assembled_expression_cassettes"] = {
            "status": "assembled",
            "source": "assemble_expression_cassette",
            "assembled_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "cassette_count": len(cassette_records),
            "cassette_files": cassette_records,
        }
        clear_downstream_fields(manifest, "assembled_expression_cassettes")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_path, manifest)

        monitor.report_end(tool_name, {"cassette_count": len(cassette_records), "cassette_dir": str(cassette_dir)})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": str(manifest_path),
                "revision": manifest["revision"],
                "cassette_count": len(cassette_records),
                "cassette_files": cassette_records,
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    expression_cassette_dir = (
        Path(__file__).resolve().parents[3]
        / "agent_workspace"
        / "users"
        / "admin"
        / "sessions"
        / "0e09e4680bb24cb4a364f3d4c6c316a5"
        / "outputs"
        / "expression_cassettes"
    )
    print(assemble_expression_cassette.invoke({
        "expression_cassette_dir": str(expression_cassette_dir),
    }))
