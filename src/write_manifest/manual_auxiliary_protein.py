"""Import user-provided auxiliary amino-acid or CDS sequences."""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
import uuid
from collections.abc import Mapping
from io import StringIO
from pathlib import Path
from typing import Any, Literal, cast

from Bio import SeqIO

from src.pathway_analyze.target_id import validate_target_compound_id
from src.write_manifest.store import read_design_manifest, update_design_manifest


MANUAL_AUXILIARY_PROTEIN_SCHEMA_VERSION = (
    "manual_auxiliary_protein_selection.v1"
)
MANUAL_AUXILIARY_DOWNSTREAM_SECTIONS = (
    "protein_selection",
    "enzyme_system_selection",
    "cds_selection",
    "expression_box_selection",
    "expression_cassette_assembly",
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
)
SequenceType = Literal["protein", "cds"]


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少有效的 {field_name}")
    return value


def _normalize_sequence_type(value: Any) -> SequenceType:
    normalized = str(value or "").strip().lower()
    if normalized not in {"protein", "cds"}:
        raise ValueError("--sequence-type 必须是 protein 或 cds")
    return cast(SequenceType, normalized)


def _safe_protein_id(value: str, fallback_index: int) -> str:
    normalized = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    normalized = normalized.strip("._-").upper()
    return normalized or f"AUXILIARY_{fallback_index}"


def _parse_uploaded_sequences(path: Path) -> list[dict[str, str]]:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"无法读取辅助序列文件：{path}") from exc
    if not text.strip():
        raise ValueError(f"辅助序列文件为空：{path}")

    content = text.lstrip("\ufeff \t\r\n")
    parsed: list[dict[str, str]] = []
    if content.startswith(">"):
        try:
            records = list(SeqIO.parse(StringIO(content), "fasta"))
        except Exception as exc:
            raise ValueError(f"无法解析辅助序列 FASTA：{path}") from exc
        if not records:
            raise ValueError(f"辅助序列 FASTA 中没有记录：{path}")
        for index, record in enumerate(records, start=1):
            sequence = re.sub(r"\s+", "", str(record.seq)).upper()
            if not sequence:
                raise ValueError(f"辅助序列 FASTA 第 {index} 条记录为空：{path}")
            original_id = str(record.id or "").strip()
            parsed.append(
                {
                    "accession": _safe_protein_id(original_id, index),
                    "original_id": original_id,
                    "protein_name": str(record.description or original_id).strip(),
                    "sequence": sequence,
                }
            )
    else:
        sequence = re.sub(r"\s+", "", content).upper()
        if not sequence:
            raise ValueError(f"辅助序列文本中没有内容：{path}")
        original_id = path.stem
        parsed.append(
            {
                "accession": _safe_protein_id(original_id, 1),
                "original_id": original_id,
                "protein_name": original_id,
                "sequence": sequence,
            }
        )

    # The last record wins when multiple headers normalize to the same ID.
    deduplicated: dict[str, dict[str, str]] = {}
    for item in parsed:
        deduplicated[item["accession"]] = item
    return list(deduplicated.values())


def _canonical_fasta(accession: str, protein_name: str, sequence: str) -> str:
    label = re.sub(r"\s+", " ", protein_name).strip()
    header = accession if not label or label == accession else f"{accession} {label}"
    return f">{header}\n{sequence}\n"


def _write_text_atomic(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _relative(project_root: Path, path: Path) -> str:
    return path.resolve().relative_to(project_root.resolve()).as_posix()


def _read_snapshot_sequence(
    project_root: Path,
    protein: Mapping[str, Any],
) -> str | None:
    sequence_file = protein.get("sequence_file")
    if not isinstance(sequence_file, Mapping):
        return None
    relative_path = str(sequence_file.get("path") or "").strip()
    if not relative_path:
        return None
    path = (project_root / relative_path).resolve()
    try:
        path.relative_to(project_root.resolve())
    except ValueError:
        return None
    if not path.is_file():
        return None
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError):
        return None
    return "".join(
        re.sub(r"\s+", "", line)
        for line in lines
        if line.strip() and not line.lstrip().startswith(">")
    ).upper()


def _existing_manual_selection(
    value: Any,
) -> tuple[list[dict[str, Any]], list[str]]:
    if not isinstance(value, Mapping):
        return [], []
    if (
        str(value.get("schema_version") or "")
        != MANUAL_AUXILIARY_PROTEIN_SCHEMA_VERSION
    ):
        return [], []
    raw_proteins = value.get("proteins")
    proteins = (
        [dict(item) for item in raw_proteins if isinstance(item, Mapping)]
        if isinstance(raw_proteins, list)
        else []
    )
    source = value.get("source")
    raw_files = source.get("input_files") if isinstance(source, Mapping) else []
    input_files = (
        [str(item) for item in raw_files if str(item).strip()]
        if isinstance(raw_files, list)
        else []
    )
    return proteins, list(dict.fromkeys(input_files))


def _main_enzyme_context(manifest: Mapping[str, Any]) -> dict[str, Any]:
    selection = _mapping(
        manifest.get("main_enzyme_selection"),
        "main_enzyme_selection；请先执行 write --main-enzyme-set N",
    )
    status = str(selection.get("selection_status") or "").strip()
    if status not in {"user_selected", "user_selected_pending_review"}:
        raise ValueError("manifest 中的主酶组合尚未由用户确认")
    coverage = _mapping(selection.get("coverage"), "main_enzyme_selection.coverage")
    if coverage.get("complete") is not True:
        raise ValueError("主酶组合尚未完整覆盖路线")
    raw_proteins = selection.get("proteins")
    if not isinstance(raw_proteins, list) or not raw_proteins:
        raise ValueError("main_enzyme_selection.proteins 不能为空")
    step_indexes = sorted({
        int(step)
        for raw in raw_proteins
        if isinstance(raw, Mapping)
        for step in (raw.get("assigned_step_indexes") or [])
    })
    if not step_indexes:
        raise ValueError("主酶组合缺少 assigned_step_indexes")
    return {
        "selected_set_id": selection.get("selected_set_id"),
        "selected_set_fingerprint": str(
            selection.get("selected_set_fingerprint") or ""
        ),
        "assigned_step_indexes": step_indexes,
    }


def _payload(
    *,
    proteins: list[dict[str, Any]],
    input_files: list[str],
    main_context: Mapping[str, Any],
) -> dict[str, Any]:
    normalized_proteins = [
        {
            **protein,
            "roles": ["auxiliary_protein"],
            "assigned_step_indexes": list(
                main_context["assigned_step_indexes"]
            ),
            "required_by_main_accessions": [],
        }
        for protein in proteins
    ]
    return {
        "schema_version": MANUAL_AUXILIARY_PROTEIN_SCHEMA_VERSION,
        "selection_mode": "manual_upload",
        "selection_status": "user_selected",
        "can_advance": True,
        "source": {
            "type": "inputs_sequence_file",
            "input_files": input_files,
            "selected_set_id": main_context["selected_set_id"],
            "selected_set_fingerprint": main_context[
                "selected_set_fingerprint"
            ],
        },
        "auxiliary_proteins_to_introduce": [
            protein["accession"] for protein in normalized_proteins
        ],
        "proteins": normalized_proteins,
    }


def add_manual_auxiliary_protein(config: Any) -> dict[str, Any]:
    """Append one uploaded file to the current manual auxiliary selection."""

    target_compound = validate_target_compound_id(config.target_name)
    sequence_type = _normalize_sequence_type(
        getattr(config, "sequence_type", None)
    )
    filename = str(getattr(config, "protein_file", "") or "").strip()
    if not filename:
        raise ValueError("缺少 --protein-file")
    inputs_root = Path(config.inputs_dir).expanduser().resolve()
    input_path = (inputs_root / filename).resolve()
    try:
        relative_input = input_path.relative_to(inputs_root).as_posix()
    except ValueError as exc:
        raise ValueError("辅助序列文件必须位于 inputs 目录内") from exc
    if not input_path.is_file():
        raise FileNotFoundError(f"未找到辅助序列文件：{input_path}")

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if recorded_target != target_compound:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target or '空'}，"
            f"与当前输入目标 {target_compound} 不一致"
        )
    try:
        current_revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision 必须是整数") from exc
    main_context = _main_enzyme_context(manifest)
    uploaded = _parse_uploaded_sequences(input_path)

    existing, input_files = _existing_manual_selection(
        manifest.get("auxiliary_protein_selection")
    )
    by_accession = {
        str(item.get("accession") or "").strip().upper(): item
        for item in existing
        if str(item.get("accession") or "").strip()
    }
    added: list[str] = []
    replaced: list[str] = []
    unchanged: list[str] = []
    staged_files: list[tuple[Path, str]] = []
    snapshot_dir = (
        project_root
        / "protein_to_cds"
        / "uploaded_sequences"
        / f"manifest_revision_{current_revision + 1:06d}"
    )
    for item in uploaded:
        accession = item["accession"]
        previous = by_accession.get(accession)
        previous_sequence = (
            _read_snapshot_sequence(project_root, previous)
            if previous is not None
            else None
        )
        if (
            previous is not None
            and str(previous.get("sequence_type") or "") == sequence_type
            and str(previous.get("protein_name") or "") == item["protein_name"]
            and previous_sequence == item["sequence"]
        ):
            unchanged.append(accession)
            continue

        snapshot_path = snapshot_dir / f"{accession}.{sequence_type}.fasta"
        staged_files.append(
            (
                snapshot_path,
                _canonical_fasta(
                    accession,
                    item["protein_name"],
                    item["sequence"],
                ),
            )
        )
        by_accession[accession] = {
            "accession": accession,
            "original_id": item["original_id"],
            "protein_name": item["protein_name"],
            "organism_name": "",
            "roles": ["auxiliary_protein"],
            "assigned_step_indexes": list(
                main_context["assigned_step_indexes"]
            ),
            "required_by_main_accessions": [],
            "sequence_type": sequence_type,
            "source_type": "user_uploaded",
            "source_input_file": relative_input,
            "sequence_file": {
                "path": _relative(project_root, snapshot_path),
                "format": "fasta",
            },
        }
        (replaced if previous is not None else added).append(accession)

    if relative_input not in input_files:
        input_files.append(relative_input)
    proteins = list(by_accession.values())
    payload = _payload(
        proteins=proteins,
        input_files=input_files,
        main_context=main_context,
    )
    current_selection = manifest.get("auxiliary_protein_selection")
    changed = not (
        isinstance(current_selection, Mapping)
        and dict(current_selection) == payload
    )
    if changed:
        for path, fasta_text in staged_files:
            _write_text_atomic(path, fasta_text)
        updated = update_design_manifest(
            manifest_path,
            target_compound_id=target_compound,
            sections={"auxiliary_protein_selection": payload},
            discard_sections=MANUAL_AUXILIARY_DOWNSTREAM_SECTIONS,
            expected_revision=current_revision,
        )
    else:
        updated = manifest

    return {
        "运行成功": True,
        "目标化合物": target_compound,
        "序列类型": sequence_type,
        "输入文件": relative_input,
        "新增辅助蛋白": added,
        "替换辅助蛋白": replaced,
        "未变化辅助蛋白": unchanged,
        "当前辅助蛋白": [item["accession"] for item in proteins],
        "清单文件": str(manifest_path),
        "清单版本": updated["revision"],
        "清单是否更新": changed,
    }


def _requested_protein_ids(value: Any) -> list[str]:
    raw_values = [value] if isinstance(value, str) else list(value or [])
    normalized: list[str] = []
    for index, raw in enumerate(raw_values, start=1):
        text = str(raw or "").strip()
        if not text:
            raise ValueError("--protein-id 不能为空")
        normalized.append(_safe_protein_id(text, index))
    normalized = list(dict.fromkeys(normalized))
    if not normalized:
        raise ValueError("至少需要一个 --protein-id")
    return normalized


def _uploaded_sequence_root(project_root: Path) -> Path:
    return (
        project_root
        / "protein_to_cds"
        / "uploaded_sequences"
    ).resolve()


def _assert_within(path: Path, root: Path, label: str) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError(f"{label} 路径超出项目上传序列目录：{path}") from exc
    return resolved


def _snapshot_paths_for_removal(
    *,
    project_root: Path,
    proteins: list[Mapping[str, Any]],
    accessions: list[str],
) -> list[Path]:
    uploaded_root = _uploaded_sequence_root(project_root)
    targets: set[Path] = set()
    requested = set(accessions)

    for protein in proteins:
        accession = str(protein.get("accession") or "").strip().upper()
        if accession not in requested:
            continue
        sequence_file = protein.get("sequence_file")
        if not isinstance(sequence_file, Mapping):
            continue
        relative_path = str(sequence_file.get("path") or "").strip()
        if not relative_path:
            continue
        current_path = _assert_within(
            project_root / relative_path,
            uploaded_root,
            f"辅助蛋白 {accession}",
        )
        if current_path.is_file():
            targets.add(current_path)

    if uploaded_root.is_dir():
        for revision_dir in uploaded_root.iterdir():
            if (
                not revision_dir.is_dir()
                or not revision_dir.name.startswith("manifest_revision_")
            ):
                continue
            safe_revision_dir = _assert_within(
                revision_dir,
                uploaded_root,
                "历史上传目录",
            )
            for accession in accessions:
                for sequence_type in ("protein", "cds"):
                    candidate = _assert_within(
                        safe_revision_dir
                        / f"{accession}.{sequence_type}.fasta",
                        uploaded_root,
                        f"辅助蛋白 {accession}",
                    )
                    if candidate.is_file():
                        targets.add(candidate)
    return sorted(targets, key=lambda item: str(item).lower())


def _restore_quarantined_files(
    moved: list[tuple[Path, Path]],
    quarantine_root: Path,
) -> list[str]:
    errors: list[str] = []
    for original, quarantined in reversed(moved):
        if not quarantined.exists():
            continue
        try:
            original.parent.mkdir(parents=True, exist_ok=True)
            os.replace(quarantined, original)
        except OSError as exc:
            errors.append(
                f"{quarantined} -> {original}: {exc}"
            )
    if quarantine_root.exists() and not errors:
        shutil.rmtree(quarantine_root, ignore_errors=True)
    return errors


def _cleanup_empty_snapshot_directories(
    uploaded_root: Path,
    original_paths: list[Path],
) -> None:
    for directory in sorted(
        {path.parent for path in original_paths},
        key=lambda item: len(item.parts),
        reverse=True,
    ):
        if directory.parent != uploaded_root:
            continue
        try:
            directory.rmdir()
        except OSError:
            pass
    quarantine_parent = uploaded_root / ".delete_quarantine"
    try:
        quarantine_parent.rmdir()
    except OSError:
        pass


def remove_manual_auxiliary_proteins(config: Any) -> dict[str, Any]:
    """Delete selected manual auxiliary records and project snapshots."""

    target_compound = validate_target_compound_id(config.target_name)
    requested = _requested_protein_ids(getattr(config, "protein_id", None))
    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    recorded_target = str(manifest.get("target_compound_id") or "").strip()
    if recorded_target != target_compound:
        raise ValueError(
            f"manifest 目标化合物为 {recorded_target or '空'}，"
            f"与当前输入目标 {target_compound} 不一致"
        )
    selection = _mapping(
        manifest.get("auxiliary_protein_selection"),
        "auxiliary_protein_selection；当前没有手动导入的辅助蛋白",
    )
    if (
        str(selection.get("schema_version") or "")
        != MANUAL_AUXILIARY_PROTEIN_SCHEMA_VERSION
    ):
        raise ValueError("当前辅助蛋白选择不是手动导入格式，不能使用删除命令")
    raw_proteins = selection.get("proteins")
    if not isinstance(raw_proteins, list):
        raise ValueError("auxiliary_protein_selection.proteins 必须是列表")
    proteins = [
        dict(item) for item in raw_proteins if isinstance(item, Mapping)
    ]
    by_accession = {
        str(item.get("accession") or "").strip().upper(): item
        for item in proteins
        if str(item.get("accession") or "").strip()
    }
    missing = [accession for accession in requested if accession not in by_accession]
    if missing:
        available = sorted(by_accession)
        raise ValueError(
            f"未找到手动辅助蛋白：{missing}；当前可删除 ID：{available}"
        )
    try:
        current_revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision 必须是整数") from exc

    remaining = [
        item
        for item in proteins
        if str(item.get("accession") or "").strip().upper() not in set(requested)
    ]
    referenced_inputs = list(dict.fromkeys(
        str(item.get("source_input_file") or "").strip()
        for item in remaining
        if str(item.get("source_input_file") or "").strip()
    ))
    if remaining:
        payload = _payload(
            proteins=remaining,
            input_files=referenced_inputs,
            main_context=_main_enzyme_context(manifest),
        )
        sections: Mapping[str, Any] = {
            "auxiliary_protein_selection": payload,
        }
        discard_sections = MANUAL_AUXILIARY_DOWNSTREAM_SECTIONS
    else:
        sections = {}
        discard_sections = (
            "auxiliary_protein_selection",
            *MANUAL_AUXILIARY_DOWNSTREAM_SECTIONS,
        )

    snapshots = _snapshot_paths_for_removal(
        project_root=project_root,
        proteins=proteins,
        accessions=requested,
    )
    uploaded_root = _uploaded_sequence_root(project_root)
    quarantine_root = (
        uploaded_root
        / ".delete_quarantine"
        / uuid.uuid4().hex
    )
    moved: list[tuple[Path, Path]] = []
    try:
        for original in snapshots:
            relative_path = original.relative_to(uploaded_root)
            quarantined = quarantine_root / relative_path
            quarantined.parent.mkdir(parents=True, exist_ok=True)
            os.replace(original, quarantined)
            moved.append((original, quarantined))
    except OSError as exc:
        rollback_errors = _restore_quarantined_files(moved, quarantine_root)
        detail = f"；回滚失败：{rollback_errors}" if rollback_errors else ""
        raise OSError(f"无法隔离待删除辅助序列：{exc}{detail}") from exc

    try:
        updated = update_design_manifest(
            manifest_path,
            target_compound_id=target_compound,
            sections=sections,
            discard_sections=discard_sections,
            expected_revision=current_revision,
        )
    except Exception as exc:
        rollback_errors = _restore_quarantined_files(moved, quarantine_root)
        if rollback_errors:
            raise RuntimeError(
                "manifest 更新失败，且辅助序列快照回滚不完整："
                + "；".join(rollback_errors)
            ) from exc
        raise

    cleanup_warnings: list[str] = []
    if quarantine_root.exists():
        try:
            shutil.rmtree(quarantine_root)
        except OSError as exc:
            cleanup_warnings.append(
                f"隔离文件未能彻底清理：{quarantine_root}: {exc}"
            )
    _cleanup_empty_snapshot_directories(uploaded_root, snapshots)
    deleted = [
        {
            "accession": accession,
            "sequence_type": str(
                by_accession[accession].get("sequence_type") or ""
            ),
        }
        for accession in requested
    ]
    return {
        "运行成功": not cleanup_warnings,
        "目标化合物": target_compound,
        "已删除辅助蛋白": deleted,
        "已删除项目快照": [
            _relative(project_root, path) for path in snapshots
        ],
        "保留 inputs 原文件": True,
        "剩余辅助蛋白": [
            str(item.get("accession") or "") for item in remaining
        ],
        "清单文件": str(manifest_path),
        "清单版本": updated["revision"],
        "清理警告": cleanup_warnings,
    }


def run_add_manual_auxiliary_protein(config: Any) -> dict[str, Any]:
    result = add_manual_auxiliary_protein(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def run_remove_manual_auxiliary_proteins(config: Any) -> dict[str, Any]:
    result = remove_manual_auxiliary_proteins(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "MANUAL_AUXILIARY_PROTEIN_SCHEMA_VERSION",
    "add_manual_auxiliary_protein",
    "remove_manual_auxiliary_proteins",
    "run_add_manual_auxiliary_protein",
    "run_remove_manual_auxiliary_proteins",
]
