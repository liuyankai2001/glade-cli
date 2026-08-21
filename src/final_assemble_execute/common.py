"""Strict sequence, plan, annotation, and output helpers for final assembly."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import warnings
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

from Bio import Restriction, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from src.final_assemble_plan.common import (
    circular_segment,
    enzyme_cut_positions,
    stable_json_hash,
)
from src.final_assemble_plan.config import (
    FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION,
    SUPPORTED_METHODS,
)
from src.final_assemble_plan.get_final_assembly_context import (
    load_final_assembly_context,
)
from src.final_assemble_execute.models import (
    FinalAssemblyExecutionContext,
    SequenceAssemblyResult,
)
from src.write_manifest.store import read_design_manifest


DNA_ALPHABET = frozenset("ACGT")


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
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
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def relative_project_path(project_root: Path, path: Path) -> str:
    return path.resolve().relative_to(project_root.resolve()).as_posix()


def resolve_project_file(project_root: Path, value: Any) -> Path:
    text = str(value or "").strip()
    if not text:
        raise ValueError("sequence file path is empty")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = project_root / path
    path = path.resolve()
    try:
        path.relative_to(project_root.resolve())
    except ValueError as exc:
        raise ValueError(f"file is outside project output: {path}") from exc
    return path


def read_genbank_strict(path: Path) -> tuple[SeqRecord, str, list[str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            record = SeqIO.read(path, "genbank")
    except Exception as exc:
        raise ValueError(f"could not parse GenBank file: {path}") from exc
    sequence = str(record.seq).upper()
    if not sequence or set(sequence) - DNA_ALPHABET:
        raise ValueError(f"GenBank sequence must contain only A/C/G/T: {path}")
    return record, sequence, [str(item.message) for item in caught]


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少 {name}")
    return value


def load_execution_context(config: Any) -> FinalAssemblyExecutionContext:
    """Load and authenticate the selected 12-plan execution bundle."""

    source_context = load_final_assembly_context(config)
    manifest = read_design_manifest(source_context.manifest_path)
    revision = int(manifest.get("revision") or 0)
    if revision != source_context.manifest_revision:
        raise ValueError("manifest revision 在读取执行上下文时发生变化，请重试")

    section = _mapping(manifest.get("final_assembly_plan"), "final_assembly_plan")
    if section.get("schema_version") != FINAL_ASSEMBLY_PLAN_SCHEMA_VERSION:
        raise ValueError("最终组装计划版本无效，请重新运行并写入 assembly --plan")
    if section.get("status") != "selected":
        raise ValueError("请先使用 write --assembly-plan 接受完整组装计划")
    recorded_selection = str(section.get("selection_fingerprint") or "")
    unsigned_selection = dict(section)
    unsigned_selection.pop("selection_fingerprint", None)
    if recorded_selection != stable_json_hash(unsigned_selection):
        raise ValueError("final_assembly_plan.selection_fingerprint 无效")

    source = _mapping(section.get("source"), "final_assembly_plan.source")
    expected_source = {
        "context_input_fingerprint": source_context.input_fingerprint,
        "parts_selection_fingerprint": source_context.parts_selection_fingerprint,
        "assembled_constructs_fingerprint": (
            source_context.assembled_constructs_fingerprint
        ),
        "plasmid_selection_fingerprint": (
            source_context.plasmid_selection_fingerprint
        ),
    }
    for key, expected in expected_source.items():
        if source.get(key) != expected:
            raise ValueError(
                "最终组装计划与当前表达构建或质粒骨架不一致，请重新生成计划"
            )

    raw_plans = section.get("design_plans")
    if not isinstance(raw_plans, list) or not raw_plans:
        raise ValueError("final_assembly_plan.design_plans 为空")
    plans = [dict(item) for item in raw_plans if isinstance(item, Mapping)]
    if len(plans) != len(raw_plans):
        raise ValueError("final_assembly_plan.design_plans 包含无效条目")
    expected_ids = {item.design_id for item in source_context.constructs}
    actual_ids = [int(item.get("parts_design_id") or 0) for item in plans]
    if (
        len(plans) != len(source_context.constructs)
        or len(actual_ids) != len(set(actual_ids))
        or set(actual_ids) != expected_ids
        or int(section.get("design_count") or 0) != len(plans)
    ):
        raise ValueError("最终组装计划没有完整且唯一地覆盖全部表达构建")
    if source.get("plan_set_fingerprint") != stable_json_hash(plans):
        raise ValueError("final_assembly_plan 的 plan_set_fingerprint 无效")

    construct_by_id = {item.design_id: item for item in source_context.constructs}
    for plan in plans:
        design_id = int(plan["parts_design_id"])
        construct = construct_by_id[design_id]
        if plan.get("assembly_method") not in SUPPORTED_METHODS:
            raise ValueError(f"design {design_id} 包含不支持的组装方法")
        recorded_plan = str(plan.get("plan_fingerprint") or "")
        unsigned_plan = dict(plan)
        unsigned_plan.pop("plan_fingerprint", None)
        if recorded_plan != stable_json_hash(unsigned_plan):
            raise ValueError(f"design {design_id} 的 plan_fingerprint 无效")
        insert = _mapping(plan.get("insert"), f"design {design_id}.insert")
        backbone = _mapping(plan.get("backbone"), f"design {design_id}.backbone")
        if (
            insert.get("file_sha256") != construct.file_sha256
            or insert.get("sequence_sha256") != construct.sequence_sha256
            or int(insert.get("length_bp") or 0) != construct.length_bp
            or backbone.get("file_sha256") != source_context.backbone.file_sha256
            or backbone.get("sequence_content_sha256")
            != source_context.backbone.sequence_content_sha256
            or int(backbone.get("length_bp") or 0)
            != source_context.backbone.length_bp
        ):
            raise ValueError(f"design {design_id} 的输入文件指纹与当前文件不一致")

    backbone_record, backbone_sequence, _ = read_genbank_strict(
        source_context.backbone.path
    )
    if backbone_sequence != source_context.backbone.sequence:
        raise ValueError("质粒骨架序列在执行上下文读取期间发生变化")
    insert_records: dict[int, SeqRecord] = {}
    for construct in source_context.constructs:
        record, sequence, _ = read_genbank_strict(construct.path)
        if sequence != construct.sequence:
            raise ValueError(
                f"design {construct.design_id} 序列在执行上下文读取期间发生变化"
            )
        insert_records[construct.design_id] = record

    return FinalAssemblyExecutionContext(
        manifest=manifest,
        manifest_path=source_context.manifest_path,
        project_output_path=source_context.project_output_path,
        manifest_revision=revision,
        source_context=source_context,
        plan_selection_fingerprint=recorded_selection,
        plans=tuple(sorted(plans, key=lambda item: int(item["parts_design_id"]))),
        backbone_record=backbone_record,
        insert_records=insert_records,
    )


def _target_indexes(target: Mapping[str, Any], vector_length: int) -> tuple[int, int]:
    mode = str(target.get("mode") or "")
    if mode == "insert_after":
        cut = int(target.get("insert_after_bp") or 0)
        if not 1 <= cut <= vector_length:
            raise ValueError("insert_after_bp 超出骨架范围")
        return cut, cut
    if mode == "replace":
        start = int(target.get("replace_start_bp") or 0)
        end = int(target.get("replace_end_bp") or 0)
        if start < 1 or end < start or end > vector_length:
            raise ValueError("replace 坐标超出骨架范围")
        return start - 1, end
    raise ValueError(f"不支持的 target.mode：{mode}")


def _enzyme_from_record(record: Mapping[str, Any]) -> Any:
    name = str(record.get("name") or "")
    enzyme = getattr(Restriction, name, None)
    if enzyme is None:
        raise ValueError(f"Biopython 不支持计划中的限制酶：{name}")
    if str(enzyme.site) != str(record.get("recognition_site") or "").upper():
        raise ValueError(f"限制酶 {name} 的识别序列与 Biopython 不一致")
    return enzyme


def _validate_linearization_enzyme(
    record: Mapping[str, Any],
    vector_sequence: str,
) -> Any:
    enzyme = _enzyme_from_record(record)
    cuts = enzyme_cut_positions(enzyme, vector_sequence, circular=True)
    if cuts != [int(record.get("cut_after_bp"))]:
        raise ValueError(
            f"限制酶 {record.get('name')} 的骨架切割坐标与计划不一致"
        )
    start = int(record.get("site_start_bp") or 0)
    end = int(record.get("site_end_bp") or 0)
    site = str(record.get("recognition_site") or "").upper()
    if start < 1 or end != start + len(site) - 1:
        raise ValueError(f"限制酶 {record.get('name')} 的识别位点坐标无效")
    if vector_sequence[start - 1 : end] != site:
        raise ValueError(f"限制酶 {record.get('name')} 的识别位点与骨架不一致")
    return enzyme


def _assemble_restriction(
    vector_sequence: str,
    insert_sequence: str,
    plan: Mapping[str, Any],
) -> SequenceAssemblyResult:
    target = _mapping(plan.get("target"), "restriction target")
    start_index, end_index = _target_indexes(target, len(vector_sequence))
    if target.get("mode") != "replace":
        raise ValueError("restriction 计划必须使用 replace target")
    restriction = _mapping(plan.get("restriction"), "restriction parameters")
    if restriction.get("restriction_site_retention") != "retain":
        raise ValueError("当前执行器仅接受 retain restriction sites 的计划")
    left_site = str(restriction.get("left_site") or "").upper()
    right_site = str(restriction.get("right_site") or "").upper()
    if not left_site or not right_site:
        raise ValueError("restriction 计划缺少左右识别序列")
    if vector_sequence[start_index : start_index + len(left_site)] != left_site:
        raise ValueError("左限制酶位点与 target.replace_start_bp 不一致")
    if vector_sequence[end_index - len(right_site) : end_index] != right_site:
        raise ValueError("右限制酶位点与 target.replace_end_bp 不一致")

    linearization = _mapping(
        plan.get("backbone_linearization"), "backbone_linearization"
    )
    enzymes = linearization.get("restriction_enzymes")
    if linearization.get("mode") != "restriction" or not isinstance(enzymes, list):
        raise ValueError("restriction 计划缺少骨架酶切线性化记录")
    by_role = {
        str(item.get("role")): item
        for item in enzymes
        if isinstance(item, Mapping)
    }
    left_record = _mapping(by_role.get("left"), "left restriction enzyme")
    right_record = _mapping(by_role.get("right"), "right restriction enzyme")
    if (
        left_record.get("name") != restriction.get("left_enzyme")
        or right_record.get("name") != restriction.get("right_enzyme")
    ):
        raise ValueError("左右限制酶名称在计划字段间不一致")
    for record in (left_record, right_record):
        enzyme = _validate_linearization_enzyme(record, vector_sequence)
        if enzyme_cut_positions(enzyme, insert_sequence, circular=False):
            raise ValueError(
                f"限制酶 {record.get('name')} 在完整表达构建内部存在切点"
            )

    payload = left_site + insert_sequence + right_site
    final_sequence = (
        vector_sequence[:start_index] + payload + vector_sequence[end_index:]
    )
    inserted_start = start_index + len(left_site) + 1
    inserted_end = inserted_start + len(insert_sequence) - 1
    return SequenceAssemblyResult(
        sequence=final_sequence,
        inserted_start_bp=inserted_start,
        inserted_end_bp=inserted_end,
        vector_replace_start_index=start_index,
        vector_replace_end_index=end_index,
        replacement_payload_length_bp=len(payload),
        target_audit={
            "mode": "restriction_replace_retain_sites",
            "replace_start_bp": start_index + 1,
            "replace_end_bp": end_index,
            "inserted_start_bp": inserted_start,
            "inserted_end_bp": inserted_end,
            "left_enzyme": restriction.get("left_enzyme"),
            "right_enzyme": restriction.get("right_enzyme"),
            "left_site": left_site,
            "right_site": right_site,
            "retained_sites": {
                "left": {
                    "start_bp": start_index + 1,
                    "end_bp": start_index + len(left_site),
                },
                "right": {
                    "start_bp": inserted_end + 1,
                    "end_bp": inserted_end + len(right_site),
                },
            },
        },
    )


def _validate_gibson_linearization(
    plan: Mapping[str, Any],
    vector_sequence: str,
) -> None:
    linearization = _mapping(
        plan.get("backbone_linearization"), "backbone_linearization"
    )
    enzymes = linearization.get("restriction_enzymes")
    if not isinstance(enzymes, list):
        raise ValueError("Gibson 计划缺少 restriction_enzymes 列表")
    mode = str(linearization.get("mode") or "")
    if mode == "pcr":
        if enzymes or linearization.get("enzyme_summary") != "none (PCR linearization)":
            raise ValueError("PCR 线性化记录必须明确为无限制酶")
        return
    if mode != "restriction" or not enzymes:
        raise ValueError("Gibson 骨架线性化方式无效")
    for item in enzymes:
        if not isinstance(item, Mapping):
            raise ValueError("Gibson 限制酶线性化记录无效")
        _validate_linearization_enzyme(item, vector_sequence)


def _assemble_gibson(
    vector_sequence: str,
    insert_sequence: str,
    plan: Mapping[str, Any],
) -> SequenceAssemblyResult:
    target = _mapping(plan.get("target"), "Gibson target")
    start_index, end_index = _target_indexes(target, len(vector_sequence))
    gibson = _mapping(plan.get("gibson"), "Gibson parameters")
    arm_length = int(gibson.get("homology_arm_length") or 0)
    left_arm = str(gibson.get("left_homology") or "").upper()
    right_arm = str(gibson.get("right_homology") or "").upper()
    if arm_length < 1 or len(left_arm) != arm_length or len(right_arm) != arm_length:
        raise ValueError("Gibson 同源臂长度无效")
    if target.get("mode") == "insert_after":
        cut = int(target["insert_after_bp"])
        left_start = cut - arm_length + 1
        right_start = cut + 1
    else:
        left_start = int(target["replace_start_bp"]) - arm_length
        right_start = int(target["replace_end_bp"]) + 1
    expected_left = circular_segment(vector_sequence, left_start, arm_length)
    expected_right = circular_segment(vector_sequence, right_start, arm_length)
    if left_arm != expected_left or right_arm != expected_right:
        raise ValueError("Gibson 同源臂与计划插入位置附近的骨架序列不一致")
    _validate_gibson_linearization(plan, vector_sequence)

    final_sequence = (
        vector_sequence[:start_index] + insert_sequence + vector_sequence[end_index:]
    )
    inserted_start = start_index + 1
    inserted_end = inserted_start + len(insert_sequence) - 1
    return SequenceAssemblyResult(
        sequence=final_sequence,
        inserted_start_bp=inserted_start,
        inserted_end_bp=inserted_end,
        vector_replace_start_index=start_index,
        vector_replace_end_index=end_index,
        replacement_payload_length_bp=len(insert_sequence),
        target_audit={
            "mode": str(target.get("mode")),
            "insert_after_bp": target.get("insert_after_bp"),
            "replace_start_bp": target.get("replace_start_bp"),
            "replace_end_bp": target.get("replace_end_bp"),
            "inserted_start_bp": inserted_start,
            "inserted_end_bp": inserted_end,
            "homology_arm_length": arm_length,
            "left_homology": left_arm,
            "right_homology": right_arm,
            "linearization": deepcopy(plan.get("backbone_linearization")),
        },
    )


def assemble_sequence(
    vector_sequence: str,
    insert_sequence: str,
    plan: Mapping[str, Any],
) -> SequenceAssemblyResult:
    method = str(plan.get("assembly_method") or "")
    if method == "restriction":
        result = _assemble_restriction(vector_sequence, insert_sequence, plan)
    elif method == "gibson":
        result = _assemble_gibson(vector_sequence, insert_sequence, plan)
    else:
        raise ValueError(f"不支持的组装方法：{method}")
    expected_length = int(plan.get("estimated_final_length_bp") or 0)
    if expected_length != len(result.sequence):
        raise ValueError(
            f"预计最终长度 {expected_length} bp 与执行结果 {len(result.sequence)} bp 不一致"
        )
    return result


def _mapped_simple_parts(
    location: Any,
    *,
    replace_start: int,
    replace_end: int,
    payload_length: int,
) -> tuple[list[SimpleLocation], bool]:
    start = int(location.start)
    end = int(location.end)
    strand = location.strand
    delta = payload_length - (replace_end - replace_start)
    changed = False
    if replace_start == replace_end:
        cut = replace_start
        if end <= cut:
            return [SimpleLocation(start, end, strand=strand)], False
        if start >= cut:
            return [SimpleLocation(start + payload_length, end + payload_length, strand=strand)], True
        changed = True
        parts = []
        if start < cut:
            parts.append(SimpleLocation(start, cut, strand=strand))
        if cut < end:
            parts.append(
                SimpleLocation(cut + payload_length, end + payload_length, strand=strand)
            )
        return parts, changed

    if end <= replace_start:
        return [SimpleLocation(start, end, strand=strand)], False
    if start >= replace_end:
        return [SimpleLocation(start + delta, end + delta, strand=strand)], True
    changed = True
    parts: list[SimpleLocation] = []
    if start < replace_start:
        parts.append(SimpleLocation(start, replace_start, strand=strand))
    if end > replace_end:
        parts.append(
            SimpleLocation(
                replace_start + payload_length,
                end + delta,
                strand=strand,
            )
        )
    return parts, changed


def _map_vector_feature(
    feature: SeqFeature,
    *,
    replace_start: int,
    replace_end: int,
    payload_length: int,
) -> tuple[SeqFeature | None, str | None]:
    if feature.type == "source" or feature.location is None:
        return None, None
    original_parts = (
        list(feature.location.parts)
        if isinstance(feature.location, CompoundLocation)
        else [feature.location]
    )
    mapped_parts: list[SimpleLocation] = []
    changed = False
    for part in original_parts:
        pieces, part_changed = _mapped_simple_parts(
            part,
            replace_start=replace_start,
            replace_end=replace_end,
            payload_length=payload_length,
        )
        mapped_parts.extend(pieces)
        changed = changed or part_changed
    label = _feature_label(feature)
    if not mapped_parts:
        return None, f"dropped vector feature inside replaced region: {label}"
    location: Any
    if len(mapped_parts) == 1:
        location = mapped_parts[0]
    else:
        operator = (
            feature.location.operator
            if isinstance(feature.location, CompoundLocation)
            else "join"
        )
        location = CompoundLocation(mapped_parts, operator=operator)
    mapped = SeqFeature(
        location=location,
        type=feature.type,
        id=feature.id,
        qualifiers=deepcopy(feature.qualifiers),
    )
    warning = f"remapped vector feature across assembly target: {label}" if changed else None
    return mapped, warning


def _feature_label(feature: SeqFeature) -> str:
    for key in ("label", "gene", "product", "note"):
        values = feature.qualifiers.get(key, [])
        if values:
            return str(values[0])
    return feature.type


def _shift_feature(feature: SeqFeature, offset: int) -> SeqFeature | None:
    if feature.type == "source" or feature.location is None:
        return None
    parts = (
        list(feature.location.parts)
        if isinstance(feature.location, CompoundLocation)
        else [feature.location]
    )
    shifted = [
        SimpleLocation(
            int(part.start) + offset,
            int(part.end) + offset,
            strand=part.strand,
        )
        for part in parts
    ]
    if len(shifted) == 1:
        location: Any = shifted[0]
    else:
        location = CompoundLocation(
            shifted,
            operator=(
                feature.location.operator
                if isinstance(feature.location, CompoundLocation)
                else "join"
            ),
        )
    qualifiers = deepcopy(feature.qualifiers)
    qualifiers.setdefault("assembly_source", ["complete_expression_construct"])
    return SeqFeature(
        location=location,
        type=feature.type,
        id=feature.id,
        qualifiers=qualifiers,
    )


def build_final_record(
    *,
    record_id: str,
    backbone_record: SeqRecord,
    insert_record: SeqRecord,
    result: SequenceAssemblyResult,
    plan: Mapping[str, Any],
) -> tuple[SeqRecord, list[str]]:
    warnings_list: list[str] = []
    features: list[SeqFeature] = [
        SeqFeature(
            SimpleLocation(0, len(result.sequence), strand=1),
            type="source",
            qualifiers={
                "organism": ["synthetic DNA construct"],
                "mol_type": ["other DNA"],
            },
        )
    ]
    for feature in backbone_record.features:
        mapped, warning = _map_vector_feature(
            feature,
            replace_start=result.vector_replace_start_index,
            replace_end=result.vector_replace_end_index,
            payload_length=result.replacement_payload_length_bp,
        )
        if mapped is not None:
            mapped.qualifiers.setdefault("assembly_source", ["plasmid_backbone"])
            features.append(mapped)
        if warning:
            warnings_list.append(warning)

    insert_offset = result.inserted_start_bp - 1
    features.append(
        SeqFeature(
            SimpleLocation(
                insert_offset,
                result.inserted_end_bp,
                strand=1,
            ),
            type="misc_feature",
            qualifiers={
                "label": [f"expression_design_{int(plan['parts_design_id']):03d}"],
                "note": ["Complete expression construct inserted in forward orientation"],
                "assembly_source": ["complete_expression_construct"],
            },
        )
    )
    for feature in insert_record.features:
        shifted = _shift_feature(feature, insert_offset)
        if shifted is not None:
            features.append(shifted)

    retained = result.target_audit.get("retained_sites")
    if isinstance(retained, Mapping):
        for role in ("left", "right"):
            site = retained.get(role)
            if not isinstance(site, Mapping):
                continue
            enzyme = result.target_audit.get(f"{role}_enzyme")
            recognition = result.target_audit.get(f"{role}_site")
            features.append(
                SeqFeature(
                    SimpleLocation(
                        int(site["start_bp"]) - 1,
                        int(site["end_bp"]),
                        strand=1,
                    ),
                    type="misc_feature",
                    qualifiers={
                        "label": [f"{enzyme}_retained_site"],
                        "note": [f"Retained {enzyme} recognition site {recognition}"],
                        "assembly_source": ["final_assembly_plan"],
                    },
                )
            )

    def feature_key(item: SeqFeature) -> tuple[int, int, str]:
        return (int(item.location.start), int(item.location.end), item.type)

    source_feature = features[0]
    remaining = sorted(features[1:], key=feature_key)
    final = SeqRecord(
        Seq(result.sequence),
        id=record_id,
        name=record_id[:16],
        description=(
            f"In-silico assembled final plasmid for expression parts design "
            f"{int(plan['parts_design_id'])}"
        ),
    )
    final.annotations["molecule_type"] = "DNA"
    final.annotations["topology"] = "circular"
    final.annotations["data_file_division"] = "SYN"
    final.features = [source_feature, *remaining]
    return final, warnings_list


def write_genbank(path: Path, record: SeqRecord) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(record, path, "genbank")


def write_fasta(path: Path, record_id: str, sequence: str) -> None:
    record = SeqRecord(Seq(sequence), id=record_id, description="in-silico final plasmid")
    SeqIO.write(record, path, "fasta")


def validate_written_outputs(
    *,
    genbank_path: Path,
    fasta_path: Path,
    expected_sequence: str,
    inserted_start_bp: int,
    inserted_end_bp: int,
    expected_insert: str,
) -> dict[str, Any]:
    genbank_record = SeqIO.read(genbank_path, "genbank")
    fasta_record = SeqIO.read(fasta_path, "fasta")
    genbank_sequence = str(genbank_record.seq).upper()
    fasta_sequence = str(fasta_record.seq).upper()
    if genbank_sequence != expected_sequence:
        raise ValueError("写出的 GenBank 序列与理论组装结果不一致")
    if fasta_sequence != expected_sequence:
        raise ValueError("写出的 FASTA 序列与理论组装结果不一致")
    if str(genbank_record.annotations.get("topology") or "").lower() != "circular":
        raise ValueError("写出的 GenBank 未保留 circular topology")
    if genbank_sequence[inserted_start_bp - 1 : inserted_end_bp] != expected_insert:
        raise ValueError("最终 GenBank 的插入片段与完整表达构建不一致")
    invalid_features = [
        feature.type
        for feature in genbank_record.features
        if int(feature.location.start) < 0
        or int(feature.location.end) > len(expected_sequence)
        or int(feature.location.end) < int(feature.location.start)
    ]
    if invalid_features:
        raise ValueError(f"最终 GenBank 包含越界 feature：{invalid_features}")
    return {
        "status": "PASS",
        "genbank_parse_ok": True,
        "fasta_parse_ok": True,
        "genbank_sequence_matches": True,
        "fasta_sequence_matches": True,
        "insert_sequence_matches": True,
        "topology_circular": True,
        "feature_coordinates_valid": True,
        "sequence_sha256": sha256_sequence(expected_sequence),
    }


def file_summary(project_root: Path, path: Path, *, file_format: str) -> dict[str, Any]:
    return {
        "path": relative_project_path(project_root, path),
        "format": file_format,
        "file_sha256": sha256_file(path),
    }


def inspect_file(project_root: Path, info: Mapping[str, Any]) -> dict[str, Any]:
    try:
        path = resolve_project_file(project_root, info.get("path"))
    except Exception as exc:
        return {"path": str(info.get("path") or ""), "exists": False, "error": str(exc)}
    exists = path.is_file()
    actual = sha256_file(path) if exists else ""
    expected = str(info.get("file_sha256") or "")
    return {
        "path": relative_project_path(project_root, path),
        "exists": exists,
        "file_sha256_matches": bool(exists and expected and actual == expected),
        "actual_file_sha256": actual,
    }


__all__ = [
    "assemble_sequence",
    "build_final_record",
    "file_summary",
    "inspect_file",
    "load_execution_context",
    "read_genbank_strict",
    "relative_project_path",
    "resolve_project_file",
    "sha256_file",
    "sha256_sequence",
    "stable_json_hash",
    "validate_written_outputs",
    "write_fasta",
    "write_genbank",
    "write_json_atomic",
]
