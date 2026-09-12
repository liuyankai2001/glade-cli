"""Read-only expression-box view with optional expression-part overlays."""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from src.write_manifest.store import read_design_manifest


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} 必须是 JSON 对象")
    return value


def _text(value: Any, field: str) -> str:
    normalized = str(value or "").strip()
    if not normalized:
        raise ValueError(f"{field} 不能为空")
    return normalized


def _positive_int(value: Any, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} 必须是正整数")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} 必须是正整数") from exc
    if parsed < 1:
        raise ValueError(f"{field} 必须是正整数")
    return parsed


def _stable_json_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _project_file(project_root: Path, path_value: Any, field: str) -> Path:
    value = _text(path_value, field)
    path = Path(value)
    if not path.is_absolute():
        path = project_root / path
    resolved = path.expanduser().resolve()
    if resolved == project_root or project_root not in resolved.parents:
        raise ValueError(f"{field} 超出当前项目输出目录")
    return resolved


def _cds_by_accession(manifest: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    selection = _mapping(manifest.get("cds_selection"), "cds_selection")
    proteins = selection.get("proteins")
    if not isinstance(proteins, list) or not proteins:
        raise ValueError("cds_selection.proteins 必须是非空列表")
    result: dict[str, dict[str, Any]] = {}
    for index, raw in enumerate(proteins):
        item = _mapping(raw, f"cds_selection.proteins[{index}]")
        accession = _text(
            item.get("accession"),
            f"cds_selection.proteins[{index}].accession",
        ).upper()
        if accession in result:
            raise ValueError(f"cds_selection 包含重复蛋白：{accession}")
        optimized = _mapping(
            item.get("optimized_cds"),
            f"cds_selection.proteins[{index}].optimized_cds",
        )
        length_nt = _positive_int(
            optimized.get("length_nt"),
            f"cds_selection.proteins[{index}].optimized_cds.length_nt",
        )
        sequence_sha256 = _text(
            optimized.get("sequence_sha256"),
            f"cds_selection.proteins[{index}].optimized_cds.sequence_sha256",
        ).lower()
        if _SHA256_RE.fullmatch(sequence_sha256) is None:
            raise ValueError(f"{accession} 的 CDS SHA-256 无效")
        roles = item.get("roles")
        if not isinstance(roles, list) or not roles:
            raise ValueError(f"{accession} 的 roles 必须是非空列表")
        result[accession] = {
            "长度_nt": length_nt,
            "序列SHA256": sequence_sha256,
            "优化已跳过": optimized.get("optimization_skipped") is True,
            "路径": str(optimized.get("path") or ""),
            "角色": [str(role) for role in roles],
        }
    return result


def _base_cassettes(
    selection: Mapping[str, Any],
    cds_by_accession: Mapping[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    raw_cassettes = selection.get("cassettes")
    if not isinstance(raw_cassettes, list) or not raw_cassettes:
        raise ValueError("expression_box_selection.cassettes 必须是非空列表")
    assigned: list[str] = []
    result: list[dict[str, Any]] = []
    for expected_index, raw in enumerate(raw_cassettes, start=1):
        cassette = _mapping(raw, f"expression_box_selection.cassettes[{expected_index - 1}]")
        cassette_index = _positive_int(
            cassette.get("cassette_index"),
            f"expression_box_selection.cassettes[{expected_index - 1}].cassette_index",
        )
        if cassette_index != expected_index:
            raise ValueError("expression_box_selection 的表达盒编号必须连续")
        raw_accessions = cassette.get("protein_accessions")
        if not isinstance(raw_accessions, list) or not raw_accessions:
            raise ValueError(f"表达盒 {cassette_index} 没有蛋白")
        accessions = [str(value or "").strip().upper() for value in raw_accessions]
        if any(not value for value in accessions):
            raise ValueError(f"表达盒 {cassette_index} 包含空蛋白 ID")
        unknown = sorted(set(accessions) - set(cds_by_accession))
        if unknown:
            raise ValueError("表达盒引用未知 CDS：" + ", ".join(unknown))
        assigned.extend(accessions)
        result.append(
            {
                "表达盒编号": cassette_index,
                "状态": "表达元件未设置",
                "Promoter": None,
                "基因列表": [
                    {
                        "顺序": order,
                        "蛋白ID": accession,
                        "角色": list(cds_by_accession[accession]["角色"]),
                        "CDS": {
                            key: value
                            for key, value in cds_by_accession[accession].items()
                            if key != "角色"
                        },
                        "RBS": None,
                        "OSTIR": None,
                    }
                    for order, accession in enumerate(accessions, start=1)
                ],
                "Terminator": None,
                "序列检查": None,
            }
        )
    if len(assigned) != len(set(assigned)):
        raise ValueError("expression_box_selection 重复分配了蛋白")
    if set(assigned) != set(cds_by_accession):
        missing = sorted(set(cds_by_accession) - set(assigned))
        raise ValueError("expression_box_selection 遗漏蛋白：" + ", ".join(missing))
    return result


def _part_summary(part: Any, expected_role: str) -> dict[str, Any]:
    item = _mapping(part, expected_role)
    role = _text(item.get("role"), f"{expected_role}.role").lower()
    if role != expected_role:
        raise ValueError(f"表达元件角色不匹配：预期 {expected_role}，实际 {role}")
    part_id = _text(item.get("part_id"), f"{expected_role}.part_id")
    length_bp = _positive_int(item.get("length_bp"), f"{part_id}.length_bp")
    sequence_sha256 = _text(
        item.get("sequence_sha256"),
        f"{part_id}.sequence_sha256",
    ).lower()
    if _SHA256_RE.fullmatch(sequence_sha256) is None:
        raise ValueError(f"{part_id} 的序列 SHA-256 无效")
    sequence_file = item.get("sequence_file")
    if isinstance(sequence_file, Mapping):
        sequence_file = sequence_file.get("path")
    return {
        "元件ID": part_id,
        "角色": role,
        "来源": str(item.get("source") or item.get("source_type") or "未记录"),
        "长度_bp": length_bp,
        "序列SHA256": sequence_sha256,
        "序列文件": str(sequence_file or ""),
        "强度": str(item.get("strength") or ""),
        "调控类型": str(item.get("regulation") or ""),
        "宿主匹配": str(item.get("host_match_kind") or ""),
    }


def _ostir_summary(value: Any) -> dict[str, Any]:
    item = _mapping(value, "ostir")
    try:
        expression = float(item.get("translation_initiation_rate"))
        energy = float(item.get("d_g_total"))
        intended = int(item.get("intended_start_position"))
        unintended = int(item.get("unintended_start_count"))
    except (TypeError, ValueError) as exc:
        raise ValueError("OSTIR 结果包含无效数值") from exc
    return {
        "翻译起始率": expression,
        "总自由能": energy,
        "预期起始位置": intended,
        "非预期起始位点数": unintended,
        "上下文SHA256": str(item.get("context_sha256") or ""),
    }


def _audit_summary(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    item = _mapping(value, "sequence_audit")
    return {
        "状态": str(item.get("gate_status") or "未知"),
        "长度_nt": item.get("length_nt"),
        "GC百分比": item.get("gc_percent"),
        "局部GC最小百分比": item.get("local_gc_min_percent"),
        "局部GC最大百分比": item.get("local_gc_max_percent"),
        "禁止位点": item.get("forbidden_site_hits", {}),
        "最长同聚物": item.get("max_homopolymer"),
        "失败检查": list(item.get("failed_checks") or []),
    }


def _read_artifact(project_root: Path, source: Mapping[str, Any]) -> dict[str, Any]:
    path = _project_file(project_root, source.get("artifact"), "parts_selection.source.artifact")
    if not path.is_file():
        raise FileNotFoundError(f"表达元件方案文件不存在：{path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"表达元件方案文件不是有效 JSON：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("表达元件方案文件根节点必须是 JSON 对象")
    payload["_resolved_path"] = str(path)
    return payload


def _design_references(selection: Mapping[str, Any]) -> dict[int, Mapping[str, Any]]:
    raw = selection.get("design_references")
    if not isinstance(raw, list) or not raw:
        raise ValueError("parts_selection.design_references 必须是非空列表")
    result: dict[int, Mapping[str, Any]] = {}
    for item_value in raw:
        item = _mapping(item_value, "parts_selection.design_references[]")
        design_id = _positive_int(item.get("design_id"), "design_reference.design_id")
        if design_id in result:
            raise ValueError(f"parts_selection 重复引用方案 {design_id}")
        result[design_id] = item
    return result


def _selected_design_ids(selection: Mapping[str, Any]) -> list[int]:
    raw = selection.get("selected_design_ids")
    if not isinstance(raw, list) or not raw:
        raise ValueError("parts_selection.selected_design_ids 必须是非空列表")
    values = [_positive_int(value, "selected_design_ids") for value in raw]
    if len(values) != len(set(values)):
        raise ValueError("parts_selection.selected_design_ids 包含重复编号")
    return values


def _overlay_design(
    cassettes: list[dict[str, Any]],
    design: Mapping[str, Any],
) -> None:
    raw_cassettes = design.get("cassettes")
    if not isinstance(raw_cassettes, list):
        raise ValueError("表达元件方案缺少 cassettes")
    by_index = {
        _positive_int(item.get("cassette_index"), "cassette_index"): item
        for item in raw_cassettes
        if isinstance(item, Mapping)
    }
    expected_indexes = {item["表达盒编号"] for item in cassettes}
    if set(by_index) != expected_indexes or len(raw_cassettes) != len(expected_indexes):
        raise ValueError("表达元件方案与当前表达盒数量不一致")

    for cassette in cassettes:
        cassette_index = cassette["表达盒编号"]
        raw = by_index[cassette_index]
        genes = raw.get("genes")
        if not isinstance(genes, list) or len(genes) != len(cassette["基因列表"]):
            raise ValueError(f"表达盒 {cassette_index} 的 RBS 数量不匹配")
        cassette["Promoter"] = _part_summary(raw.get("promoter"), "promoter")
        cassette["Terminator"] = _part_summary(raw.get("terminator"), "terminator")
        cassette["序列检查"] = _audit_summary(raw.get("sequence_audit"))
        for raw_gene_value, gene in zip(genes, cassette["基因列表"], strict=True):
            raw_gene = _mapping(raw_gene_value, "gene")
            if str(raw_gene.get("accession") or "").strip().upper() != gene["蛋白ID"]:
                raise ValueError(f"表达盒 {cassette_index} 的蛋白顺序不匹配")
            if raw_gene.get("cds_sequence_sha256") != gene["CDS"]["序列SHA256"]:
                raise ValueError(f"{gene['蛋白ID']} 的 CDS 已发生变化")
            gene["RBS"] = _part_summary(raw_gene.get("rbs"), "rbs")
            gene["OSTIR"] = _ostir_summary(raw_gene.get("ostir"))
        cassette["状态"] = "表达元件已设置"


def _construct_summary(
    manifest: Mapping[str, Any],
    project_root: Path,
    selection: Mapping[str, Any],
    design_id: int,
    warnings: list[str],
) -> dict[str, Any]:
    empty = {
        "状态": "未生成",
        "路径": "",
        "文件存在": False,
        "文件哈希匹配": False,
        "长度_bp": None,
        "序列SHA256": "",
        "文件SHA256": "",
        "序列检查": None,
    }
    assembled = manifest.get("assembled_expression_constructs")
    if not isinstance(assembled, Mapping):
        return empty
    if assembled.get("source_parts_selection_fingerprint") != selection.get(
        "selection_fingerprint"
    ):
        warnings.append("完整表达构建与当前表达元件选择不一致")
        return {**empty, "状态": "信息不可用"}
    constructs = assembled.get("constructs")
    if not isinstance(constructs, list):
        warnings.append("assembled_expression_constructs.constructs 无效")
        return {**empty, "状态": "信息不可用"}
    matches = [
        item
        for item in constructs
        if isinstance(item, Mapping) and item.get("parts_design_id") == design_id
    ]
    if len(matches) != 1:
        warnings.append(f"未找到表达元件方案 {design_id} 对应的完整表达构建")
        return empty
    item = matches[0]
    try:
        path = _project_file(project_root, item.get("path"), "construct.path")
    except ValueError as exc:
        warnings.append(str(exc))
        return {**empty, "状态": "信息不可用"}
    exists = path.is_file()
    expected_hash = str(item.get("file_sha256") or "")
    hash_matches = exists and bool(expected_hash) and _sha256_file(path) == expected_hash
    if not exists:
        warnings.append(f"完整表达构建文件不存在：{path}")
    elif not hash_matches:
        warnings.append(f"完整表达构建文件哈希不匹配：{path}")
    return {
        "状态": "已生成" if hash_matches else "文件校验失败",
        "路径": str(path),
        "文件存在": exists,
        "文件哈希匹配": hash_matches,
        "长度_bp": item.get("length_bp"),
        "序列SHA256": str(item.get("sequence_sha256") or ""),
        "文件SHA256": expected_hash,
        "序列检查": _audit_summary(item.get("sequence_audit")),
    }


def _parts_unavailable(result: dict[str, Any], message: str) -> None:
    result["表达元件状态"] = "信息不可用"
    result["警告"].append(message)


def _draft_part_summary(
    part: Any,
    project_root: Path,
) -> dict[str, Any]:
    item = _mapping(part, "expression_parts_draft.promoter")
    summary = _part_summary(item, "promoter")
    sequence_file = _mapping(item.get("sequence_file"), "promoter.sequence_file")
    path = _project_file(project_root, sequence_file.get("path"), "promoter.sequence_file.path")
    if not path.is_file():
        raise FileNotFoundError(f"启动子快照文件不存在：{path}")
    expected_hash = _text(
        sequence_file.get("file_sha256"),
        "promoter.sequence_file.file_sha256",
    )
    if _sha256_file(path) != expected_hash:
        raise ValueError(f"启动子快照文件哈希不匹配：{path}")
    summary["序列文件"] = str(path)
    return summary


def _apply_draft(
    result: dict[str, Any],
    manifest: Mapping[str, Any],
    project_root: Path,
) -> None:
    raw_draft = manifest.get("expression_parts_draft")
    if raw_draft is None:
        return
    try:
        draft = _mapping(raw_draft, "expression_parts_draft")
        if draft.get("schema_version") != "expression_parts_draft.v1":
            raise ValueError("不支持的 expression_parts_draft schema_version")
        if str(draft.get("target_compound_id") or "") != result["目标化合物"]:
            raise ValueError("表达元件草稿与当前目标化合物不一致")
        source = _mapping(draft.get("source"), "expression_parts_draft.source")
        if source.get("expression_box_selection_fingerprint") != result[
            "表达盒方案"
        ]["方案Fingerprint"]:
            raise ValueError("表达元件草稿与当前表达盒方案不一致")
        content = {
            "target_compound_id": draft.get("target_compound_id"),
            "source": dict(source),
            "cassettes": draft.get("cassettes"),
        }
        if draft.get("draft_fingerprint") != _stable_json_hash(content):
            raise ValueError("表达元件草稿 fingerprint 校验失败")
        raw_cassettes = draft.get("cassettes")
        if not isinstance(raw_cassettes, list):
            raise ValueError("expression_parts_draft.cassettes 必须是列表")
        by_index = {
            _positive_int(item.get("cassette_index"), "draft.cassette_index"): item
            for item in raw_cassettes
            if isinstance(item, Mapping)
        }
        expected_indexes = {
            cassette["表达盒编号"] for cassette in result["表达盒列表"]
        }
        if set(by_index) != expected_indexes or len(raw_cassettes) != len(
            expected_indexes
        ):
            raise ValueError("表达元件草稿与当前表达盒数量不一致")
        promoter_count = 0
        for cassette in result["表达盒列表"]:
            raw = by_index[cassette["表达盒编号"]]
            expected_accessions = [
                gene["蛋白ID"] for gene in cassette["基因列表"]
            ]
            if raw.get("protein_accessions") != expected_accessions:
                raise ValueError(
                    f"表达盒 {cassette['表达盒编号']} 的草稿蛋白顺序不匹配"
                )
            promoter = raw.get("promoter")
            if promoter is None:
                continue
            cassette["Promoter"] = _draft_part_summary(promoter, project_root)
            cassette["状态"] = "表达元件部分设置"
            promoter_count += 1
        result["表达元件草稿"] = {
            "状态": str(draft.get("status") or "partial"),
            "已上传启动子数": promoter_count,
            "表达盒总数": len(result["表达盒列表"]),
            "草稿Fingerprint": str(draft.get("draft_fingerprint") or ""),
        }
        if promoter_count:
            result["表达元件状态"] = "部分设置"
    except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
        _parts_unavailable(result, str(exc))


def _apply_parts(
    result: dict[str, Any],
    manifest: Mapping[str, Any],
    project_root: Path,
    requested_design_id: int | None,
) -> None:
    raw_selection = manifest.get("parts_selection")
    if raw_selection is None:
        if requested_design_id is not None:
            raise ValueError("尚未选择表达元件方案，不能使用 --parts-design")
        _apply_draft(result, manifest, project_root)
        return
    try:
        selection = _mapping(raw_selection, "parts_selection")
        selected_ids = _selected_design_ids(selection)
        primary_id = _positive_int(selection.get("primary_design_id"), "primary_design_id")
        if primary_id not in selected_ids:
            raise ValueError("parts_selection.primary_design_id 不在已选方案中")
        references = _design_references(selection)
        if set(references) != set(selected_ids):
            raise ValueError("parts_selection 的方案引用与已选编号不一致")
    except (TypeError, ValueError) as exc:
        _parts_unavailable(result, str(exc))
        return

    result["已选表达元件方案"] = selected_ids
    result["主表达元件方案"] = primary_id
    display_id = requested_design_id if requested_design_id is not None else primary_id
    if display_id not in selected_ids:
        raise ValueError(
            f"表达元件方案 {display_id} 未被选中；可查看方案：{selected_ids}"
        )
    result["当前显示方案"] = display_id

    source = selection.get("source")
    if not isinstance(source, Mapping):
        _parts_unavailable(result, "parts_selection.source 无效")
        return
    box_fingerprint = result["表达盒方案"]["方案Fingerprint"]
    if source.get("expression_box_selection_fingerprint") != box_fingerprint:
        _parts_unavailable(result, "表达元件选择与当前表达盒方案不一致")
        return

    try:
        artifact = _read_artifact(project_root, source)
        if str(artifact.get("target_compound_id") or "") != result["目标化合物"]:
            raise ValueError("表达元件方案文件与当前目标化合物不一致")
        artifact_source = _mapping(artifact.get("source"), "expression-parts artifact.source")
        if artifact_source.get("expression_box_selection_fingerprint") != box_fingerprint:
            raise ValueError("表达元件方案文件与当前表达盒方案不一致")
        raw_designs = artifact.get("designs")
        if not isinstance(raw_designs, list):
            raise ValueError("表达元件方案文件缺少 designs")
        by_id = {
            _positive_int(item.get("design_id"), "design.design_id"): item
            for item in raw_designs
            if isinstance(item, Mapping)
        }
        if display_id not in by_id:
            raise ValueError(f"表达元件方案文件中不存在方案 {display_id}")
        design = by_id[display_id]
        expected_design_hash = str(references[display_id].get("design_fingerprint") or "")
        if not expected_design_hash or _stable_json_hash(design) != expected_design_hash:
            raise ValueError(f"表达元件方案 {display_id} 的内容 fingerprint 不匹配")
        _overlay_design(result["表达盒列表"], design)
    except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
        _parts_unavailable(result, str(exc))
        return

    reference = references[display_id]
    result["表达元件状态"] = "已设置"
    result["表达元件方案"] = {
        "方案编号": display_id,
        "排名": reference.get("rank"),
        "名称": str(design.get("name") or ""),
        "策略": str(design.get("strategy") or ""),
        "表达模式": str(reference.get("expression_regime") or ""),
        "成功评分": reference.get("expression_success_score"),
        "表达负担": reference.get("expression_burden"),
        "是否系统推荐": reference.get("system_recommended") is True,
        "来源文件": artifact.get("_resolved_path"),
    }
    result["完整表达构建"] = _construct_summary(
        manifest,
        project_root,
        selection,
        display_id,
        result["警告"],
    )


def get_expression_box_info(config: Any) -> dict[str, Any]:
    """Return the selected expression-box layout with one optional parts design."""

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    manifest = read_design_manifest(manifest_path)
    target = str(manifest.get("target_compound_id") or "").strip()
    if target and target != config.target_name:
        raise ValueError(
            f"manifest 目标化合物为 {target}，与当前输入目标 {config.target_name} 不一致"
        )
    requested = getattr(config, "parts_design", None)
    requested_design_id = (
        None if requested is None else _positive_int(requested, "--parts-design")
    )
    result: dict[str, Any] = {
        "运行成功": True,
        "状态": "当前无表达盒方案",
        "目标化合物": target or str(config.target_name),
        "表达盒方案": None,
        "表达元件状态": "未设置",
        "已选表达元件方案": [],
        "主表达元件方案": None,
        "当前显示方案": None,
        "表达元件方案": None,
        "表达元件草稿": None,
        "表达盒列表": [],
        "完整表达构建": None,
        "警告": [],
    }
    raw_selection = manifest.get("expression_box_selection")
    if raw_selection is None:
        if requested_design_id is not None:
            raise ValueError("尚未选择表达盒，不能使用 --parts-design")
        result["提示"] = (
            "当前无表达盒方案，请先运行 expression --design --box 并选择方案，"
            "或使用 --custom 直接写入。"
        )
        return result

    selection = _mapping(raw_selection, "expression_box_selection")
    if selection.get("schema_version") != "expression_box_selection.v1":
        raise ValueError("不支持的 expression_box_selection schema_version")
    fingerprint = _text(
        selection.get("selected_design_fingerprint"),
        "expression_box_selection.selected_design_fingerprint",
    )
    if _SHA256_RE.fullmatch(fingerprint) is None:
        raise ValueError("expression_box_selection 的方案 fingerprint 无效")
    cds_by_accession = _cds_by_accession(manifest)
    result["状态"] = "已选择表达盒方案"
    result["表达盒方案"] = {
        "选择状态": str(selection.get("selection_status") or ""),
        "方案编号": selection.get("selected_design_id"),
        "名称": str(selection.get("name") or ""),
        "策略": str(selection.get("strategy") or ""),
        "是否系统推荐": selection.get("recommended") is True,
        "方案Fingerprint": fingerprint,
    }
    result["表达盒列表"] = _base_cassettes(selection, cds_by_accession)
    _apply_parts(result, manifest, project_root, requested_design_id)
    return result


def _cell(value: Any) -> str:
    return " ".join(str(value).split())


def _width(value: str) -> int:
    return sum(2 if unicodedata.east_asian_width(char) in "WF" else 1 for char in value)


def _table(headers: list[str], rows: list[list[str]]) -> list[str]:
    widths = [
        max(_width(row[index]) for row in [headers, *rows])
        for index in range(len(headers))
    ]

    def render(row: list[str]) -> str:
        return " | ".join(
            cell + " " * (width - _width(cell))
            for cell, width in zip(row, widths, strict=True)
        ).rstrip()

    return [render(headers), "-+-".join("-" * width for width in widths), *map(render, rows)]


def _length(value: Any, unit: str) -> str:
    return "-" if value is None else f"{value} {unit}"


def _part_row(
    order: int,
    role: str,
    accession: str,
    part: Mapping[str, Any] | None,
) -> list[str]:
    return [
        str(order),
        role,
        accession or "-",
        _cell(part["元件ID"]) if part else "-",
        _length(part["长度_bp"], "bp") if part else "-",
        _cell(part["来源"]) if part else "-",
        "已设置" if part else "未设置",
    ]


def format_expression_box_info(result: Mapping[str, Any]) -> str:
    """Render a stable expression-box layout before or after part selection."""

    if "提示" in result:
        return str(result["提示"])
    scheme = _mapping(result.get("表达盒方案"), "表达盒方案")
    lines = [
        f"目标化合物：{result['目标化合物']}",
        f"表达盒方案：{scheme.get('名称') or scheme.get('方案编号') or '未命名'}",
        f"分组策略：{scheme.get('策略') or '未记录'}",
        f"表达元件状态：{result['表达元件状态']}",
        f"表达盒数量：{len(result['表达盒列表'])}",
    ]
    if result["已选表达元件方案"]:
        lines.append(
            "已选表达元件方案："
            + ", ".join(str(value) for value in result["已选表达元件方案"])
        )
        lines.append(f"主方案：{result['主表达元件方案']}")
        lines.append(f"当前显示方案：{result['当前显示方案']}")
    part_design = result.get("表达元件方案")
    if isinstance(part_design, Mapping):
        lines.append(
            f"成功评分：{part_design.get('成功评分')}｜"
            f"表达模式：{part_design.get('表达模式') or '未记录'}"
        )
    draft = result.get("表达元件草稿")
    if isinstance(draft, Mapping):
        lines.append(
            f"已上传启动子：{draft.get('已上传启动子数', 0)}/"
            f"{draft.get('表达盒总数', 0)}"
        )

    headers = ["组装顺序", "组件", "对应蛋白", "元件 ID", "长度", "来源", "状态"]
    for cassette in result["表达盒列表"]:
        lines.extend(["", f"表达盒 {cassette['表达盒编号']}"])
        rows: list[list[str]] = []
        order = 1
        rows.append(_part_row(order, "Promoter", "", cassette["Promoter"]))
        order += 1
        for gene in cassette["基因列表"]:
            rows.append(_part_row(order, "RBS", gene["蛋白ID"], gene["RBS"]))
            order += 1
            rows.append(
                [
                    str(order),
                    "CDS",
                    gene["蛋白ID"],
                    gene["蛋白ID"],
                    _length(gene["CDS"]["长度_nt"], "nt"),
                    "CDS选择",
                    "已就绪",
                ]
            )
            order += 1
        rows.append(_part_row(order, "Terminator", "", cassette["Terminator"]))
        lines.extend(_table(headers, rows))

        rbs_rows = [
            [
                gene["蛋白ID"],
                gene["RBS"]["元件ID"],
                str(gene["OSTIR"]["翻译起始率"]),
                str(gene["OSTIR"]["非预期起始位点数"]),
            ]
            for gene in cassette["基因列表"]
            if gene["RBS"] is not None and gene["OSTIR"] is not None
        ]
        if rbs_rows:
            lines.extend(["", "RBS 检查"])
            lines.extend(
                _table(
                    ["蛋白", "RBS ID", "翻译起始率", "非预期起始位点"],
                    rbs_rows,
                )
            )
        audit = cassette.get("序列检查")
        if isinstance(audit, Mapping):
            lines.append(
                "完整表达盒检查："
                f"{audit.get('状态')}｜长度 {_length(audit.get('长度_nt'), 'bp')}｜"
                f"GC {audit.get('GC百分比') if audit.get('GC百分比') is not None else '-'}%"
            )

    construct = result.get("完整表达构建")
    if isinstance(construct, Mapping):
        lines.extend(
            [
                "",
                f"完整表达构建：{construct.get('状态')}",
                f"GenBank：{construct.get('路径') or '未生成'}",
            ]
        )
    if result["警告"]:
        lines.extend(["", "警告："])
        lines.extend(f"- {warning}" for warning in result["警告"])
    return "\n".join(lines)


def run_expression_box_info(config: Any) -> dict[str, Any]:
    result = get_expression_box_info(config)
    print(format_expression_box_info(result))
    return result


__all__ = [
    "format_expression_box_info",
    "get_expression_box_info",
    "run_expression_box_info",
]
