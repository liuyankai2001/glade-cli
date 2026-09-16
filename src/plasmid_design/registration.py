"""Prepare an independently versioned export bundle in existing manifest shapes."""

from __future__ import annotations

from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path

from src.final_assemble_execute.common import (
    sha256_file,
    sha256_sequence,
    validate_written_outputs,
    write_fasta,
    write_genbank,
    write_json_atomic,
)
from src.final_assemble_plan.config import (
    ASSEMBLY_PLAN_ALGORITHM_VERSION,
    ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION,
)
from src.final_assemble_plan.get_final_assembly_context import (
    build_final_assembly_context,
)
from src.plasmid_selection.get_plasmid_context import stable_json_hash

DESIGN_VERSION = "component_plasmid_design.v4"
DESIGN_WARNING = "这些文件是计算拼接设计，实际质粒需制备并验证。"
PROTECTIVE_WARNING = (
    "制备片段两端已加入 6 bp 保护碱基，实验前请核对所选酶的末端切割条件。"
)


def _signed(payload: dict) -> dict:
    payload["selection_fingerprint"] = stable_json_hash(payload)
    return payload


def _file_info(
    path: Path, root: Path, label: str, file_format: str, sequence: str | None = None
) -> dict:
    info = {
        "path": path.relative_to(root).as_posix(),
        "format": file_format,
        "file_sha256": sha256_file(path),
        "label": label,
    }
    if sequence is not None:
        info.update(
            {
                "length_bp": len(sequence),
                "sequence_sha256": sha256_sequence(sequence),
                "sequence_content_sha256": sha256_sequence(sequence),
            }
        )
    return info


def export_design(
    snapshot, design, directory: Path, generation_id: str
) -> dict[str, dict]:
    """Write artifacts and construct all payloads, without updating the manifest."""
    source, root = snapshot.source, snapshot.source.project_output_path
    components = design.preview["components"]
    resistances = [
        (c, snapshot.catalog.get_resistance(c["module_id"]))
        for c in components
        if c["component_type"] == "resistance"
    ]
    resistance = resistances[0][1]
    replication = snapshot.catalog.get_replication(
        next(c["module_id"] for c in components if c["component_type"] == "replication")
    )
    expression_segment = next(
        s for s in design.preview["segments"] if s["kind"] == "expression"
    )
    resistance_description = " + ".join(
        f"{module.name}（{module.antibiotic}）" for _, module in resistances
    )
    design_id = source.constructs[0].design_id
    generated_at = datetime.now(UTC).isoformat()
    warnings = list(
        dict.fromkeys([DESIGN_WARNING, PROTECTIVE_WARNING, *design.preview["warnings"]])
    )
    preparation_files = {}
    backbone_path = directory / "backbone.gb"
    write_genbank(backbone_path, design.backbone_record)
    backbone_file = _file_info(
        backbone_path, root, "骨架 GenBank", "genbank", str(design.backbone_record.seq)
    )
    for key, label in [("backbone", "骨架制备片段"), ("insert", "表达插入制备片段")]:
        record = design.preparation_records[key]
        prefix = "backbone_preparation" if key == "backbone" else "insert_preparation"
        for file_format, suffix in [("fasta", ".fasta"), ("genbank", ".gb")]:
            path = directory / (prefix + suffix)
            if file_format == "fasta":
                write_fasta(path, record.id, str(record.seq))
            else:
                write_genbank(path, record)
            file_id = prefix if file_format == "fasta" else prefix + "_genbank"
            preparation_files[file_id] = _file_info(
                path,
                root,
                label + (" FASTA" if file_format == "fasta" else " GenBank"),
                file_format,
                str(record.seq),
            )

    region = deepcopy(design.plan["target"]["insertion_region"])
    protected = []
    for start, end in (
        (1, region["start_bp"] - 1),
        (region["end_bp"] + 1, len(design.backbone_record)),
    ):
        if start <= end:
            protected.append(
                {
                    "label": "immutable_component_backbone",
                    "start_bp": start,
                    "end_bp": end,
                }
            )
    source_provenance = dict(snapshot.catalog.scaffold["provenance"])
    # Mapping proxies are thawed via JSON-safe source fields, never mutated.
    source_provenance = {
        k: dict(v) if hasattr(v, "items") else v for k, v in source_provenance.items()
    }
    component_design = {
        "generation_id": generation_id,
        "algorithm_version": DESIGN_VERSION,
        "source_fingerprint": snapshot.fingerprint,
        "catalog_fingerprint": snapshot.catalog.fingerprint,
        "resistance_id": resistance.id,
        "replication_id": replication.id,
        "t0_id": design.preview["t0_id"],
        "t1_id": design.preview["t1_id"],
        "terminator_warnings": list(design.preview["terminator_warnings"]),
        "component_order": list(design.preview["component_order"]),
        "components": deepcopy(components),
        "scaffold_version": snapshot.catalog.scaffold.get(
            "version", "basic_seva_fixed.v1"
        ),
        "preparation_files": preparation_files,
    }
    selection = _signed(
        {
            "schema_version": "plasmid_selection.v2",
            "status": "selected",
            "selection_status": "user_selected",
            "source": "web_component_design",
            "backbone_policy": "one_backbone_for_all_selected_constructs",
            "covered_parts_design_ids": [design_id],
            "covered_design_count": 1,
            "vector": {
                "plasmid_id": f"component_{generation_id}",
                "name": " + ".join(
                    [*(module.name for _, module in resistances), replication.name]
                ),
                "description": "Component-derived backbone for restriction assembly",
                "vector_type": "component_derived",
                "length_bp": len(design.backbone_record),
                "topology": "circular",
                "copy_number_class": replication.copy_number,
                "assembly_policy": "insert_into_mcs",
                "cargo_type": "temporary_landing_pad",
                "requires_cargo_replacement": True,
                "host_compatibility": [source.host_key],
                "origins": [
                    {
                        "module_id": replication.id,
                        "name": replication.name,
                        "host_range": replication.host_range,
                    }
                ],
                "selection_markers": [
                    {
                        "instance_id": component["instance_id"],
                        "module_id": module.id,
                        "name": module.name,
                        "antibiotic": module.antibiotic,
                        "gene": module.resistance_gene,
                    }
                    for component, module in resistances
                ],
                "insertion_regions": [region],
                "protected_features": protected,
                "selected_sequence_file": backbone_file,
                "source": {
                    "source": "curated_BASIC_SEVA_modules",
                    "source_provenance": source_provenance,
                    "experimental_validation": "not_demonstrated_for_this_new_design",
                },
            },
            "provenance": {
                "context_input_fingerprint": source.input_fingerprint,
                "parts_selection_fingerprint": source.parts_selection_fingerprint,
                "assembled_constructs_fingerprint": source.assembled_constructs_fingerprint,
            },
            "component_design": component_design,
            "warnings": warnings,
        }
    )
    final_context = build_final_assembly_context(source, selection)
    construct = final_context.constructs[0]
    plan = deepcopy(design.plan)
    plan.update(
        {
            "backbone": {
                "path": backbone_file["path"],
                "length_bp": final_context.backbone.length_bp,
                "file_sha256": final_context.backbone.file_sha256,
                "sequence_content_sha256": final_context.backbone.sequence_content_sha256,
            },
            "insert": {
                "path": construct.path.relative_to(root).as_posix(),
                "length_bp": construct.length_bp,
                "file_sha256": construct.file_sha256,
                "sequence_sha256": construct.sequence_sha256,
            },
            "selection_reason": "user_selected_modules_and_manifest_enzymes",
            "warnings": warnings,
        }
    )
    plan["plan_fingerprint"] = stable_json_hash(plan)
    plan_source = {
        "context_input_fingerprint": final_context.input_fingerprint,
        "parts_selection_fingerprint": source.parts_selection_fingerprint,
        "assembled_constructs_fingerprint": source.assembled_constructs_fingerprint,
        "plasmid_selection_fingerprint": selection["selection_fingerprint"],
        "request_fingerprint": stable_json_hash(
            {
                "source": snapshot.fingerprint,
                "resistance": resistance.id,
                "replication": replication.id,
                "component_order": design.preview["component_order"],
                "components": components,
                "t0_id": design.preview["t0_id"],
                "t1_id": design.preview["t1_id"],
                "enzymes": design.preview["enzymes"],
            }
        ),
        "plan_set_fingerprint": stable_json_hash([plan]),
    }
    plan_path = directory / "assembly_plan.json"
    plan_artifact = {
        "schema_version": ASSEMBLY_PLAN_RECOMMENDATIONS_SCHEMA_VERSION,
        "algorithm_version": ASSEMBLY_PLAN_ALGORITHM_VERSION,
        "source_algorithm_version": DESIGN_VERSION,
        "status": "complete",
        "target_compound_id": source.target_compound_id,
        "selection_mode": "user_selected_components",
        "requested_method": "restriction",
        "design_count": 1,
        "plans": [plan],
        "plan_set_fingerprint": plan_source["plan_set_fingerprint"],
        "source": plan_source,
        "warnings": warnings,
    }
    write_json_atomic(plan_path, plan_artifact)
    plan_file = _file_info(plan_path, root, "组装计划 JSON", "json")
    plan_selection = _signed(
        {
            "schema_version": "final_assembly_plan.v2",
            "status": "selected",
            "selection_mode": "user_selected_components",
            "requested_method": "restriction",
            "design_count": 1,
            "design_plans": [plan],
            "source": {
                **plan_source,
                "artifact": plan_file["path"],
                "artifact_schema_version": plan_artifact["schema_version"],
                "algorithm_version": ASSEMBLY_PLAN_ALGORITHM_VERSION,
                "artifact_file": plan_file,
            },
            "warnings": warnings,
        }
    )
    final_sequence = str(design.final_record.seq)
    final_genbank = directory / "final_plasmid.gb"
    final_fasta = directory / "final_plasmid.fasta"
    write_genbank(final_genbank, design.final_record)
    write_fasta(final_fasta, design.final_record.id, final_sequence)
    inserted_start = expression_segment["start_bp"]
    validation = validate_written_outputs(
        genbank_path=final_genbank,
        fasta_path=final_fasta,
        expected_sequence=final_sequence,
        inserted_start_bp=inserted_start,
        inserted_end_bp=inserted_start + len(snapshot.insert_record) - 1,
        expected_insert=str(snapshot.insert_record.seq),
    )
    final_files = {
        "genbank": _file_info(
            final_genbank, root, "最终质粒 GenBank", "genbank", final_sequence
        ),
        "fasta": _file_info(
            final_fasta, root, "最终质粒 FASTA", "fasta", final_sequence
        ),
    }
    assembly_path = directory / "assembly.json"
    write_json_atomic(
        assembly_path,
        {
            "schema_version": "final_assembly_design_report.v1",
            "status": "assembled_in_silico",
            "algorithm_version": DESIGN_VERSION,
            "target_compound_id": source.target_compound_id,
            "parts_design_id": design_id,
            "source": {
                "final_assembly_plan_selection_fingerprint": plan_selection[
                    "selection_fingerprint"
                ],
                "plan_fingerprint": plan["plan_fingerprint"],
                "backbone_file_sha256": final_context.backbone.file_sha256,
                "insert_file_sha256": construct.file_sha256,
            },
            "assembly_method": "restriction",
            "backbone_linearization": plan["backbone_linearization"],
            "target": plan["target"],
            "final_construct": {
                "length_bp": len(final_sequence),
                "topology": "circular",
                "sequence_sha256": sha256_sequence(final_sequence),
                "files": final_files,
            },
            "validation": validation,
            "warnings": warnings,
        },
    )
    final_files["assembly_report"] = _file_info(
        assembly_path, root, "组装记录 JSON", "json"
    )
    final_construct = {
        "parts_design_id": design_id,
        "status": "assembled_in_silico",
        "assembly_method": "restriction",
        "enzyme_summary": plan["backbone_linearization"]["enzyme_summary"],
        "plan_fingerprint": plan["plan_fingerprint"],
        "length_bp": len(final_sequence),
        "topology": "circular",
        "sequence_sha256": sha256_sequence(final_sequence),
        "feature_count": len(design.final_record.features),
        "files": final_files,
        "validation": validation,
        "warnings": warnings,
    }
    bundle_fingerprint = stable_json_hash(
        {
            "plan_selection_fingerprint": plan_selection["selection_fingerprint"],
            "constructs": [final_construct],
            "algorithm_version": DESIGN_VERSION,
        }
    )
    report_path = directory / "design_report_zh.md"
    rows = "\n".join(
        f"| {e['role']} | {e['name']} | {e['site']} | {e['overhang']} |"
        for e in design.preview["enzymes"]
    )
    order_labels = {
        "resistance": "抗性模块",
        "replication": "复制模块",
        "expression": "完整表达构建",
        "t0": "T0",
        "t1": "T1",
    }
    order_description = " → ".join(
        f"{order_labels[c['component_type']]}（{c.get('module_id', '当前项目')}；{c['instance_id']}）"
        for c in components
    )
    terminator_description = (
        "；".join(
            f"{c['component_type'].upper()} = {c['module_id']}（{c['instance_id']}）"
            for c in components
            if c["component_type"] in ("t0", "t1")
        )
        or "未选择"
    )
    report_path.write_text(
        f"# 质粒设计：{source.target_compound_id}\n\n抗性模块：{resistance_description}\n\n复制模块：{replication.name}；拷贝数类别：{replication.copy_number}\n\n表达构建：方案 {design_id}，{construct.length_bp} bp\n\n目标质粒：{len(final_sequence)} bp，GC {design.preview['gc_percent']:.2f}%，环状。\n\n保留来源中的模块间隔；终止子仅包含用户明确选入的 T0/T1。所选模块和表达构建序列完整保留。\n\n| 边界 | 酶 | 识别序列 | 末端 |\n| --- | --- | --- | --- |\n{rows}\n\n制备骨架和表达插入两个 DNA 片段，再用指定双酶切割并连接。骨架制备片段中的边界酶顺序为右酶→骨架→左酶，插入片段为左酶→表达构建→右酶。保护碱基在切割时去除；骨架占位序列不保留在目标质粒中。\n\n禁止位点检查覆盖各模块、表达构建、连接处及环状闭合处；只保留两个指定边界位点。序列、导出文件和插入区域已核对。\n\n"
        + f"组件顺序：{order_description}\n\n"
        + f"终止子引用：{terminator_description}。\n\n"
        + "\n".join(f"- {w}" for w in warnings)
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    report_file = _file_info(report_path, root, "中文设计报告", "markdown")
    summary_path = directory / "run_summary.json"
    write_json_atomic(
        summary_path,
        {
            "schema_version": "final_assembly_run_summary.v1",
            "status": "complete",
            "generated_at": generated_at,
            "target_compound_id": source.target_compound_id,
            "planned_design_count": 1,
            "succeeded_count": 1,
            "failed_count": 0,
            "method_counts": {"restriction": 1},
            "constructs": [final_construct],
            "failures": [],
            "bundle_fingerprint": bundle_fingerprint,
            "warnings": warnings,
        },
    )
    summary_file = _file_info(summary_path, root, "生成汇总 JSON", "json")
    final_section = _signed(
        {
            "schema_version": "final_assembly.v2",
            "status": "complete",
            "result_kind": "in_silico_theoretical_assembly",
            "source": "web_component_design",
            "algorithm_version": DESIGN_VERSION,
            "generated_at": generated_at,
            "target_compound_id": source.target_compound_id,
            "planned_design_count": 1,
            "succeeded_count": 1,
            "failed_count": 0,
            "method_counts": {"restriction": 1},
            "output_dir": directory.relative_to(root).as_posix(),
            "source_fingerprints": {
                "final_assembly_plan_selection_fingerprint": plan_selection[
                    "selection_fingerprint"
                ],
                "plasmid_selection_fingerprint": selection["selection_fingerprint"],
                "assembled_constructs_fingerprint": source.assembled_constructs_fingerprint,
            },
            "constructs": [final_construct],
            "failures": [],
            "bundle_fingerprint": bundle_fingerprint,
            "warnings": warnings,
            "run_summary_file": summary_file,
            "report_file": report_file,
        }
    )
    design_report = _signed(
        {
            "schema_version": "final_design_report.v2",
            "status": "exported",
            "source": "web_component_design",
            "generated_at": generated_at,
            "language": "zh-CN",
            "generated_by": "deterministic_template",
            "report_file": report_file,
            "source_final_assembly_bundle_fingerprint": bundle_fingerprint,
            "warnings": [],
        }
    )
    return {
        "plasmid_selection": selection,
        "final_assembly_plan": plan_selection,
        "final_assembly": final_section,
        "final_design_report": design_report,
    }


def committed_files(manifest: dict) -> dict[str, dict]:
    selection = manifest.get("plasmid_selection", {})
    final = manifest.get("final_assembly", {})
    files = {}
    backbone_file = selection.get("vector", {}).get("selected_sequence_file")
    if isinstance(backbone_file, dict):
        files["backbone_genbank"] = backbone_file
    files.update(selection.get("component_design", {}).get("preparation_files", {}))
    constructs = final.get("constructs", [])
    if len(constructs) == 1:
        for key, info in constructs[0].get("files", {}).items():
            files[
                {
                    "genbank": "final_genbank",
                    "fasta": "final_fasta",
                    "assembly_report": "assembly_json",
                }.get(key, key)
            ] = info
    for file_id, key in [
        ("report", "report_file"),
        ("run_summary", "run_summary_file"),
    ]:
        if isinstance(final.get(key), dict):
            files[file_id] = final[key]
    plan_file = (
        manifest.get("final_assembly_plan", {}).get("source", {}).get("artifact_file")
    )
    if isinstance(plan_file, dict):
        files["plan_json"] = plan_file
    return files
