"""Render a deterministic Chinese report from manifest assembly facts."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _items(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _split_ids(value: Any) -> list[str]:
    if isinstance(value, list):
        raw = value
    else:
        raw = str(value or "").replace(",", ";").split(";")
    return [str(item).strip() for item in raw if str(item).strip()]


def _table(value: Any) -> str:
    text = str(value if value not in (None, "") else "—")
    return text.replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def _join(value: Any) -> str:
    values = _split_ids(value)
    return ", ".join(values) if values else "—"


def build_report_payload(
    *,
    manifest: Mapping[str, Any],
    final_assembly: Mapping[str, Any],
    run_summary: Mapping[str, Any],
) -> dict[str, Any]:
    solution = _mapping(manifest.get("solution"))
    solution_summary = _mapping(solution.get("summary"))
    main_selection = _mapping(manifest.get("main_enzyme_selection"))
    cds_selection = _mapping(manifest.get("cds_selection"))
    expression_box = _mapping(manifest.get("expression_box_selection"))
    parts_selection = _mapping(manifest.get("parts_selection"))
    plasmid_selection = _mapping(manifest.get("plasmid_selection"))
    vector = _mapping(plasmid_selection.get("vector"))

    route_steps = []
    for item in _items(solution.get("steps")):
        if not isinstance(item, Mapping):
            continue
        route_steps.append(
            {
                "step_index": item.get("step_index"),
                "reaction_id": item.get("reaction_id"),
                "reaction_name": item.get("reaction_name"),
                "status": item.get("status"),
                "ko_ids": _split_ids(item.get("ko_ids")),
                "enzyme_ecs": _split_ids(item.get("enzyme_ecs")),
            }
        )
    route_steps.sort(key=lambda item: int(item.get("step_index") or 0))

    main_proteins = []
    for item in _items(main_selection.get("proteins")):
        if isinstance(item, Mapping):
            main_proteins.append(
                {
                    "accession": item.get("accession"),
                    "protein_name": item.get("protein_name"),
                    "organism_name": item.get("organism_name"),
                    "reviewed": item.get("reviewed"),
                    "assigned_step_indexes": _items(
                        item.get("assigned_step_indexes")
                    ),
                }
            )
    main_proteins.sort(
        key=lambda item: (
            min(item["assigned_step_indexes"] or [999]),
            str(item.get("accession") or ""),
        )
    )

    optimized_cds = []
    for item in _items(cds_selection.get("proteins")):
        if not isinstance(item, Mapping):
            continue
        optimized = _mapping(item.get("optimized_cds"))
        final_metrics = _mapping(_mapping(optimized.get("metrics")).get("final"))
        optimized_cds.append(
            {
                "accession": item.get("accession"),
                "protein_name": item.get("protein_name"),
                "assigned_step_indexes": _items(item.get("assigned_step_indexes")),
                "length_nt": optimized.get("length_nt"),
                "path": optimized.get("path"),
                "gc_percent": final_metrics.get("gc_percent"),
                "cai": final_metrics.get("cai"),
                "gate_status": final_metrics.get("gate_status"),
            }
        )
    optimized_cds.sort(
        key=lambda item: (
            min(item["assigned_step_indexes"] or [999]),
            str(item.get("accession") or ""),
        )
    )

    cassettes = [
        {
            "cassette_index": item.get("cassette_index"),
            "protein_accessions": _items(item.get("protein_accessions")),
            "protein_count": item.get("protein_count"),
            "total_cds_length_nt": item.get("total_cds_length_nt"),
        }
        for item in _items(expression_box.get("cassettes"))
        if isinstance(item, Mapping)
    ]
    design_references = [
        {
            "design_id": item.get("design_id"),
            "rank": item.get("rank"),
            "expression_success_score": item.get("expression_success_score"),
            "expression_regime": item.get("expression_regime"),
            "system_recommended": item.get("system_recommended"),
            "expression_burden": dict(_mapping(item.get("expression_burden"))),
        }
        for item in _items(parts_selection.get("design_references"))
        if isinstance(item, Mapping)
    ]
    design_references.sort(key=lambda item: int(item.get("design_id") or 0))

    constructs = [
        {
            "parts_design_id": item.get("parts_design_id"),
            "assembly_method": item.get("assembly_method"),
            "enzyme_summary": item.get("enzyme_summary"),
            "target": dict(_mapping(item.get("target"))),
            "length_bp": item.get("length_bp"),
            "sequence_sha256": item.get("sequence_sha256"),
            "files": dict(_mapping(item.get("files"))),
            "validation": dict(_mapping(item.get("validation"))),
            "warnings": [str(value) for value in _items(item.get("warnings"))],
        }
        for item in _items(final_assembly.get("constructs"))
        if isinstance(item, Mapping)
    ]
    constructs.sort(key=lambda item: int(item.get("parts_design_id") or 0))

    return {
        "target_compound_id": manifest.get("target_compound_id"),
        "solution_summary": dict(solution_summary),
        "route_steps": route_steps,
        "host": dict(_mapping(cds_selection.get("host"))),
        "host_gene_knockouts": [],
        "main_enzyme_selection": {
            "selection_status": main_selection.get("selection_status"),
            "proteins": main_proteins,
            "unresolved_reviews": [
                dict(item)
                for item in _items(main_selection.get("unresolved_reviews"))
                if isinstance(item, Mapping)
            ],
            "warnings": [
                str(item) for item in _items(main_selection.get("warnings"))
            ],
        },
        "cds_selection": {
            "status": cds_selection.get("status"),
            "proteins": optimized_cds,
            "warnings": [str(item) for item in _items(cds_selection.get("warnings"))],
        },
        "expression_box": {
            "selection_status": expression_box.get("selection_status"),
            "strategy": expression_box.get("strategy"),
            "name": expression_box.get("name"),
            "summary": dict(_mapping(expression_box.get("summary"))),
            "cassettes": cassettes,
        },
        "parts_selection": {
            "primary_design_id": parts_selection.get("primary_design_id"),
            "design_count": parts_selection.get("design_count"),
            "design_references": design_references,
            "warnings": [
                str(item) for item in _items(parts_selection.get("warnings"))
            ],
        },
        "plasmid": {
            "plasmid_id": vector.get("plasmid_id"),
            "name": vector.get("name"),
            "length_bp": vector.get("length_bp"),
            "topology": vector.get("topology"),
            "replicon_family": vector.get("replicon_family"),
            "copy_number": vector.get("copy_number"),
            "copy_number_class": vector.get("copy_number_class"),
            "marker": vector.get("marker"),
            "bacterial_resistance": vector.get("bacterial_resistance"),
            "assembly_policy": vector.get("assembly_policy"),
            "warnings": [
                str(item) for item in _items(plasmid_selection.get("warnings"))
            ],
        },
        "execution": {
            "status": final_assembly.get("status"),
            "planned_design_count": final_assembly.get("planned_design_count"),
            "succeeded_count": final_assembly.get("succeeded_count"),
            "failed_count": final_assembly.get("failed_count"),
            "method_counts": dict(_mapping(run_summary.get("method_counts"))),
            "output_dir": final_assembly.get("output_dir"),
        },
        "constructs": constructs,
        "failures": [
            dict(item)
            for item in _items(final_assembly.get("failures"))
            if isinstance(item, Mapping)
        ],
        "warnings": [
            str(item) for item in _items(final_assembly.get("warnings"))
        ],
    }


_REVIEW_LABELS = {
    "reaction_fit": "目标反应匹配仍需复核",
    "reaction_direction": "生理反应方向证据仍需复核",
    "reaction_specificity": "底物或产物特异性仍需复核",
    "electron_assessment": "电子供给或再生机制仍需补充",
}


def _target_text(target: Mapping[str, Any]) -> str:
    mode = str(target.get("mode") or "")
    if mode == "restriction_replace_retain_sites":
        return (
            f"替换 {target.get('replace_start_bp')}–{target.get('replace_end_bp')} bp；"
            f"insert {target.get('inserted_start_bp')}–{target.get('inserted_end_bp')} bp；"
            "保留酶切位点"
        )
    if mode == "insert_after":
        return (
            f"在 {target.get('insert_after_bp')} bp 后插入；"
            f"insert {target.get('inserted_start_bp')}–{target.get('inserted_end_bp')} bp"
        )
    if mode == "replace":
        return (
            f"替换 {target.get('replace_start_bp')}–{target.get('replace_end_bp')} bp；"
            f"insert {target.get('inserted_start_bp')}–{target.get('inserted_end_bp')} bp"
        )
    return mode or "—"


def _file_path(files: Mapping[str, Any], key: str) -> str:
    return str(_mapping(files.get(key)).get("path") or "")


def _render_report(payload: Mapping[str, Any]) -> str:
    summary = _mapping(payload.get("solution_summary"))
    main = _mapping(payload.get("main_enzyme_selection"))
    cds = _mapping(payload.get("cds_selection"))
    expression_box = _mapping(payload.get("expression_box"))
    parts = _mapping(payload.get("parts_selection"))
    plasmid = _mapping(payload.get("plasmid"))
    execution = _mapping(payload.get("execution"))
    target_name = summary.get("target_compound_name") or payload.get(
        "target_compound_id"
    )
    lines = [
        "# GLADE 最终理论组装报告",
        "",
        "> 本报告由程序根据 design manifest 和组装输出固定生成。报告描述计算设计，"
        "不代表湿实验组装、表达或目标产物合成已经成功。",
        "",
        "## 1. 设计与通路摘要",
        "",
        "| 项目 | 内容 |",
        "|---|---|",
        f"| 目标化合物 | {_table(target_name)} ({_table(payload.get('target_compound_id'))}) |",
        f"| 通路总步骤 | {_table(summary.get('total_steps'))} |",
        f"| 异源步骤 | {_table(summary.get('heterologous_steps'))} |",
        f"| 可达锚点 | {_table(summary.get('reachable_anchor_labels'))} |",
        f"| 反应记录解析状态 | {_table(summary.get('reaction_resolution_status'))} |",
        f"| 阻断反应数 | {_table(summary.get('blocking_reaction_count'))} |",
        "",
        "> KO 在本报告中始终表示 KEGG Orthology，不表示 gene knockout（基因敲除）。",
        "",
        "| 步骤 | 反应 | 反应名称 | 状态 | 异源酶 KEGG Orthology | EC |",
        "|---:|---|---|---|---|---|",
    ]
    for item in _items(payload.get("route_steps")):
        if isinstance(item, Mapping):
            lines.append(
                f"| {_table(item.get('step_index'))} | {_table(item.get('reaction_id'))} | "
                f"{_table(item.get('reaction_name'))} | {_table(item.get('status'))} | "
                f"{_table(_join(item.get('ko_ids')))} | {_table(_join(item.get('enzyme_ecs')))} |"
            )

    host = _mapping(payload.get("host"))
    lines.extend(
        [
            "",
            "## 2. 宿主与基因敲除",
            "",
            f"- 当前设计宿主：{_table(host.get('name'))}",
            f"- 宿主标识：{_table(host.get('chassis_key'))}",
            "- 宿主基因敲除：当前流程未进行宿主基因敲除分析，暂无敲除建议。",
            "",
            "## 3. 主酶选择与待复核风险",
            "",
            f"- 选择状态：`{_table(main.get('selection_status'))}`",
            "",
            "| 步骤 | UniProt | 蛋白 | 来源生物 | Reviewed |",
            "|---|---|---|---|---|",
        ]
    )
    for item in _items(main.get("proteins")):
        if isinstance(item, Mapping):
            reviewed = "是" if item.get("reviewed") else "否"
            lines.append(
                f"| {_table(_join(item.get('assigned_step_indexes')))} | "
                f"{_table(item.get('accession'))} | {_table(item.get('protein_name'))} | "
                f"{_table(item.get('organism_name'))} | {reviewed} |"
            )
    lines.extend(["", "### 未解决的复核项", ""])
    reviews = _items(main.get("unresolved_reviews"))
    if reviews:
        for item in reviews:
            if not isinstance(item, Mapping):
                continue
            label = _REVIEW_LABELS.get(
                str(item.get("review_type") or ""),
                str(item.get("review_type") or "未分类复核项"),
            )
            reason = str(item.get("reason") or "")
            if item.get("review_type") != "electron_assessment":
                reason = label
            lines.append(
                f"- 步骤 {_join(item.get('step_indexes'))}：{label}；"
                f"状态 `{_table(item.get('status'))}`"
                + (f"；{reason}" if reason and reason != label else "")
            )
    else:
        lines.append("- 无。")
    main_warnings = _items(main.get("warnings"))
    if main_warnings:
        lines.extend(["", "### 来源记录中的具体警告", ""])
        lines.extend(f"- {item}" for item in main_warnings)

    lines.extend(
        [
            "",
            "## 4. 已生成的优化 CDS",
            "",
            f"- CDS 处理状态：`{_table(cds.get('status'))}`",
            "- 下列优化序列已经生成并进入表达构建；实验前应复核并据此下单合成，"
            "无需重新从零进行密码子优化。",
            "",
            "| 步骤 | UniProt | 蛋白 | 长度 (nt) | GC (%) | CAI | Gate | 文件 |",
            "|---|---|---|---:|---:|---:|---|---|",
        ]
    )
    for item in _items(cds.get("proteins")):
        if isinstance(item, Mapping):
            lines.append(
                f"| {_table(_join(item.get('assigned_step_indexes')))} | "
                f"{_table(item.get('accession'))} | {_table(item.get('protein_name'))} | "
                f"{_table(item.get('length_nt'))} | {_table(item.get('gc_percent'))} | "
                f"{_table(item.get('cai'))} | {_table(item.get('gate_status'))} | "
                f"`{_table(item.get('path'))}` |"
            )

    box_summary = _mapping(expression_box.get("summary"))
    lines.extend(
        [
            "",
            "## 5. 表达盒分组与表达方案",
            "",
            f"- 分组策略：`{_table(expression_box.get('strategy'))}`",
            f"- 表达盒数量：{_table(box_summary.get('cassette_count'))}",
            f"- 蛋白数量：{_table(box_summary.get('protein_count'))}",
            "",
            "| 表达盒 | 蛋白 accession | 蛋白数 | CDS 总长度 (nt) |",
            "|---:|---|---:|---:|",
        ]
    )
    for item in _items(expression_box.get("cassettes")):
        if isinstance(item, Mapping):
            lines.append(
                f"| {_table(item.get('cassette_index'))} | "
                f"{_table(_join(item.get('protein_accessions')))} | "
                f"{_table(item.get('protein_count'))} | "
                f"{_table(item.get('total_cds_length_nt'))} |"
            )
    lines.extend(
        [
            "",
            "12 个 design 表示不同表达元件组合；评分是工程启发式指标，不是实验成功率。",
            "",
            "| Design | 排名 | 表达评分 | 表达模式 | 负担分数 | 负担等级 | 系统主推荐 |",
            "|---:|---:|---:|---|---:|---|---|",
        ]
    )
    for item in _items(parts.get("design_references")):
        if not isinstance(item, Mapping):
            continue
        burden = _mapping(item.get("expression_burden"))
        lines.append(
            f"| {_table(item.get('design_id'))} | {_table(item.get('rank'))} | "
            f"{_table(item.get('expression_success_score'))} | "
            f"{_table(item.get('expression_regime'))} | {_table(burden.get('score'))} | "
            f"{_table(burden.get('level'))} | "
            f"{'是' if item.get('system_recommended') else '否'} |"
        )

    lines.extend(
        [
            "",
            "## 6. 质粒骨架",
            "",
            "| 项目 | 内容 |",
            "|---|---|",
            f"| 名称 | {_table(plasmid.get('name'))} |",
            f"| ID | {_table(plasmid.get('plasmid_id'))} |",
            f"| 原始长度 | {_table(plasmid.get('length_bp'))} bp |",
            f"| 拓扑 | {_table(plasmid.get('topology'))} |",
            f"| 复制子家族 | {_table(plasmid.get('replicon_family'))} |",
            f"| 拷贝数等级 | {_table(plasmid.get('copy_number_class'))} |",
            f"| 筛选标记 | {_table(plasmid.get('marker'))} |",
            f"| 组装策略 | {_table(plasmid.get('assembly_policy'))} |",
            "",
            "## 7. 最终理论组装结果",
            "",
            f"- 状态：`{_table(execution.get('status'))}`",
            f"- 计划/成功/失败：{_table(execution.get('planned_design_count'))} / "
            f"{_table(execution.get('succeeded_count'))} / "
            f"{_table(execution.get('failed_count'))}",
            "",
            "| Design | 方法 | 酶/线性化 | 插入或替换位置 | 最终长度 | 验证 | 输出文件 |",
            "|---:|---|---|---|---:|---|---|",
        ]
    )
    for item in _items(payload.get("constructs")):
        if not isinstance(item, Mapping):
            continue
        files = _mapping(item.get("files"))
        file_text = "<br>".join(
            f"`{path}`"
            for path in (
                _file_path(files, "genbank"),
                _file_path(files, "fasta"),
                _file_path(files, "assembly_report"),
            )
            if path
        )
        validation = _mapping(item.get("validation"))
        lines.append(
            f"| {_table(item.get('parts_design_id'))} | "
            f"{_table(item.get('assembly_method'))} | "
            f"{_table(item.get('enzyme_summary'))} | "
            f"{_table(_target_text(_mapping(item.get('target'))))} | "
            f"{_table(item.get('length_bp'))} | {_table(validation.get('status'))} | "
            f"{file_text or '—'} |"
        )

    lines.extend(["", "## 8. 失败项和计算验证边界", ""])
    failures = _items(payload.get("failures"))
    if failures:
        for item in failures:
            if isinstance(item, Mapping):
                lines.append(
                    f"- Design {_table(item.get('parts_design_id'))}："
                    f"{_table(item.get('error_type'))}: {_table(item.get('message'))}"
                )
    else:
        lines.append("- 本次无组装执行失败项。")
    lines.extend(
        [
            "- PASS 仅表示输出序列、文件格式、环状拓扑、feature 坐标和 insert 一致性通过计算检查。",
            "- PASS 不代表酶具有目标活性、蛋白能够表达、质粒已经构建或宿主能够合成目标化合物。",
            "",
            "## 9. 实验复核重点",
            "",
            "- 复核主酶与目标反应的底物、产物、方向和特异性；优先处理上述 pending review。",
            "- 处理蛋白定位肽、辅因子和外部电子再生需求后，再确定最终表达版本。",
            "- 使用已生成的优化 CDS 文件进行合成前复核，不要把 KEGG Orthology 当作敲除基因。",
            "- 查阅 BamHI、SalI、AvrII 等实际采用酶的说明书，确认缓冲液和反应条件。",
            "- 通过测序、表达检测、酶活和产物检测完成湿实验验证。",
            "",
        ]
    )
    return "\n".join(lines)


def generate_final_design_report(
    payload: Mapping[str, Any],
) -> tuple[str, str, None]:
    """Return deterministic Markdown without invoking any language model."""

    return _render_report(payload).rstrip() + "\n", "system_template", None


__all__ = [
    "build_report_payload",
    "generate_final_design_report",
]
