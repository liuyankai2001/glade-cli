"""Immutable module assembly and independent circular restriction-site checks."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from itertools import product
from typing import Any

from Bio import Restriction
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from src.final_assemble_execute.common import (
    assemble_sequence,
    build_final_record,
    sha256_sequence,
)
from src.final_assemble_execute.models import SequenceAssemblyResult
from src.final_assemble_plan.common import enzyme_cut_positions, gc_percent
from src.plasmid_design.errors import DesignError
from src.plasmid_design.components import legacy_components, normalize_components

LEGACY_COMPONENT_ORDER = ("resistance", "replication", "expression")
DEFAULT_COMPONENT_ORDER = ("resistance", "replication", "t1", "expression", "t0")


def normalize_component_order(order=DEFAULT_COMPONENT_ORDER) -> list[str]:
    if (
        isinstance(order, (list, tuple))
        and len(order) == 3
        and all(isinstance(item, str) for item in order)
        and set(order) == set(LEGACY_COMPONENT_ORDER)
    ):
        return [
            piece
            for kind in order
            for piece in (
                ("t1", "expression", "t0") if kind == "expression" else (kind,)
            )
        ]
    if (
        not isinstance(order, (list, tuple))
        or len(order) != 5
        or any(not isinstance(item, str) for item in order)
        or set(order) != set(DEFAULT_COMPONENT_ORDER)
    ):
        raise DesignError(
            "组件顺序必须包含抗性、复制、完整表达构建、T0、T1 各一次。",
            code="invalid_component_order",
        )
    return list(order)


@dataclass(frozen=True, slots=True)
class EnzymeChoice:
    name: str
    site: str
    role: str
    enzyme: Any

    def to_payload(self) -> dict:
        prime = "5′" if self.enzyme.is_5overhang() else "3′"
        return {
            "name": self.name,
            "site": self.site,
            "role": self.role,
            "overhang": f"{prime} {self.enzyme.ovhgseq}",
        }


def resolve_enzymes(names) -> tuple[EnzymeChoice, EnzymeChoice]:
    index = {name.lower(): name for name in rest_dict}
    if not isinstance(names, (list, tuple)):
        raise DesignError(
            "请先用 optimize --enzyme 指定两种组装酶。", code="enzyme_count"
        )
    canonical = set()
    for name in names:
        actual = index.get(name.strip().lower()) if isinstance(name, str) else None
        if actual is None:
            raise DesignError(
                f"无法识别限制酶：{name}；请检查前面的酶设置。", code="unknown_enzyme"
            )
        canonical.add(actual)
    if len(canonical) != 2:
        raise DesignError(
            "首版使用双酶组装，请先用 optimize --enzyme 设置恰好两种不同的酶。",
            code="enzyme_count",
        )
    choices = []
    for role, name in zip(
        ("left", "right"), sorted(canonical, key=str.lower), strict=True
    ):
        enzyme = getattr(Restriction, name)
        site = str(enzyme.site).upper()
        n = len(site)
        if (
            not site
            or set(site) - set("ACGT")
            or not enzyme.is_palindromic()
            or not enzyme.is_defined()
            or enzyme.is_blunt()
            or enzyme.scd5 is not None
            or enzyme.scd3 is not None
            or not 0 <= enzyme.fst5 <= n
            or not 0 <= n + enzyme.fst3 <= n
        ):
            raise DesignError(
                f"{name} 不适用于首版常规双酶组装，请更换为切点在识别序列内、末端确定的黏性末端酶。",
                code="unsupported_enzyme",
            )
        choices.append(EnzymeChoice(name, site, role, enzyme))
    if choices[0].enzyme % choices[1].enzyme:
        raise DesignError(
            f"{choices[0].name} 与 {choices[1].name} 产生兼容末端，不能保证插入方向，请更换其中一种酶。",
            code="compatible_ends",
        )
    return tuple(choices)


def audit_sequence(
    sequence: str,
    enzymes,
    *,
    circular: bool,
    allowed: dict[str, set[int]] | None = None,
) -> list[dict]:
    """Look for overlapping motifs on both strands, including circular closure."""
    sequence = sequence.upper()
    allowed = allowed or {}
    issues = []
    for choice in enzymes:
        motifs = {choice.site, str(Seq(choice.site).reverse_complement())}
        starts = set()
        for motif in motifs:
            extended = sequence + sequence[: len(motif) - 1] if circular else sequence
            start = extended.find(motif)
            while start >= 0:
                if start < len(sequence):
                    starts.add(start + 1)
                start = extended.find(motif, start + 1)
        for start in sorted(starts - allowed.get(choice.name, set())):
            end = start + len(choice.site) - 1
            issues.append(
                {
                    "code": "restriction_conflict",
                    "message": f"{choice.name} 存在非边界识别位点（{start}–{end} bp）。",
                    "enzyme": choice.name,
                    "start_bp": start,
                    "end_bp": end,
                }
            )
    return issues


def feature_payloads(record: SeqRecord) -> list[dict]:
    result = []
    for feature in record.features:
        if feature.type == "source" or feature.location is None:
            continue
        label = next(
            (
                str(feature.qualifiers[k][0])
                for k in ("label", "gene", "product")
                if feature.qualifiers.get(k)
            ),
            feature.type,
        )
        kind = str(feature.qualifiers.get("web_kind", [""])[0])
        if not kind:
            if feature.type == "CDS":
                kind = "gene"
            elif feature.type == "terminator" or "SEVA_T" in label:
                kind = "terminator"
            elif "retained_site" in label:
                kind = "restriction"
            else:
                kind = (
                    "expression"
                    if feature.qualifiers.get("assembly_source")
                    == ["complete_expression_construct"]
                    else "linker"
                )
        result.append(
            {
                "label": label,
                "type": feature.type,
                "start_bp": int(feature.location.start) + 1,
                "end_bp": int(feature.location.end),
                "strand": int(feature.location.strand or 0),
                "kind": kind,
                **(
                    {"instance_id": "".join(feature.qualifiers["web_instance_id"])}
                    if feature.qualifiers.get("web_instance_id")
                    else {}
                ),
            }
        )
    return result


def _annotation(start, end, label, kind, *, feature_type="misc_feature"):
    return SeqFeature(
        SimpleLocation(start, end, strand=1),
        type=feature_type,
        qualifiers={"label": [label], "web_kind": [kind]},
    )


def _module_features(module, offset: int, instance_id: str) -> list[SeqFeature]:
    features = [
        _annotation(offset, offset + module.length_bp, module.name, module.type)
    ]
    for item in module.features:
        parts = [
            SimpleLocation(
                int(p["start"]) + offset, int(p["end"]) + offset, strand=p.get("strand")
            )
            for p in item["parts"]
        ]
        location = (
            parts[0]
            if len(parts) == 1
            else CompoundLocation(parts, operator=item.get("operator") or "join")
        )
        qualifiers = {
            key: list(values) for key, values in item.get("qualifiers", {}).items()
        }
        features.append(SeqFeature(location, type=item["type"], qualifiers=qualifiers))
    # Repeat the qualifier in short ordered chunks to avoid spaces inserted by
    # GenBank readers when wrapping long qualifier words.
    chunks = [
        instance_id[index : index + 36] for index in range(0, len(instance_id), 36)
    ]
    for feature in features:
        feature.qualifiers["web_instance_id"] = list(chunks)
    return features


@dataclass(slots=True)
class MolecularDesign:
    backbone_record: SeqRecord
    final_record: SeqRecord
    plan: dict
    preview: dict
    preparation_records: dict[str, SeqRecord]


def _locate_issues(issues, segments, size):
    for issue in issues:
        start, end = issue["start_bp"], issue["end_bp"]
        label = "环状闭合边界" if end > size else "模块连接边界"
        for segment in segments:
            if segment["start_bp"] <= start <= end <= segment["end_bp"]:
                label = segment["label"]
                break
        issue["module"] = label
        issue["message"] = (
            f"{issue['enzyme']} 在{label}存在额外位点（质粒坐标 {start}–{end} bp）；请更换酶或模块。"
        )


def _prepare_fragment(
    payload: str, template: SeqRecord, offset: int, enzymes, record_id: str
) -> SeqRecord:
    # These bases are outside retained modules and are removed by digestion.
    for bases in product("ACGT", repeat=6):
        guard = "".join(bases)
        if not 2 <= sum(b in "GC" for b in guard) <= 4:
            continue
        sequence = guard + payload + guard
        allowed = {e.name: {payload.find(e.site) + 7} for e in enzymes}
        if audit_sequence(sequence, enzymes, circular=False, allowed=allowed):
            continue
        if any(
            len(enzyme_cut_positions(e.enzyme, sequence, circular=False)) != 1
            for e in enzymes
        ):
            continue
        record = SeqRecord(
            Seq(sequence),
            id=record_id,
            name=record_id[:16],
            description="Linear preparation fragment with six protective bases at each end",
        )
        record.annotations = {
            "molecule_type": "DNA",
            "topology": "linear",
            "data_file_division": "SYN",
        }
        for feature in template.features:
            if feature.type == "source" or feature.location is None:
                continue
            copied = deepcopy(feature)
            copied.location = copied.location + offset + 6
            record.features.append(copied)
        record.features.extend(
            [
                _annotation(0, 6, "5prime_protective_bases", "linker"),
                _annotation(
                    len(sequence) - 6,
                    len(sequence),
                    "3prime_protective_bases",
                    "linker",
                ),
            ]
        )
        for e in enzymes:
            start = payload.find(e.site) + 6
            record.features.append(
                _annotation(
                    start,
                    start + len(e.site),
                    f"{e.name}_preparation_site",
                    "restriction",
                )
            )
        return record
    raise DesignError(
        "无法为制备片段找到不引入额外酶位点的保护碱基，请更换酶。",
        code="preparation_conflict",
    )


def build_design(
    catalog,
    insert_record: SeqRecord,
    enzyme_names,
    resistance_id: str | None = None,
    replication_id: str | None = None,
    design_id: int = 1,
    component_order=DEFAULT_COMPONENT_ORDER,
    t0_id: str | None = None,
    t1_id: str | None = None,
    components: list[dict] | None = None,
) -> MolecularDesign:
    legacy_order = (
        normalize_component_order(component_order) if components is None else None
    )
    components = normalize_components(
        (
            legacy_components(
                {
                    "resistance_id": resistance_id,
                    "replication_id": replication_id,
                    "t0_id": t0_id,
                    "t1_id": t1_id,
                },
                legacy_order,
            )
            if components is None
            else components
        ),
        catalog,
    )
    order = legacy_order or [c["instance_id"] for c in components]
    enzymes = resolve_enzymes(enzyme_names)
    modules = {}
    for component in components:
        kind = component["component_type"]
        if kind == "expression":
            continue
        identifier = component["module_id"]
        modules[component["instance_id"]] = (
            catalog.get_terminator(identifier)
            if kind in ("t0", "t1")
            else (
                catalog.get_resistance(identifier)
                if kind == "resistance"
                else catalog.get_replication(identifier)
            )
        )
    resistances = [
        modules[c["instance_id"]]
        for c in components
        if c["component_type"] == "resistance"
    ]
    resistance = resistances[0]
    replication = next(
        modules[c["instance_id"]]
        for c in components
        if c["component_type"] == "replication"
    )
    t0_id = next(
        (c["module_id"] for c in components if c["component_type"] == "t0"), None
    )
    t1_id = next(
        (c["module_id"] for c in components if c["component_type"] == "t1"), None
    )
    terminator_warnings = []
    actual_order = [c["component_type"] for c in components]
    expression_index = actual_order.index("expression")
    for role, relative in (("t1", -1), ("t0", 1)):
        boundary = "上游" if relative < 0 else "下游"
        if role not in actual_order:
            terminator_warnings.append(
                f"当前设计缺少 {role.upper()} 终止子，仍可继续生成。"
            )
        elif actual_order[(expression_index + relative) % len(actual_order)] != role:
            terminator_warnings.append(
                f"{role.upper()} 未位于完整表达构建的{boundary}边界，仍可继续生成。"
            )
    insert = str(insert_record.seq).upper()
    if not insert or set(insert) - set("ACGT"):
        raise DesignError(
            "完整表达构建必须包含连续的 A/C/G/T 序列。", code="invalid_sequence"
        )
    left, right = enzymes
    spacer = catalog.scaffold["landing_pad_spacer"]

    def layout(expression_dna, *, placeholder=False):
        expression_pieces = [
            ("left_site", left.name, "restriction", left.site),
            (
                "landing_pad_spacer" if placeholder else "expression",
                "landing_pad_spacer" if placeholder else "完整表达构建",
                "linker" if placeholder else "expression",
                expression_dna,
            ),
            ("right_site", right.name, "restriction", right.site),
        ]
        pieces = []
        for component in components:
            owner, kind = component["instance_id"], component["component_type"]
            if kind == "expression":
                block = expression_pieces
            else:
                module = modules[owner]
                block = [
                    (
                        module.id,
                        module.name,
                        "terminator" if kind in ("t0", "t1") else kind,
                        module.sequence,
                    )
                ]
                if kind == "resistance":
                    block.append(
                        (
                            "module_interval",
                            "固定模块间隔",
                            "linker",
                            catalog.scaffold["resistance_to_replication"],
                        )
                    )
                elif kind == "replication":
                    block.append(
                        (
                            "t1_interval",
                            "复制模块尾部间隔",
                            "linker",
                            catalog.scaffold["replication_to_t1"],
                        )
                    )
            pieces.extend((owner, kind, *piece) for piece in block)
        segments, features, cursor = [], [], 0
        for owner, component_type, identifier, label, kind, dna in pieces:
            segments.append(
                {
                    "id": identifier,
                    "label": label,
                    "kind": kind,
                    "component_type": component_type,
                    "instance_id": owner,
                    "start_bp": cursor + 1,
                    "end_bp": cursor + len(dna),
                    "length_bp": len(dna),
                }
            )
            if kind in ("resistance", "replication"):
                features.extend(_module_features(modules[owner], cursor, owner))
            elif kind == "terminator":
                # One source annotation per standalone terminator; no duplicate wrapper.
                features.extend(_module_features(modules[owner], cursor, owner)[1:])
            elif kind not in ("restriction", "expression"):
                features.append(
                    _annotation(
                        cursor,
                        cursor + len(dna),
                        label,
                        kind,
                        feature_type=(
                            "terminator" if kind == "terminator" else "misc_feature"
                        ),
                    )
                )
            cursor += len(dna)
        return "".join(p[5] for p in pieces), segments, features

    backbone_sequence, bone_segments, features = layout(spacer, placeholder=True)
    final_sequence, segments, _ = layout(insert)
    bone_by_id = {s["id"]: s for s in bone_segments}
    final_by_id = {s["id"]: s for s in segments}
    replace_start = bone_by_id["left_site"]["start_bp"] - 1
    replace_end = bone_by_id["right_site"]["end_bp"]
    backbone = SeqRecord(
        Seq(backbone_sequence),
        id="designed_backbone",
        name="designed_bb",
        description=f"Component-derived {resistance.name} / {replication.name} backbone with restriction landing pad",
    )
    backbone.annotations = {
        "molecule_type": "DNA",
        "topology": "circular",
        "data_file_division": "SYN",
    }
    backbone.features = features
    site_starts = {
        left.name: bone_by_id["left_site"]["start_bp"],
        right.name: bone_by_id["right_site"]["start_bp"],
    }
    plan = {
        "parts_design_id": design_id,
        "assembly_method": "restriction",
        "component_order": order,
        "components": components,
        "t0_id": t0_id,
        "t1_id": t1_id,
        "estimated_final_length_bp": len(final_sequence),
        "target": {
            "mode": "replace",
            "replace_start_bp": replace_start + 1,
            "replace_end_bp": replace_end,
            "insertion_region": {
                "label": "designed_landing_pad",
                "start_bp": replace_start + 1,
                "end_bp": replace_end,
            },
        },
        "backbone_linearization": {
            "mode": "restriction",
            "restriction_enzymes": [],
            "enzyme_summary": f"{left.name}/{right.name}",
        },
        "restriction": {
            "left_enzyme": left.name,
            "right_enzyme": right.name,
            "left_site": left.site,
            "right_site": right.site,
            "restriction_site_retention": "retain",
        },
        "warnings": [],
    }
    bone_issues = audit_sequence(
        backbone_sequence,
        enzymes,
        circular=True,
        allowed={name: {start} for name, start in site_starts.items()},
    )
    _locate_issues(bone_issues, bone_segments, len(backbone_sequence))
    issues = audit_sequence(
        final_sequence,
        enzymes,
        circular=True,
        allowed={
            left.name: {final_by_id["left_site"]["start_bp"]},
            right.name: {final_by_id["right_site"]["start_bp"]},
        },
    )
    _locate_issues(issues, segments, len(final_sequence))
    for issue in bone_issues:
        if not any(
            (i["enzyme"], i["start_bp"], i["end_bp"])
            == (issue["enzyme"], issue["start_bp"], issue["end_bp"])
            for i in issues
        ):
            issues.append(issue)
    for e in enzymes:
        cuts = enzyme_cut_positions(e.enzyme, backbone_sequence, circular=True)
        start = site_starts[e.name]
        plan["backbone_linearization"]["restriction_enzymes"].append(
            {
                "role": e.role,
                "name": e.name,
                "recognition_site": e.site,
                "site_start_bp": start,
                "site_end_bp": start + len(e.site) - 1,
                "cut_after_bp": (start - 1 + e.enzyme.fst5) % len(backbone_sequence),
                "overhang": "five_prime" if e.enzyme.is_5overhang() else "three_prime",
            }
        )
        if not issues and len(cuts) != 1:
            issues.append(
                {
                    "code": "unexpected_cuts",
                    "message": f"{e.name} 的骨架切割结果不是一个唯一切点。",
                    "enzyme": e.name,
                }
            )
    inserted_start = final_by_id["expression"]["start_bp"]
    if issues:
        assembled = SequenceAssemblyResult(
            final_sequence,
            inserted_start,
            inserted_start + len(insert) - 1,
            replace_start,
            replace_end,
            len(left.site) + len(insert) + len(right.site),
            {
                "retained_sites": {
                    "left": {
                        "start_bp": final_by_id["left_site"]["start_bp"],
                        "end_bp": final_by_id["left_site"]["end_bp"],
                    },
                    "right": {
                        "start_bp": final_by_id["right_site"]["start_bp"],
                        "end_bp": final_by_id["right_site"]["end_bp"],
                    },
                },
                "left_enzyme": left.name,
                "right_enzyme": right.name,
                "left_site": left.site,
                "right_site": right.site,
            },
        )
    else:
        assembled = assemble_sequence(backbone_sequence, insert, plan)
    final, mapping_warnings = build_final_record(
        record_id=f"expression_D{design_id:03d}_FIN",
        backbone_record=backbone,
        insert_record=insert_record,
        result=assembled,
        plan=plan,
    )
    warnings = [w for module in modules.values() for w in module.notes]
    warnings.extend(terminator_warnings)
    warnings.extend(w for w in mapping_warnings if "landing_pad_spacer" not in w)
    preparations = {}
    if not issues:
        # Digested backbone runs from the right site around to the left site.
        # SeqRecord slicing/concatenation also remaps every retained annotation.
        core_record = backbone[replace_end:] + backbone[:replace_start]
        core = str(core_record.seq)
        try:
            preparations = {
                "backbone": _prepare_fragment(
                    right.site + core + left.site,
                    core_record,
                    len(right.site),
                    enzymes,
                    "backbone_prep",
                ),
                "insert": _prepare_fragment(
                    left.site + insert + right.site,
                    insert_record,
                    len(left.site),
                    enzymes,
                    "expression_prep",
                ),
            }
        except DesignError as exc:
            issues.append({"code": exc.code, "message": str(exc)})
    preview = {
        "component_order": order,
        "components": components,
        "t0_id": t0_id,
        "t1_id": t1_id,
        "terminator_warnings": terminator_warnings,
        "valid": not issues,
        "issues": issues,
        "warnings": warnings,
        "sequence": str(final.seq),
        "length_bp": len(final),
        "gc_percent": gc_percent(str(final.seq)),
        "sequence_sha256": sha256_sequence(str(final.seq)),
        "segments": segments,
        "features": feature_payloads(final),
        "enzymes": [e.to_payload() for e in enzymes],
    }
    return MolecularDesign(backbone, final, plan, preview, preparations)
