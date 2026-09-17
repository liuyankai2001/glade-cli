"""Reproduce the reviewed gap library from pinned author sources (run from repo root)."""

import ast
import csv
import hashlib
import io
import json
from copy import deepcopy
from pathlib import Path

import httpx
from Bio import SeqIO

from src.plasmid_design.gap_curation import extract_reviewed_gaps

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "src/plasmid_design/source_features.json"
SCAFFOLD = ROOT / "src/plasmid_design/scaffold.json"
VERSION = "cad60c243aa2b12c6ed52ce36be454da3f357a6e"
BASE = f"https://raw.githubusercontent.com/LondonBiofoundry/basicsynbio/{VERSION}/"


def sha(value):
    return hashlib.sha256(
        value.encode() if isinstance(value, str) else value
    ).hexdigest()


def curate():
    meta = json.loads(META.read_text(encoding="utf-8"))
    scaffold = json.loads(SCAFFOLD.read_text(encoding="utf-8"))
    archive = ROOT / "data/gap_sources"
    archive.mkdir(exist_ok=True)
    contents = {}
    additions = {
        "standard_linkers": "basicsynbio/parts_linkers/biolegio_linkers.py",
        "pkd46_core": "sequences/genbank_files/BASIC_SEVA_collection/misc_seqs/BS_x7x_pKD46_core_seq.gb",
    }
    for key, path in additions.items():
        meta["_sources"].setdefault(
            key, {"source_url": BASE + path, "source_version": VERSION}
        )
    for key, source in meta["_sources"].items():
        target = archive / (key + (".py.txt" if key == "standard_linkers" else ".gb"))
        if target.exists():
            raw = target.read_bytes()
        else:
            for attempt in range(3):
                try:
                    raw = (
                        httpx.get(source["source_url"], timeout=30)
                        .raise_for_status()
                        .content
                    )
                    break
                except httpx.TransportError:
                    if attempt == 2:
                        raise
        if source.get("source_sha256") and sha(raw) != source["source_sha256"]:
            raise ValueError(f"pinned source hash mismatch: {key}")
        source["source_sha256"] = sha(raw)
        target.write_bytes(raw)
        contents[key] = raw
    records = {
        key: {r.name: r for r in SeqIO.parse(io.StringIO(raw.decode()), "genbank")}
        for key, raw in contents.items()
        if key != "standard_linkers"
    }
    gaps = {}

    def provenance(key, record, start, end):
        return {
            **meta["_sources"][key],
            "source_key": key,
            "source_record": record,
            "start_1based": start,
            "end_1based": end,
        }

    def add(sequence, name, alias, purpose, source, evidence, note):
        identifier = "gap_" + sha(sequence)
        if not sequence or set(sequence) - set("ACGT"):
            raise ValueError("gap must be unambiguous retained DNA")
        item = gaps.setdefault(
            identifier,
            {
                "id": identifier,
                "name": name,
                "sequence": sequence,
                "type": "gap",
                "sha256": sha(sequence),
                "aliases": [],
                "purpose": purpose,
                "evidence_status": "source_sequence_verified",
                "notes": [],
                "provenance": source,
                "sources": [],
                "features": [],
            },
        )
        if alias not in item["aliases"]:
            item["aliases"].append(alias)
        if note not in item["notes"]:
            item["notes"].append(note)
        occurrence = {**source, "evidence": evidence}
        if occurrence not in item["sources"]:
            item["sources"].append(occurrence)
        return identifier

    tree = ast.parse(contents["standard_linkers"].decode())
    values = next(
        ast.literal_eval(n.value)
        for n in tree.body
        if isinstance(n, ast.Assign)
        and any(
            isinstance(t, ast.Name) and t.id == "STANDARD_LINKERS" for t in n.targets
        )
    )
    for alias in [f"L{i}" for i in range(1, 7)]:
        sequence = values[alias]
        add(
            sequence,
            f"BASIC {alias}",
            alias,
            "通用中性 linker 候选",
            provenance(
                "standard_linkers", f"STANDARD_LINKERS.{alias}", 1, len(sequence)
            ),
            "Author STANDARD_LINKERS literal; exclude the prepended GG oligo adapter.",
            "来源定义的 BASIC linker；收录保留 DNA，不含制备用 GG 接头。新位置组合尚未实验验证。",
        )
    record = records["collection"]["BASIC_SEVA_15a.10"]
    annotated = next(
        f for f in record.features if f.qualifiers.get("label") == ["BSEVA_L1"]
    )
    add(
        str(annotated.extract(record.seq)).upper(),
        "BSEVA L1",
        "BSEVA_L1",
        "抗性与复制模块之间的中性 linker 候选",
        provenance(
            "collection",
            record.name,
            int(annotated.location.start) + 1,
            int(annotated.location.end),
        ),
        "Author BSEVA_L1 annotation; BASIC SEVA paper DOI 10.1093/synbio/ysac023.",
        "作者为抗性与复制模块间设计的中性 linker；其他位置组合尚未实验验证。",
    )
    native_ids = {}
    purposes = {
        "resistance_to_replication": "来源原生模块接口（76 bp）",
        "replication_to_t1": "来源原生复制模块末端接口（14 bp）",
        "landing_pad_spacer": "临时骨架占位片段（24 bp）",
    }
    names = {
        "resistance_to_replication": "SEVA 模块间隔 76",
        "replication_to_t1": "SEVA 复制接口 14",
        "landing_pad_spacer": "骨架占位 24",
    }
    for alias in purposes:
        detail = scaffold["provenance"][alias]
        source = provenance(
            "collection", record.name, detail["start_1based"], detail["end_1based"]
        )
        sequence = str(
            record.seq[detail["start_1based"] - 1 : detail["end_1based"]]
        ).upper()
        assert sha(sequence) == detail["sha256"]
        native_ids[alias] = add(
            sequence,
            names[alias],
            alias,
            purposes[alias],
            source,
            "Pinned native scaffold coordinates, exact sequence hash verified.",
            (
                "来源原生接口，不表示在任意上下文中都具有生物学中性。"
                if alias != "landing_pad_spacer"
                else "仅用于临时骨架占位；最终质粒只在用户明确选入时包含该片段。"
            ),
        )
    audit = []
    for typ in ("resistance", "replication", "terminator"):
        path = ROOT / f"data/{typ}_modules.csv"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            columns, rows = reader.fieldnames, list(reader)
        before_rows = deepcopy(rows)
        for row in rows:
            item = meta[row["id"]]
            prior = item.get("gap_extraction")
            original = prior["original_sequence"] if prior else row["sequence"]
            features = deepcopy(
                prior["original_features"] if prior else item["features"]
            )
            p = item["provenance"]
            source_record = records[p["source_key"]][p["source_record"]]
            assert (
                str(source_record.seq[p["start_1based"] - 1 : p["end_1based"]]).upper()
                == original
            )
            spans, pending = [], []
            if typ == "resistance":
                prefix_len = (
                    34 if row["id"] in ("basic_seva_sm_sp", "basic_seva_tet_5a") else 26
                )
                assert (
                    min(part["start"] for f in features for part in f["parts"])
                    == prefix_len
                )
                evidence = "Pinned source whole resistance cassette annotation starts after this native boundary interval."
                span = {
                    "start_1based": 1,
                    "end_1based": prefix_len,
                    "evidence": evidence,
                }
                source = {**p, "end_1based": p["start_1based"] + prefix_len - 1}
                span["gap_id"] = add(
                    original[:prefix_len],
                    f"SEVA 抗性接口 {prefix_len}",
                    f"resistance_prefix_{prefix_len}",
                    "来源原生抗性模块前部接口",
                    source,
                    evidence,
                    "位于作者定义的完整抗性功能模块之外；保留来源接口用途，不作为通用中性 spacer。",
                )
                spans.append(span)
                pending.append(
                    {
                        "start_1based": prefix_len + 1,
                        "end_1based": len(original),
                        "reason": "Source defines this entire resistance cassette as functional; uncovered promoter/terminator/regulatory DNA is retained, not inferred to be gap.",
                    }
                )
            elif row["id"] == "basic_seva_psc101_pkd46_ts":
                core = next(iter(records["pkd46_core"].values()))
                assert str(core.seq[: len(original)]).upper() == original
                assert any(
                    int(f.location.start) == 1551
                    for f in core.features
                    if f.qualifiers.get("label") == ["SEVA_T1"]
                )
                assert (
                    original[-14:] == gaps[native_ids["replication_to_t1"]]["sequence"]
                )
                evidence = "Identical 14 bp native interface immediately preceding author-annotated SEVA_T1 in both collection and independent pKD46 core."
                source = {**p, "start_1based": p["end_1based"] - 13}
                identifier = add(
                    original[-14:],
                    names["replication_to_t1"],
                    "pkd46_tail_14",
                    purposes["replication_to_t1"],
                    source,
                    evidence,
                    "亦从 pKD46 来源记录尾部提取；复制功能注释与其余调控区保留。",
                )
                add(
                    original[-14:],
                    names["replication_to_t1"],
                    "pkd46_tail_14",
                    purposes["replication_to_t1"],
                    provenance("pkd46_core", core.name, 1538, 1551),
                    evidence,
                    "pKD46 独立核心来源的同一接口。",
                )
                spans.append(
                    {
                        "start_1based": 1538,
                        "end_1based": 1551,
                        "gap_id": identifier,
                        "evidence": evidence,
                    }
                )
                pending = [
                    {
                        "start_1based": a,
                        "end_1based": b,
                        "reason": "Unannotated DNA may regulate replication; no evidence it is a removable gap.",
                    }
                    for a, b in ((1, 61), (1013, 1059), (1283, 1537))
                ]
            clean, mapped, retained = extract_reviewed_gaps(original, features, spans)
            for segment in retained:
                segment.update(
                    source_start_1based=p["start_1based"]
                    + segment["original_start_1based"]
                    - 1,
                    source_end_1based=p["start_1based"]
                    + segment["original_end_1based"]
                    - 1,
                )
            item.update(
                sha256=sha(clean),
                features=mapped,
                gap_extraction={
                    "status": (
                        "confirmed_removed" if spans else "no_confirmed_internal_gap"
                    ),
                    "original_sequence": original,
                    "original_sha256": sha(original),
                    "original_features": features,
                    "removed_gaps": spans,
                    "retained_segments": retained,
                    "pending_regions": pending,
                },
            )
            row["sequence"] = clean
            audit.append(
                f"| {row['id']} | {len(original)} | {len(clean)} | {', '.join(str(s['end_1based']-s['start_1based']+1) for s in spans) or '0'} | {len(pending)} |"
            )
        if rows != before_rows:
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    for identifier, item in gaps.items():
        meta[identifier] = {
            k: v for k, v in item.items() if k not in ("id", "name", "sequence")
        }
    meta["_schema_version"] = 3
    META.write_text(
        json.dumps(meta, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    with (ROOT / "data/gap_modules.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["id", "name", "sequence"], lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            {k: item[k] for k in ("id", "name", "sequence")} for item in gaps.values()
        )
    SCAFFOLD.write_text(
        json.dumps(
            {
                "_schema_version": 3,
                "version": "basic_seva_gap_references.v2",
                "gap_ids": native_ids,
                "provenance": scaffold["provenance"],
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    (ROOT / "docs/gap-extraction-audit.md").write_text(
        "# Gap 来源及提取审计\n\n范围：当前 7 个抗性、5 个复制、2 个终止子和固定接口；不扫描完整表达构建。\n\n"
        "只拆出有来源依据且不覆盖功能注释的边界片段。未确认的非编码/调控区保留在模块中。所有原始 DNA、哈希、注释、删除范围、保留坐标映射及待核查区域存于 source_features.json。\n\n"
        "| 模块 | 原始 bp | 当前 bp | 已拆出 bp | 保留待核查区域数 |\n| --- | ---: | ---: | --- | ---: |\n"
        + "\n".join(audit)
        + "\n\n"
        "抗性前部26/34 bp依据作者的完整功能模块注释边界；pKD46尾部14 bp在集合及独立核心来源中紧邻SEVA_T1，且与原生复制接口完全一致。Ap末端含终止作用的DNA仍属于功能模块，未拆除。其他复制模块及T0/T1没有可确认的内部gap。pKD46局部1–61、1013–1059、1283–1537 bp可能参与调控，保留待核查。\n\n"
        f"归库 {len(gaps)} 条完整DNA去重序列：BASIC L1–L6、BSEVA_L1、原生76/14 bp接口、24 bp占位及抗性26/34 bp接口。53 bp linker不含制备用GG接头；占位片段不自动加入最终质粒。相同DNA的全部本次提取来源保存在sources数组。\n\n"
        "来源：[作者固定版本仓库](https://github.com/LondonBiofoundry/basicsynbio/tree/"
        + VERSION
        + ")；[BASIC SEVA论文](https://doi.org/10.1093/synbio/ysac023)。原始文件保存在data/gap_sources，逐文件SHA256校验；重现：`python -B -m data.curate_gap_modules`。来源定义的中性linker是候选条目，不表示新位置组合已经实验验证。\n",
        encoding="utf-8",
        newline="\n",
    )
    (ROOT / "data/gap_modules.README.md").write_text(
        "# Gap 序列库\n\nCSV字段严格为id,name,sequence；DNA为大写ACGT。ID=gap_+完整序列SHA256；完全相同DNA合并，保留aliases及所有sources。\n\n"
        "当前12条：6个BASIC标准linker、BSEVA_L1、76/14bp原生接口、24bp占位、26/34bp抗性接口。用途、证据及注释存于src/plasmid_design/source_features.json；模块原始快照及删除映射也存于该文件。详见docs/gap-extraction-audit.md。\n\n"
        "新装配schema2只包含明确选入的gap实例；旧装配仅将原76/14bp外部间隔显式迁移，不恢复已经拆出的内部边界。可重复选入同一条DNA。增加数据时需同步固定来源、特征快照及用途；不要把任意非编码DNA推断为gap。\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"Curated {len(gaps)} gap sequences; audited {len(audit)} modules.")


if __name__ == "__main__":
    curate()
