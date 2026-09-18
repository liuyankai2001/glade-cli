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
    # Native boundary removals remain audit evidence, not selectable Gap modules.
    native_replication_tail = str(
        records["collection"]["BASIC_SEVA_15a.10"].seq[4978:4992]
    ).upper()
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
                span["gap_id"] = "gap_" + sha(original[:prefix_len])
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
                assert original[-14:] == native_replication_tail
                evidence = "Identical 14 bp native interface immediately preceding author-annotated SEVA_T1 in both collection and independent pKD46 core."
                assert str(core.seq[1537:1551]).upper() == native_replication_tail
                identifier = "gap_" + sha(original[-14:])
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
    for identifier in list(meta):
        if isinstance(meta[identifier], dict) and meta[identifier].get("type") == "gap":
            del meta[identifier]
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
                "version": "basic_linker_landing_pad.v1",
                "gap_ids": {"landing_pad_spacer": "gap_" + sha(values["L1"])},
                "provenance": {
                    **provenance(
                        "standard_linkers", "STANDARD_LINKERS.L1", 1, len(values["L1"])
                    ),
                    "landing_pad_spacer": {
                        "start_1based": 1,
                        "end_1based": len(values["L1"]),
                        "sha256": sha(values["L1"]),
                    },
                },
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    (ROOT / "docs/gap-extraction-audit.md").write_text(
        "# Gap 来源及提取审计\n\n范围：当前 7 个抗性、5 个复制、2 个终止子；不扫描完整表达构建。\n\n"
        "只拆出有来源依据且不覆盖功能注释的边界片段。未确认的非编码/调控区保留在模块中。所有原始 DNA、哈希、注释、删除范围、保留坐标映射及待核查区域存于 source_features.json。\n\n"
        "| 模块 | 原始 bp | 当前 bp | 已拆出 bp | 保留待核查区域数 |\n| --- | ---: | ---: | --- | ---: |\n"
        + "\n".join(audit)
        + "\n\n"
        "抗性前部26/34 bp依据作者的完整功能模块注释边界；pKD46尾部14 bp在集合及独立核心来源中紧邻SEVA_T1，且与原生复制接口完全一致。Ap末端含终止作用的DNA仍属于功能模块，未拆除。其他复制模块及T0/T1没有可确认的内部gap。pKD46局部1–61、1013–1059、1283–1537 bp可能参与调控，保留待核查。\n\n"
        f"可选 Gap 库仅保留 {len(gaps)} 条 BASIC L1–L6，均为53 bp，不含制备用GG接头。BSEVA_L1、76/14 bp接口、24 bp占位及26/34 bp抗性接口均不再归库。已拆出的边界以内容哈希留存审计，校验不依赖可选 Gap 库。临时骨架占位改用BASIC L1，生成最终质粒时由完整表达构建替换。新旧请求均只插入用户明确选择的Gap。\n\n"
        "来源：[作者固定版本仓库](https://github.com/LondonBiofoundry/basicsynbio/tree/"
        + VERSION
        + ")；[BASIC SEVA论文](https://doi.org/10.1093/synbio/ysac023)。原始文件保存在data/gap_sources，逐文件SHA256校验；重现：`python -B -m data.curate_gap_modules`。来源定义的中性linker是候选条目，不表示新位置组合已经实验验证。\n",
        encoding="utf-8",
        newline="\n",
    )
    (ROOT / "data/gap_modules.README.md").write_text(
        "# Gap 序列库\n\nCSV字段严格为id,name,sequence；DNA为大写ACGT。ID=gap_+完整序列SHA256；完全相同DNA合并，保留aliases及所有sources。\n\n"
        "当前仅6条：BASIC L1–L6，均为53 bp通用中性linker候选，不含制备用GG接头。用途、证据及注释存于src/plasmid_design/source_features.json；模块原始快照及删除映射也存于该文件。已拆出的原生边界仅留存来源审计，不再作为可选Gap。详见docs/gap-extraction-audit.md。\n\n"
        "新旧请求均只包含明确选入的gap实例，不自动补入旧76/14 bp间隔。临时骨架占位使用BASIC L1，最终由完整表达构建替换，不自动残留。纯序列拼接可重复选入同一条DNA；实际BASIC装配一轮内应使用不同linker。候选用于独立功能模块边界，不保证任意组合均具有生物学中性。重新运行归库脚本仍只生成L1–L6。\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"Curated {len(gaps)} gap sequences; audited {len(audit)} modules.")


if __name__ == "__main__":
    curate()
