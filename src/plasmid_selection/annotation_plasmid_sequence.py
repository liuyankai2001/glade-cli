from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field
from src.runtime.monitor import monitor
from src.tools.common.session_paths import inputs_dir as resolve_inputs_dir
from src.tools.common.session_paths import plasmid_annotation_dir as resolve_plasmid_annotation_dir


class AnnotatePlasmidSequenceArgs(BaseModel):
    vector_id: str = Field(default="", description="可选；vector 名称或编号")
    linear: bool = Field(default=False, description="False 表示按环状质粒注释；True 表示按线性序列注释")
    db: str = Field(default="addgene", description='pLannotate 数据库："addgene"、"fpbase" 或 "snapgene"')
    identity_threshold: float = Field(default=95.0, ge=0, le=100, description="高置信 feature 的 identity 阈值")
    coverage_threshold: float = Field(default=90.0, ge=0, le=100, description="高置信 feature 的 query coverage 阈值")


VALID_BASES = set("ATGCNRYSWKMBDHV")
SEQUENCE_EXTENSIONS = {".fa", ".fasta", ".fna", ".gb", ".gbk", ".txt"}


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _default_input_dir() -> Path:
    return resolve_inputs_dir()


def _resolve_input_path(input_value: str | None, vector_id: str) -> Path:
    input_path = Path(input_value).resolve() if input_value else _default_input_dir()
    if input_path.is_file():
        return input_path
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if not input_path.is_dir():
        raise ValueError(f"input_dir must be a file or directory: {input_path}")

    candidates = [
        path
        for path in input_path.iterdir()
        if path.is_file() and path.suffix.lower() in SEQUENCE_EXTENSIONS
    ]
    if vector_id:
        matched = [
            path
            for path in candidates
            if path.stem.lower() == vector_id.lower()
        ]
        if len(matched) == 1:
            return matched[0].resolve()
        if len(matched) > 1:
            raise ValueError(f"multiple sequence files matched vector_id {vector_id}: {matched}")

    if len(candidates) != 1:
        raise ValueError(
            f"input directory must contain exactly one sequence file, found {len(candidates)}: {candidates}"
        )
    return candidates[0].resolve()


def _sanitize_name(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return "vector"
    return "".join(char if char.isalnum() or char in {"_", "-", "."} else "_" for char in text)


def _load_sequence(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix in {".gb", ".gbk"}:
        return _load_genbank(path), path.stem
    if suffix in {".fa", ".fasta", ".fna"}:
        return _load_fasta(path), path.stem
    return _load_plain_or_fasta(path), path.stem


def _load_fasta(path: Path) -> str:
    try:
        from Bio import SeqIO

        record = next(SeqIO.parse(str(path), "fasta"))
        return str(record.seq).upper()
    except ImportError:
        return _load_plain_or_fasta(path)


def _load_genbank(path: Path) -> str:
    try:
        from Bio import SeqIO

        record = next(SeqIO.parse(str(path), "genbank"))
        return str(record.seq).upper()
    except ImportError:
        return _load_genbank_origin(path)


def _load_plain_or_fasta(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    return "".join(lines).upper()


def _load_genbank_origin(path: Path) -> str:
    in_origin = False
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.upper() == "ORIGIN":
            in_origin = True
            continue
        if in_origin and stripped == "//":
            break
        if in_origin:
            chunks.append("".join(char for char in stripped if char.isalpha()))
    return "".join(chunks).upper()


def _validate_sequence(sequence: str, vector_id: str) -> dict[str, Any]:
    sequence = sequence.upper()
    invalid_bases = sorted(set(sequence) - VALID_BASES)
    if invalid_bases:
        raise ValueError(f"invalid bases in {vector_id}: {invalid_bases}")
    if len(sequence) < 100:
        raise ValueError(f"sequence too short for plasmid annotation: {len(sequence)} bp")
    gc_percent = round((sequence.count("G") + sequence.count("C")) * 100 / len(sequence), 2)
    return {
        "length_bp": len(sequence),
        "gc_percent": gc_percent,
    }


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relpath_or_abs(base_dir: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _feature_value(row: Any, key: str) -> str:
    value = row.get(key, "") if hasattr(row, "get") else ""
    if value is None:
        return ""
    return str(value)


def _feature_name(row: Any) -> str:
    for key in ("Feature", "feature", "label", "name"):
        value = _feature_value(row, key).strip()
        if value:
            return value
    return ""


def _feature_type(row: Any) -> str:
    for key in ("Feature_type", "feature_type", "type"):
        value = _feature_value(row, key).strip()
        if value:
            return value
    return ""


def _contains_any(text: str, keywords: tuple[str, ...]) -> bool:
    lowered = text.lower()
    return any(keyword in lowered for keyword in keywords)


def _summarize_detected_features(results: Any) -> dict[str, list[str]]:
    detected = {
        "origins": [],
        "resistance_markers": [],
        "promoters": [],
        "terminators": [],
        "mcs": [],
        "reporters": [],
        "tags": [],
    }
    seen = {key: set() for key in detected}

    for _, row in results.iterrows():
        name = _feature_name(row)
        feature_type = _feature_type(row)
        haystack = f"{name} {feature_type}"
        buckets = []
        if _contains_any(haystack, ("origin", "ori", "replicon")):
            buckets.append("origins")
        if _contains_any(haystack, ("ampr", "ampicillin", "kanr", "kanamycin", "cmr", "chloramphenicol", "tetr", "tetracycline", "spectinomycin", "zeocin", "hygromycin", "resistance", "resistant", "beta-lactamase", "bla")):
            buckets.append("resistance_markers")
        if _contains_any(haystack, ("promoter", "t7", "lac", "tac", "trc", "ara", "tet")):
            buckets.append("promoters")
        if _contains_any(haystack, ("terminator", "term")):
            buckets.append("terminators")
        if _contains_any(haystack, ("mcs", "multiple cloning site", "polylinker")):
            buckets.append("mcs")
        if _contains_any(haystack, ("gfp", "rfp", "cfp", "yfp", "mcherry", "fluorescent")):
            buckets.append("reporters")
        if _contains_any(haystack, ("his-tag", "flag", "ha tag", "myc", "strep", "tag")):
            buckets.append("tags")

        for bucket in buckets:
            if name and name not in seen[bucket]:
                seen[bucket].add(name)
                detected[bucket].append(name)
    return detected


def _feature_rows(results: Any) -> list[dict[str, Any]]:
    columns = [
        "Feature",
        "Feature_type",
        "start",
        "end",
        "strand",
        "pct_identity",
        "pct_query_cov",
        "database",
    ]
    rows = []
    for _, row in results.iterrows():
        item = {}
        for column in columns:
            if column in results.columns:
                value = row[column]
                if hasattr(value, "item"):
                    value = value.item()
                item[column] = value
        rows.append(item)
    return rows


def _high_confidence(results: Any, identity_threshold: float, coverage_threshold: float) -> Any:
    if "pct_identity" not in results.columns or "pct_query_cov" not in results.columns:
        return results.iloc[0:0].copy()
    return results[
        (results["pct_identity"] >= identity_threshold)
        & (results["pct_query_cov"] >= coverage_threshold)
    ].copy()


@tool(args_schema=AnnotatePlasmidSequenceArgs)
def annotate_plasmid_sequence(
    vector_id: str = "",
    linear: bool = False,
    db: str = "addgene",
    identity_threshold: float = 95.0,
    coverage_threshold: float = 90.0,
) -> str:
    """
    注释用户提供的质粒序列，生成 vector 选择可用的结构摘要。
    
    调用时机：用户上传自定义质粒，或模板库没有合适 vector。
    输入：vector_id、序列文件、topology 和阈值。
    返回：ok、序列文件摘要、检测到的 features、注释文件路径和 warnings。
    限制：不写 plasmid_selection；用户确认后调用 write_plasmid_selection。
    """

    tool_name = "annotate_plasmid_sequence"
    monitor.report_start(tool_name, {"vector_id": vector_id, "db": db})
    try:
        input_path = _resolve_input_path(None, vector_id)
        output_path = resolve_plasmid_annotation_dir()
        output_root = output_path.parent
        output_path.mkdir(parents=True, exist_ok=True)

        sequence, inferred_id = _load_sequence(input_path)
        resolved_vector_id = _sanitize_name(vector_id or inferred_id)
        sequence_info = _validate_sequence(sequence, resolved_vector_id)

        try:
            from plannotate import annotate, create_bokeh_chart, write_genbank
        except ImportError as exc:
            raise RuntimeError(
                "plannotate is not installed or not importable. Install plannotate and BLAST+ before using this tool."
            ) from exc

        monitor.report_running(tool_name, "正在运行 pLannotate 注释质粒序列...", progress=0.45)
        results = annotate(sequence, linear=linear, db=db)
        all_features_csv = output_path / f"{resolved_vector_id}_features.csv"
        high_confidence_csv = output_path / f"{resolved_vector_id}_high_confidence_features.csv"
        genbank_path = output_path / f"{resolved_vector_id}_annotated.gb"
        html_path = output_path / f"{resolved_vector_id}_map.html"
        annotation_result_json = output_path / f"{resolved_vector_id}_annotation_result.json"

        results.to_csv(all_features_csv, index=False)
        high_conf = _high_confidence(results, identity_threshold, coverage_threshold)
        high_conf.to_csv(high_confidence_csv, index=False)
        write_genbank(sequence, results, output_file=str(genbank_path))
        create_bokeh_chart(sequence, results, output_file=str(html_path))

        detected_features = _summarize_detected_features(results)
        warnings = []
        if not detected_features["origins"]:
            warnings.append("未识别到明确的复制起点 origin/ori。")
        if not detected_features["resistance_markers"]:
            warnings.append("未识别到明确的抗性标记。")
        if len(high_conf) < len(results):
            warnings.append(
                f"存在 {len(results) - len(high_conf)} 个低置信或部分匹配 feature，建议人工复核。"
            )

        result = {
            "ok": True,
            "vector_id": resolved_vector_id,
            "input_sequence_file": {
                "path": _relpath_or_abs(output_root, input_path),
                "sha256": _sha256_file(input_path),
            },
            "length_bp": sequence_info["length_bp"],
            "gc_percent": sequence_info["gc_percent"],
            "topology": "linear" if linear else "circular",
            "database": db,
            "feature_count": len(results),
            "high_confidence_feature_count": len(high_conf),
            "annotation_files": {
                "genbank": _relpath_or_abs(output_root, genbank_path),
                "html_map": _relpath_or_abs(output_root, html_path),
                "features_csv": _relpath_or_abs(output_root, all_features_csv),
                "high_confidence_features_csv": _relpath_or_abs(output_root, high_confidence_csv),
                "annotation_result_json": _relpath_or_abs(output_root, annotation_result_json),
            },
            "detected_features": detected_features,
            "features_preview": _feature_rows(results.head(20)),
            "warnings": warnings,
        }
        _write_json(annotation_result_json, result)
        monitor.report_end(tool_name, {"vector_id": resolved_vector_id, "feature_count": len(results)})
        return json.dumps(result, ensure_ascii=False, default=str)
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(annotate_plasmid_sequence.invoke({}))
