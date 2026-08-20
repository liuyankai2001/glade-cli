from __future__ import annotations
import json
from pathlib import Path
from typing import Any

from langchain.tools import tool
from pydantic import BaseModel, Field, field_validator

from src.runtime.monitor import monitor
from src.tools.common.session_paths import (
    design_manifest_file,
    outputs_dir as resolve_outputs_dir,
    session_dir as resolve_session_dir,
)


class RecommendGibsonInsertionSitesArgs(BaseModel):
    cassette_order: list[int] | None = Field(
        default=None,
        description="表达盒顺序；为空时按 cassette_index 升序拼接 insert",
    )
    linker_sequence: str = Field(default="", description="多个 cassette 之间的连接序列，可为空")
    homology_arm_length: int = Field(default=30, ge=15, le=80, description="Gibson 同源臂长度")
    max_recommendations: int = Field(default=3, ge=1, le=10, description="最多返回多少个候选插入位点")
    preferred_feature_keywords: list[str] | None = Field(
        default=None,
        description="可选；优先匹配的 feature 关键词，例如 ['MCS', 'lacZ']",
    )

    @field_validator("cassette_order", mode="before")
    @classmethod
    def normalize_cassette_order(cls, value: Any) -> list[int] | None:
        return _normalize_cassette_order(value)

    @field_validator("preferred_feature_keywords", mode="before")
    @classmethod
    def normalize_preferred_feature_keywords(cls, value: Any) -> list[str] | None:
        return _normalize_string_list_parameter(value)


PROTECTED_FEATURE_TYPES = {
    "CDS",
    "promoter",
    "terminator",
    "rep_origin",
    "origin_of_replication",
    "oriT",
}


def _json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def _session_root(
    session_dir: str | None,
    output_dir: str | None,
    manifest_path: str | None,
) -> Path:
    if session_dir:
        return Path(session_dir).resolve()

    if output_dir:
        output_path = Path(output_dir).resolve()
        return output_path.parent.resolve() if output_path.name == "outputs" else output_path.parent.resolve()

    if manifest_path:
        path = Path(manifest_path).resolve()
        if path.name == "design_manifest.json" and path.parent.name == "outputs":
            return path.parent.parent.resolve()

    return resolve_session_dir()


def _output_dir(session_root: Path, output_dir: str | None) -> Path:
    if output_dir:
        return Path(output_dir).resolve()

    return resolve_outputs_dir()


def _manifest_file(session_root: Path, output_dir: Path, manifest_path: str | None) -> Path:
    if manifest_path:
        path = Path(manifest_path)
        if not path.is_absolute():
            output_relative = output_dir / path
            session_relative = session_root / path
            path = output_relative if output_relative.exists() else session_relative
        return path.resolve()

    return design_manifest_file()


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"design_manifest.json root must be a JSON object: {path}")
    return data


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _normalize_cassette_order(value: Any) -> list[int] | None:
    if value is None or value == "":
        return None
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            value = json.loads(text)
        except json.JSONDecodeError:
            value = [
                part.strip()
                for part in text.replace("|", ",").replace(";", ",").split(",")
                if part.strip()
            ]
    if not isinstance(value, list):
        raise ValueError("cassette_order must be a list of positive cassette indexes, for example [1, 2, 3]")

    order = []
    for item in value:
        if isinstance(item, bool):
            raise ValueError("cassette_order values must be positive integers, not booleans")
        try:
            index = int(item)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"cassette_order contains a non-integer value: {item!r}") from exc
        if index <= 0:
            raise ValueError("cassette_order values must be positive cassette indexes")
        order.append(index)
    return order or None


def _normalize_string_list_parameter(
    value: Any,
) -> list[str] | None:
    if value is None or value == "":
        return None
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            value = json.loads(text)
        except json.JSONDecodeError:
            value = [
                part.strip()
                for part in text.replace("|", ",").replace(";", ",").split(",")
                if part.strip()
            ]
    if not isinstance(value, list):
        raise ValueError("expected a list of non-empty strings")

    result = []
    seen = set()
    for item in value:
        if not isinstance(item, str):
            raise ValueError("string list values must be strings")
        text = item.strip()
        if not text:
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        result.append(text)
    return result or None


def _relpath_or_abs(base_dir: Path, path_value: Any) -> str:
    if not path_value:
        return ""
    path = Path(str(path_value))
    if not path.is_absolute():
        path = base_dir / path
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _resolve_file_path(output_dir: Path, path_value: Any) -> Path:
    if not path_value:
        raise ValueError("sequence file path is empty")

    path = Path(str(path_value))
    if path.is_absolute():
        return path.resolve()

    output_relative = (output_dir / path).resolve()
    if output_relative.exists():
        return output_relative

    session_relative = (output_dir.parent / path).resolve()
    if session_relative.exists():
        return session_relative

    return output_relative


def _clean_dna(value: Any) -> str:
    return "".join(base for base in str(value or "").upper() if base in "ACGTN")


def _read_sequence(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".gb", ".gbk", ".gbj"}:
        sequence = _read_genbank_sequence(path)
    else:
        sequence = _read_plain_or_fasta_sequence(path)

    sequence = _clean_dna(sequence)
    if not sequence:
        raise ValueError(f"sequence file is empty or unsupported: {path}")
    return sequence


def _read_plain_or_fasta_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    return "".join(
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    )


def _read_genbank_sequence(path: Path) -> str:
    try:
        from Bio import SeqIO

        record = next(SeqIO.parse(str(path), "genbank"))
        return str(record.seq)
    except Exception:
        return _read_genbank_origin(path)


def _read_genbank_origin(path: Path) -> str:
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
    return "".join(chunks)


def _feature_labels(feature: Any) -> list[str]:
    labels = []
    for key in ("label", "gene", "product", "note", "bound_moiety"):
        labels.extend(str(value) for value in feature.qualifiers.get(key, []))
    return [label for label in labels if label]


def _genbank_features(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() not in {".gb", ".gbk", ".gbj"}:
        return []

    try:
        from Bio import SeqIO
    except Exception:
        return []

    try:
        record = next(SeqIO.parse(str(path), "genbank"))
    except Exception:
        return []

    features = []
    for feature in record.features:
        labels = _feature_labels(feature)
        start_bp = int(feature.location.start) + 1
        end_bp = int(feature.location.end)
        features.append({
            "type": str(feature.type),
            "label": labels[0] if labels else str(feature.type),
            "labels": labels,
            "start_bp": start_bp,
            "end_bp": end_bp,
            "length_bp": end_bp - start_bp + 1,
            "haystack": f"{feature.type} {' '.join(labels)}".lower(),
        })
    return features


def _ordered_cassettes(manifest: dict[str, Any], cassette_order: list[int] | None) -> list[dict[str, Any]]:
    cassette_files = [
        item
        for item in _as_list(_as_dict(manifest.get("assembled_expression_cassettes")).get("cassette_files"))
        if isinstance(item, dict)
    ]
    cassette_by_index = {
        int(item.get("cassette_index")): item
        for item in cassette_files
    }
    order = cassette_order or sorted(cassette_by_index)

    missing = [int(index) for index in order if int(index) not in cassette_by_index]
    if missing:
        raise ValueError(f"cassette_order contains unknown cassette_index values: {missing}")

    return [cassette_by_index[int(index)] for index in order]


def _build_insert_sequence(
    output_dir: Path,
    cassettes: list[dict[str, Any]],
    linker_sequence: str,
) -> tuple[str, list[dict[str, Any]], str]:
    sequences = []
    files = []
    linker = _clean_dna(linker_sequence)

    for index, cassette in enumerate(cassettes):
        path = _resolve_file_path(output_dir, cassette.get("path"))
        sequence = _read_sequence(path)
        if index > 0 and linker:
            sequences.append(linker)
        sequences.append(sequence)
        files.append({
            "cassette_index": cassette.get("cassette_index"),
            "path": _relpath_or_abs(output_dir, path),
            "length_bp": len(sequence),
            "sha256": cassette.get("sha256", ""),
        })

    return "".join(sequences), files, linker


def _candidate_score(feature: dict[str, Any], preferred_keywords: list[str] | None) -> tuple[int, str]:
    haystack = str(feature.get("haystack") or "")
    score = 0
    category = ""

    if any(keyword in haystack for keyword in ("mcs", "multiple cloning", "polylinker")):
        score = 100
        category = "mcs"
    elif any(keyword in haystack for keyword in ("cloning site", "insertion site")):
        score = 90
        category = "cloning_site"
    elif any(keyword in haystack for keyword in ("insert", "stuffer")):
        score = 70
        category = "insert_like"
    elif "lacz" in haystack:
        score = 50
        category = "lacz_region"

    for keyword in preferred_keywords or []:
        cleaned = str(keyword or "").strip().lower()
        if cleaned and cleaned in haystack:
            score += 20
            category = category or "preferred_feature"

    return score, category


def _candidate_regions(features: list[dict[str, Any]], preferred_keywords: list[str] | None) -> list[dict[str, Any]]:
    regions = []
    seen = set()
    for feature in features:
        if feature.get("type") == "source":
            continue

        score, category = _candidate_score(feature, preferred_keywords)
        if score <= 0:
            continue

        key = (int(feature["start_bp"]), int(feature["end_bp"]), str(feature.get("label") or ""))
        if key in seen:
            continue
        seen.add(key)

        regions.append({
            "label": feature.get("label", ""),
            "type": feature.get("type", ""),
            "category": category,
            "start_bp": int(feature["start_bp"]),
            "end_bp": int(feature["end_bp"]),
            "length_bp": int(feature["length_bp"]),
            "score": score,
        })

    regions.sort(key=lambda item: (-int(item["score"]), int(item["length_bp"]), int(item["start_bp"])))
    return regions


def _audited_insertion_regions(vector: dict[str, Any]) -> list[dict[str, Any]]:
    regions = []
    for item in _as_list(vector.get("insertion_regions")):
        if not isinstance(item, dict):
            continue
        start_bp = int(item.get("start_bp") or 0)
        end_bp = int(item.get("end_bp") or 0)
        if start_bp <= 0 or end_bp < start_bp:
            continue
        is_seva_cargo = (
            vector.get("assembly_policy") == "replace_seva_cargo_paci_spei"
            and item.get("type") == "seva_cargo"
        )
        regions.append(
            {
                "label": item.get("name") or item.get("type") or "audited insertion region",
                "type": item.get("type") or "audited_insertion_region",
                "category": "audited_seva_cargo" if is_seva_cargo else "audited_insertion_region",
                "start_bp": start_bp,
                "end_bp": end_bp,
                "length_bp": end_bp - start_bp + 1,
                "score": 1000 if is_seva_cargo else 500,
                "audit_source": "plasmid_selection.vector.insertion_regions",
            }
        )
    return regions


def _extract_segment(
    sequence: str,
    start_bp: int,
    length: int,
    *,
    circular: bool,
) -> dict[str, Any] | None:
    if length <= 0:
        raise ValueError("homology_arm_length must be positive")

    sequence_length = len(sequence)
    end_bp = start_bp + length - 1

    if not circular and (start_bp < 1 or end_bp > sequence_length):
        return None

    if circular:
        bases = [
            sequence[(start_bp - 1 + offset) % sequence_length]
            for offset in range(length)
        ]
        normalized_start = ((start_bp - 1) % sequence_length) + 1
        normalized_end = ((end_bp - 1) % sequence_length) + 1
        wraps_origin = start_bp < 1 or end_bp > sequence_length or normalized_start > normalized_end
    else:
        bases = list(sequence[start_bp - 1:end_bp])
        normalized_start = start_bp
        normalized_end = end_bp
        wraps_origin = False

    return {
        "sequence": "".join(bases),
        "start_bp": normalized_start,
        "end_bp": normalized_end,
        "wraps_origin": wraps_origin,
    }


def _count_occurrences(sequence: str, query: str, *, circular: bool = False) -> int:
    if not sequence or not query or len(query) > len(sequence):
        return 0

    search_space = sequence
    if circular:
        search_space = sequence + sequence[: len(query) - 1]

    count = 0
    for index in range(0, len(search_space) - len(query) + 1):
        if circular and index >= len(sequence):
            break
        if search_space[index:index + len(query)] == query:
            count += 1
    return count


def _overlaps(left: dict[str, Any], right: dict[str, Any]) -> bool:
    return int(left["start_bp"]) <= int(right["end_bp"]) and int(right["start_bp"]) <= int(left["end_bp"])


def _protected_overlaps(region: dict[str, Any], features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    overlaps = []
    region_bounds = {
        "start_bp": int(region["start_bp"]),
        "end_bp": int(region["end_bp"]),
    }
    region_label = str(region.get("label") or "").lower()

    for feature in features:
        feature_type = str(feature.get("type") or "")
        feature_label = str(feature.get("label") or "")
        haystack = str(feature.get("haystack") or "")

        if feature_type in {"source", "primer_bind"}:
            continue
        if feature_label.lower() == region_label and _overlaps(region_bounds, feature):
            continue
        if any(keyword in haystack for keyword in ("mcs", "multiple cloning", "polylinker")):
            continue

        protected = feature_type in PROTECTED_FEATURE_TYPES or any(
            keyword in haystack
            for keyword in ("resistance", "origin", "ori", "promoter", "terminator")
        )
        if protected and _overlaps(region_bounds, feature):
            overlaps.append({
                "type": feature_type,
                "label": feature_label,
                "start_bp": feature.get("start_bp"),
                "end_bp": feature.get("end_bp"),
            })

    return overlaps


def _recommendation_for_region(
    *,
    region: dict[str, Any],
    features: list[dict[str, Any]],
    vector_sequence: str,
    vector_is_circular: bool,
    insert_sequence: str,
    homology_arm_length: int,
    cassette_order: list[int],
    linker_sequence: str,
) -> dict[str, Any] | None:
    start_bp = int(region["start_bp"])
    end_bp = int(region["end_bp"])
    left_arm = _extract_segment(
        vector_sequence,
        start_bp - homology_arm_length,
        homology_arm_length,
        circular=vector_is_circular,
    )
    right_arm = _extract_segment(
        vector_sequence,
        end_bp + 1,
        homology_arm_length,
        circular=vector_is_circular,
    )

    if left_arm is None or right_arm is None:
        return None

    left_homology = str(left_arm["sequence"])
    right_homology = str(right_arm["sequence"])
    overlaps = _protected_overlaps(region, features)
    if region.get("category") == "audited_seva_cargo":
        overlaps = [
            overlap
            for overlap in overlaps
            if not any(
                token in str(overlap.get("label") or "").lower()
                for token in ("lacz", "lac promoter", "cloning cargo")
            )
        ]

    warnings = []
    if overlaps:
        warnings.append("候选替换区间与关键 feature 重叠，执行前需要人工复核。")
    if _count_occurrences(vector_sequence, left_homology, circular=vector_is_circular) != 1:
        warnings.append("左侧同源臂在 vector 中不是唯一序列。")
    if _count_occurrences(vector_sequence, right_homology, circular=vector_is_circular) != 1:
        warnings.append("右侧同源臂在 vector 中不是唯一序列。")
    if _count_occurrences(insert_sequence, left_homology) > 0:
        warnings.append("左侧同源臂序列也出现在 insert 中，可能造成非预期重组。")
    if _count_occurrences(insert_sequence, right_homology) > 0:
        warnings.append("右侧同源臂序列也出现在 insert 中，可能造成非预期重组。")

    rationale = [
        f"候选区间来自 vector feature: {region.get('label') or region.get('type')}",
        "推荐用表达盒 insert 替换该区间，并使用区间两侧序列作为 Gibson 同源臂。",
    ]
    if region.get("category") == "mcs":
        rationale.append("该区间被识别为 MCS/multiple cloning site，适合作为插入位点。")
    if region.get("category") == "audited_seva_cargo":
        rationale.append("该区间是载体审计字段指定的 PacI–SpeI cargo，必须整体替换。")

    plan_args = {
        "assembly_method": "gibson",
        "cassette_order": cassette_order,
        "linker_sequence": linker_sequence,
        "insert_after_bp": None,
        "replace_start_bp": start_bp,
        "replace_end_bp": end_bp,
        "homology_arm_length": homology_arm_length,
        "left_homology": left_homology,
        "right_homology": right_homology,
    }

    return {
        "label": region.get("label", ""),
        "feature_type": region.get("type", ""),
        "mode": "replace",
        "replace_start_bp": start_bp,
        "replace_end_bp": end_bp,
        "replaced_length_bp": end_bp - start_bp + 1,
        "insert_after_bp": None,
        "homology_arm_length": homology_arm_length,
        "left_homology": left_homology,
        "right_homology": right_homology,
        "left_homology_region": {
            "start_bp": left_arm["start_bp"],
            "end_bp": left_arm["end_bp"],
            "wraps_origin": left_arm["wraps_origin"],
        },
        "right_homology_region": {
            "start_bp": right_arm["start_bp"],
            "end_bp": right_arm["end_bp"],
            "wraps_origin": right_arm["wraps_origin"],
        },
        "estimated_final_length_bp": len(vector_sequence) - (end_bp - start_bp + 1) + len(insert_sequence),
        "protected_feature_overlaps": overlaps,
        "warnings": warnings,
        "rationale": rationale,
        "plan_final_assembly_args": plan_args,
    }


@tool(args_schema=RecommendGibsonInsertionSitesArgs)
def recommend_gibson_insertion_sites(
    cassette_order: list[int] | str | None = None,
    linker_sequence: str = "",
    homology_arm_length: int = 30,
    max_recommendations: int = 3,
    preferred_feature_keywords: list[str] | str | None = None,
) -> str:
    """
    为当前 vector 和 assembled insert 推荐 Gibson 插入/替换位点及同源臂。
    
    调用时机：用户选择 Gibson，但缺少插入坐标、替换区间或 homology arms。
    返回：ok、vector/insert 摘要、candidate_regions、recommendations、plan_final_assembly_args 和 warnings。
    限制：只读；不写 manifest；用户确认候选后才调用 plan_final_assembly。
    """

    tool_name = "recommend_gibson_insertion_sites"
    monitor.report_start(tool_name, {"homology_arm_length": homology_arm_length, "max_recommendations": max_recommendations})
    try:
        cassette_order = _normalize_cassette_order(cassette_order)
        preferred_feature_keywords = _normalize_string_list_parameter(preferred_feature_keywords)

        session_root = _session_root(None, None, None)
        outputs = _output_dir(session_root, None)
        manifest_file = _manifest_file(session_root, outputs, None)
        manifest = _read_json(manifest_file)

        vector = _as_dict(_as_dict(manifest.get("plasmid_selection")).get("vector"))
        vector_sequence_file = _as_dict(vector.get("sequence_file"))
        vector_path = _resolve_file_path(outputs, vector_sequence_file.get("path"))
        vector_sequence = _read_sequence(vector_path)
        vector_is_circular = str(vector.get("topology") or "").strip().lower() == "circular"
        if len(vector_sequence) < homology_arm_length * 2:
            raise ValueError("vector sequence is shorter than two homology arms")

        cassettes = _ordered_cassettes(manifest, cassette_order)
        resolved_order = [int(item.get("cassette_index")) for item in cassettes]
        insert_sequence, cassette_files, cleaned_linker = _build_insert_sequence(outputs, cassettes, linker_sequence)

        features = _genbank_features(vector_path)
        audited_regions = _audited_insertion_regions(vector)
        if vector.get("assembly_policy") == "replace_seva_cargo_paci_spei":
            candidate_regions = [
                region
                for region in audited_regions
                if region.get("category") == "audited_seva_cargo"
            ]
            if len(candidate_regions) != 1:
                raise ValueError(
                    "selected SEVA lacZ-alpha vector lacks one audited PacI-SpeI cargo region"
                )
        else:
            candidate_regions = audited_regions + _candidate_regions(
                features,
                preferred_feature_keywords,
            )

        recommendations = []
        for region in candidate_regions:
            recommendation = _recommendation_for_region(
                region=region,
                features=features,
                vector_sequence=vector_sequence,
                vector_is_circular=vector_is_circular,
                insert_sequence=insert_sequence,
                homology_arm_length=homology_arm_length,
                cassette_order=resolved_order,
                linker_sequence=cleaned_linker,
            )
            if recommendation is not None:
                recommendations.append(recommendation)
            if len(recommendations) >= max_recommendations:
                break

        for rank, recommendation in enumerate(recommendations, start=1):
            recommendation["rank"] = rank

        monitor.report_end(tool_name, {"recommendation_count": len(recommendations)})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": _relpath_or_abs(session_root, manifest_file),
                "vector": {
                    "name": vector.get("name", ""),
                    "path": _relpath_or_abs(outputs, vector_path),
                    "length_bp": len(vector_sequence),
                    "topology": vector.get("topology", ""),
                    "assembly_policy": vector.get("assembly_policy", ""),
                },
                "insert": {
                    "cassette_order": resolved_order,
                    "cassette_files": cassette_files,
                    "linker_sequence": cleaned_linker,
                    "length_bp": len(insert_sequence),
                },
                "homology_arm_length": homology_arm_length,
                "candidate_regions": candidate_regions,
                "recommendations": recommendations,
                "recommendation_count": len(recommendations),
                "next_action": (
                    "用户确认某个 Gibson 候选后，调用 plan_final_assembly，assembly_method='gibson'，"
                    "并传入该候选的 plan_final_assembly_args。"
                ),
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(recommend_gibson_insertion_sites.invoke({}))
