from __future__ import annotations
import json
import re
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


class RecommendRestrictionSitesArgs(BaseModel):
    cassette_order: list[int] | None = Field(
        default=None,
        description="表达盒顺序；为空时按 cassette_index 升序拼接 insert",
    )
    linker_sequence: str = Field(default="", description="多个 cassette 之间的连接序列，可为空")
    max_recommendations: int = Field(default=3, ge=1, le=10, description="最多返回多少对酶")
    enzyme_names: list[str] | None = Field(
        default=None,
        description="可选；限定候选酶名称，例如 ['EcoRI', 'HindIII']",
    )

    @field_validator("cassette_order", mode="before")
    @classmethod
    def normalize_cassette_order(cls, value: Any) -> list[int] | None:
        return _normalize_cassette_order(value)

    @field_validator("enzyme_names", mode="before")
    @classmethod
    def normalize_enzyme_names(cls, value: Any) -> list[str] | None:
        return _normalize_enzyme_names(value)


ENZYME_CATALOG: list[dict[str, Any]] = [
    {"name": "EcoRI", "site": "GAATTC", "priority": 95},
    {"name": "BamHI", "site": "GGATCC", "priority": 95},
    {"name": "HindIII", "site": "AAGCTT", "priority": 95},
    {"name": "XhoI", "site": "CTCGAG", "priority": 92},
    {"name": "XbaI", "site": "TCTAGA", "priority": 90},
    {"name": "SpeI", "site": "ACTAGT", "priority": 88},
    {"name": "NheI", "site": "GCTAGC", "priority": 88},
    {"name": "PstI", "site": "CTGCAG", "priority": 90},
    {"name": "KpnI", "site": "GGTACC", "priority": 88},
    {"name": "SacI", "site": "GAGCTC", "priority": 88},
    {"name": "SalI", "site": "GTCGAC", "priority": 86},
    {"name": "NotI", "site": "GCGGCCGC", "priority": 86},
    {"name": "SphI", "site": "GCATGC", "priority": 84},
    {"name": "BglII", "site": "AGATCT", "priority": 84},
    {"name": "AgeI", "site": "ACCGGT", "priority": 82},
    {"name": "AvrII", "site": "CCTAGG", "priority": 82},
    {"name": "BstBI", "site": "TTCGAA", "priority": 80},
    {"name": "ClaI", "site": "ATCGAT", "priority": 80},
    {"name": "MluI", "site": "ACGCGT", "priority": 80},
    {"name": "NcoI", "site": "CCATGG", "priority": 78},
    {"name": "NdeI", "site": "CATATG", "priority": 78},
    {"name": "PacI", "site": "TTAATTAA", "priority": 76},
    {"name": "AscI", "site": "GGCGCGCC", "priority": 76},
    {"name": "SbfI", "site": "CCTGCAGG", "priority": 76},
    {"name": "SmaI", "site": "CCCGGG", "priority": 72},
    {"name": "ApaI", "site": "GGGCCC", "priority": 72},
    {"name": "EagI", "site": "CGGCCG", "priority": 72},
]


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


def _normalize_enzyme_names(value: Any) -> list[str] | None:
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
        raise ValueError('enzyme_names must be a list of enzyme names, for example ["BamHI", "EcoRI"]')

    names = []
    seen = set()
    for item in value:
        if not isinstance(item, str):
            raise ValueError("enzyme_names values must be strings")
        name = item.strip()
        if not name:
            continue
        key = name.lower()
        if key in seen:
            continue
        seen.add(key)
        names.append(name)
    return names or None


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


def _read_sequence(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".gb", ".gbk", ".gbj"}:
        sequence = _read_genbank_sequence(path)
    else:
        sequence = _read_plain_or_fasta_sequence(path)

    sequence = "".join(base for base in sequence.upper() if base in "ACGTN")
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


def _genbank_mcs_regions(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() not in {".gb", ".gbk", ".gbj"}:
        return []

    try:
        from Bio import SeqIO
    except Exception:
        return []

    regions = []
    try:
        record = next(SeqIO.parse(str(path), "genbank"))
    except Exception:
        return []

    for feature in record.features:
        labels = []
        for key in ("label", "gene", "product", "note"):
            labels.extend(str(value) for value in feature.qualifiers.get(key, []))
        haystack = f"{feature.type} {' '.join(labels)}".lower()
        if "mcs" not in haystack and "multiple cloning" not in haystack and "polylinker" not in haystack:
            continue
        regions.append({
            "label": labels[0] if labels else feature.type,
            "start_bp": int(feature.location.start) + 1,
            "end_bp": int(feature.location.end),
        })
    return regions


def _site_positions(sequence: str, site: str, circular: bool = False) -> list[int]:
    sequence = sequence.upper()
    site = site.upper()
    if not sequence or not site:
        return []

    search_space = sequence
    if circular and len(sequence) >= len(site):
        search_space = sequence + sequence[: len(site) - 1]

    positions = []
    for match in re.finditer(f"(?={re.escape(site)})", search_space):
        index = match.start()
        if index < len(sequence):
            positions.append(index + 1)
    return positions


def _in_regions(position: int, regions: list[dict[str, Any]]) -> bool:
    return any(
        int(region.get("start_bp") or 0) <= position <= int(region.get("end_bp") or 0)
        for region in regions
    )


def _restriction_catalog(enzyme_names: list[str] | None) -> list[dict[str, Any]]:
    if not enzyme_names:
        return ENZYME_CATALOG
    allowed = {name.strip().lower() for name in enzyme_names if name.strip()}
    return [enzyme for enzyme in ENZYME_CATALOG if enzyme["name"].lower() in allowed]


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
    return [cassette_by_index[int(index)] for index in order]


def _build_insert_sequence(
    output_dir: Path,
    cassettes: list[dict[str, Any]],
    linker_sequence: str,
) -> tuple[str, list[dict[str, Any]]]:
    sequences = []
    files = []
    linker = "".join(base for base in str(linker_sequence or "").upper() if base in "ACGTN")

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

    return "".join(sequences), files


def _single_cutters(
    *,
    vector_sequence: str,
    insert_sequence: str,
    vector_is_circular: bool,
    mcs_regions: list[dict[str, Any]],
    enzyme_names: list[str] | None,
) -> list[dict[str, Any]]:
    cutters = []
    for enzyme in _restriction_catalog(enzyme_names):
        site = enzyme["site"]
        vector_positions = _site_positions(vector_sequence, site, circular=vector_is_circular)
        insert_positions = _site_positions(insert_sequence, site, circular=False)
        if len(vector_positions) != 1 or insert_positions:
            continue

        cut_position = vector_positions[0]
        cutters.append({
            "enzyme": enzyme["name"],
            "site": site,
            "vector_position_bp": cut_position,
            "insert_site_count": len(insert_positions),
            "in_mcs": _in_regions(cut_position, mcs_regions),
            "priority": enzyme.get("priority", 0),
        })
    return cutters


def _pair_score(left: dict[str, Any], right: dict[str, Any]) -> tuple[int, int, int, int]:
    mcs_score = int(bool(left.get("in_mcs"))) + int(bool(right.get("in_mcs")))
    distance = abs(int(left["vector_position_bp"]) - int(right["vector_position_bp"]))
    close_score = 1 if distance <= 250 else 0
    priority = int(left.get("priority") or 0) + int(right.get("priority") or 0)
    return (mcs_score, close_score, priority, -distance)


def _recommend_pairs(cutters: list[dict[str, Any]], max_recommendations: int) -> list[dict[str, Any]]:
    pairs = []
    for left_index, left in enumerate(cutters):
        for right in cutters[left_index + 1:]:
            ordered = sorted([left, right], key=lambda item: int(item["vector_position_bp"]))
            pair = {
                "left_enzyme": ordered[0]["enzyme"],
                "right_enzyme": ordered[1]["enzyme"],
                "left_site": ordered[0]["site"],
                "right_site": ordered[1]["site"],
                "left_vector_position_bp": ordered[0]["vector_position_bp"],
                "right_vector_position_bp": ordered[1]["vector_position_bp"],
                "both_in_mcs": bool(ordered[0]["in_mcs"] and ordered[1]["in_mcs"]),
                "rationale": [
                    "Both enzymes cut the vector exactly once.",
                    "Neither recognition site occurs in the assembled insert sequence.",
                ],
            }
            if pair["both_in_mcs"]:
                pair["rationale"].append("Both unique cut sites fall inside the vector MCS region.")
            pair["_score"] = _pair_score(ordered[0], ordered[1])
            pairs.append(pair)

    pairs.sort(key=lambda item: item["_score"], reverse=True)
    result = []
    for pair in pairs[:max_recommendations]:
        pair.pop("_score", None)
        result.append(pair)
    return result


def _seva_cargo_restriction_recommendation(
    *,
    vector: dict[str, Any],
    vector_sequence: str,
    insert_sequence: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    regions = [
        region
        for region in _as_list(vector.get("insertion_regions"))
        if isinstance(region, dict) and region.get("type") == "seva_cargo"
    ]
    if len(regions) != 1:
        raise ValueError("SEVA cargo replacement requires exactly one audited cargo region")
    region = regions[0]
    mcs_regions = [
        {
            "label": region.get("name") or "audited SEVA cargo",
            "start_bp": int(region.get("start_bp") or 0),
            "end_bp": int(region.get("end_bp") or 0),
        }
    ]
    cutters = []
    for enzyme, site in (("PacI", "TTAATTAA"), ("SpeI", "ACTAGT")):
        vector_positions = _site_positions(vector_sequence, site, circular=True)
        insert_positions = _site_positions(insert_sequence, site, circular=False)
        if len(vector_positions) != 1 or insert_positions:
            continue
        cutters.append(
            {
                "enzyme": enzyme,
                "site": site,
                "vector_position_bp": vector_positions[0],
                "insert_site_count": len(insert_positions),
                "in_mcs": True,
                "priority": 1000,
            }
        )
    recommendations = []
    if len(cutters) == 2:
        recommendations.append(
            {
                "left_enzyme": "PacI",
                "right_enzyme": "SpeI",
                "left_site": "TTAATTAA",
                "right_site": "ACTAGT",
                "left_vector_position_bp": next(
                    item["vector_position_bp"] for item in cutters if item["enzyme"] == "PacI"
                ),
                "right_vector_position_bp": next(
                    item["vector_position_bp"] for item in cutters if item["enzyme"] == "SpeI"
                ),
                "both_in_mcs": True,
                "rationale": [
                    "The audited vector policy requires complete PacI-SpeI cargo replacement.",
                    "PacI and SpeI each cut the circular vector once.",
                    "Neither recognition site occurs in the assembled insert sequence.",
                ],
            }
        )
    return mcs_regions, cutters, recommendations


@tool(args_schema=RecommendRestrictionSitesArgs)
def recommend_restriction_sites(
    cassette_order: list[int] | str | None = None,
    linker_sequence: str = "",
    max_recommendations: int = 3,
    enzyme_names: list[str] | str | None = None,
) -> str:
    """
    为当前 vector 和 assembled insert 推荐 restriction cloning 可用酶对。
    
    调用时机：用户选择 restriction，但未指定 left/right enzyme 和 site。
    返回：ok、vector/insert 摘要、single_cutters、recommendations 和 next_action。
    限制：只读；只推荐 vector 唯一切且不切 insert 的酶；用户确认后才调用 plan_final_assembly。
    """

    tool_name = "recommend_restriction_sites"
    monitor.report_start(tool_name, {"max_recommendations": max_recommendations, "enzyme_names": enzyme_names})
    try:
        cassette_order = _normalize_cassette_order(cassette_order)
        enzyme_names = _normalize_enzyme_names(enzyme_names)

        session_root = _session_root(None, None, None)
        outputs = _output_dir(session_root, None)
        manifest_file = _manifest_file(session_root, outputs, None)
        manifest = _read_json(manifest_file)

        vector = _as_dict(_as_dict(manifest.get("plasmid_selection")).get("vector"))
        vector_sequence_file = _as_dict(vector.get("sequence_file"))
        vector_path = _resolve_file_path(outputs, vector_sequence_file.get("path"))
        vector_sequence = _read_sequence(vector_path)
        vector_is_circular = str(vector.get("topology") or "").strip().lower() == "circular"

        cassettes = _ordered_cassettes(manifest, cassette_order)
        insert_sequence, cassette_files = _build_insert_sequence(outputs, cassettes, linker_sequence)
        if vector.get("assembly_policy") == "replace_seva_cargo_paci_spei":
            mcs_regions, cutters, recommendations = _seva_cargo_restriction_recommendation(
                vector=vector,
                vector_sequence=vector_sequence,
                insert_sequence=insert_sequence,
            )
        else:
            mcs_regions = _genbank_mcs_regions(vector_path)
            cutters = _single_cutters(
                vector_sequence=vector_sequence,
                insert_sequence=insert_sequence,
                vector_is_circular=vector_is_circular,
                mcs_regions=mcs_regions,
                enzyme_names=enzyme_names,
            )
            recommendations = _recommend_pairs(cutters, max_recommendations)

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
                    "cassette_order": [int(item.get("cassette_index")) for item in cassettes],
                    "cassette_files": cassette_files,
                    "linker_sequence": linker_sequence,
                    "length_bp": len(insert_sequence),
                },
                "mcs_regions": mcs_regions,
                "single_cutters": cutters,
                "recommendations": recommendations,
                "recommendation_count": len(recommendations),
                "next_action": (
                    "用户确认某一对酶后，调用 plan_final_assembly，assembly_method='restriction'，"
                    "并传入 left_enzyme/right_enzyme/left_site/right_site。"
                ),
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(recommend_restriction_sites.invoke({}))
