from __future__ import annotations

import hashlib
import json
import tempfile
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any

from src.tools.common.session_paths import (
    design_manifest_file,
    final_assembly_dir as resolve_context_final_assembly_dir,
    outputs_dir as resolve_context_outputs_dir,
    session_dir as resolve_context_session_dir,
)

GENBANK_EXTENSIONS = {".gb", ".gbk", ".gbj"}


def json_error(error: str, **details: Any) -> str:
    return json.dumps({"ok": False, "error": error, **details}, ensure_ascii=False, default=str)


def session_root(session_dir: str | None, output_dir: str | None, manifest_path: str | None) -> Path:
    if session_dir:
        return Path(session_dir).resolve()

    if output_dir:
        output_path = Path(output_dir).resolve()
        return output_path.parent.resolve() if output_path.name == "outputs" else output_path.parent.resolve()

    if manifest_path:
        path = Path(manifest_path).resolve()
        if path.name == "design_manifest.json" and path.parent.name == "outputs":
            return path.parent.parent.resolve()

    return resolve_context_session_dir()


def outputs_dir(root: Path, output_dir: str | None) -> Path:
    if output_dir:
        return Path(output_dir).resolve()

    return resolve_context_outputs_dir()


def manifest_file(root: Path, output_dir: Path, manifest_path: str | None) -> Path:
    if manifest_path:
        path = Path(manifest_path)
        if not path.is_absolute():
            output_relative = output_dir / path
            session_relative = root / path
            path = output_relative if output_relative.exists() else session_relative
        return path.resolve()

    return design_manifest_file()


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as file:
        data = json.load(file)
    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return data


def write_json_atomic(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
        suffix=".tmp",
    ) as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        tmp_path = Path(handle.name)
    tmp_path.replace(path)


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def clean_dna(value: Any) -> str:
    return "".join(base for base in str(value or "").upper() if base in "ACGTN")


def sanitize_name(value: Any) -> str:
    text = str(value or "").strip() or "final_construct"
    return "".join(char if char.isalnum() or char in {"_", "-", "."} else "_" for char in text)


def relpath_or_abs(base_dir: Path, path_value: Any) -> str:
    if not path_value:
        return ""
    path = Path(str(path_value))
    if not path.is_absolute():
        path = base_dir / path
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def resolve_file_path(output_dir: Path, path_value: Any) -> Path:
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_genbank_origin(path: Path) -> str:
    in_origin = False
    chunks: list[str] = []
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


def read_genbank(path: Path) -> tuple[Any | None, str, list[str]]:
    if path.suffix.lower() not in GENBANK_EXTENSIONS:
        raise ValueError(f"expected GenBank file (.gb/.gbk/.gbj): {path}")

    parse_warnings: list[str] = []
    try:
        from Bio import SeqIO

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with path.open("r", encoding="utf-8") as handle:
                record = next(SeqIO.parse(handle, "genbank"))
        parse_warnings = [str(item.message) for item in caught]
        sequence = clean_dna(str(record.seq))
        if not sequence:
            raise ValueError("parsed GenBank sequence is empty")
        return record, sequence, parse_warnings
    except Exception as exc:
        sequence = clean_dna(read_genbank_origin(path))
        if not sequence:
            raise ValueError(f"could not parse GenBank sequence: {path}") from exc
        parse_warnings.append(f"Biopython parse failed; recovered sequence from ORIGIN: {exc}")
        return None, sequence, parse_warnings


def feature_label(feature: Any) -> str:
    for key in ("label", "gene", "product", "note"):
        values = getattr(feature, "qualifiers", {}).get(key, [])
        if values:
            return str(values[0])
    return str(getattr(feature, "type", "misc_feature"))


def feature_note(feature: Any) -> str:
    values = getattr(feature, "qualifiers", {}).get("note", [])
    return str(values[0]) if values else feature_label(feature)


def file_summary(output_dir: Path, path: Path, *, sequence: str, file_format: str) -> dict[str, Any]:
    return {
        "path": relpath_or_abs(output_dir, path),
        "sha256": sha256_file(path),
        "sequence_sha256": sha256_text(sequence),
        "length_bp": len(sequence),
        "format": file_format,
    }


def validate_file_digest(path: Path, file_info: dict[str, Any], sequence: str) -> list[str]:
    warnings_list = []
    expected_sha = str(file_info.get("sha256") or "").strip()
    if expected_sha and expected_sha != sha256_file(path):
        warnings_list.append(f"sha256 mismatch for {path.name}")

    expected_sequence_sha = str(file_info.get("sequence_sha256") or "").strip()
    if expected_sequence_sha and expected_sequence_sha != sha256_text(sequence):
        warnings_list.append(f"sequence_sha256 mismatch for {path.name}")

    expected_length = file_info.get("length_bp")
    if expected_length is not None and int(expected_length) != len(sequence):
        warnings_list.append(f"length mismatch for {path.name}: expected {expected_length}, actual {len(sequence)}")

    return warnings_list


def site_positions(sequence: str, site: str) -> list[int]:
    site = clean_dna(site)
    if not sequence or not site:
        return []
    return [
        index + 1
        for index in range(0, len(sequence) - len(site) + 1)
        if sequence[index:index + len(site)] == site
    ]


def sequence_segment(sequence: str, start_bp: int, length: int, *, circular: bool) -> str:
    if length <= 0:
        raise ValueError("segment length must be positive")
    if not circular and (start_bp < 1 or start_bp + length - 1 > len(sequence)):
        raise ValueError("segment is outside linear vector bounds")
    if circular:
        return "".join(sequence[(start_bp - 1 + offset) % len(sequence)] for offset in range(length))
    return sequence[start_bp - 1:start_bp + length - 1]


def get_plan(manifest: dict[str, Any]) -> dict[str, Any]:
    plan = as_dict(manifest.get("final_assembly_plan"))
    if plan.get("status") != "planned":
        raise ValueError('final_assembly_plan.status must be "planned"')
    strategy = as_dict(plan.get("strategy"))
    method = strategy.get("assembly_method")
    if method not in {"direct", "restriction", "gibson"}:
        raise ValueError("final_assembly_plan.strategy.assembly_method must be direct, restriction, or gibson")
    return plan


def assembled_cassette_by_index(manifest: dict[str, Any]) -> dict[int, dict[str, Any]]:
    assembled = as_dict(manifest.get("assembled_expression_cassettes"))
    if assembled.get("status") != "assembled":
        raise ValueError('assembled_expression_cassettes.status must be "assembled"')

    cassette_files = [item for item in as_list(assembled.get("cassette_files")) if isinstance(item, dict)]
    if not cassette_files:
        raise ValueError("assembled_expression_cassettes.cassette_files is empty")
    return {int(item.get("cassette_index")): item for item in cassette_files}


def planned_cassettes(manifest: dict[str, Any], plan: dict[str, Any]) -> list[dict[str, Any]]:
    cassette_by_index = assembled_cassette_by_index(manifest)
    order = [int(index) for index in as_list(as_dict(plan.get("insert")).get("cassette_order"))]
    if not order:
        order = sorted(cassette_by_index)

    missing = [index for index in order if index not in cassette_by_index]
    if missing:
        raise ValueError(f"final_assembly_plan references missing cassette_index values: {missing}")
    return [cassette_by_index[index] for index in order]


def load_vector(manifest: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    plasmid_selection = as_dict(manifest.get("plasmid_selection"))
    if plasmid_selection.get("status") != "selected":
        raise ValueError('plasmid_selection.status must be "selected"')

    vector = as_dict(plasmid_selection.get("vector"))
    sequence_file = as_dict(vector.get("sequence_file"))
    vector_path = resolve_file_path(output_dir, sequence_file.get("path"))
    if not vector_path.exists():
        raise FileNotFoundError(vector_path)

    record, sequence, parse_warnings = read_genbank(vector_path)
    digest_warnings = validate_file_digest(vector_path, sequence_file, sequence)
    topology = str(vector.get("topology") or getattr(record, "annotations", {}).get("topology", "") or "").lower()

    return {
        "manifest_vector": vector,
        "path": vector_path,
        "record": record,
        "sequence": sequence,
        "topology": topology or "circular",
        "warnings": parse_warnings + digest_warnings,
    }


def load_insert(manifest: dict[str, Any], plan: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    cassettes = planned_cassettes(manifest, plan)
    linker_sequence = clean_dna(as_dict(plan.get("insert")).get("linker_sequence"))

    cursor = 1
    chunks: list[str] = []
    files: list[dict[str, Any]] = []
    insert_components: list[dict[str, Any]] = []

    for index, cassette in enumerate(cassettes):
        if index > 0 and linker_sequence:
            linker_start = cursor
            chunks.append(linker_sequence)
            cursor += len(linker_sequence)
            insert_components.append({
                "type": "linker",
                "label": f"linker_before_cassette_{cassette.get('cassette_index')}",
                "start_bp": linker_start,
                "end_bp": cursor - 1,
                "source": "final_assembly_plan.insert.linker_sequence",
            })

        cassette_path = resolve_file_path(output_dir, cassette.get("path"))
        if not cassette_path.exists():
            raise FileNotFoundError(cassette_path)

        _, sequence, parse_warnings = read_genbank(cassette_path)
        digest_warnings = validate_file_digest(cassette_path, cassette, sequence)
        cassette_start = cursor
        chunks.append(sequence)
        cursor += len(sequence)
        cassette_end = cursor - 1

        files.append({
            "cassette_index": cassette.get("cassette_index"),
            "path": relpath_or_abs(output_dir, cassette_path),
            "length_bp": len(sequence),
            "sha256": sha256_file(cassette_path),
            "sequence_sha256": sha256_text(sequence),
            "parse_warnings": parse_warnings,
            "validation_warnings": digest_warnings,
        })
        insert_components.append({
            "type": "expression_cassette",
            "label": f"cassette_{cassette.get('cassette_index')}",
            "cassette_index": cassette.get("cassette_index"),
            "start_bp": cassette_start,
            "end_bp": cassette_end,
            "source": relpath_or_abs(output_dir, cassette_path),
        })

        for component in as_list(cassette.get("components")):
            if not isinstance(component, dict):
                continue
            component_start = int(component.get("start_bp") or 0)
            component_end = int(component.get("end_bp") or 0)
            if component_start <= 0 or component_end <= 0:
                continue
            insert_components.append({
                "type": component.get("type", "misc_feature"),
                "label": component.get("label") or component.get("source_id") or component.get("cds_id") or "component",
                "source_id": component.get("source_id", ""),
                "cds_id": component.get("cds_id", ""),
                "protein_name": component.get("protein_name", ""),
                "cassette_index": cassette.get("cassette_index"),
                "start_bp": cassette_start + component_start - 1,
                "end_bp": cassette_start + component_end - 1,
                "source": as_dict(component.get("sequence_file")).get("path", ""),
            })

    sequence = "".join(chunks)
    if not sequence:
        raise ValueError("assembled insert sequence is empty")

    expected_length = as_dict(plan.get("insert")).get("total_length_bp")
    warnings_list = []
    if expected_length is not None and int(expected_length) != len(sequence):
        warnings_list.append(f"insert length mismatch: expected {expected_length}, actual {len(sequence)}")

    for item in files:
        warnings_list.extend(item["parse_warnings"])
        warnings_list.extend(item["validation_warnings"])

    return {
        "sequence": sequence,
        "cassette_files": files,
        "components": insert_components,
        "linker_sequence": linker_sequence,
        "warnings": warnings_list,
    }


def validate_target_bounds(vector_length: int, start_bp: int, end_bp: int) -> None:
    if start_bp < 1 or end_bp < 1:
        raise ValueError("target coordinates must be >= 1")
    if start_bp > end_bp:
        raise ValueError("replace_start_bp must be <= replace_end_bp")
    if end_bp > vector_length:
        raise ValueError("target coordinates exceed vector length")


def assembly_target(
    *,
    vector_sequence: str,
    insert_sequence: str,
    strategy: dict[str, Any],
) -> dict[str, Any]:
    method = strategy.get("assembly_method")
    vector_length = len(vector_sequence)

    if method in {"direct", "gibson"}:
        insert_after_bp = strategy.get("insert_after_bp")
        replace_start_bp = strategy.get("replace_start_bp")
        replace_end_bp = strategy.get("replace_end_bp")

        if insert_after_bp is not None:
            cut = int(insert_after_bp)
            if cut < 0 or cut > vector_length:
                raise ValueError("insert_after_bp must be between 0 and vector length")
            inserted_start = cut + 1
            return {
                "mode": "insert_after",
                "insert_after_bp": cut,
                "replace_start_bp": None,
                "replace_end_bp": None,
                "replaced_length_bp": 0,
                "inserted_start_bp": inserted_start,
                "inserted_end_bp": inserted_start + len(insert_sequence) - 1,
                "final_sequence": vector_sequence[:cut] + insert_sequence + vector_sequence[cut:],
            }

        if replace_start_bp is None or replace_end_bp is None:
            raise ValueError(f"{method} strategy requires insert_after_bp or replace_start_bp/replace_end_bp")

        start_bp = int(replace_start_bp)
        end_bp = int(replace_end_bp)
        validate_target_bounds(vector_length, start_bp, end_bp)
        inserted_start = start_bp
        return {
            "mode": "replace",
            "insert_after_bp": None,
            "replace_start_bp": start_bp,
            "replace_end_bp": end_bp,
            "replaced_length_bp": end_bp - start_bp + 1,
            "inserted_start_bp": inserted_start,
            "inserted_end_bp": inserted_start + len(insert_sequence) - 1,
            "final_sequence": vector_sequence[:start_bp - 1] + insert_sequence + vector_sequence[end_bp:],
        }

    left_site = clean_dna(strategy.get("left_site"))
    right_site = clean_dna(strategy.get("right_site"))
    if not left_site or not right_site:
        raise ValueError("restriction strategy requires left_site and right_site")
    if site_positions(insert_sequence, left_site) or site_positions(insert_sequence, right_site):
        raise ValueError("restriction site occurs inside assembled insert sequence")

    left_positions = site_positions(vector_sequence, left_site)
    right_positions = site_positions(vector_sequence, right_site)
    if len(left_positions) != 1:
        raise ValueError(f"left restriction site must occur exactly once in vector, found {len(left_positions)}")
    if len(right_positions) != 1:
        raise ValueError(f"right restriction site must occur exactly once in vector, found {len(right_positions)}")

    left_start = left_positions[0]
    right_start = right_positions[0]
    if left_start >= right_start:
        raise ValueError("restriction execution requires left_site to occur before right_site in the vector sequence")

    replace_start_bp = left_start
    replace_end_bp = right_start + len(right_site) - 1
    validate_target_bounds(vector_length, replace_start_bp, replace_end_bp)
    retention = str(strategy.get("restriction_site_retention") or "remove").strip().lower()
    if retention not in {"retain", "remove"}:
        raise ValueError("restriction_site_retention must be 'retain' or 'remove'")

    if retention == "retain":
        replacement_payload = left_site + insert_sequence + right_site
        inserted_start = replace_start_bp + len(left_site)
        left_site_retained_start = replace_start_bp
        right_site_retained_start = inserted_start + len(insert_sequence)
        mode = "restriction_replace_retain_sites"
    else:
        replacement_payload = insert_sequence
        inserted_start = replace_start_bp
        left_site_retained_start = None
        right_site_retained_start = None
        mode = "restriction_replace"

    junction_flank_length = min(20, len(insert_sequence))
    return {
        "mode": mode,
        "left_enzyme": strategy.get("left_enzyme"),
        "right_enzyme": strategy.get("right_enzyme"),
        "left_site": left_site,
        "right_site": right_site,
        "restriction_site_retention": retention,
        "left_site_position_bp": left_start,
        "right_site_position_bp": right_start,
        "left_site_retained_start_bp": left_site_retained_start,
        "left_site_retained_end_bp": (
            left_site_retained_start + len(left_site) - 1
            if left_site_retained_start is not None
            else None
        ),
        "right_site_retained_start_bp": right_site_retained_start,
        "right_site_retained_end_bp": (
            right_site_retained_start + len(right_site) - 1
            if right_site_retained_start is not None
            else None
        ),
        "restriction_junctions": {
            "left": {
                "enzyme": strategy.get("left_enzyme"),
                "recognition_site": left_site,
                "sequence": left_site + insert_sequence[:junction_flank_length]
                if retention == "retain"
                else insert_sequence[:junction_flank_length],
            },
            "right": {
                "enzyme": strategy.get("right_enzyme"),
                "recognition_site": right_site,
                "sequence": insert_sequence[-junction_flank_length:] + right_site
                if retention == "retain"
                else insert_sequence[-junction_flank_length:],
            },
        },
        "primer_tail_requirements": {
            "left_forward": {
                "enzyme": strategy.get("left_enzyme"),
                "recognition_site_tail": left_site,
                "protective_clamp": "not_specified",
            },
            "right_reverse": {
                "enzyme": strategy.get("right_enzyme"),
                "recognition_site_tail": right_site,
                "protective_clamp": "not_specified",
            },
        },
        "insert_after_bp": None,
        "replace_start_bp": replace_start_bp,
        "replace_end_bp": replace_end_bp,
        "replaced_length_bp": replace_end_bp - replace_start_bp + 1,
        "replacement_payload_length_bp": len(replacement_payload),
        "inserted_start_bp": inserted_start,
        "inserted_end_bp": inserted_start + len(insert_sequence) - 1,
        "final_sequence": (
            vector_sequence[:replace_start_bp - 1]
            + replacement_payload
            + vector_sequence[replace_end_bp:]
        ),
    }


def validate_gibson_homology(vector: dict[str, Any], strategy: dict[str, Any], target: dict[str, Any]) -> None:
    if strategy.get("assembly_method") != "gibson":
        return

    homology_arm_length = int(strategy.get("homology_arm_length") or 0)
    left_homology = clean_dna(strategy.get("left_homology"))
    right_homology = clean_dna(strategy.get("right_homology"))
    if homology_arm_length <= 0 or not left_homology or not right_homology:
        raise ValueError("gibson strategy requires homology_arm_length, left_homology, and right_homology")
    if len(left_homology) != homology_arm_length or len(right_homology) != homology_arm_length:
        raise ValueError("gibson homology arm lengths do not match homology_arm_length")

    vector_sequence = vector["sequence"]
    circular = str(vector.get("topology") or "").lower() == "circular"
    if target["mode"] == "insert_after":
        insert_after = int(target["insert_after_bp"])
        left_start = insert_after - homology_arm_length + 1
        right_start = insert_after + 1
    else:
        left_start = int(target["replace_start_bp"]) - homology_arm_length
        right_start = int(target["replace_end_bp"]) + 1

    expected_left = sequence_segment(vector_sequence, left_start, homology_arm_length, circular=circular)
    expected_right = sequence_segment(vector_sequence, right_start, homology_arm_length, circular=circular)
    if expected_left != left_homology:
        raise ValueError("left_homology does not match vector sequence around planned insertion site")
    if expected_right != right_homology:
        raise ValueError("right_homology does not match vector sequence around planned insertion site")


def map_vector_feature(feature: Any, target: dict[str, Any], insert_length: int) -> tuple[dict[str, Any] | None, str | None]:
    if getattr(feature, "location", None) is None or getattr(feature, "type", "") == "source":
        return None, None

    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    replacement_payload_length = int(target.get("replacement_payload_length_bp", insert_length))
    shift = replacement_payload_length - int(target["replaced_length_bp"])
    replace_start = target.get("replace_start_bp")
    replace_end = target.get("replace_end_bp")
    insert_after = target.get("insert_after_bp")

    if replace_start is not None and replace_end is not None:
        replace_start = int(replace_start)
        replace_end = int(replace_end)
        if end < replace_start:
            mapped_start, mapped_end = start, end
        elif start > replace_end:
            mapped_start, mapped_end = start + shift, end + shift
        elif start >= replace_start and end <= replace_end:
            return None, f"dropped vector feature inside replaced region: {feature_label(feature)}"
        else:
            return None, f"dropped vector feature overlapping replaced region: {feature_label(feature)}"
    else:
        cut = int(insert_after)
        if end <= cut:
            mapped_start, mapped_end = start, end
        elif start > cut:
            mapped_start, mapped_end = start + insert_length, end + insert_length
        else:
            mapped_start, mapped_end = start, end + insert_length
            return {
                "type": getattr(feature, "type", "misc_feature"),
                "label": feature_label(feature),
                "note": f"{feature_note(feature)}; feature spans inserted region",
                "start_bp": mapped_start,
                "end_bp": mapped_end,
                "strand": getattr(feature.location, "strand", None),
                "source": "vector",
            }, f"vector feature spans insertion site and was expanded: {feature_label(feature)}"

    return {
        "type": getattr(feature, "type", "misc_feature"),
        "label": feature_label(feature),
        "note": feature_note(feature),
        "start_bp": mapped_start,
        "end_bp": mapped_end,
        "strand": getattr(feature.location, "strand", None),
        "source": "vector",
    }, None


def final_features(vector: dict[str, Any], insert: dict[str, Any], target: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    features: list[dict[str, Any]] = []
    warnings_list: list[str] = []
    insert_length = len(insert["sequence"])
    record = vector.get("record")

    if record is not None:
        for feature in getattr(record, "features", []):
            mapped, warning = map_vector_feature(feature, target, insert_length)
            if mapped is not None:
                features.append(mapped)
            if warning:
                warnings_list.append(warning)

    inserted_start = int(target["inserted_start_bp"])
    if target.get("restriction_site_retention") == "retain":
        for side in ("left", "right"):
            start = target.get(f"{side}_site_retained_start_bp")
            end = target.get(f"{side}_site_retained_end_bp")
            if start is None or end is None:
                continue
            enzyme = target.get(f"{side}_enzyme") or side
            features.append({
                "type": "misc_feature",
                "label": f"{enzyme}_retained_restriction_site",
                "note": (
                    f"Retained {enzyme} recognition site "
                    f"({target.get(f'{side}_site', '')}) at the insert junction"
                ),
                "start_bp": int(start),
                "end_bp": int(end),
                "strand": 1,
                "source": "final_assembly_plan.restriction_site",
            })

    for component in insert["components"]:
        if component.get("type") == "expression_cassette":
            continue

        start = inserted_start + int(component["start_bp"]) - 1
        end = inserted_start + int(component["end_bp"]) - 1
        features.append({
            **component,
            "note": component.get("note")
            or component.get("protein_name")
            or component.get("source_id")
            or component.get("label", ""),
            "start_bp": start,
            "end_bp": end,
            "source": component.get("source", "insert"),
        })

    features = [
        feature
        for feature in features
        if int(feature.get("start_bp") or 0) >= 1
        and int(feature.get("end_bp") or 0) >= int(feature.get("start_bp") or 0)
    ]
    features.sort(key=lambda item: (int(item["start_bp"]), int(item["end_bp"])))
    return features, warnings_list


def biopython_feature(feature: dict[str, Any]) -> Any:
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    start = int(feature["start_bp"]) - 1
    end = int(feature["end_bp"])
    strand = feature.get("strand")
    feature_type = str(feature.get("type") or "misc_feature")
    if feature_type.lower() == "rbs":
        feature_type = "RBS"
    elif feature_type.lower() == "cds":
        feature_type = "CDS"
    if len(feature_type) > 15 or feature_type in {"expression_cassette", "linker"}:
        feature_type = "misc_feature"

    qualifiers = {"label": [str(feature.get("label") or feature.get("type") or "feature")]}
    note = str(feature.get("note") or "").strip()
    if note:
        qualifiers["note"] = [note]
    if feature.get("cds_id"):
        qualifiers["gene"] = [str(feature["cds_id"])]
    if feature.get("protein_name"):
        qualifiers["product"] = [str(feature["protein_name"])]

    return SeqFeature(
        SimpleLocation(start, end, strand=strand if strand in {-1, 1} else None),
        type=feature_type,
        qualifiers=qualifiers,
    )


def write_fasta(path: Path, record_id: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f">{record_id}"]
    for index in range(0, len(sequence), 80):
        lines.append(sequence[index:index + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_genbank(path: Path, record_id: str, sequence: str, topology: str, features: list[dict[str, Any]]) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(
        Seq(sequence),
        id=record_id[:16],
        name=record_id[:16],
        description=f"{record_id} final assembled construct",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular" if topology == "circular" else "linear"
    record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
    record.features = [biopython_feature(feature) for feature in features]
    SeqIO.write(record, str(path), "genbank")


def read_written_sequence(path: Path, file_format: str) -> tuple[str, list[str]]:
    from Bio import SeqIO

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with path.open("r", encoding="utf-8") as handle:
            record = next(SeqIO.parse(handle, file_format))
    return clean_dna(str(record.seq)), [str(item.message) for item in caught]


def validate_written_outputs(genbank_path: Path, fasta_path: Path, expected_sequence: str) -> dict[str, Any]:
    genbank_sequence, genbank_warnings = read_written_sequence(genbank_path, "genbank")
    fasta_sequence, fasta_warnings = read_written_sequence(fasta_path, "fasta")
    if genbank_sequence != expected_sequence:
        raise ValueError("written GenBank sequence does not match assembled final sequence")
    if fasta_sequence != expected_sequence:
        raise ValueError("written FASTA sequence does not match assembled final sequence")
    return {
        "genbank_parse_ok": True,
        "fasta_parse_ok": True,
        "genbank_sequence_matches": True,
        "fasta_sequence_matches": True,
        "genbank_parse_warnings": genbank_warnings,
        "fasta_parse_warnings": fasta_warnings,
    }


def prepare_context(output_dir: str | None, session_dir: str | None, manifest_path: str | None) -> dict[str, Any]:
    root = session_root(session_dir, output_dir, manifest_path)
    outputs = outputs_dir(root, output_dir)
    manifest = manifest_file(root, outputs, manifest_path)
    manifest_data = read_json(manifest)
    plan = get_plan(manifest_data)
    vector = load_vector(manifest_data, outputs)
    insert = load_insert(manifest_data, plan, outputs)
    strategy = as_dict(plan.get("strategy"))
    target = assembly_target(
        vector_sequence=vector["sequence"],
        insert_sequence=insert["sequence"],
        strategy=strategy,
    )
    validate_gibson_homology(vector, strategy, target)
    warnings_list = vector["warnings"] + insert["warnings"]
    return {
        "session_root": root,
        "outputs": outputs,
        "manifest_file": manifest,
        "manifest": manifest_data,
        "plan": plan,
        "vector": vector,
        "insert": insert,
        "strategy": strategy,
        "target": target,
        "warnings": warnings_list,
    }


def final_assembly_dir(outputs: Path, final_assembly_dir_value: str | None) -> Path:
    return Path(final_assembly_dir_value).resolve() if final_assembly_dir_value else resolve_context_final_assembly_dir()


def final_assembly_summary(manifest: dict[str, Any], outputs: Path) -> dict[str, Any]:
    final_assembly = as_dict(manifest.get("final_assembly"))
    if not final_assembly:
        return {"available": False}

    construct_files = as_dict(final_assembly.get("construct_files"))
    resolved_files = {}
    file_warnings = []
    for key, info in construct_files.items():
        file_info = info if isinstance(info, dict) else {"path": info}
        path_value = file_info.get("path")
        if not path_value:
            resolved_files[key] = {"exists": False, "path": ""}
            continue
        path = resolve_file_path(outputs, path_value)
        exists = path.exists()
        actual_sha = sha256_file(path) if exists else ""
        expected_sha = str(file_info.get("sha256") or "")
        if exists and expected_sha and expected_sha != actual_sha:
            file_warnings.append(f"{key} sha256 mismatch")
        resolved_files[key] = {
            **file_info,
            "path": relpath_or_abs(outputs, path),
            "exists": exists,
            "sha256_matches": bool(exists and (not expected_sha or expected_sha == actual_sha)),
        }

    return {
        "available": True,
        "status": final_assembly.get("status"),
        "source": final_assembly.get("source"),
        "assembled_at": final_assembly.get("assembled_at"),
        "construct_name": final_assembly.get("construct_name"),
        "assembly_method": as_dict(final_assembly.get("strategy")).get("assembly_method"),
        "length_bp": final_assembly.get("length_bp"),
        "topology": final_assembly.get("topology"),
        "component_count": final_assembly.get("component_count"),
        "construct_files": resolved_files,
        "validation": as_dict(final_assembly.get("validation")),
        "warnings": as_list(final_assembly.get("warnings")) + file_warnings,
    }
