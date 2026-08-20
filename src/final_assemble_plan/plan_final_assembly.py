from __future__ import annotations

import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Literal

from langchain.tools import tool
from pydantic import BaseModel, Field, field_validator

from src.runtime.monitor import monitor
from src.tools.common.manifest import clear_downstream_fields
from src.tools.common.session_paths import (
    design_manifest_file,
    outputs_dir as resolve_outputs_dir,
    session_dir as resolve_session_dir,
)


class PlanFinalAssemblyArgs(BaseModel):
    construct_name: str = Field(default="final_construct", description="最终 construct 名称")
    cassette_order: list[int] | None = Field(
        default=None,
        description="表达盒组装顺序；为空时默认按 cassette_index 升序",
    )
    linker_sequence: str = Field(default="", description="多个表达盒之间的连接序列，可为空")
    notes: list[str] = Field(default_factory=list, description="用户确认的备注")

    assembly_method: Literal["direct", "restriction", "gibson"] = Field(
        default="direct",
        description="最终组装方法：direct、restriction 或 gibson",
    )

    insert_after_bp: int | None = Field(default=None, description="direct/gibson：在 vector 指定 bp 后插入")
    replace_start_bp: int | None = Field(default=None, description="direct/gibson：替换区域起点")
    replace_end_bp: int | None = Field(default=None, description="direct/gibson：替换区域终点")

    left_enzyme: str | None = Field(default=None, description="restriction：左侧限制性内切酶名称")
    right_enzyme: str | None = Field(default=None, description="restriction：右侧限制性内切酶名称")
    left_site: str | None = Field(default=None, description="restriction：左侧识别序列")
    right_site: str | None = Field(default=None, description="restriction：右侧识别序列")
    restriction_site_retention: Literal["retain", "remove"] = Field(
        default="retain",
        description="restriction：最终构建中保留或移除两端识别位点；新计划默认保留",
    )

    homology_arm_length: int = Field(default=30, description="gibson：同源臂长度")
    left_homology: str = Field(default="", description="gibson：左侧同源臂；应由 recommend_gibson_insertion_sites 生成")
    right_homology: str = Field(default="", description="gibson：右侧同源臂；应由 recommend_gibson_insertion_sites 生成")

    @field_validator("cassette_order", mode="before")
    @classmethod
    def normalize_cassette_order(cls, value: Any) -> list[int] | None:
        return _normalize_cassette_order(value)


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


def _write_json_atomic(path: Path, data: dict[str, Any]) -> None:
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
        raise ValueError("file path is missing")

    path = Path(str(path_value))
    if path.is_absolute():
        return path.resolve()

    output_relative = (output_dir / path).resolve()
    if output_relative.exists():
        return output_relative

    return (output_dir.parent / path).resolve()


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


def _clean_dna(value: Any) -> str:
    return "".join(base for base in str(value or "").upper() if base in "ACGTN")


def _cassette_summary(output_dir: Path, cassette: dict[str, Any]) -> dict[str, Any]:
    return {
        "cassette_index": cassette.get("cassette_index"),
        "format": cassette.get("format", ""),
        "path": _relpath_or_abs(output_dir, cassette.get("path")),
        "sha256": cassette.get("sha256", ""),
        "sequence_sha256": cassette.get("sequence_sha256", ""),
        "length_bp": cassette.get("length_bp"),
        "component_count": cassette.get("component_count"),
    }


def _vector_summary(output_dir: Path, vector: dict[str, Any]) -> dict[str, Any]:
    sequence_file = _as_dict(vector.get("sequence_file"))
    if sequence_file:
        sequence_file = {
            **sequence_file,
            "path": _relpath_or_abs(output_dir, sequence_file.get("path")),
        }

    return {
        "plasmid_id": vector.get("plasmid_id", ""),
        "addgene_id": vector.get("addgene_id", ""),
        "name": vector.get("name", ""),
        "replicon_family": vector.get("replicon_family", ""),
        "cargo_type": vector.get("cargo_type", ""),
        "assembly_policy": vector.get("assembly_policy", ""),
        "insertion_regions": _as_list(vector.get("insertion_regions")),
        "audit": _as_dict(vector.get("audit")),
        "topology": vector.get("topology", ""),
        "length_bp": vector.get("length_bp"),
        "sequence_file": sequence_file,
    }


def _strategy_payload(
    *,
    assembly_method: Literal["direct", "restriction", "gibson"],
    insert_after_bp: int | None,
    replace_start_bp: int | None,
    replace_end_bp: int | None,
    left_enzyme: str | None,
    right_enzyme: str | None,
    left_site: str | None,
    right_site: str | None,
    restriction_site_retention: Literal["retain", "remove"],
    homology_arm_length: int,
    left_homology: str,
    right_homology: str,
) -> dict[str, Any]:
    if assembly_method == "direct":
        return {
            "assembly_method": "direct",
            "insert_after_bp": insert_after_bp,
            "replace_start_bp": replace_start_bp,
            "replace_end_bp": replace_end_bp,
        }

    if assembly_method == "restriction":
        return {
            "assembly_method": "restriction",
            "left_enzyme": left_enzyme,
            "right_enzyme": right_enzyme,
            "left_site": left_site,
            "right_site": right_site,
            "restriction_site_retention": restriction_site_retention,
        }

    target_mode = "replace" if replace_start_bp is not None or replace_end_bp is not None else "insert_after"
    return {
        "assembly_method": "gibson",
        "target_mode": target_mode,
        "insert_after_bp": insert_after_bp,
        "replace_start_bp": replace_start_bp,
        "replace_end_bp": replace_end_bp,
        "homology_arm_length": homology_arm_length,
        "left_homology": left_homology,
        "right_homology": right_homology,
    }


def _validate_gibson_strategy(
    *,
    insert_after_bp: int | None,
    replace_start_bp: int | None,
    replace_end_bp: int | None,
    homology_arm_length: int,
    left_homology: str,
    right_homology: str,
) -> None:
    has_insert_after = insert_after_bp is not None
    has_replace_value = replace_start_bp is not None or replace_end_bp is not None

    if has_insert_after and has_replace_value:
        raise ValueError("gibson strategy must use either insert_after_bp or replace_start_bp/replace_end_bp, not both")
    if not has_insert_after and not has_replace_value:
        raise ValueError("gibson strategy requires insert_after_bp or replace_start_bp/replace_end_bp")

    if has_insert_after and int(insert_after_bp) < 0:
        raise ValueError("gibson insert_after_bp must be >= 0")

    if has_replace_value:
        if replace_start_bp is None or replace_end_bp is None:
            raise ValueError("gibson replacement requires both replace_start_bp and replace_end_bp")
        if int(replace_start_bp) < 1 or int(replace_end_bp) < 1:
            raise ValueError("gibson replace_start_bp and replace_end_bp must be >= 1")
        if int(replace_start_bp) > int(replace_end_bp):
            raise ValueError("gibson replace_start_bp must be <= replace_end_bp")

    if homology_arm_length < 15:
        raise ValueError("gibson homology_arm_length must be at least 15 bp")
    if not left_homology or not right_homology:
        raise ValueError("gibson strategy requires left_homology and right_homology")
    if len(left_homology) != homology_arm_length:
        raise ValueError("left_homology length must match homology_arm_length")
    if len(right_homology) != homology_arm_length:
        raise ValueError("right_homology length must match homology_arm_length")


def _validate_direct_strategy(
    *,
    insert_after_bp: int | None,
    replace_start_bp: int | None,
    replace_end_bp: int | None,
) -> None:
    has_insert_after = insert_after_bp is not None
    has_replace_value = replace_start_bp is not None or replace_end_bp is not None

    if has_insert_after and has_replace_value:
        raise ValueError("direct strategy must use either insert_after_bp or replace_start_bp/replace_end_bp, not both")
    if not has_insert_after and not has_replace_value:
        raise ValueError("direct strategy requires insert_after_bp or replace_start_bp/replace_end_bp")

    if has_insert_after and int(insert_after_bp) < 0:
        raise ValueError("direct insert_after_bp must be >= 0")

    if has_replace_value:
        if replace_start_bp is None or replace_end_bp is None:
            raise ValueError("direct replacement requires both replace_start_bp and replace_end_bp")
        if int(replace_start_bp) < 1 or int(replace_end_bp) < 1:
            raise ValueError("direct replace_start_bp and replace_end_bp must be >= 1")
        if int(replace_start_bp) > int(replace_end_bp):
            raise ValueError("direct replace_start_bp must be <= replace_end_bp")


def _validate_restriction_strategy(
    *,
    left_enzyme: str | None,
    right_enzyme: str | None,
    left_site: str | None,
    right_site: str | None,
) -> None:
    missing = [
        name
        for name, value in [
            ("left_enzyme", left_enzyme),
            ("right_enzyme", right_enzyme),
            ("left_site", left_site),
            ("right_site", right_site),
        ]
        if not str(value or "").strip()
    ]
    if missing:
        raise ValueError(f"restriction strategy missing required fields: {', '.join(missing)}")

    cleaned_left_site = _clean_dna(left_site)
    cleaned_right_site = _clean_dna(right_site)
    if len(cleaned_left_site) != len(str(left_site or "").strip().upper()):
        raise ValueError("left_site must contain only DNA bases A/C/G/T/N")
    if len(cleaned_right_site) != len(str(right_site or "").strip().upper()):
        raise ValueError("right_site must contain only DNA bases A/C/G/T/N")
    if not cleaned_left_site or not cleaned_right_site:
        raise ValueError("restriction sites must not be empty")


def _validate_ready_for_final_assembly(manifest: dict[str, Any], output_dir: Path) -> None:
    assembled = _as_dict(manifest.get("assembled_expression_cassettes"))
    if assembled.get("status") != "assembled":
        raise ValueError('assembled_expression_cassettes.status must be "assembled"')

    cassette_files = [
        item
        for item in _as_list(assembled.get("cassette_files"))
        if isinstance(item, dict)
    ]
    if not cassette_files:
        raise ValueError("assembled_expression_cassettes.cassette_files is empty")

    for item in cassette_files:
        cassette_index = item.get("cassette_index")
        path_value = item.get("path")
        if not path_value:
            raise ValueError(f"expression cassette file path is missing for cassette_index={cassette_index}")
        cassette_path = _resolve_file_path(output_dir, path_value)
        if not cassette_path.exists():
            raise FileNotFoundError(f"expression cassette file does not exist: {cassette_path}")
        if int(item.get("length_bp") or 0) <= 0:
            raise ValueError(f"expression cassette length_bp must be > 0 for cassette_index={cassette_index}")

    plasmid_selection = _as_dict(manifest.get("plasmid_selection"))
    if plasmid_selection.get("status") != "selected":
        raise ValueError('plasmid_selection.status must be "selected"')

    vector = _as_dict(plasmid_selection.get("vector"))
    if not vector:
        raise ValueError("plasmid_selection.vector is missing")

    sequence_file = _as_dict(vector.get("sequence_file"))
    vector_path = _resolve_file_path(output_dir, sequence_file.get("path"))
    if not vector_path.exists():
        raise FileNotFoundError(f"vector sequence_file does not exist: {vector_path}")
    if int(vector.get("length_bp") or sequence_file.get("length_bp") or 0) <= 0:
        raise ValueError("plasmid_selection.vector length_bp must be > 0")


def _validate_vector_assembly_policy(
    vector: dict[str, Any],
    *,
    assembly_method: Literal["direct", "restriction", "gibson"],
    insert_after_bp: int | None,
    replace_start_bp: int | None,
    replace_end_bp: int | None,
    left_enzyme: str | None,
    right_enzyme: str | None,
    left_site: str | None,
    right_site: str | None,
) -> None:
    policy = str(vector.get("assembly_policy") or "")
    if vector.get("has_expression_cargo") is True and (
        replace_start_bp is None or replace_end_bp is None
    ):
        raise ValueError(
            "selected vector contains expression cargo; final assembly must explicitly replace it"
        )
    if policy != "replace_seva_cargo_paci_spei":
        return

    regions = [
        region
        for region in _as_list(vector.get("insertion_regions"))
        if isinstance(region, dict) and region.get("type") == "seva_cargo"
    ]
    if len(regions) != 1:
        raise ValueError(
            "replace_seva_cargo_paci_spei requires exactly one audited SEVA cargo region"
        )
    region = regions[0]
    expected_start = int(region.get("start_bp") or 0)
    expected_end = int(region.get("end_bp") or 0)
    if expected_start <= 0 or expected_end < expected_start:
        raise ValueError("audited SEVA cargo replacement coordinates are invalid")

    if assembly_method in {"direct", "gibson"}:
        if insert_after_bp is not None:
            raise ValueError(
                "SEVA lacZ-alpha vectors must replace the PacI-SpeI cargo, not use insert_after_bp"
            )
        if replace_start_bp != expected_start or replace_end_bp != expected_end:
            raise ValueError(
                "SEVA cargo replacement coordinates must match the audited insertion region: "
                f"{expected_start}-{expected_end}"
            )
        return

    if (
        str(left_enzyme or "").strip().lower() != "paci"
        or str(right_enzyme or "").strip().lower() != "spei"
        or _clean_dna(left_site) != "TTAATTAA"
        or _clean_dna(right_site) != "ACTAGT"
    ):
        raise ValueError(
            "restriction assembly for this vector must replace the cargo with PacI/SpeI"
        )


@tool(args_schema=PlanFinalAssemblyArgs)
def plan_final_assembly(
    construct_name: str = "final_construct",
    cassette_order: list[int] | str | None = None,
    linker_sequence: str = "",
    notes: list[str] | None = None,
    assembly_method: Literal["direct", "restriction", "gibson"] = "direct",
    insert_after_bp: int | None = None,
    replace_start_bp: int | None = None,
    replace_end_bp: int | None = None,
    left_enzyme: str | None = None,
    right_enzyme: str | None = None,
    left_site: str | None = None,
    right_site: str | None = None,
    restriction_site_retention: Literal["retain", "remove"] = "retain",
    homology_arm_length: int = 30,
    left_homology: str = "",
    right_homology: str = "",
) -> str:
    """
    把用户确认的最终组装方案写入 design_manifest.json 的 final_assembly_plan。
    
    调用时机：用户已确认 direct、restriction 或 Gibson 的必要参数。
    返回：ok、manifest_path、revision、construct_name、assembly_method、cassette_order 和长度摘要。
    限制：只写计划；不推荐位点/酶对，不生成最终质粒；已有 final_assembly 会被标记为需重新执行。
    """

    tool_name = "plan_final_assembly"
    monitor.report_start(tool_name, {"construct_name": construct_name, "assembly_method": assembly_method})
    try:
        cassette_order = _normalize_cassette_order(cassette_order)

        if assembly_method == "direct":
            _validate_direct_strategy(
                insert_after_bp=insert_after_bp,
                replace_start_bp=replace_start_bp,
                replace_end_bp=replace_end_bp,
            )
        elif assembly_method == "restriction":
            _validate_restriction_strategy(
                left_enzyme=left_enzyme,
                right_enzyme=right_enzyme,
                left_site=left_site,
                right_site=right_site,
            )
            left_enzyme = str(left_enzyme or "").strip()
            right_enzyme = str(right_enzyme or "").strip()
            left_site = _clean_dna(left_site)
            right_site = _clean_dna(right_site)
        elif assembly_method == "gibson":
            left_homology = _clean_dna(left_homology)
            right_homology = _clean_dna(right_homology)
            _validate_gibson_strategy(
                insert_after_bp=insert_after_bp,
                replace_start_bp=replace_start_bp,
                replace_end_bp=replace_end_bp,
                homology_arm_length=homology_arm_length,
                left_homology=left_homology,
                right_homology=right_homology,
            )

        session_root = _session_root(None, None, None)
        outputs = _output_dir(session_root, None)
        manifest_file = _manifest_file(session_root, outputs, None)
        manifest = _read_json(manifest_file)
        _validate_ready_for_final_assembly(manifest, outputs)

        cassette_files = [
            item
            for item in _as_list(_as_dict(manifest.get("assembled_expression_cassettes")).get("cassette_files"))
            if isinstance(item, dict)
        ]
        cassette_by_index = {
            int(item.get("cassette_index")): item
            for item in cassette_files
        }
        resolved_order = cassette_order or sorted(cassette_by_index)
        if not resolved_order:
            raise ValueError("cassette_order is empty because no expression cassette files are available")
        missing_indexes = [
            int(index)
            for index in resolved_order
            if int(index) not in cassette_by_index
        ]
        if missing_indexes:
            raise ValueError(f"cassette_order references missing cassette_index values: {missing_indexes}")

        ordered_cassettes = [
            _cassette_summary(outputs, cassette_by_index[int(index)])
            for index in resolved_order
        ]

        vector = _as_dict(_as_dict(manifest.get("plasmid_selection")).get("vector"))
        _validate_vector_assembly_policy(
            vector,
            assembly_method=assembly_method,
            insert_after_bp=insert_after_bp,
            replace_start_bp=replace_start_bp,
            replace_end_bp=replace_end_bp,
            left_enzyme=left_enzyme,
            right_enzyme=right_enzyme,
            left_site=left_site,
            right_site=right_site,
        )
        vector_plan = _vector_summary(outputs, vector)
        insert_total_length = sum(int(item.get("length_bp") or 0) for item in ordered_cassettes)
        if insert_total_length <= 0:
            raise ValueError("final assembly insert total length must be > 0")

        plan = {
            "status": "planned",
            "source": "plan_final_assembly",
            "planned_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "construct_name": construct_name,
            "vector": vector_plan,
            "insert": {
                "cassette_order": [int(index) for index in resolved_order],
                "cassette_files": ordered_cassettes,
                "linker_sequence": linker_sequence,
                "total_length_bp": insert_total_length,
            },
            "strategy": _strategy_payload(
                assembly_method=assembly_method,
                insert_after_bp=insert_after_bp,
                replace_start_bp=replace_start_bp,
                replace_end_bp=replace_end_bp,
                left_enzyme=left_enzyme,
                right_enzyme=right_enzyme,
                left_site=left_site,
                right_site=right_site,
                restriction_site_retention=restriction_site_retention,
                homology_arm_length=homology_arm_length,
                left_homology=left_homology,
                right_homology=right_homology,
            ),
            "notes": notes or [],
        }

        current_revision = int(manifest.get("revision", 0))
        manifest["final_assembly_plan"] = plan
        clear_downstream_fields(manifest, "final_assembly_plan")
        manifest["revision"] = current_revision + 1
        _write_json_atomic(manifest_file, manifest)

        monitor.report_end(tool_name, {"revision": manifest["revision"], "assembly_method": assembly_method})
        return json.dumps(
            {
                "ok": True,
                "manifest_path": _relpath_or_abs(session_root, manifest_file),
                "revision": manifest["revision"],
                "construct_name": construct_name,
                "assembly_method": assembly_method,
                "cassette_order": plan["insert"]["cassette_order"],
                "cassette_count": len(ordered_cassettes),
                "insert_total_length_bp": insert_total_length,
                "vector_name": vector_plan.get("name", ""),
            },
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        return _json_error(type(exc).__name__, message=str(exc))


if __name__ == "__main__":
    print(plan_final_assembly.invoke({}))
