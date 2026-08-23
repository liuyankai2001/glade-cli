"""Build auxiliary-protein research units from the selected design manifest."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from src.main_protein_selection.provenance import (
    solution_fingerprint,
    stable_json_hash,
)
from src.protein_selection.state import (
    AssignedReactionStep,
    MainEnzymeResearchUnit,
    WholeReactionContext,
)
from src.write_manifest.main_enzyme import (
    MAIN_ENZYME_MANIFEST_SCHEMA_VERSION,
)
from src.write_manifest.store import (
    SCHEMA_VERSION as DESIGN_MANIFEST_SCHEMA_VERSION,
    read_design_manifest,
)


_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_ALLOWED_SELECTION_STATUSES = {
    "user_selected",
    "user_selected_pending_review",
}


def _string(value: Any) -> str:
    return str(value or "").strip()


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest 缺少有效的 {field_name}")
    return value


def _mapping_list(value: Any, field_name: str) -> list[Mapping[str, Any]]:
    if not isinstance(value, list) or any(
        not isinstance(item, Mapping) for item in value
    ):
        raise ValueError(f"manifest 中的 {field_name} 必须是对象列表")
    return value


def _integer(value: Any, field_name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ValueError(f"manifest 中的 {field_name} 必须是整数")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"manifest 中的 {field_name} 必须是整数") from exc
    if result < minimum:
        raise ValueError(
            f"manifest 中的 {field_name} 必须大于等于 {minimum}"
        )
    return result


def _sorted_unique_step_indexes(value: Any, field_name: str) -> list[int]:
    if not isinstance(value, list):
        raise ValueError(f"manifest 中的 {field_name} 必须是 Step 列表")
    result = [
        _integer(item, field_name, minimum=1)
        for item in value
    ]
    if result != sorted(set(result)):
        raise ValueError(
            f"manifest 中的 {field_name} 必须按升序排列且不能重复"
        )
    return result


def _split_values(value: Any) -> list[str]:
    if value is None:
        return []
    raw_values: list[Any]
    if isinstance(value, Sequence) and not isinstance(value, str):
        raw_values = list(value)
    else:
        raw_values = re.split(r"\s*[;|]\s*", str(value))
    normalized = [
        str(item).strip()
        for item in raw_values
        if str(item).strip()
    ]
    return list(dict.fromkeys(normalized))


def _route_step_index(
    steps: list[Mapping[str, Any]],
) -> dict[int, Mapping[str, Any]]:
    result: dict[int, Mapping[str, Any]] = {}
    for step in steps:
        step_index = _integer(
            step.get("step_index"),
            "solution.steps[].step_index",
            minimum=1,
        )
        if step_index in result:
            raise ValueError(f"manifest 路线包含重复 Step {step_index}")
        result[step_index] = step
    if not result:
        raise ValueError("manifest 中的路线没有反应步骤")
    return result


def _validate_route_lock(
    solution: Mapping[str, Any],
    selection: Mapping[str, Any],
    route_steps: list[Mapping[str, Any]],
) -> str:
    solution_id = _integer(
        solution.get("solution_id"),
        "solution.solution_id",
        minimum=1,
    )
    expansion_depth = _integer(
        solution.get("expansion_depth", 0),
        "solution.expansion_depth",
    )
    if _integer(
        selection.get("selected_solution_id"),
        "main_enzyme_selection.selected_solution_id",
        minimum=1,
    ) != solution_id:
        raise ValueError(
            "主酶组合对应的路径与 manifest 当前路线不一致，"
            "请重新运行主酶选择并写入组合"
        )
    if _integer(
        selection.get("expansion_depth", 0),
        "main_enzyme_selection.expansion_depth",
    ) != expansion_depth:
        raise ValueError(
            "主酶组合对应的扩展深度与 manifest 当前路线不一致，"
            "请重新运行主酶选择并写入组合"
        )

    current_fingerprint = solution_fingerprint(
        solution_id,
        [dict(step) for step in route_steps],
    )
    source = _mapping(
        selection.get("source"),
        "main_enzyme_selection.source",
    )
    stored_fingerprint = _string(source.get("solution_fingerprint"))
    if stored_fingerprint != current_fingerprint:
        raise ValueError(
            "manifest 路线已在主酶组合写入后发生变化，请依次重新运行 "
            "main-enzyme、main-enzyme-sets 和 write --main-enzyme-set N"
        )
    return current_fingerprint


def _normalized_assignments(
    selection: Mapping[str, Any],
    route_by_index: Mapping[int, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    assignments: list[dict[str, Any]] = []
    seen_steps: set[int] = set()
    for raw in _mapping_list(
        selection.get("step_assignments"),
        "main_enzyme_selection.step_assignments",
    ):
        step_index = _integer(
            raw.get("step_index"),
            "step_assignments[].step_index",
            minimum=1,
        )
        if step_index in seen_steps:
            raise ValueError(f"主酶组合对 Step {step_index} 重复分配")
        seen_steps.add(step_index)
        route_step = route_by_index.get(step_index)
        if route_step is None:
            raise ValueError(
                f"主酶组合引用了 manifest 路线中不存在的 Step {step_index}"
            )

        reaction_id = _string(raw.get("reaction_id")).upper()
        route_reaction_id = _string(route_step.get("reaction_id")).upper()
        if not reaction_id or reaction_id != route_reaction_id:
            raise ValueError(
                f"Step {step_index} 的反应已变化：组合记录为 "
                f"{reaction_id or '空'}，路线当前为 "
                f"{route_reaction_id or '空'}"
            )
        accession = _string(raw.get("accession")).upper()
        if not accession:
            raise ValueError(f"Step {step_index} 缺少主酶 accession")
        assignments.append(
            {
                "step_index": step_index,
                "reaction_id": reaction_id,
                "accession": accession,
                "candidate_rank": _integer(
                    raw.get("candidate_rank"),
                    "step_assignments[].candidate_rank",
                    minimum=1,
                ),
            }
        )
    assignments.sort(key=lambda item: item["step_index"])
    return assignments


def _validate_coverage(
    selection: Mapping[str, Any],
    assignments: list[dict[str, Any]],
) -> list[int]:
    coverage = _mapping(
        selection.get("coverage"),
        "main_enzyme_selection.coverage",
    )
    if coverage.get("complete") is not True:
        raise ValueError("主酶组合没有完整覆盖全部必需 Step")
    required = _sorted_unique_step_indexes(
        coverage.get("required_step_indexes"),
        "coverage.required_step_indexes",
    )
    covered = _sorted_unique_step_indexes(
        coverage.get("covered_step_indexes"),
        "coverage.covered_step_indexes",
    )
    assigned = [item["step_index"] for item in assignments]
    if not required or required != covered or required != assigned:
        raise ValueError(
            "主酶组合的 required、covered 和 assigned Step 不一致"
        )
    return required


def _validate_set_fingerprint(
    selection: Mapping[str, Any],
    accessions: list[str],
    assignments: list[dict[str, Any]],
) -> None:
    source = _mapping(
        selection.get("source"),
        "main_enzyme_selection.source",
    )
    candidate_pool_fingerprint = _string(
        source.get("candidate_pool_fingerprint")
    ).lower()
    if not _SHA256_PATTERN.fullmatch(candidate_pool_fingerprint):
        raise ValueError("主酶组合缺少有效的 candidate_pool_fingerprint")
    expected = stable_json_hash(
        {
            "candidate_pool_fingerprint": candidate_pool_fingerprint,
            "accessions": sorted(accessions),
            "assignments": [
                {
                    "step_index": assignment["step_index"],
                    "accession": assignment["accession"],
                    "candidate_rank": assignment["candidate_rank"],
                }
                for assignment in assignments
            ],
        }
    )
    recorded = _string(
        selection.get("selected_set_fingerprint")
    ).lower()
    if recorded != expected:
        raise ValueError(
            "主酶组合指纹校验失败，manifest 可能已损坏或过期；"
            "请重新运行 main-enzyme-sets 并写入组合"
        )


def _assigned_reaction_step(
    step_index: int,
    route_step: Mapping[str, Any],
) -> AssignedReactionStep:
    produced_compound_id = _string(
        route_step.get("produced_compound_id")
    ).upper()
    if not produced_compound_id:
        raise ValueError(
            f"manifest 路线 Step {step_index} 缺少 produced_compound_id"
        )
    return {
        "step_index": step_index,
        "reaction_id": _string(route_step.get("reaction_id")).upper(),
        "reaction_name": _string(route_step.get("reaction_name")),
        "equation": _string(route_step.get("equation")),
        "direction": _string(route_step.get("direction")),
        "precursor_compound_ids": [
            value.upper()
            for value in _split_values(
                route_step.get("precursor_compound_ids")
            )
        ],
        "produced_compound_id": produced_compound_id,
        "produced_compound_name": _string(
            route_step.get("produced_compound_name")
        ),
        "rhea_ids": [
            value.upper()
            for value in _split_values(route_step.get("rhea_ids"))
        ],
    }


def _whole_reaction_context(
    reaction_steps: list[AssignedReactionStep],
) -> WholeReactionContext:
    produced = {
        step["produced_compound_id"]
        for step in reaction_steps
        if step["produced_compound_id"]
    }
    precursors = {
        compound_id
        for step in reaction_steps
        for compound_id in step["precursor_compound_ids"]
    }
    return {
        "equation": None,
        "rhea_ids": [],
        "start_compound_ids": sorted(precursors - produced),
        "end_compound_ids": sorted(produced - precursors),
        "intermediate_compound_ids": sorted(produced & precursors),
        "evidence_status": "unavailable",
        "evidence": [],
    }


def build_main_enzyme_research_units(
    manifest: Mapping[str, Any],
) -> list[MainEnzymeResearchUnit]:
    """Validate one selected manifest and build deterministic research units."""

    if _string(manifest.get("schema_version")) != DESIGN_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            "不支持当前 design manifest 格式，请重新写入路线"
        )
    solution = _mapping(manifest.get("solution"), "solution")
    selection = _mapping(
        manifest.get("main_enzyme_selection"),
        "main_enzyme_selection；请先执行 write --main-enzyme-set N",
    )
    if (
        _string(selection.get("schema_version"))
        != MAIN_ENZYME_MANIFEST_SCHEMA_VERSION
    ):
        raise ValueError(
            "不支持当前主酶选择格式，请重新运行 "
            "write --main-enzyme-set N"
        )
    selection_status = _string(selection.get("selection_status"))
    if selection_status not in _ALLOWED_SELECTION_STATUSES:
        raise ValueError(
            "manifest 中的主酶组合尚未由用户确认，请先执行 "
            "write --main-enzyme-set N"
        )
    _integer(
        selection.get("selected_set_id"),
        "main_enzyme_selection.selected_set_id",
        minimum=1,
    )

    route_steps = _mapping_list(solution.get("steps"), "solution.steps")
    route_by_index = _route_step_index(route_steps)
    _validate_route_lock(solution, selection, route_steps)
    assignments = _normalized_assignments(selection, route_by_index)
    _validate_coverage(selection, assignments)

    assignments_by_accession: dict[str, list[dict[str, Any]]] = {}
    for assignment in assignments:
        assignments_by_accession.setdefault(
            assignment["accession"], []
        ).append(assignment)

    proteins: dict[str, Mapping[str, Any]] = {}
    for protein in _mapping_list(
        selection.get("proteins"),
        "main_enzyme_selection.proteins",
    ):
        accession = _string(protein.get("accession")).upper()
        if not accession:
            raise ValueError("主酶组合包含缺少 accession 的蛋白")
        if accession in proteins:
            raise ValueError(f"主酶组合包含重复 accession：{accession}")
        assigned_steps = _sorted_unique_step_indexes(
            protein.get("assigned_step_indexes"),
            f"{accession}.assigned_step_indexes",
        )
        expected_steps = [
            item["step_index"]
            for item in assignments_by_accession.get(accession, [])
        ]
        if assigned_steps != expected_steps:
            raise ValueError(
                f"主酶 {accession} 的 assigned_step_indexes 与 "
                "step_assignments 不一致"
            )
        proteins[accession] = protein

    assignment_accessions = set(assignments_by_accession)
    if set(proteins) != assignment_accessions:
        missing = sorted(assignment_accessions - set(proteins))
        unused = sorted(set(proteins) - assignment_accessions)
        details: list[str] = []
        if missing:
            details.append(f"缺少蛋白记录：{','.join(missing)}")
        if unused:
            details.append(f"存在未分配蛋白：{','.join(unused)}")
        raise ValueError("主酶组合成员与步骤分配不一致；" + "；".join(details))

    _validate_set_fingerprint(
        selection,
        list(proteins),
        assignments,
    )

    units: list[MainEnzymeResearchUnit] = []
    for accession, protein in proteins.items():
        protein_assignments = assignments_by_accession[accession]
        assigned_step_indexes = [
            item["step_index"] for item in protein_assignments
        ]
        reaction_steps = [
            _assigned_reaction_step(
                step_index,
                route_by_index[step_index],
            )
            for step_index in assigned_step_indexes
        ]
        unit: MainEnzymeResearchUnit = {
            "accession": accession,
            "reaction_scope": (
                "single_step"
                if len(reaction_steps) == 1
                else "multi_step"
            ),
            "assigned_step_indexes": assigned_step_indexes,
            "reaction_steps": reaction_steps,
        }
        protein_name = _string(protein.get("protein_name"))
        organism_name = _string(protein.get("organism_name"))
        if protein_name:
            unit["protein_name"] = protein_name
        if organism_name:
            unit["organism_name"] = organism_name
        if len(reaction_steps) > 1:
            unit["whole_reaction"] = _whole_reaction_context(
                reaction_steps
            )
        units.append(unit)

    units.sort(
        key=lambda item: (
            min(item["assigned_step_indexes"]),
            item["accession"],
        )
    )
    return units


def load_main_enzyme_research_units(
    manifest_path: str | Path,
) -> list[MainEnzymeResearchUnit]:
    """Read a design manifest and return one unit per selected main enzyme."""

    path = Path(manifest_path).expanduser()
    if not path.is_file():
        raise FileNotFoundError(
            f"未找到 design manifest，请先写入路线：{path}"
        )
    manifest = read_design_manifest(path)
    return build_main_enzyme_research_units(manifest)


__all__ = [
    "build_main_enzyme_research_units",
    "load_main_enzyme_research_units",
]
