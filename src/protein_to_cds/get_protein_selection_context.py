"""Adapt the current design manifest into proteins requiring CDS design."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from src.write_manifest.store import read_design_manifest
from src.main_protein_selection.genome_ec_filter import (
    UNSPECIFIED_GENOME, read_current_genome_selection, validate_genome_ec_filter,
)

@dataclass(frozen=True, slots=True)
class SelectedProteinForCds:
    accession: str
    roles: tuple[str, ...]
    protein_name: str
    organism_name: str
    assigned_step_indexes: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ProteinToCdsContext:
    manifest_path: Path
    manifest_revision: int
    target_compound_id: str
    chassis_key: str
    selected_solution_id: int
    selected_set_id: int
    source_fingerprint: str
    proteins: tuple[SelectedProteinForCds, ...]
    warnings: tuple[str, ...]


def _mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest missing object: {field_name}")
    return value


def _nonempty_text(value: Any, field_name: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"manifest field must not be empty: {field_name}")
    return text


def _positive_int(value: Any, field_name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"manifest field must be an integer: {field_name}")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str) and re.fullmatch(r"[0-9]+", value.strip()):
        result = int(value.strip())
    else:
        raise ValueError(f"manifest field must be an integer: {field_name}")
    if result < 1:
        raise ValueError(f"manifest field must be positive: {field_name}")
    return result


def _step_indexes(value: Any, field_name: str) -> tuple[int, ...]:
    if not isinstance(value, list):
        raise ValueError(f"manifest field must be a list: {field_name}")
    indexes = tuple(sorted({_positive_int(item, field_name) for item in value}))
    if not indexes:
        raise ValueError(f"manifest field must not be empty: {field_name}")
    return indexes


def _stable_fingerprint(payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _main_proteins(
    selection: Mapping[str, Any],
) -> dict[str, SelectedProteinForCds]:
    raw_proteins = selection.get("proteins")
    if not isinstance(raw_proteins, list) or not raw_proteins:
        raise ValueError("main_enzyme_selection.proteins must be a non-empty list")

    proteins: dict[str, SelectedProteinForCds] = {}
    for index, raw in enumerate(raw_proteins):
        item = _mapping(raw, f"main_enzyme_selection.proteins[{index}]")
        accession = _nonempty_text(
            item.get("accession"),
            f"main_enzyme_selection.proteins[{index}].accession",
        ).upper()
        if accession in proteins:
            raise ValueError(f"duplicate main-enzyme accession: {accession}")
        proteins[accession] = SelectedProteinForCds(
            accession=accession,
            roles=("main_enzyme",),
            protein_name=str(item.get("protein_name") or "").strip(),
            organism_name=str(item.get("organism_name") or "").strip(),
            assigned_step_indexes=_step_indexes(
                item.get("assigned_step_indexes"),
                f"main_enzyme_selection.proteins[{index}].assigned_step_indexes",
            ),
        )
    return proteins


def get_proteins_for_cds(
    manifest_path: str | Path, *, genome_path: Any = UNSPECIFIED_GENOME,
) -> ProteinToCdsContext:
    """Return the unique manifest-selected proteins that require CDS design."""

    from src.write_manifest.main_enzyme import (
        MAIN_ENZYME_MANIFEST_SCHEMA_VERSION,
    )

    path = Path(manifest_path).expanduser().resolve()
    manifest = read_design_manifest(path)
    target_compound_id = _nonempty_text(
        manifest.get("target_compound_id"), "target_compound_id"
    )
    try:
        revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision must be an integer") from exc

    selection = _mapping(manifest.get("main_enzyme_selection"), "main_enzyme_selection")
    if selection.get("schema_version") != MAIN_ENZYME_MANIFEST_SCHEMA_VERSION:
        raise ValueError("主酶选择结果已过期，请重新运行 main-enzyme、main-enzyme-sets 和 write --main-enzyme-set")
    selection_status = _nonempty_text(
        selection.get("selection_status"),
        "main_enzyme_selection.selection_status",
    )
    if selection_status not in {"user_selected", "user_selected_pending_review"}:
        raise ValueError("main enzyme set is not user selected: " + selection_status)
    coverage = _mapping(selection.get("coverage"), "main_enzyme_selection.coverage")
    if coverage.get("complete") is not True:
        raise ValueError("main enzyme selection does not completely cover the route")

    chassis_key = _nonempty_text(
        selection.get("chassis_key"), "main_enzyme_selection.chassis_key"
    )
    selected_solution_id = _positive_int(
        selection.get("selected_solution_id"),
        "main_enzyme_selection.selected_solution_id",
    )
    selected_set_id = _positive_int(
        selection.get("selected_set_id"),
        "main_enzyme_selection.selected_set_id",
    )
    warnings: list[str] = []
    if selection_status == "user_selected_pending_review":
        warnings.append(
            "main enzyme selection contains unresolved review items; CDS design continued"
        )

    proteins = _main_proteins(selection)
    from src.main_protein_selection.common import get_solution_steps, heterologous_requirements
    from src.main_protein_selection.provenance import solution_fingerprint

    solution_id, steps = get_solution_steps(manifest)
    genome_filter = validate_genome_ec_filter(
        selection.get("genome_ec_filter"), heterologous_requirements(steps), chassis_key, genome_path,
    )
    if genome_filter is not None:
        source = _mapping(selection.get("source"), "main_enzyme_selection.source")
        if source.get("solution_fingerprint") != solution_fingerprint(solution_id, steps):
            raise ValueError("Selected route changed; rerun main-enzyme and write --main-enzyme-set")
        assigned = {step for protein in proteins.values() for step in protein.assigned_step_indexes}
        if (assigned != set(genome_filter.external_required_step_indexes)
                or coverage.get("required_step_indexes") != genome_filter.required_step_indexes
                or coverage.get("covered_step_indexes") != genome_filter.required_step_indexes):
            raise ValueError("Genome and external enzyme coverage are inconsistent; rerun main-enzyme")
        current_path = path.parent / "main_protein_selection" / "main_enzyme_selection.json"
        if current_path.is_file():
            current = read_current_genome_selection(path, genome_path)
            from src.main_protein_selection.build_main_enzyme_sets import main_enzyme_selection_fingerprint
            if (current is None or current.genome_ec_filter != genome_filter
                    or source.get("candidate_selection_fingerprint") != main_enzyme_selection_fingerprint(current)):
                raise ValueError("Genome candidate results changed; rerun main-enzyme-sets and write --main-enzyme-set")
    ordered = tuple(proteins[accession] for accession in sorted(proteins))
    fingerprint_payload = {
        "schema_version": "protein_to_cds.input.v3",
        "target_compound_id": target_compound_id,
        "selected_solution_id": selected_solution_id,
        "selected_set_id": selected_set_id,
        "selected_set_fingerprint": str(
            selection.get("selected_set_fingerprint") or ""
        ),
        "chassis_key": chassis_key,
        "proteins": [
            {
                "accession": item.accession,
                "roles": item.roles,
                "protein_name": item.protein_name,
                "organism_name": item.organism_name,
                "assigned_step_indexes": item.assigned_step_indexes,
            }
            for item in ordered
        ],
        **({"genome_ec_filter_fingerprint": genome_filter.fingerprint} if genome_filter else {}),
    }
    return ProteinToCdsContext(
        manifest_path=path,
        manifest_revision=revision,
        target_compound_id=target_compound_id,
        chassis_key=chassis_key,
        selected_solution_id=selected_solution_id,
        selected_set_id=selected_set_id,
        source_fingerprint=_stable_fingerprint(fingerprint_payload),
        proteins=ordered,
        warnings=tuple(warnings),
    )


# Transitional name retained for callers that imported the old module-level
# function.  The return type is now the typed current-manifest context.
get_protein_selection_context = get_proteins_for_cds


__all__ = [
    "ProteinToCdsContext",
    "SelectedProteinForCds",
    "get_protein_selection_context",
    "get_proteins_for_cds",
]
