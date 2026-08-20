"""Adapt the current design manifest into proteins requiring CDS design."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from src.write_manifest.store import read_design_manifest

_SHA256_RE = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True, slots=True)
class SelectedProteinForCds:
    accession: str
    roles: tuple[str, ...]
    protein_name: str
    organism_name: str
    assigned_step_indexes: tuple[int, ...]
    required_by_main_accessions: tuple[str, ...]
    expected_sequence_sha256: str | None


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


def _sha256_or_none(value: Any, field_name: str) -> str | None:
    if value is None or str(value).strip() == "":
        return None
    normalized = str(value).strip().lower()
    if _SHA256_RE.fullmatch(normalized) is None:
        raise ValueError(f"manifest field is not a SHA-256 digest: {field_name}")
    return normalized


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
            required_by_main_accessions=(),
            expected_sequence_sha256=_sha256_or_none(
                item.get("sequence_sha256"),
                f"main_enzyme_selection.proteins[{index}].sequence_sha256",
            ),
        )
    return proteins


def _auxiliary_details(
    auxiliary: Mapping[str, Any],
    main_proteins: Mapping[str, SelectedProteinForCds],
) -> dict[str, SelectedProteinForCds]:
    raw_introduce = auxiliary.get("auxiliary_proteins_to_introduce")
    if not isinstance(raw_introduce, list):
        raise ValueError(
            "auxiliary_protein_selection.auxiliary_proteins_to_introduce must be a list"
        )
    introduce = {str(value or "").strip().upper() for value in raw_introduce}
    introduce.discard("")
    if not introduce:
        return {}

    raw_main_entries = auxiliary.get("main_enzymes")
    if not isinstance(raw_main_entries, list):
        raise ValueError("auxiliary_protein_selection.main_enzymes must be a list")

    details: dict[str, SelectedProteinForCds] = {}
    for main_index, raw_main in enumerate(raw_main_entries):
        main_entry = _mapping(
            raw_main,
            f"auxiliary_protein_selection.main_enzymes[{main_index}]",
        )
        main_accession = _nonempty_text(
            main_entry.get("accession"),
            f"auxiliary_protein_selection.main_enzymes[{main_index}].accession",
        ).upper()
        main = main_proteins.get(main_accession)
        if main is None:
            raise ValueError(
                "auxiliary protein research references an unselected main enzyme: "
                f"{main_accession}"
            )
        confirmed = main_entry.get("confirmed_auxiliary_proteins")
        if not isinstance(confirmed, list):
            raise ValueError(
                "auxiliary_protein_selection.main_enzymes"
                f"[{main_index}].confirmed_auxiliary_proteins must be a list"
            )
        for aux_index, raw_aux in enumerate(confirmed):
            aux = _mapping(
                raw_aux,
                "auxiliary_protein_selection.main_enzymes"
                f"[{main_index}].confirmed_auxiliary_proteins[{aux_index}]",
            )
            accession = _nonempty_text(
                aux.get("accession"), "auxiliary accession"
            ).upper()
            if accession not in introduce:
                continue
            previous = details.get(accession)
            if previous is None:
                details[accession] = SelectedProteinForCds(
                    accession=accession,
                    roles=("auxiliary_protein",),
                    protein_name=str(aux.get("protein_name") or "").strip(),
                    organism_name=str(aux.get("organism_name") or "").strip(),
                    assigned_step_indexes=main.assigned_step_indexes,
                    required_by_main_accessions=(main_accession,),
                    expected_sequence_sha256=None,
                )
            else:
                details[accession] = replace(
                    previous,
                    assigned_step_indexes=tuple(
                        sorted(
                            set(previous.assigned_step_indexes)
                            | set(main.assigned_step_indexes)
                        )
                    ),
                    required_by_main_accessions=tuple(
                        sorted(
                            set(previous.required_by_main_accessions) | {main_accession}
                        )
                    ),
                )

    missing = sorted(introduce - set(details))
    if missing:
        raise ValueError(
            "auxiliary_proteins_to_introduce lacks confirmed protein details: "
            + ", ".join(missing)
        )
    return details


def get_proteins_for_cds(manifest_path: str | Path) -> ProteinToCdsContext:
    """Return the unique manifest-selected proteins that require CDS design."""

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
    auxiliary_raw = manifest.get("auxiliary_protein_selection")
    auxiliary_source_fingerprint = ""
    if auxiliary_raw is None:
        warnings.append(
            "auxiliary_protein_selection is absent; only selected main enzymes were processed"
        )
    else:
        auxiliary = _mapping(auxiliary_raw, "auxiliary_protein_selection")
        if auxiliary.get("can_advance") is not True:
            raise ValueError("auxiliary protein selection cannot advance to CDS design")
        source = auxiliary.get("source")
        if isinstance(source, Mapping):
            auxiliary_source_fingerprint = str(
                source.get("result_fingerprint") or ""
            ).strip()
        for accession, auxiliary_protein in _auxiliary_details(
            auxiliary, proteins
        ).items():
            if accession in proteins:
                main = proteins[accession]
                proteins[accession] = replace(
                    main,
                    roles=tuple(sorted(set(main.roles) | set(auxiliary_protein.roles))),
                    assigned_step_indexes=tuple(
                        sorted(
                            set(main.assigned_step_indexes)
                            | set(auxiliary_protein.assigned_step_indexes)
                        )
                    ),
                    required_by_main_accessions=auxiliary_protein.required_by_main_accessions,
                )
            else:
                proteins[accession] = auxiliary_protein

    ordered = tuple(proteins[accession] for accession in sorted(proteins))
    fingerprint_payload = {
        "target_compound_id": target_compound_id,
        "selected_solution_id": selected_solution_id,
        "selected_set_id": selected_set_id,
        "selected_set_fingerprint": str(
            selection.get("selected_set_fingerprint") or ""
        ),
        "auxiliary_source_fingerprint": auxiliary_source_fingerprint,
        "chassis_key": chassis_key,
        "proteins": [
            {
                "accession": item.accession,
                "roles": item.roles,
                "protein_name": item.protein_name,
                "organism_name": item.organism_name,
                "assigned_step_indexes": item.assigned_step_indexes,
                "required_by_main_accessions": (item.required_by_main_accessions),
                "expected_sequence_sha256": item.expected_sequence_sha256,
            }
            for item in ordered
        ],
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
