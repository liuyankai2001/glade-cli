"""Validate a protein-to-CDS batch and commit it to the design manifest."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq

from src.protein_to_cds.codon_optimization import CdsOptimizationResult
from src.protein_to_cds.config import HostProfile
from src.protein_to_cds.get_protein_selection_context import (
    ProteinToCdsContext,
    SelectedProteinForCds,
)
from src.protein_to_cds.search_protein_sequence import ProteinSequenceRecord
from src.write_manifest.store import read_design_manifest, update_design_manifest

CDS_SELECTION_SCHEMA_VERSION = "protein_to_cds.selection.v2"
CDS_SELECTION_DOWNSTREAM_SECTIONS = (
    "expression_box_selection",
    "expression_cassette_assembly",
    "parts_selection",
    "assembled_expression_cassettes",
    "assembled_expression_constructs",
    "plasmid_selection",
    "final_assembly_plan",
    "final_assembly",
)


@dataclass(frozen=True, slots=True)
class CompletedProteinCds:
    """All source and output artifacts for one successful protein."""

    selected: SelectedProteinForCds
    protein: ProteinSequenceRecord
    optimization: CdsOptimizationResult


@dataclass(frozen=True, slots=True)
class CompletedDirectCds:
    """One user-supplied CDS that intentionally skipped optimization."""

    selected: SelectedProteinForCds
    cds_path: Path
    sequence: str


@dataclass(frozen=True, slots=True)
class FailedProteinCds:
    """Serializable failure associated with one selected protein."""

    selected: SelectedProteinForCds
    error_type: str
    message: str


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _read_single_fasta(
    path: Path,
    alphabet: frozenset[str],
) -> tuple[str, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"could not read FASTA: {path}") from exc
    if len(records) != 1:
        raise ValueError(
            f"FASTA must contain exactly one record; found {len(records)}: {path}"
        )
    record = records[0]
    sequence = str(record.seq).upper()
    if not sequence or not set(sequence).issubset(alphabet):
        raise ValueError(f"invalid FASTA sequence: {path}")
    return str(record.id), sequence


def _relative_path(project_output_path: Path, path: Path) -> str:
    resolved_root = project_output_path.resolve()
    resolved_path = path.resolve()
    try:
        return resolved_path.relative_to(resolved_root).as_posix()
    except ValueError as exc:
        raise ValueError(
            f"protein-to-CDS artifact is outside project output: {resolved_path}"
        ) from exc


def _translated_protein(cds: str) -> str:
    if len(cds) % 3:
        raise ValueError("CDS length is not a multiple of three")
    if not cds.startswith("ATG"):
        raise ValueError("CDS does not start with ATG")
    if cds[-3:] not in {"TAA", "TAG", "TGA"}:
        raise ValueError("CDS does not end with a bacterial stop codon")
    translated = str(Seq(cds).translate(table=11, to_stop=False))
    if translated.count("*") != 1 or not translated.endswith("*"):
        raise ValueError("CDS contains an internal stop codon")
    return translated[:-1]


def _validated_success_payload(
    success: CompletedProteinCds,
    *,
    host: HostProfile,
    project_output_path: Path,
) -> dict[str, Any]:
    selected = success.selected
    protein = success.protein
    optimization = success.optimization
    accession = selected.accession

    if (
        protein.requested_accession != accession
        or protein.primary_accession != accession
    ):
        raise ValueError(
            f"protein artifact accession does not match manifest selection: {accession}"
        )
    protein_record_id, protein_sequence = _read_single_fasta(
        protein.fasta_path,
        frozenset("ACDEFGHIKLMNPQRSTVWY"),
    )
    if protein_record_id != accession:
        raise ValueError(f"protein FASTA header accession mismatch: {accession}")
    protein_hash = _sha256_text(protein_sequence)
    if protein_sequence != protein.sequence or protein_hash != protein.sequence_sha256:
        raise ValueError(f"protein FASTA changed after validation: {accession}")

    raw_record_id, raw_sequence = _read_single_fasta(
        optimization.raw_fasta_path,
        frozenset("ACGT"),
    )
    final_record_id, final_sequence = _read_single_fasta(
        optimization.optimized_fasta_path,
        frozenset("ACGT"),
    )
    if raw_record_id != accession or final_record_id != accession:
        raise ValueError(f"CDS FASTA header accession mismatch: {accession}")
    if optimization.accession != accession:
        raise ValueError(f"optimization result accession mismatch: {accession}")
    if optimization.protein_sequence_sha256 != protein_hash:
        raise ValueError(f"optimization result protein hash mismatch: {accession}")
    if optimization.codon_transformer_organism_id != host.codon_transformer_organism_id:
        raise ValueError(f"optimization result host mismatch: {accession}")
    if raw_sequence != optimization.raw_sequence:
        raise ValueError(f"raw CDS FASTA changed after optimization: {accession}")
    if final_sequence != optimization.final_sequence:
        raise ValueError(f"optimized CDS FASTA changed after optimization: {accession}")
    if _sha256_text(raw_sequence) != optimization.raw_sequence_sha256:
        raise ValueError(f"raw CDS hash mismatch: {accession}")
    if _sha256_text(final_sequence) != optimization.final_sequence_sha256:
        raise ValueError(f"optimized CDS hash mismatch: {accession}")
    if _translated_protein(raw_sequence) != protein_sequence:
        raise ValueError(f"raw CDS does not encode the selected protein: {accession}")
    if _translated_protein(final_sequence) != protein_sequence:
        raise ValueError(
            f"optimized CDS does not encode the selected protein: {accession}"
        )

    try:
        report = json.loads(optimization.report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"optimization report is invalid JSON: {accession}") from exc
    if not isinstance(report, dict):
        raise ValueError(f"optimization report root is not an object: {accession}")
    if report != optimization.report:
        raise ValueError(f"optimization report changed after optimization: {accession}")
    if (
        report.get("status") != "PASS"
        or report.get("final", {}).get("gate_status") != "PASS"
    ):
        raise ValueError(f"optimization constraint gate did not pass: {accession}")
    if report.get("protein", {}).get("sequence_sha256") != protein_hash:
        raise ValueError(f"optimization report protein hash mismatch: {accession}")
    if report.get("input_fingerprint") != optimization.input_fingerprint:
        raise ValueError(f"optimization report input fingerprint mismatch: {accession}")
    if report.get("raw", {}).get("sequence_sha256") != optimization.raw_sequence_sha256:
        raise ValueError(f"optimization report raw CDS hash mismatch: {accession}")
    if (
        report.get("final", {}).get("sequence_sha256")
        != optimization.final_sequence_sha256
    ):
        raise ValueError(f"optimization report final CDS hash mismatch: {accession}")
    if (
        report.get("host", {}).get("codon_transformer_organism_id")
        != host.codon_transformer_organism_id
    ):
        raise ValueError(f"optimization report host mismatch: {accession}")

    return {
        "accession": accession,
        "roles": list(selected.roles),
        "protein_name": selected.protein_name,
        "organism_name": selected.organism_name,
        "assigned_step_indexes": list(selected.assigned_step_indexes),
        "required_by_main_accessions": list(selected.required_by_main_accessions),
        "sequence_input": {
            "type": "protein",
            "source": (
                "user_uploaded"
                if selected.source_sequence_path is not None
                else "uniprot"
            ),
        },
        "protein_sequence": {
            "path": _relative_path(project_output_path, protein.fasta_path),
            "file_sha256": _sha256_file(protein.fasta_path),
            "length_aa": len(protein_sequence),
            "source_url": protein.source_url,
            "reused_existing": protein.reused_existing,
        },
        "raw_cds": {
            "path": _relative_path(project_output_path, optimization.raw_fasta_path),
            "file_sha256": _sha256_file(optimization.raw_fasta_path),
            "sequence_sha256": optimization.raw_sequence_sha256,
            "length_nt": len(raw_sequence),
        },
        "optimized_cds": {
            "path": _relative_path(
                project_output_path,
                optimization.optimized_fasta_path,
            ),
            "file_sha256": _sha256_file(optimization.optimized_fasta_path),
            "sequence_sha256": optimization.final_sequence_sha256,
            "length_nt": len(final_sequence),
            "report": {
                "path": _relative_path(project_output_path, optimization.report_path),
                "file_sha256": _sha256_file(optimization.report_path),
            },
            "input_fingerprint": optimization.input_fingerprint,
            "constraint_policy_version": report.get("constraint_policy_version"),
            "metrics": {
                "raw": {
                    key: value
                    for key, value in report.get("raw", {}).items()
                    if key != "sequence_path"
                },
                "final": {
                    key: value
                    for key, value in report.get("final", {}).items()
                    if key != "sequence_path"
                },
                "changes": report.get("changes", {}),
            },
            "reused_existing": optimization.reused_existing,
            "optimization_skipped": False,
        },
    }


def _read_uploaded_cds(path: Path) -> str:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"could not read uploaded CDS: {path}") from exc
    sequence = "".join(
        re.sub(r"\s+", "", line)
        for line in lines
        if line.strip() and not line.lstrip().startswith(">")
    ).upper()
    if not sequence:
        raise ValueError(f"uploaded CDS is empty: {path}")
    return sequence


def _validated_direct_cds_payload(
    direct: CompletedDirectCds,
    *,
    project_output_path: Path,
) -> dict[str, Any]:
    selected = direct.selected
    path = direct.cds_path.resolve()
    sequence = _read_uploaded_cds(path)
    if sequence != direct.sequence:
        raise ValueError(
            f"uploaded CDS changed during protein-to-CDS processing: {selected.accession}"
        )
    return {
        "accession": selected.accession,
        "roles": list(selected.roles),
        "protein_name": selected.protein_name,
        "organism_name": selected.organism_name,
        "assigned_step_indexes": list(selected.assigned_step_indexes),
        "required_by_main_accessions": list(
            selected.required_by_main_accessions
        ),
        "sequence_input": {
            "type": "cds",
            "source": "user_uploaded",
        },
        "optimized_cds": {
            "path": _relative_path(project_output_path, path),
            "file_sha256": _sha256_file(path),
            "sequence_sha256": _sha256_text(sequence),
            "length_nt": len(sequence),
            "source_type": "user_uploaded",
            "optimization_skipped": True,
        },
    }


def _failure_payload(failure: FailedProteinCds) -> dict[str, Any]:
    return {
        "accession": failure.selected.accession,
        "roles": list(failure.selected.roles),
        "assigned_step_indexes": list(failure.selected.assigned_step_indexes),
        "required_by_main_accessions": list(
            failure.selected.required_by_main_accessions
        ),
        "error_type": failure.error_type,
        "message": failure.message,
    }


def write_cds_selection_to_manifest(
    *,
    context: ProteinToCdsContext,
    host: HostProfile,
    project_output_path: str | Path,
    successes: Iterable[CompletedProteinCds],
    failures: Iterable[FailedProteinCds],
    warnings: Iterable[str],
    run_summary_path: str | Path,
    direct_cds: Iterable[CompletedDirectCds] = (),
) -> dict[str, Any]:
    """Validate the full batch and atomically replace ``cds_selection``."""

    project_root = Path(project_output_path).expanduser().resolve()
    summary_path = Path(run_summary_path).expanduser().resolve()
    manifest = read_design_manifest(context.manifest_path)
    try:
        current_revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision must be an integer") from exc
    if current_revision != context.manifest_revision:
        raise ValueError(
            "manifest changed during protein-to-CDS processing: "
            f"expected revision {context.manifest_revision}, current {current_revision}"
        )
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)

    success_items = [
        _validated_success_payload(
            success,
            host=host,
            project_output_path=project_root,
        )
        for success in successes
    ]
    direct_items = [
        _validated_direct_cds_payload(
            item,
            project_output_path=project_root,
        )
        for item in direct_cds
    ]
    success_items.extend(direct_items)
    failure_items = [_failure_payload(failure) for failure in failures]
    selected_accessions = {item.accession for item in context.proteins}
    recorded_accessions = {item["accession"] for item in success_items + failure_items}
    if recorded_accessions != selected_accessions:
        raise ValueError(
            "CDS batch does not account for every selected protein; missing/extra: "
            + ", ".join(sorted(selected_accessions ^ recorded_accessions))
        )
    if len(success_items) + len(failure_items) != len(context.proteins):
        raise ValueError("CDS batch contains duplicate protein results")

    if success_items and not failure_items:
        status = "complete"
    elif success_items:
        status = "partial"
    else:
        status = "failed"
    payload = {
        "schema_version": CDS_SELECTION_SCHEMA_VERSION,
        "status": status,
        "generated_at": _utc_now(),
        "source_manifest_revision": context.manifest_revision,
        "source_fingerprint": context.source_fingerprint,
        "source_selection": {
            "selected_solution_id": context.selected_solution_id,
            "selected_set_id": context.selected_set_id,
        },
        "host": {
            "chassis_key": host.chassis_key,
            "name": host.host_name,
            "codon_transformer_organism_id": (host.codon_transformer_organism_id),
        },
        "counts": {
            "selected": len(context.proteins),
            "succeeded": len(success_items),
            "failed": len(failure_items),
        },
        "proteins": success_items,
        "failures": failure_items,
        "warnings": list(dict.fromkeys(str(value) for value in warnings if str(value))),
        "run_summary": {
            "path": _relative_path(project_root, summary_path),
        },
    }
    updated = update_design_manifest(
        context.manifest_path,
        target_compound_id=context.target_compound_id,
        sections={"cds_selection": payload},
        discard_sections=CDS_SELECTION_DOWNSTREAM_SECTIONS,
        expected_revision=context.manifest_revision,
    )
    return {
        "manifest_path": str(context.manifest_path),
        "manifest_revision": updated["revision"],
        "cds_selection": payload,
    }


__all__ = [
    "CDS_SELECTION_DOWNSTREAM_SECTIONS",
    "CDS_SELECTION_SCHEMA_VERSION",
    "CompletedDirectCds",
    "CompletedProteinCds",
    "FailedProteinCds",
    "write_cds_selection_to_manifest",
]
