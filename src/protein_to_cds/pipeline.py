"""Manifest-driven batch workflow from selected proteins to optimized CDSs."""

from __future__ import annotations

import json
import os
import re
import tempfile
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from src.protein_to_cds.codon_optimization import (
    normalize_forbidden_motifs,
    optimize_protein_cds,
)
from src.protein_to_cds.config import host_profile_for_chassis
from src.protein_to_cds.get_protein_selection_context import get_proteins_for_cds
from src.protein_to_cds.search_protein_sequence import (
    load_uploaded_protein_sequence,
    search_protein_sequence,
)
from src.protein_to_cds.write_protein_to_manifest import (
    CompletedDirectCds,
    CompletedProteinCds,
    FailedProteinCds,
    write_cds_selection_to_manifest,
)

RUN_SUMMARY_SCHEMA_VERSION = "protein_to_cds.run_summary.v2"


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _relative(project_root: Path, path: Path) -> str:
    return path.resolve().relative_to(project_root.resolve()).as_posix()


def _status(success_count: int, failure_count: int) -> str:
    if success_count and not failure_count:
        return "complete"
    if success_count:
        return "partial"
    return "failed"


def _uploaded_sequence_path(
    project_root: Path,
    relative_path: str | None,
) -> Path:
    value = str(relative_path or "").strip()
    if not value:
        raise ValueError("uploaded sequence path is missing from the manifest")
    path = (project_root / value).resolve()
    try:
        path.relative_to(project_root.resolve())
    except ValueError as exc:
        raise ValueError("uploaded sequence path is outside the project output") from exc
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def _read_direct_cds(path: Path) -> str:
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


def run_protein_to_cds(
    config: Any,
    *,
    device: str | None = None,
    additional_forbidden_motifs: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Run the complete selected-protein to optimized-CDS workflow."""

    manifest_path = Path(config.manifest_output_path).expanduser().resolve()
    project_root = Path(config.project_output_path).expanduser().resolve()
    context = get_proteins_for_cds(manifest_path)
    host = host_profile_for_chassis(context.chassis_key)
    requested_device = (
        str(
            (device if device is not None else getattr(config, "device", "auto"))
            or "auto"
        )
        .strip()
        .lower()
    )
    if requested_device not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be one of: auto, cpu, cuda")
    motif_values = (
        additional_forbidden_motifs
        if additional_forbidden_motifs is not None
        else getattr(config, "forbidden_motif", ())
    ) or ()
    motif_items = (
        (motif_values,) if isinstance(motif_values, str) else tuple(motif_values)
    )
    motifs = normalize_forbidden_motifs(motif_items)

    output_root = project_root / "protein_to_cds"
    protein_sequence_dir = output_root / "protein_sequences"
    for directory in (
        protein_sequence_dir,
        output_root / "raw_cds",
        output_root / "optimized_cds",
        output_root / "reports",
    ):
        directory.mkdir(parents=True, exist_ok=True)

    successes: list[CompletedProteinCds] = []
    direct_cds: list[CompletedDirectCds] = []
    failures: list[FailedProteinCds] = []
    for selected in context.proteins:
        try:
            if selected.sequence_type == "cds":
                cds_path = _uploaded_sequence_path(
                    project_root,
                    selected.source_sequence_path,
                )
                direct_cds.append(
                    CompletedDirectCds(
                        selected=selected,
                        cds_path=cds_path,
                        sequence=_read_direct_cds(cds_path),
                    )
                )
                continue
            if selected.source_sequence_path is not None:
                protein = load_uploaded_protein_sequence(
                    selected.accession,
                    _uploaded_sequence_path(
                        project_root,
                        selected.source_sequence_path,
                    ),
                )
            else:
                protein = search_protein_sequence(
                    selected.accession,
                    protein_sequence_dir,
                )
            if protein.primary_accession != selected.accession:
                raise ValueError(
                    "UniProt response accession does not match the manifest: "
                    f"{selected.accession} -> {protein.primary_accession}"
                )
            optimization = optimize_protein_cds(
                protein,
                host,
                output_root,
                device=requested_device,
                additional_forbidden_motifs=motifs,
            )
            successes.append(
                CompletedProteinCds(
                    selected=selected,
                    protein=protein,
                    optimization=optimization,
                )
            )
        except Exception as exc:
            failures.append(
                FailedProteinCds(
                    selected=selected,
                    error_type=type(exc).__name__,
                    message=str(exc),
                )
            )

    success_count = len(successes) + len(direct_cds)
    status = _status(success_count, len(failures))
    summary_path = output_root / "run_summary.json"
    summary: dict[str, Any] = {
        "schema_version": RUN_SUMMARY_SCHEMA_VERSION,
        "status": status,
        "generated_at": _utc_now(),
        "target_compound_id": context.target_compound_id,
        "source_manifest": str(context.manifest_path),
        "source_manifest_revision": context.manifest_revision,
        "source_fingerprint": context.source_fingerprint,
        "host": {
            "chassis_key": host.chassis_key,
            "name": host.host_name,
            "codon_transformer_organism_id": (host.codon_transformer_organism_id),
        },
        "device_request": requested_device,
        "additional_forbidden_motifs": list(motifs),
        "counts": {
            "selected": len(context.proteins),
            "succeeded": success_count,
            "failed": len(failures),
        },
        "successes": [
            {
                "accession": success.selected.accession,
                "roles": list(success.selected.roles),
                "processing_mode": "optimized_from_protein",
                "protein_sequence": _relative(
                    project_root,
                    success.protein.fasta_path,
                ),
                "optimized_cds": _relative(
                    project_root,
                    success.optimization.optimized_fasta_path,
                ),
                "optimization_report": _relative(
                    project_root,
                    success.optimization.report_path,
                ),
                "protein_reused": success.protein.reused_existing,
                "optimization_reused": success.optimization.reused_existing,
            }
            for success in successes
        ]
        + [
            {
                "accession": item.selected.accession,
                "roles": list(item.selected.roles),
                "processing_mode": "user_uploaded_cds",
                "uploaded_cds": _relative(project_root, item.cds_path),
                "optimization_skipped": True,
            }
            for item in direct_cds
        ],
        "failures": [
            {
                "accession": failure.selected.accession,
                "roles": list(failure.selected.roles),
                "error_type": failure.error_type,
                "message": failure.message,
            }
            for failure in failures
        ],
        "warnings": list(context.warnings),
        "manifest_written": False,
    }
    _write_json_atomic(summary_path, summary)

    try:
        manifest_result = write_cds_selection_to_manifest(
            context=context,
            host=host,
            project_output_path=project_root,
            successes=successes,
            failures=failures,
            warnings=context.warnings,
            run_summary_path=summary_path,
            direct_cds=direct_cds,
        )
    except Exception as exc:
        summary["manifest_error"] = {
            "error_type": type(exc).__name__,
            "message": str(exc),
        }
        _write_json_atomic(summary_path, summary)
        raise

    summary["manifest_written"] = True
    summary["manifest_revision"] = manifest_result["manifest_revision"]
    _write_json_atomic(summary_path, summary)
    return {
        "ok": status == "complete",
        "status": status,
        "target_compound_id": context.target_compound_id,
        "counts": summary["counts"],
        "warnings": summary["warnings"],
        "failures": summary["failures"],
        "output_dir": str(output_root),
        "run_summary": str(summary_path),
        "manifest_path": str(context.manifest_path),
        "manifest_revision": manifest_result["manifest_revision"],
    }


__all__ = ["RUN_SUMMARY_SCHEMA_VERSION", "run_protein_to_cds"]
