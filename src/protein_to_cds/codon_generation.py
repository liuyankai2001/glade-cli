"""CodonTransformer-only generation, separate from the retained repair workflow."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from src.protein_to_cds.codon_optimization import (
    CdsOptimizationError,
    CdsOptimizationResult,
    _dna_fasta,
    _model_identity,
    _read_dna_fasta,
    _stable_fingerprint,
    _write_atomic,
    _write_json_atomic,
    normalize_forbidden_motifs,
    predict_cds_sequence,
)
from src.protein_to_cds.config import CDS_CONSTRAINT_CONFIG, HostProfile
from src.protein_to_cds.search_protein_sequence import ProteinSequenceRecord
from src.protein_to_cds.sequence_constraints import CdsConstraintError, assess_generated_cds, sha256_text

GENERATION_SCHEMA_VERSION = "protein_to_cds.generation.v1"
GENERATION_MODE = "codon_transformer_only"


def generate_protein_cds(
    protein: ProteinSequenceRecord,
    host: HostProfile,
    output_dir: str | Path,
    *,
    device: str = "auto",
    additional_forbidden_motifs: Iterable[str] = (),
) -> CdsOptimizationResult:
    """Generate an unchanged CDS; quality metrics are advisory, not acceptance gates.

    The legacy result's ``final``/``optimized`` fields point to the raw output so
    downstream readers can continue using the current CDS-selection schema.
    Existing repair reports and optimized FASTA files are left untouched.
    """
    if protein.primary_accession != protein.requested_accession:
        raise CdsOptimizationError("UniProt redirected the manifest accession")
    if len(protein.sequence) > 2046:
        raise CdsOptimizationError(
            f"protein {protein.primary_accession} exceeds the model limit of 2046 aa"
        )
    motifs = normalize_forbidden_motifs(additional_forbidden_motifs)
    model_identity = _model_identity()
    fingerprint = _stable_fingerprint({
        "schema_version": GENERATION_SCHEMA_VERSION,
        "processing_mode": GENERATION_MODE,
        "protein_sequence_sha256": protein.sequence_sha256,
        "host_organism_id": host.codon_transformer_organism_id,
        "assessment_policy_version": CDS_CONSTRAINT_CONFIG.policy_version,
        "additional_forbidden_motifs": motifs,
        "model": model_identity,
    })
    root = Path(output_dir).expanduser().resolve()
    raw_path = root / "raw_cds" / f"{protein.primary_accession}.raw.fasta"
    report_path = root / "reports" / f"{protein.primary_accession}.generation.json"
    report: dict[str, Any] | None = None
    sequence = ""
    if raw_path.is_file() and report_path.is_file():
        try:
            candidate = json.loads(report_path.read_text(encoding="utf-8"))
            sequence = _read_dna_fasta(raw_path)
            audit = assess_generated_cds(
                sequence, protein.sequence, host.codon_transformer_organism_id, motifs
            )
            if (
                isinstance(candidate, dict)
                and candidate.get("schema_version") == GENERATION_SCHEMA_VERSION
                and candidate.get("processing_mode") == GENERATION_MODE
                and candidate.get("status") == "PASS"
                and candidate.get("constraint_repair_applied") is False
                and candidate.get("quality_checks_enforced") is False
                and candidate.get("input_fingerprint") == fingerprint
                and candidate.get("changes", {}).get("codon_change_count") == 0
                and candidate.get("changes", {}).get("nucleotide_change_count") == 0
                and candidate.get("protein", {}).get("sequence_sha256") == protein.sequence_sha256
                and candidate.get("host", {}).get("codon_transformer_organism_id") == host.codon_transformer_organism_id
                and all(
                    candidate.get(stage, {}).get(key) == value
                    for stage in ("raw", "final") for key, value in audit.items()
                )
            ):
                report = candidate
        except (OSError, ValueError, TypeError, AttributeError, CdsConstraintError):
            pass

    reused = report is not None
    if report is None:
        try:
            sequence, device_name = predict_cds_sequence(protein, host, device)
            _write_atomic(
                raw_path, _dna_fasta(protein.primary_accession, "codon_transformer_raw", sequence)
            )
            audit = assess_generated_cds(
                sequence, protein.sequence, host.codon_transformer_organism_id, motifs
            )
        except Exception as exc:
            _write_json_atomic(report_path, {
                "schema_version": GENERATION_SCHEMA_VERSION,
                "processing_mode": GENERATION_MODE,
                "status": "FAIL",
                "input_fingerprint": fingerprint,
                "protein": {"accession": protein.primary_accession, "sequence_sha256": protein.sequence_sha256},
                "error": f"{type(exc).__name__}: {exc}",
            })
            raise CdsOptimizationError(f"CDS generation failed for {protein.primary_accession}: {exc}") from exc
        metrics = {**audit, "sequence_path": raw_path.relative_to(root).as_posix()}
        report = {
            "schema_version": GENERATION_SCHEMA_VERSION,
            "processing_mode": GENERATION_MODE,
            "status": "PASS",
            "constraint_repair_applied": False,
            "quality_checks_enforced": False,
            "input_fingerprint": fingerprint,
            "protein": {
                "accession": protein.primary_accession,
                "sequence_sha256": protein.sequence_sha256,
                "length_aa": protein.length_aa,
            },
            "host": {
                "name": host.host_name,
                "chassis_key": host.chassis_key,
                "codon_transformer_organism_id": host.codon_transformer_organism_id,
            },
            "model": {
                **model_identity, "deterministic": True,
                "match_protein": True, "device": device_name,
            },
            "assessment_policy_version": CDS_CONSTRAINT_CONFIG.policy_version,
            "additional_forbidden_motifs": list(motifs),
            "raw": metrics,
            "final": dict(metrics),
            "changes": {
                "nucleotide_change_count": 0, "nucleotide_change_fraction": 0.0,
                "codon_change_count": 0, "codon_change_fraction": 0.0,
            },
        }
        _write_json_atomic(report_path, report)
    return CdsOptimizationResult(
        accession=protein.primary_accession,
        protein_sequence_sha256=protein.sequence_sha256,
        input_fingerprint=fingerprint,
        host_name=host.host_name,
        codon_transformer_organism_id=host.codon_transformer_organism_id,
        raw_sequence=sequence,
        final_sequence=sequence,
        raw_sequence_sha256=sha256_text(sequence),
        final_sequence_sha256=sha256_text(sequence),
        raw_fasta_path=raw_path,
        optimized_fasta_path=raw_path,
        report_path=report_path,
        report=report,
        reused_existing=reused,
    )


__all__ = ["GENERATION_SCHEMA_VERSION", "GENERATION_MODE", "generate_protein_cds"]
