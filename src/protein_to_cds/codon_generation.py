"""CodonTransformer generation with read-only CDS quality assessment."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from src.protein_to_cds.codon_optimization import (
    CdsOptimizationError,
    CdsOptimizationResult,
    _dna_fasta,
    _model_identity,
    _write_atomic,
    _write_json_atomic,
    normalize_forbidden_motifs,
    predict_cds_sequence,
)
from src.protein_to_cds.config import CDS_CONSTRAINT_CONFIG, HostProfile
from src.protein_to_cds.search_protein_sequence import ProteinSequenceRecord
from src.protein_to_cds.sequence_constraints import assess_generated_cds
from src.protein_to_cds.artifacts import ArtifactTransaction

GENERATION_SCHEMA_VERSION = "protein_to_cds.generation.v3"
GENERATION_MODE = "codon_transformer_only"
GENERATION_VALIDATION_METHOD = "codon_transformer_match_protein"


def generate_protein_cds(
    protein: ProteinSequenceRecord,
    host: HostProfile,
    output_dir: str | Path,
    *,
    device: str = "auto",
    additional_forbidden_motifs: Iterable[str] = (),
) -> CdsOptimizationResult:
    """Generate an unchanged CDS; quality metrics are advisory, not acceptance gates.

    Every successful call resets a separate optimized working file, including
    repeated calls for the same protein.
    """
    if protein.primary_accession != protein.requested_accession:
        raise CdsOptimizationError("UniProt redirected the manifest accession")
    if len(protein.sequence) > 2046:
        raise CdsOptimizationError(
            f"protein {protein.primary_accession} exceeds the model limit of 2046 aa"
        )
    motifs = normalize_forbidden_motifs(additional_forbidden_motifs)
    model_identity = _model_identity()
    root = Path(output_dir).expanduser().resolve()
    raw_path = root / "raw_cds" / f"{protein.primary_accession}.raw.fasta"
    optimized_path = root / "optimized_cds" / f"{protein.primary_accession}.fasta"
    report_path = root / "reports" / f"{protein.primary_accession}.generation.json"
    try:
        sequence, device_name = predict_cds_sequence(protein, host, device)
        audit = assess_generated_cds(
            sequence, protein.sequence, host.codon_transformer_organism_id, motifs,
            host_taxon_id=host.taxon_id,
            verify_translation=False,
        )
    except Exception as exc:
        failure_path = report_path.with_name(f"{protein.primary_accession}.generation.failed.json") if report_path.exists() else report_path
        _write_json_atomic(failure_path, {
            "schema_version": GENERATION_SCHEMA_VERSION,
            "processing_mode": GENERATION_MODE,
            "status": "FAIL",
            "protein": {"accession": protein.primary_accession},
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
        "encoding_validation": GENERATION_VALIDATION_METHOD,
        "independent_translation_performed": False,
        "protein": {
            "accession": protein.primary_accession,
            "length_aa": protein.length_aa,
        },
        "host": host.to_dict(),
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
    with ArtifactTransaction((raw_path, optimized_path, report_path)) as transaction:
        _write_atomic(raw_path, _dna_fasta(protein.primary_accession, "codon_transformer_raw", sequence))
        report["raw_cds"] = {
            "path": f"protein_to_cds/raw_cds/{protein.primary_accession}.raw.fasta",
            "length_nt": len(sequence),
        }
        report["final"]["sequence_path"] = optimized_path.relative_to(root).as_posix()
        _write_atomic(optimized_path, _dna_fasta(protein.primary_accession, "optimized_cds", sequence))
        _write_json_atomic(report_path, report)
        transaction.commit()
    return CdsOptimizationResult(
        accession=protein.primary_accession,
        host_name=host.host_name,
        codon_transformer_organism_id=host.codon_transformer_organism_id,
        raw_sequence=sequence,
        final_sequence=sequence,
        raw_fasta_path=raw_path,
        optimized_fasta_path=optimized_path,
        report_path=report_path,
        report=report,
    )


__all__ = ["GENERATION_SCHEMA_VERSION", "GENERATION_MODE", "GENERATION_VALIDATION_METHOD", "generate_protein_cds"]
