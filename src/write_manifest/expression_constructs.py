"""Build and transactionally install complete expression-construct GenBank files."""

from __future__ import annotations

import hashlib
import shutil
import tempfile
import uuid
from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib.metadata import version
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from dnachisel import AvoidPattern, DnaOptimizationProblem, EnforceGCContent

from src.expression_box.parts_models import ExpressionPartsContext
from src.protein_to_cds.sequence_constraints import (
    DEFAULT_FORBIDDEN_MOTIFS,
    DNA_ALPHABET,
    HOMOPOLYMER_LIMIT,
    LOCAL_GC_MAX,
    LOCAL_GC_MIN,
    LOCAL_GC_WINDOW_NT,
    gc_fraction,
    local_gc_values,
    max_homopolymer_length,
    motif_hits,
    reverse_complement,
)


ASSEMBLED_EXPRESSION_CONSTRUCTS_SCHEMA_VERSION = (
    "assembled_expression_constructs.v1"
)
EXPRESSION_CONSTRUCTS_DIRNAME = "expression_constructs"
JUNCTION_POLICY = "direct_concatenation"


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_expression_construct_sequence(sequence: str) -> dict[str, Any]:
    """Apply the same hard sequence gates used by expression-parts design."""

    normalized = str(sequence or "").strip().upper()
    valid_alphabet = bool(normalized) and set(normalized).issubset(DNA_ALPHABET)
    if not valid_alphabet:
        return {
            "engine": "DNA Chisel",
            "engine_version": version("dnachisel"),
            "gate_status": "FAIL",
            "checks": {"valid_alphabet": False},
            "failed_checks": ["valid_alphabet"],
        }

    constrained_motifs: set[str] = set()
    for motif in DEFAULT_FORBIDDEN_MOTIFS.values():
        constrained_motifs.add(motif)
        constrained_motifs.add(reverse_complement(motif))
    constraints: list[Any] = [
        EnforceGCContent(mini=0.30, maxi=0.70),
        EnforceGCContent(
            mini=LOCAL_GC_MIN,
            maxi=LOCAL_GC_MAX,
            window=LOCAL_GC_WINDOW_NT,
        ),
        *(AvoidPattern(motif) for motif in sorted(constrained_motifs)),
        *(AvoidPattern(base * HOMOPOLYMER_LIMIT) for base in "ACGT"),
    ]
    problem = DnaOptimizationProblem(
        normalized,
        constraints=constraints,
        objectives=[],
        logger=None,
    )
    local_values = local_gc_values(normalized)
    hits = motif_hits(normalized, DEFAULT_FORBIDDEN_MOTIFS)
    checks = {
        "valid_alphabet": True,
        "global_gc_pass": 0.30 <= gc_fraction(normalized) <= 0.70,
        "local_gc_pass": all(
            LOCAL_GC_MIN <= value <= LOCAL_GC_MAX for value in local_values
        ),
        "forbidden_motif_pass": not hits,
        "homopolymer_pass": (
            max_homopolymer_length(normalized) < HOMOPOLYMER_LIMIT
        ),
        "dnachisel_constraints_pass": problem.all_constraints_pass(),
    }
    return {
        "engine": "DNA Chisel",
        "engine_version": version("dnachisel"),
        "gate_status": "PASS" if all(checks.values()) else "FAIL",
        "sequence_sha256": _sha256_text(normalized),
        "length_nt": len(normalized),
        "gc_percent": round(100.0 * gc_fraction(normalized), 8),
        "local_gc_min_percent": (
            round(100.0 * min(local_values), 8) if local_values else None
        ),
        "local_gc_max_percent": (
            round(100.0 * max(local_values), 8) if local_values else None
        ),
        "forbidden_site_hits": hits,
        "max_homopolymer": max_homopolymer_length(normalized),
        "checks": checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
    }


def _part_sequence(part: Mapping[str, Any], role: str) -> str:
    actual_role = str(part.get("role") or "").strip().lower()
    sequence = str(part.get("sequence") or "").strip().upper()
    if actual_role != role or not sequence or set(sequence) - DNA_ALPHABET:
        raise ValueError(f"invalid {role} sequence in expression-parts design")
    if part.get("sequence_sha256") != _sha256_text(sequence):
        raise ValueError(f"{role} sequence hash does not match expression-parts design")
    if int(part.get("length_bp") or 0) != len(sequence):
        raise ValueError(f"{role} sequence length does not match expression-parts design")
    return sequence


def _qualifiers(**values: Any) -> dict[str, list[str]]:
    return {
        key: [str(value)]
        for key, value in values.items()
        if value is not None and str(value) != ""
    }


def _part_qualifiers(
    part: Mapping[str, Any],
    *,
    label: str,
    parts_design_id: int,
    cassette_index: int,
) -> dict[str, list[str]]:
    return _qualifiers(
        label=label,
        part_id=part.get("part_id"),
        role=part.get("role"),
        strength=part.get("strength"),
        regulation=part.get("regulation"),
        direction=part.get("direction"),
        host_match_kind=part.get("host_match_kind"),
        evidence_grade=part.get("evidence_grade"),
        sequence_sha256=part.get("sequence_sha256"),
        parts_design_id=parts_design_id,
        cassette_index=cassette_index,
    )


def _append_feature(
    features: list[SeqFeature],
    *,
    feature_type: str,
    start: int,
    end: int,
    qualifiers: Mapping[str, list[str]],
) -> None:
    features.append(
        SeqFeature(
            FeatureLocation(start, end, strand=1),
            type=feature_type,
            qualifiers=dict(qualifiers),
        )
    )


def _expected_cassette_sequence(
    raw_cassette: Mapping[str, Any],
    sequence: str,
    design_id: int,
    cassette_index: int,
) -> None:
    assembled = raw_cassette.get("assembled_sequence")
    if not isinstance(assembled, Mapping):
        raise ValueError(
            f"expression-parts design {design_id} cassette {cassette_index} "
            "is missing assembled_sequence"
        )
    expected_length = int(assembled.get("length_nt") or 0)
    expected_hash = str(assembled.get("sequence_sha256") or "")
    if expected_length != len(sequence) or expected_hash != _sha256_text(sequence):
        raise ValueError(
            f"expression-parts design {design_id} cassette {cassette_index} "
            "cannot be reconstructed from the selected parts and optimized CDS"
        )
    audit = raw_cassette.get("sequence_audit")
    if (
        not isinstance(audit, Mapping)
        or audit.get("gate_status") != "PASS"
        or audit.get("sequence_sha256") != expected_hash
    ):
        raise ValueError(
            f"expression-parts design {design_id} cassette {cassette_index} "
            "does not have a matching PASS sequence audit"
        )


def _build_record(
    design: Mapping[str, Any],
    context: ExpressionPartsContext,
    selection_fingerprint: str,
) -> tuple[SeqRecord, dict[str, Any]]:
    design_id = int(design["design_id"])
    rank = int(design["rank"])
    score = float(design["expression_success_score"])
    raw_cassettes = design.get("cassettes")
    if not isinstance(raw_cassettes, list):
        raise ValueError(f"expression-parts design {design_id} is missing cassettes")
    by_index = {
        int(item.get("cassette_index") or 0): item
        for item in raw_cassettes
        if isinstance(item, Mapping)
    }

    sequence_chunks: list[str] = []
    features: list[SeqFeature] = []
    cassette_ranges: list[dict[str, Any]] = []
    component_count = 0
    cursor = 0

    for cassette in sorted(context.cassettes, key=lambda item: item.cassette_index):
        cassette_index = cassette.cassette_index
        raw = by_index.get(cassette_index)
        if raw is None:
            raise ValueError(
                f"expression-parts design {design_id} is missing cassette {cassette_index}"
            )
        promoter = raw.get("promoter")
        terminator = raw.get("terminator")
        genes = raw.get("genes")
        if not isinstance(promoter, Mapping) or not isinstance(terminator, Mapping):
            raise ValueError(f"expression-parts design {design_id} has invalid parts")
        if not isinstance(genes, list) or len(genes) != len(cassette.cds):
            raise ValueError(f"expression-parts design {design_id} has invalid genes")

        cassette_start = cursor
        promoter_sequence = _part_sequence(promoter, "promoter")
        sequence_chunks.append(promoter_sequence)
        _append_feature(
            features,
            feature_type="promoter",
            start=cursor,
            end=cursor + len(promoter_sequence),
            qualifiers=_part_qualifiers(
                promoter,
                label=f"cassette_{cassette_index}_promoter",
                parts_design_id=design_id,
                cassette_index=cassette_index,
            ),
        )
        cursor += len(promoter_sequence)
        component_count += 1

        cassette_sequence_chunks = [promoter_sequence]
        for gene_number, (raw_gene, cds) in enumerate(
            zip(genes, cassette.cds, strict=True),
            start=1,
        ):
            if not isinstance(raw_gene, Mapping):
                raise ValueError(f"expression-parts design {design_id} has invalid gene")
            if str(raw_gene.get("accession") or "") != cds.accession:
                raise ValueError(f"expression-parts design {design_id} CDS order changed")
            if raw_gene.get("cds_sequence_sha256") != cds.sequence_sha256:
                raise ValueError(f"expression-parts design {design_id} CDS hash changed")

            rbs = raw_gene.get("rbs")
            ostir = raw_gene.get("ostir")
            if not isinstance(rbs, Mapping) or not isinstance(ostir, Mapping):
                raise ValueError(f"expression-parts design {design_id} has invalid RBS")
            rbs_sequence = _part_sequence(rbs, "rbs")
            sequence_chunks.append(rbs_sequence)
            cassette_sequence_chunks.append(rbs_sequence)
            rbs_qualifiers = _part_qualifiers(
                rbs,
                label=f"cassette_{cassette_index}_{cds.accession}_RBS",
                parts_design_id=design_id,
                cassette_index=cassette_index,
            )
            rbs_qualifiers.update(
                _qualifiers(
                    target_accession=cds.accession,
                    initiation_rate=ostir.get(
                        "translation_initiation_rate"
                    ),
                    d_g_total=ostir.get("d_g_total"),
                    start_position=ostir.get("intended_start_position"),
                    off_target_starts=ostir.get("unintended_start_count"),
                    ostir_context_sha256=ostir.get("context_sha256"),
                )
            )
            _append_feature(
                features,
                feature_type="RBS",
                start=cursor,
                end=cursor + len(rbs_sequence),
                qualifiers=rbs_qualifiers,
            )
            cursor += len(rbs_sequence)
            component_count += 1

            sequence_chunks.append(cds.sequence)
            cassette_sequence_chunks.append(cds.sequence)
            _append_feature(
                features,
                feature_type="CDS",
                start=cursor,
                end=cursor + len(cds.sequence),
                qualifiers=_qualifiers(
                    label=f"cassette_{cassette_index}_{cds.accession}_CDS",
                    accession=cds.accession,
                    gene_order=gene_number,
                    cds_sequence_sha256=cds.sequence_sha256,
                    parts_design_id=design_id,
                    cassette_index=cassette_index,
                    codon_start=1,
                    transl_table=11,
                ),
            )
            cursor += len(cds.sequence)
            component_count += 1

        terminator_sequence = _part_sequence(terminator, "terminator")
        sequence_chunks.append(terminator_sequence)
        cassette_sequence_chunks.append(terminator_sequence)
        _append_feature(
            features,
            feature_type="terminator",
            start=cursor,
            end=cursor + len(terminator_sequence),
            qualifiers=_part_qualifiers(
                terminator,
                label=f"cassette_{cassette_index}_terminator",
                parts_design_id=design_id,
                cassette_index=cassette_index,
            ),
        )
        cursor += len(terminator_sequence)
        component_count += 1

        cassette_sequence = "".join(cassette_sequence_chunks)
        _expected_cassette_sequence(
            raw,
            cassette_sequence,
            design_id,
            cassette_index,
        )
        cassette_ranges.append(
            {
                "cassette_index": cassette_index,
                "start_1based": cassette_start + 1,
                "end_1based": cursor,
                "length_bp": cursor - cassette_start,
                "protein_accessions": [item.accession for item in cassette.cds],
                "sequence_sha256": _sha256_text(cassette_sequence),
            }
        )

    whole_sequence = "".join(sequence_chunks)
    audit = audit_expression_construct_sequence(whole_sequence)
    if audit.get("gate_status") != "PASS":
        failed = ", ".join(str(item) for item in audit.get("failed_checks", []))
        raise ValueError(
            f"expression-parts design {design_id} complete construct failed "
            f"sequence safety checks: {failed or 'unknown'}"
        )

    record_id = f"{context.target_compound_id}_D{design_id:03d}"[:16]
    record = SeqRecord(
        Seq(whole_sequence),
        id=record_id,
        name=record_id,
        description=(
            f"complete concatenated expression construct for parts design {design_id}"
        ),
    )
    record.annotations = {
        "molecule_type": "DNA",
        "topology": "linear",
        "data_file_division": "SYN",
        "date": "01-JAN-1980",
        "source": "synthetic DNA construct",
        "organism": "synthetic DNA construct",
        "keywords": ["synthetic biology", "expression construct"],
        "comment": (
            f"parts_design_id={design_id}; rank={rank}; "
            f"expression_success_score={score}; junction_policy={JUNCTION_POLICY}; "
            "components are concatenated without linker sequence"
        ),
    }
    source_feature = SeqFeature(
        FeatureLocation(0, len(whole_sequence), strand=1),
        type="source",
        qualifiers=_qualifiers(
            organism="synthetic DNA construct",
            mol_type="other DNA",
            target_compound_id=context.target_compound_id,
            parts_design_id=design_id,
            expression_score=score,
            selection_sha256=selection_fingerprint,
            junction_policy=JUNCTION_POLICY,
        ),
    )
    cassette_features = [
        SeqFeature(
            FeatureLocation(
                int(item["start_1based"]) - 1,
                int(item["end_1based"]),
                strand=1,
            ),
            type="misc_feature",
            qualifiers=_qualifiers(
                label=f"expression_cassette_{item['cassette_index']}",
                cassette_index=item["cassette_index"],
                parts_design_id=design_id,
                protein_accessions=",".join(item["protein_accessions"]),
                sequence_sha256=item["sequence_sha256"],
            ),
        )
        for item in cassette_ranges
    ]
    record.features = [source_feature, *cassette_features, *features]

    metadata = {
        "parts_design_id": design_id,
        "rank": rank,
        "expression_success_score": score,
        "expression_regime": str(design.get("expression_regime") or ""),
        "record_id": record_id,
        "length_bp": len(whole_sequence),
        "sequence_sha256": audit["sequence_sha256"],
        "cassette_count": len(cassette_ranges),
        "component_count": component_count,
        "cassette_ranges": cassette_ranges,
        "sequence_audit": audit,
    }
    return record, metadata


def _write_and_validate_record(
    path: Path,
    record: SeqRecord,
    metadata: Mapping[str, Any],
) -> str:
    written = SeqIO.write(record, path, "genbank")
    if written != 1:
        raise ValueError(f"failed to write one GenBank record: {path}")
    try:
        parsed = SeqIO.read(path, "genbank")
    except Exception as exc:
        raise ValueError(f"could not read back generated GenBank: {path}") from exc
    if parsed.id != record.id or str(parsed.seq).upper() != str(record.seq).upper():
        raise ValueError(f"generated GenBank sequence changed during round trip: {path}")
    if _sha256_text(str(parsed.seq).upper()) != metadata["sequence_sha256"]:
        raise ValueError(f"generated GenBank sequence hash is invalid: {path}")
    if len(parsed.features) != len(record.features):
        raise ValueError(f"generated GenBank feature table is incomplete: {path}")
    expected_features = [
        (item.type, int(item.location.start), int(item.location.end))
        for item in record.features
    ]
    parsed_features = [
        (item.type, int(item.location.start), int(item.location.end))
        for item in parsed.features
    ]
    if parsed_features != expected_features:
        raise ValueError(f"generated GenBank feature coordinates changed: {path}")
    return _sha256_file(path)


def _is_safe_child(path: Path, root: Path) -> bool:
    resolved = path.resolve()
    return resolved != root and root in resolved.parents


def _remove_transaction_path(path: Path, root: Path, allowed_prefix: str) -> None:
    if not path.exists():
        return
    if not _is_safe_child(path, root) or not path.name.startswith(allowed_prefix):
        raise ValueError(f"refusing to remove unsafe construct transaction path: {path}")
    shutil.rmtree(path)


def _current_artifacts_match(
    current: Any,
    desired: Mapping[str, Any],
    target_dir: Path,
) -> bool:
    if not isinstance(current, Mapping) or dict(current) != dict(desired):
        return False
    expected_names = {
        Path(str(item["path"])).name for item in desired.get("constructs", [])
    }
    if not target_dir.is_dir():
        return False
    entries = list(target_dir.iterdir())
    actual_names = {path.name for path in entries}
    if actual_names != expected_names:
        return False
    if any(not path.is_file() for path in entries):
        return False
    for item in desired.get("constructs", []):
        path = target_dir / Path(str(item["path"])).name
        if not path.is_file() or _sha256_file(path) != item["file_sha256"]:
            return False
        try:
            parsed = SeqIO.read(path, "genbank")
        except Exception:
            return False
        if (
            len(parsed.seq) != int(item["length_bp"])
            or _sha256_text(str(parsed.seq).upper()) != item["sequence_sha256"]
        ):
            return False
    return True


@dataclass(slots=True)
class ExpressionConstructTransaction:
    """Staged directory replacement that can be rolled back after manifest failure."""

    project_root: Path
    target_dir: Path
    staging_dir: Path
    section: dict[str, Any]
    needs_install: bool
    is_repair: bool
    backup_dir: Path | None = field(default=None, init=False)
    installed: bool = field(default=False, init=False)

    def install(self) -> None:
        if not self.needs_install or self.installed:
            return
        if not _is_safe_child(self.target_dir, self.project_root):
            raise ValueError(f"unsafe expression construct target: {self.target_dir}")
        if self.target_dir.name != EXPRESSION_CONSTRUCTS_DIRNAME:
            raise ValueError(f"unexpected expression construct target: {self.target_dir}")
        if self.target_dir.exists():
            self.backup_dir = self.project_root / (
                f".{EXPRESSION_CONSTRUCTS_DIRNAME}.backup-{uuid.uuid4().hex}"
            )
            self.target_dir.rename(self.backup_dir)
        try:
            self.staging_dir.rename(self.target_dir)
        except Exception:
            if self.backup_dir is not None and self.backup_dir.exists():
                self.backup_dir.rename(self.target_dir)
            raise
        self.installed = True

    def rollback(self) -> None:
        if not self.installed:
            self.cleanup_staging()
            return
        failed_dir = self.project_root / (
            f".{EXPRESSION_CONSTRUCTS_DIRNAME}.failed-{uuid.uuid4().hex}"
        )
        if self.target_dir.exists():
            self.target_dir.rename(failed_dir)
        if self.backup_dir is not None and self.backup_dir.exists():
            self.backup_dir.rename(self.target_dir)
        _remove_transaction_path(
            failed_dir,
            self.project_root,
            f".{EXPRESSION_CONSTRUCTS_DIRNAME}.failed-",
        )
        self.installed = False
        self.backup_dir = None
        self.cleanup_staging()

    def finalize(self) -> str | None:
        """Remove backups after commit; return a non-fatal cleanup warning."""

        warning: str | None = None
        if self.backup_dir is not None and self.backup_dir.exists():
            try:
                _remove_transaction_path(
                    self.backup_dir,
                    self.project_root,
                    f".{EXPRESSION_CONSTRUCTS_DIRNAME}.backup-",
                )
            except OSError as exc:
                warning = f"could not remove old expression-construct backup: {exc}"
        self.backup_dir = None
        self.cleanup_staging()
        return warning

    def cleanup_staging(self) -> None:
        if self.staging_dir.exists():
            _remove_transaction_path(
                self.staging_dir,
                self.project_root,
                f".{EXPRESSION_CONSTRUCTS_DIRNAME}.staging-",
            )


def prepare_expression_constructs(
    *,
    selected_designs: list[dict[str, Any]],
    selection_payload: Mapping[str, Any],
    context: ExpressionPartsContext,
    current_section: Any,
) -> ExpressionConstructTransaction:
    """Build every selected design in staging and describe the desired manifest section."""

    root = context.project_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    target_dir = root / EXPRESSION_CONSTRUCTS_DIRNAME
    staging_dir = Path(
        tempfile.mkdtemp(
            prefix=f".{EXPRESSION_CONSTRUCTS_DIRNAME}.staging-",
            dir=root,
        )
    ).resolve()
    try:
        constructs: list[dict[str, Any]] = []
        for design in selected_designs:
            record, metadata = _build_record(
                design,
                context,
                str(selection_payload["selection_fingerprint"]),
            )
            filename = f"design_{int(design['design_id']):03d}.gb"
            staged_path = staging_dir / filename
            file_hash = _write_and_validate_record(staged_path, record, metadata)
            constructs.append(
                {
                    **metadata,
                    "format": "genbank",
                    "path": f"{EXPRESSION_CONSTRUCTS_DIRNAME}/{filename}",
                    "file_sha256": file_hash,
                }
            )
        section = {
            "schema_version": ASSEMBLED_EXPRESSION_CONSTRUCTS_SCHEMA_VERSION,
            "status": "assembled",
            "source": "write_expression_parts_selection",
            "source_parts_selection_fingerprint": selection_payload[
                "selection_fingerprint"
            ],
            "target_compound_id": context.target_compound_id,
            "junction_policy": JUNCTION_POLICY,
            "linker_sequence": "",
            "output_dir": EXPRESSION_CONSTRUCTS_DIRNAME,
            "design_count": len(constructs),
            "constructs": constructs,
        }
        artifacts_match = _current_artifacts_match(
            current_section,
            section,
            target_dir,
        )
        is_repair = (
            isinstance(current_section, Mapping)
            and dict(current_section) == section
            and not artifacts_match
        )
        transaction = ExpressionConstructTransaction(
            project_root=root,
            target_dir=target_dir,
            staging_dir=staging_dir,
            section=section,
            needs_install=not artifacts_match,
            is_repair=is_repair,
        )
        if artifacts_match:
            transaction.cleanup_staging()
        return transaction
    except Exception:
        _remove_transaction_path(
            staging_dir,
            root,
            f".{EXPRESSION_CONSTRUCTS_DIRNAME}.staging-",
        )
        raise


__all__ = [
    "ASSEMBLED_EXPRESSION_CONSTRUCTS_SCHEMA_VERSION",
    "EXPRESSION_CONSTRUCTS_DIRNAME",
    "ExpressionConstructTransaction",
    "JUNCTION_POLICY",
    "audit_expression_construct_sequence",
    "prepare_expression_constructs",
]
