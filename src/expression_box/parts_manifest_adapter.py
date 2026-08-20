"""Build a validated expression-parts context from the current manifest."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio import SeqIO

from src.expression_box.expression_host_context import resolve_expression_host
from src.expression_box.parts_models import (
    ExpressionCds,
    ExpressionPartsCassette,
    ExpressionPartsContext,
)
from src.write_manifest.store import read_design_manifest


_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_DNA_RE = re.compile(r"[ACGT]+")


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"manifest missing object: {field}")
    return value


def _text(value: Any, field: str) -> str:
    normalized = str(value or "").strip()
    if not normalized:
        raise ValueError(f"manifest field must not be empty: {field}")
    return normalized


def _positive_int(value: Any, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"manifest field must be a positive integer: {field}")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"manifest field must be a positive integer: {field}"
        ) from exc
    if parsed < 1:
        raise ValueError(f"manifest field must be a positive integer: {field}")
    return parsed


def _sha(value: Any, field: str) -> str:
    normalized = str(value or "").strip().lower()
    if _SHA256_RE.fullmatch(normalized) is None:
        raise ValueError(f"manifest field is not a SHA-256 digest: {field}")
    return normalized


def _stable_hash(payload: Any) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _safe_project_path(project_root: Path, path_value: Any, field: str) -> Path:
    value = _text(path_value, field)
    path = Path(value)
    if not path.is_absolute():
        path = project_root / path
    path = path.expanduser().resolve()
    if path != project_root and project_root not in path.parents:
        raise ValueError(f"manifest path escapes project output: {field}")
    return path


def _read_cds(path: Path, expected_sha: str, expected_length: int) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"could not read optimized CDS FASTA: {path}") from exc
    if len(records) != 1:
        raise ValueError(
            "optimized CDS FASTA must contain exactly one record: " + str(path)
        )
    sequence = re.sub(r"\s+", "", str(records[0].seq)).upper()
    if _DNA_RE.fullmatch(sequence) is None:
        raise ValueError(f"optimized CDS contains invalid DNA: {path}")
    if len(sequence) != expected_length:
        raise ValueError(f"optimized CDS length does not match manifest: {path}")
    if hashlib.sha256(sequence.encode("utf-8")).hexdigest() != expected_sha:
        raise ValueError(f"optimized CDS sequence hash does not match manifest: {path}")
    return sequence


def _cds_by_accession(
    cds_selection: Mapping[str, Any],
    project_root: Path,
) -> dict[str, ExpressionCds]:
    raw_proteins = cds_selection.get("proteins")
    if not isinstance(raw_proteins, list) or not raw_proteins:
        raise ValueError("cds_selection.proteins must be a non-empty list")
    result: dict[str, ExpressionCds] = {}
    for index, raw in enumerate(raw_proteins):
        field = f"cds_selection.proteins[{index}]"
        item = _mapping(raw, field)
        accession = _text(item.get("accession"), f"{field}.accession").upper()
        if accession in result:
            raise ValueError(f"cds_selection contains duplicate accession: {accession}")
        optimized = _mapping(item.get("optimized_cds"), f"{field}.optimized_cds")
        sequence_sha = _sha(
            optimized.get("sequence_sha256"),
            f"{field}.optimized_cds.sequence_sha256",
        )
        length_nt = _positive_int(
            optimized.get("length_nt"),
            f"{field}.optimized_cds.length_nt",
        )
        path = _safe_project_path(
            project_root,
            optimized.get("path"),
            f"{field}.optimized_cds.path",
        )
        result[accession] = ExpressionCds(
            accession=accession,
            sequence=_read_cds(path, sequence_sha, length_nt),
            sequence_sha256=sequence_sha,
            path=path,
        )
    return result


def _cassettes(
    selection: Mapping[str, Any],
    cds_by_accession: Mapping[str, ExpressionCds],
) -> tuple[ExpressionPartsCassette, ...]:
    raw_cassettes = selection.get("cassettes")
    if not isinstance(raw_cassettes, list) or not raw_cassettes:
        raise ValueError("expression_box_selection.cassettes must be a non-empty list")
    cassettes: list[ExpressionPartsCassette] = []
    assigned: list[str] = []
    for index, raw in enumerate(raw_cassettes, start=1):
        field = f"expression_box_selection.cassettes[{index - 1}]"
        item = _mapping(raw, field)
        cassette_index = _positive_int(item.get("cassette_index"), f"{field}.cassette_index")
        if cassette_index != index:
            raise ValueError("expression_box_selection cassette indexes must be sequential")
        raw_accessions = item.get("protein_accessions")
        if not isinstance(raw_accessions, list) or not raw_accessions:
            raise ValueError(f"{field}.protein_accessions must be a non-empty list")
        accessions = [str(value or "").strip().upper() for value in raw_accessions]
        if any(not value for value in accessions):
            raise ValueError(f"{field}.protein_accessions contains an empty value")
        unknown = sorted(set(accessions) - set(cds_by_accession))
        if unknown:
            raise ValueError(
                "expression-box selection references unknown CDS accessions: "
                + ", ".join(unknown)
            )
        assigned.extend(accessions)
        cassettes.append(
            ExpressionPartsCassette(
                cassette_index=cassette_index,
                cds=tuple(cds_by_accession[value] for value in accessions),
            )
        )
    if len(assigned) != len(set(assigned)):
        raise ValueError("expression-box selection assigns a CDS more than once")
    if set(assigned) != set(cds_by_accession):
        raise ValueError("expression-box selection does not assign every optimized CDS")
    return tuple(cassettes)


def load_expression_parts_context(
    manifest_path: str | Path,
    project_output_path: str | Path,
) -> ExpressionPartsContext:
    path = Path(manifest_path).expanduser().resolve()
    project_root = Path(project_output_path).expanduser().resolve()
    manifest = read_design_manifest(path)
    target = _text(manifest.get("target_compound_id"), "target_compound_id")
    cds_selection = _mapping(manifest.get("cds_selection"), "cds_selection")
    if cds_selection.get("status") != "complete":
        raise ValueError("cds_selection must be complete before expression-parts design")
    selection = _mapping(
        manifest.get("expression_box_selection"),
        "expression_box_selection",
    )
    if selection.get("schema_version") != "expression_box_selection.v1":
        raise ValueError("unsupported expression_box_selection schema_version")
    if selection.get("selection_status") != "user_selected":
        raise ValueError("expression-box design must be user selected")
    selection_fingerprint = _sha(
        selection.get("selected_design_fingerprint"),
        "expression_box_selection.selected_design_fingerprint",
    )
    cds_source = _text(
        cds_selection.get("source_fingerprint"),
        "cds_selection.source_fingerprint",
    )
    selection_source = _mapping(
        selection.get("source"),
        "expression_box_selection.source",
    )
    if selection_source.get("cds_selection_source_fingerprint") != cds_source:
        raise ValueError(
            "expression-box selection does not match the current CDS selection; "
            "rerun expression --design --box and write --expression-box"
        )
    cds_map = _cds_by_accession(cds_selection, project_root)
    cassettes = _cassettes(selection, cds_map)
    host = resolve_expression_host(dict(manifest))
    try:
        revision = int(manifest.get("revision", 0))
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest revision must be an integer") from exc
    fingerprint_payload = {
        "target_compound_id": target,
        "expression_box_selection_fingerprint": selection_fingerprint,
        "cds_selection_source_fingerprint": cds_source,
        "host": host.to_dict(),
        "cassettes": [
            {
                "cassette_index": cassette.cassette_index,
                "cds": [
                    {
                        "accession": cds.accession,
                        "sequence_sha256": cds.sequence_sha256,
                    }
                    for cds in cassette.cds
                ],
            }
            for cassette in cassettes
        ],
    }
    return ExpressionPartsContext(
        manifest_path=path,
        manifest_revision=revision,
        project_root=project_root,
        target_compound_id=target,
        expression_box_selection_fingerprint=selection_fingerprint,
        cds_selection_source_fingerprint=cds_source,
        host_name=host.codon_transformer_name,
        host_key=host.expression_host_key,
        host_labels=host.milvus_host_labels,
        input_fingerprint=_stable_hash(fingerprint_payload),
        cassettes=cassettes,
    )


__all__ = ["load_expression_parts_context"]
