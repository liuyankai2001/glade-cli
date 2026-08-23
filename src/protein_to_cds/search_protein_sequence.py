"""Fetch and persist one validated protein sequence from UniProtKB."""

from __future__ import annotations

import hashlib
import os
import re
import tempfile
import threading
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from Bio import SeqIO

UNIPROT_FASTA_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
DEFAULT_TIMEOUT_SECONDS = 30.0
USER_AGENT = "glade/0.1.0"
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# UniProt accessions use either the historical six-character form or the
# newer ten-character form.  Isoform accessions add a numeric suffix.
_ACCESSION_RE = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:-[1-9][0-9]*)?"
)
_WRITE_LOCK = threading.RLock()


class ProteinSequenceError(Exception):
    """Base error raised by the UniProt protein-sequence layer."""


class InvalidUniProtAccessionError(ProteinSequenceError):
    """Raised when a requested UniProt accession has an invalid shape."""


class UniProtNotFoundError(ProteinSequenceError):
    """Raised when UniProt does not contain the requested accession."""


class UniProtServiceError(ProteinSequenceError):
    """Raised when UniProt cannot provide a usable response."""


class InvalidProteinSequenceError(ProteinSequenceError):
    """Raised when a fetched FASTA record cannot be used for CDS design."""


class ProteinSequenceConflictError(ProteinSequenceError):
    """Raised when an accession would overwrite a different local sequence."""


@dataclass(frozen=True, slots=True)
class ProteinSequenceRecord:
    """One validated and locally persisted UniProt protein sequence."""

    requested_accession: str
    primary_accession: str
    entry_name: str
    sequence: str
    length_aa: int
    sequence_sha256: str
    source_url: str
    fasta_path: Path
    reused_existing: bool


def normalize_uniprot_accession(accession: str) -> str:
    """Return a normalized UniProt accession or raise a validation error."""

    normalized = str(accession or "").strip().upper()
    if not normalized:
        raise InvalidUniProtAccessionError("UniProt accession must not be empty")
    if _ACCESSION_RE.fullmatch(normalized) is None:
        raise InvalidUniProtAccessionError(f"invalid UniProt accession: {accession!r}")
    return normalized


def _sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _parse_record_identity(
    record_id: str,
    record_description: str = "",
) -> tuple[str, str]:
    tokens = str(record_id or "").strip().split("|")
    if len(tokens) == 3 and tokens[0].lower() in {"sp", "tr"}:
        primary_accession = normalize_uniprot_accession(tokens[1])
        entry_name = tokens[2].strip()
    else:
        primary_accession = normalize_uniprot_accession(tokens[0])
        description_tokens = str(record_description or "").strip().split()
        entry_name = (
            description_tokens[1] if len(description_tokens) > 1 else primary_accession
        )

    if not entry_name or any(character.isspace() for character in entry_name):
        raise InvalidProteinSequenceError(
            f"invalid UniProt FASTA record identifier: {record_id!r}"
        )
    return primary_accession, entry_name


def _normalize_protein_sequence(sequence: str) -> str:
    normalized = re.sub(r"\s+", "", str(sequence or "")).upper()
    if not normalized:
        raise InvalidProteinSequenceError("UniProt protein sequence is empty")

    invalid = sorted(set(normalized) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise InvalidProteinSequenceError(
            "UniProt protein sequence contains unsupported amino-acid "
            f"characters: {', '.join(invalid)}"
        )
    return normalized


def _parse_single_fasta(fasta_text: str) -> tuple[str, str, str]:
    if not isinstance(fasta_text, str) or not fasta_text.lstrip().startswith(">"):
        raise InvalidProteinSequenceError("UniProt response is not FASTA")

    try:
        records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    except Exception as exc:
        raise InvalidProteinSequenceError(
            "UniProt response could not be parsed as FASTA"
        ) from exc

    if len(records) != 1:
        raise InvalidProteinSequenceError(
            "UniProt response must contain exactly one FASTA record; "
            f"found {len(records)}"
        )

    record = records[0]
    primary_accession, entry_name = _parse_record_identity(
        record.id,
        record.description,
    )
    sequence = _normalize_protein_sequence(str(record.seq))
    return primary_accession, entry_name, sequence


def _download_uniprot_fasta(
    accession: str,
    *,
    timeout_seconds: float,
) -> tuple[str, str]:
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be greater than zero")

    source_url = UNIPROT_FASTA_URL_TEMPLATE.format(accession=accession)
    request = Request(
        source_url,
        headers={
            "Accept": "text/x-fasta",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            payload = response.read()
            get_url = getattr(response, "geturl", None)
            if callable(get_url):
                resolved_url = str(get_url() or source_url)
            else:
                resolved_url = source_url
    except HTTPError as exc:
        if exc.code == 404:
            raise UniProtNotFoundError(
                f"UniProt accession {accession} was not found"
            ) from exc
        raise UniProtServiceError(
            f"UniProt returned HTTP {exc.code} for accession {accession}"
        ) from exc
    except (URLError, TimeoutError, OSError) as exc:
        detail = getattr(exc, "reason", None) or str(exc)
        raise UniProtServiceError(
            f"UniProt request failed for accession {accession}: {detail}"
        ) from exc

    try:
        fasta_text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidProteinSequenceError(
            "UniProt FASTA response is not valid UTF-8"
        ) from exc
    return fasta_text, resolved_url


def _canonical_fasta(
    primary_accession: str,
    entry_name: str,
    sequence: str,
) -> str:
    return f">{primary_accession} {entry_name}\n{sequence}\n"


def _read_existing_sequence(path: Path) -> str:
    try:
        text = path.read_text(encoding="utf-8")
        _, _, sequence = _parse_single_fasta(text)
    except (OSError, UnicodeError, ProteinSequenceError) as exc:
        raise ProteinSequenceConflictError(
            f"existing protein FASTA is invalid and will not be overwritten: {path}"
        ) from exc
    return sequence


def _write_text_atomic(path: Path, value: str) -> None:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
            temporary_path = Path(handle.name)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def search_protein_sequence(
    accession: str,
    output_dir: str | Path,
    *,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
) -> ProteinSequenceRecord:
    """Download, validate, and persist one UniProtKB protein sequence.

    The function is intentionally independent of CLI, session, model, and
    manifest state.  Callers must provide the directory in which the canonical
    FASTA file should be stored.
    """

    requested_accession = normalize_uniprot_accession(accession)
    destination_dir = Path(output_dir).expanduser().resolve()
    try:
        destination_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise ProteinSequenceError(
            f"could not create protein-sequence directory: {destination_dir}"
        ) from exc
    if not destination_dir.is_dir():
        raise ProteinSequenceError(
            f"protein-sequence output path is not a directory: {destination_dir}"
        )

    existing_path = destination_dir / f"{requested_accession}.fasta"
    if existing_path.exists():
        try:
            existing_text = existing_path.read_text(encoding="utf-8")
            primary_accession, entry_name, sequence = _parse_single_fasta(existing_text)
        except (OSError, UnicodeError, ProteinSequenceError) as exc:
            raise ProteinSequenceConflictError(
                "existing protein FASTA is invalid and will not be overwritten: "
                f"{existing_path}"
            ) from exc
        if primary_accession != requested_accession:
            raise ProteinSequenceConflictError(
                f"existing FASTA {existing_path} identifies {primary_accession}, "
                f"not {requested_accession}"
            )
        return ProteinSequenceRecord(
            requested_accession=requested_accession,
            primary_accession=primary_accession,
            entry_name=entry_name,
            sequence=sequence,
            length_aa=len(sequence),
            sequence_sha256=_sequence_sha256(sequence),
            source_url=UNIPROT_FASTA_URL_TEMPLATE.format(accession=requested_accession),
            fasta_path=existing_path,
            reused_existing=True,
        )

    fasta_text, source_url = _download_uniprot_fasta(
        requested_accession,
        timeout_seconds=float(timeout_seconds),
    )
    primary_accession, entry_name, sequence = _parse_single_fasta(fasta_text)
    fasta_path = destination_dir / f"{primary_accession}.fasta"
    sequence_hash = _sequence_sha256(sequence)
    canonical_fasta = _canonical_fasta(primary_accession, entry_name, sequence)

    with _WRITE_LOCK:
        if fasta_path.exists():
            existing_sequence = _read_existing_sequence(fasta_path)
            if _sequence_sha256(existing_sequence) != sequence_hash:
                raise ProteinSequenceConflictError(
                    "refusing to overwrite a different protein sequence for "
                    f"UniProt accession {primary_accession}: {fasta_path}"
                )
        else:
            try:
                _write_text_atomic(fasta_path, canonical_fasta)
            except OSError as exc:
                raise ProteinSequenceError(
                    f"could not write protein FASTA: {fasta_path}"
                ) from exc

    return ProteinSequenceRecord(
        requested_accession=requested_accession,
        primary_accession=primary_accession,
        entry_name=entry_name,
        sequence=sequence,
        length_aa=len(sequence),
        sequence_sha256=sequence_hash,
        source_url=source_url,
        fasta_path=fasta_path,
        reused_existing=False,
    )


def load_uploaded_protein_sequence(
    accession: str,
    path: str | Path,
) -> ProteinSequenceRecord:
    """Load one user-uploaded amino-acid FASTA without a UniProt lookup."""

    requested_accession = str(accession or "").strip().upper()
    if not requested_accession:
        raise InvalidProteinSequenceError("uploaded protein ID must not be empty")
    fasta_path = Path(path).expanduser().resolve()
    if not fasta_path.is_file():
        raise FileNotFoundError(fasta_path)
    try:
        with fasta_path.open("r", encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except (OSError, UnicodeError) as exc:
        raise InvalidProteinSequenceError(
            f"uploaded protein FASTA could not be read: {fasta_path}"
        ) from exc
    if len(records) != 1:
        raise InvalidProteinSequenceError(
            "uploaded protein FASTA must contain exactly one record; "
            f"found {len(records)}: {fasta_path}"
        )
    record = records[0]
    record_accession = str(record.id or "").strip().upper()
    if record_accession != requested_accession:
        raise InvalidProteinSequenceError(
            "uploaded protein FASTA header does not match the manifest: "
            f"{record_accession or 'empty'} != {requested_accession}"
        )
    sequence = _normalize_protein_sequence(str(record.seq))
    description = str(record.description or "").strip()
    entry_name = (
        description.split(maxsplit=1)[1]
        if len(description.split(maxsplit=1)) == 2
        else requested_accession
    )
    return ProteinSequenceRecord(
        requested_accession=requested_accession,
        primary_accession=requested_accession,
        entry_name=entry_name,
        sequence=sequence,
        length_aa=len(sequence),
        sequence_sha256=_sequence_sha256(sequence),
        source_url=f"uploaded:{fasta_path.name}",
        fasta_path=fasta_path,
        reused_existing=True,
    )


__all__ = [
    "DEFAULT_TIMEOUT_SECONDS",
    "InvalidProteinSequenceError",
    "InvalidUniProtAccessionError",
    "load_uploaded_protein_sequence",
    "ProteinSequenceConflictError",
    "ProteinSequenceError",
    "ProteinSequenceRecord",
    "UniProtNotFoundError",
    "UniProtServiceError",
    "normalize_uniprot_accession",
    "search_protein_sequence",
]
