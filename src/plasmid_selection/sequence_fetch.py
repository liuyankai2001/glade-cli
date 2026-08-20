"""Download and audit authoritative plasmid GenBank sequences."""

from __future__ import annotations

import hashlib
import io
import re
import warnings
from collections.abc import Callable, Mapping
from typing import Any
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen

from Bio import SeqIO
from Bio import BiopythonParserWarning

from src.plasmid_selection.config import DOWNLOAD_TIMEOUT_SECONDS
from src.plasmid_selection.models import DownloadedSequence


_ACCESSION_RE = re.compile(r"^[A-Z]{1,6}[A-Z0-9_]*\d+(?:\.\d+)?$", re.IGNORECASE)


class PlasmidSequenceError(RuntimeError):
    """A candidate source sequence is unavailable or fails its audit."""


def _minimum_rotation(sequence: str) -> str:
    """Return the lexicographically minimal circular rotation (Booth)."""

    if not sequence:
        return sequence
    doubled = sequence + sequence
    length = len(sequence)
    left, right, offset = 0, 1, 0
    while left < length and right < length and offset < length:
        a = doubled[left + offset]
        b = doubled[right + offset]
        if a == b:
            offset += 1
            continue
        if a > b:
            left = left + offset + 1
            if left <= right:
                left = right + 1
        else:
            right = right + offset + 1
            if right <= left:
                right = left + 1
        offset = 0
    start = min(left, right)
    return doubled[start : start + length]


def canonical_circular_sequence(sequence: str) -> str:
    sequence = sequence.upper()
    reverse_complement = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(_minimum_rotation(sequence), _minimum_rotation(reverse_complement))


def canonical_circular_sha256(sequence: str) -> str:
    canonical = canonical_circular_sequence(sequence)
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def _source_url(template: Mapping[str, Any]) -> str:
    source = str(template.get("source") or "").strip().lower()
    accession = str(template.get("source_record_id") or "").strip()
    source_url = str(template.get("source_url") or "").strip()
    sequence_url = str(template.get("sequence_url") or "").strip()
    if (
        accession
        and _ACCESSION_RE.fullmatch(accession)
        and (source == "ncbi" or "ncbi.nlm.nih.gov" in source_url.lower())
    ):
        return (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={quote(accession)}&rettype=gbwithparts&retmode=text"
        )
    for value in (sequence_url, source_url):
        parsed = urlparse(value)
        lowered_path = parsed.path.lower()
        if parsed.scheme in {"http", "https"} and (
            lowered_path.endswith((".gb", ".gbk", ".gbff", ".genbank"))
            or "download" in lowered_path
        ):
            return value
    if source == "addgene" or "addgene.org" in sequence_url.lower():
        raise PlasmidSequenceError(
            "Addgene sequence page is not a raw GenBank download"
        )
    raise PlasmidSequenceError("no supported raw GenBank source is available")


def _read_response_bytes(response: Any) -> bytes:
    content = response.read()
    if not isinstance(content, bytes):
        content = bytes(content)
    return content


def _open_url(
    url: str,
    *,
    opener: Callable[..., Any],
    timeout: int,
) -> bytes:
    request = Request(
        url,
        headers={
            "User-Agent": "GLADE/0.1 plasmid-selection (GenBank validation)",
            "Accept": "text/plain, application/octet-stream;q=0.9, */*;q=0.1",
        },
    )
    try:
        response = opener(request, timeout=timeout)
        if hasattr(response, "__enter__"):
            with response as handle:
                return _read_response_bytes(handle)
        try:
            return _read_response_bytes(response)
        finally:
            close = getattr(response, "close", None)
            if callable(close):
                close()
    except PlasmidSequenceError:
        raise
    except Exception as exc:
        raise PlasmidSequenceError(
            f"sequence download failed: {type(exc).__name__}"
        ) from exc


def validate_genbank_bytes(
    content: bytes,
    template: Mapping[str, Any],
    *,
    source_download_url: str,
) -> DownloadedSequence:
    if not content.strip():
        raise PlasmidSequenceError("downloaded GenBank file is empty")
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = content.decode("latin-1")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    except Exception as exc:
        raise PlasmidSequenceError("downloaded file is not valid GenBank") from exc
    if len(records) != 1:
        raise PlasmidSequenceError(
            f"expected one GenBank record, found {len(records)}"
        )
    record = records[0]
    sequence = str(record.seq).upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise PlasmidSequenceError(
            "plasmid sequence must contain only A/C/G/T"
        )
    topology = str(record.annotations.get("topology") or "").lower()
    if topology != "circular":
        raise PlasmidSequenceError("downloaded GenBank record is not circular")
    try:
        expected_length = int(template.get("length_bp") or 0)
    except (TypeError, ValueError) as exc:
        raise PlasmidSequenceError("template length_bp is invalid") from exc
    if expected_length < 1 or len(sequence) != expected_length:
        raise PlasmidSequenceError(
            f"sequence length mismatch: expected {expected_length}, got {len(sequence)}"
        )

    content_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    canonical_hash = canonical_circular_sha256(sequence)
    expected_content_hash = str(
        template.get("sequence_content_sha256") or ""
    ).strip().lower()
    expected_canonical_hash = str(
        template.get("canonical_sequence_sha256") or ""
    ).strip().lower()
    content_matches = content_hash == expected_content_hash
    canonical_matches = canonical_hash == expected_canonical_hash
    if not content_matches and not canonical_matches:
        raise PlasmidSequenceError(
            "downloaded sequence does not match the audited content/canonical SHA-256"
        )
    return DownloadedSequence(
        content=content,
        length_bp=len(sequence),
        sequence_content_sha256=content_hash,
        canonical_sequence_sha256=canonical_hash,
        file_sha256=hashlib.sha256(content).hexdigest(),
        source_download_url=source_download_url,
    )


def download_and_validate_template(
    template: Mapping[str, Any],
    *,
    opener: Callable[..., Any] = urlopen,
    timeout: int = DOWNLOAD_TIMEOUT_SECONDS,
) -> DownloadedSequence:
    url = _source_url(template)
    content = _open_url(url, opener=opener, timeout=timeout)
    return validate_genbank_bytes(
        content,
        template,
        source_download_url=url,
    )


__all__ = [
    "PlasmidSequenceError",
    "canonical_circular_sequence",
    "canonical_circular_sha256",
    "download_and_validate_template",
    "validate_genbank_bytes",
]
