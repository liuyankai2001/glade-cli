"""KEGG MOL retrieval and deterministic structure normalization for RetroPath."""

from __future__ import annotations

import hashlib
import json
import math
import re
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional, Protocol
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import rdkit
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdMolDescriptors

from src.pathway_analyze.retropath_identity import structure_identity
from src.pathway_analyze.retropath_models import PredictedCompound


KEGG_REST_BASE_URL = "https://rest.kegg.jp"
KEGG_COMPOUND_ID_PATTERN = re.compile(r"^C\d{5}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
STRUCTURE_CACHE_SCHEMA_VERSION = 1
DEFAULT_HTTP_TIMEOUT_SECONDS = 30.0
DEFAULT_HTTP_RETRIES = 3
DEFAULT_REQUEST_SLEEP_SECONDS = 0.2

FetchText = Callable[[str, float], str]


class StructureProvider(Protocol):
    """Resolve one KEGG compound into the shared P1 compound model."""

    def resolve(
        self,
        compound_id: str,
        *,
        minimum_depth: Optional[int] = None,
    ) -> PredictedCompound:
        """Return one validated structure or raise ``StructureResolutionError``."""


class StructureResolutionError(ValueError):
    """A stable, auditable structure resolution failure."""

    def __init__(self, code: str, compound_id: str, detail: str) -> None:
        self.code = str(code).strip()
        self.compound_id = str(compound_id).strip()
        self.detail = str(detail).strip()
        super().__init__(f"{self.code}: {self.compound_id}: {self.detail}")


@dataclass(frozen=True)
class KeggMolCacheRecord:
    """Validated MOL cache content and its deterministic provenance."""

    compound_id: str
    mol_block: str
    mol_sha256: str
    source_url: str


def _normalize_kegg_compound_id(value: str) -> str:
    compound_id = str(value).strip().upper()
    if not KEGG_COMPOUND_ID_PATTERN.fullmatch(compound_id):
        raise StructureResolutionError(
            "invalid_kegg_compound_id",
            compound_id or str(value),
            "expected a KEGG Compound identifier in Cxxxxx format",
        )
    return compound_id


def _normalize_minimum_depth(value: Optional[int], compound_id: str) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise StructureResolutionError(
            "invalid_minimum_depth",
            compound_id,
            "minimum_depth must be null or an integer greater than or equal to 0",
        )
    return value


def _canonicalize_mol_block(value: str, compound_id: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise StructureResolutionError(
            "empty_mol_response",
            compound_id,
            "KEGG returned an empty MOL response",
        )
    # MOL headers may intentionally start with an empty title line.  Removing
    # leading newlines shifts the V2000 counts line and makes a valid file fail
    # strict parsing, so only normalize line endings and trailing newlines.
    normalized = value.replace("\r\n", "\n").replace("\r", "\n").rstrip("\n")
    normalized_message = normalized.strip().lower()
    if normalized_message in {"404", "not found"} or (
        "no such data was found" in normalized_message
    ):
        raise StructureResolutionError(
            "kegg_structure_not_found",
            compound_id,
            "KEGG did not return a MOL structure",
        )
    return normalized + "\n"


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(text)
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _structure_error(
    code: str,
    compound_id: str,
    detail: str,
    exc: Optional[BaseException] = None,
) -> StructureResolutionError:
    error = StructureResolutionError(code, compound_id, detail)
    if exc is not None:
        error.__cause__ = exc
    return error


def compound_from_kegg_mol(
    compound_id: str,
    mol_block: str,
    *,
    minimum_depth: Optional[int] = None,
    source_url: Optional[str] = None,
    mol_sha256: Optional[str] = None,
) -> PredictedCompound:
    """Convert one KEGG MOL block into a fully populated P1 compound model."""

    normalized_id = _normalize_kegg_compound_id(compound_id)
    normalized_depth = _normalize_minimum_depth(minimum_depth, normalized_id)
    normalized_mol = _canonicalize_mol_block(mol_block, normalized_id)
    observed_sha256 = _sha256_text(normalized_mol)
    if mol_sha256 is not None:
        expected_sha256 = str(mol_sha256).strip().lower()
        if not SHA256_PATTERN.fullmatch(expected_sha256):
            raise StructureResolutionError(
                "invalid_mol_sha256",
                normalized_id,
                "mol_sha256 must be a full SHA-256 digest",
            )
        if expected_sha256 != observed_sha256:
            raise StructureResolutionError(
                "mol_checksum_mismatch",
                normalized_id,
                "the supplied MOL checksum does not match the normalized MOL block",
            )

    try:
        molecule = Chem.MolFromMolBlock(
            normalized_mol,
            sanitize=False,
            removeHs=False,
            strictParsing=True,
        )
    except (RuntimeError, ValueError) as exc:
        raise _structure_error(
            "mol_parse_failed",
            normalized_id,
            str(exc) or "RDKit could not parse the MOL block",
            exc,
        )
    if molecule is None:
        raise StructureResolutionError(
            "mol_parse_failed",
            normalized_id,
            "RDKit could not parse the MOL block",
        )

    try:
        Chem.SanitizeMol(molecule)
        molecule = Chem.RemoveHs(molecule, sanitize=True)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    except (RuntimeError, ValueError) as exc:
        raise _structure_error(
            "mol_sanitize_failed",
            normalized_id,
            str(exc) or "RDKit could not sanitize the molecule",
            exc,
        )

    try:
        standard_inchi = rd_inchi.MolToInchi(molecule).strip()
    except (RuntimeError, ValueError) as exc:
        raise _structure_error(
            "inchi_generation_failed",
            normalized_id,
            str(exc) or "RDKit could not generate InChI",
            exc,
        )
    if not standard_inchi.startswith("InChI=1S/"):
        raise StructureResolutionError(
            "inchi_generation_failed",
            normalized_id,
            "RDKit did not generate a standard InChI v1 value",
        )

    try:
        inchikey = rd_inchi.InchiToInchiKey(standard_inchi).strip().upper()
    except (RuntimeError, ValueError) as exc:
        raise _structure_error(
            "inchikey_generation_failed",
            normalized_id,
            str(exc) or "RDKit could not generate InChIKey",
            exc,
        )
    if not inchikey:
        raise StructureResolutionError(
            "inchikey_generation_failed",
            normalized_id,
            "RDKit returned an empty InChIKey",
        )

    try:
        smiles = Chem.MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=True,
        ).strip()
        formula = rdMolDescriptors.CalcMolFormula(molecule).strip()
        charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    except (RuntimeError, ValueError) as exc:
        raise _structure_error(
            "structure_fields_incomplete",
            normalized_id,
            str(exc) or "RDKit could not generate all required structure fields",
            exc,
        )
    if not smiles or not formula or not isinstance(charge, int) or not math.isfinite(charge):
        raise StructureResolutionError(
            "structure_fields_incomplete",
            normalized_id,
            "SMILES, formula, or formal charge is missing",
        )

    resolved_url = (
        str(source_url).strip()
        if source_url is not None and str(source_url).strip()
        else f"{KEGG_REST_BASE_URL}/get/cpd:{normalized_id}/mol"
    )
    provenance = (
        "inchi:standard-v1",
        f"kegg_mol_sha256:{observed_sha256}",
        f"kegg_rest:{resolved_url}",
        f"rdkit:{rdkit.__version__}",
    )
    try:
        identity = structure_identity(standard_inchi)
        return PredictedCompound.create(
            compound_id=normalized_id,
            inchi=standard_inchi,
            inchikey=inchikey,
            stereo_stripped_inchikey=identity.stereo_stripped_inchikey,
            stereo_specified=identity.stereo_specified,
            isomeric_smiles=smiles,
            formula=formula,
            charge=charge,
            kegg_ids=(normalized_id,),
            minimum_depth=normalized_depth,
            structure_provenance=provenance,
        )
    except ValueError as exc:
        raise _structure_error(
            "structure_model_invalid",
            normalized_id,
            str(exc),
            exc,
        )


class KeggMolStructureProvider:
    """Resolve KEGG structures through a checksummed disk cache and REST fallback."""

    def __init__(
        self,
        cache_dir: str | Path,
        *,
        timeout_seconds: float = DEFAULT_HTTP_TIMEOUT_SECONDS,
        retries: int = DEFAULT_HTTP_RETRIES,
        request_sleep_seconds: float = DEFAULT_REQUEST_SLEEP_SECONDS,
        fetch_text: Optional[FetchText] = None,
    ) -> None:
        self.cache_dir = Path(cache_dir).expanduser().resolve() / "mol"
        self.timeout_seconds = float(timeout_seconds)
        self.retries = int(retries)
        self.request_sleep_seconds = float(request_sleep_seconds)
        self.fetch_text = fetch_text or self._fetch_text
        if not math.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be a positive finite number")
        if self.retries < 1:
            raise ValueError("retries must be greater than or equal to 1")
        if (
            not math.isfinite(self.request_sleep_seconds)
            or self.request_sleep_seconds < 0
        ):
            raise ValueError("request_sleep_seconds must be a non-negative finite number")

    @staticmethod
    def _fetch_text(url: str, timeout_seconds: float) -> str:
        request = Request(url, headers={"User-Agent": "GLADE/0.1 RetroPath-input"})
        with urlopen(request, timeout=timeout_seconds) as response:
            return response.read().decode("utf-8")

    def _paths(self, compound_id: str) -> tuple[Path, Path]:
        return (
            self.cache_dir / f"{compound_id}.mol",
            self.cache_dir / f"{compound_id}.json",
        )

    def _read_cache(self, compound_id: str, source_url: str) -> Optional[KeggMolCacheRecord]:
        mol_path, metadata_path = self._paths(compound_id)
        if not mol_path.is_file() or not metadata_path.is_file():
            return None
        try:
            mol_block = _canonicalize_mol_block(
                mol_path.read_text(encoding="utf-8"),
                compound_id,
            )
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, StructureResolutionError):
            return None
        if not isinstance(metadata, dict):
            return None
        observed_sha256 = _sha256_text(mol_block)
        if (
            metadata.get("schema_version") != STRUCTURE_CACHE_SCHEMA_VERSION
            or metadata.get("kegg_id") != compound_id
            or metadata.get("source_url") != source_url
            or metadata.get("mol_sha256") != observed_sha256
        ):
            return None
        return KeggMolCacheRecord(
            compound_id=compound_id,
            mol_block=mol_block,
            mol_sha256=observed_sha256,
            source_url=source_url,
        )

    def _fetch_with_retries(self, compound_id: str, source_url: str) -> str:
        last_error: Optional[BaseException] = None
        for attempt in range(1, self.retries + 1):
            try:
                text = self.fetch_text(source_url, self.timeout_seconds)
                mol_block = _canonicalize_mol_block(text, compound_id)
                if self.request_sleep_seconds > 0:
                    time.sleep(self.request_sleep_seconds)
                return mol_block
            except StructureResolutionError:
                raise
            except HTTPError as exc:
                if exc.code == 404:
                    raise _structure_error(
                        "kegg_structure_not_found",
                        compound_id,
                        "KEGG returned HTTP 404 for the MOL structure",
                        exc,
                    )
                last_error = exc
                if attempt < self.retries:
                    time.sleep(float(attempt))
            except (URLError, TimeoutError, OSError, UnicodeError) as exc:
                last_error = exc
                if attempt < self.retries:
                    time.sleep(float(attempt))
        detail = str(last_error) if last_error is not None else "unknown KEGG error"
        raise _structure_error(
            "kegg_fetch_failed",
            compound_id,
            detail,
            last_error,
        )

    def _write_cache(
        self,
        compound_id: str,
        source_url: str,
        mol_block: str,
    ) -> KeggMolCacheRecord:
        mol_path, metadata_path = self._paths(compound_id)
        mol_sha256 = _sha256_text(mol_block)
        metadata = {
            "schema_version": STRUCTURE_CACHE_SCHEMA_VERSION,
            "kegg_id": compound_id,
            "source_url": source_url,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "mol_sha256": mol_sha256,
        }
        _atomic_write_text(mol_path, mol_block)
        _atomic_write_text(
            metadata_path,
            json.dumps(
                metadata,
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
            )
            + "\n",
        )
        return KeggMolCacheRecord(
            compound_id=compound_id,
            mol_block=mol_block,
            mol_sha256=mol_sha256,
            source_url=source_url,
        )

    def resolve(
        self,
        compound_id: str,
        *,
        minimum_depth: Optional[int] = None,
    ) -> PredictedCompound:
        normalized_id = _normalize_kegg_compound_id(compound_id)
        normalized_depth = _normalize_minimum_depth(minimum_depth, normalized_id)
        source_url = f"{KEGG_REST_BASE_URL}/get/cpd:{normalized_id}/mol"
        record = self._read_cache(normalized_id, source_url)
        if record is not None:
            try:
                return compound_from_kegg_mol(
                    normalized_id,
                    record.mol_block,
                    minimum_depth=normalized_depth,
                    source_url=record.source_url,
                    mol_sha256=record.mol_sha256,
                )
            except StructureResolutionError:
                # A checksum-valid cache may still be chemically invalid after
                # an RDKit upgrade. Fetch it once more before reporting failure.
                pass

        mol_block = self._fetch_with_retries(normalized_id, source_url)
        mol_sha256 = _sha256_text(mol_block)
        compound = compound_from_kegg_mol(
            normalized_id,
            mol_block,
            minimum_depth=normalized_depth,
            source_url=source_url,
            mol_sha256=mol_sha256,
        )
        self._write_cache(normalized_id, source_url, mol_block)
        return compound


__all__ = [
    "DEFAULT_HTTP_RETRIES",
    "DEFAULT_HTTP_TIMEOUT_SECONDS",
    "DEFAULT_REQUEST_SLEEP_SECONDS",
    "KEGG_REST_BASE_URL",
    "KeggMolCacheRecord",
    "KeggMolStructureProvider",
    "STRUCTURE_CACHE_SCHEMA_VERSION",
    "StructureProvider",
    "StructureResolutionError",
    "compound_from_kegg_mol",
]
