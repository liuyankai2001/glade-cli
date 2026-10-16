"""Exact-structure names for RetroPath results, independent of pathway identity.

Run with ``python -m src.pathway_analyze.retropath_compound_annotation --help``.
Names are copied from database responses; candidate files are only ever read.
"""
from __future__ import annotations

import argparse
import base64
import copy
import csv
import hashlib
import io
import json
import re
import sqlite3
import tempfile
import time
import uuid
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Mapping

import httpx
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi

SCHEMA = "retropath_compound_annotations.v1"
ANNOTATIONS_FILE = "compound_annotations.json"
ANNOTATIONS_CSV = "compound_annotations.csv"
_KEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_HASH = re.compile(r"^[0-9a-f]{64}$")
_IDENTIFIER = re.compile(
    r"^(?:\d{2,7}-\d{2}-\d|(?:CID|SID|CHEBI|MNXM)\s*:?\s*\d+|SCHEMBL\d+|Q\d+)$",
    re.IGNORECASE,
)
_PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
_CHEBI = "https://www.ebi.ac.uk/chebi/backend/api/public"
_INPUT_FILES = ("candidate_routes.csv", "candidate_steps.csv")


class StructureMismatch(ValueError):
    """Returned or supplied structure differs from the complete query key."""


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _file_sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        temporary.replace(path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _write_json(path: Path, value: Any) -> None:
    _atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def _usable_name(value: Any) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    if value.strip().lower() in {"na", "none", "null"}:
        return None
    if _IDENTIFIER.fullmatch(value.strip()) or _KEY.fullmatch(value.strip()):
        return None
    return value  # Preserve the actual database string, including punctuation.


def _verify_structure(key: str, inchi: Any) -> None:
    if not isinstance(inchi, str) or not inchi.startswith("InChI=1S/"):
        raise ValueError("Missing standard InChI")
    if rd_inchi.InchiToInchiKey(inchi) != key:
        raise StructureMismatch("InChI does not reproduce the complete InChIKey")


def _input_identity(compound_id: str, inchi: str | None, smiles: str | None) -> tuple[str, str | None]:
    if not isinstance(compound_id, str) or not compound_id.startswith("RP2CPD:"):
        raise ValueError("Expected RP2CPD:<complete InChIKey or SHA-256>")
    suffix = compound_id.split(":", 1)[1]
    if inchi:
        _verify_structure(rd_inchi.InchiToInchiKey(inchi) or "", inchi)
    elif smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid source SMILES")
        inchi = rd_inchi.MolToInchi(mol)
    if _KEY.fullmatch(suffix):
        if inchi:
            _verify_structure(suffix, inchi)
        return suffix, inchi
    if _HASH.fullmatch(suffix):
        if not inchi:
            raise LookupError("SHA-256 identifiers require the original structure")
        from src.pathway_analyze.retropath_models import PredictedCompound
        if PredictedCompound.create(inchi=inchi).compound_id != compound_id:
            raise ValueError("Source structure does not reproduce the SHA-256 identifier")
        key = rd_inchi.InchiToInchiKey(inchi)
        if not key or not _KEY.fullmatch(key):
            raise ValueError("Unable to generate a complete InChIKey")
        return key, inchi
    raise ValueError("Invalid RP2CPD identifier")


class CompoundAnnotator:
    """Reusable, injectable HTTP session with exact matching and response evidence."""

    def __init__(self, cache_dir: str | Path | None = None, *,
                 evidence_dir: str | Path | None = None, offline: bool = False,
                 force: bool = False, client: httpx.Client | None = None,
                 timeout: float = 10.0, attempts: int = 3, interval: float = 0.5,
                 sleep=time.sleep, monotonic=time.monotonic):
        if offline and force:
            raise ValueError("offline and force cannot be combined")
        if attempts < 1 or timeout <= 0 or interval < 0:
            raise ValueError("Invalid request limits")
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.evidence_dir = Path(evidence_dir) if evidence_dir is not None else None
        self.offline, self.force = offline, force
        self.client, self._owns_client = client, client is None
        self.timeout, self.attempts, self.interval = timeout, attempts, interval
        self.sleep, self.monotonic = sleep, monotonic
        self._last_request = float("-inf")
        self.responses: list[dict[str, Any]] = []
        self._memory: dict[str, dict[str, Any]] = {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self) -> None:
        if self._owns_client and self.client is not None:
            self.client.close()
            self.client = None
        if self.evidence_dir is not None:
            records = [{k: v for k, v in row.items() if k != "body_base64"} for row in self.responses]
            _write_json(self.evidence_dir / "requests.json", records)

    def _capture(self, record: Mapping[str, Any]) -> None:
        row = copy.deepcopy(dict(record))
        row["record_id"] = len(self.responses) + 1
        if "body_base64" in row:
            body = base64.b64decode(row["body_base64"], validate=True)
            if hashlib.sha256(body).hexdigest() != row["body_sha256"]:
                raise ValueError("Response evidence checksum mismatch")
            if self.evidence_dir is not None:
                self.evidence_dir.mkdir(parents=True, exist_ok=True)
                row["body_file"] = f"{row['record_id']:06d}.body"
                (self.evidence_dir / row["body_file"]).write_bytes(body)
        self.responses.append(row)

    def _request(self, key: str, method: str, url: str, **kwargs) -> tuple[Any, str]:
        if self.client is None:
            self.client = httpx.Client(follow_redirects=True)
        for attempt in range(1, self.attempts + 1):
            delay = self.interval - (self.monotonic() - self._last_request)
            if delay > 0:
                self.sleep(delay)
            self._last_request = self.monotonic()
            row = {"inchikey": key, "method": method, "url": url, "arguments": kwargs,
                   "source": "network", "attempt": attempt, "retrieved_at": _now()}
            try:
                response = self.client.request(method, url, timeout=self.timeout, **kwargs)
                row.update(http_status=response.status_code, response_url=str(response.url),
                           headers=dict(response.headers),
                           body_base64=base64.b64encode(response.content).decode("ascii"),
                           body_sha256=hashlib.sha256(response.content).hexdigest())
                self._capture(row)
                if response.status_code in {429, 500, 502, 503, 504} and attempt < self.attempts:
                    self.sleep(2 ** (attempt - 1))
                    continue
                if response.status_code == 404:
                    return None, "not_found"
                if response.status_code != 200:
                    return None, "query_failed"
                try:
                    payload = response.json()
                    if not isinstance(payload, dict):
                        return None, "invalid_response"
                    return payload, "ok"
                except ValueError:
                    return None, "invalid_response"
            except httpx.HTTPError as exc:
                row.update(error_type=type(exc).__name__, error=str(exc))
                self._capture(row)
                if attempt == self.attempts:
                    return None, "query_failed"
                self.sleep(2 ** (attempt - 1))
        raise AssertionError("Unreachable retry loop")

    def _pubchem(self, key: str, inchi: str | None, smiles: str | None) -> dict[str, Any]:
        properties = "IUPACName,InChI,InChIKey,CanonicalSMILES"
        queries = []
        if inchi:
            queries.append(("POST", f"{_PUBCHEM}/inchi/property/{properties}/JSON", {"data": {"inchi": inchi}}))
        queries.append(("GET", f"{_PUBCHEM}/inchikey/{key}/property/{properties}/JSON", {}))
        if smiles:
            queries.append(("POST", f"{_PUBCHEM}/smiles/property/{properties}/JSON", {"data": {"smiles": smiles}}))
        failures = []
        for method, url, arguments in queries:
            payload, status = self._request(key, method, url, **arguments)
            if status != "ok":
                failures.append(status)
                continue
            try:
                records = payload["PropertyTable"]["Properties"]
                if not isinstance(records, list) or not records:
                    raise ValueError("Missing PubChem properties")
                # Structure property queries can return {"CID": 0} with HTTP
                # 200 when no registered compound exists. Keep trying the
                # other identifiers; this placeholder cannot supply a name.
                records = [row for row in records if not (
                    isinstance(row, dict) and set(row) == {"CID"}
                    and type(row["CID"]) is int and row["CID"] == 0)]
                if not records:
                    failures.append("not_found")
                    continue
                matches = []
                for props in records:
                    if props["InChIKey"] != key:
                        raise StructureMismatch("Returned complete InChIKey differs")
                    _verify_structure(key, props["InChI"])
                    cid = props["CID"]
                    if isinstance(cid, bool) or not isinstance(cid, int) or cid <= 0:
                        raise ValueError("Invalid PubChem CID")
                    name = _usable_name(props.get("IUPACName"))
                    if name:
                        matches.append({"name": name, "name_kind": "iupac", "database": "PubChem",
                            "database_id": str(cid), "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                            "inchikey": key, "inchi": props["InChI"],
                            "smiles": props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")})
                names = {record["name"] for record in matches}
                return {"status": "exact" if len(names) == 1 else "ambiguous" if names else "no_name",
                        "matches": sorted(matches, key=lambda row: int(row["database_id"]))}
            except StructureMismatch:
                failures.append("structure_mismatch")
            except (ValueError, KeyError, TypeError, AttributeError, RuntimeError):
                failures.append("invalid_response")
        status = next((s for s in ("query_failed", "structure_mismatch", "invalid_response") if s in failures), "not_found")
        return {"status": status, "matches": []}

    def _chebi(self, key: str) -> dict[str, Any]:
        payload, status = self._request(key, "GET", f"{_CHEBI}/es_search/",
                                       params={"term": key, "size": 50, "page": 1})
        if status != "ok":
            return {"status": status, "matches": []}
        try:
            hits = payload["results"]
            if not isinstance(hits, list):
                raise ValueError("Missing ChEBI search results")
            if payload.get("number_pages", 1) > 1:
                return {"status": "ambiguous", "matches": []}
            identifiers = sorted({str(hit["_id"]) for hit in hits
                                  if hit["_source"].get("inchikey") == key})
            if not identifiers:
                return {"status": "structure_mismatch" if hits else "not_found", "matches": []}
            matches, errors = [], []
            for identifier in identifiers:
                if not identifier.isdigit():
                    errors.append("invalid_response")
                    continue
                detail, result = self._request(key, "GET", f"{_CHEBI}/compound/{identifier}/")
                if result != "ok":
                    errors.append(result)
                    continue
                structure = detail["default_structure"]
                if structure["standard_inchi_key"] != key:
                    errors.append("structure_mismatch")
                    continue
                _verify_structure(key, structure["standard_inchi"])
                if str(detail["id"]) != identifier or detail.get("is_released") is False:
                    raise ValueError("Invalid ChEBI accession")
                name = _usable_name(detail.get("name"))
                if name:
                    matches.append({"name": name, "name_kind": "database_preferred", "database": "ChEBI",
                        "database_id": f"CHEBI:{identifier}", "url": f"https://www.ebi.ac.uk/chebi/CHEBI:{identifier}",
                        "inchikey": key, "inchi": structure["standard_inchi"], "smiles": structure.get("smiles")})
            names = {record["name"] for record in matches}
            status = errors[0] if errors else "exact" if len(names) == 1 else "ambiguous" if names else "no_name"
            return {"status": status, "matches": matches}
        except StructureMismatch:
            return {"status": "structure_mismatch", "matches": []}
        except (ValueError, KeyError, TypeError, AttributeError, RuntimeError):
            return {"status": "invalid_response", "matches": []}

    def _cached(self, key: str) -> dict[str, Any] | None:
        if self.cache_dir is None or self.force:
            return None
        try:
            wrapper = json.loads((self.cache_dir / f"{key}.json").read_text(encoding="utf-8"))
            data = wrapper["data"]
            if wrapper["sha256"] != _digest(data) or data["schema_version"] != SCHEMA:
                return None
            result = data["annotation"]
            if result["inchikey"] != key or result["status"] != "resolved":
                return None
            _validate_resolved(result)
            for row in data["responses"]:
                if ("body_base64" in row and
                        hashlib.sha256(base64.b64decode(row["body_base64"], validate=True)).hexdigest() != row["body_sha256"]):
                    return None
            _verify_name_provenance(result, data["responses"],
                                    lambda row: base64.b64decode(row["body_base64"], validate=True))
            for row in data["responses"]:
                self._capture({**row, "source": "cache", "reused_at": _now()})
            return copy.deepcopy(result)
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            return None

    def resolve(self, compound_id: str, *, inchi: str | None = None, smiles: str | None = None) -> dict[str, Any]:
        result = {"compound_id": compound_id, "inchikey": None, "status": "invalid_input",
                  "display_name": compound_id, "name": None, "name_kind": None,
                  "iupac_name": None, "source": None, "lookups": {}, "retrieved_at": _now(),
                  "cache_hit": False, "evidence_request_ids": [], "input_inchi": inchi, "input_smiles": smiles}
        try:
            key, inchi = _input_identity(compound_id, inchi, smiles)
            result.update(inchikey=key, input_inchi=inchi)
        except LookupError as exc:
            return {**result, "status": "missing_structure", "detail": str(exc)}
        except (ValueError, TypeError, RuntimeError) as exc:
            return {**result, "detail": str(exc)}
        start = len(self.responses)
        cached = copy.deepcopy(self._memory.get(key)) or self._cached(key)
        if cached is not None:
            cached.update(compound_id=compound_id, input_inchi=inchi, input_smiles=smiles, cache_hit=True)
            if not cached.get("name"):
                cached["display_name"] = compound_id
            if len(self.responses) > start:
                cached["evidence_request_ids"] = list(range(start + 1, len(self.responses) + 1))
            return cached
        if self.offline:
            return {**result, "status": "offline", "detail": "No verified cached annotation"}
        pubchem = self._pubchem(key, inchi, smiles)
        chebi = self._chebi(key)
        result["lookups"] = {"PubChem": pubchem, "ChEBI": chebi}
        result["evidence_request_ids"] = list(range(start + 1, len(self.responses) + 1))
        if pubchem["status"] == "exact":
            result["iupac_name"] = pubchem["matches"][0]["name"]
        chosen = chebi if chebi["status"] == "exact" else pubchem if pubchem["status"] == "exact" else None
        if chosen:
            source = chosen["matches"][0]
            result.update(status="resolved", name=source["name"], display_name=source["name"],
                          name_kind=source["name_kind"], source=source)
        else:
            statuses = {pubchem["status"], chebi["status"]}
            result["status"] = next((s for s in ("query_failed", "invalid_response", "structure_mismatch", "ambiguous")
                                     if s in statuses), "not_found")
        self._memory[key] = copy.deepcopy(result)
        complete = all(r["status"] in {"exact", "not_found", "no_name"} for r in (pubchem, chebi))
        if self.cache_dir is not None and chosen and complete:
            data = {"schema_version": SCHEMA, "annotation": result, "responses": self.responses[start:]}
            _write_json(self.cache_dir / f"{key}.json", {"data": data, "sha256": _digest(data)})
        return result


def _validate_resolved(row: Mapping[str, Any]) -> None:
    key, _ = _input_identity(row["compound_id"], row.get("input_inchi"), row.get("input_smiles"))
    source = row["source"]
    if key != row["inchikey"] or source["inchikey"] != key:
        raise ValueError("Annotation structure differs from input")
    _verify_structure(key, source["inchi"])
    if not _usable_name(row["name"]) or row["name"] != source["name"] or row["display_name"] != row["name"]:
        raise ValueError("Invalid annotation name")


def _verify_name_provenance(annotation, records, read_body):
    source, key = annotation["source"], annotation["inchikey"]
    name_found, iupac_found = False, annotation.get("iupac_name") is None
    for record in records:
        if record.get("inchikey") != key or record.get("http_status") != 200:
            continue
        try:
            payload = json.loads(read_body(record))
        except (ValueError, UnicodeDecodeError):
            continue
        if not isinstance(payload, dict):
            continue
        for props in payload.get("PropertyTable", {}).get("Properties", []):
            if props.get("InChIKey") != key:
                continue
            if props.get("IUPACName") == annotation.get("iupac_name"):
                iupac_found = True
            if (source["database"] == "PubChem" and str(props.get("CID")) == source["database_id"]
                    and props.get("IUPACName") == annotation["name"]):
                name_found = True
        structure = payload.get("default_structure") or {}
        if (source["database"] == "ChEBI" and f"CHEBI:{payload.get('id')}" == source["database_id"]
                and structure.get("standard_inchi_key") == key and payload.get("name") == annotation["name"]):
            name_found = True
    if not name_found or not iupac_found:
        raise ValueError("Annotation name is absent from its original database responses")


def candidate_compound_ids(rows: Iterable[Mapping[str, Any]]) -> list[str]:
    ids = set()
    for row in rows:
        for field in ("substrate_stoichiometry_json", "product_stoichiometry_json"):
            values = json.loads(row[field])
            if not isinstance(values, list):
                raise ValueError(f"Invalid {field}")
            for identifier, _ in values:
                if str(identifier).startswith("RP2CPD:"):
                    ids.add(str(identifier))
    return sorted(ids)


def _input_hashes(directory: Path) -> dict[str, str]:
    return {name: _file_sha(directory / name) for name in _INPUT_FILES}


def annotate_retropath_results(retropath_dir: str | Path, *, compounds: Iterable[Any] = (),
                               output_dir: str | Path | None = None, cache_dir: str | Path | None = None,
                               offline: bool = False, force: bool = False,
                               client: httpx.Client | None = None, **request_options) -> dict[str, Any]:
    """Annotate final candidates; write an independent, input-bound sidecar only."""
    directory = Path(retropath_dir).resolve()
    destination = Path(output_dir).resolve() if output_dir is not None else directory
    inputs = _input_hashes(directory)
    with (directory / "candidate_steps.csv").open(encoding="utf-8-sig", newline="") as stream:
        ids = candidate_compound_ids(csv.DictReader(stream))
    structures = {}
    for compound in compounds:
        data = dict(compound) if isinstance(compound, Mapping) else compound.to_dict()
        structures[data["compound_id"]] = data
    if cache_dir is None:
        from src.config.config import CACHE_DIR
        cache_dir = CACHE_DIR / "retropath" / "compound_annotation"
    token = datetime.now(UTC).strftime("%Y%m%dT%H%M%S_%fZ") + "_" + uuid.uuid4().hex[:8]
    evidence_dir = destination / "compound_annotation_evidence" / token
    with CompoundAnnotator(cache_dir, evidence_dir=evidence_dir, offline=offline, force=force,
                           client=client, **request_options) as annotator:
        annotations = [annotator.resolve(identifier, inchi=structures.get(identifier, {}).get("inchi"),
                                         smiles=structures.get(identifier, {}).get("isomeric_smiles")) for identifier in ids]
    if _input_hashes(directory) != inputs:
        raise ValueError("Candidate inputs changed during annotation")
    counts = dict(Counter(row["status"] for row in annotations))
    payload = {"schema_version": SCHEMA, "generated_at": _now(), "inputs": inputs,
               "annotations": annotations, "annotations_sha256": _digest(annotations),
               "evidence_directory": evidence_dir.relative_to(destination).as_posix(),
               "requests_sha256": _file_sha(evidence_dir / "requests.json"),
               "summary": {"total": len(ids), "resolved": counts.get("resolved", 0), "statuses": counts}}
    buffer = io.StringIO(newline="")
    fields = ["compound_id", "inchikey", "status", "display_name", "iupac_name", "name_kind",
              "source_database", "source_id", "source_url", "retrieved_at"]
    writer = csv.DictWriter(buffer, fields)
    writer.writeheader()
    for row in annotations:
        source = row.get("source") or {}
        writer.writerow({**{key: row.get(key) for key in fields if not key.startswith("source_")},
                         "source_database": source.get("database"), "source_id": source.get("database_id"),
                         "source_url": source.get("url")})
    _atomic_text(destination / ANNOTATIONS_CSV, buffer.getvalue())
    payload["csv_sha256"] = _file_sha(destination / ANNOTATIONS_CSV)
    _write_json(destination / ANNOTATIONS_FILE, payload)
    return {**payload["summary"], "state": "completed", "json_path": str(destination / ANNOTATIONS_FILE),
            "csv_path": str(destination / ANNOTATIONS_CSV), "sha256": _file_sha(destination / ANNOTATIONS_FILE)}


def load_compound_annotations(retropath_dir: str | Path, *, annotation_dir: str | Path | None = None) -> dict[str, Any]:
    """Read and validate saved evidence without creating files or using HTTP."""
    input_directory = Path(retropath_dir).resolve()
    directory = Path(annotation_dir).resolve() if annotation_dir is not None else input_directory
    path = directory / ANNOTATIONS_FILE
    if not path.exists():
        return {"state": "missing", "annotations": {}, "summary": {}}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload["schema_version"] != SCHEMA:
            raise ValueError("Unsupported annotation schema")
        if payload["inputs"] != _input_hashes(input_directory):
            return {"state": "stale", "annotations": {}, "summary": {}}
        if payload["annotations_sha256"] != _digest(payload["annotations"]):
            raise ValueError("Annotation table checksum mismatch")
        if "csv_sha256" in payload and _file_sha(directory / ANNOTATIONS_CSV) != payload["csv_sha256"]:
            raise ValueError("Annotation CSV checksum mismatch")
        evidence = (directory / payload["evidence_directory"]).resolve()
        if not evidence.is_relative_to(directory / "compound_annotation_evidence"):
            raise ValueError("Evidence directory escapes annotation output")
        if _file_sha(evidence / "requests.json") != payload["requests_sha256"]:
            raise ValueError("Request index checksum mismatch")
        records = json.loads((evidence / "requests.json").read_text(encoding="utf-8"))
        indexed = {row["record_id"]: row for row in records}
        for record in records:
            if "body_file" in record:
                body = (evidence / record["body_file"]).resolve()
                if body.parent != evidence or _file_sha(body) != record["body_sha256"]:
                    raise ValueError("Original response checksum mismatch")
        annotations = {}
        for row in payload["annotations"]:
            if row["compound_id"] in annotations:
                raise ValueError("Duplicate annotation identifier")
            if row["status"] == "resolved":
                _validate_resolved(row)
                if not row["evidence_request_ids"] or any(indexed[i]["inchikey"] != row["inchikey"]
                                                        for i in row["evidence_request_ids"]):
                    raise ValueError("Missing response provenance")
                _verify_name_provenance(row, [indexed[i] for i in row["evidence_request_ids"]],
                                        lambda record: (evidence / record["body_file"]).read_bytes())
            annotations[row["compound_id"]] = row
        with (input_directory / "candidate_steps.csv").open(encoding="utf-8-sig", newline="") as stream:
            expected = set(candidate_compound_ids(csv.DictReader(stream)))
        if set(annotations) != expected:
            raise ValueError("Annotation identifiers differ from candidate compounds")
        counts = dict(Counter(row["status"] for row in annotations.values()))
        if payload["summary"] != {"total": len(expected), "resolved": counts.get("resolved", 0), "statuses": counts}:
            raise ValueError("Annotation summary differs from the saved records")
        return {"state": "current", "annotations": annotations, "summary": payload["summary"]}
    except (OSError, ValueError, KeyError, TypeError, AttributeError, RuntimeError) as exc:
        return {"state": "invalid", "annotations": {}, "summary": {}, "detail": str(exc)}


def legacy_exact_annotations(compound_ids: Iterable[str], *, mnxref_dir: str | Path | None = None,
                             rules_path: str | Path | None = None) -> dict[str, Any]:
    """Read-only compatibility for old results; never promote a connectivity match."""
    from src.pathway_analyze.retropath_mnxref import MnxrefIndex, MnxrefIndexError
    keys = {identifier: identifier.split(":", 1)[1] for identifier in compound_ids
            if identifier.startswith("RP2CPD:") and _KEY.fullmatch(identifier.split(":", 1)[1])}
    if not keys:
        return {}
    try:
        with MnxrefIndex(mnxref_dir, rules_path) as index:
            matches = index.chemicals_by_connectivity_keys(key[:14] for key in keys.values())
            result = {}
            for identifier, key in keys.items():
                exact = [item for item in matches.get(key[:14], ())
                         if item.inchikey == key and _usable_name(item.name)]
                if len({item.name for item in exact}) != 1:
                    continue
                item = exact[0]
                _verify_structure(key, item.inchi)
                result[identifier] = {"compound_id": identifier, "inchikey": key, "status": "resolved",
                    "name": item.name, "display_name": item.name, "iupac_name": None,
                    "source": {"database": "MetaNetX", "database_id": item.mnxm_id,
                               "url": f"https://www.metanetx.org/chem_info/{item.mnxm_id}"}}
            return result
    except (OSError, ValueError, TypeError, RuntimeError, sqlite3.Error, MnxrefIndexError):
        return {}


def annotation_view(retropath_dir: str | Path, *, compound_ids: Iterable[str] = (), config=None) -> dict[str, Any]:
    view = load_compound_annotations(retropath_dir)
    if view["state"] == "missing":
        data_dir = getattr(config, "data_dir", None)
        local = legacy_exact_annotations(compound_ids,
            mnxref_dir=Path(data_dir) / "retropath/mnxref/3.0" if data_dir is not None else None,
            rules_path=getattr(config, "retropath_rules_path", None))
        if local:
            return {"state": "legacy_exact", "annotations": local, "summary": {"resolved": len(local)}}
    return view


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retropath-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    switches = parser.add_mutually_exclusive_group()
    switches.add_argument("--offline", action="store_true")
    switches.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = annotate_retropath_results(**vars(args))
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (OSError, ValueError) as exc:
        print(json.dumps({"state": "failed", "detail": str(exc)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
