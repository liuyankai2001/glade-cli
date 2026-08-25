"""Pinned MNXref 3.0 subset installer used by RetroPath stoichiometry repair."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import re
import shutil
import sqlite3
import tempfile
import time
import urllib.request
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


MNXREF_VERSION = "3.0"
MNXREF_RELEASE = "2017/05/04"
MNXREF_BASE_URL = f"https://www.metanetx.org/ftp/{MNXREF_VERSION}"
MNXREF_FILES = (
    "reac_prop.tsv",
    "reac_xref.tsv",
    "chem_prop.tsv",
    "chem_xref.tsv",
)
MNXREF_INDEX_SCHEMA = "mnxref_rr02_subset.v1"
MNXREF_MANIFEST_SCHEMA = "mnxref_rr02_subset_manifest.v1"
MNXREF_INDEX_FILE_NAME = "mnxref_rr02_subset.sqlite3"
MNXREF_MANIFEST_FILE_NAME = "mnxref_rr02_subset_manifest.json"
DEFAULT_RULES_RELATIVE_PATH = Path(
    "data/retropath/rules/rr02/retrorules_rr02_rp2_flat_retro.csv"
)
DEFAULT_INSTALL_RELATIVE_PATH = Path("data/retropath/mnxref/3.0")

_RULE_REQUIRED_COLUMNS = {
    "Rule ID",
    "Legacy ID",
    "Reaction direction",
    "Rule relative direction",
    "Rule usage",
}
_MNXR_PATTERN = re.compile(r"MNXR\d+")
_MNXM_PATTERN = re.compile(r"MNXM\d+")
_TERM_PATTERN = re.compile(
    r"^\s*(?P<coefficient>\d+(?:\.\d+)?)\s+"
    r"(?P<compound>MNXM\d+)(?:@(?P<compartment>MNX[CD]\d+))?\s*$"
)

DownloadFile = Callable[[str, Path], None]


@dataclass(frozen=True)
class MnxrefChemical:
    mnxm_id: str
    name: str
    formula: str
    charge: int | None
    mass: float | None
    inchi: str
    smiles: str
    reference: str
    inchikey: str
    xrefs: tuple[str, ...] = tuple()


@dataclass(frozen=True)
class MnxrefReactionTerm:
    side: str
    coefficient: float
    mnxm_id: str
    compartment: str
    ordinal: int


@dataclass(frozen=True)
class MnxrefReactionTemplate:
    rule_id: str
    mnxr_id: str
    main_mnxm_id: str
    reaction_direction: str
    rule_relative_direction: str
    rule_usage: str
    equation: str
    balanced: bool
    transport: bool
    reference: str
    parse_status: str
    terms: tuple[MnxrefReactionTerm, ...]
    reaction_xrefs: tuple[str, ...]


class MnxrefIndexError(RuntimeError):
    """Stable failure for an absent, stale, or inconsistent MNXref index."""


def default_project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_rules_path() -> Path:
    return default_project_root() / DEFAULT_RULES_RELATIVE_PATH


def default_install_dir() -> Path:
    return default_project_root() / DEFAULT_INSTALL_RELATIVE_PATH


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
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
            json.dump(
                payload,
                handle,
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
                allow_nan=False,
            )
            handle.write("\n")
            temporary_path = Path(handle.name)
        temporary_path.replace(path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _download_file(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    user_agent = "GLADE/0.1 P8"
    probe = urllib.request.Request(
        url,
        headers={"User-Agent": user_agent, "Range": "bytes=0-0"},
    )
    with urllib.request.urlopen(probe, timeout=120) as response:
        content_range = str(response.headers.get("Content-Range") or "")
        match = re.fullmatch(r"bytes\s+0-0/(\d+)", content_range)
        if match is None:
            with path.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
            return
        total_size = int(match.group(1))
    chunk_size = 8 * 1024 * 1024
    ranges = [
        (start, min(total_size - 1, start + chunk_size - 1))
        for start in range(0, total_size, chunk_size)
    ]

    def fetch_range(start: int, end: int) -> tuple[int, bytes]:
        last_error: Exception | None = None
        for attempt in range(3):
            try:
                request = urllib.request.Request(
                    url,
                    headers={
                        "User-Agent": user_agent,
                        "Range": f"bytes={start}-{end}",
                    },
                )
                with urllib.request.urlopen(request, timeout=120) as response:
                    if getattr(response, "status", 206) != 206:
                        raise MnxrefIndexError(
                            f"server ignored byte range for {url}"
                        )
                    block = response.read()
                expected = end - start + 1
                if len(block) != expected:
                    raise MnxrefIndexError(
                        f"short byte range for {url}: {len(block)} != {expected}"
                    )
                return start, block
            except Exception as exc:
                last_error = exc
                if attempt == 2:
                    raise MnxrefIndexError(
                        f"failed to download byte range {start}-{end}: {url}"
                    ) from last_error
                time.sleep(0.5 * (2**attempt))
        raise AssertionError("unreachable download retry state")

    completed_bytes = 0
    with path.open("w+b") as handle:
        handle.truncate(total_size)
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(fetch_range, start, end): (start, end)
                for start, end in ranges
            }
            for future in concurrent.futures.as_completed(futures):
                start, end = futures[future]
                block_start, block = future.result()
                handle.seek(block_start)
                handle.write(block)
                completed_bytes += len(block)
            percent = 100.0 * completed_bytes / total_size
            print(
                f"[INFO] downloaded {path.name}: {percent:.1f}% "
                f"({completed_bytes}/{total_size} bytes)",
                flush=True,
            )


def _read_official_md5(path: Path) -> str:
    token = path.read_text(encoding="ascii").strip().split()[0].lower()
    if not re.fullmatch(r"[0-9a-f]{32}", token):
        raise MnxrefIndexError(f"invalid official MD5 file: {path}")
    return token


def _rule_links(
    rules_path: Path,
) -> tuple[list[tuple[str, str, str, str, str, str]], set[str], set[str]]:
    rows: set[tuple[str, str, str, str, str, str]] = set()
    reaction_ids: set[str] = set()
    main_compound_ids: set[str] = set()
    with rules_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(_RULE_REQUIRED_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise MnxrefIndexError(f"RR02 file is missing columns: {missing}")
        for raw in reader:
            rule_id = str(raw.get("Rule ID") or "").strip()
            legacy = str(raw.get("Legacy ID") or "").strip()
            mnxr_ids = _MNXR_PATTERN.findall(legacy)
            mnxm_ids = _MNXM_PATTERN.findall(legacy)
            if not rule_id or not mnxr_ids or not mnxm_ids:
                continue
            for mnxr_id in mnxr_ids:
                row = (
                    rule_id,
                    mnxr_id,
                    mnxm_ids[0],
                    str(raw.get("Reaction direction") or "").strip(),
                    str(raw.get("Rule relative direction") or "").strip(),
                    str(raw.get("Rule usage") or "").strip(),
                )
                rows.add(row)
                reaction_ids.add(mnxr_id)
                main_compound_ids.add(mnxm_ids[0])
    if not rows:
        raise MnxrefIndexError("RR02 file contains no MNXR/MNXM legacy links")
    return sorted(rows), reaction_ids, main_compound_ids


def _data_rows(path: Path) -> Iterable[list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or not row[0] or row[0].startswith("#"):
                continue
            yield row


def _as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _parse_equation(
    equation: str,
) -> tuple[str, list[tuple[str, float, str, str, int]]]:
    if equation.count("=") != 1:
        return "equation_invalid", []
    left, right = equation.split("=", 1)
    terms: list[tuple[str, float, str, str, int]] = []
    for side_name, raw_side in (("left", left), ("right", right)):
        raw_terms = [item.strip() for item in raw_side.split("+") if item.strip()]
        if not raw_terms:
            return "equation_invalid", []
        for ordinal, raw_term in enumerate(raw_terms, start=1):
            match = _TERM_PATTERN.fullmatch(raw_term)
            if match is None:
                return "nonnumeric_or_invalid_term", []
            terms.append(
                (
                    side_name,
                    float(match.group("coefficient")),
                    match.group("compound"),
                    match.group("compartment") or "MNXD1",
                    ordinal,
                )
            )
    return "ok", terms


def _create_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE rule_templates (
            rule_id TEXT NOT NULL,
            mnxr_id TEXT NOT NULL,
            main_mnxm_id TEXT NOT NULL,
            reaction_direction TEXT NOT NULL,
            rule_relative_direction TEXT NOT NULL,
            rule_usage TEXT NOT NULL,
            PRIMARY KEY (rule_id, mnxr_id, main_mnxm_id)
        );
        CREATE TABLE reactions (
            mnxr_id TEXT PRIMARY KEY,
            equation TEXT NOT NULL,
            human_equation TEXT NOT NULL,
            balanced INTEGER NOT NULL,
            transport INTEGER NOT NULL,
            reference TEXT NOT NULL,
            parse_status TEXT NOT NULL
        );
        CREATE TABLE reaction_terms (
            mnxr_id TEXT NOT NULL,
            side TEXT NOT NULL,
            coefficient REAL NOT NULL,
            mnxm_id TEXT NOT NULL,
            compartment TEXT NOT NULL,
            ordinal INTEGER NOT NULL,
            PRIMARY KEY (mnxr_id, side, ordinal)
        );
        CREATE TABLE chemicals (
            mnxm_id TEXT PRIMARY KEY,
            name TEXT NOT NULL,
            formula TEXT NOT NULL,
            charge INTEGER,
            mass REAL,
            inchi TEXT NOT NULL,
            smiles TEXT NOT NULL,
            reference TEXT NOT NULL,
            inchikey TEXT NOT NULL
        );
        CREATE TABLE chemical_xrefs (
            xref TEXT NOT NULL,
            mnxm_id TEXT NOT NULL,
            evidence TEXT NOT NULL,
            description TEXT NOT NULL,
            PRIMARY KEY (xref, mnxm_id)
        );
        CREATE TABLE reaction_xrefs (
            xref TEXT NOT NULL,
            mnxr_id TEXT NOT NULL,
            PRIMARY KEY (xref, mnxr_id)
        );
        CREATE INDEX idx_rule_templates_rule ON rule_templates(rule_id);
        CREATE INDEX idx_terms_compound ON reaction_terms(mnxm_id);
        CREATE INDEX idx_chemical_xrefs_mnxm ON chemical_xrefs(mnxm_id);
        CREATE INDEX idx_reaction_xrefs_mnxr ON reaction_xrefs(mnxr_id);
        """
    )


def build_mnxref_subset(
    *,
    rules_path: str | Path,
    source_paths: Mapping[str, str | Path],
    output_dir: str | Path,
    source_metadata: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Build an atomic SQLite subset from already downloaded MNXref 3.0 TSVs."""

    resolved_rules = Path(rules_path).expanduser().resolve()
    resolved_output = Path(output_dir).expanduser().resolve()
    if not resolved_rules.is_file():
        raise FileNotFoundError(f"RR02 rules file not found: {resolved_rules}")
    normalized_sources = {
        name: Path(source_paths[name]).expanduser().resolve() for name in MNXREF_FILES
    }
    missing_sources = [str(path) for path in normalized_sources.values() if not path.is_file()]
    if missing_sources:
        raise FileNotFoundError(f"MNXref source files not found: {missing_sources}")
    rule_rows, requested_reactions, requested_main_compounds = _rule_links(
        resolved_rules
    )
    resolved_output.mkdir(parents=True, exist_ok=True)
    temporary_db: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=resolved_output,
            prefix=f".{MNXREF_INDEX_FILE_NAME}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_db = Path(handle.name)
        connection = sqlite3.connect(temporary_db)
        try:
            _create_schema(connection)
            connection.executemany(
                "INSERT INTO rule_templates VALUES (?, ?, ?, ?, ?, ?)",
                rule_rows,
            )
            required_compounds = set(requested_main_compounds)
            reaction_count = 0
            term_count = 0
            for row in _data_rows(normalized_sources["reac_prop.tsv"]):
                if len(row) < 6 or row[0] not in requested_reactions:
                    continue
                mnxr_id, equation, human_equation = row[:3]
                balanced = _as_bool(row[3])
                transport = _as_bool(row[4])
                reference = row[5].strip()
                parse_status, terms = _parse_equation(equation)
                connection.execute(
                    "INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        mnxr_id,
                        equation,
                        human_equation,
                        int(balanced),
                        int(transport),
                        reference,
                        parse_status,
                    ),
                )
                connection.executemany(
                    "INSERT INTO reaction_terms VALUES (?, ?, ?, ?, ?, ?)",
                    ((mnxr_id, *term) for term in terms),
                )
                required_compounds.update(term[2] for term in terms)
                reaction_count += 1
                term_count += len(terms)

            chemical_count = 0
            found_compounds: set[str] = set()
            for row in _data_rows(normalized_sources["chem_prop.tsv"]):
                if len(row) < 9 or row[0] not in required_compounds:
                    continue
                charge = None
                mass = None
                try:
                    charge = int(row[3]) if row[3].strip() else None
                except ValueError:
                    pass
                try:
                    mass = float(row[4]) if row[4].strip() else None
                except ValueError:
                    pass
                connection.execute(
                    "INSERT INTO chemicals VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (
                        row[0],
                        row[1],
                        row[2],
                        charge,
                        mass,
                        row[5],
                        row[6],
                        row[7],
                        row[8],
                    ),
                )
                found_compounds.add(row[0])
                chemical_count += 1

            chemical_xref_count = 0
            for row in _data_rows(normalized_sources["chem_xref.tsv"]):
                if len(row) < 2 or row[1] not in required_compounds:
                    continue
                padded = [*row, "", ""]
                connection.execute(
                    "INSERT OR IGNORE INTO chemical_xrefs VALUES (?, ?, ?, ?)",
                    (padded[0], padded[1], padded[2], padded[3]),
                )
                chemical_xref_count += 1

            reaction_xref_count = 0
            for row in _data_rows(normalized_sources["reac_xref.tsv"]):
                if len(row) < 2 or row[1] not in requested_reactions:
                    continue
                connection.execute(
                    "INSERT OR IGNORE INTO reaction_xrefs VALUES (?, ?)",
                    (row[0], row[1]),
                )
                reaction_xref_count += 1

            connection.executemany(
                "INSERT INTO metadata VALUES (?, ?)",
                (
                    ("schema_version", MNXREF_INDEX_SCHEMA),
                    ("mnxref_version", MNXREF_VERSION),
                    ("mnxref_release", MNXREF_RELEASE),
                    ("rr02_sha256", _sha256_file(resolved_rules)),
                ),
            )
            connection.commit()
        finally:
            connection.close()

        index_path = resolved_output / MNXREF_INDEX_FILE_NAME
        temporary_db.replace(index_path)
        temporary_db = None
        sources = {}
        for name, path in normalized_sources.items():
            metadata = dict((source_metadata or {}).get(name, {}))
            sources[name] = {
                **metadata,
                "size_bytes": path.stat().st_size,
                "sha256": _sha256_file(path),
                "md5": _md5_file(path),
            }
        manifest = {
            "schema_version": MNXREF_MANIFEST_SCHEMA,
            "mnxref_version": MNXREF_VERSION,
            "mnxref_release": MNXREF_RELEASE,
            "installed_at": datetime.now(UTC).isoformat(),
            "rr02_path": str(resolved_rules),
            "rr02_sha256": _sha256_file(resolved_rules),
            "index_path": str(index_path),
            "index_sha256": _sha256_file(index_path),
            "sources": sources,
            "counts": {
                "rule_template_links": len(rule_rows),
                "requested_reactions": len(requested_reactions),
                "indexed_reactions": reaction_count,
                "reaction_terms": term_count,
                "requested_chemicals": len(required_compounds),
                "indexed_chemicals": chemical_count,
                "missing_chemicals": len(required_compounds - found_compounds),
                "chemical_xrefs": chemical_xref_count,
                "reaction_xrefs": reaction_xref_count,
            },
        }
        _atomic_write_json(resolved_output / MNXREF_MANIFEST_FILE_NAME, manifest)
        return manifest
    finally:
        if temporary_db is not None:
            temporary_db.unlink(missing_ok=True)


def install_mnxref_subset(
    *,
    rules_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    downloader: DownloadFile = _download_file,
    force: bool = False,
) -> dict[str, Any]:
    """Download pinned files, verify official MD5 values, and build the subset."""

    resolved_rules = Path(rules_path or default_rules_path()).expanduser().resolve()
    resolved_output = Path(output_dir or default_install_dir()).expanduser().resolve()
    manifest_path = resolved_output / MNXREF_MANIFEST_FILE_NAME
    resolved_output.mkdir(parents=True, exist_ok=True)
    for stale_path in resolved_output.iterdir():
        if (
            stale_path.name.startswith(".mnxref-install-")
            and stale_path.is_dir()
            and stale_path.resolve().parent == resolved_output
        ):
            shutil.rmtree(stale_path, ignore_errors=True)
    if not force and manifest_path.is_file():
        return validate_mnxref_index(resolved_output, resolved_rules)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=".mnxref-install-", dir=resolved_output)
    ).resolve()
    try:
        source_paths: dict[str, Path] = {}
        source_metadata: dict[str, dict[str, Any]] = {}
        for name in MNXREF_FILES:
            url = f"{MNXREF_BASE_URL}/{name}"
            md5_url = f"{url}.md5"
            source_path = staging_dir / name
            md5_path = staging_dir / f"{name}.md5"
            print(f"[INFO] downloading MNXref {MNXREF_VERSION}: {name}")
            downloader(url, source_path)
            downloader(md5_url, md5_path)
            official_md5 = _read_official_md5(md5_path)
            observed_md5 = _md5_file(source_path)
            if observed_md5 != official_md5:
                raise MnxrefIndexError(
                    f"MD5 mismatch for {name}: {observed_md5} != {official_md5}"
                )
            source_paths[name] = source_path
            source_metadata[name] = {
                "url": url,
                "official_md5_url": md5_url,
                "official_md5": official_md5,
            }
        return build_mnxref_subset(
            rules_path=resolved_rules,
            source_paths=source_paths,
            output_dir=resolved_output,
            source_metadata=source_metadata,
        )
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)


def validate_mnxref_index(
    output_dir: str | Path | None = None,
    rules_path: str | Path | None = None,
) -> dict[str, Any]:
    resolved_output = Path(output_dir or default_install_dir()).expanduser().resolve()
    resolved_rules = Path(rules_path or default_rules_path()).expanduser().resolve()
    manifest_path = resolved_output / MNXREF_MANIFEST_FILE_NAME
    index_path = resolved_output / MNXREF_INDEX_FILE_NAME
    install_command = "uv run python -m src.pathway_analyze.retropath_mnxref install"
    if not manifest_path.is_file() or not index_path.is_file():
        raise MnxrefIndexError(
            f"MNXref {MNXREF_VERSION} index is not installed; run: {install_command}"
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise MnxrefIndexError(f"invalid MNXref manifest: {manifest_path}") from exc
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != MNXREF_MANIFEST_SCHEMA
        or manifest.get("mnxref_version") != MNXREF_VERSION
    ):
        raise MnxrefIndexError(f"unsupported MNXref manifest: {manifest_path}")
    if not resolved_rules.is_file():
        raise MnxrefIndexError(f"RR02 rules file not found: {resolved_rules}")
    if manifest.get("rr02_sha256") != _sha256_file(resolved_rules):
        raise MnxrefIndexError("MNXref subset was built for a different RR02 file")
    if manifest.get("index_sha256") != _sha256_file(index_path):
        raise MnxrefIndexError("MNXref SQLite index checksum mismatch")
    connection = sqlite3.connect(f"file:{index_path.as_posix()}?mode=ro", uri=True)
    try:
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    finally:
        connection.close()
    if (
        metadata.get("schema_version") != MNXREF_INDEX_SCHEMA
        or metadata.get("mnxref_version") != MNXREF_VERSION
        or metadata.get("rr02_sha256") != manifest.get("rr02_sha256")
    ):
        raise MnxrefIndexError("MNXref SQLite metadata does not match its manifest")
    return manifest


class MnxrefIndex:
    """Read-only query facade over one validated subset index."""

    def __init__(
        self,
        output_dir: str | Path | None = None,
        rules_path: str | Path | None = None,
    ) -> None:
        self.output_dir = Path(output_dir or default_install_dir()).expanduser().resolve()
        self.rules_path = Path(rules_path or default_rules_path()).expanduser().resolve()
        self.manifest = validate_mnxref_index(self.output_dir, self.rules_path)
        self.index_path = self.output_dir / MNXREF_INDEX_FILE_NAME
        self._connection = sqlite3.connect(
            f"file:{self.index_path.as_posix()}?mode=ro",
            uri=True,
        )
        self._connection.row_factory = sqlite3.Row

    def close(self) -> None:
        self._connection.close()

    def __enter__(self) -> "MnxrefIndex":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def templates_for_rules(
        self,
        rule_ids: Sequence[str],
    ) -> tuple[MnxrefReactionTemplate, ...]:
        normalized = sorted({str(value).strip() for value in rule_ids if str(value).strip()})
        if not normalized:
            return tuple()
        placeholders = ",".join("?" for _ in normalized)
        rows = self._connection.execute(
            f"""
            SELECT rt.*, r.equation, r.balanced, r.transport, r.reference,
                   r.parse_status
            FROM rule_templates rt
            LEFT JOIN reactions r ON r.mnxr_id = rt.mnxr_id
            WHERE rt.rule_id IN ({placeholders})
            ORDER BY rt.rule_id, rt.mnxr_id, rt.main_mnxm_id
            """,
            normalized,
        ).fetchall()
        templates: list[MnxrefReactionTemplate] = []
        for row in rows:
            terms = tuple(
                MnxrefReactionTerm(
                    side=item["side"],
                    coefficient=float(item["coefficient"]),
                    mnxm_id=item["mnxm_id"],
                    compartment=item["compartment"],
                    ordinal=int(item["ordinal"]),
                )
                for item in self._connection.execute(
                    """
                    SELECT side, coefficient, mnxm_id, compartment, ordinal
                    FROM reaction_terms WHERE mnxr_id = ?
                    ORDER BY CASE side WHEN 'left' THEN 0 ELSE 1 END, ordinal
                    """,
                    (row["mnxr_id"],),
                )
            )
            reaction_xrefs = tuple(
                item[0]
                for item in self._connection.execute(
                    "SELECT xref FROM reaction_xrefs WHERE mnxr_id = ? ORDER BY xref",
                    (row["mnxr_id"],),
                )
            )
            templates.append(
                MnxrefReactionTemplate(
                    rule_id=row["rule_id"],
                    mnxr_id=row["mnxr_id"],
                    main_mnxm_id=row["main_mnxm_id"],
                    reaction_direction=row["reaction_direction"],
                    rule_relative_direction=row["rule_relative_direction"],
                    rule_usage=row["rule_usage"],
                    equation=row["equation"] or "",
                    balanced=bool(row["balanced"] or False),
                    transport=bool(row["transport"] or False),
                    reference=row["reference"] or "",
                    parse_status=row["parse_status"] or "reaction_missing",
                    terms=terms,
                    reaction_xrefs=reaction_xrefs,
                )
            )
        return tuple(templates)

    def chemicals(
        self,
        mnxm_ids: Iterable[str],
    ) -> dict[str, MnxrefChemical]:
        normalized = sorted({str(value).strip() for value in mnxm_ids if str(value).strip()})
        if not normalized:
            return {}
        placeholders = ",".join("?" for _ in normalized)
        rows = self._connection.execute(
            f"SELECT * FROM chemicals WHERE mnxm_id IN ({placeholders})",
            normalized,
        ).fetchall()
        result: dict[str, MnxrefChemical] = {}
        for row in rows:
            xrefs = tuple(
                item[0]
                for item in self._connection.execute(
                    "SELECT xref FROM chemical_xrefs WHERE mnxm_id = ? ORDER BY xref",
                    (row["mnxm_id"],),
                )
            )
            result[row["mnxm_id"]] = MnxrefChemical(
                mnxm_id=row["mnxm_id"],
                name=row["name"],
                formula=row["formula"],
                charge=row["charge"],
                mass=row["mass"],
                inchi=row["inchi"],
                smiles=row["smiles"],
                reference=row["reference"],
                inchikey=row["inchikey"],
                xrefs=xrefs,
            )
        return result


def _main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Install the RR02 MNXref 3.0 subset")
    subparsers = parser.add_subparsers(dest="command", required=True)
    install = subparsers.add_parser("install")
    install.add_argument("--rules", type=Path, default=default_rules_path())
    install.add_argument("--output", type=Path, default=default_install_dir())
    install.add_argument("--force", action="store_true")
    status = subparsers.add_parser("status")
    status.add_argument("--rules", type=Path, default=default_rules_path())
    status.add_argument("--output", type=Path, default=default_install_dir())
    args = parser.parse_args(argv)
    if args.command == "install":
        result = install_mnxref_subset(
            rules_path=args.rules,
            output_dir=args.output,
            force=args.force,
        )
    else:
        result = validate_mnxref_index(args.output, args.rules)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())


__all__ = [
    "MNXREF_INDEX_FILE_NAME",
    "MNXREF_MANIFEST_FILE_NAME",
    "MNXREF_VERSION",
    "MnxrefChemical",
    "MnxrefIndex",
    "MnxrefIndexError",
    "MnxrefReactionTemplate",
    "MnxrefReactionTerm",
    "build_mnxref_subset",
    "default_install_dir",
    "default_rules_path",
    "install_mnxref_subset",
    "validate_mnxref_index",
]
