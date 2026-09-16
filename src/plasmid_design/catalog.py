"""Validated, offline catalogue of the curated BASIC SEVA modules."""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

_DNA = frozenset("ACGT")
_EXPECTED_COLUMNS = {
    "resistance": ("id", "name", "antibiotic", "resistance_gene", "sequence"),
    "replication": ("id", "name", "host_range", "copy_number", "sequence"),
    "terminator": ("id", "name", "sequence"),
}
_NOTES = {
    "basic_seva_gm": (
        "原始庆大霉素模块。",
        "该版本不是作者针对 p15A/pBR322 组合公布的 66.11/69.11 变体；如需该组合的已公布变体，请选择 basic_seva_gm_11。",
    ),
    "basic_seva_gm_11": (
        "66.11/69.11 庆大霉素变体，作者公布用于 p15A/pBR322 组合。",
        "其与其他复制模块的组合尚未验证。",
    ),
    "basic_seva_psc101_pkd46_ts": (
        "温敏 pKD46/pSC101 模块：通常在 30°C 维持，在 37°C 不稳定。",
    ),
}


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


@dataclass(frozen=True, slots=True)
class Module:
    id: str
    name: str
    type: str
    sequence: str
    antibiotic: str | None = None
    resistance_gene: str | None = None
    host_range: str | None = None
    copy_number: str | None = None
    role: str | None = None
    notes: tuple[str, ...] = ()
    features: tuple[Mapping[str, Any], ...] = ()

    @property
    def length_bp(self) -> int:
        return len(self.sequence)

    @property
    def gc_percent(self) -> float:
        return (
            100 * (self.sequence.count("G") + self.sequence.count("C")) / self.length_bp
        )

    def to_payload(self) -> dict[str, Any]:
        payload = {
            "id": self.id,
            "name": self.name,
            "type": self.type,
            "sequence": self.sequence,
            "notes": list(self.notes),
            "length_bp": self.length_bp,
            "gc_percent": self.gc_percent,
        }
        if self.type == "resistance":
            payload.update(
                antibiotic=self.antibiotic, resistance_gene=self.resistance_gene
            )
        elif self.type == "replication":
            payload.update(host_range=self.host_range, copy_number=self.copy_number)
        elif self.type == "terminator":
            payload["role"] = self.role
        return payload


class ModuleCatalog:
    """Load modules from a data directory and validate their curated annotations."""

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)
        package_dir = Path(__file__).parent
        self._feature_snapshot = self._read_json(package_dir / "source_features.json")
        if self._feature_snapshot.get("_schema_version") != 2:
            raise ValueError("unsupported source feature schema version")
        self.scaffold = self._load_scaffold(package_dir / "scaffold.json")
        self._resistance = MappingProxyType(self._load_type("resistance"))
        self._replication = MappingProxyType(self._load_type("replication"))
        self._terminator = MappingProxyType(self._load_type("terminator"))
        self._validate_derivations()
        self._fingerprint = self._compute_fingerprint()

    @property
    def fingerprint(self) -> str:
        return self._fingerprint

    @property
    def resistance(self) -> Mapping[str, Module]:
        return self._resistance

    @property
    def replication(self) -> Mapping[str, Module]:
        return self._replication

    @property
    def terminator(self) -> Mapping[str, Module]:
        return self._terminator

    def get_resistance(self, module_id: str) -> Module:
        return self._get(self.resistance, module_id, "resistance")

    def get_replication(self, module_id: str) -> Module:
        return self._get(self.replication, module_id, "replication")

    def get_terminator(self, module_id: str) -> Module:
        return self._get(self.terminator, module_id, "terminator")

    @staticmethod
    def _get(modules: Mapping[str, Module], module_id: str, module_type: str) -> Module:
        try:
            return modules[module_id]
        except KeyError as exc:
            raise KeyError(f"unknown {module_type} module: {module_id}") from exc

    @staticmethod
    def _read_json(path: Path) -> dict[str, Any]:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"invalid catalog metadata {path}: {exc}") from exc

    def _load_scaffold(self, path: Path) -> Mapping[str, Any]:
        scaffold = self._read_json(path)
        if scaffold.get("_schema_version") != 2:
            raise ValueError("unsupported scaffold schema version")
        provenance = scaffold.get("provenance")
        if not isinstance(provenance, dict):
            raise ValueError("scaffold provenance is required")  # noqa: TRY004 - invalid persisted JSON data
        for key in (
            "resistance_to_replication",
            "replication_to_t1",
            "landing_pad_spacer",
        ):
            sequence = scaffold.get(key)
            details = provenance.get(key)
            if not isinstance(sequence, str) or not sequence or set(sequence) - _DNA:
                raise ValueError(f"invalid scaffold sequence: {key}")
            if not isinstance(details, dict) or hashlib.sha256(
                sequence.encode()
            ).hexdigest() != details.get("sha256"):
                raise ValueError(f"scaffold source hash mismatch: {key}")
            self._validate_provenance(key, sequence, {**provenance, **details})
        return _freeze(scaffold)

    def _load_type(self, module_type: str) -> dict[str, Module]:
        path = self.data_dir / f"{module_type}_modules.csv"
        try:
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                if tuple(reader.fieldnames or ()) != _EXPECTED_COLUMNS[module_type]:
                    raise ValueError(
                        f"unexpected columns in {path.name}: {reader.fieldnames}"
                    )
                rows = list(reader)
        except OSError as exc:
            raise ValueError(f"unable to load {path}: {exc}") from exc
        modules: dict[str, Module] = {}
        for row in rows:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed CSV row in {path.name}")
            for existing_type in ("resistance", "replication", "terminator"):
                existing = getattr(self, f"_{existing_type}", {})
                if row["id"] in existing:
                    raise ValueError(
                        f"module identifiers must be unique across types: {row['id']}"
                    )
            module = self._module_from_row(module_type, row)
            if module.id in modules:
                raise ValueError(f"duplicate {module_type} module id: {module.id}")
            modules[module.id] = module
        if not modules:
            raise ValueError(f"no {module_type} modules found")
        return modules

    def _module_from_row(
        self, module_type: str, row: Mapping[str, str | None]
    ) -> Module:
        module_id = row["id"] or ""
        name = row["name"] or ""
        sequence = row["sequence"] or ""
        required = _EXPECTED_COLUMNS[module_type]
        if (
            any(not (row.get(field) or "").strip() for field in required)
            or set(sequence) - _DNA
        ):
            raise ValueError(
                f"invalid {module_type} module row: {module_id or '<missing id>'}"
            )
        snapshot = self._feature_snapshot.get(module_id)
        if not isinstance(snapshot, dict) or hashlib.sha256(
            sequence.encode()
        ).hexdigest() != snapshot.get("sha256"):
            raise ValueError(
                f"source feature snapshot mismatch for module: {module_id}"
            )
        if snapshot.get("type") != module_type:
            raise ValueError(f"source module type mismatch: {module_id}")
        self._validate_provenance(module_id, sequence, snapshot.get("provenance"))
        features = snapshot.get("features")
        if not isinstance(features, list):
            raise ValueError(f"invalid feature snapshot for module: {module_id}")  # noqa: TRY004 - catalog data validation
        for feature in features:
            self._validate_feature(module_id, sequence, feature)
        common = {
            "id": module_id,
            "name": name,
            "type": module_type,
            "sequence": sequence,
            "notes": _NOTES.get(module_id, ()),
            "features": tuple(_freeze(feature) for feature in features),
        }
        if module_type == "resistance":
            return Module(
                **common,
                antibiotic=row["antibiotic"],
                resistance_gene=row["resistance_gene"],
            )
        if module_type == "terminator":
            role = snapshot.get("role")
            if role not in ("t0", "t1"):
                raise ValueError(f"invalid terminator role for module: {module_id}")
            return Module(**common, role=role)
        return Module(
            **common, host_range=row["host_range"], copy_number=row["copy_number"]
        )

    def _validate_provenance(
        self, module_id: str, sequence: str, provenance: Any
    ) -> None:
        if not isinstance(provenance, dict):
            raise ValueError(f"source provenance is required: {module_id}")  # noqa: TRY004 - invalid persisted JSON data
        sources = self._feature_snapshot.get("_sources")
        source_key = provenance.get("source_key")
        source = (
            sources.get(source_key)
            if isinstance(sources, dict) and isinstance(source_key, str)
            else None
        )
        if not isinstance(source, dict):
            raise ValueError(f"unknown pinned source for module: {module_id}")  # noqa: TRY004 - invalid persisted JSON data
        for key, length in (("source_sha256", 64), ("source_version", 40)):
            value = provenance.get(key)
            if (
                not isinstance(value, str)
                or len(value) != length
                or set(value) - frozenset("0123456789abcdef")
                or value != source.get(key)
            ):
                raise ValueError(f"invalid pinned source {key}: {module_id}")
        url = provenance.get("source_url")
        record = provenance.get("source_record")
        if (
            not isinstance(url, str)
            or url != source.get("source_url")
            or provenance["source_version"] not in url
            or not isinstance(record, str)
            or not record.strip()
        ):
            raise ValueError(f"invalid pinned source reference: {module_id}")
        start, end = provenance.get("start_1based"), provenance.get("end_1based")
        if (
            type(start) is not int
            or type(end) is not int
            or start < 1
            or end - start + 1 != len(sequence)
        ):
            raise ValueError(f"invalid source coordinates for module: {module_id}")

    def _validate_derivations(self) -> None:
        """Verify recorded extraction history without transforming CSV sequences."""
        for module in self.resistance.values():
            snapshot = self._feature_snapshot[module.id]
            derivation = snapshot.get("derivation")
            if derivation is None:
                continue
            if not isinstance(derivation, dict):
                raise ValueError(f"invalid source derivation: {module.id}")  # noqa: TRY004 - invalid persisted JSON data
            prefix = self.terminator.get(derivation.get("removed_prefix_module_id"))
            provenance = snapshot["provenance"]
            if (
                prefix is None
                or prefix.role != "t0"
                or derivation.get("operation") != "remove_verified_prefix_once"
                or derivation.get("removed_prefix_bp") != prefix.length_bp
                or derivation.get("removed_prefix_sha256")
                != hashlib.sha256(prefix.sequence.encode()).hexdigest()
                or derivation.get("parent_sha256")
                != hashlib.sha256(
                    (prefix.sequence + module.sequence).encode()
                ).hexdigest()
                or derivation.get("parent_start_1based")
                != provenance["start_1based"] - prefix.length_bp
                or derivation.get("parent_end_1based") != provenance["end_1based"]
            ):
                raise ValueError(f"source derivation mismatch: {module.id}")

    @staticmethod
    def _validate_feature(module_id: str, sequence: str, feature: Any) -> None:
        if (
            not isinstance(feature, dict)
            or not isinstance(feature.get("type"), str)
            or not isinstance(feature.get("parts"), list)
        ):
            raise ValueError(f"invalid feature metadata for module: {module_id}")  # noqa: TRY004 - catalog data validation
        for part in feature["parts"]:
            if not isinstance(part, dict) or not all(
                isinstance(part.get(key), int) for key in ("start", "end")
            ):
                raise ValueError(f"invalid feature range for module: {module_id}")
            if not 0 <= part["start"] < part["end"] <= len(sequence):
                raise ValueError(f"feature outside module sequence: {module_id}")

    def _compute_fingerprint(self) -> str:
        digest = hashlib.sha256()
        for path in (
            self.data_dir / "resistance_modules.csv",
            self.data_dir / "replication_modules.csv",
            self.data_dir / "terminator_modules.csv",
            Path(__file__).parent / "scaffold.json",
            Path(__file__).parent / "source_features.json",
        ):
            digest.update(path.name.encode())
            digest.update(path.read_bytes())
        return digest.hexdigest()
