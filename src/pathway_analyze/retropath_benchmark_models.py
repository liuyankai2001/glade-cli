"""Versioned data contracts for the P11.1 RetroPath benchmark.

The benchmark deliberately separates the information visible to the search
pipeline from the KEGG/MNXref gold standard used only during scoring.  This
module has no network or Docker dependency and is shared by the runner,
report generator, and offline tests.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


BENCHMARK_DATASET_SCHEMA = "retropath_benchmark_cases.v1"
BENCHMARK_RUN_SCHEMA = "retropath_benchmark_run.v1"
BENCHMARK_TASK_SCHEMA = "retropath_benchmark_task.v1"

KEGG_COMPOUND_RE = re.compile(r"^C\d{5}$")
KEGG_REACTION_RE = re.compile(r"^R\d{5}$")
MNXR_RE = re.compile(r"^MNXR\d+$")
EC_RE = re.compile(r"^\d+\.(?:\d+|-)\.(?:\d+|-|[A-Za-z]+)\.(?:\d+|-|[A-Za-z]+)$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    text = json.dumps(
        value,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be an object")
    return value


def _sequence(value: Any, field: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{field} must be an array")
    return value


def _text(value: Any, field: str) -> str:
    result = str(value or "").strip()
    if not result:
        raise ValueError(f"{field} must be a non-empty string")
    return result


def _unique_strings(
    value: Any,
    field: str,
    *,
    pattern: re.Pattern[str] | None = None,
    allow_empty: bool = False,
) -> tuple[str, ...]:
    raw = _sequence(value, field)
    values = tuple(dict.fromkeys(str(item).strip() for item in raw if str(item).strip()))
    if not values and not allow_empty:
        raise ValueError(f"{field} must contain at least one value")
    if pattern is not None:
        invalid = [item for item in values if pattern.fullmatch(item) is None]
        if invalid:
            raise ValueError(f"{field} contains invalid values: {invalid}")
    return values


def _positive_int(value: Any, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a positive integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a positive integer") from exc
    if result < 1:
        raise ValueError(f"{field} must be a positive integer")
    return result


@dataclass(frozen=True)
class BenchmarkResource:
    name: str
    path: Path
    sha256: str
    compatible_sha256: tuple[str, ...] = ()

    @classmethod
    def from_mapping(
        cls,
        name: str,
        value: Any,
        *,
        root: Path,
    ) -> "BenchmarkResource":
        payload = _mapping(value, f"resources.{name}")
        raw_path = Path(_text(payload.get("path"), f"resources.{name}.path"))
        path = raw_path if raw_path.is_absolute() else root / raw_path
        digest = _text(payload.get("sha256"), f"resources.{name}.sha256").lower()
        if SHA256_RE.fullmatch(digest) is None:
            raise ValueError(f"resources.{name}.sha256 must be a SHA-256 hex digest")
        compatible = _unique_strings(
            payload.get("compatible_sha256", []),
            f"resources.{name}.compatible_sha256",
            pattern=SHA256_RE,
            allow_empty=True,
        )
        return cls(
            name=name,
            path=path.resolve(),
            sha256=digest,
            compatible_sha256=compatible,
        )

    def verify(self) -> None:
        if not self.path.is_file():
            raise FileNotFoundError(f"benchmark resource is missing: {self.path}")
        observed = sha256_file(self.path)
        if observed != self.sha256:
            raise ValueError(
                f"benchmark resource checksum mismatch for {self.name}: "
                f"expected {self.sha256}, observed {observed}"
            )

    def to_dict(self) -> dict[str, str]:
        return {
            "name": self.name,
            "path": str(self.path),
            "sha256": self.sha256,
            "compatible_sha256": list(self.compatible_sha256),
        }


@dataclass(frozen=True)
class BenchmarkGoldStep:
    reaction_id: str
    mnxr_id: str
    ec_numbers: tuple[str, ...]
    uniprot_accessions: tuple[str, ...]

    @classmethod
    def from_mapping(cls, value: Any, *, field: str) -> "BenchmarkGoldStep":
        payload = _mapping(value, field)
        reaction_id = _text(payload.get("reaction_id"), f"{field}.reaction_id").upper()
        mnxr_id = _text(payload.get("mnxr_id"), f"{field}.mnxr_id").upper()
        if KEGG_REACTION_RE.fullmatch(reaction_id) is None:
            raise ValueError(f"{field}.reaction_id must be Rxxxxx")
        if MNXR_RE.fullmatch(mnxr_id) is None:
            raise ValueError(f"{field}.mnxr_id must be MNXR followed by digits")
        ecs = _unique_strings(payload.get("ec_numbers"), f"{field}.ec_numbers")
        invalid_ecs = [item for item in ecs if EC_RE.fullmatch(item) is None]
        if invalid_ecs:
            raise ValueError(f"{field}.ec_numbers contains invalid values: {invalid_ecs}")
        accessions = _unique_strings(
            payload.get("uniprot_accessions", []),
            f"{field}.uniprot_accessions",
            allow_empty=True,
        )
        return cls(reaction_id, mnxr_id, ecs, accessions)

    def to_dict(self) -> dict[str, Any]:
        return {
            "reaction_id": self.reaction_id,
            "mnxr_id": self.mnxr_id,
            "ec_numbers": list(self.ec_numbers),
            "uniprot_accessions": list(self.uniprot_accessions),
        }


@dataclass(frozen=True)
class BenchmarkCase:
    case_id: str
    ec_class: int
    target_kegg_id: str
    target_mnxm_id: str
    controlled_sink_kegg_ids: tuple[str, ...]
    controlled_sink_mnxm_ids: tuple[str, ...]
    gold_steps: tuple[BenchmarkGoldStep, ...]
    notes: str

    @classmethod
    def from_mapping(cls, value: Any, *, index: int) -> "BenchmarkCase":
        field = f"cases[{index}]"
        payload = _mapping(value, field)
        case_id = _text(payload.get("case_id"), f"{field}.case_id")
        if re.fullmatch(r"[a-z0-9][a-z0-9_-]*", case_id) is None:
            raise ValueError(f"{field}.case_id must use lowercase letters, digits, _ or -")
        ec_class = _positive_int(payload.get("ec_class"), f"{field}.ec_class")
        if ec_class not in range(1, 7):
            raise ValueError(f"{field}.ec_class must be between 1 and 6")
        search = _mapping(payload.get("search"), f"{field}.search")
        target = _text(
            search.get("target_kegg_id"), f"{field}.search.target_kegg_id"
        ).upper()
        target_mnxm = _text(
            search.get("target_mnxm_id"), f"{field}.search.target_mnxm_id"
        ).upper()
        if KEGG_COMPOUND_RE.fullmatch(target) is None:
            raise ValueError(f"{field}.search.target_kegg_id must be Cxxxxx")
        if re.fullmatch(r"^MNXM\d+$", target_mnxm) is None:
            raise ValueError(f"{field}.search.target_mnxm_id must be MNXM followed by digits")
        sinks = _unique_strings(
            search.get("controlled_sink_kegg_ids"),
            f"{field}.search.controlled_sink_kegg_ids",
            pattern=KEGG_COMPOUND_RE,
        )
        sink_mnxm = _unique_strings(
            search.get("controlled_sink_mnxm_ids"),
            f"{field}.search.controlled_sink_mnxm_ids",
            pattern=re.compile(r"^MNXM\d+$"),
        )
        if target in sinks:
            raise ValueError(f"{field} leaks the target into the controlled sink")
        gold = _mapping(payload.get("gold"), f"{field}.gold")
        steps = tuple(
            BenchmarkGoldStep.from_mapping(item, field=f"{field}.gold.steps[{offset}]")
            for offset, item in enumerate(
                _sequence(gold.get("steps"), f"{field}.gold.steps")
            )
        )
        if not steps:
            raise ValueError(f"{field}.gold.steps must not be empty")
        if any(int(step.ec_numbers[0].split(".", 1)[0]) != ec_class for step in steps):
            raise ValueError(f"{field}.ec_class does not match its gold EC numbers")
        return cls(
            case_id=case_id,
            ec_class=ec_class,
            target_kegg_id=target,
            target_mnxm_id=target_mnxm,
            controlled_sink_kegg_ids=sinks,
            controlled_sink_mnxm_ids=sink_mnxm,
            gold_steps=steps,
            notes=str(payload.get("notes") or "").strip(),
        )

    @property
    def gold_reaction_ids(self) -> tuple[str, ...]:
        return tuple(step.reaction_id for step in self.gold_steps)

    @property
    def gold_mnxr_ids(self) -> tuple[str, ...]:
        return tuple(step.mnxr_id for step in self.gold_steps)

    @property
    def gold_ec_numbers(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(ec for step in self.gold_steps for ec in step.ec_numbers))

    @property
    def gold_uniprot_accessions(self) -> tuple[str, ...]:
        return tuple(
            dict.fromkeys(
                accession
                for step in self.gold_steps
                for accession in step.uniprot_accessions
            )
        )

    def search_dict(self) -> dict[str, Any]:
        """Return only information that may be exposed to RetroPath."""

        return {
            "target_kegg_id": self.target_kegg_id,
            "target_mnxm_id": self.target_mnxm_id,
            "controlled_sink_kegg_ids": list(self.controlled_sink_kegg_ids),
            "controlled_sink_mnxm_ids": list(self.controlled_sink_mnxm_ids),
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "ec_class": self.ec_class,
            "search": self.search_dict(),
            "gold": {"steps": [step.to_dict() for step in self.gold_steps]},
            "notes": self.notes,
        }


@dataclass(frozen=True)
class BenchmarkDefaults:
    max_steps: int
    topx: int
    dmin: int
    dmax: int
    max_candidates: int
    top_k: tuple[int, ...]
    run_fva: bool
    request_timeout_seconds: int
    get_attempts: int
    wait_timeout_seconds: int

    @classmethod
    def from_mapping(cls, value: Any) -> "BenchmarkDefaults":
        payload = _mapping(value, "defaults")
        top_k = tuple(
            sorted(
                set(
                    _positive_int(item, "defaults.top_k item")
                    for item in _sequence(payload.get("top_k"), "defaults.top_k")
                )
            )
        )
        dmin = _positive_int(payload.get("dmin"), "defaults.dmin")
        dmax = _positive_int(payload.get("dmax"), "defaults.dmax")
        if dmax < dmin:
            raise ValueError("defaults.dmax must be greater than or equal to dmin")
        max_candidates = _positive_int(
            payload.get("max_candidates"), "defaults.max_candidates"
        )
        if top_k[-1] > max_candidates:
            raise ValueError("defaults.top_k cannot exceed max_candidates")
        run_fva = payload.get("run_fva", True)
        if not isinstance(run_fva, bool):
            raise ValueError("defaults.run_fva must be boolean")
        return cls(
            max_steps=_positive_int(payload.get("max_steps"), "defaults.max_steps"),
            topx=_positive_int(payload.get("topx"), "defaults.topx"),
            dmin=dmin,
            dmax=dmax,
            max_candidates=max_candidates,
            top_k=top_k,
            run_fva=run_fva,
            request_timeout_seconds=_positive_int(
                payload.get("request_timeout_seconds"),
                "defaults.request_timeout_seconds",
            ),
            get_attempts=_positive_int(
                payload.get("get_attempts"), "defaults.get_attempts"
            ),
            wait_timeout_seconds=_positive_int(
                payload.get("wait_timeout_seconds"),
                "defaults.wait_timeout_seconds",
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "max_steps": self.max_steps,
            "topx": self.topx,
            "dmin": self.dmin,
            "dmax": self.dmax,
            "max_candidates": self.max_candidates,
            "top_k": list(self.top_k),
            "run_fva": self.run_fva,
            "request_timeout_seconds": self.request_timeout_seconds,
            "get_attempts": self.get_attempts,
            "wait_timeout_seconds": self.wait_timeout_seconds,
        }


@dataclass(frozen=True)
class BenchmarkDataset:
    path: Path
    benchmark_id: str
    description: str
    resources: Mapping[str, BenchmarkResource]
    defaults: BenchmarkDefaults
    cases: tuple[BenchmarkCase, ...]
    dataset_sha256: str

    def case_by_id(self, case_id: str) -> BenchmarkCase:
        for case in self.cases:
            if case.case_id == case_id:
                return case
        raise KeyError(case_id)

    def verify_resources(self) -> None:
        for resource in self.resources.values():
            resource.verify()

    def validate_pilot_shape(self) -> None:
        if len(self.cases) != 12:
            raise ValueError("P11.1 pilot dataset must contain exactly 12 cases")
        counts = {ec_class: 0 for ec_class in range(1, 7)}
        for case in self.cases:
            counts[case.ec_class] += 1
        if any(count != 2 for count in counts.values()):
            raise ValueError(
                "P11.1 pilot dataset must contain exactly two cases per EC class"
            )

    def validate_a0_snapshot(self) -> None:
        snapshot = self.resources["a0_snapshot"].path
        with snapshot.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "kegg_id" not in reader.fieldnames:
                raise ValueError("A0 snapshot must contain a kegg_id column")
            a0 = {
                str(row.get("kegg_id") or "").strip().upper()
                for row in reader
                if KEGG_COMPOUND_RE.fullmatch(
                    str(row.get("kegg_id") or "").strip().upper()
                )
            }
        if not a0:
            raise ValueError("A0 snapshot contains no KEGG compounds")
        for case in self.cases:
            if case.target_kegg_id in a0:
                raise ValueError(f"{case.case_id} target is already present in A0")
            missing = sorted(set(case.controlled_sink_kegg_ids) - a0)
            if missing:
                raise ValueError(
                    f"{case.case_id} controlled sinks are absent from A0: {missing}"
                )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": BENCHMARK_DATASET_SCHEMA,
            "benchmark_id": self.benchmark_id,
            "description": self.description,
            "resources": {
                key: value.to_dict() for key, value in sorted(self.resources.items())
            },
            "defaults": self.defaults.to_dict(),
            "cases": [case.to_dict() for case in self.cases],
            "dataset_sha256": self.dataset_sha256,
        }


REQUIRED_RESOURCES = (
    "model",
    "medium",
    "rr02",
    "mnxref_index",
    "a0_snapshot",
)


def load_benchmark_dataset(
    path: str | Path,
    *,
    root: str | Path | None = None,
    verify_resources: bool = True,
    require_pilot_shape: bool = True,
) -> BenchmarkDataset:
    resolved_path = Path(path).expanduser().resolve()
    try:
        payload = json.loads(resolved_path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid benchmark dataset: {resolved_path}") from exc
    document = _mapping(payload, "dataset")
    if document.get("schema_version") != BENCHMARK_DATASET_SCHEMA:
        raise ValueError(
            f"unsupported benchmark dataset schema: {document.get('schema_version')!r}"
        )
    resolved_root = Path(root or repository_root()).expanduser().resolve()
    resource_payload = _mapping(document.get("resources"), "resources")
    missing = [name for name in REQUIRED_RESOURCES if name not in resource_payload]
    if missing:
        raise ValueError(f"benchmark dataset is missing resources: {missing}")
    resources = {
        name: BenchmarkResource.from_mapping(
            name,
            resource_payload[name],
            root=resolved_root,
        )
        for name in REQUIRED_RESOURCES
    }
    cases = tuple(
        BenchmarkCase.from_mapping(item, index=index)
        for index, item in enumerate(_sequence(document.get("cases"), "cases"))
    )
    ids = [case.case_id for case in cases]
    if len(ids) != len(set(ids)):
        raise ValueError("benchmark case_id values must be unique")
    dataset = BenchmarkDataset(
        path=resolved_path,
        benchmark_id=_text(document.get("benchmark_id"), "benchmark_id"),
        description=str(document.get("description") or "").strip(),
        resources=resources,
        defaults=BenchmarkDefaults.from_mapping(document.get("defaults")),
        cases=cases,
        dataset_sha256=sha256_file(resolved_path),
    )
    if require_pilot_shape:
        dataset.validate_pilot_shape()
    if verify_resources:
        dataset.verify_resources()
        dataset.validate_a0_snapshot()
    return dataset


def finite_float(value: Any, default: float | None = None) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def split_ids(value: Any) -> tuple[str, ...]:
    if value is None:
        return tuple()
    if isinstance(value, (list, tuple, set)):
        values: Iterable[Any] = value
    else:
        values = re.split(r"[;,|]", str(value))
    return tuple(dict.fromkeys(str(item).strip() for item in values if str(item).strip()))


__all__ = [
    "BENCHMARK_DATASET_SCHEMA",
    "BENCHMARK_RUN_SCHEMA",
    "BENCHMARK_TASK_SCHEMA",
    "BenchmarkCase",
    "BenchmarkDataset",
    "BenchmarkDefaults",
    "BenchmarkGoldStep",
    "BenchmarkResource",
    "canonical_sha256",
    "finite_float",
    "load_benchmark_dataset",
    "repository_root",
    "sha256_file",
    "split_ids",
]
