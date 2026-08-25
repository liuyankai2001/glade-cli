"""Build deterministic RetroPath source/sink CSV files from KEGG structures."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

from src.pathway_analyze.expand_chassis_metabolites import ExpansionBundle
from src.pathway_analyze.retropath_models import PredictedCompound
from src.pathway_analyze.retropath_structure import (
    StructureProvider,
    StructureResolutionError,
)


TARGET_SOURCE_FILE_NAME = "target_source.csv"
CHASSIS_SINK_FILE_NAME = "chassis_sink.csv"
COMPOUND_MAPPING_FILE_NAME = "compound_mapping.csv"
REJECTED_COMPOUNDS_FILE_NAME = "rejected_compounds.csv"

KEGG_COMPOUND_ID_PATTERN = re.compile(r"^C\d{5}$")
ROLE_ORDER = {"target": 0, "sink": 1}

RETROPATH_COMPOUND_COLUMNS = ("Name", "InChI")
COMPOUND_MAPPING_COLUMNS = (
    "role",
    "kegg_id",
    "representative_kegg_id",
    "is_representative",
    "minimum_depth",
    "inchi",
    "inchikey",
    "isomeric_smiles",
    "formula",
    "charge",
    "structure_provenance",
)
REJECTED_COMPOUND_COLUMNS = (
    "role",
    "kegg_id",
    "minimum_depth",
    "reason_code",
    "reason_detail",
)


@dataclass(frozen=True)
class StructureRejection:
    """One target or sink compound excluded before RetroPath execution."""

    role: str
    kegg_id: str
    minimum_depth: Optional[int]
    reason_code: str
    reason_detail: str

    def to_row(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "kegg_id": self.kegg_id,
            "minimum_depth": (
                "" if self.minimum_depth is None else self.minimum_depth
            ),
            "reason_code": self.reason_code,
            "reason_detail": self.reason_detail,
        }


@dataclass(frozen=True)
class CompoundMapping:
    """One KEGG ID to target/sink representative structure mapping."""

    role: str
    compound: PredictedCompound
    representative_kegg_id: str
    is_representative: bool

    def to_row(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "kegg_id": self.compound.compound_id,
            "representative_kegg_id": self.representative_kegg_id,
            "is_representative": str(self.is_representative).lower(),
            "minimum_depth": (
                ""
                if self.compound.minimum_depth is None
                else self.compound.minimum_depth
            ),
            "inchi": self.compound.inchi,
            "inchikey": self.compound.inchikey or "",
            "isomeric_smiles": self.compound.isomeric_smiles or "",
            "formula": self.compound.formula or "",
            "charge": "" if self.compound.charge is None else self.compound.charge,
            "structure_provenance": json.dumps(
                list(self.compound.structure_provenance),
                ensure_ascii=False,
                separators=(",", ":"),
            ),
        }


@dataclass(frozen=True)
class RetroPathInputBundle:
    """Paths, hashes, structures, and audit records for one P2 build."""

    expansion_depth: int
    reachable_compound_count: int
    target_compound: PredictedCompound
    sink_compounds: Tuple[PredictedCompound, ...]
    mappings: Tuple[CompoundMapping, ...]
    rejected_compounds: Tuple[StructureRejection, ...]
    target_source_path: Path
    chassis_sink_path: Path
    compound_mapping_path: Path
    rejected_compounds_path: Path
    target_source_sha256: str
    chassis_sink_sha256: str

    @property
    def sink_structure_count(self) -> int:
        return len(self.sink_compounds)

    @property
    def rejected_compound_count(self) -> int:
        return len(self.rejected_compounds)


class RetroPathInputBuildError(ValueError):
    """Input generation failed after writing deterministic audit artifacts."""

    def __init__(
        self,
        code: str,
        detail: str,
        *,
        output_dir: Path,
        rejected_compounds: Sequence[StructureRejection],
    ) -> None:
        self.code = str(code).strip()
        self.detail = str(detail).strip()
        self.output_dir = output_dir
        self.rejected_compounds = tuple(rejected_compounds)
        super().__init__(f"{self.code}: {self.detail}; audit directory: {output_dir}")


def _normalize_target_compound_id(value: str) -> str:
    compound_id = str(value).strip().upper()
    if not KEGG_COMPOUND_ID_PATTERN.fullmatch(compound_id):
        raise ValueError("target_compound_id must be a KEGG Cxxxxx identifier")
    return compound_id


def _normalized_bundle_depth(bundle: ExpansionBundle) -> int:
    if not isinstance(bundle, ExpansionBundle):
        raise ValueError("expansion_bundle must be an ExpansionBundle")
    depth = bundle.depth
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 0:
        raise ValueError("expansion_bundle.depth must be greater than or equal to 0")
    return depth


def _complete_structure_or_error(
    compound: PredictedCompound,
    *,
    expected_id: str,
    expected_depth: Optional[int],
) -> PredictedCompound:
    if not isinstance(compound, PredictedCompound):
        raise StructureResolutionError(
            "structure_model_invalid",
            expected_id,
            "structure provider did not return PredictedCompound",
        )
    if compound.compound_id != expected_id or expected_id not in compound.kegg_ids:
        raise StructureResolutionError(
            "structure_identity_mismatch",
            expected_id,
            f"provider returned identity {compound.compound_id!r}",
        )
    if compound.minimum_depth != expected_depth:
        raise StructureResolutionError(
            "minimum_depth_mismatch",
            expected_id,
            (
                f"provider returned depth {compound.minimum_depth!r}; "
                f"expected {expected_depth!r}"
            ),
        )
    missing = []
    if not compound.inchikey:
        missing.append("inchikey")
    if not compound.isomeric_smiles:
        missing.append("isomeric_smiles")
    if not compound.formula:
        missing.append("formula")
    if compound.charge is None:
        missing.append("charge")
    if missing:
        raise StructureResolutionError(
            "structure_fields_incomplete",
            expected_id,
            f"missing required P2 fields: {', '.join(missing)}",
        )
    return compound


def _rejection(
    role: str,
    kegg_id: str,
    minimum_depth: Optional[int],
    error: StructureResolutionError,
) -> StructureRejection:
    return StructureRejection(
        role=role,
        kegg_id=kegg_id,
        minimum_depth=minimum_depth,
        reason_code=error.code,
        reason_detail=error.detail,
    )


def _resolve_compound(
    provider: StructureProvider,
    *,
    role: str,
    compound_id: str,
    minimum_depth: Optional[int],
) -> tuple[Optional[PredictedCompound], Optional[StructureRejection]]:
    try:
        compound = provider.resolve(compound_id, minimum_depth=minimum_depth)
        return (
            _complete_structure_or_error(
                compound,
                expected_id=compound_id,
                expected_depth=minimum_depth,
            ),
            None,
        )
    except StructureResolutionError as exc:
        return None, _rejection(role, compound_id, minimum_depth, exc)


def _validate_sink_depth(
    compound_id: str,
    depth_by_compound: Mapping[str, int],
    maximum_depth: int,
) -> tuple[Optional[int], Optional[StructureRejection]]:
    raw_depth = depth_by_compound.get(compound_id)
    if (
        isinstance(raw_depth, bool)
        or not isinstance(raw_depth, int)
        or raw_depth < 0
    ):
        return None, StructureRejection(
            role="sink",
            kegg_id=compound_id,
            minimum_depth=None,
            reason_code="missing_minimum_depth",
            reason_detail="the cumulative expansion does not contain a valid minimum depth",
        )
    if raw_depth > maximum_depth:
        return None, StructureRejection(
            role="sink",
            kegg_id=compound_id,
            minimum_depth=raw_depth,
            reason_code="minimum_depth_exceeds_requested_depth",
            reason_detail=(
                f"minimum depth {raw_depth} exceeds requested cumulative "
                f"depth {maximum_depth}"
            ),
        )
    return raw_depth, None


def _merge_sink_group(
    compounds: Sequence[PredictedCompound],
) -> PredictedCompound:
    representative = min(
        compounds,
        key=lambda item: (
            item.minimum_depth if item.minimum_depth is not None else 2**31,
            item.compound_id,
        ),
    )
    all_kegg_ids = tuple(sorted(item.compound_id for item in compounds))
    all_provenance = tuple(
        sorted(
            {
                value
                for compound in compounds
                for value in compound.structure_provenance
            }
        )
    )
    return PredictedCompound.create(
        compound_id=representative.compound_id,
        name=representative.name,
        inchi=representative.inchi,
        inchikey=representative.inchikey,
        isomeric_smiles=representative.isomeric_smiles,
        formula=representative.formula,
        charge=representative.charge,
        kegg_ids=all_kegg_ids,
        minimum_depth=representative.minimum_depth,
        structure_provenance=all_provenance,
    )


def _group_sink_structures(
    compounds: Iterable[PredictedCompound],
) -> tuple[
    Tuple[PredictedCompound, ...],
    Tuple[CompoundMapping, ...],
    Tuple[StructureRejection, ...],
]:
    by_inchikey: dict[str, list[PredictedCompound]] = {}
    for compound in compounds:
        if compound.inchikey is None:
            raise ValueError("P2 sink compound unexpectedly has no InChIKey")
        by_inchikey.setdefault(compound.inchikey, []).append(compound)

    representatives: list[PredictedCompound] = []
    mappings: list[CompoundMapping] = []
    rejections: list[StructureRejection] = []
    for inchikey in sorted(by_inchikey):
        group = sorted(
            by_inchikey[inchikey],
            key=lambda item: (
                item.minimum_depth if item.minimum_depth is not None else 2**31,
                item.compound_id,
            ),
        )
        observed_inchis = {item.inchi for item in group}
        if len(observed_inchis) != 1:
            ids = ",".join(item.compound_id for item in group)
            for item in group:
                rejections.append(
                    StructureRejection(
                        role="sink",
                        kegg_id=item.compound_id,
                        minimum_depth=item.minimum_depth,
                        reason_code="structure_identity_conflict",
                        reason_detail=(
                            f"InChIKey {inchikey} maps to multiple standard InChI "
                            f"values across: {ids}"
                        ),
                    )
                )
            continue
        representative = _merge_sink_group(group)
        representatives.append(representative)
        mappings.extend(
            CompoundMapping(
                role="sink",
                compound=item,
                representative_kegg_id=representative.compound_id,
                is_representative=item.compound_id == representative.compound_id,
            )
            for item in group
        )

    representatives.sort(key=lambda item: item.compound_id)
    mappings.sort(
        key=lambda item: (
            ROLE_ORDER[item.role],
            item.representative_kegg_id,
            item.compound.compound_id,
        )
    )
    rejections.sort(key=_rejection_sort_key)
    return tuple(representatives), tuple(mappings), tuple(rejections)


def _rejection_sort_key(item: StructureRejection) -> tuple[Any, ...]:
    return (
        ROLE_ORDER.get(item.role, 99),
        item.minimum_depth if item.minimum_depth is not None else -1,
        item.kegg_id,
        item.reason_code,
        item.reason_detail,
    )


def _render_csv(
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


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


def _write_artifacts(
    output_dir: Path,
    *,
    target_compound: Optional[PredictedCompound],
    sink_compounds: Sequence[PredictedCompound],
    mappings: Sequence[CompoundMapping],
    rejections: Sequence[StructureRejection],
) -> tuple[Path, Path, Path, Path, str, str]:
    target_source_path = output_dir / TARGET_SOURCE_FILE_NAME
    chassis_sink_path = output_dir / CHASSIS_SINK_FILE_NAME
    compound_mapping_path = output_dir / COMPOUND_MAPPING_FILE_NAME
    rejected_compounds_path = output_dir / REJECTED_COMPOUNDS_FILE_NAME

    target_text = _render_csv(
        RETROPATH_COMPOUND_COLUMNS,
        (
            []
            if target_compound is None
            else [
                {
                    "Name": target_compound.compound_id,
                    "InChI": target_compound.inchi,
                }
            ]
        ),
    )
    sink_text = _render_csv(
        RETROPATH_COMPOUND_COLUMNS,
        (
            {"Name": compound.compound_id, "InChI": compound.inchi}
            for compound in sink_compounds
        ),
    )
    mapping_text = _render_csv(
        COMPOUND_MAPPING_COLUMNS,
        (item.to_row() for item in mappings),
    )
    rejected_text = _render_csv(
        REJECTED_COMPOUND_COLUMNS,
        (item.to_row() for item in sorted(rejections, key=_rejection_sort_key)),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    for path, text in (
        (target_source_path, target_text),
        (chassis_sink_path, sink_text),
        (compound_mapping_path, mapping_text),
        (rejected_compounds_path, rejected_text),
    ):
        _atomic_write_text(path, text)
    return (
        target_source_path,
        chassis_sink_path,
        compound_mapping_path,
        rejected_compounds_path,
        hashlib.sha256(target_text.encode("utf-8")).hexdigest(),
        hashlib.sha256(sink_text.encode("utf-8")).hexdigest(),
    )


def build_retropath_inputs(
    target_compound_id: str,
    expansion_bundle: ExpansionBundle,
    structure_provider: StructureProvider,
    output_dir: str | Path,
) -> RetroPathInputBundle:
    """Build source, cumulative sink, mapping, and rejection CSV artifacts."""

    target_id = _normalize_target_compound_id(target_compound_id)
    expansion_depth = _normalized_bundle_depth(expansion_bundle)
    resolved_output_dir = Path(output_dir).expanduser().resolve()
    reachable_ids = tuple(
        sorted(str(value).strip().upper() for value in expansion_bundle.reachable_compounds)
    )

    target, target_rejection = _resolve_compound(
        structure_provider,
        role="target",
        compound_id=target_id,
        minimum_depth=None,
    )
    if target is None:
        rejections = tuple(
            item for item in (target_rejection,) if item is not None
        )
        _write_artifacts(
            resolved_output_dir,
            target_compound=None,
            sink_compounds=tuple(),
            mappings=tuple(),
            rejections=rejections,
        )
        raise RetroPathInputBuildError(
            "target_structure_invalid",
            f"target {target_id} has no usable structure",
            output_dir=resolved_output_dir,
            rejected_compounds=rejections,
        )

    resolved_sink: list[PredictedCompound] = []
    rejections: list[StructureRejection] = []
    for compound_id in reachable_ids:
        if not KEGG_COMPOUND_ID_PATTERN.fullmatch(compound_id):
            rejections.append(
                StructureRejection(
                    role="sink",
                    kegg_id=compound_id,
                    minimum_depth=None,
                    reason_code="invalid_kegg_compound_id",
                    reason_detail="expected a KEGG Compound identifier in Cxxxxx format",
                )
            )
            continue
        minimum_depth, depth_rejection = _validate_sink_depth(
            compound_id,
            expansion_bundle.depth_by_compound,
            expansion_depth,
        )
        if depth_rejection is not None:
            rejections.append(depth_rejection)
            continue
        compound, structure_rejection = _resolve_compound(
            structure_provider,
            role="sink",
            compound_id=compound_id,
            minimum_depth=minimum_depth,
        )
        if structure_rejection is not None:
            rejections.append(structure_rejection)
            continue
        if compound is not None:
            resolved_sink.append(compound)

    sink_compounds, sink_mappings, conflict_rejections = _group_sink_structures(
        resolved_sink
    )
    rejections.extend(conflict_rejections)
    rejections.sort(key=_rejection_sort_key)
    target_mapping = CompoundMapping(
        role="target",
        compound=target,
        representative_kegg_id=target.compound_id,
        is_representative=True,
    )
    mappings = (target_mapping, *sink_mappings)
    artifacts = _write_artifacts(
        resolved_output_dir,
        target_compound=target,
        sink_compounds=sink_compounds,
        mappings=mappings,
        rejections=rejections,
    )
    if not sink_compounds:
        raise RetroPathInputBuildError(
            "sink_structure_empty",
            "no cumulative chassis compound has a usable unique structure",
            output_dir=resolved_output_dir,
            rejected_compounds=rejections,
        )

    return RetroPathInputBundle(
        expansion_depth=expansion_depth,
        reachable_compound_count=len(reachable_ids),
        target_compound=target,
        sink_compounds=sink_compounds,
        mappings=mappings,
        rejected_compounds=tuple(rejections),
        target_source_path=artifacts[0],
        chassis_sink_path=artifacts[1],
        compound_mapping_path=artifacts[2],
        rejected_compounds_path=artifacts[3],
        target_source_sha256=artifacts[4],
        chassis_sink_sha256=artifacts[5],
    )


__all__ = [
    "CHASSIS_SINK_FILE_NAME",
    "COMPOUND_MAPPING_COLUMNS",
    "COMPOUND_MAPPING_FILE_NAME",
    "CompoundMapping",
    "REJECTED_COMPOUND_COLUMNS",
    "REJECTED_COMPOUNDS_FILE_NAME",
    "RETROPATH_COMPOUND_COLUMNS",
    "RetroPathInputBuildError",
    "RetroPathInputBundle",
    "StructureRejection",
    "TARGET_SOURCE_FILE_NAME",
    "build_retropath_inputs",
]
