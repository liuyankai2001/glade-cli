"""Local genome annotation evidence for skipping external enzyme construction."""

from __future__ import annotations

import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .host_registry import get_host_preset
from .genome_identity import GenomeIdentity, read_genome_identity


UNSPECIFIED_GENOME = object()
_COMPLETE_EC = re.compile(r"[1-7]\.\d+\.\d+\.\d+")
_ROOT = Path(__file__).resolve().parents[2]


def resolve_genome_path(value: Any, root_dir: Path = _ROOT) -> Path | None:
    """Resolve an optional input value without opening or validating its file."""
    if value is None:
        return None
    if not isinstance(value, (str, Path)):
        raise ValueError("genome_file must be a path string or null")
    if not str(value).strip():
        return None
    path = Path(str(value).strip()).expanduser()
    return (path if path.is_absolute() else root_dir / path).resolve()


class GenomeGeneEvidence(BaseModel):
    model_config = ConfigDict(extra="forbid")
    record_id: str
    locus_tag: str = ""
    gene: str = ""
    protein_id: str = ""
    product: str = ""


class GenomeEcMatch(BaseModel):
    model_config = ConfigDict(extra="forbid")
    ec_number: str = Field(pattern=r"^[1-7]\.\d+\.\d+\.\d+$")
    source_fields: list[str]
    genes: list[GenomeGeneEvidence] = Field(min_length=1)


class GenomeSupportedStep(BaseModel):
    model_config = ConfigDict(extra="forbid")
    step_index: int = Field(ge=1)
    reaction_id: str
    step_source: str = ""
    rule_ids: list[str] = Field(default_factory=list)
    matches: list[GenomeEcMatch] = Field(min_length=1)


class GenomeEcFilter(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal["genome_ec_filter.v2"] = "genome_ec_filter.v2"
    policy_version: Literal["complete_ec_any_source.v1"] = "complete_ec_any_source.v1"
    genome_file: str
    chassis_key: str
    record_ids: list[str]
    required_step_indexes: list[int]
    external_required_step_indexes: list[int]
    matched_steps: list[GenomeSupportedStep]

    @model_validator(mode="after")
    def validate_partition(self) -> "GenomeEcFilter":
        supported = [item.step_index for item in self.matched_steps]
        for values in (self.required_step_indexes, self.external_required_step_indexes, supported):
            if values != sorted(set(values)) or any(value < 1 for value in values):
                raise ValueError("genome EC step indexes must be positive, sorted and unique")
        if set(supported) & set(self.external_required_step_indexes) or sorted(
            supported + self.external_required_step_indexes
        ) != self.required_step_indexes:
            raise ValueError("genome EC matches must partition the construction requirements")
        return self

    @property
    def supported_step_indexes(self) -> list[int]:
        return [item.step_index for item in self.matched_steps]

    @property
    def no_external_enzymes_required(self) -> bool:
        return bool(self.matched_steps) and not self.external_required_step_indexes


def _values(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return sorted({part for item in value for part in _values(item)})
    return sorted({part.strip() for part in re.split(r"[;|,]", str(value or "")) if part.strip()})


def _genome_index(path: Path, chassis_key: str, identity: GenomeIdentity) -> tuple[list[str], dict[str, list[GenomeGeneEvidence]]]:
    preset = get_host_preset(chassis_key)
    accepted_taxa = {preset["taxon_id"]}
    if chassis_key == "ecoli":
        accepted_taxa.add(511145)
    if identity.taxon_id not in accepted_taxa:
        raise ValueError(f"genome_file 菌株与 {chassis_key} 不匹配: taxon={identity.taxon_id}")
    index: dict[str, list[GenomeGeneEvidence]] = {}
    for record in identity.records:
        pseudo_loci = {locus for feature in record.features if feature.type == "gene"
                       and ("pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers)
                       for locus in feature.qualifiers.get("locus_tag", [])}
        for feature in record.features:
            qualifiers = feature.qualifiers
            if (feature.type != "CDS" or "pseudo" in qualifiers or "pseudogene" in qualifiers
                    or pseudo_loci.intersection(qualifiers.get("locus_tag", []))):
                continue
            gene = GenomeGeneEvidence(record_id=record.id, **{
                name: str((qualifiers.get(name) or [""])[0])
                for name in ("locus_tag", "gene", "protein_id", "product")
            })
            for ec in sorted(set(qualifiers.get("EC_number", []))):
                if _COMPLETE_EC.fullmatch(ec):
                    index.setdefault(ec, []).append(gene)
    return sorted(record.id for record in identity.records), index


def build_genome_ec_filter(
    genome_path: str | Path, requirements: Sequence[Mapping[str, Any]], chassis_key: str,
    *, genome_identity: GenomeIdentity | None = None,
) -> GenomeEcFilter:
    path = resolve_genome_path(genome_path)
    if path is None:
        raise ValueError("genome_file must not be empty when filtering")
    identity = genome_identity if genome_identity is not None else read_genome_identity(path)
    if identity.path != path:
        raise ValueError("genome_file identity belongs to a different path")
    record_ids, index = _genome_index(path, chassis_key, identity)
    matched: list[GenomeSupportedStep] = []
    input_requirements: list[dict[str, Any]] = []
    required = sorted({int(item["step_index"]) for item in requirements})
    for item in sorted(requirements, key=lambda value: int(value["step_index"])):
        fields = {
            "enzyme_ecs": item.get("enzyme_ecs", item.get("ec_numbers")),
            "locked_enzyme_ecs": item.get("locked_enzyme_ecs", item.get("locked_ec_numbers")),
        }
        provenance = item.get("prediction_provenance") or {}
        if not isinstance(provenance, Mapping):
            provenance = {}
        if item.get("step_source") == "retropath":
            fields.update({"source_ec_numbers": item.get("source_ec_numbers"),
                           "prediction_provenance.source_ec_numbers": provenance.get("source_ec_numbers")})
        rule_ids = sorted(set(_values(item.get("retropath_rule_id")) + _values(provenance.get("rule_ids"))
                              + _values(provenance.get("rule_id"))))
        # Unmatched alternatives still affect the input binding even when the
        # same genomic EC continues to support this step.
        input_requirements.append({
            "step_index": int(item["step_index"]), "reaction_id": str(item.get("reaction_id") or ""),
            "step_source": str(item.get("step_source") or ""), "rule_ids": rule_ids,
            "ec_sources": {name: _values(value) for name, value in fields.items()},
        })
        sources: dict[str, list[str]] = {}
        for name, value in fields.items():
            for ec in _values(value):
                if _COMPLETE_EC.fullmatch(ec) and ec in index:
                    sources.setdefault(ec, []).append(name)
        if sources:
            matched.append(GenomeSupportedStep(
                step_index=int(item["step_index"]), reaction_id=str(item.get("reaction_id") or ""),
                step_source=str(item.get("step_source") or ""),
                rule_ids=rule_ids,
                matches=[GenomeEcMatch(ec_number=ec, source_fields=sorted(names), genes=index[ec])
                         for ec, names in sorted(sources.items())],
            ))
    supported = {item.step_index for item in matched}
    payload = {"schema_version": "genome_ec_filter.v2", "policy_version": "complete_ec_any_source.v1",
               "genome_file": str(path), "chassis_key": chassis_key,
               "record_ids": record_ids, "required_step_indexes": required,
               "external_required_step_indexes": sorted(set(required) - supported),
               "matched_steps": [item.model_dump(mode="json") for item in matched]}
    return GenomeEcFilter(**payload)


def validate_genome_ec_filter(
    record: GenomeEcFilter | Mapping[str, Any] | None,
    requirements: Sequence[Mapping[str, Any]], chassis_key: str,
    genome_path: Any = UNSPECIFIED_GENOME,
) -> GenomeEcFilter | None:
    """Recompute live evidence so changed configuration or files cannot skip work."""
    expected = resolve_genome_path(genome_path) if genome_path is not UNSPECIFIED_GENOME else UNSPECIFIED_GENOME
    if record is None:
        if expected is not UNSPECIFIED_GENOME and expected is not None:
            raise ValueError("genome_file configuration changed; rerun main-enzyme")
        return None
    saved = record if isinstance(record, GenomeEcFilter) else GenomeEcFilter.model_validate(record)
    if expected is not UNSPECIFIED_GENOME and expected != Path(saved.genome_file).resolve():
        raise ValueError("genome_file configuration changed; rerun main-enzyme")
    try:
        live = build_genome_ec_filter(saved.genome_file, requirements, chassis_key)
    except (OSError, ValueError) as exc:
        raise ValueError(f"genome_file is stale; rerun main-enzyme: {exc}") from exc
    if live.external_required_step_indexes != saved.external_required_step_indexes:
        raise ValueError("genome EC coverage changed; rerun main-enzyme")
    return live


def read_current_genome_selection(manifest_path: str | Path, genome_path: Any = UNSPECIFIED_GENOME):
    """Read the current selection only when genome filtering is relevant."""
    from .common import get_solution_steps, heterologous_requirements, read_manifest
    from .models import MainEnzymeSelectionResult

    manifest_path = Path(manifest_path)
    path = manifest_path.parent / "main_protein_selection" / "main_enzyme_selection.json"
    configured = genome_path is not UNSPECIFIED_GENOME and resolve_genome_path(genome_path) is not None
    if not path.is_file():
        if configured:
            raise ValueError("genome_file requires current candidate results; rerun main-enzyme")
        return None
    raw = json.loads(path.read_text(encoding="utf-8-sig"))
    if raw.get("genome_ec_filter") is None and not configured:
        return None
    selection = MainEnzymeSelectionResult.model_validate(raw)
    manifest = read_manifest(manifest_path)
    solution_id, steps = get_solution_steps(manifest)
    if (selection.selected_solution_id != solution_id
            or selection.expansion_depth != int(manifest["solution"].get("expansion_depth") or 0)):
        raise ValueError("Selected route changed; rerun main-enzyme")
    validate_genome_ec_filter(selection.genome_ec_filter, heterologous_requirements(steps),
                             selection.chassis_key, genome_path)
    return selection
