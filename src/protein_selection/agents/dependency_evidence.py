"""Atomic protein-dependency evidence and deterministic synthesis gates."""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any, Literal, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from src.protein_selection.reaction_scope import (
    ReactionScopedModel,
    normalize_reaction_ids,
    require_reaction_subset,
)


DependencyNecessity = Literal["required", "enhancing", "uncertain"]
SynthesisNecessity = Literal["required", "enhancing"]
ExperimentType = Literal[
    "subunit_omission",
    "gene_knockout",
    "protein_depletion",
    "gene_loss_of_function",
    "addback_or_complementation",
    "purified_component_comparison",
    "residue_mutagenesis",
    "complex_purification",
    "physical_interaction",
    "other",
]
ActivityStatus = Literal[
    "none",
    "conditional_loss",
    "detectable_reduced",
    "detectable",
    "comparable",
    "increased",
    "full_restored",
    "not_reported",
]
CandidateScope = Literal[
    "whole_protein",
    "residue_or_domain",
    "unspecified",
]
EvidenceFact = Literal[
    "input_identity",
    "candidate_identity",
    "experimental_organism",
    "reaction_context",
    "candidate_absence",
    "candidate_presence",
    "candidate_addition",
    "input_alone",
    "activity_without_candidate",
    "activity_with_candidate",
    "rescue",
    "candidate_component_function",
    "input_component_function",
    "coupled_target_activity",
    "complex_composition",
    "physical_interaction",
    "genetic_loss_of_function",
]
EvidenceSpanFact = Literal[
    "input_identity",
    "candidate_identity",
    "experimental_organism",
    "reaction_context",
    "candidate_absence",
    "candidate_presence",
    "candidate_addition",
    "input_alone",
    "activity_without_candidate",
    "activity_with_candidate",
    "rescue",
    "candidate_component_function",
    "input_component_function",
    "coupled_target_activity",
    "complex_composition",
    "physical_interaction",
    "genetic_loss_of_function",
    "curated_dependency_assertion",
]
DependencyDecisionPath = Literal[
    "loss_and_reconstitution",
    "residual_activity_and_enhancement",
    "reaction_component_coupling",
    "genetic_loss_of_function",
]
AssertionType = Literal[
    "explicit_dependency",
    "complex_composition",
    "gpr",
    "physical_interaction",
    "functional_association",
    "other",
]


class EvidenceSpan(BaseModel):
    """One exact source span supporting named relation facts."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    span_id: str = Field(min_length=1)
    raw_evidence_id: str = Field(min_length=1)
    source_locator: str = Field(min_length=1)
    text: str = Field(min_length=1)
    supports: list[EvidenceSpanFact] = Field(min_length=1)


class DependencyEvidenceAtom(ReactionScopedModel):
    """One source-grounded observation about a protein dependency."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    atom_id: str = Field(min_length=1)
    paper_id: str = Field(min_length=1)
    evidence_id: str = Field(min_length=1)
    input_uniprot_id: str = Field(min_length=1)
    candidate_uniprot_id: str = Field(min_length=1)
    input_protein_mentions: list[str] = Field(min_length=1)
    candidate_protein_mentions: list[str] = Field(min_length=1)
    experimental_organism: str | None = None
    experimental_taxon_id: int | None = None
    rhea_family: str | None = None
    experiment_type: ExperimentType
    candidate_scope: CandidateScope
    fact: EvidenceFact
    activity_status: ActivityStatus = "not_reported"
    evidence_span: EvidenceSpan
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_span_support(self) -> Self:
        if self.fact not in self.evidence_span.supports:
            raise ValueError("atom fact must be named by evidence_span.supports")
        return self


class DependencyEvidenceSynthesis(ReactionScopedModel):
    """A deterministic decision proposal over one or more evidence atoms."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    synthesis_id: str = Field(min_length=1)
    input_uniprot_id: str = Field(min_length=1)
    candidate_uniprot_id: str = Field(min_length=1)
    organism_name: str | None = None
    taxon_id: int | None = None
    rhea_family: str | None = None
    necessity: SynthesisNecessity
    decision_path: DependencyDecisionPath
    atom_ids: list[str] = Field(min_length=1)
    paper_ids: list[str] = Field(min_length=1)
    evidence_ids: list[str] = Field(min_length=1)
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_local_references(self) -> Self:
        for field_name in ("atom_ids", "paper_ids", "evidence_ids"):
            values = getattr(self, field_name)
            if len(values) != len(set(values)):
                raise ValueError(
                    f"dependency synthesis {field_name} values must be unique"
                )
        return self


class CuratedDependencyAssertion(ReactionScopedModel):
    """A protein-level dependency assertion from one curated lineage."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    assertion_id: str = Field(min_length=1)
    evidence_ids: list[str] = Field(min_length=1)
    database: str = Field(min_length=1)
    record_id: str = Field(min_length=1)
    lineage_key: str = Field(min_length=1)
    input_protein_mentions: list[str] = Field(min_length=1)
    candidate_protein_mentions: list[str] = Field(min_length=1)
    candidate_uniprot_id: str | None = None
    organism_name: str | None = None
    taxon_id: int | None = None
    assertion_type: AssertionType
    whole_protein_scope: bool
    necessity: DependencyNecessity
    evidence_spans: list[EvidenceSpan] = Field(min_length=1)
    limitations: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_local_references(self) -> Self:
        span_ids = [span.span_id for span in self.evidence_spans]
        if len(span_ids) != len(set(span_ids)):
            raise ValueError("curated assertion span_id values must be unique")
        return self


class DependencyEvidenceExtractionResult(BaseModel):
    """Small per-paper atom extraction; synthesis remains deterministic."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    paper_id: str = Field(min_length=1)
    raw_evidence_id: str = Field(min_length=1)
    input_uniprot_id: str = Field(min_length=1)
    candidate_uniprot_id: str = Field(min_length=1)
    atoms: list[DependencyEvidenceAtom] = Field(default_factory=list)
    unresolved_reasons: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_scope(self) -> Self:
        atom_ids = [atom.atom_id for atom in self.atoms]
        if len(atom_ids) != len(set(atom_ids)):
            raise ValueError("dependency atom IDs must be unique")
        for atom in self.atoms:
            if atom.paper_id != self.paper_id:
                raise ValueError("atom must reference the requested paper")
            if atom.input_uniprot_id != self.input_uniprot_id:
                raise ValueError("atom input protein must match the request")
            if atom.candidate_uniprot_id != self.candidate_uniprot_id:
                raise ValueError("atom candidate must match the request")
            if atom.evidence_span.raw_evidence_id != self.raw_evidence_id:
                raise ValueError("atom span must reference the requested raw record")
        return self


class CuratedDependencyAssertionExtractionResult(BaseModel):
    """Small database-model output containing only explicit assertions."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    assertions: list[CuratedDependencyAssertion] = Field(default_factory=list)
    unresolved_reasons: list[str] = Field(default_factory=list)


@dataclass(frozen=True, slots=True)
class DependencyEvaluation:
    """Independent validation result consumed by final evidence gates."""

    valid: bool
    necessity: DependencyNecessity
    failures: tuple[str, ...]


def evaluate_dependency_synthesis(
    synthesis: DependencyEvidenceSynthesis,
    atoms: list[DependencyEvidenceAtom],
    *,
    expected_input_uniprot_id: str,
    expected_candidate_uniprot_id: str,
    expected_organism: str | None,
    expected_taxon_id: int | None,
    expected_reaction_ids: Sequence[str],
    expected_rhea_family: str | None = None,
) -> DependencyEvaluation:
    """Validate atom references, scope and one of four decision paths."""

    failures: list[str] = []
    if synthesis.input_uniprot_id != expected_input_uniprot_id:
        failures.append("综合结论的输入蛋白与请求不一致")
    if synthesis.candidate_uniprot_id != expected_candidate_uniprot_id:
        failures.append("综合结论的候选蛋白与请求不一致")
    if not _same_organism(
        synthesis.organism_name,
        synthesis.taxon_id,
        expected_organism,
        expected_taxon_id,
    ):
        failures.append("综合结论的物种与候选来源物种不一致或未核验")
    try:
        require_reaction_subset(
            synthesis.reaction_ids,
            expected_reaction_ids,
            label="dependency synthesis",
        )
    except ValueError:
        failures.append("综合结论未对应请求范围内的目标反应")
    if (
        expected_rhea_family is not None
        and synthesis.rhea_family != expected_rhea_family
    ):
        failures.append("综合结论的 Rhea 反应家族不一致或未核验")

    atom_ids = [atom.atom_id for atom in atoms]
    duplicate_atom_ids = sorted(
        atom_id for atom_id in set(atom_ids) if atom_ids.count(atom_id) > 1
    )
    if duplicate_atom_ids:
        failures.append(
            "输入证据原子的 ID 重复: " + ", ".join(duplicate_atom_ids)
        )
    atoms_by_id = {atom.atom_id: atom for atom in atoms}
    unknown_atom_ids = sorted(set(synthesis.atom_ids) - set(atoms_by_id))
    if unknown_atom_ids:
        failures.append(
            "综合结论引用未知证据原子: " + ", ".join(unknown_atom_ids)
        )
    selected_atoms = [
        atoms_by_id[atom_id]
        for atom_id in synthesis.atom_ids
        if atom_id in atoms_by_id
    ]

    referenced_paper_ids = {atom.paper_id for atom in selected_atoms}
    if set(synthesis.paper_ids) != referenced_paper_ids:
        failures.append("综合结论的 paper_ids 与引用原子集合不一致")
    referenced_evidence_ids = {atom.evidence_id for atom in selected_atoms}
    if set(synthesis.evidence_ids) != referenced_evidence_ids:
        failures.append("综合结论的 evidence_ids 与引用原子集合不一致")

    for atom in selected_atoms:
        if atom.input_uniprot_id != expected_input_uniprot_id:
            failures.append(f"原子 {atom.atom_id} 的输入蛋白不一致")
        if atom.candidate_uniprot_id != expected_candidate_uniprot_id:
            failures.append(f"原子 {atom.atom_id} 的候选蛋白不一致")
        if atom.input_uniprot_id != synthesis.input_uniprot_id:
            failures.append(f"原子 {atom.atom_id} 与综合输入蛋白不一致")
        if atom.candidate_uniprot_id != synthesis.candidate_uniprot_id:
            failures.append(f"原子 {atom.atom_id} 与综合候选蛋白不一致")
        if not _same_organism(
            atom.experimental_organism,
            atom.experimental_taxon_id,
            expected_organism,
            expected_taxon_id,
        ):
            failures.append(f"原子 {atom.atom_id} 的实验物种不一致或未核验")
        try:
            require_reaction_subset(
                atom.reaction_ids,
                expected_reaction_ids,
                label=f"atom {atom.atom_id}",
            )
        except ValueError:
            failures.append(f"原子 {atom.atom_id} 未对应请求范围内的目标反应")
        if normalize_reaction_ids(atom.reaction_ids) != normalize_reaction_ids(
            synthesis.reaction_ids
        ):
            failures.append(f"原子 {atom.atom_id} 与综合结论的反应范围不一致")
        if (
            expected_rhea_family is not None
            and atom.rhea_family is not None
            and atom.rhea_family != expected_rhea_family
        ):
            failures.append(f"原子 {atom.atom_id} 的 Rhea 反应家族不一致")

    facts = {atom.fact for atom in selected_atoms}
    if "input_identity" not in facts or "candidate_identity" not in facts:
        failures.append("缺少输入蛋白或候选蛋白的逐字身份证据原子")
    if "reaction_context" not in facts:
        failures.append("缺少目标反应的逐字证据原子")

    if synthesis.decision_path == "loss_and_reconstitution":
        _evaluate_loss_and_reconstitution(synthesis, selected_atoms, failures)
    elif synthesis.decision_path == "residual_activity_and_enhancement":
        _evaluate_residual_activity(synthesis, selected_atoms, failures)
    elif synthesis.decision_path == "reaction_component_coupling":
        _evaluate_component_coupling(synthesis, selected_atoms, failures)
    else:
        _evaluate_genetic_loss(synthesis, selected_atoms, failures)

    return DependencyEvaluation(
        valid=not failures,
        necessity=synthesis.necessity,
        failures=tuple(failures),
    )


def _evaluate_loss_and_reconstitution(
    synthesis: DependencyEvidenceSynthesis,
    atoms: list[DependencyEvidenceAtom],
    failures: list[str],
) -> None:
    if synthesis.necessity != "required":
        failures.append("缺失/重构路径只能证明 required")
    if not _has_whole_protein_fact(atoms, "candidate_absence"):
        failures.append("required 缺少完整候选蛋白缺失证据")
    if not _has_activity(atoms, "activity_without_candidate", {"none"}):
        failures.append("候选缺失时并非无可检测目标活性")
    if not _has_whole_protein_fact(atoms, "candidate_addition"):
        failures.append("required 缺少完整候选蛋白补回证据")
    if not _has_activity(
        atoms,
        "activity_with_candidate",
        {"detectable", "comparable", "increased", "full_restored"},
    ):
        failures.append("补回候选蛋白后的目标活性结果不足")
    if not _has_fact(atoms, "rescue"):
        failures.append("required 缺少补回后的活性恢复证据")


def _evaluate_residual_activity(
    synthesis: DependencyEvidenceSynthesis,
    atoms: list[DependencyEvidenceAtom],
    failures: list[str],
) -> None:
    if synthesis.necessity != "enhancing":
        failures.append("残余活性/增强路径只能证明 enhancing")
    if not _has_fact(atoms, "input_alone"):
        failures.append("未明确测定输入蛋白单独状态")
    if not _has_activity(
        atoms,
        "activity_without_candidate",
        {"detectable_reduced", "detectable"},
    ):
        failures.append("输入蛋白单独活性未得到明确测量")
    if not _has_whole_protein_fact(atoms, "candidate_addition"):
        failures.append("未明确加入完整候选蛋白")
    if not _has_activity(
        atoms,
        "activity_with_candidate",
        {"increased", "full_restored"},
    ):
        failures.append("加入候选后未显示明确目标活性提高")


def _evaluate_component_coupling(
    synthesis: DependencyEvidenceSynthesis,
    atoms: list[DependencyEvidenceAtom],
    failures: list[str],
) -> None:
    if synthesis.necessity != "required":
        failures.append("反应组分耦联路径只能证明 required")
    if not _has_whole_protein_fact(atoms, "candidate_component_function"):
        failures.append("缺少完整候选蛋白的反应组分功能证据")
    if not _has_fact(atoms, "coupled_target_activity"):
        failures.append("缺少与候选组分耦联的目标活性证据")


def _evaluate_genetic_loss(
    synthesis: DependencyEvidenceSynthesis,
    atoms: list[DependencyEvidenceAtom],
    failures: list[str],
) -> None:
    if synthesis.necessity != "required":
        failures.append("遗传功能缺失路径只能证明 required")
    genetic_atoms = [
        atom
        for atom in atoms
        if atom.fact == "genetic_loss_of_function"
        and atom.candidate_scope == "whole_protein"
        and atom.experiment_type
        in {"gene_knockout", "protein_depletion", "gene_loss_of_function"}
    ]
    if not genetic_atoms:
        failures.append("缺少完整候选蛋白的遗传功能缺失证据")
    if not _has_activity(
        atoms,
        "activity_without_candidate",
        {"none", "conditional_loss"},
    ):
        failures.append("遗传缺失未导致目标活性消失或条件性丧失")


def _has_fact(
    atoms: list[DependencyEvidenceAtom],
    fact: EvidenceFact,
) -> bool:
    return any(atom.fact == fact for atom in atoms)


def _has_whole_protein_fact(
    atoms: list[DependencyEvidenceAtom],
    fact: EvidenceFact,
) -> bool:
    return any(
        atom.fact == fact and atom.candidate_scope == "whole_protein"
        for atom in atoms
    )


def _has_activity(
    atoms: list[DependencyEvidenceAtom],
    fact: EvidenceFact,
    accepted_statuses: set[ActivityStatus],
) -> bool:
    return any(
        atom.fact == fact and atom.activity_status in accepted_statuses
        for atom in atoms
    )


def evaluate_curated_assertion(
    assertion: CuratedDependencyAssertion,
    necessity: Literal["required", "enhancing"],
    *,
    expected_organism: str | None,
    expected_taxon_id: int | None,
    reaction_ids: Sequence[str],
) -> DependencyEvaluation:
    """Accept only an explicit, same-organism, whole-protein assertion."""

    failures: list[str] = []
    if assertion.necessity != necessity:
        failures.append("整理断言的依赖等级与候选等级不一致")
    if assertion.assertion_type != "explicit_dependency":
        failures.append("复合物组成、GPR 或互作不能替代依赖断言")
    if not assertion.whole_protein_scope:
        failures.append("整理断言不是完整蛋白层级")
    try:
        require_reaction_subset(
            assertion.reaction_ids,
            reaction_ids,
            label="curated dependency assertion",
        )
    except ValueError:
        failures.append("整理断言未对应请求范围内的目标反应")
    if not _same_organism(
        assertion.organism_name,
        assertion.taxon_id,
        expected_organism,
        expected_taxon_id,
    ):
        failures.append("整理断言的物种不一致或未核验")
    facts = {
        fact
        for span in assertion.evidence_spans
        for fact in span.supports
    }
    if "curated_dependency_assertion" not in facts:
        failures.append("缺少逐字蛋白依赖断言")
    return DependencyEvaluation(
        valid=not failures,
        necessity=necessity,
        failures=tuple(failures),
    )


def validate_relation_span_grounding(
    spans: list[EvidenceSpan],
    successful_content: dict[str, str],
) -> list[str]:
    """Require every relation span to occur in its named raw record."""

    failures: list[str] = []
    for span in spans:
        raw = successful_content.get(span.raw_evidence_id)
        if raw is None:
            failures.append(
                f"span {span.span_id} cites unavailable raw evidence"
            )
            continue
        if not source_text_occurs(span.text, raw):
            failures.append(f"span {span.span_id} was absent from raw evidence")
    return failures


def source_text_occurs(value: str, raw_content: str) -> bool:
    """Match source text against raw JSON and its decoded string values."""

    needle = _normalize_text(value)
    if not needle:
        return False
    if needle in _normalize_text(raw_content):
        return True
    try:
        payload = json.loads(raw_content)
    except (TypeError, ValueError):
        return False
    return any(
        needle in _normalize_text(candidate)
        for candidate in _iter_source_strings(payload)
    )


def _iter_source_strings(value: Any):
    if isinstance(value, str):
        yield value
    elif isinstance(value, Mapping):
        for key, child in value.items():
            yield str(key)
            yield from _iter_source_strings(child)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for child in value:
            yield from _iter_source_strings(child)


def dependency_claim(
    candidate_name: str,
    necessity: Literal["required", "enhancing"],
) -> str:
    if necessity == "required":
        return (
            f"{candidate_name} is required: validated atomic evidence "
            "establishes complete-partner dependence for target activity."
        )
    return (
        f"{candidate_name} enhances activity: the input remains active alone "
        "and addition of the complete partner increases activity."
    )


def _same_organism(
    observed_name: str | None,
    observed_taxon: int | None,
    expected_name: str | None,
    expected_taxon: int | None,
) -> bool:
    if observed_taxon is not None and expected_taxon is not None:
        if observed_taxon == expected_taxon:
            return True
    if observed_name and expected_name:
        return _normalize_species(observed_name) == _normalize_species(
            expected_name
        )
    return False


def _normalize_species(value: str) -> str:
    without_strain = re.sub(
        r"\s*\([^)]*(?:strain|substrain)[^)]*\)\s*",
        " ",
        value,
        flags=re.IGNORECASE,
    )
    return " ".join(without_strain.casefold().split())


def _normalize_text(value: str) -> str:
    return " ".join(value.split()).casefold()
