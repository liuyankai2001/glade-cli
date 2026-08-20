"""Deterministic system recommendations for protein-to-cassette grouping."""

from __future__ import annotations

from collections import defaultdict

from src.expression_box.config import (
    BALANCED_MAX_MAIN_UNITS_PER_CASSETTE,
    COMPACT_MAX_CDS_LENGTH_NT,
    COMPACT_MAX_PROTEINS,
)
from src.expression_box.models import (
    ExpressionCassette,
    ExpressionGroupingDesign,
    ExpressionProtein,
    GroupingGenerationResult,
    ProteinUnit,
    SkippedGroupingStrategy,
)


def _protein_order(protein: ExpressionProtein) -> tuple[int, str]:
    return protein.first_step_index, protein.accession


def _build_units(
    proteins: tuple[ExpressionProtein, ...],
) -> tuple[tuple[ProteinUnit, ...], tuple[ExpressionProtein, ...], tuple[str, ...]]:
    mains = {
        protein.accession: protein for protein in proteins if protein.is_main_enzyme
    }
    dedicated: dict[str, list[ExpressionProtein]] = defaultdict(list)
    standalone: list[ExpressionProtein] = []
    warnings: list[str] = []

    for protein in proteins:
        if protein.is_main_enzyme:
            if protein.is_auxiliary_protein:
                warnings.append(
                    f"{protein.accession} 同时承担主酶和辅助蛋白角色；"
                    "按主酶分组且不重复放置"
                )
            continue
        if len(protein.required_by_main_accessions) == 1:
            dedicated[protein.required_by_main_accessions[0]].append(protein)
        else:
            standalone.append(protein)
            if not protein.required_by_main_accessions:
                warnings.append(
                    f"辅助蛋白 {protein.accession} 缺少明确主酶依赖；"
                    "系统将其安排为独立表达盒"
                )

    units: list[ProteinUnit] = []
    for main in sorted(mains.values(), key=_protein_order):
        auxiliaries = tuple(
            sorted(dedicated.get(main.accession, []), key=_protein_order)
        )
        units.append(
            ProteinUnit(
                main_accession=main.accession,
                proteins=(main, *auxiliaries),
                first_step_index=main.first_step_index,
            )
        )
    standalone.sort(key=_protein_order)
    return tuple(units), tuple(standalone), tuple(dict.fromkeys(warnings))


def _independent_design(
    units: tuple[ProteinUnit, ...],
    standalone: tuple[ExpressionProtein, ...],
    warnings: tuple[str, ...],
) -> ExpressionGroupingDesign:
    cassettes = [
        ExpressionCassette(
            proteins=unit.proteins,
            reason=("主酶独立调控；专属辅助蛋白与对应主酶保持在同一表达盒"),
        )
        for unit in units
    ]
    cassettes.extend(
        ExpressionCassette(
            proteins=(protein,),
            reason=("共享或无明确归属的辅助蛋白独立表达，避免重复放置"),
        )
        for protein in standalone
    )
    return ExpressionGroupingDesign(
        strategy="independent_control",
        name="独立调控方案",
        recommended=False,
        cassettes=tuple(cassettes),
        warnings=warnings,
    )


def _balanced_design(
    units: tuple[ProteinUnit, ...],
    standalone: tuple[ExpressionProtein, ...],
    warnings: tuple[str, ...],
) -> ExpressionGroupingDesign:
    cassettes: list[ExpressionCassette] = []
    for start in range(0, len(units), BALANCED_MAX_MAIN_UNITS_PER_CASSETTE):
        group = units[start : start + BALANCED_MAX_MAIN_UNITS_PER_CASSETTE]
        members = tuple(protein for unit in group for protein in unit.proteins)
        reason = (
            "合并相邻的两个主酶单元，以平衡独立调控能力和表达盒数量"
            if len(group) > 1
            else "该主酶单元保持独立；专属辅助蛋白随主酶表达"
        )
        cassettes.append(ExpressionCassette(proteins=members, reason=reason))
    cassettes.extend(
        ExpressionCassette(
            proteins=(protein,),
            reason=("共享或无明确归属的辅助蛋白独立表达，避免重复放置"),
        )
        for protein in standalone
    )
    return ExpressionGroupingDesign(
        strategy="balanced_dependency",
        name="均衡方案",
        recommended=True,
        cassettes=tuple(cassettes),
        warnings=warnings,
    )


def _compact_design(
    units: tuple[ProteinUnit, ...],
    standalone: tuple[ExpressionProtein, ...],
    warnings: tuple[str, ...],
) -> ExpressionGroupingDesign:
    members = tuple(
        [protein for unit in units for protein in unit.proteins] + list(standalone)
    )
    return ExpressionGroupingDesign(
        strategy="compact_operon",
        name="紧凑表达方案",
        recommended=False,
        cassettes=(
            ExpressionCassette(
                proteins=members,
                reason="将所有蛋白放入同一表达盒，以减少调控元件数量",
            ),
        ),
        warnings=warnings,
    )


def _validate_design(
    design: ExpressionGroupingDesign,
    proteins: tuple[ExpressionProtein, ...],
) -> None:
    expected = {protein.accession for protein in proteins}
    assigned = [
        protein.accession
        for cassette in design.cassettes
        for protein in cassette.proteins
    ]
    assigned_set = set(assigned)
    if assigned_set != expected:
        difference = sorted(assigned_set ^ expected)
        raise ValueError(
            f"grouping strategy {design.strategy} has missing or extra proteins: "
            + ", ".join(difference)
        )
    if len(assigned) != len(assigned_set):
        duplicates = sorted(
            {accession for accession in assigned if assigned.count(accession) > 1}
        )
        raise ValueError(
            f"grouping strategy {design.strategy} duplicates proteins: "
            + ", ".join(duplicates)
        )
    if not design.cassettes or any(
        not cassette.proteins for cassette in design.cassettes
    ):
        raise ValueError(
            f"grouping strategy {design.strategy} contains an empty design"
        )


def generate_grouping_designs(
    proteins: tuple[ExpressionProtein, ...],
) -> GroupingGenerationResult:
    """Generate balanced, independent, and eligible compact recommendations."""

    if not proteins:
        raise ValueError("expression grouping requires at least one protein")
    units, standalone, warnings = _build_units(proteins)
    if not units:
        raise ValueError("expression grouping requires at least one main-enzyme unit")

    candidates = [
        _balanced_design(units, standalone, warnings),
        _independent_design(units, standalone, warnings),
    ]
    skipped: list[SkippedGroupingStrategy] = []
    total_length = sum(item.optimized_cds_length_nt for item in proteins)
    if (
        len(proteins) <= COMPACT_MAX_PROTEINS
        and total_length <= COMPACT_MAX_CDS_LENGTH_NT
    ):
        candidates.append(_compact_design(units, standalone, warnings))
    else:
        reasons: list[str] = []
        if len(proteins) > COMPACT_MAX_PROTEINS:
            reasons.append(f"蛋白数量 {len(proteins)} 超过上限 {COMPACT_MAX_PROTEINS}")
        if total_length > COMPACT_MAX_CDS_LENGTH_NT:
            reasons.append(
                f"CDS 总长度 {total_length} nt 超过上限 {COMPACT_MAX_CDS_LENGTH_NT} nt"
            )
        skipped.append(
            SkippedGroupingStrategy(
                strategy="compact_operon",
                reason="；".join(reasons),
            )
        )

    designs: list[ExpressionGroupingDesign] = []
    signature_owner: dict[tuple[tuple[str, ...], ...], str] = {}
    for candidate in candidates:
        _validate_design(candidate, proteins)
        owner = signature_owner.get(candidate.signature)
        if owner is not None:
            skipped.append(
                SkippedGroupingStrategy(
                    strategy=candidate.strategy,
                    reason=f"与 {owner} 产生相同分组，已自动去重",
                )
            )
            continue
        signature_owner[candidate.signature] = candidate.strategy
        designs.append(candidate)

    if not designs:
        raise ValueError("no valid expression-box grouping design was generated")
    return GroupingGenerationResult(
        designs=tuple(designs),
        skipped_strategies=tuple(skipped),
        warnings=warnings,
    )


__all__ = ["generate_grouping_designs"]
