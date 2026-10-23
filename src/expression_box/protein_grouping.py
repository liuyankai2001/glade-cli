"""Deterministic system recommendations for protein-to-cassette grouping."""

from __future__ import annotations


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
) -> tuple[ProteinUnit, ...]:
    if any(protein.roles != ("main_enzyme",) for protein in proteins):
        raise ValueError("expression grouping accepts selected main enzymes only")
    return tuple(
        ProteinUnit(
            main_accession=protein.accession,
            proteins=(protein,),
            first_step_index=protein.first_step_index,
        )
        for protein in sorted(proteins, key=_protein_order)
    )


def _independent_design(
    units: tuple[ProteinUnit, ...],
) -> ExpressionGroupingDesign:
    cassettes = [
        ExpressionCassette(
            proteins=unit.proteins,
            reason="每个主酶独立调控",
        )
        for unit in units
    ]
    return ExpressionGroupingDesign(
        strategy="independent_control",
        name="独立调控方案",
        recommended=False,
        cassettes=tuple(cassettes),
    )


def _balanced_design(
    units: tuple[ProteinUnit, ...],
) -> ExpressionGroupingDesign:
    cassettes: list[ExpressionCassette] = []
    for start in range(0, len(units), BALANCED_MAX_MAIN_UNITS_PER_CASSETTE):
        group = units[start : start + BALANCED_MAX_MAIN_UNITS_PER_CASSETTE]
        members = tuple(protein for unit in group for protein in unit.proteins)
        reason = (
            "合并相邻的两个主酶单元，以平衡独立调控能力和表达盒数量"
            if len(group) > 1
            else "该主酶保持独立调控"
        )
        cassettes.append(ExpressionCassette(proteins=members, reason=reason))
    return ExpressionGroupingDesign(
        strategy="balanced_main_enzymes",
        name="均衡方案",
        recommended=True,
        cassettes=tuple(cassettes),
    )


def _compact_design(
    units: tuple[ProteinUnit, ...],
) -> ExpressionGroupingDesign:
    members = tuple(protein for unit in units for protein in unit.proteins)
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
    units = _build_units(proteins)
    if not units:
        raise ValueError("expression grouping requires at least one main-enzyme unit")

    candidates = [
        _balanced_design(units),
        _independent_design(units),
    ]
    skipped: list[SkippedGroupingStrategy] = []
    total_length = sum(item.optimized_cds_length_nt for item in proteins)
    if (
        len(proteins) <= COMPACT_MAX_PROTEINS
        and total_length <= COMPACT_MAX_CDS_LENGTH_NT
    ):
        candidates.append(_compact_design(units))
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
        warnings=(),
    )


__all__ = ["generate_grouping_designs"]
