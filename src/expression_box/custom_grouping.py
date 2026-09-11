"""Parse and validate user-defined expression-box groupings."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence

from src.expression_box.models import (
    ExpressionCassette,
    ExpressionGroupingDesign,
    ExpressionProtein,
)


def parse_custom_expression_groups(
    tokens: Sequence[str],
) -> tuple[tuple[str, ...], ...]:
    """Parse ``[P1 P2] [P3]`` tokens into ordered cassette groups."""

    if isinstance(tokens, (str, bytes)) or not tokens:
        raise ValueError("--custom 后必须提供至少一个 [蛋白 accession ...] 分组")

    raw = " ".join(str(token) for token in tokens).strip()
    if not raw:
        raise ValueError("--custom 后必须提供至少一个 [蛋白 accession ...] 分组")

    groups: list[tuple[str, ...]] = []
    cursor = 0
    while cursor < len(raw):
        while cursor < len(raw) and raw[cursor].isspace():
            cursor += 1
        if cursor >= len(raw):
            break
        if raw[cursor] != "[":
            raise ValueError("--custom 分组格式错误：每个表达盒必须使用 [...] 包围")

        close = raw.find("]", cursor + 1)
        if close < 0:
            raise ValueError("--custom 分组格式错误：缺少右方括号 ]")
        content = raw[cursor + 1 : close]
        if "[" in content:
            raise ValueError("--custom 分组格式错误：不允许嵌套方括号")

        members = tuple(value.upper() for value in content.split())
        if not members:
            raise ValueError("--custom 不允许空表达盒 []")
        groups.append(members)
        cursor = close + 1

    return tuple(groups)


def build_custom_grouping_design(
    tokens: Sequence[str],
    proteins: tuple[ExpressionProtein, ...],
) -> ExpressionGroupingDesign:
    """Build a complete, ordered custom grouping from the current CDS selection."""

    groups = parse_custom_expression_groups(tokens)
    proteins_by_accession = {protein.accession: protein for protein in proteins}
    assigned = [accession for group in groups for accession in group]
    counts = Counter(assigned)

    unknown = sorted(set(assigned) - set(proteins_by_accession))
    duplicates = sorted(accession for accession, count in counts.items() if count > 1)
    missing = sorted(set(proteins_by_accession) - set(assigned))

    problems: list[str] = []
    if unknown:
        problems.append("未知蛋白：" + ", ".join(unknown))
    if duplicates:
        problems.append("重复蛋白：" + ", ".join(duplicates))
    if missing:
        problems.append("遗漏蛋白：" + ", ".join(missing))
    if problems:
        raise ValueError("自定义表达盒分组无效：" + "；".join(problems))

    cassettes = tuple(
        ExpressionCassette(
            proteins=tuple(proteins_by_accession[accession] for accession in group),
            reason="用户自定义表达盒分组",
        )
        for group in groups
    )
    return ExpressionGroupingDesign(
        strategy="custom",
        name="用户自定义方案",
        recommended=False,
        cassettes=cassettes,
    )


__all__ = [
    "build_custom_grouping_design",
    "parse_custom_expression_groups",
]
