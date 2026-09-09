"""Read-only pFBA tracing from medium substrates to gap-route anchors."""

from __future__ import annotations

import csv
import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import cobra
from cobra.flux_analysis import pfba
from cobra.util.solver import linear_reaction_coefficients

from src.pathway_analyze.gem_validation import split_equation


DEFAULT_GROWTH_FRACTION = 0.1
DEFAULT_FLUX_THRESHOLD = 1e-8
DEFAULT_NODE_LIMIT = 500
KEGG_ANNOTATION_KEYS = ("kegg.compound", "kegg.drug", "kegg.glycan")

# Carbon-containing molecules that are inorganic inputs rather than organic
# pathway roots.  Both KEGG and common BiGG identifiers are checked.
INORGANIC_CARBON_KEGG_IDS = frozenset({
    "C00011",  # CO2
    "C00237",  # CO
    "C00288",  # bicarbonate
    "C01353",  # carbonate
})
INORGANIC_CARBON_BIGG_IDS = frozenset({"co", "co2", "h2co3", "hco3"})

# These metabolites support energy, redox or group transfer.  They are shown
# as side dependencies and are deliberately not used as main-carbon shortcuts.
CURRENCY_METABOLITE_BIGG_IDS = frozenset({
    "adp",
    "amp",
    "atp",
    "cdp",
    "cmp",
    "coa",
    "ctp",
    "fad",
    "fadh2",
    "flxr",
    "flxso",
    "gdp",
    "gmp",
    "gtp",
    "h",
    "h2o",
    "mql8",
    "mqn8",
    "nad",
    "nadh",
    "nadp",
    "nadph",
    "pi",
    "ppi",
    "q8",
    "q8h2",
    "udp",
    "ump",
    "utp",
})


def _split_values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [
            item.strip()
            for item in re.split(r"[;,|]", value)
            if item.strip()
        ]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    return [text] if text else []


def _base_metabolite_id(metabolite_id: str) -> str:
    value = str(metabolite_id or "").strip().lower()
    return re.sub(r"_(?:c|e|p|m|n|r|x|v|g|u|l)$", "", value)


def _annotation_values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [str(item) for item in value if item is not None]
    return [str(value)]


def _kegg_ids(metabolite: cobra.Metabolite) -> list[str]:
    annotation = getattr(metabolite, "annotation", {}) or {}
    values: list[str] = []
    for key in KEGG_ANNOTATION_KEYS:
        values.extend(_annotation_values(annotation.get(key)))
    return list(dict.fromkeys(
        str(value).strip().upper().split(":")[-1]
        for value in values
        if str(value).strip()
    ))


def _reaction_kegg_ids(reaction: cobra.Reaction) -> list[str]:
    annotation = getattr(reaction, "annotation", {}) or {}
    return list(dict.fromkeys(
        str(value).strip().upper().split(":")[-1]
        for value in _annotation_values(annotation.get("kegg.reaction"))
        if str(value).strip()
    ))


def _elements(metabolite: cobra.Metabolite) -> dict[str, float] | None:
    if not str(getattr(metabolite, "formula", "") or "").strip():
        return None
    try:
        elements = metabolite.elements or {}
    except (TypeError, ValueError):
        return None
    return {str(key): float(value) for key, value in elements.items()}


def _has_carbon(metabolite: cobra.Metabolite) -> bool:
    elements = _elements(metabolite)
    if elements is not None:
        return float(elements.get("C", 0.0)) > 0
    # Missing formulae must remain traceable.  They are classified as unknown
    # at a boundary instead of silently terminating an otherwise valid route.
    return True


def _is_inorganic_carbon(metabolite: cobra.Metabolite) -> bool:
    return bool(
        set(_kegg_ids(metabolite)) & INORGANIC_CARBON_KEGG_IDS
        or _base_metabolite_id(metabolite.id) in INORGANIC_CARBON_BIGG_IDS
    )


def _is_currency_metabolite(metabolite: cobra.Metabolite) -> bool:
    return _base_metabolite_id(metabolite.id) in CURRENCY_METABOLITE_BIGG_IDS


def _is_organic(metabolite: cobra.Metabolite) -> bool:
    elements = _elements(metabolite)
    return bool(
        elements is not None
        and float(elements.get("C", 0.0)) > 0
        and not _is_inorganic_carbon(metabolite)
    )


def _metabolite_payload(
    metabolite: cobra.Metabolite,
    *,
    exchange_reaction_id: str = "",
    uptake_flux: float | None = None,
) -> dict[str, Any]:
    return {
        "model_metabolite_id": metabolite.id,
        "name": metabolite.name or metabolite.id,
        "compartment": metabolite.compartment,
        "formula": metabolite.formula or "",
        "kegg_ids": _kegg_ids(metabolite),
        "exchange_reaction_id": exchange_reaction_id,
        "uptake_flux": uptake_flux,
    }


def _oriented_sides(
    reaction: cobra.Reaction,
    flux: float,
) -> tuple[list[tuple[cobra.Metabolite, float]], list[tuple[cobra.Metabolite, float]]]:
    if flux >= 0:
        substrates = [
            (metabolite, abs(float(coefficient)))
            for metabolite, coefficient in reaction.metabolites.items()
            if coefficient < 0
        ]
        products = [
            (metabolite, abs(float(coefficient)))
            for metabolite, coefficient in reaction.metabolites.items()
            if coefficient > 0
        ]
    else:
        substrates = [
            (metabolite, abs(float(coefficient)))
            for metabolite, coefficient in reaction.metabolites.items()
            if coefficient > 0
        ]
        products = [
            (metabolite, abs(float(coefficient)))
            for metabolite, coefficient in reaction.metabolites.items()
            if coefficient < 0
        ]
    return substrates, products


def _format_term(metabolite: cobra.Metabolite, coefficient: float) -> str:
    prefix = "" if math.isclose(coefficient, 1.0, abs_tol=1e-12) else f"{coefficient:g} "
    return f"{prefix}{metabolite.id}"


def _oriented_equation(
    substrates: Sequence[tuple[cobra.Metabolite, float]],
    products: Sequence[tuple[cobra.Metabolite, float]],
) -> str:
    left = " + ".join(_format_term(*item) for item in substrates) or "∅"
    right = " + ".join(_format_term(*item) for item in products) or "∅"
    return f"{left} -> {right}"


def _canonical_compound_id(
    value: Any,
    compound_aliases: Mapping[str, str],
) -> str:
    compound_id = str(value or "").strip()
    alias = str(compound_aliases.get(compound_id) or "").strip()
    match = re.match(r"^([CDG]\d{5})(?:\b|\s|\()", alias, flags=re.IGNORECASE)
    return match.group(1).upper() if match else compound_id.upper()


def _parse_generic_equation_side(
    side: str,
    compound_aliases: Mapping[str, str],
) -> tuple[tuple[str, float], ...]:
    rows: list[tuple[str, float]] = []
    for raw_part in str(side or "").split("+"):
        tokens = raw_part.strip().split()
        if not tokens:
            continue
        amount = 1.0
        if len(tokens) > 1 and re.fullmatch(r"\d+(?:\.\d+)?", tokens[0]):
            amount = float(tokens[0])
        compound_id = _canonical_compound_id(tokens[-1], compound_aliases)
        if compound_id:
            rows.append((compound_id, amount))
    return tuple(rows)


def _directional_row_stoichiometry(
    row: Mapping[str, Any],
    compound_aliases: Mapping[str, str],
) -> tuple[tuple[tuple[str, float], ...], tuple[tuple[str, float], ...]]:
    equation = str(row.get("equation") or "").strip()
    left, right = split_equation(equation)
    left_values = _parse_generic_equation_side(left, compound_aliases)
    right_values = _parse_generic_equation_side(right, compound_aliases)
    direction = str(row.get("direction") or "").strip().lower()
    if direction in {"reverse", "right_to_left", "backward", "rtl"}:
        return right_values, left_values
    return left_values, right_values


def derive_anchor_requirements(
    steps: Sequence[Mapping[str, Any]],
    *,
    target_compound: str,
    anchor_compound_ids: Sequence[str],
    compound_aliases: Mapping[str, str] | None = None,
) -> dict[str, float]:
    """Derive anchor demand ratios required for one unit of route target."""

    aliases = dict(compound_aliases or {})
    target = _canonical_compound_id(target_compound, aliases)
    anchors = tuple(dict.fromkeys(
        _canonical_compound_id(value, aliases)
        for value in anchor_compound_ids
        if str(value or "").strip()
    ))
    if not anchors:
        raise ValueError("路线没有可用于上游追踪的底盘锚点")

    rows_by_product: dict[str, Mapping[str, Any]] = {}
    for row in steps:
        product = _canonical_compound_id(row.get("produced_compound_id"), aliases)
        if not product:
            continue
        if product in rows_by_product:
            raise ValueError(f"路线包含多个生成 {product} 的步骤，无法确定计量关系")
        rows_by_product[product] = row

    requirements = {anchor: 0.0 for anchor in anchors}
    active: set[str] = set()

    def require(compound_id: str, amount: float) -> None:
        compound = _canonical_compound_id(compound_id, aliases)
        if compound in requirements:
            requirements[compound] += amount
            return
        if compound in active:
            raise ValueError(f"gap 路线在 {compound} 处存在计量循环")
        row = rows_by_product.get(compound)
        if row is None:
            raise ValueError(
                f"无法从目标追溯到锚点：{compound} 缺少对应路线步骤"
            )

        reaction_id = str(row.get("reaction_id") or "").strip()
        consumed, produced = _directional_row_stoichiometry(row, aliases)
        produced_amount = next(
            (value for item, value in produced if item.upper() == compound),
            None,
        )
        if produced_amount is None or produced_amount <= 0:
            raise ValueError(
                f"路线步骤 {reaction_id or '<unknown>'} 的方程无法确认 "
                f"{compound} 的生成计量"
            )
        extent = amount / produced_amount
        precursor_ids = _split_values(row.get("precursor_compound_ids"))
        if not precursor_ids:
            raise ValueError(
                f"路线步骤 {reaction_id or '<unknown>'} 缺少主要前体"
            )
        consumed_by_id = {item.upper(): value for item, value in consumed}
        active.add(compound)
        try:
            for precursor in precursor_ids:
                precursor_id = _canonical_compound_id(precursor, aliases)
                precursor_amount = consumed_by_id.get(precursor_id)
                if precursor_amount is None or precursor_amount <= 0:
                    raise ValueError(
                        f"路线步骤 {reaction_id or '<unknown>'} 的方程无法确认 "
                        f"前体 {precursor_id} 的消耗计量"
                    )
                require(precursor_id, extent * precursor_amount)
        finally:
            active.remove(compound)

    require(target, 1.0)
    missing = [anchor for anchor, amount in requirements.items() if amount <= 0]
    if missing:
        raise ValueError(f"路线锚点未连接到目标：{', '.join(missing)}")
    return requirements


def _read_medium(path: Path) -> dict[str, float]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"无法读取培养基配置：{path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"培养基配置根节点必须是对象：{path}")
    try:
        return {str(key): float(value) for key, value in payload.items()}
    except (TypeError, ValueError) as exc:
        raise ValueError(f"培养基摄取上限必须是数值：{path}") from exc


def _read_chassis_settings(path: Path) -> tuple[float, float]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            values = {
                str(row.get("item") or "").strip(): str(row.get("value") or "").strip()
                for row in csv.DictReader(handle)
            }
    except OSError as exc:
        raise ValueError(f"无法读取底盘分析摘要，请先运行 chassis：{path}") from exc
    try:
        growth_fraction = float(
            values.get("growth_fraction", DEFAULT_GROWTH_FRACTION)
        )
        flux_threshold = float(
            values.get("flux_threshold", DEFAULT_FLUX_THRESHOLD)
        )
    except ValueError as exc:
        raise ValueError(f"底盘分析摘要中的通量参数无效：{path}") from exc
    if not 0 < growth_fraction <= 1 or flux_threshold <= 0:
        raise ValueError(f"底盘分析摘要中的通量参数超出有效范围：{path}")
    return growth_fraction, flux_threshold


def _build_model_kegg_index(
    model: cobra.Model,
) -> dict[str, list[cobra.Metabolite]]:
    result: dict[str, list[cobra.Metabolite]] = {}
    for metabolite in model.metabolites:
        for kegg_id in _kegg_ids(metabolite):
            result.setdefault(kegg_id, []).append(metabolite)
    return result


def _map_anchor_metabolites(
    model: cobra.Model,
    anchor_requirements: Mapping[str, float],
    producible_csv: Path,
) -> dict[str, cobra.Metabolite]:
    try:
        with producible_csv.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError as exc:
        raise ValueError(
            f"无法读取底盘可生成代谢物，请先运行 chassis：{producible_csv}"
        ) from exc

    model_index = _build_model_kegg_index(model)
    mapped: dict[str, cobra.Metabolite] = {}
    for anchor in anchor_requirements:
        csv_ids = [
            str(row.get("met_id") or "").strip()
            for row in rows
            if str(row.get("kegg_id") or "").strip().upper() == anchor
            and str(row.get("met_id") or "").strip() in model.metabolites
        ]
        candidates = [model.metabolites.get_by_id(value) for value in csv_ids]
        candidates.extend(
            metabolite
            for metabolite in model_index.get(anchor, [])
            if metabolite not in candidates
        )
        candidates.sort(key=lambda item: (
            0 if item.compartment == "c" else 1,
            item.compartment,
            item.id,
        ))
        if not candidates:
            raise ValueError(f"底盘 GEM 中无法映射路线锚点 {anchor}")
        mapped[anchor] = candidates[0]
    return mapped


def _objective_reaction_id(model: cobra.Model) -> str:
    coefficients = linear_reaction_coefficients(model)
    if not coefficients:
        raise ValueError("底盘 GEM 没有可识别的生长目标反应")
    reaction, _ = max(coefficients.items(), key=lambda item: abs(item[1]))
    return reaction.id


def _active_uptakes(
    model: cobra.Model,
    fluxes: Any,
    flux_threshold: float,
) -> dict[str, tuple[cobra.Metabolite, float]]:
    result: dict[str, tuple[cobra.Metabolite, float]] = {}
    for reaction in model.exchanges:
        if reaction.id not in fluxes.index:
            continue
        flux = float(fluxes.loc[reaction.id])
        if abs(flux) <= flux_threshold:
            continue
        _, products = _oriented_sides(reaction, flux)
        for metabolite, coefficient in products:
            uptake = abs(flux) * coefficient
            if uptake > flux_threshold:
                result[reaction.id] = (metabolite, uptake)
    return result


def _trace_active_network(
    model: cobra.Model,
    fluxes: Any,
    anchor_metabolites: Mapping[str, cobra.Metabolite],
    anchor_requirements: Mapping[str, float],
    anchor_sink_flux: float,
    *,
    anchor_sink_id: str,
    flux_threshold: float,
    node_limit: int,
) -> dict[str, Any]:
    exchange_ids = {reaction.id for reaction in model.exchanges}
    production_index: dict[
        str,
        list[tuple[cobra.Reaction, float, float]],
    ] = {}
    oriented_by_reaction: dict[
        str,
        tuple[
            list[tuple[cobra.Metabolite, float]],
            list[tuple[cobra.Metabolite, float]],
            float,
        ],
    ] = {}
    for reaction in model.reactions:
        if reaction.id == anchor_sink_id or reaction.id not in fluxes.index:
            continue
        flux = float(fluxes.loc[reaction.id])
        if abs(flux) <= flux_threshold:
            continue
        substrates, products = _oriented_sides(reaction, flux)
        oriented_by_reaction[reaction.id] = (substrates, products, flux)
        for metabolite, coefficient in products:
            production_index.setdefault(metabolite.id, []).append(
                (reaction, coefficient, abs(flux) * coefficient)
            )
    for producers in production_index.values():
        producers.sort(key=lambda item: (-item[2], item[0].id))

    selected_extent: dict[str, float] = {}
    selected_distance: dict[str, int] = {}
    allocated_production: dict[tuple[str, str], float] = {}
    organic_sources: dict[str, dict[str, Any]] = {}
    side_dependencies: dict[str, dict[str, Any]] = {}
    unresolved: set[str] = set()
    warnings: set[str] = set()
    active_stack: set[str] = set()
    truncated = False

    def remember_side(metabolite: cobra.Metabolite, reaction_id: str) -> None:
        payload = side_dependencies.setdefault(
            metabolite.id,
            {
                **_metabolite_payload(metabolite),
                "used_by_reaction_ids": [],
            },
        )
        if reaction_id not in payload["used_by_reaction_ids"]:
            payload["used_by_reaction_ids"].append(reaction_id)

    def visit(metabolite: cobra.Metabolite, required_rate: float, distance: int) -> None:
        nonlocal truncated
        if required_rate <= flux_threshold or truncated:
            return
        if metabolite.id in active_stack:
            warnings.add(f"检测到循环并在 {metabolite.id} 截断")
            return
        if _is_currency_metabolite(metabolite):
            remember_side(metabolite, "")
            return

        active_stack.add(metabolite.id)
        remaining = required_rate
        try:
            producers = production_index.get(metabolite.id, [])
            if not producers:
                unresolved.add(metabolite.id)
                return
            for reaction, product_coefficient, capacity in producers:
                key = (metabolite.id, reaction.id)
                available = capacity - allocated_production.get(key, 0.0)
                if available <= flux_threshold:
                    continue
                used = min(remaining, available)
                allocated_production[key] = allocated_production.get(key, 0.0) + used
                required_extent = allocated_production[key] / product_coefficient
                previous_extent = selected_extent.get(reaction.id, 0.0)
                if reaction.id not in selected_extent and len(selected_extent) >= node_limit:
                    truncated = True
                    warnings.add(f"内源路线超过 {node_limit} 个反应节点，结果已截断")
                    return
                selected_extent[reaction.id] = max(previous_extent, required_extent)
                selected_distance[reaction.id] = max(
                    selected_distance.get(reaction.id, 0),
                    distance,
                )
                if reaction.id in exchange_ids:
                    if _is_organic(metabolite):
                        organic_sources[metabolite.id] = _metabolite_payload(
                            metabolite,
                            exchange_reaction_id=reaction.id,
                            uptake_flux=capacity,
                        )
                    else:
                        if _elements(metabolite) is None:
                            warnings.add(
                                f"边界代谢物 {metabolite.id} 缺少可靠分子式，"
                                "未判定为有机起点"
                            )
                        remember_side(metabolite, reaction.id)
                elif reaction.boundary:
                    unresolved.add(metabolite.id)
                    warnings.add(
                        f"{metabolite.id} 由非 exchange 边界反应 {reaction.id} 提供"
                    )
                elif required_extent > previous_extent + flux_threshold:
                    substrates, _, _ = oriented_by_reaction[reaction.id]
                    delta_extent = required_extent - previous_extent
                    for substrate, coefficient in substrates:
                        if _is_currency_metabolite(substrate) or not _has_carbon(substrate):
                            remember_side(substrate, reaction.id)
                        elif _is_inorganic_carbon(substrate):
                            remember_side(substrate, reaction.id)
                        else:
                            visit(
                                substrate,
                                delta_extent * coefficient,
                                distance + 1,
                            )
                remaining -= used
                if remaining <= flux_threshold:
                    break
            if remaining > flux_threshold:
                unresolved.add(metabolite.id)
        finally:
            active_stack.remove(metabolite.id)

    for anchor, metabolite in anchor_metabolites.items():
        visit(
            metabolite,
            anchor_requirements[anchor] * anchor_sink_flux,
            0,
        )

    reaction_rows: list[dict[str, Any]] = []
    for reaction_id in selected_extent:
        reaction = model.reactions.get_by_id(reaction_id)
        substrates, products, flux = oriented_by_reaction[reaction_id]
        reaction_rows.append({
            "model_reaction_id": reaction.id,
            "reaction_name": reaction.name or reaction.id,
            "direction": "forward" if flux >= 0 else "reverse",
            "signed_pfba_flux": flux,
            "pfba_flux": abs(flux),
            "attributed_flux": selected_extent[reaction_id],
            "equation": _oriented_equation(substrates, products),
            "substrates": [
                {**_metabolite_payload(item), "coefficient": coefficient}
                for item, coefficient in substrates
            ],
            "products": [
                {**_metabolite_payload(item), "coefficient": coefficient}
                for item, coefficient in products
            ],
            "kegg_reaction_ids": _reaction_kegg_ids(reaction),
            "gene_ids": sorted(gene.id for gene in reaction.genes),
            "is_exchange": reaction.id in exchange_ids,
            "distance_to_anchor": selected_distance[reaction_id],
        })
    reaction_rows.sort(key=lambda row: (
        0 if row["is_exchange"] else 1,
        -int(row["distance_to_anchor"]),
        str(row["model_reaction_id"]),
    ))

    active_uptakes = _active_uptakes(model, fluxes, flux_threshold)
    auxiliary_inputs = [
        _metabolite_payload(
            metabolite,
            exchange_reaction_id=reaction_id,
            uptake_flux=uptake,
        )
        for reaction_id, (metabolite, uptake) in active_uptakes.items()
        if not _is_organic(metabolite)
    ]
    other_organic_inputs = [
        _metabolite_payload(
            metabolite,
            exchange_reaction_id=reaction_id,
            uptake_flux=uptake,
        )
        for reaction_id, (metabolite, uptake) in active_uptakes.items()
        if _is_organic(metabolite) and metabolite.id not in organic_sources
    ]
    if other_organic_inputs:
        warnings.add(
            "存在只用于生长或未分配到锚点路径的其他有机培养基摄取："
            + ", ".join(item["model_metabolite_id"] for item in other_organic_inputs)
        )

    return {
        "trace_complete": not truncated and not unresolved,
        "truncated": truncated,
        "organic_sources": sorted(
            organic_sources.values(),
            key=lambda item: (item["name"], item["model_metabolite_id"]),
        ),
        "auxiliary_medium_inputs": sorted(
            auxiliary_inputs,
            key=lambda item: item["model_metabolite_id"],
        ),
        "other_active_organic_inputs": sorted(
            other_organic_inputs,
            key=lambda item: item["model_metabolite_id"],
        ),
        "side_dependencies": sorted(
            side_dependencies.values(),
            key=lambda item: item["model_metabolite_id"],
        ),
        "reactions": reaction_rows,
        "unresolved_carbon_metabolites": sorted(unresolved),
        "warnings": sorted(warnings),
    }


def trace_upstream_pathway(
    *,
    model_path: str | Path,
    medium_path: str | Path,
    chassis_producible_csv: str | Path,
    chassis_summary_csv: str | Path,
    steps: Sequence[Mapping[str, Any]],
    target_compound: str,
    anchor_compound_ids: Sequence[str],
    compound_aliases: Mapping[str, str] | None = None,
    node_limit: int = DEFAULT_NODE_LIMIT,
) -> dict[str, Any]:
    """Return one parsimonious native pathway from active medium roots to anchors."""

    if node_limit < 1:
        raise ValueError("上游路线节点上限必须大于等于 1")
    model_file = Path(model_path).expanduser().resolve()
    medium_file = Path(medium_path).expanduser().resolve()
    producible_file = Path(chassis_producible_csv).expanduser().resolve()
    summary_file = Path(chassis_summary_csv).expanduser().resolve()
    if not model_file.is_file():
        raise FileNotFoundError(f"未找到底盘 GEM：{model_file}")
    if model_file.suffix.lower() != ".json":
        raise ValueError(f"完整路线当前只支持 JSON GEM：{model_file}")
    if not medium_file.is_file():
        raise FileNotFoundError(f"未找到培养基配置：{medium_file}")

    anchor_requirements = derive_anchor_requirements(
        steps,
        target_compound=target_compound,
        anchor_compound_ids=anchor_compound_ids,
        compound_aliases=compound_aliases,
    )
    growth_fraction, flux_threshold = _read_chassis_settings(summary_file)
    model = cobra.io.load_json_model(str(model_file))
    model.medium = _read_medium(medium_file)
    biomass_reaction_id = _objective_reaction_id(model)
    baseline_growth = model.slim_optimize()
    if baseline_growth is None or not math.isfinite(float(baseline_growth)):
        raise ValueError("当前培养基下底盘 GEM 无法求出生长通量")
    baseline_growth = float(baseline_growth)
    if baseline_growth <= flux_threshold:
        raise ValueError("当前培养基下底盘 GEM 没有有效生长通量")
    required_growth = baseline_growth * growth_fraction
    biomass = model.reactions.get_by_id(biomass_reaction_id)
    biomass.lower_bound = max(float(biomass.lower_bound), required_growth)

    anchor_metabolites = _map_anchor_metabolites(
        model,
        anchor_requirements,
        producible_file,
    )
    pathway_anchor_requirements = {
        anchor: requirement
        for anchor, requirement in anchor_requirements.items()
        if _is_organic(anchor_metabolites[anchor])
        and not _is_currency_metabolite(anchor_metabolites[anchor])
    }
    if not pathway_anchor_requirements:
        raise ValueError("路线没有可用于追踪有机底物来源的主碳骨架锚点")
    sink_id = "DM_INFO_ALL_ROUTE_ANCHORS"
    if sink_id in model.reactions:
        raise ValueError(f"底盘 GEM 已存在保留反应 ID：{sink_id}")
    sink = cobra.Reaction(sink_id)
    sink.name = "Coupled demand for info --solution --all"
    sink.lower_bound = 0.0
    sink.upper_bound = 1000.0
    sink.add_metabolites({
        anchor_metabolites[anchor]: -float(requirement)
        for anchor, requirement in anchor_requirements.items()
    })
    model.add_reactions([sink])
    model.objective = sink
    model.objective_direction = "max"

    fba_solution = model.optimize()
    if fba_solution.status != "optimal":
        raise ValueError(f"锚点联合供给 FBA 求解失败：{fba_solution.status}")
    fba_anchor_flux = float(fba_solution.objective_value or 0.0)
    if fba_anchor_flux <= flux_threshold:
        raise ValueError("当前培养基和生长约束下无法联合供给路线锚点")
    try:
        pfba_solution = pfba(model, fraction_of_optimum=1.0)
    except Exception as exc:
        raise ValueError(f"锚点联合供给 pFBA 求解失败：{exc}") from exc
    if pfba_solution.status != "optimal":
        raise ValueError(f"锚点联合供给 pFBA 求解失败：{pfba_solution.status}")
    anchor_sink_flux = float(pfba_solution.fluxes.loc[sink_id])
    if anchor_sink_flux <= flux_threshold:
        raise ValueError("pFBA 解没有有效的锚点联合供给通量")

    trace = _trace_active_network(
        model,
        pfba_solution.fluxes,
        {
            anchor: anchor_metabolites[anchor]
            for anchor in pathway_anchor_requirements
        },
        pathway_anchor_requirements,
        anchor_sink_flux,
        anchor_sink_id=sink_id,
        flux_threshold=flux_threshold,
        node_limit=node_limit,
    )
    return {
        "schema_version": "upstream_pathway.v1",
        "trace_complete": trace["trace_complete"],
        "model_path": str(model_file),
        "medium_path": str(medium_file),
        "target_compound": str(target_compound).strip().upper(),
        "anchors": [
            {
                "kegg_id": anchor,
                "relative_requirement": requirement,
                "role": (
                    "main_carbon_anchor"
                    if anchor in pathway_anchor_requirements
                    else "route_side_dependency"
                ),
                **_metabolite_payload(anchor_metabolites[anchor]),
            }
            for anchor, requirement in anchor_requirements.items()
        ],
        "organic_sources": trace["organic_sources"],
        "auxiliary_medium_inputs": trace["auxiliary_medium_inputs"],
        "other_active_organic_inputs": trace["other_active_organic_inputs"],
        "side_dependencies": trace["side_dependencies"],
        "reactions": trace["reactions"],
        "unresolved_carbon_metabolites": trace[
            "unresolved_carbon_metabolites"
        ],
        "warnings": trace["warnings"],
        "pfba_summary": {
            "baseline_growth": baseline_growth,
            "growth_fraction": growth_fraction,
            "required_growth": required_growth,
            "biomass_reaction_id": biomass_reaction_id,
            "flux_threshold": flux_threshold,
            "fba_anchor_flux": fba_anchor_flux,
            "pfba_anchor_flux": anchor_sink_flux,
            "active_reaction_count": sum(
                abs(float(value)) > flux_threshold
                for value in pfba_solution.fluxes
            ),
            "traced_reaction_count": len(trace["reactions"]),
            "node_limit": node_limit,
            "truncated": trace["truncated"],
        },
    }


__all__ = [
    "DEFAULT_NODE_LIMIT",
    "derive_anchor_requirements",
    "trace_upstream_pathway",
]
