from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import cobra
from cobra.util.solver import fix_objective_as_constraint


DEFAULT_GROWTH_FRACTION = 0.1
DEFAULT_FLUX_THRESHOLD = 1e-8
DEFAULT_TEST_COMPARTMENTS = ("c",)


def _test_compartments(config: Any) -> set[str] | None:
    value = getattr(config, "test_compartments", DEFAULT_TEST_COMPARTMENTS)
    if value is None:
        return None
    if isinstance(value, str):
        value = value.replace(",", " ").split()
    compartments = {str(item).strip() for item in value if str(item).strip()}
    return compartments or None


def _load_medium(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Medium JSON root must be an object: {path}")

    try:
        return {str(reaction_id): float(bound) for reaction_id, bound in payload.items()}
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Medium bounds must be numeric: {path}") from exc


def _kegg_ids(metabolite: cobra.Metabolite) -> list[str]:
    annotation = getattr(metabolite, "annotation", {}) or {}
    values: list[str] = []
    for key in ("kegg.compound", "kegg.drug", "kegg.glycan"):
        annotation_values = annotation.get(key, [])
        if not isinstance(annotation_values, (list, tuple, set)):
            annotation_values = [annotation_values]
        values.extend(str(value) for value in annotation_values if value is not None)
    return list(dict.fromkeys(values))


def _analyze_producibility(
    model: cobra.Model,
    growth_fraction: float,
    flux_threshold: float,
    compartments: set[str] | None,
) -> tuple[list[dict[str, str]], float, int]:
    """在维持最低生长率时，逐个检测胞内代谢物的最大 demand flux。"""

    baseline_growth = model.slim_optimize()
    if baseline_growth is None or baseline_growth <= flux_threshold:
        raise ValueError("当前培养基下模型几乎不生长，请检查培养基设置")

    external_counts = Counter()
    for reaction in model.exchanges:
        exchange_metabolites = list(reaction.metabolites)
        if len(exchange_metabolites) == 1:
            external_counts[exchange_metabolites[0].compartment] += 1
    external = (
        {external_counts.most_common(1)[0][0]} if external_counts else {"e"}
    )
    metabolites = [
        metabolite
        for metabolite in model.metabolites
        if metabolite.compartment not in external
        and (compartments is None or metabolite.compartment in compartments)
    ]
    kegg_rows: list[dict[str, str]] = []
    producible_count = 0

    for index, metabolite in enumerate(metabolites, start=1):
        with model:
            fix_objective_as_constraint(model, fraction=growth_fraction)
            demand = model.add_boundary(metabolite, type="demand")
            model.objective = demand
            model.objective_direction = "max"
            flux = model.slim_optimize()

        if flux is not None and flux > flux_threshold:
            producible_count += 1
            for kegg_id in _kegg_ids(metabolite):
                kegg_rows.append(
                    {
                        "source": "producible",
                        "met_id": metabolite.id,
                        "met_name": metabolite.name,
                        "compartment": metabolite.compartment,
                        "kegg_id": kegg_id,
                    }
                )

        if index % 500 == 0:
            print(f"[INFO] tested {index}/{len(metabolites)} metabolites")

    kegg_rows.sort(key=lambda row: (row["kegg_id"], row["met_id"]))
    return kegg_rows, float(baseline_growth), producible_count


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def analyze_chassis_metabolites(config: Any) -> dict[str, Any]:
    model_path = Path(config.model_path).expanduser().resolve()
    medium_path = Path(config.medium_path).expanduser().resolve()
    output_dir = Path(config.chassis_output_path).expanduser().resolve()
    producible_csv = Path(config.chassis_producible_csv).expanduser().resolve()
    summary_csv = Path(config.chassis_metabolites_summary_csv).expanduser().resolve()
    growth_fraction = float(
        getattr(config, "growth_fraction", DEFAULT_GROWTH_FRACTION)
    )
    flux_threshold = float(getattr(config, "flux_threshold", DEFAULT_FLUX_THRESHOLD))
    compartments = _test_compartments(config)

    if not model_path.is_file():
        raise FileNotFoundError(f"GEM model file not found: {model_path}")
    if model_path.suffix.lower() != ".json":
        raise ValueError(f"Only JSON GEM models are supported: {model_path}")
    if not medium_path.is_file():
        raise FileNotFoundError(f"Medium file not found: {medium_path}")
    if not 0 < growth_fraction <= 1:
        raise ValueError("growth_fraction must be in the interval (0, 1]")
    if flux_threshold <= 0:
        raise ValueError("flux_threshold must be greater than 0")

    output_dir.mkdir(parents=True, exist_ok=True)
    model = cobra.io.load_json_model(str(model_path))
    model.medium = _load_medium(medium_path)
    kegg_rows, baseline_growth, producible_count = _analyze_producibility(
        model,
        growth_fraction,
        flux_threshold,
        compartments,
    )
    _write_csv(
        producible_csv,
        ["source", "met_id", "met_name", "compartment", "kegg_id"],
        kegg_rows,
    )

    summary_rows = [
        {"item": "baseline_growth", "value": baseline_growth},
        {"item": "growth_fraction", "value": growth_fraction},
        {"item": "flux_threshold", "value": flux_threshold},
        {
            "item": "test_compartments",
            "value": ";".join(sorted(compartments or ())) or "non_external",
        },
        {"item": "producible_metabolites", "value": producible_count},
        {"item": "producible_kegg_compounds", "value": len(kegg_rows)},
        {"item": "reachable_kegg_compounds", "value": len(kegg_rows)},
    ]
    _write_csv(summary_csv, ["item", "value"], summary_rows)

    return {
        "ok": True,
        "baseline_growth": baseline_growth,
        "growth_fraction": growth_fraction,
        "flux_threshold": flux_threshold,
        "test_compartments": sorted(compartments or ()),
        "producible_metabolites": producible_count,
        "producible_kegg_compounds": len(kegg_rows),
        "reachable_kegg_file": str(producible_csv),
        "summary_file": str(summary_csv),
    }


def run_chassis(config: Any) -> dict[str, Any]:
    """运行入口；JSON 读取和 RunConfig 构造由 ``main.py`` 负责。"""

    result = analyze_chassis_metabolites(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result
