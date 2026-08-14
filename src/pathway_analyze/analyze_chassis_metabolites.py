from pathlib import Path
import json
import pandas as pd
from typing import Any, Counter, Dict, List
import cobra
from cobra.util.solver import fix_objective_as_constraint
from langchain.tools import tool

from src.config.pathway_config import GEM_VALIDATION_CONFIG
from src.runtime.monitor import monitor
from src.tools.common.session_paths import (
    chassis_reachability_dir,
    default_medium_file,
    gem_model_file,
    outputs_dir as resolve_outputs_dir,
)


CHASSIS_ANALYSIS_RUNNING_MARKER = ".analyze_chassis_metabolites.running"


# 保留旧字段名，便于实验脚本记录；阈值含义统一见 pathway_config.py。
GROWTH_FRACTION = GEM_VALIDATION_CONFIG.growth_fraction
FLUX_THRESHOLD = GEM_VALIDATION_CONFIG.flux_threshold
TEST_COMPARTMENTS = set(GEM_VALIDATION_CONFIG.test_compartments)

def safe_mkdir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_table(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, encoding="utf-8-sig")


def guess_external_compartments(model: cobra.Model) -> set[str]:
    """
    从 exchange reactions 猜测外部区室。
    常见结果是 {'e'}。
    """
    counter = Counter()

    for rxn in model.exchanges:
        mets = list(rxn.metabolites.keys())
        if len(mets) == 1:
            counter[mets[0].compartment] += 1

    if not counter:
        return {"e"}

    most_common_compartment, _ = counter.most_common(1)[0]
    return {most_common_compartment}


def choose_test_metabolites(
    model: cobra.Model,
    external_compartments: set[str],
    compartments: set[str] | None = None,
) -> List[cobra.Metabolite]:
    selected = []

    for met in model.metabolites:
        if met.compartment in external_compartments:
            continue
        if compartments is not None and met.compartment not in compartments:
            continue
        selected.append(met)

    return selected

def normalize_annotation_value(value: Any) -> List[str]:
    """
    将接收到的字符串、列表、元组、集合转换为列表
    """
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple, set)):
        out = []
        for x in value:
            if x is None:
                continue
            out.append(str(x))
        return out
    return [str(value)]
def extract_kegg_ids(met: cobra.Metabolite) -> List[str]:
    """
    从 metabolite.annotation 里提取 KEGG 相关注释。
    常见键可能包括:
    - kegg.compound
    - kegg.drug
    - kegg.glycan

    若模型本身没有这些注释，则返回空列表。
    """
    ann = getattr(met, "annotation", {}) or {}
    kegg_ids = []

    for key in ("kegg.compound", "kegg.drug", "kegg.glycan"):
        if key in ann:
            kegg_ids.extend(normalize_annotation_value(ann[key]))

    # 去重且保持顺序
    seen = set()
    uniq = []
    for x in kegg_ids:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq

def metabolite_basic_record(met: cobra.Metabolite) -> Dict[str, Any]:
    kegg_ids = extract_kegg_ids(met)

    return {
        "met_id": met.id,
        "met_name": met.name,
        "compartment": met.compartment,
        "formula": met.formula,
        "charge": met.charge,
        "kegg_ids": ";".join(kegg_ids),
        "kegg_count": len(kegg_ids),
    }


def get_producible_metabolites(
    model: cobra.Model,
    growth_fraction: float = 0.1,
    flux_threshold: float = 1e-8,
    compartments: set[str] | None = None,
) -> pd.DataFrame:
    """
    在当前培养基下，测试哪些内部代谢物可被网络生产。
    做法：
    1. 先求最大生长率
    2. 固定至少 growth_fraction * max_growth 的生长
    3. 对每个代谢物临时加 demand 反应并最大化它
    4. 若 demand 最优值 > 阈值，则认为该代谢物可生产
    """
    baseline_growth = model.slim_optimize()
    if baseline_growth is None or baseline_growth <= flux_threshold:
        raise ValueError(
            "当前培养基下模型几乎不生长。请先检查培养基设置是否合理。"
        )

    external_compartments = guess_external_compartments(model)
    test_mets = choose_test_metabolites(model, external_compartments, compartments)

    rows = []

    for i, met in enumerate(test_mets, start=1):
        with model:
            # 先把“当前目标函数”（通常是 biomass）固定住一部分
            fix_objective_as_constraint(model, fraction=growth_fraction)

            # 临时加 demand reaction，测试该代谢物最大可生成量
            dm = model.add_boundary(met, type="demand")
            model.objective = dm
            model.objective_direction = "max"

            value = model.slim_optimize()

            if value is not None and value > flux_threshold:
                rec = metabolite_basic_record(met)
                rec.update(
                    {
                        "max_demand_flux": float(value),
                        "baseline_growth": float(baseline_growth),
                        "required_growth_fraction": growth_fraction,
                    }
                )
                rows.append(rec)

        if i % 500 == 0:
            print(f"[INFO] tested {i}/{len(test_mets)} metabolites")

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["compartment", "met_id"]).reset_index(drop=True)
    return df

def dataframe_to_kegg_table(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    """
    把带 kegg_ids 列的表，展开成一行一个 KEGG ID。
    """
    rows = []

    if df.empty or "kegg_ids" not in df.columns:
        return pd.DataFrame(columns=["source", "met_id", "met_name", "kegg_id"])

    for _, row in df.iterrows():
        kegg_ids = str(row["kegg_ids"]).split(";") if pd.notna(row["kegg_ids"]) else []
        kegg_ids = [x.strip() for x in kegg_ids if x.strip()]

        for kid in kegg_ids:
            rows.append(
                {
                    "source": source_name,
                    "met_id": row["met_id"],
                    "met_name": row["met_name"],
                    "compartment": row.get("compartment", ""),
                    "kegg_id": kid,
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.drop_duplicates().sort_values(["kegg_id", "met_id"]).reset_index(drop=True)
    return out

@tool
def analyze_chassis_metabolites() -> dict[str, Any]:
    """
    分析底盘 GEM 在给定培养基下可产生的 KEGG 化合物集合。
    
    调用时机：通路 gap 分析前需要刷新 chassis reachability。
    返回：成功状态、基线生长率、分析阈值、可产物数量及两个输出文件路径。
    限制：只读模型并写分析输出；不写 design_manifest。
    """

    tool_name = "analyze_chassis_metabolites"
    monitor.report_start(tool_name)
    try:
        model_path_obj = gem_model_file()
        medium_path_obj = default_medium_file()
        output_root = resolve_outputs_dir()
        outdir = safe_mkdir(chassis_reachability_dir())
        running_marker = outdir / CHASSIS_ANALYSIS_RUNNING_MARKER
        running_marker.write_text(
            json.dumps(
                {
                    "tool": tool_name,
                    "model_path": str(model_path_obj),
                    "medium_path": str(medium_path_obj),
                    "output_dir": str(output_root),
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        monitor.report_running(tool_name, "正在加载 GEM 模型和培养基...", progress=0.1)
        model = cobra.io.load_json_model(str(model_path_obj))
        with medium_path_obj.open(mode='r', encoding="utf-8") as f:
            medium = json.load(f)
        model.medium = dict(medium)
        growth = model.slim_optimize()
        monitor.report_running(tool_name, "正在计算底盘可生产代谢物集合...", progress=0.35)
        producible_df = get_producible_metabolites(
            model,
            growth_fraction=GROWTH_FRACTION,
            flux_threshold=FLUX_THRESHOLD,
            compartments=TEST_COMPARTMENTS,
        )
        monitor.report_running(tool_name, "正在整理 KEGG 可达化合物表...", progress=0.8)
        producible_kegg_df = dataframe_to_kegg_table(producible_df, "producible")
        reachable_kegg_df = producible_kegg_df.copy()
        reachable_kegg_file = (outdir / "producible_kegg_compounds.csv").resolve()
        summary_file = (outdir / "analyze_chassis_metabolites_summary.csv").resolve()
        save_table(producible_kegg_df, reachable_kegg_file)
        analyze_chassis_metabolites_summary = pd.DataFrame(
            [
                {"item": "baseline_growth", "value": float(growth)},
                {"item": "growth_fraction", "value": float(GROWTH_FRACTION)},
                {"item": "flux_threshold", "value": float(FLUX_THRESHOLD)},
                {"item": "test_compartments", "value": ";".join(sorted(TEST_COMPARTMENTS or [])) or "non_external"},
                {"item": "producible_metabolites", "value": int(len(producible_df))},
                {"item": "producible_kegg_compounds", "value": int(len(producible_kegg_df))},
                {"item": "reachable_kegg_compounds", "value": int(len(reachable_kegg_df))},
            ]
        )
        save_table(analyze_chassis_metabolites_summary, summary_file)

        result = {
            "ok": True,
            "baseline_growth": float(growth),
            "growth_fraction": float(GROWTH_FRACTION),
            "flux_threshold": float(FLUX_THRESHOLD),
            "test_compartments": sorted(TEST_COMPARTMENTS or []),
            "producible_metabolites": int(len(producible_df)),
            "producible_kegg_compounds": int(len(producible_kegg_df)),
            "reachable_kegg_file": str(reachable_kegg_file),
            "summary_file": str(summary_file),
        }
        monitor.report_end(tool_name, result)
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise
    finally:
        if "running_marker" in locals():
            running_marker.unlink(missing_ok=True)


if  __name__ == "__main__":
    analyze_chassis_metabolites.invoke({})
