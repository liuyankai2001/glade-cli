from __future__ import annotations

import ast
import json
import runpy
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

import cobra
import pandas as pd
from cobra.util.solver import fix_objective_as_constraint


CHASSIS_ANALYSIS_RUNNING_MARKER = ".analyze_chassis_metabolites.running"
PROJECT_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_MODEL_PATH = PROJECT_ROOT / "data" / "gem_models" / "iML1515.json"
DEFAULT_MEDIUM_PATH = PROJECT_ROOT / "data" / "mediums" / "default_medium.json"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs" / "chassis_result"
DEFAULT_GROWTH_FRACTION = 0.1
DEFAULT_FLUX_THRESHOLD = 1e-8
DEFAULT_TEST_COMPARTMENTS = ("c",)


@dataclass(frozen=True)
class ChassisAnalysisConfig:
    model_path: Path
    medium_path: Path
    output_dir: Path
    producible_csv: Path
    summary_csv: Path
    growth_fraction: float = DEFAULT_GROWTH_FRACTION
    flux_threshold: float = DEFAULT_FLUX_THRESHOLD
    test_compartments: tuple[str, ...] | None = DEFAULT_TEST_COMPARTMENTS


def _config_value(source: Any, *names: str, default: Any = None) -> Any:
    for name in names:
        if isinstance(source, Mapping) and name in source:
            return source[name]
        if hasattr(source, name):
            return getattr(source, name)
    return default


def _config_value_or_default(source: Any, *names: str, default: Any) -> Any:
    value = _config_value(source, *names, default=default)
    return default if value is None else value


def _resolve_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


def _normalize_compartments(value: Any) -> tuple[str, ...] | None:
    if value is None:
        return None
    if isinstance(value, str):
        raw_values = value.replace(",", " ").split()
    else:
        raw_values = value
    compartments = tuple(
        dict.fromkeys(str(item).strip() for item in raw_values if str(item).strip())
    )
    return compartments or None


def _coerce_config(source: Any) -> ChassisAnalysisConfig:
    model_path = _resolve_path(
        _config_value_or_default(source, "model_path", default=DEFAULT_MODEL_PATH)
    )
    medium_path = _resolve_path(
        _config_value_or_default(source, "medium_path", default=DEFAULT_MEDIUM_PATH)
    )
    target_name = _config_value(source, "target_name")
    default_output_dir = DEFAULT_OUTPUT_DIR
    if target_name:
        default_output_dir = (
            PROJECT_ROOT / "outputs" / str(target_name) / "chassis_result"
        )
    output_dir = _resolve_path(
        _config_value_or_default(
            source,
            "chassis_output_path",
            "output_dir",
            default=default_output_dir,
        )
    )
    producible_csv = _resolve_path(
        _config_value_or_default(
            source,
            "chassis_producible_csv",
            "producible_csv",
            default=output_dir / "producible_kegg_compounds.csv",
        )
    )
    summary_csv = _resolve_path(
        _config_value_or_default(
            source,
            "chassis_metabolites_summary_csv",
            "summary_csv",
            default=output_dir / "analyze_chassis_metabolites_summary.csv",
        )
    )
    return ChassisAnalysisConfig(
        model_path=model_path,
        medium_path=medium_path,
        output_dir=output_dir,
        producible_csv=producible_csv,
        summary_csv=summary_csv,
        growth_fraction=float(
            _config_value_or_default(
                source,
                "growth_fraction",
                default=DEFAULT_GROWTH_FRACTION,
            )
        ),
        flux_threshold=float(
            _config_value_or_default(
                source,
                "flux_threshold",
                default=DEFAULT_FLUX_THRESHOLD,
            )
        ),
        test_compartments=_normalize_compartments(
            _config_value(
                source,
                "test_compartments",
                "compartments",
                default=DEFAULT_TEST_COMPARTMENTS,
            )
        ),
    )


def _load_python_config(path: Path) -> Any:
    try:
        namespace = runpy.run_path(str(path))
    except Exception as execution_error:
        # 兼容只包含 RunConfig(...) 和简单属性赋值的声明式配置文件。
        # 当配置模块本身暂时不可导入时，仍可读取其中的字面量参数。
        try:
            return _load_static_python_config(path)
        except Exception:
            raise RuntimeError(f"Failed to load Python config: {path}") from execution_error

    preferred_names = ("CONFIG", "run_config", "config")
    for name in preferred_names:
        candidate = namespace.get(name)
        if candidate is not None and _config_value(candidate, "model_path") is not None:
            return candidate

    candidates = [
        value
        for name, value in namespace.items()
        if not name.startswith("__")
        and _config_value(value, "model_path") is not None
        and _config_value(value, "medium_path") is not None
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError(f"No chassis configuration object found in: {path}")
    raise ValueError(
        f"Multiple chassis configuration objects found in {path}; name the intended one CONFIG."
    )


def _literal_value(node: ast.AST) -> Any:
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError) as exc:
        raise ValueError("Only literal values are supported in static Python configs") from exc


def _load_static_python_config(path: Path) -> dict[str, Any]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    configs: dict[str, dict[str, Any]] = {}

    for statement in tree.body:
        if isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Call):
            function = statement.value.func
            function_name = function.id if isinstance(function, ast.Name) else None
            if function_name != "RunConfig" or len(statement.targets) != 1:
                continue
            target = statement.targets[0]
            if not isinstance(target, ast.Name):
                continue
            configs[target.id] = {
                keyword.arg: _literal_value(keyword.value)
                for keyword in statement.value.keywords
                if keyword.arg is not None
            }
            continue

        if isinstance(statement, ast.Assign) and len(statement.targets) == 1:
            target = statement.targets[0]
            if not isinstance(target, ast.Attribute) or not isinstance(target.value, ast.Name):
                continue
            config_values = configs.get(target.value.id)
            if config_values is not None:
                config_values[target.attr] = _literal_value(statement.value)

    if "CONFIG" in configs:
        return configs["CONFIG"]
    if len(configs) == 1:
        return next(iter(configs.values()))
    if not configs:
        raise ValueError(f"No RunConfig declaration found in: {path}")
    raise ValueError(
        f"Multiple RunConfig declarations found in {path}; name the intended one CONFIG."
    )


def load_chassis_config(source: Any = None) -> ChassisAnalysisConfig:
    """把 argparse Namespace、配置对象或配置文件转换为底盘分析配置。"""

    if isinstance(source, ChassisAnalysisConfig):
        return source
    if source is None:
        return _coerce_config({})

    # 支持直接传入 argparse Namespace；其 input 可以为空、JSON 或 Python 配置文件。
    if not isinstance(source, (str, Path, Mapping)) and hasattr(source, "input"):
        input_path = getattr(source, "input")
        if input_path is not None:
            return load_chassis_config(input_path)
        return _coerce_config(source)

    if isinstance(source, (str, Path)):
        path = _resolve_path(source)
        if not path.is_file():
            raise FileNotFoundError(f"Input configuration file not found: {path}")
        if path.suffix.lower() == ".json":
            with path.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
            if not isinstance(payload, dict):
                raise ValueError(f"Configuration JSON root must be an object: {path}")
            return _coerce_config(payload)
        if path.suffix.lower() == ".py":
            return _coerce_config(_load_python_config(path))
        raise ValueError("Input configuration must be a .json or .py file")

    return _coerce_config(source)

def safe_mkdir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def save_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
        return pd.DataFrame(
            columns=["source", "met_id", "met_name", "compartment", "kegg_id"]
        )

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

def _validate_config(config: ChassisAnalysisConfig) -> None:
    if not config.model_path.is_file():
        raise FileNotFoundError(f"GEM model file not found: {config.model_path}")
    if config.model_path.suffix.lower() != ".json":
        raise ValueError(
            f"Only JSON GEM models are currently supported: {config.model_path}"
        )
    if not config.medium_path.is_file():
        raise FileNotFoundError(f"Medium file not found: {config.medium_path}")
    if not 0 < config.growth_fraction <= 1:
        raise ValueError("growth_fraction must be greater than 0 and at most 1")
    if config.flux_threshold <= 0:
        raise ValueError("flux_threshold must be greater than 0")


def _load_medium(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Medium JSON root must be an object: {path}")

    medium: dict[str, float] = {}
    for reaction_id, bound in payload.items():
        try:
            medium[str(reaction_id)] = float(bound)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Medium bound for {reaction_id!r} must be numeric, got: {bound!r}"
            ) from exc
    return medium


def analyze_chassis_metabolites(config_source: Any = None) -> dict[str, Any]:
    """
    分析底盘 GEM 在给定培养基下可产生的 KEGG 化合物集合。
    
    调用时机：通路 gap 分析前需要刷新 chassis reachability。
    返回：成功状态、基线生长率、分析阈值、可产物数量及两个输出文件路径。
    限制：只读模型并写分析输出；不写 design_manifest。
    """

    tool_name = "analyze_chassis_metabolites"
    try:
        config = load_chassis_config(config_source)
        _validate_config(config)
        outdir = safe_mkdir(config.output_dir)
        running_marker = outdir / CHASSIS_ANALYSIS_RUNNING_MARKER
        running_marker.write_text(
            json.dumps(
                {
                    "tool": tool_name,
                    "model_path": str(config.model_path),
                    "medium_path": str(config.medium_path),
                    "output_dir": str(config.output_dir),
                    "growth_fraction": config.growth_fraction,
                    "flux_threshold": config.flux_threshold,
                    "test_compartments": list(config.test_compartments or ()),
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        model = cobra.io.load_json_model(str(config.model_path))
        model.medium = _load_medium(config.medium_path)
        growth = model.slim_optimize()
        producible_df = get_producible_metabolites(
            model,
            growth_fraction=config.growth_fraction,
            flux_threshold=config.flux_threshold,
            compartments=(
                set(config.test_compartments)
                if config.test_compartments is not None
                else None
            ),
        )
        producible_kegg_df = dataframe_to_kegg_table(producible_df, "producible")
        reachable_kegg_df = producible_kegg_df.copy()
        reachable_kegg_file = config.producible_csv
        summary_file = config.summary_csv
        save_table(producible_kegg_df, reachable_kegg_file)
        analyze_chassis_metabolites_summary = pd.DataFrame(
            [
                {"item": "baseline_growth", "value": float(growth)},
                {"item": "growth_fraction", "value": config.growth_fraction},
                {"item": "flux_threshold", "value": config.flux_threshold},
                {"item": "test_compartments", "value": ";".join(config.test_compartments or ()) or "non_external"},
                {"item": "producible_metabolites", "value": int(len(producible_df))},
                {"item": "producible_kegg_compounds", "value": int(len(producible_kegg_df))},
                {"item": "reachable_kegg_compounds", "value": int(len(reachable_kegg_df))},
            ]
        )
        save_table(analyze_chassis_metabolites_summary, summary_file)

        result = {
            "ok": True,
            "baseline_growth": float(growth),
            "growth_fraction": config.growth_fraction,
            "flux_threshold": config.flux_threshold,
            "test_compartments": list(config.test_compartments or ()),
            "producible_metabolites": int(len(producible_df)),
            "producible_kegg_compounds": int(len(producible_kegg_df)),
            "reachable_kegg_file": str(reachable_kegg_file),
            "summary_file": str(summary_file),
        }
        return result
    finally:
        if "running_marker" in locals():
            running_marker.unlink(missing_ok=True)


def run_chassis(args: Any) -> dict[str, Any]:
    """供 ``main.py`` 的 argparse handler 直接调用。"""

    result = analyze_chassis_metabolites(args)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result
