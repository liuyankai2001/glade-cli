import sys
import os
from pathlib import Path
from dataclasses import dataclass
from src.config import config

@dataclass
class RunConfig:
    # 基础配置
    root_dir:Path = config.ROOT
    data_dir:Path = config.DATA_DIR
    gem_dir:Path = config.GEM_DIR
    medium_dir:Path = config.MEDIUM_DIR
    inputs_dir:Path = config.INPUTS_DIR
    outputs_dir:Path = config.OUTPUTS_DIR
    cache_dir:Path = config.CACHE_DIR

    model_path:Path = gem_dir / 'iML1515.json'
    medium_path:Path = medium_dir / 'default_medium.json'

    target_name:str = None
    project_output_path:Path = outputs_dir / target_name

    # 底盘细胞分析输出文件
    chassis_output_path:Path = project_output_path / 'chassis_result'
    chassis_producible_csv:Path = chassis_output_path / 'producible_kegg_compounds.csv'
    chassis_metabolites_summary_csv:Path = chassis_output_path / 'analyze_chassis_metabolites_summary.csv'

    