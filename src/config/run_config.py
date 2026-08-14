import sys
import os
from pathlib import Path
from dataclasses import dataclass
from src.config import config

class RunConfig:
    def __init__(self, target_name: str = None):
        # 基础配置
        self.root_dir: Path = config.ROOT
        self.data_dir: Path = config.DATA_DIR
        self.gem_dir: Path = config.GEM_DIR
        self.medium_dir: Path = config.MEDIUM_DIR
        self.inputs_dir: Path = config.INPUTS_DIR
        self.outputs_dir: Path = config.OUTPUTS_DIR
        self.cache_dir: Path = config.CACHE_DIR

        self.model_path: Path = self.gem_dir / 'iML1515.json'
        self.medium_path: Path = self.medium_dir / 'default_medium.json'

        self.target_name: str = target_name
        self.project_output_path: Path = self.outputs_dir / self.target_name

        # 底盘细胞分析输出文件
        self.chassis_output_path: Path = self.project_output_path / 'chassis_result'
        self.chassis_producible_csv: Path = self.chassis_output_path / 'producible_kegg_compounds.csv'
        self.chassis_metabolites_summary_csv: Path = self.chassis_output_path / 'analyze_chassis_metabolites_summary.csv'

