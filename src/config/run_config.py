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

        # 底盘细胞代谢通路分析参数
        self.max_total_steps = 20
        self.max_new_enzymes = self.max_total_steps
        self.max_solutions = 10
        self.max_reactions_per_compound = 15 # 每次反向扩展一个化合物时最多考虑多少反应
        self.max_major_precursors = 4 # 一个候选反应最多允许包含多少种需要继续向前追溯的“主要前体”。
        self.max_routes_per_state = 3 #  相同待解决化合物情况下，允许保留的不同历史路线数
        self.max_module_chain_depth = 20
        self.module_filter_mode = 'prefer' # prefer、strict 或 off
        self.electron_avoidance_mode = 'strict_with_fallback' # off、prefer、strict 或 strict_with_fallback
        self.reaction_resolution_mode = 'strict' # strict 或 audit
        self.gap_output_path = self.project_output_path / f'kegg_gap_{self.target_name}'