import json
from argparse import ArgumentParser
from src.pathway_analyze.analyze_chassis_metabolites import run_chassis
from src.pathway_analyze.kegg_gap_analyze import run_gap
from src.pathway_analyze.gem_validation import run_validation
from src.config.run_config import RunConfig
from src.config.config import INPUTS_DIR

def load_input(path):
    with open(path,'r',encoding='utf-8') as f:
        target_name = json.load(f).get('target_name')

    config = RunConfig(target_name=target_name)
    return config
    
# 主parser
parser = ArgumentParser()
sub_parser = parser.add_subparsers(dest='command',required=True)

# 底盘细胞可提供代谢物
chassis_parser = sub_parser.add_parser('chassis',help='底盘细胞可提供代谢物分析')
chassis_parser.add_argument('-i','--input',type=str,help='输入配置文件')
chassis_parser.set_defaults(func=run_chassis)

# gap分析
gap_parser = sub_parser.add_parser('gap',help='底盘细胞通路分析')
gap_parser.add_argument('-i','--input',type=str,help='输入配置文件')
gap_parser.set_defaults(func=run_gap)

# gem通量验证
validation_parser = sub_parser.add_parser('validate')
validation_parser.add_argument('-i','--input',type=str,help='输入配置文件')
validation_parser.add_argument('-s','--solutions',type=int,nargs='+',help='指定需要验证的solution ID')
validation_parser.add_argument(
    '-m',
    "--mode",
    dest="validation_mode",
    choices=("per", "pooled", "both"),
    default="per",
    help="路径验证模式",
)

validation_parser.add_argument(
    '-c',
    "--cofactor-mode",
    dest="validation_cofactor_mode",
    choices=("strict", "relaxed"),
    default="strict",
    help="辅因子模式",
)
validation_parser.set_defaults(func=run_validation)


def main():
    args = parser.parse_args()
    config_path = INPUTS_DIR / args.input
    config = load_input(config_path)

    if args.command == "validate":
        config.solutions = args.solutions
        config.validation_mode = args.validation_mode
        config.validation_cofactor_mode = args.validation_cofactor_mode

    args.func(config)

if __name__ == '__main__':
    main()