import json
from argparse import ArgumentParser
from src.pathway_analyze.analyze_chassis_metabolites import run_chassis
from src.pathway_analyze.kegg_gap_analyze import run_gap
from src.config.run_config import RunConfig
from src.config.config import INPUTS_DIR

def load_input(path):
    with open(path,'r',encoding='utf-8') as f:
        target_name = json.load(f).get('target_name')

    config = RunConfig(target_name=target_name)
    return config
    

parser = ArgumentParser()
sub_parser = parser.add_subparsers(dest='command',required=True)

chassis_parser = sub_parser.add_parser('chassis',help='底盘细胞可提供代谢物分析')
chassis_parser.add_argument('-i','--input',type=str,help='输入配置文件')
chassis_parser.set_defaults(func=run_chassis)

gap_parser = sub_parser.add_parser('gap',help='底盘细胞通路分析')
gap_parser.add_argument('-i','--input',type=str,help='输入配置文件')
gap_parser.set_defaults(func=run_gap)

def main():
    args = parser.parse_args()
    config_path = INPUTS_DIR / args.input
    config = load_input(config_path)
    args.func(config)

if __name__ == '__main__':
    main()