from src.pathway_analyze.analyze_chassis_metabolites import run_chassis

def register(subparsers):
    p = subparsers.add_parser('chassis', help='底盘细胞可提供代谢物分析')
    p.add_argument('-i', '--input',required=True, type=str, help='输入配置文件')
    p.set_defaults(func=run_chassis)
    return p