from src.pathway_analyze.expand_chassis_metabolites import run_expand

def register(subparsers):
    p = subparsers.add_parser('expand', help='根据 KEGG 反应分层扩展底盘可提供代谢物集合')
    p.add_argument('-d', '--depth', type=int, help='扩展深度，必须大于等于 1')
    p.add_argument('-i', '--input', required=True, help='inputs 目录下的输入配置文件名')
    p.set_defaults(func=run_expand)
    return p