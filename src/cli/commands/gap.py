from src.pathway_analyze.kegg_gap_analyze import run_gap

def register(subparsers):
    p = subparsers.add_parser('gap', help='底盘细胞通路分析')
    p.add_argument('-i', '--input',required=True, type=str, help='输入配置文件')
    p.add_argument('-d', '--depth', type=int, default=0, help='使用的底盘扩展深度；0 表示原始 A0')
    p.set_defaults(func=run_gap)
    return p