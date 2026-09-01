from src.pathway_analyze.retropath_pipeline import run_gap_command

def register(subparsers):
    p = subparsers.add_parser('gap', help='底盘细胞通路分析')
    p.add_argument('-i', '--input',required=True, type=str, help='输入配置文件')
    p.add_argument('-d', '--depth', type=int, default=0, help='使用的底盘扩展深度；0 表示原始 A0')
    p.add_argument(
        '--retropath',
        action='store_true',
        default=None,
        help='显式启用 RetroPath 预测搜索；默认仍使用原始 KEGG 搜索',
    )
    p.add_argument(
        '--step',
        dest='retropath_max_steps',
        type=int,
        choices=range(1, 11),
        default=None,
        metavar='N',
        help='RetroPath 最大逆合成步数（1–10；默认 3；仅与 --retropath 一起使用）',
    )
    p.set_defaults(func=run_gap_command)
    return p
