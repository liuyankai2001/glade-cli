from src.main_protein_selection import run_main_protein_selection

def register(subparsers):
    p = subparsers.add_parser('main-enzyme', help='为 manifest 中已选路线检索和排序主酶候选')
    p.add_argument('-i', '--input', required=True, help='inputs 目录下的输入配置文件名')
    p.add_argument('-n', '--top-n', type=int, default=5, help='每个反应最终保留的主酶候选数量，默认 5')
    p.add_argument(
        '--literature-search',
        action='store_true',
        default=False,
        help='为标准数据库未覆盖的步骤在线检索论文实验酶活证据',
    )
    p.set_defaults(func=run_main_protein_selection)
    return p
