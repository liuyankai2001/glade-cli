from src.main_protein_selection import run_main_enzyme_sets


def register(subparsers):
    p = subparsers.add_parser(
        "main-enzyme-sets",
        help="根据主酶候选生成并排序主酶组合",
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的输入配置文件名",
    )
    p.add_argument(
        "--max-sets",
        type=int,
        default=20,
        help="最多输出多少个主酶组合，默认 20",
    )
    p.add_argument(
        "--max-search-nodes",
        type=int,
        default=1_000_000,
        help="组合搜索节点上限，通常无需修改",
    )
    p.set_defaults(func=run_main_enzyme_sets)
    return p