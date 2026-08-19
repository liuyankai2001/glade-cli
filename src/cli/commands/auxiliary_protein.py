from src.protein_selection import run_auxiliary_protein_research

def register(subparsers):
    p = subparsers.add_parser(
        "auxiliary-protein",
        help="研究已选主酶组合可能需要的辅助蛋白",
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的输入配置文件名",
    )
    p.add_argument(
        "--research-mode",
        choices=("balanced", "deep"),
        default="balanced",
        help="研究深度，默认 balanced；deep 更全面但耗时更长",
    )
    p.set_defaults(func=run_auxiliary_protein_research)
    return p