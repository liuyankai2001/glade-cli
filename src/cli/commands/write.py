from src.write_manifest.auxiliary_protein import (
    run_write_auxiliary_protein_research,
)
from src.write_manifest.main_enzyme import run_write_main_enzyme_set
from src.write_manifest.solution import run_write_solution


def run_write(config):
    if bool(getattr(config, "auxiliary_protein", False)):
        return run_write_auxiliary_protein_research(config)

    if getattr(config, "main_enzyme_set", None) is not None:
        return run_write_main_enzyme_set(config)

    if getattr(config, "solution", None) is not None:
        return run_write_solution(config)

    raise ValueError("未指定需要写入的内容")


def register(subparsers):
    p = subparsers.add_parser(
        "write",
        help="将所选信息写入 design manifest",
    )

    action = p.add_mutually_exclusive_group(required=True)

    action.add_argument("--solution", type=int, metavar="N")
    action.add_argument("--main-enzyme-set", type=int, metavar="N")
    action.add_argument(
        "--auxiliary-protein",
        action="store_true",
        help="将当前主酶组合的全部辅助蛋白研究结果写入manifest",
    )

    p.add_argument("-i", "--input", required=True)
    p.add_argument("-d", "--depth", type=int, default=0)

    p.set_defaults(func=run_write)
    return p
