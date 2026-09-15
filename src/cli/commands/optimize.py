"""CLI entry point for explicitly requested CDS edits."""

from src.protein_to_cds.restriction_optimization import run_cds_optimization


def register(subparsers):
    parser = subparsers.add_parser("optimize", help="调整单条 CDS 的 GC，或批量消除用户指定的限制酶位点")
    parser.add_argument("-i", "--input", required=True, help="inputs 目录下的输入配置文件名")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--cds", metavar="ID", help="GC 模式要优化的蛋白编号，例如 P21683")
    mode.add_argument("--enzyme", nargs="+", metavar="NAME", help="批量消除所有 CDS 的指定酶位点，忽略大小写；none 清空选择")
    parser.add_argument("--gc-min", metavar="PERCENT", help="GC 模式必填下限百分比；传 --window 时用于局部窗口")
    parser.add_argument("--gc-max", metavar="PERCENT", help="GC 模式必填上限百分比；传 --window 时用于局部窗口")
    parser.add_argument("--window", type=int, default=None, metavar="NT", help="局部 GC 窗口长度，例如 50；不传则调整整体 GC")
    parser.set_defaults(func=run_cds_optimization)
    return parser
