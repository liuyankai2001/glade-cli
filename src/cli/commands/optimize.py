"""CLI entry point for explicitly requested CDS edits."""

from src.protein_to_cds.gc_optimization import run_gc_optimization


def register(subparsers):
    parser = subparsers.add_parser("optimize", help="按指定范围调整单条 CDS 的整体或局部 GC")
    parser.add_argument("-i", "--input", required=True, help="inputs 目录下的输入配置文件名")
    parser.add_argument("--cds", required=True, metavar="ID", help="要优化的蛋白编号，例如 P21683")
    parser.add_argument("--gc-min", required=True, metavar="PERCENT", help="GC 下限百分比；传 --window 时用于局部窗口")
    parser.add_argument("--gc-max", required=True, metavar="PERCENT", help="GC 上限百分比；传 --window 时用于局部窗口")
    parser.add_argument("--window", type=int, default=None, metavar="NT", help="局部 GC 窗口长度，例如 50；不传则调整整体 GC")
    parser.set_defaults(func=run_gc_optimization)
    return parser
