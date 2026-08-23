"""CLI registration for deleting manually uploaded auxiliary sequences."""

from src.write_manifest.manual_auxiliary_protein import (
    run_remove_manual_auxiliary_proteins,
)


def register(subparsers):
    parser = subparsers.add_parser(
        "remove-auxiliary-protein",
        help="删除手动导入的辅助蛋白及其项目序列快照",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的项目输入配置文件名",
    )
    parser.add_argument(
        "--protein-id",
        action="append",
        required=True,
        help="要删除的手动辅助蛋白 ID；可重复传入",
    )
    parser.set_defaults(func=run_remove_manual_auxiliary_proteins)
    return parser


__all__ = ["register"]
