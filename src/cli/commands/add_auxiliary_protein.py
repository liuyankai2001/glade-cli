"""CLI registration for manually uploaded auxiliary sequences."""

from src.write_manifest.manual_auxiliary_protein import (
    run_add_manual_auxiliary_protein,
)


def register(subparsers):
    parser = subparsers.add_parser(
        "add-auxiliary-protein",
        help="从 inputs 导入用户提供的辅助蛋白氨基酸或 CDS 序列",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的项目输入配置文件名",
    )
    parser.add_argument(
        "--protein-file",
        required=True,
        help="inputs 目录下的 FASTA、FAA 或纯文本序列文件名",
    )
    parser.add_argument(
        "--sequence-type",
        required=True,
        choices=("protein", "cds"),
        help="protein 执行密码子优化；cds 直接使用并跳过优化",
    )
    parser.set_defaults(func=run_add_manual_auxiliary_protein)
    return parser


__all__ = ["register"]
