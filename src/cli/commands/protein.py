"""CLI entry point for local protein sequence downloads."""
from __future__ import annotations

import json
import sys
from argparse import ArgumentTypeError
from pathlib import PureWindowsPath

from src.protein_selection import download_protein_sequence


def _input_filename(value: str) -> str:
    if (not value.strip() or value in {".", ".."} or "/" in value or "\\" in value
            or PureWindowsPath(value).drive or "\x00" in value):
        raise ArgumentTypeError("-i 只接受 inputs 目录下的配置文件名")
    return value


def _run_download(config):
    try:
        result = download_protein_sequence(
            config, uniprot=getattr(config, "uniprot", None),
            name=getattr(config, "name", None), force=getattr(config, "force", False),
        )
    except (ValueError, OSError) as exc:
        print(json.dumps({"ok": False, "message": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(2) from exc
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def register(subparsers):
    parser = subparsers.add_parser("protein", help="获取本地蛋白氨基酸序列")
    actions = parser.add_subparsers(dest="protein_action", required=True)
    download = actions.add_parser("download", help="按 UniProt 编号或名称下载 FASTA 和 JSON")
    download.add_argument("-i", "--input", required=True, type=_input_filename, help="inputs 目录下的配置文件名")
    source = download.add_mutually_exclusive_group(required=True)
    source.add_argument("--uniprot", help="指定 UniProt accession，例如 P22873")
    source.add_argument("--name", help="基因名或蛋白名称；按当前底盘分类接近度自动选择")
    download.add_argument("--force", action="store_true", help="重新获取并替换同编号的 FASTA 和 JSON")
    download.set_defaults(func=_run_download)
    return parser
