"""CLI entry points for protein selection, sequence downloads and CDS optimization."""
from __future__ import annotations

import json
import sys
from argparse import ArgumentTypeError
from pathlib import PureWindowsPath

from src.cli.commands import optimize
from src.protein_selection import (
    add_protein_to_reaction, download_protein_sequence, remove_protein_from_reaction, replace_protein_in_reaction,
)


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
    parser = subparsers.add_parser("protein", help="下载蛋白序列，添加、删除、替换已选蛋白及优化 CDS")
    actions = parser.add_subparsers(dest="protein_action", required=True)
    download = actions.add_parser("download", help="按 UniProt 编号或名称下载 FASTA 和 JSON")
    download.add_argument("-i", "--input", required=True, type=_input_filename, help="inputs 目录下的配置文件名")
    source = download.add_mutually_exclusive_group(required=True)
    source.add_argument("--uniprot", help="指定 UniProt accession，例如 P22873")
    source.add_argument("--name", help="基因名或蛋白名称；按当前底盘分类接近度自动选择")
    download.add_argument("--force", action="store_true", help="重新获取并替换同编号的 FASTA 和 JSON")
    download.set_defaults(func=_run_download)
    add = actions.add_parser("add", help="为当前路线反应追加本地已选蛋白")
    add.add_argument("-i", "--input", required=True, type=_input_filename, help="inputs 下的配置文件名")
    add.add_argument("-r", "--reaction", dest="reaction_id", required=True, help="当前路线中的反应编号")
    add.add_argument("-f", "--file", dest="sequence_file", required=True, help="inputs/protein_sequences 下的 .fasta 文件名")
    add.set_defaults(func=_run_add)
    remove = actions.add_parser("remove", help="删除当前路线反应中的一条已选蛋白关联，保留序列文件")
    remove.add_argument("-i", "--input", required=True, type=_input_filename, help="inputs 下的配置文件名")
    remove.add_argument("-r", "--reaction", dest="reaction_id", required=True, help="当前路线中的反应编号")
    remove.add_argument("-p", "--protein", dest="accession", required=True, help="UniProt 或本地蛋白编号")
    remove.set_defaults(func=_run_remove)
    replace = actions.add_parser("replace", help="替换反应中的一个已选蛋白，保留其他蛋白及原文件")
    replace.add_argument("-i", "--input", required=True, type=_input_filename, help="inputs 下的配置文件名")
    replace.add_argument("-r", "--reaction", dest="reaction_id", required=True, help="当前路线中的反应编号")
    replace.add_argument("-p", "--protein", dest="old_accession", required=True, help="该反应已选中的旧蛋白编号")
    replace.add_argument("-f", "--file", dest="sequence_file", required=True, help="inputs/protein_sequences 下的新 .fasta 文件名")
    replace.set_defaults(func=_run_replace)
    optimize.register(actions)
    return parser


def _run_add(config):
    try:
        result = add_protein_to_reaction(config, reaction_id=config.reaction_id, sequence_file=config.sequence_file)
    except (ValueError, OSError) as exc:
        print(json.dumps({"ok": False, "message": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(2) from exc
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def _run_remove(config):
    try:
        result = remove_protein_from_reaction(config, reaction_id=config.reaction_id, accession=config.accession)
    except (ValueError, OSError) as exc:
        print(json.dumps({"ok": False, "message": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(2) from exc
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def _run_replace(config):
    try:
        result = replace_protein_in_reaction(config, reaction_id=config.reaction_id,
                                            old_accession=config.old_accession, sequence_file=config.sequence_file)
    except (ValueError, OSError) as exc:
        print(json.dumps({"ok": False, "message": str(exc)}, ensure_ascii=False), file=sys.stderr)
        raise SystemExit(2) from exc
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result
