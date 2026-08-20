"""CLI registration for the manifest-driven protein-to-CDS workflow."""

from __future__ import annotations

import json
from typing import Any

from src.protein_to_cds import run_protein_to_cds


def _run(config: Any) -> dict[str, Any]:
    result = run_protein_to_cds(
        config,
        device=getattr(config, "device", "auto"),
        additional_forbidden_motifs=getattr(config, "forbidden_motif", ()),
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if result["status"] != "complete":
        raise SystemExit(2)
    return result


def register(subparsers):
    parser = subparsers.add_parser(
        "protein-to-cds",
        help="读取 manifest 中的已选蛋白并生成密码子优化 CDS",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的输入配置文件名",
    )
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        default="auto",
        help="模型推理设备，默认自动选择",
    )
    parser.add_argument(
        "--forbidden-motif",
        action="append",
        default=[],
        help="额外禁止的 DNA motif；可重复传入",
    )
    parser.set_defaults(func=_run)
    return parser


__all__ = ["register"]
