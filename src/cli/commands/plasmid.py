"""CLI registration for system plasmid-backbone recommendation."""

from __future__ import annotations

import json
from argparse import ArgumentTypeError
from typing import Any

from src.plasmid_selection.config import (
    DEFAULT_CANDIDATE_COUNT,
    MAX_CANDIDATE_COUNT,
    MIN_CANDIDATE_COUNT,
    SUPPORTED_PRIORITIES,
)
from src.plasmid_selection.pipeline import run_plasmid_recommendation


def _candidate_count(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ArgumentTypeError("--n-candidates must be an integer") from exc
    if not MIN_CANDIDATE_COUNT <= parsed <= MAX_CANDIDATE_COUNT:
        raise ArgumentTypeError(
            "--n-candidates must be between "
            f"{MIN_CANDIDATE_COUNT} and {MAX_CANDIDATE_COUNT}"
        )
    return parsed


def _run(config: Any) -> dict[str, Any]:
    result = run_plasmid_recommendation(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def register(subparsers):
    parser = subparsers.add_parser(
        "plasmid",
        help="从远端 Milvus 推荐可覆盖全部表达构建的质粒骨架",
    )
    parser.add_argument(
        "--recommend",
        action="store_true",
        required=True,
        help="生成系统质粒候选",
    )
    parser.add_argument(
        "--n-candidates",
        type=_candidate_count,
        default=DEFAULT_CANDIDATE_COUNT,
        metavar="N",
        help=f"候选数量（默认 {DEFAULT_CANDIDATE_COUNT}）",
    )
    parser.add_argument(
        "--priority",
        choices=SUPPORTED_PRIORITIES,
        default="stability",
        help="骨架选择优先级（默认 stability）",
    )
    parser.add_argument(
        "--preferred-resistance",
        default=None,
        metavar="MARKER",
        help="优先抗性标记，例如 kanamycin",
    )
    parser.add_argument(
        "--exclude-resistance",
        nargs="+",
        default=None,
        metavar="MARKER",
        help="排除一个或多个抗性标记",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的输入配置文件名",
    )
    parser.set_defaults(func=_run)
    return parser


__all__ = ["register"]
