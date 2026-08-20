"""CLI registration for expression-box design."""

from __future__ import annotations

import json
from typing import Any

from src.expression_box import run_expression_box_design


def _run(config: Any) -> dict[str, Any]:
    result = run_expression_box_design(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def register(subparsers):
    parser = subparsers.add_parser(
        "expression",
        help="设计表达盒或表达元件方案",
    )
    parser.add_argument(
        "--design",
        action="store_true",
        required=True,
        help="进入表达设计流程",
    )
    parser.add_argument(
        "--box",
        action="store_true",
        required=True,
        help="生成系统推荐的蛋白表达盒分组方案",
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
