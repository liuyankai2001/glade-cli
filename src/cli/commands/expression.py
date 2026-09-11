"""CLI registration for expression-box design."""

from __future__ import annotations

import json
from argparse import ArgumentTypeError
from typing import Any

from src.expression_box.config import (
    MAX_EXPRESSION_PARTS_DESIGN_COUNT,
    MIN_EXPRESSION_PARTS_DESIGN_COUNT,
)
from src.expression_box import (
    run_expression_box_design,
    run_expression_parts_design,
)
from src.write_manifest.expression_box import write_custom_expression_box_selection


def _parts_design_count(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ArgumentTypeError("--n-designs must be an integer") from exc
    if not MIN_EXPRESSION_PARTS_DESIGN_COUNT <= parsed <= MAX_EXPRESSION_PARTS_DESIGN_COUNT:
        raise ArgumentTypeError(
            "--n-designs must be between "
            f"{MIN_EXPRESSION_PARTS_DESIGN_COUNT} and "
            f"{MAX_EXPRESSION_PARTS_DESIGN_COUNT}"
        )
    return parsed


def _run(config: Any) -> dict[str, Any]:
    custom = getattr(config, "custom", None)
    if bool(getattr(config, "parts", False)):
        if custom is not None:
            raise ValueError("--custom can only be used with --box")
        result = run_expression_parts_design(config)
    else:
        if getattr(config, "n_designs", None) is not None:
            raise ValueError("--n-designs can only be used with --parts")
        if custom is not None:
            result = write_custom_expression_box_selection(config)
        else:
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
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--box",
        action="store_true",
        help="生成系统推荐的蛋白表达盒分组方案",
    )
    mode.add_argument(
        "--parts",
        action="store_true",
        help="从远端Milvus生成系统推荐的表达元件方案",
    )
    parser.add_argument(
        "--custom",
        nargs="+",
        default=None,
        metavar="GROUP",
        help=(
            "自定义表达盒分组并直接写入manifest，例如 "
            "--custom [P00001 P00002] [P00003]"
        ),
    )
    parser.add_argument(
        "--n-designs",
        type=_parts_design_count,
        default=None,
        metavar="N",
        help="稳定表达元件方案数量（默认12，范围3-96，仅用于--parts）",
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
