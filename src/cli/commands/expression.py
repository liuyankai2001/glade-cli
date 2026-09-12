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
from src.write_manifest.expression_parts_draft import upload_expression_promoter


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
    promoter = getattr(config, "promoter", None)
    custom = getattr(config, "custom", None)
    if promoter is not None:
        if bool(getattr(config, "box", False)) or bool(
            getattr(config, "parts", False)
        ):
            raise ValueError("--promoter 不能与 --box 或 --parts 同时使用")
        if custom is not None or getattr(config, "n_designs", None) is not None:
            raise ValueError("--promoter 不能与 --custom 或 --n-designs 同时使用")
        result = upload_expression_promoter(config)
    elif bool(getattr(config, "parts", False)):
        if custom is not None:
            raise ValueError("--custom can only be used with --box")
        result = run_expression_parts_design(config)
    elif bool(getattr(config, "box", False)):
        if getattr(config, "n_designs", None) is not None:
            raise ValueError("--n-designs can only be used with --parts")
        if custom is not None:
            result = write_custom_expression_box_selection(config)
        else:
            result = run_expression_box_design(config)
    else:
        raise ValueError("--design 必须与 --box 或 --parts 一起使用")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def register(subparsers):
    parser = subparsers.add_parser(
        "expression",
        help="设计表达盒或表达元件方案",
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "--design",
        action="store_true",
        help="进入表达设计流程",
    )
    action.add_argument(
        "--promoter",
        nargs=2,
        metavar=("BOX", "FILE"),
        help="为指定表达盒上传 inputs/parts 下的启动子文件",
    )
    mode = parser.add_mutually_exclusive_group(required=False)
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
