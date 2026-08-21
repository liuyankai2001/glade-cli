"""CLI registration for final-assembly planning and execution."""

from __future__ import annotations

import json
from typing import Any

from src.final_assemble_execute import run_final_assembly_execute
from src.final_assemble_plan.config import SUPPORTED_METHODS
from src.final_assemble_plan.plan_final_assembly import run_final_assembly_plan


def _run(config: Any) -> dict[str, Any]:
    if bool(getattr(config, "assembly_execute", False)):
        if getattr(config, "assembly_method", None) is not None:
            raise ValueError("--method 只能与 assembly --plan 一起使用")
        return run_final_assembly_execute(config)
    result = run_final_assembly_plan(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def register(subparsers):
    parser = subparsers.add_parser(
        "assembly",
        help="生成或执行最终组装计划",
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "--plan",
        action="store_true",
        help="生成整套最终组装计划",
    )
    action.add_argument(
        "--execute",
        dest="assembly_execute",
        action="store_true",
        help="执行manifest中已接受的整套最终组装计划",
    )
    parser.add_argument(
        "--method",
        dest="assembly_method",
        choices=SUPPORTED_METHODS,
        default=None,
        help="统一指定组装方法；省略时逐 design 自动推荐",
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
