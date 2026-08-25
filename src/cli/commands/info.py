from src.info_show import run_info

def register(subparsers):
    p = subparsers.add_parser('info', help='信息查看')
    p.add_argument('-i', '--input', required=True, help='inputs 目录下的输入配置文件名')
    view_group = p.add_mutually_exclusive_group(required=True)
    view_group.add_argument('--chassis', action='store_true', help='查看底盘模型、培养基和可生成代谢物摘要')
    view_group.add_argument(
        "--gap",
        action="store_true",
        help="查看 gap 候选路线摘要",
    )
    view_group.add_argument(
        "--solution",
        type=int,
        help="查看指定路线的完整步骤",
    )
    view_group.add_argument(
        "--retropath",
        action="store_true",
        help="查看 RetroPath 运行状态和预测候选摘要",
    )
    view_group.add_argument(
        "--retropath-candidate",
        type=int,
        metavar="N",
        help="查看排名第 N 的 RetroPath 预测候选详情",
    )

    view_group.add_argument(
        "--main-enzyme-candidates",
        action="store_true",
        help="按步骤查看当前已选路线的主酶候选",
    )

    view_group.add_argument(
        "--main-enzyme-candidate",
        type=int,
        metavar="N",
        help="查看指定 Step 中排名第 N 的主酶候选详情",
    )


    view_group.add_argument(
        "--main-enzyme-sets",
        action="store_true",
        help="查看当前路线的主酶组合列表",
    )

    view_group.add_argument(
        "--main-enzyme-set",
        type=int,
        metavar="N",
        help="查看排名第 N 的主酶组合详情",
    )
    view_group.add_argument(
        "--proteins",
        action="store_true",
        help="查看当前 manifest 中的主酶和辅助蛋白",
    )
    view_group.add_argument(
        "--protein",
        metavar="ID",
        help="查看当前 manifest 中指定蛋白的详情",
    )

    p.add_argument('-d', '--depth',type=int, default=0, help='指定信息查看深度，具体含义由查看类型决定')
    p.add_argument(
        "--step",
        type=int,
        metavar="N",
        default=None,
        help=(
            "查看路线或候选中的第 N 步，支持 --solution、"
            "--retropath-candidate 和主酶候选视图"
        ),
    )
    p.set_defaults(func=run_info)
    return p
