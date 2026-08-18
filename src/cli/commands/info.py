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

    p.add_argument('-d', '--depth',type=int, default=0, help='指定信息查看深度，具体含义由查看类型决定')
    p.add_argument(
        "--step",
        type=int,
        metavar="N",
        default=None,
        help="查看路线中的第 N 步，必须与 --solution 一起使用",
    )
    p.set_defaults(func=run_info)
    return p