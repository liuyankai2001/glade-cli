from src.write_manifest.solution import run_write_solution


def register(subparsers):
    p = subparsers.add_parser(
        "write",
        help="将所选信息写入 design manifest",
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="inputs 目录下的输入配置文件名",
    )
    p.add_argument(
        '-s',
        "--solution",
        type=int,
        required=True,
        metavar="N",
        help="写入第 N 条候选路线",
    )
    p.add_argument(
        "-d",
        "--depth",
        type=int,
        default=0,
        metavar="N",
        help="该路线所属的 gap 深度，默认 0",
    )
    p.set_defaults(func=run_write_solution)
    return p
