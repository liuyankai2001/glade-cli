from src.info_show import run_info

def register(subparsers):
    p = subparsers.add_parser('info', help='信息查看')
    p.add_argument('-i', '--input', required=True, help='inputs 目录下的输入配置文件名')
    p.add_argument('--chassis',dest="show_chassis", action='store_true', help='查看底盘模型、培养基和可生成代谢物摘要')
    p.add_argument('-d', '--depth', dest="info_depth",type=int, default=None, help='指定信息查看深度，具体含义由查看类型决定')
    p.set_defaults(func=run_info)
    return p