from src.pathway_analyze.list_solution_steps import run_solution_steps

def register(subparsers):
    p = subparsers.add_parser('solution', help='查看某个具体路线详情')
    p.add_argument('-s', '--solution',dest='solution_id',type=int, help='需要查看的solution ID', required=True)
    p.add_argument('-i', '--input', required=True, help='输入配置文件')
    p.set_defaults(func=run_solution_steps)
    return p