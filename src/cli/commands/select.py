from src.pathway_analyze.write_solution_to_manifest import run_select_solution

def register(subparsers):
    p = subparsers.add_parser('select', help='将指定 solution 写入 manifest')
    p.add_argument('-i', '--input', required=True)
    p.add_argument('-s', '--solution',dest="solution_id", type=int, required=True, help='需要保存的solution ID')
    p.add_argument('-d', '--depth', type=int, default=0, help='写入选择的深度')
    p.set_defaults(func=run_select_solution)
    return p