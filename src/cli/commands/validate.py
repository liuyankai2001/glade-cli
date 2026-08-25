from src.pathway_analyze.gem_validation import run_validation

def register(subparsers):
    p = subparsers.add_parser('validate', help='通路通量分析')
    p.add_argument('-i', '--input', required=True,type=str, help='输入配置文件')
    selection_group = p.add_mutually_exclusive_group()
    selection_group.add_argument(
        '-s',
        '--solutions',
        type=int,
        nargs='+',
        help='指定需要验证的 KEGG solution ID',
    )
    selection_group.add_argument(
        '--retropath-candidates',
        type=int,
        nargs='*',
        metavar='N',
        default=None,
        help='验证 RetroPath 候选；不提供编号时验证全部候选',
    )
    p.add_argument('-m', '--mode', dest='validation_mode',
                   choices=('per', 'pooled', 'both'), default='per', help='路径验证模式')
    p.add_argument('-c', '--cofactor-mode', dest='validation_cofactor_mode',
                   choices=('strict', 'relaxed'), default='strict', help='辅因子模式')
    p.add_argument('-d', '--depth', type=int, default=0, help='需要验证的 gap 结果深度，默认 0')
    p.set_defaults(func=run_validation)
    return p
