"""Build a GEM through the pre-deployed Bactabolize HTTP service."""


def run_gem_build(config):
    from src.pathway_analyze.gem_build.pipeline import run_gem_build as run

    return run(config)


def register(subparsers):
    parser = subparsers.add_parser('gem-build', help='通过 Bactabolize 服务重建底盘 GEM')
    parser.add_argument('-i', '--input', required=True, help='inputs 目录下的输入配置文件名')
    parser.add_argument('--force', action='store_true', help='重新提交重建任务，保留历史结果')
    parser.set_defaults(func=run_gem_build)
    return parser
