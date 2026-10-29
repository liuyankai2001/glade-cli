import json
from src.config.run_config import RunConfig
from src.config.config import INPUTS_DIR

def load_input(path):
    from src.main_protein_selection.genome_ec_filter import resolve_genome_path

    with open(path, 'r', encoding='utf-8') as f:
        payload = json.load(f)
    config = RunConfig(target_name=payload.get('target_name'))
    config.genome_file = resolve_genome_path(payload.get('genome_file'), config.root_dir)
    return config


def apply_args_to_config(config, args):
    ignored_fields = {
        "command",
        "func",
        "input",
    }
    for name, value in vars(args).items():
        if name in ignored_fields:
            continue

        # 不用 None 覆盖 RunConfig 中已有的默认值
        if value is None:
            continue
        setattr(config, name, value)

    return config
