import json
from src.config.run_config import RunConfig
from src.config.config import INPUTS_DIR

def load_input(path):
    with open(path, 'r', encoding='utf-8') as f:
        target_name = json.load(f).get('target_name')
    return RunConfig(target_name=target_name)


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