import json
from src.config.run_config import RunConfig
from src.config.config import INPUTS_DIR
from src.config.config import ROOT
from pathlib import Path

def load_input(path, *, command=None):
    with open(path, 'r', encoding='utf-8') as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError('input JSON must be an object')
    if command == 'gem-build':
        from src.pathway_analyze.gem_build.config import GemBuildConfig

        return GemBuildConfig.from_payload(payload, ROOT)

    from src.main_protein_selection.genome_ec_filter import resolve_genome_path

    value = payload.get('model_file')
    model_path = None
    if value is not None:
        if not isinstance(value, str):
            raise ValueError('model_file must be a string, null, or omitted')
        if value.strip():
            candidate = Path(value.strip()).expanduser()
            model_path = (candidate if candidate.is_absolute() else ROOT / candidate).resolve()
            if not model_path.is_file():
                raise FileNotFoundError(f'model_file does not exist: {model_path}')
    config = RunConfig(target_name=payload.get('target_name'), model_path=model_path)
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
