from pathlib import Path

ROOT = Path('__file__').parent.parent.parent
DATA_DIR = ROOT / 'data'
GEM_DIR = DATA_DIR / 'gem_models'
MEDIUM_DIR = DATA_DIR / 'mediums'

INPUTS_DIR = ROOT / 'inputs'
OUTPUTS_DIR = ROOT / 'outputs'
CACHE_DIR = ROOT / 'cache'

if __name__ == "__main__":
    print(DATA_DIR)