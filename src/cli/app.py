from argparse import ArgumentParser

from src.cli.common import load_input, apply_args_to_config
from src.config.config import INPUTS_DIR
from src.cli.commands import (
    chassis, expand, gap, validate, select, main_enzyme, info,
)

COMMAND_MODULES = [chassis, expand, gap, validate, select, main_enzyme, info]

def build_parser():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='command', required=True)
    for module in COMMAND_MODULES:
        module.register(subparsers)
    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()

    config_path = INPUTS_DIR / args.input
    config = load_input(config_path)
    config = apply_args_to_config(config, args)

    args.func(config)
