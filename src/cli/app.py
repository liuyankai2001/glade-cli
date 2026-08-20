from argparse import ArgumentParser

from src.cli.commands import (
    auxiliary_protein,
    chassis,
    expression,
    expand,
    gap,
    info,
    main_enzyme,
    main_enzyme_sets,
    protein_to_cds,
    validate,
    write,
)
from src.cli.common import apply_args_to_config, load_input
from src.config.config import INPUTS_DIR

COMMAND_MODULES = [
    chassis,
    expression,
    expand,
    gap,
    validate,
    write,
    main_enzyme,
    info,
    main_enzyme_sets,
    auxiliary_protein,
    protein_to_cds,
]

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
