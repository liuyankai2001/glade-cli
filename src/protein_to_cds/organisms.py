"""Resolve integer generation categories against the CSV and installed package."""

from __future__ import annotations

import ast
import csv
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path


ORGANISMS_PATH = Path(__file__).resolve().parents[2] / 'data/codontransformer_organisms.csv'
_COLUMNS = ['codon_transformer_organism_id', 'organism_name']


def _installed_organisms() -> dict[int, str]:
    # Reading literal constants avoids loading torch or any model for validation.
    try:
        package = distribution('CodonTransformer')
        source = Path(package.locate_file('CodonTransformer/CodonUtils.py'))
        tree = ast.parse(source.read_text(encoding='utf-8-sig'))
        mapping = next(
            ast.literal_eval(node.value) for node in tree.body
            if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
            and node.target.id == 'ORGANISM2ID'
        )
    except (PackageNotFoundError, OSError, UnicodeError, SyntaxError, ValueError, StopIteration) as exc:
        raise ValueError('Cannot read installed CodonTransformer organism definitions') from exc
    if (not isinstance(mapping, dict) or not mapping
            or any(not isinstance(name, str) or not name.strip()
                   or type(identifier) is not int or identifier < 0
                   for name, identifier in mapping.items())
            or len(set(mapping.values())) != len(mapping)):
        raise ValueError('Invalid installed CodonTransformer organism definitions')
    return {identifier: name for name, identifier in mapping.items()}


def _catalog() -> dict[int, str]:
    path = Path(ORGANISMS_PATH)
    result: dict[int, str] = {}
    names: set[str] = set()
    try:
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle, strict=True)
            if reader.fieldnames != _COLUMNS:
                raise ValueError('unexpected headers')
            for row in reader:
                raw_id, name = row.get(_COLUMNS[0]), row.get(_COLUMNS[1])
                if (set(row) != set(_COLUMNS) or not isinstance(raw_id, str)
                        or not raw_id.isascii() or not raw_id.isdecimal()
                        or not isinstance(name, str) or not name.strip() or name != name.strip()):
                    raise ValueError('invalid row')
                identifier = int(raw_id)
                if identifier in result or name in names:
                    raise ValueError('duplicate ID or name')
                result[identifier] = name
                names.add(name)
        if not result:
            raise ValueError('empty catalog')
    except FileNotFoundError as exc:
        raise FileNotFoundError(f'CodonTransformer organism catalog not found: {path}') from exc
    except (ValueError, UnicodeError, csv.Error) as exc:
        raise ValueError(f'Invalid CodonTransformer organism catalog {path}: {exc}') from exc
    return result


def resolve_organism(value: int) -> tuple[int, str]:
    """Accept JSON integers only, including zero; return the official name."""
    if type(value) is not int:
        raise ValueError('codon_transformer_organism must be an integer model ID or null')
    catalog = _catalog()
    if value not in catalog:
        raise ValueError(f'codon_transformer_organism unsupported ID: {value}; see {ORGANISMS_PATH}')
    if catalog != _installed_organisms():
        raise ValueError('Organism catalog differs from installed CodonTransformer; use matching package and CSV versions')
    return value, catalog[value]
