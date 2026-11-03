"""Strict host identity and immutable input binding for one GenBank file."""
from __future__ import annotations

import io
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqRecord import SeqRecord

from .host_registry import get_host_preset, register_host


@dataclass(frozen=True)
class GenomeIdentity:
    path: Path
    name: str
    taxon_id: int
    records: tuple[SeqRecord, ...] = field(repr=False, compare=False)



def read_genome_identity(path: str | Path) -> GenomeIdentity:
    path = Path(path).expanduser().resolve()
    data = path.read_bytes()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error', BiopythonParserWarning)
            records = tuple(SeqIO.parse(io.StringIO(data.decode('utf-8')), 'genbank'))
        if not records:
            raise ValueError('no GenBank records')
        identities = set()
        names = []
        for record in records:
            sources = [feature for feature in record.features if feature.type == 'source']
            if not sources:
                raise ValueError(f'missing source annotation: {record.id}')
            for source in sources:
                values = source.qualifiers.get('organism', [])
                if len(values) != 1 or not isinstance(values[0], str) or not values[0].strip():
                    raise ValueError(f'missing or ambiguous organism: {record.id}')
                name = ' '.join(values[0].split())
                taxa = set()
                for value in source.qualifiers.get('db_xref', []):
                    if value.startswith('taxon:'):
                        token = value.removeprefix('taxon:')
                        if not re.fullmatch(r'[0-9]+', token) or int(token) <= 0:
                            raise ValueError(f'invalid taxon ID: {record.id}')
                        taxa.add(int(token))
                if len(taxa) != 1:
                    raise ValueError(f'missing or ambiguous taxon ID: {record.id}')
                names.append(name)
                identities.add((name.casefold(), next(iter(taxa))))
        if len(identities) != 1:
            raise ValueError('conflicting source names or taxon IDs')
    except (ValueError, BiopythonParserWarning) as exc:
        raise ValueError(f'Invalid GenBank genome_file {path}: {exc}') from exc
    _, taxon_id = next(iter(identities))
    identity = GenomeIdentity(path, names[0], taxon_id, records)
    return identity


def resolve_genome_chassis(identity: GenomeIdentity, chassis_key: str | None = None) -> str:
    """Candidate generation only: explicit keys constrain identity; otherwise register."""
    if chassis_key is not None:
        preset = get_host_preset(chassis_key)
        if preset['taxon_id'] != identity.taxon_id:
            raise ValueError(f'genome_file host mismatch / 菌株不匹配: {chassis_key}, taxon={identity.taxon_id}')
        return chassis_key
    key = register_host(identity.name, identity.taxon_id)
    return key
