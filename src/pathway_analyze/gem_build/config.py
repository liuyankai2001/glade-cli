"""Independent configuration for model construction (no target compound required)."""
from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse


def validate_managed_paths(output, sources):
    """The sole published JSON must never overwrite an input or medium file."""
    for source in sources:
        if source.resolve() == output.resolve() or (source.is_file() and output.is_file() and source.samefile(output)):
            raise ValueError(f'output must not alias any input or medium file: {source}')


@dataclass
class GemBuildConfig:
    root_dir: Path
    model_file: Path
    reference_model_file: Path
    reference_genome_file: Path
    target_genome_file: Path
    service_url: str
    parameters: dict
    medium_path: Path
    force: bool = False
    request_timeout_seconds: float = 30.0
    wait_timeout_seconds: float = 1900.0
    poll_interval_seconds: float = 1.0
    get_attempts: int = 3

    @classmethod
    def from_payload(cls, payload, root_dir):
        if not isinstance(payload, dict):
            raise ValueError('gem-build payload must be an object')
        root = Path(root_dir).expanduser().resolve()
        options = payload.get('bactabolize')
        if not isinstance(options, dict):
            raise ValueError('bactabolize must be an object')

        def path(value, name):
            if not isinstance(value, (str, Path)) or not str(value).strip():
                raise ValueError(f'{name} must be an explicit file path')
            result = Path(value).expanduser()
            return (root / result).resolve() if not result.is_absolute() else result.resolve()

        output = path(payload.get('model_file'), 'model_file')
        if output.suffix.lower() != '.json' or output.is_dir():
            raise ValueError('model_file must be a JSON output file')
        inputs = {}
        for field in ('reference_model_file', 'reference_genome_file', 'target_genome_file'):
            source = path(options.get(field), field)
            if source.suffix.lower() != ('.json' if field == 'reference_model_file' else '.gbff'):
                raise ValueError(f'{field} requires ' + ('JSON' if field == 'reference_model_file' else 'GBFF'))
            if not source.is_file():
                raise FileNotFoundError(f'{field}: {source}')
            with source.open('rb') as handle:
                handle.read(1)
            if source == output or (output.exists() and source.samefile(output)):
                raise ValueError('output must not alias any input file')
            inputs[field] = source
        parameters = {'min_pident': 80.0, 'min_coverage': 25.0,
                      'media_type': 'm9', 'atmosphere_type': 'aerobic', 'biomass_reaction_id': None}
        for field in parameters:
            if field in options:
                parameters[field] = options[field]
        for field in ('min_pident', 'min_coverage'):
            value = parameters[field]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value <= 100:
                raise ValueError(f'{field} must be a finite percentage in (0, 100]')
            parameters[field] = float(value)
        media = {'bg11', 'cdm_mendoza', 'lb_carveme', 'lb', 'lb_SEED', 'm9', 'm9_SEED', 'nutrient',
                 'pmm5_mendoza', 'pmm7_mendoza', 'tsa', 'tsa_sheep_blood'}
        for field, allowed in [('media_type', media), ('atmosphere_type', {'aerobic', 'anaerobic'})]:
            if not isinstance(parameters[field], str) or parameters[field] not in allowed:
                raise ValueError(f'{field} must be one of {sorted(allowed)}')
        biomass = parameters['biomass_reaction_id']
        if biomass is not None and (not isinstance(biomass, str) or not biomass.strip()):
            raise ValueError('biomass_reaction_id must be a nonempty string or null')
        if options.get('no_reannotation', True) is not True:
            raise ValueError('no_reannotation is fixed true; supply GBFF annotations')
        url = options.get('service_url', os.environ.get('BACTABOLIZE_SERVICE_URL', 'http://127.0.0.1:8766'))
        if not isinstance(url, str):
            raise ValueError('service_url must be an HTTP URL')
        parsed = urlparse(url)
        if parsed.scheme not in {'http', 'https'} or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment:
            raise ValueError('service_url must be an HTTP URL without credentials/query/fragment')
        medium = path(payload.get('medium_path', 'data/mediums/default_medium.json'), 'medium_path')
        validate_managed_paths(output, [*inputs.values(), medium])
        return cls(root_dir=root, model_file=output, **inputs, service_url=url.rstrip('/'),
                   parameters=parameters, medium_path=medium)
