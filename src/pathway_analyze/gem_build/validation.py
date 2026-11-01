"""Validate the original JSON before COBRA can normalize or merge its objects."""
from __future__ import annotations

import json
import math
import warnings

from .client import GemBuildError


def _finite(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def validate_model(model_path, medium_path):
    import cobra
    from cobra.core import GPR
    try:
        original = json.loads(model_path.read_text(encoding='utf-8'))
        ids = {}
        for group in ('metabolites', 'reactions', 'genes'):
            rows = original[group]
            if not isinstance(rows, list) or (not rows and group != 'genes'):
                raise ValueError(f'{group} must be a list of objects')
            values = [row['id'] for row in rows]
            if any(not isinstance(value, str) or not value.strip() for value in values) or len(values) != len(set(values)):
                raise ValueError(f'duplicate or invalid {group} IDs')
            ids[group] = set(values)
        objective = []
        for reaction in original['reactions']:
            lower, upper = reaction['lower_bound'], reaction['upper_bound']
            if not _finite(lower) or not _finite(upper) or lower > upper:
                raise ValueError(f'invalid bounds: {reaction["id"]}')
            metabolites = reaction['metabolites']
            if not isinstance(metabolites, dict) or set(metabolites) - ids['metabolites'] or any(not _finite(v) for v in metabolites.values()):
                raise ValueError(f'invalid metabolite references: {reaction["id"]}')
            with warnings.catch_warnings():
                warnings.simplefilter('error', SyntaxWarning)
                genes = GPR.from_string(reaction.get('gene_reaction_rule', '')).genes
            if set(genes) - ids['genes']:
                raise ValueError(f'undefined gene reference: {reaction["id"]}')
            coefficient = reaction.get('objective_coefficient', 0)
            if not _finite(coefficient):
                raise ValueError('invalid objective coefficient')
            if coefficient:
                objective.append(reaction['id'])
        if not objective:
            raise ValueError('model has no objective')
        model = cobra.io.load_json_model(str(model_path))
        for group in ('metabolites', 'reactions', 'genes'):
            if {item.id for item in getattr(model, group)} != ids[group]:
                raise ValueError(f'JSON {group} changed while loading')
    except Exception as exc:
        raise GemBuildError(f'invalid_model: {exc}') from exc
    report = {'model_readable': True, 'cobra_version': cobra.__version__,
              'counts': {group: len(getattr(model, group)) for group in ids},
              'compartments': dict(model.compartments), 'objective_reaction_ids': objective,
              'kegg_annotation_counts': {group: sum(any(str(k).startswith('kegg') and bool(v)
                    for k, v in item.annotation.items()) for item in getattr(model, group)) for group in ids}}
    growth = {'medium': str(medium_path), 'status': 'not_validated', 'biomass_flux': None, 'passed': False}
    try:
        medium = json.loads(medium_path.read_text(encoding='utf-8'))
        if not isinstance(medium, dict) or any(not _finite(v) or v < 0 for v in medium.values()):
            raise ValueError('medium must map exchange IDs to finite nonnegative uptake rates')
        model.medium = medium
        model.solver.configuration.timeout = 30
        solution = model.optimize()
        flux = solution.objective_value
        growth.update(status=solution.status, biomass_flux=float(flux) if flux is not None and math.isfinite(flux) else None,
                      passed=solution.status == 'optimal' and flux is not None and math.isfinite(flux) and flux > 1e-8)
    except Exception as exc:
        growth.update(status='medium_or_solver_error', error=str(exc))
    report['growth'] = growth
    return report
