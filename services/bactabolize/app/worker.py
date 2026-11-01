"""Everything scientific runs inside this isolated subprocess's process group."""
import json
import math
import subprocess
import sys
from pathlib import Path

def validate_model(job_dir, parameters, native_return_code):
    import cobra
    from bactabolize import package_data
    path = job_dir / 'raw/draft_model.json'
    model = cobra.io.load_json_model(str(path))
    # Ensure both formal native outputs were fully written and are readable.
    cobra.io.read_sbml_model(str(job_dir / 'raw/draft_model.xml'))
    if not model.reactions or not model.metabolites:
        raise ValueError('generated model is empty')
    growth = dict(medium=parameters['media_type'],atmosphere=parameters['atmosphere_type'],status='error',biomass_flux=None,passed=False)
    try:
        with model:
            for reaction in model.exchanges:
                reaction.lower_bound = 0
            medium = package_data.get_data('media_definitions', parameters['media_type'])
            missing = []
            for reaction_id, lower_bound in medium['exchanges'].items():
                if reaction_id in model.reactions:
                    model.reactions.get_by_id(reaction_id).lower_bound = lower_bound
                else:
                    missing.append(reaction_id)
            if 'EX_o2_e' in model.reactions:
                model.reactions.get_by_id('EX_o2_e').lower_bound = -20 if parameters['atmosphere_type'] == 'aerobic' else 0
            # The prepared reference already has the selected biomass objective.
            # Reassigning it is unnecessary and triggers stale GLPK column
            # indices with legacy optlang after native model reconstruction.
            from cobra.util.solver import linear_reaction_coefficients
            coefficients = linear_reaction_coefficients(model)
            if len(coefficients) != 1 or next(iter(coefficients)).id != parameters['biomass_reaction_id'] or next(iter(coefficients.values())) <= 0:
                raise ValueError('draft biomass objective differs from resolved parameters')
            solution = model.optimize()
            flux = solution.objective_value
            flux = float(flux) if flux is not None and math.isfinite(flux) else None
            growth.update(status=solution.status,biomass_flux=flux,passed=solution.status == 'optimal' and flux is not None and flux >= 1e-4, missing_exchanges=missing)
    except Exception as exc:
        growth['error'] = str(exc)
    report = dict(model_readable=True,growth=growth,native_return_code=native_return_code, gapfill_applied=False)
    (job_dir / 'service_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')


def main():
    if sys.argv[1] in ('prepare','validate'):
        stage = sys.argv[1]
        job_dir = Path(sys.argv[2]).resolve()
        if stage == 'prepare':
            from .preparation import prepare
            parameters = json.loads((job_dir/'parameters.json').read_text(encoding='utf-8'))
            prepare(job_dir,parameters)
        else:
            resolved = json.loads((job_dir/'resolved_parameters.json').read_text(encoding='utf-8'))
            validate_model(job_dir,resolved,int(sys.argv[3]))
        return
    job_dir = Path(sys.argv[1]).resolve()
    stage_command = [sys.executable,'-B','-m',__package__ + '.worker']
    # No scientific objects live in this orchestrator. Each stage receives a
    # fresh interpreter; all children inherit the outer worker's process group.
    preparation = subprocess.run(stage_command + ['prepare',str(job_dir)])
    if preparation.returncode:
        raise SystemExit(preparation.returncode)
    resolved = json.loads((job_dir/'resolved_parameters.json').read_text(encoding='utf-8'))
    raw = job_dir / 'raw'
    raw.mkdir(exist_ok=True)
    work = job_dir / 'work'
    args = [sys.executable,'-B','-m','bactabolize','draft_model','--assembly_fp',str(work/'target.gbff'), '--ref_genes_fp',str(work/'genes.fna'), '--ref_proteins_fp',str(work/'proteins.faa'), '--ref_model_fp',str(work/'reference_model.json'), '--output_fp',str(raw/'draft'), '--no_reannotation']
    for key in ('min_pident','min_coverage','media_type','atmosphere_type','biomass_reaction_id'):
        args.extend(['--' + key,str(resolved[key])])
    native_code = subprocess.run(args).returncode
    if native_code not in (0,101):
        raise SystemExit(native_code)
    # Native 101 denotes growth failure after exporting the draft. Validate only
    # once the native process has completely exited and closed all output files.
    validation = subprocess.run(stage_command + ['validate',str(job_dir),str(native_code)])
    if validation.returncode:
        raise SystemExit(validation.returncode)


if __name__ == '__main__':
    main()
