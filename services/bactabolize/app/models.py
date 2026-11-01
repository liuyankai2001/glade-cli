import math

TERMINAL = {'succeeded', 'failed', 'timed_out', 'cancelled'}
INPUTS = {'reference_model_file':'reference_model.json', 'reference_genome_file':'reference.gbff', 'target_genome_file':'target.gbff'}


def validate_parameters(data):
    defaults = dict(min_pident=80.0, min_coverage=25.0, media_type='m9', atmosphere_type='aerobic', biomass_reaction_id=None)
    if not isinstance(data, dict) or set(data) - set(defaults) - {'no_reannotation'}:
        raise ValueError('unknown parameters')
    if data.get('no_reannotation', True) is not True:
        raise ValueError('no_reannotation must be true')
    defaults.update({k:v for k,v in data.items() if k in defaults})
    for key in ('min_pident','min_coverage'):
        value = defaults[key]
        if isinstance(value, bool) or not isinstance(value, (int,float)) or not math.isfinite(value) or not 0 < value <= 100:
            raise ValueError(key + ' must be finite and in (0,100]')
    if defaults['atmosphere_type'] not in ('aerobic','anaerobic'):
        raise ValueError('invalid atmosphere_type')
    if defaults['media_type'] not in ('bg11','cdm_mendoza','lb_carveme','lb','lb_SEED','m9','m9_SEED','nutrient','pmm5_mendoza','pmm7_mendoza','tsa','tsa_sheep_blood'):
        raise ValueError('unsupported media_type')
    for key in ('media_type','biomass_reaction_id'):
        value = defaults[key]
        if key == 'biomass_reaction_id' and value is None:
            continue
        if not isinstance(value, str) or not value.strip() or len(value) > 200:
            raise ValueError('invalid ' + key)
    return defaults


def public_job(job):
    return {k:job.get(k) for k in ('job_id','status','parameters','created_at','started_at','finished_at','return_code','failure_code','error')}
