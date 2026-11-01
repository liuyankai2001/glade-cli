"""Prepare writable model/GBFF/FASTA copies without changing source inputs."""
import copy
import json
from pathlib import Path


def valid_translation(feature):
    values = feature.qualifiers.get('translation', [])
    protein = values[0].rstrip('*') if values else ''
    return protein if protein and '*' not in protein and protein.isascii() and protein.isalpha() else ''


def clean_records(records):
    cleaned = copy.deepcopy(records)
    stats = dict(removed_pseudo_cds=0, kept_cds=0, alternative_cds_removed=0, invalid_translation_cds=0)
    for record in cleaned:
        pseudo = lambda f: any(k in f.qualifiers for k in ('pseudo','pseudogene'))
        parent = {tag for f in record.features if f.type == 'gene' and pseudo(f) for tag in f.qualifiers.get('locus_tag',[])}
        kept, selected = [], {}
        for feature in record.features:
            if feature.type != 'CDS':
                kept.append(feature)
                continue
            tags = feature.qualifiers.get('locus_tag',[])
            if pseudo(feature) or parent.intersection(tags):
                stats['removed_pseudo_cds'] += 1
                continue
            if len(tags) != 1 or not tags[0].isascii() or any(c.isspace() for c in tags[0]):
                raise ValueError('CDS requires one ASCII locus_tag: ' + str(tags))
            protein = valid_translation(feature)
            if not protein:
                stats['invalid_translation_cds'] += 1
                continue
            tag = tags[0]
            if tag in selected:
                stats['alternative_cds_removed'] += 1
                if (len(protein),len(feature)) <= (len(valid_translation(selected[tag])),len(selected[tag])):
                    continue
                kept.remove(selected[tag])
            selected[tag] = feature
            kept.append(feature)
        stats['kept_cds'] += len(selected)
        record.features = kept
    return cleaned, stats


def prepare(job_dir, parameters):
    import cobra
    from Bio import SeqIO
    from cobra.util.solver import linear_reaction_coefficients
    job_dir = Path(job_dir)
    work = job_dir / 'work'
    work.mkdir(exist_ok=True)
    model = cobra.io.load_json_model(str(job_dir / 'input/reference_model.json'))
    resolved = dict(parameters)
    biomass = parameters.get('biomass_reaction_id')
    if biomass is None:
        coefficients = linear_reaction_coefficients(model)
        if len(coefficients) != 1 or next(iter(coefficients.values())) <= 0 or model.objective.direction != 'max':
            raise ValueError('reference objective is ambiguous; provide biomass_reaction_id')
        biomass = next(iter(coefficients)).id
    if biomass not in model.reactions:
        raise ValueError('biomass_reaction_id not in reference model')
    model.objective = model.reactions.get_by_id(biomass)
    model.objective.direction = 'max'
    resolved.update(biomass_reaction_id=biomass, no_reannotation=True)
    (job_dir / 'resolved_parameters.json').write_text(json.dumps(resolved),encoding='utf-8')
    ref, ref_stats = clean_records(list(SeqIO.parse(str(job_dir / 'input/reference.gbff'),'genbank')))
    target, target_stats = clean_records(list(SeqIO.parse(str(job_dir / 'input/target.gbff'),'genbank')))
    if not ref or not target or not target_stats['kept_cds']:
        raise ValueError('GBFF inputs require annotated nonpseudo translated CDS')
    SeqIO.write(ref, str(work / 'reference.gbff'), 'genbank')
    SeqIO.write(target, str(work / 'target.gbff'), 'genbank')
    genes, proteins = {}, {}
    model_genes = {g.id for g in model.genes}
    for record in ref:
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            tag = feature.qualifiers['locus_tag'][0]
            if tag not in model_genes:
                continue
            if tag in genes:
                raise ValueError('duplicate reference locus_tag: ' + tag)
            genes[tag] = str(feature.extract(record.seq))
            proteins[tag] = valid_translation(feature)
    if not genes:
        raise ValueError('reference FASTA has no model gene matches')
    missing = sorted(model_genes - set(genes) - {'s0001','KPN_SPONT'})
    for name, sequences in [('genes.fna',genes),('proteins.faa',proteins)]:
        with (work / name).open('w',encoding='ascii',newline='\n') as stream:
            for tag, sequence in sorted(sequences.items()):
                stream.write('>' + tag + '\n' + sequence + '\n')
    mapping = {}
    if 's0001' in model.genes:
        if 'KPN_SPONT' in model.genes:
            raise ValueError('protected gene KPN_SPONT already exists')
        from cobra.manipulation.modify import rename_genes
        rename_genes(model, {'s0001':'KPN_SPONT'})
        mapping = {'s0001':'KPN_SPONT'}
    cobra.io.save_json_model(model, str(work / 'reference_model.json'))
    report = dict(reference_filter=ref_stats,target_filter=target_stats, reference_missing_biological_genes=missing,
                  missing_reference_affected_reactions=[dict(id=r.id,gpr=r.gene_reaction_rule) for r in model.reactions if {g.id for g in r.genes}.intersection(missing)],
                  spontaneous_label_mapping=mapping, reference_sequence_genes=len(genes), resolved_parameters=resolved)
    (job_dir / 'input_manifest.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return resolved
