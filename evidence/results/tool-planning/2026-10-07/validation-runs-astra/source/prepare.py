"""Build execution/truth adapters using frozen metadata only, before outcomes exist."""
import hashlib
from collections import Counter

from common import ART, DESIGN, DOWNLOAD, EXPERT, FREEZE, PINS, ROOT, RUN, json_write, load_module, read, sha, verify_object, write


def coded(accession, layer):
    return 'VAL-' + layer.upper() + '-' + hashlib.sha256(
        ('Enterobacter-validation-execution-20261007-v1|' + accession).encode()).hexdigest()[:12].upper()


def truth_base(genome, layer, digest, block):
    return dict(genome_id=genome, layer=layer, truth_kind='', label_source='', truth_species_id='',
                truth_subspecies_id='', ani_species_passing='', subspecies_structural='',
                subspecies_label_source='', expected_species_ids='', expected_subspecies_id='',
                block=block, component_id=block, query_sha256=digest)


def main():
    verification = verify_object()
    if sha(DOWNLOAD) != '2bedcd18d50dc240658ce2bbaba78ed2a7cdc15da598f1eba95034ab42e8f5a3':
        raise ValueError('download manifest differs from posted seal')
    if sha(EXPERT) != '039c544311093bd4a6eb50862748598f365bdd30bc9602d321acf27db56c584b':
        raise ValueError('expert truth differs from posted seal')
    from enterobacter_classifier.reporting import ReportingContext
    context = ReportingContext(PINS['package'][0], expected_sha256=PINS['package'][1])
    concepts = context._concepts
    species = {c['display_name']: c['concept_id'] for c in concepts.values() if c['rank'] == 'species'}
    subspecies = {c['display_name']: c['concept_id'] for c in concepts.values() if c['rank'] == 'subspecies'}
    builder = load_module(DESIGN / 'source/build_design_artifacts.py', 'frozen_artifact_builder')
    expert = {r['genome_id']: r for r in read(EXPERT)}
    layer1, truth1, truth3, held3, manifest = read(ART / 'layer1-manifest.tsv'), [], [], [], []
    for r in layer1:
        named = r['stratum'] != 'context_or_unnamed_challenge'
        accession = r['accession']
        code = coded(accession, 'l1')
        manifest.append(dict(layer='layer1', genome_id=accession, coded_id=code,
                             fasta_path=str(ROOT / r['whole_fasta']), query_sha256=r['whole_sha256']))
        t = truth_base(accession, 'layer1', r['whole_sha256'], r['clone_block'])
        t.update(truth_kind='named' if named else 'context', label_source='project_label',
                 truth_species_id=species[r['project_species']] if named else '',
                 truth_subspecies_id=subspecies[r['project_species'] + ' subsp. ' + r['project_subspecies']]
                 if named and r['project_subspecies'] else '',
                 coded_id=code, **{k: r[k] for k in ['partition', 'stratum', 'assembly_level', 'label_status',
                     'saved_pairs_ani_class', 'reference_in_clone_block', 'reference_in_split_block', 'rule_design_role']})
        truth1.append(t)
    dispositions = {r['assembly_unit']: r for r in read(DESIGN / 'outputs/exposure-audit-attempt03/reserve-exposure-dispositions.tsv')}
    crosswalk = {r['ncbi_name_normalised']: r for r in read(ART / 'ncbi-name-crosswalk.tsv')}
    downloads = read(DOWNLOAD)
    for d in downloads:
        if d['status'] != 'downloaded':
            raise ValueError('unexpected unavailable assembly; preserve and review disposition')
        accession = d['accession']
        r = dispositions[accession]
        e = expert.get(accession)
        code = e['coded_id'] if e else coded(accession, 'l3')
        if e and e['query_sha256'] != d['sha256']:
            raise ValueError('expert query hash does not match download')
        manifest.append(dict(layer='layer3', genome_id=accession, coded_id=code,
                             fasta_path=str(ROOT / d['fasta_path']), query_sha256=d['sha256']))
        name = builder.normalise(r['ncbi_names'])
        c = crosswalk[name]
        t = truth_base(accession, 'layer3' if r['disposition'] == 'effective_reserve' else 'layer3_held',
                       d['sha256'], r['component_id'])
        t.update(truth_kind='crosswalk' if c['expected_species_ids'] else 'crosswalk_none', label_source='none',
                 expected_species_ids=c['expected_species_ids'], expected_subspecies_id=c['expected_subspecies_id'],
                 coded_id=code, ncbi_name_normalised=name, crosswalk_kind=c['kind'], qc_stratum=r['qc_stratum'],
                 assembly_level=d['assembly_level'], expert_sample=str(bool(e)), disposition=r['disposition'],
                 crosswalk_decision_note=c['open_decision'])
        (truth3 if r['disposition'] == 'effective_reserve' else held3).append(t)
    assert len(truth1) == 946 and sum(t['truth_kind'] == 'named' for t in truth1) == 922
    assert len(truth3) == 4821 and len(held3) == 36 and len(expert) == 100
    assert len({r['genome_id'] for r in manifest}) == len(manifest) == 5803
    assert len({r['coded_id'] for r in manifest}) == len(manifest)
    assert set(expert) <= {t['genome_id'] for t in truth3}
    for name, rows in [('execution-manifest', manifest), ('layer1-truth', truth1),
                       ('layer3-truth', truth3), ('layer3-held-truth', held3)]:
        write(RUN / 'inputs' / (name + '.tsv'), rows)
    paths = [DOWNLOAD, EXPERT, ART / 'layer1-manifest.tsv', ART / 'ncbi-name-crosswalk.tsv',
             ART / 'layer3-pool-strata.tsv', DESIGN / 'outputs/exposure-audit-attempt03/reserve-exposure-dispositions.tsv',
             DESIGN / 'validation-design.md', FREEZE / 'package-source.tsv', FREEZE / 'freeze-hashes.tsv',
             DESIGN / 'source/score_validation.py', DESIGN / 'source/draw_disagreement_review.py',
             DESIGN / 'source/build_design_artifacts.py', *[x[0] for x in PINS.values()]]
    paths += sorted((RUN / 'inputs').glob('*-truth.tsv')) + [RUN / 'inputs/execution-manifest.tsv']
    write(RUN / 'inputs/execution-pins.tsv', [dict(path=str(p), sha256=sha(p)) for p in paths])
    json_write(RUN / 'inputs/preparation.json', dict(verification=verification, genomes=len(manifest),
                layers=dict(Counter(r['layer'] for r in manifest)), layer1_named=922, layer3_effective=4821,
                layer3_held=36, expert=100, package_positive_policy=False,
                note='Metadata adapters only; frozen crosswalk wording is retained verbatim, including historical decision notes.'))
    print('Prepared 5803 execution rows and truth adapters; no classifier run.')


if __name__ == '__main__':
    main()
