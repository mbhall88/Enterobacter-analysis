#!/usr/bin/env python3
"""Join Task 12 evidence without selecting names or turning missing ANI into zero.

Input schemas are named in config/case-evidence.json. A fresh output directory is
required. Shared archive parsing/hash helpers come from the completed Task 12 audit.
"""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import statistics
from collections import Counter, defaultdict
from audit_ani_tables import sha256_file, accession_of, genome_key
from gtdb_context import parse_newick


def read_tsv(path):
    with open(path) as f: return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    rows=list(rows)
    if not fields:
        fields=list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def passes(ani, af_a, af_b, ani_min=95, af_min=80):
    return all(v is not None for v in [ani, af_a, af_b]) and ani>=ani_min and af_a>=af_min and af_b>=af_min


def load_pairs(path, allowed):
    """Canonicalize reciprocal rows, swapping AF with its genome, reject disagreement."""
    pairs={}; counts=Counter()
    with open(path) as f:
        reader=csv.DictReader(f,delimiter='\t')
        required={'Ref_file','Query_file','ANI','Align_fraction_ref','Align_fraction_query'}
        if not required <= set(reader.fieldnames): raise ValueError('Not a skani AF table')
        for r in reader:
            a,b=[accession_of(genome_key(r[k])) for k in ['Ref_file','Query_file']]
            if a not in allowed or b not in allowed: raise ValueError(f'Unexpected input: {a} {b}')
            vals=tuple(float(r[k]) for k in ['ANI','Align_fraction_ref','Align_fraction_query'])
            if not all(math.isfinite(v) and 0<=v<=100 for v in vals): raise ValueError('Invalid percent')
            if a==b: counts['self_rows']+=1; continue
            if a>b: a,b=b,a; vals=(vals[0],vals[2],vals[1])
            if (a,b) in pairs:
                if pairs[a,b]!=vals: raise ValueError(f'Conflicting reciprocal/duplicate {a} {b}')
                counts['duplicate_or_reciprocal_rows']+=1
            else: pairs[a,b]=vals
    counts['unique_measured_pairs']=len(pairs)
    return pairs,dict(counts)


def distribution(pairs, members):
    keys=list(itertools.combinations(sorted(members),2))
    measured=[pairs[k] for k in keys if k in pairs]
    ani=sorted(v[0] for v in measured)
    n=len(ani)
    return dict(genomes=len(members),expected_pairs=len(keys),measured_pairs=n,
                missing_pairs=len(keys)-n,ani_min=min(ani) if ani else '',
                ani_q25=ani[int((n-1)*.25)] if ani else '',
                ani_median=statistics.median(ani) if ani else '',
                ani_q75=ani[int((n-1)*.75)] if ani else '',ani_max=max(ani) if ani else '',
                passing_95_80=sum(passes(*v) for v in measured),
                ani_below_95=sum(v[0]<95 for v in measured),
                coverage_below_80=sum(min(v[1:])<80 for v in measured))


def tree_evidence(config, pairs, out):
    comparison=json.loads(Path(config['tree_comparison']).read_text())
    rooted=comparison['comparison_sets']['ingroup-1076']['rooting']['per_tree']
    membership=defaultdict(lambda:defaultdict(list)); records=[]; summary=[]; conflicts=[]
    for method in ['bac120_fasttree','verticall_fastme','mashtree_consensus']:
        item=rooted[method]
        tree=parse_newick(Path(item['rootings']['both_set12_outgroups']['tree']).read_text())
        tips={tree.label[n]:n for n in tree.tips()}
        policy=comparison['trees'][method]
        for group,g in item['type_anchor_topology']['named_groups'].items():
            if not g.get('resolvable'): continue
            node=tree.mrca([tips[a] for a in g['type_anchor_tips']])
            stack=[node]; members=[]
            while stack:
                n=stack.pop()
                if tree.children[n]: stack.extend(tree.children[n])
                else: members.append(tree.label[n])
            for a in members:
                if a in pairs.allowed: membership[a][method].append(group)
            for a in sorted(members): records.append(dict(method=method,group=group,accession=a))
            summary.append(dict(group=group,method=method,tree_scope='ingroup-1076',
                clade_tips=g['mrca_clade_tips'],exclusive_of_other_types=g['exclusive'],
                type_accessions=';'.join(g['type_anchor_tips']),
                support_convention=policy['support_convention']['kind'],
                support_rule=json.dumps(policy['support_convention'],sort_keys=True),
                near_zero_threshold=policy['near_zero_branch_length'],
                interpretation='raw MRCA membership; not a retained or supported clade claim',
                **distribution(pairs,{a for a in members if a in pairs.allowed})))
        for case,spec in config['tree_conflicts'].items():
            for a,b in itertools.combinations(sorted(spec['types']),2):
                key=a+'|'+b
                conflicts.append(dict(case=case,a=a,b=b,method=method,
                    intervening_type_anchors=item['type_anchor_topology']['intervening_type_anchors_per_pair'].get(key,''),
                    disposition=spec['status'], source=config['tree_comparison']))
    write_tsv(out/'tree-group-membership.tsv',records)
    write_tsv(out/'group-ani-distributions.tsv',summary)
    write_tsv(out/'task17-pending-conflicts.tsv',conflicts)
    # Preserve the original three-state split results and all method-specific rules verbatim.
    (out/'tree-evidence.json').write_text(json.dumps(dict(
        source=config['tree_comparison'],source_sha256=sha256_file(config['tree_comparison']),
        trees=comparison['trees'],
        primary=comparison['comparison_sets']['full-1089']['comparisons'],
        type_projection=comparison['comparison_sets']['type-anchor-projection-35']['comparisons'],
        terminology='present / unresolved / contradicted; support conventions remain separate',
        baseline_source=config['tree_baselines']),indent=2)+'\n')
    return membership


def merge_measurements(primary, supplement):
    """Fill gaps, rejecting a supplement that revises any existing measurement."""
    for key in primary.keys() & supplement.keys():
        if primary[key] != supplement[key]:
            raise ValueError(f'Supplement disagrees with existing pair: {key}')
    return {**primary, **supplement},set(supplement)-set(primary)


class PairTable(dict):
    pass


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--config',required=True)
    ap.add_argument('--output',required=True); args=ap.parse_args(argv)
    config=json.loads(Path(args.config).read_text()); config_hash=sha256_file(args.config)
    if not config.get('cohort'): raise ValueError('Named cohort required')
    if config['ani_rule'] != dict(ani_min=95,af_min=80,operator='>=',scope='chromosome_to_chromosome',unit='percent'):
        raise ValueError('This implementation requires the settled DEC-007 rule')
    out=Path(args.output); out.mkdir(parents=True,exist_ok=False)
    study={r['accession_version']:r for r in read_tsv(config['study'])}
    full={r['accession_version']:r for r in read_tsv(config['full'])}
    genomes={r['accession']:r for r in read_tsv(config['genomes'])}
    types=defaultdict(list)
    for r in read_tsv(config['types']):
        if r['reference_role']=='nomenclatural_type_anchor': types[r['selected_accession_version']].append(r['taxon_name'])
    data,counts=load_pairs(config['skani'],genomes)
    primary_keys=set(data)
    supplemental_sources={}
    for source in config.get('supplemental_skani',[]):
        extra,_=load_pairs(source,genomes)
        data,added=merge_measurements(data,extra)
        supplemental_sources.update({key:source for key in added})
    counts['supplemental_new_pairs']=len(supplemental_sources)
    counts['unique_measured_pairs']=len(data)
    pairs=PairTable(data); del data
    pairs.allowed=set(genomes)
    ledger=Path(config.get('pair_ledger',str(Path(config['run'])/'skani/pairs.tsv')))
    if ledger.exists(): raise ValueError('Pair ledger already exists; select a fresh run')
    missing=[]
    def expected():
        for a,b in itertools.combinations(sorted(genomes),2):
            v=pairs.get((a,b))
            if v is None: missing.append(dict(a=a,b=b,status='not_reported_by_skani'))
            yield dict(a=a,b=b,status='measured' if v else 'not_reported_by_skani',
                ani=v[0] if v else '',af_a=v[1] if v else '',af_b=v[2] if v else '',
                passes_95_80=passes(*v) if v else '',method='skani_0.3.1',
                scope='chromosome_to_chromosome',
                measurement_source=config['skani'] if (a,b) in primary_keys else supplemental_sources.get((a,b),''))
    write_tsv(ledger,expected())
    write_tsv(out/'missing-pairs.tsv',missing,['a','b','status'])
    counts.update(expected_pairs=len(genomes)*(len(genomes)-1)//2,missing_pairs=len(missing))
    membership=tree_evidence(config,pairs,out)
    current={}
    with open(config['current_ncbi']) as f:
        for line in f:
            r=json.loads(line)
            if r['accession'] in current: raise ValueError('Duplicate NCBI accession')
            current[r['accession']]=r
    cases=[]; typehits=[]; typepairs=[]
    for a,g in genomes.items():
        candidates=[]
        for b,names in types.items():
            if a==b: continue  # Type identity is a manifest fact, not independent self evidence.
            v=pairs.get(tuple(sorted((a,b))))
            if v:
                ani,aa,ab=v if a<b else (v[0],v[2],v[1])
                candidates.append((ani,min(aa,ab),b,aa,ab))
                typehits.append(dict(accession=a,type_accession=b,type_names=';'.join(names),
                    ani=ani,af_study=aa,af_type=ab,passes_95_80=passes(ani,aa,ab)))
            else: typehits.append(dict(accession=a,type_accession=b,type_names=';'.join(names),
                                      ani='',af_study='',af_type='',passes_95_80=''))
        candidates.sort(key=lambda x:(-x[0],-x[1],x[2]))
        hits=[c for c in candidates if passes(c[0],c[3],c[4])]
        best=candidates[0] if candidates else None
        s=study.get(a,{})
        n=current.get(a,{}); organism=n.get('organism',{}).get('organism_name','')
        flags=[]
        if a in config['quality_flags']: flags.append(config['quality_flags'][a])
        if any(w in s.get('legacy_revised_subspecies','').lower() for w in ['novel','exclude','?']): flags.append('legacy_uncertain_label')
        if full.get(a,{}).get('ingroup_status')=='declared_out_of_genus_context': flags.append('out_of_genus_context')
        if s.get('usage_restrictions'): flags.append('restricted_use')
        if a in types: status='verified_type_anchor'
        elif not hits: status='no_type_passes_95_80'
        elif len(hits)>1: status='multiple_type_anchors_pass_95_80'
        else: status='one_type_anchor_passes_95_80'
        row=dict(accession=a,cohort=config['cohort'] if a in study else 'additional_type_anchors',
            in_study=int(a in study),is_type=int(a in types),formal_type_names=';'.join(types.get(a,[])),
            historical_ncbi=s.get('historical_ncbi',''),historical_kraken_chromosome=s.get('historical_kraken_chromosome',''),
            historical_kraken_whole=s.get('historical_kraken_whole',''),
            legacy_revised_species=s.get('legacy_revised_species',''),legacy_revised_subspecies=s.get('legacy_revised_subspecies',''),
            current_ncbi=organism,current_ncbi_date=config['current_ncbi_date'] if n else '',
            current_ncbi_status=n.get('assembly_info',{}).get('assembly_status','not_returned'),
            current_ncbi_name_changed=(int(organism!=s['historical_ncbi']) if organism and s else ''),
            gtdb_r232_species=full.get(a,{}).get('gtdb_species',''),
            ingroup_status=g['ingroup_status'],quality_flags=';'.join(flags),restriction=g['restriction'],
            best_type_accession=best[2] if best else '',best_type_ani=best[0] if best else '',
            best_type_af_query=best[3] if best else '',best_type_af_reference=best[4] if best else '',
            passing_type_accessions=';'.join(c[2] for c in hits),n_passing_types=len(hits),
            passing_type_names=';'.join(sorted({name for c in hits for name in types[c[2]]})),
            type_pair_coverage=len(candidates),expected_type_pairs=len(types)-int(a in types),
            evidence_status=status,naming_outcome=config['naming_policy'],
            supported_biological_correction='not_adjudicated',config_sha256=config_hash)
        for m in ['bac120_fasttree','verticall_fastme','mashtree_consensus']:
            row[m+'_raw_group_membership']=';'.join(sorted(membership[a][m]))
        cases.append(row)
    for r in read_tsv(config['exclusions']):
        if r['cohort']=='public_discovery' and r['entity_id'] not in genomes:
            cases.append(dict(accession=r['entity_id'],cohort='Task01_public_exclusions',in_study=0,is_type=0,
                evidence_status='not_in_ANI_cohort',quality_flags=r['disposition_reasons'],
                restriction=r['usage_restrictions'],supported_biological_correction='not_adjudicated',config_sha256=config_hash))
    for a,b in itertools.combinations(sorted(types),2):
        v=pairs.get((a,b))
        typepairs.append(dict(a=a,b=b,type_names_a=';'.join(types[a]),type_names_b=';'.join(types[b]),
            ani=v[0] if v else '',af_a=v[1] if v else '',af_b=v[2] if v else '',passes_95_80=passes(*v) if v else ''))
    write_tsv(out/'cases.tsv',cases)
    write_tsv(out/'type-hits.tsv',typehits)
    write_tsv(out/'type-pairs.tsv',typepairs)
    namegroups=defaultdict(set)
    for a,s in study.items(): namegroups[s['legacy_revised_subspecies'] or 'legacy_label_missing'].add(a)
    write_tsv(out/'legacy-label-distributions.tsv',[
        dict(label=label,label_source='manual_legacy_revision_not_formal_taxonomy',**distribution(pairs,ids))
        for label,ids in sorted(namegroups.items())])
    (out/'summary.json').write_text(json.dumps(dict(cohort=config['cohort'],config_sha256=config_hash,
        pairs=counts,study=len(study),union=len(genomes),cases=len(cases),
        current_ncbi_returned=len(current),current_ncbi_missing=sorted(set(study)-current.keys()),
        evidence_status=dict(Counter(r['evidence_status'] for r in cases))),indent=2)+'\n')
    sourcekeys=['study','full','types','exclusions','crosswalk','genomes','skani','current_ncbi','tree_comparison','tree_baselines']
    (out/'inputs.sha256').write_text(''.join(f'{sha256_file(config[k])}  {config[k]}\n' for k in sourcekeys)+
        ''.join(f'{sha256_file(p)}  {p}\n' for p in config.get('supplemental_skani',[]))+
        f'{config_hash}  {args.config}\n')
    print(json.dumps(counts))


if __name__=='__main__': main()
