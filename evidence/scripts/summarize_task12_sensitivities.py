#!/usr/bin/env python3
"""Summarize ANIm and Galah independently against the direct Task 12 pair ledger."""
import argparse
import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path
from audit_ani_tables import sha256_file, accession_of, genome_key
from build_case_register import read_tsv, write_tsv, passes


def matrix(path):
    with open(path) as f:
        return {r['']: {k:float(v) for k,v in r.items() if k} for r in csv.DictReader(f,delimiter='\t')}


def summarize_anim(run, pairs, output):
    selected=read_tsv(run/'anim/selected-pairs.tsv')
    result=[]; hashes=[]
    for r in selected:
        a,b=r['a'],r['b']; dest=run/'anim/pairs'/(a+'__'+b)
        record=json.loads((dest/'command.json').read_text())
        if record.get('exit_code')!=0: raise ValueError(f'Incomplete ANIm: {a} {b}')
        matrices={}
        for name in ['percentage_identity','alignment_coverage','alignment_lengths','similarity_errors']:
            path=dest/'work'/('ANIm_'+name+'.tab')
            matrices[name]=matrix(path); hashes.append((sha256_file(path),path))
        m=matrices
        ani_ab=m['percentage_identity'][a][b]*100
        ani_ba=m['percentage_identity'][b][a]*100
        af_a=m['alignment_coverage'][a][b]*100
        af_b=m['alignment_coverage'][b][a]*100
        s=pairs.get((a,b))
        row=dict(r,anim_identity_a_to_b=ani_ab,anim_identity_b_to_a=ani_ba,
            anim_coverage_a=af_a,anim_coverage_b=af_b,
            anim_aligned_bases_a=m['alignment_lengths'][a][b],
            anim_aligned_bases_b=m['alignment_lengths'][b][a],
            anim_errors_a_to_b=m['similarity_errors'][a][b],anim_errors_b_to_a=m['similarity_errors'][b][a],
            anim_passes_conservative_95_80=passes(min(ani_ab,ani_ba),af_a,af_b),
            skani_ani=s[0] if s else '',skani_af_a=s[1] if s else '',skani_af_b=s[2] if s else '',
            skani_passes_95_80=passes(*s) if s else '',
            identity_min_delta_anim_minus_skani=min(ani_ab,ani_ba)-s[0] if s else '',
            implementation='legacy_pyani_package_0.2.13.1_cli_0.3.0-alpha_ANIm_MUMmer_package_3.23',
            rule='ANIm identities retained both ways; conservative sensitivity uses minimum identity and both coverage values; separate estimator, not a tiebreaker')
        result.append(row)
    write_tsv(output/'anim-vs-skani.tsv',result)
    (output/'anim-matrix-hashes.sha256').write_text(''.join(f'{h}  {p}\n' for h,p in hashes))
    return dict(pairs=len(result),max_identity_direction_difference=max(abs(r['anim_identity_a_to_b']-r['anim_identity_b_to_a']) for r in result),
        rule_outcome_changes=sum(r['anim_passes_conservative_95_80']!=r['skani_passes_95_80'] for r in result),
        min_identity_delta=min(r['identity_min_delta_anim_minus_skani'] for r in result),
        max_identity_delta=max(r['identity_min_delta_anim_minus_skani'] for r in result))


def cluster_map(path):
    assignments={}
    with open(path) as f:
        for line in f:
            rep,member=line.rstrip().split('\t')
            rep,member=[accession_of(genome_key(x)) for x in [rep,member]]
            if not rep or not member: raise ValueError(f'Bad cluster row {line}')
            # Galah emits the reference self-row twice if --reference-genomes
            # also occurs in the input list. This is one genome, not two.
            # Never collapse conflicting assignments or duplicate member rows.
            if member in assignments:
                if rep == member == assignments[member]:
                    continue
                raise ValueError(f'Conflicting/duplicate cluster row {line}')
            assignments[member]=rep
    if not set(assignments.values())<=set(assignments): raise ValueError('Representative missing self row')
    return assignments


def audit_clusters(assignments, pairs, threshold):
    groups=defaultdict(list)
    for a,rep in assignments.items(): groups[rep].append(a)
    summary=dict(genomes=len(assignments),clusters=len(groups),within_expected=0,within_measured=0,
        within_missing=0,within_below_ani=0,within_below_reciprocal_af=0,within_fails_ani_or_af=0,
        representative_pairs=0,representative_missing=0,representative_fails_ani_or_af=0,
        eligible_cross_cluster_pairs=0,borderline_eligible_pairs=0,borderline_cross_cluster_pairs=0)
    repfail=[]
    for rep,members in groups.items():
        for a,b in itertools.combinations(sorted(members),2):
            summary['within_expected']+=1
            v=pairs.get((a,b))
            if not v: summary['within_missing']+=1; continue
            summary['within_measured']+=1
            summary['within_below_ani']+=v[0]<threshold
            summary['within_below_reciprocal_af']+=min(v[1:])<80
            summary['within_fails_ani_or_af']+=not passes(*v,ani_min=threshold)
        for a in members:
            if a==rep: continue
            summary['representative_pairs']+=1
            v=pairs.get(tuple(sorted((a,rep))))
            if not v: summary['representative_missing']+=1
            elif not passes(*v,ani_min=threshold):
                summary['representative_fails_ani_or_af']+=1
                repfail.append(dict(member=a,representative=rep,ani=v[0],af_min=min(v[1:]),af_max=max(v[1:])))
    for (a,b),v in pairs.items():
        if a not in assignments or b not in assignments or not passes(*v,ani_min=threshold): continue
        cross=assignments[a]!=assignments[b]
        summary['eligible_cross_cluster_pairs']+=cross
        if v[0]<=threshold+0.1:
            summary['borderline_eligible_pairs']+=1
            summary['borderline_cross_cluster_pairs']+=cross
    return summary,repfail


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--run',required=True)
    ap.add_argument('--output',required=True); ap.add_argument('--allow-incomplete-galah',action='store_true')
    ap.add_argument('--pairs',type=Path,help='Explicit pair ledger; defaults to the original run ledger')
    args=ap.parse_args(argv); run=Path(args.run); out=Path(args.output); out.mkdir(parents=True,exist_ok=False)
    pairs={}
    ledger=args.pairs or run/'skani/pairs.tsv'
    for r in read_tsv(ledger):
        if r['status']=='measured': pairs[r['a'],r['b']]=tuple(float(r[k]) for k in ['ani','af_a','af_b'])
    anim=summarize_anim(run,pairs,out)
    results=[]; failures=[]; maps={}; incomplete=[]
    for path in sorted((run/'galah').glob('ani*')):
        rec=json.loads((path/'command.json').read_text())
        if rec.get('exit_code')!=0: incomplete.append(path.name); continue
        assignments=cluster_map(path/'clusters.tsv')
        cutoff=float(path.name.split('-')[0][3:])
        stats,flagged=audit_clusters(assignments,pairs,cutoff)
        results.append(dict(run=path.name,threshold=cutoff,**stats,
            command_sha256=sha256_file(path/'command.json'),clusters_sha256=sha256_file(path/'clusters.tsv'),
            rule='Galah skani backend precluster 90; requested AF80 is one-sided-max; audit here requires both AF>=80'))
        failures.extend(dict(run=path.name,**r) for r in flagged)
        maps[path.name]=assignments
    expected=[f'ani{x:g}-forward' for x in [94.5,95,95.5,96,97.5,98,98.5]]
    if not args.allow_incomplete_galah and (incomplete or any(x not in maps for x in expected)):
        raise ValueError('Galah sweep incomplete')
    write_tsv(out/'galah-audit.tsv',results)
    write_tsv(out/'galah-representative-rule-failures.tsv',failures,['run','member','representative','ani','af_min','af_max'])
    dependence=[]
    for name,mapping in maps.items():
        if name.endswith('-forward') and '-type-' not in name: continue
        prefix=name.split('-')[0]
        baseline=maps.get(prefix+'-type-forward')
        if not baseline: continue
        if set(baseline)!=set(mapping): raise ValueError('Order sensitivity input mismatch')
        requested=name.split('-force-',1)[1] if '-force-' in name else ''
        dependence.append(dict(run=name,baseline=prefix+'-type-forward',genomes=len(mapping),
            requested_reference=requested,
            requested_reference_is_representative=mapping.get(requested)==requested if requested else '',
            representatives_changed=sum(mapping[a]!=baseline[a] for a in mapping),
            co_membership_pairs_changed=sum((mapping[a]==mapping[b])!=(baseline[a]==baseline[b]) for a,b in itertools.combinations(sorted(mapping),2)),
            scope='33 type genomes, diagnostic representatives only; not an operational reference choice'))
    write_tsv(out/'galah-order-representative-dependence.tsv',dependence)
    (out/'summary.json').write_text(json.dumps(dict(pair_ledger=str(ledger),pair_ledger_sha256=sha256_file(ledger),anim=anim,galah_completed=len(results),
        galah_incomplete=incomplete,galah_missing_primary=[x for x in expected if x not in maps]),indent=2)+'\n')


if __name__=='__main__': main()
