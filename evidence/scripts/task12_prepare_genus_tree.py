#!/usr/bin/env python3
"""Build explicit BAC120 genus-context panels from compatible R232 HMM coordinates."""
import argparse,csv,gzip,hashlib,json,re
from pathlib import Path
AA=set('ACDEFGHIKLMNPQRSTVWY')

def fasta(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        name=None;parts=[]
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name,''.join(parts)
                name=line[1:].split()[0];parts=[]
            else: parts.append(line.strip())
        if name is not None:yield name,''.join(parts)

def acc(s):return re.search(r'GC[AF]_\d+\.\d+',s).group()
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def seqsha(p):
    h=hashlib.sha256()
    # Compare nucleotide sequences ignoring FASTA header/line-wrap repairs.
    for _,seq in fasta(p):h.update(seq.upper().encode());h.update(b'\n')
    return h.hexdigest()
def project(seq,keep):return ''.join(seq[i] for i in keep)
def writefa(p,rows):p.write_text(''.join('>'+n+'\n'+s+'\n' for n,s in sorted(rows.items())))
def writetsv(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def occupancy(s):return sum(c in AA for c in s)/len(s)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();r=a.run
    db=Path('/scratch/opendata/genomics/GTDB/releases/release232')
    source=Path('results/full-cohort/2026-09-22/task11-opus-full/marker')
    search=Path('results/ani-cases/2026-10-01/task12-astra-gtdb-neighbours')
    tax=json.loads((r/'inputs/representative-taxonomy.json').read_text())
    q=list(csv.DictReader(open(search/'queries.tsv'),delimiter='\t'));qids={x['accession'] for x in q}
    aliases=json.loads((search/'validated-aliases.json').read_text());remove={x for v in aliases.values() for x in v}
    sources=[search/'ncbi-hits.jsonl',r/'inputs/ncbi-additional.jsonl',Path('results/ani-cases/2026-09-29/task12-astra-residual-five/inputs/ncbi-summary-initial.jsonl')]
    metadata={d['accession']:d for p in sources for d in map(json.loads,p.read_text().splitlines())}
    mask=(db/'masks/gtdb_r232_bac120.mask').read_text().strip();assert len(mask)==41084 and mask.count('1')==5010
    keep=[i for i,c in enumerate(mask) if c=='1']
    refs={};rawids={}
    for n,s in fasta(db/'msa/gtdb_r232_bac120.faa'):
        ac=acc(n)
        if ac in tax:
            assert ac not in refs and len(s)==41084
            refs[ac]=s;rawids[ac]=n
    assert refs.keys()==tax.keys()
    cohort=dict(fasta(source/'alignments/full.user_msa.fasta.gz'))
    canonical=dict(fasta(source/'alignments/canonical.user_msa.fasta.gz'))
    batch={n.removeprefix('full_'):Path(p) for p,n in csv.reader(open(source/'inputs/full.abs.batch.tsv'),delimiter='\t')}
    comparisons=[]
    for x in q:
        n=x['accession'];assert len(cohort[n])==41084
        assert project(cohort[n],keep)==canonical[n], 'Mask projection disagrees with saved canonical alignment'
        assert sha(x['staged_path'])==x['sha256']
        assert seqsha(x['staged_path'])==seqsha(batch[n]), 'Repaired marker input differs in nucleotide sequence'
        for alt in aliases[n]:
            if alt in refs:
                aa=cohort[n];bb=refs[alt]
                comparable=sum(u in AA and v in AA for u,v in zip(aa,bb))
                diffs=sum(u!=v for u,v in zip(aa,bb) if u in AA and v in AA)
                comparisons.append(dict(query=n,paired_reference=alt,compared_residues=comparable,residue_differences=diffs,query_occupancy=occupancy(aa),reference_occupancy=occupancy(bb)))
    allrows={k:v for k,v in refs.items() if k not in remove};allrows.update({n:cohort[n] for n in qids})
    # Primary sampling criterion predeclared before inference: >=50% canonical occupancy.
    # Rejected public reference rows remain in the all-candidate sensitivity alignment.
    eligible={n:s for n,s in allrows.items() if occupancy(project(s,keep))>=0.5}
    assert qids <= eligible.keys()
    filtered=[i for i in range(41084) if sum(s[i] in AA for s in eligible.values())/len(eligible)>=0.8]
    full_filtered={n:project(s,filtered) for n,s in eligible.items()}
    masked={n:project(s,keep) for n,s in eligible.items()}
    writefa(r/'inputs/candidates-full.faa',allrows)
    writefa(r/'inputs/candidates-masked.faa',{n:project(s,keep) for n,s in allrows.items()})
    writefa(r/'inputs/primary-full80.faa',full_filtered);writefa(r/'inputs/canonical-masked.faa',masked)
    rows=[]
    for n,s in sorted(allrows.items()):
        d=metadata[n];c=d.get('checkm_info',{});isquery=n in qids
        tx=tax.get(n) or next(tax[x] for x in aliases[n] if x in tax) if n in tax or any(x in tax for x in aliases.get(n,[])) else ''
        # Queries absent as representatives inherit explicit species-cluster taxonomy from previous joined results.
        if not tx:
            matches=[x for x in csv.DictReader(open(search/'all-hits.tsv'),delimiter='\t') if x['query']==n]
            tx=matches[0]['GTDB_R232_taxonomy'].rsplit(';',1)[0]+';'+matches[0]['query_R232_species']
        rows.append(dict(tip_id=n,role='focal_chromosome' if isquery else 'GTDB_R232_representative',selection='primary' if n in eligible else 'sensitivity_only_low_marker_occupancy',
          NCBI_name=d['organism']['organism_name'],strain=json.dumps(d['organism'].get('infraspecific_names',{}),sort_keys=True),
          GTDB_R232=tx,type_material=d.get('type_material',{}).get('type_label',''),
          assembly_level=d['assembly_info']['assembly_level'],completeness=c.get('completeness',''),contamination=c.get('contamination',''),
          full_marker_occupancy=occupancy(s),canonical_occupancy=occupancy(project(s,keep)),
          primary_occupancy=occupancy(full_filtered[n]) if n in eligible else '',source_row=n if isquery else rawids[n],
          source_alignment=str(source/'alignments/full.user_msa.fasta.gz') if isquery else str(db/'msa/gtdb_r232_bac120.faa'),sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
    writetsv(r/'inputs/panel.tsv',rows);writetsv(r/'tables/paired-marker-scope-check.tsv',comparisons)
    writetsv(r/'inputs/column-map.tsv',[dict(output_column=j,source_column=i+1) for j,i in enumerate(filtered,1)])
    hashpaths=[db/'msa/gtdb_r232_bac120.faa',db/'masks/gtdb_r232_bac120.mask',db/'taxonomy/gtdb_taxonomy.tsv',source/'alignments/full.user_msa.fasta.gz',source/'alignments/canonical.user_msa.fasta.gz',source/'inputs/full.abs.batch.tsv']+sources
    writetsv(r/'inputs/source-sha256.tsv',[dict(path=str(p),sha256=sha(p)) for p in hashpaths])
    summary=dict(candidate_tips=len(allrows),primary_tips=len(eligible),focal_tips=len(qids),full_columns=41084,primary_columns=len(filtered),canonical_columns=len(keep),held_references=sorted(allrows.keys()-eligible.keys()),all_query_nucleotide_inputs_match=True,all_query_canonical_projections_match=True)
    (r/'tables/preparation-validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
