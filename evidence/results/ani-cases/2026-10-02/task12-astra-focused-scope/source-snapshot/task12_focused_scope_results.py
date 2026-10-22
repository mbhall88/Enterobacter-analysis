#!/usr/bin/env python3
"""Recompute directional ANIm, conserved inputs, record context and marker evidence."""
import json,re
from collections import Counter
from Bio import SeqIO
from pyani.anim import parse_delta
import task12_focused_scope as s
R=s.RUN

def anim_results():
    raw=[];summary={}
    for j in s.rows(R/'anim-jobs.tsv'):
        d=__import__('pathlib').Path(j['reuse_source']) if j['reuse_source'] else R/'anim'/j['mode']/j['pair']
        cmd=json.loads((d/'command.json').read_text());assert cmd['exit_code']==0
        lengths={j['target']:sum(len(v[1]) for v in s.fasta(j['target_path']).values()),j['comparator']:sum(len(v[1]) for v in s.fasta(j['comparator_path']).values())}
        for side in ['target','comparator']:assert s.sha(j[side+'_path'])==j[side+'_sha256']
        fs=sorted((d/'work/nucmer_output').glob('*/*.filter'));assert len(fs)==2
        import csv
        with open(d/'work/ANIm_percentage_identity.tab') as f:
            m=list(csv.reader(f,delimiter='\t'));matrix={r[0]:dict(zip(m[0][1:],map(float,r[1:]))) for r in m[1:]}
        rr=[]
        for p in fs:
            ref,q=p.stem.split('_vs_');lr,lq,ani,errors=parse_delta(p);assert abs(ani-matrix[ref][q])<1e-10
            r=dict(pair=j['pair'],mode=j['mode'],reference=ref,query=q,ani_pct=ani*100,reference_aligned_bp=lr,query_aligned_bp=lq,reference_bp=lengths[ref],query_bp=lengths[q],reference_af=lr/lengths[ref]*100,query_af=lq/lengths[q]*100,errors=errors,raw_path=str(p),sha256=s.sha(p),reused=bool(j['reuse_source']))
            raw.append(r);rr.append(r)
        a=j['target'];ani=min(x['ani_pct'] for x in rr);af=min(x['reference_af'] if x['reference']==a else x['query_af'] for x in rr);bf=min(x['query_af'] if x['reference']==a else x['reference_af'] for x in rr)
        summary[(j['pair'],j['mode'])]=dict(ani=ani,target_af=af,comparator_af=bf,passes=ani>=95 and min(af,bf)>=80)
    s.write(R/'anim-raw-directions.tsv',raw)
    out=[]
    for pair in sorted({j['pair'] for j in s.rows(R/'anim-jobs.tsv')}):
        o,p=summary[(pair,'original')],summary[(pair,'plus')];r={'pair':pair}
        for k in o:r['original_'+k]=o[k];r['plus_'+k]=p[k]
        for k in ['ani','target_af','comparator_af']:r[k+'_change_pp']=p[k]-o[k]
        r['pass_changed']=o['passes']!=p['passes'];out.append(r)
    s.write(R/'anim-sensitivity.tsv',out)
    (R/'anim-validation.json').write_text(json.dumps(dict(pairs=len(out),raw_directions=len(raw),identity_matrices_verified=True,input_hashes_unchanged=True,pass_changes=sum(r['pass_changed'] for r in out)),indent=2))

def annotations():
    summary=[];features=[];blocks=[]
    for j in s.rows(R/'inputs.tsv'):
        a,rid=j['accession'],j['flagged_record'];rec=SeqIO.read(R/'inputs'/f'{rid}.gb','genbank')
        expected=s.fasta(j['whole_path'])[rid][1];assert str(rec.seq).upper()==expected and rec.id==rid
        counters=Counter()
        for f in rec.features:
            counters[f.type]+=1
            q=f.qualifiers;product=';'.join(q.get('product',[]));gene=';'.join(q.get('gene',[]))
            if f.type not in ['CDS','rRNA','tRNA']:continue
            mobile=bool(re.search(r'integras|transpos|conjug|mobiliz|relaxas|phage|plasmid|partition|type IV secretion',product,re.I))
            if mobile:counters['mobile_keyword_CDS']+=1
            features.append(dict(accession=a,record=rid,type=f.type,start=int(f.location.start)+1,end=int(f.location.end),strand=f.location.strand,location=str(f.location),locus_tag=';'.join(q.get('locus_tag',[])),gene=gene,product=product,mobile_keyword=mobile,pseudo='pseudo' in q or 'pseudogene' in q))
        source=next(f for f in rec.features if f.type=='source')
        matches={}
        for mode in ['neighbour','self_chromosome']:
            path=s.old.RUN/'alignments'/f'S1_{a}'/f'{mode}.coords'
            for l in path.read_text().splitlines():
                x=l.split('\t');assert len(x)==11
                if x[10]!=rid:continue
                blocks.append(dict(accession=a,mode=mode,target_record=x[9],query_record=x[10],target_start=int(x[0]),target_end=int(x[1]),query_start=int(x[2]),query_end=int(x[3]),target_block_bp=int(x[4]),query_block_bp=int(x[5]),identity_pct=float(x[6])))
            bs=[x for x in blocks if x['accession']==a and x['mode']==mode];iv=[(b['query_start'],b['query_end']) for b in bs]
            matches[mode+'_union_bp']=s.old.bp(iv);matches[mode+'_largest_block_bp']=max(b['query_block_bp'] for b in bs)
        summary.append(dict(accession=a,record=rid,bp=len(rec),gc_pct=100*(expected.count('G')+expected.count('C'))/len(expected),topology=rec.annotations.get('topology'),deposited_plasmid=';'.join(source.qualifiers.get('plasmid',[])),CDS=counters['CDS'],rRNA=counters['rRNA'],tRNA=counters['tRNA'],mobile_keyword_CDS=counters['mobile_keyword_CDS'],**matches,annotation_sequence_sha256=s.sha(R/'inputs'/f'{rid}.gb')))
    s.write(R/'record-context.tsv',summary);s.write(R/'record-features.tsv',features);s.write(R/'alignment-blocks.tsv',blocks)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':
    import sys
    {'anim':anim_results,'annotations':annotations}[sys.argv[1]]()
