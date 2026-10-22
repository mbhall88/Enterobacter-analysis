"""Quantitative local-tree sensitivity and placement tests, never species-rank tests."""
from pathlib import Path
import argparse,csv,json,os,re,subprocess,time,hashlib
import numpy as np
from Bio import SeqIO,Phylo
r=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--reference',required=True);a=p.parse_args();o=r/'local'/a.reference
assert json.loads((o/'validation.json').read_text())['status']=='PASS'
rows=list(csv.DictReader((r/'inputs/panel-qc.tsv').open(),delimiter='\t'));rows=[x for x in rows if x['local_primary']=='True'];groups={}
for x in rows:
 g='nematophilus_candidate' if x['accession']=='GCA_026344075.1' else x['group'];groups.setdefault(g,[]).append(x['accession'])
assert {k:len(v) for k,v in groups.items()}=={'bugandensis':26,'nematophilus_candidate':6,'boundary':2,'context':4}
t=o/'tests';t.mkdir(exist_ok=False);(t/'logs').mkdir();iq='${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3';genv=r.parents[3]/'.work/envs/task12-localphylo' # overwritten by fixed project path
project=Path('${PROJECT_ROOT}');genv=project/'.work/envs/task12-localphylo'
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');env['PATH']=str(genv/'bin')+':'+env['PATH']
def run(name,argv):
 b=time.monotonic();rec={'argv':list(map(str,argv)),'cwd':str(t),'job':os.getenv('SLURM_JOB_ID')}
 with (t/'logs'/(name+'.stdout')).open('w') as so,(t/'logs'/(name+'.stderr')).open('w') as se:
  status=subprocess.run(['/usr/bin/time','-v','-o',str(t/'logs'/(name+'.time'))]+list(map(str,argv)),env=env,stdout=so,stderr=se,cwd=t).returncode
 rec.update(exit_code=status,seconds=time.monotonic()-b);(t/'logs'/(name+'.command.json')).write_text(json.dumps(rec,indent=2)+'\n');assert status==0,(name,status)
 print(name,round(rec['seconds'],1),flush=True)
model=re.search(r'Best-fit model according to BIC: (\S+)',(o/'tree/primary.iqtree').read_text()).group(1)
primary=o/'primary.core.fna';primarytree=o/'tree/primary.treefile';assert len(Phylo.read(primarytree,'newick').get_terminals())==38
# Added after the primary alignment composition screen flagged 12/38 taxa.
# R/Y recoding removes within-purine/pyrimidine transition signal; it is a sensitivity,
# not a guarantee against every form of compositional heterogeneity.
with (t/'RY.core.fna').open('w') as f:
 for rec in SeqIO.parse(primary,'fasta'):
  seq=str(rec.seq);assert set(seq)<=set('ACGT');f.write('>'+rec.id+'\n'+seq.translate(str.maketrans('AGCT','0011'))+'\n')
run('RY-ml',[iq,'-s',t/'RY.core.fna','-st','BIN','-m','GTR2+FO+G4','-T','4','-B','1000','--alrt','1000','--bnni','--seed','20260927','--keep-ident','--prefix',t/'RY'])

# Site concordance across two quartet seeds, on the fixed inferred topology.
for seed in [20260927,20260928]:run('scf-'+str(seed),[iq,'-s',primary,'-t',primarytree,'--tree-fix','-m',model,'--scfl','100','--seed',str(seed),'-T','4','--keep-ident','--prefix',t/('primary.scf-'+str(seed))])
# Add either flagged 118C5 assembly on the exact same core columns; they are never independent replicate evidence.
for label in ['118C5-original','118C5-alternate']:
 run(label,[iq,'-s',o/(label+'.core.fna'),'-m',model,'-T','4','--seed','20260927','--keep-ident','--prefix',t/label])
# Recombination is a sensitivity model on coordinate-preserving primary data, not a species test.
(t/'recomb').mkdir()
# Gubbins 3.4.3 accepts ONE outgroup name (verified in common.py), unlike
# the four-tip clade used for display rooting. Isolated cwd prevents its
# basename-derived intermediate files colliding between reference jobs.
assert 'GCA_003594915.1' in groups['context']
run('gubbins',[genv/'bin/run_gubbins.py','--verbose','--prefix',t/'recomb/gubbins','--threads','4','--tree-builder','fasttree','--model','GTRGAMMA','--model-fitter','iqtree','--starting-tree',primarytree,'--filter-percentage','100','--iterations','5','--seed','20260927','--invariant-site-correction','--outgroup','GCA_003594915.1',o/'primary.full-coordinate.fna'])
gff=t/'recomb/gubbins.recombination_predictions.gff';assert gff.exists()
core=np.load(o/'strict-core-mask.npy');mask=np.zeros(len(core),dtype=bool);events=0
for l in gff.open():
 if l.startswith('#') or not l.strip():continue
 f=l.rstrip().split('\t');start,end=int(f[3])-1,int(f[4]);assert 0<=start<end<=len(mask);mask[start:end]=True;events+=1
keep=core&~mask;assert keep.sum()>10000,int(keep.sum());np.save(t/'nonrecombinant-core-mask.npy',keep)
with (t/'nonrecombinant.core.fna').open('w') as f:
 for rec in SeqIO.parse(o/'primary.full-coordinate.fna','fasta'):
  seq=np.frombuffer(str(rec.seq).encode(),dtype='S1')[keep];assert np.isin(seq,np.array(list('ACGT'),dtype='S1')).all();f.write('>'+rec.id+'\n'+seq.tobytes().decode()+'\n')
(t/'recombination-summary.json').write_text(json.dumps({'events':events,'original_core':int(core.sum()),'removed_core_bases':int((core&mask).sum()),'retained_core_bases':int(keep.sum()),'mask':'union of all inferred recombinant intervals across all taxa; common sites retained for comparability'},indent=2)+'\n')
run('masked-ml',[iq,'-s',t/'nonrecombinant.core.fna','-m',model,'-T','4','-B','1000','--alrt','1000','--bnni','--seed','20260927','--keep-ident','--prefix',t/'masked'])
run('masked-scf',[iq,'-s',t/'nonrecombinant.core.fna','-t',t/'masked.treefile','--tree-fix','-m',model,'--scfl','100','--seed','20260927','-T','4','--keep-ident','--prefix',t/'masked.scf'])
# Three prespecified alternatives: boundary pair joins nematophilus, bugandensis, or is outside both.
clade=lambda k:'('+','.join(sorted(groups[k]))+')'
B,N,X,C=map(clade,['bugandensis','nematophilus_candidate','boundary','context'])
constraints={'boundary_with_nematophilus':f'({C},({B},({N},{X})));','boundary_with_bugandensis':f'({C},({N},({B},{X})));','boundary_outside_both':f'({C},({X},({B},{N})));'}
(t/'constraints').mkdir()
for name,s in constraints.items():(t/'constraints'/(name+'.nwk')).write_text(s+'\n')
for label,aln,free_tree in [('primary',primary,primarytree),('masked',t/'nonrecombinant.core.fna',t/'masked.treefile')]:
 folder=t/('topology-'+label);folder.mkdir();trees=[free_tree.read_text().strip()];ledger=[{'index':1,'hypothesis':'unconstrained','tree':str(free_tree)}]
 for name in constraints:
  prefix=folder/name
  run(label+'-'+name,[iq,'-s',aln,'-g',t/'constraints'/(name+'.nwk'),'-m',model,'-T','4','--seed','20260927','--keep-ident','--prefix',prefix])
  path=Path(str(prefix)+'.treefile');trees.append(path.read_text().strip());ledger.append({'index':len(trees),'hypothesis':name,'tree':str(path)})
 (folder/'candidate-trees.nwk').write_text('\n'.join(trees)+'\n');(folder/'candidate-order.json').write_text(json.dumps(ledger,indent=2)+'\n')
 run(label+'-AU',[iq,'-s',aln,'-z',folder/'candidate-trees.nwk','-n','0','-m',model,'-T','4','-zb','10000','-au','--seed','20260927','--keep-ident','--prefix',folder/'AU'])
(t/'validation.json').write_text(json.dumps({'status':'PASS','primary_tips':38,'model':model,'scf_primary_seeds':[20260927,20260928],'scf_quartets_per_branch':100,'hypotheses':list(constraints),'AU_resamples':10000,'masked_and_unmasked_tested':True,'118C5_individual_sensitivities':True,'recombination':'5 iterations maximum; inspect convergence log','species_rank_test':False},indent=2)+'\n')
