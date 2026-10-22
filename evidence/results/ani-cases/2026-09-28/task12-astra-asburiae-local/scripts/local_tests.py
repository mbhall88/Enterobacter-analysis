"""Support, concordance and prespecified topology comparisons; not species-rank tests."""
from pathlib import Path
import argparse,csv,json,os,re,subprocess,time,hashlib
from Bio import SeqIO,Phylo
R=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser();p.add_argument('--reference',required=True);a=p.parse_args();o=R/'local'/a.reference
assert json.loads((o/'validation.json').read_text())['status']=='PASS'
assert json.loads((R/'outputs/test-design-validation.json').read_text())['status']=='PASS'
T=o/'tests';T.mkdir(exist_ok=False);(T/'logs').mkdir();iq='${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3'
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');assert int(env.get('SLURM_CPUS_PER_TASK','1'))>=4
names={x['accession'] for x in csv.DictReader((R/'inputs/panel.tsv').open(),delimiter='\t')};assert len(names)==93

def run(name,argv):
 argv=list(map(str,argv));b=time.monotonic();rec={'argv':argv,'cwd':str(T),'job':os.getenv('SLURM_JOB_ID')}
 with (T/'logs'/(name+'.stdout')).open('w') as so,(T/'logs'/(name+'.stderr')).open('w') as se:
  status=subprocess.run(['/usr/bin/time','-v','-o',str(T/'logs'/(name+'.time'))]+argv,env=env,stdout=so,stderr=se,cwd=T).returncode
 rec.update(exit_code=status,seconds=time.monotonic()-b);(T/'logs'/(name+'.command.json')).write_text(json.dumps(rec,indent=2)+'\n');assert status==0,(name,status)
 print(name,round(rec['seconds'],1),flush=True)

def signature(path):
 t=Phylo.read(path,'newick');assert {x.name for x in t.get_terminals()}==names
 result=set()
 for n in t.get_nonterminals():
  ss={x.name for x in n.get_terminals()};co=names-ss
  if min(len(ss),len(co))<2:continue
  result.add(min(tuple(sorted(ss)),tuple(sorted(co))))
 return frozenset(result)

primary=o/'primary.core.fna';free=o/'tree/primary.treefile';report=(o/'tree/primary.iqtree').read_text();model=re.search(r'Best-fit model according to BIC: (\S+)',report).group(1);signature(free)
for seed in [20260928,20260929]:
 run('scf-'+str(seed),[iq,'-s',primary,'-t',free,'--tree-fix','-m',model,'--scfl','100','--seed',str(seed),'-T','4','--keep-ident','--prefix',T/('primary.scf-'+str(seed))])
# All constraints have the same 93 leaves, with only the specified union constrained.
folder=T/'topology';folder.mkdir();records=[dict(hypothesis='unconstrained',tree=str(free))]
for h in json.loads((R/'inputs/hypotheses.json').read_text()):
 name=h['hypothesis'];prefix=folder/name
 assert hashlib.sha256(Path(h['constraint']).read_bytes()).hexdigest()==h['sha256']
 run(name,[iq,'-s',primary,'-g',h['constraint'],'-m',model,'-T','4','--seed','20260928','--keep-ident','--prefix',prefix])
 tree=Path(str(prefix)+'.treefile');sig=signature(tree);ingroup=set(h['clade_members']);canon=min(tuple(sorted(ingroup)),tuple(sorted(names-ingroup)));assert canon in sig,(name,'constraint not met')
 records.append(dict(hypothesis=name,tree=str(tree)))
seen={};unique=[]
for rec in records:
 sig=signature(rec['tree'])
 if sig not in seen:seen[sig]=len(unique)+1;unique.append(Path(rec['tree']).read_text().strip())
 rec['test_index']=seen[sig]
(folder/'candidate-trees.nwk').write_text('\n'.join(unique)+'\n');(folder/'candidate-order.json').write_text(json.dumps(records,indent=2)+'\n')
assert len(unique)>=2
run('AU',[iq,'-s',primary,'-z',folder/'candidate-trees.nwk','-n','0','-m',model,'-T','4','-zb','10000','-au','--seed','20260928','--keep-ident','--prefix',folder/'AU'])
# Conditional prespecified composition sensitivity; no Gubbins or dDDH.
failed=[]
for line in report.splitlines():
 if 'failed' in line:
  match=re.search(r'(GCA_[0-9]+\.[0-9]+)',line)
  if match:failed.append(match.group(1))
if failed:
 with (T/'RY.core.fna').open('w') as f:
  for rec in SeqIO.parse(primary,'fasta'):
   seq=str(rec.seq);assert set(seq)<=set('ACGT');f.write('>'+rec.id+'\n'+seq.translate(str.maketrans('AGCT','0011'))+'\n')
 run('RY',[iq,'-s',T/'RY.core.fna','-st','BIN','-m','GTR2+FO+G4','-T','4','-B','1000','--alrt','1000','--bnni','--seed','20260928','--keep-ident','--prefix',T/'RY']);signature(T/'RY.treefile')
(T/'validation.json').write_text(json.dumps(dict(status='PASS',tips=93,reference=a.reference,model=model,scf_seeds=[20260928,20260929],scf_quartets=100,AU_resamples=10000,candidates=len(records),unique_topologies=len(unique),composition_failures=failed,RY_run=bool(failed),species_rank_test=False,recombination_masking=False),indent=2)+'\n')
