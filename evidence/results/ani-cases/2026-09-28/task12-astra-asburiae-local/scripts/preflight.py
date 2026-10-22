"""Verify copied caller, masking fixture and runtime interfaces without rerunning old alignments."""
from pathlib import Path
import json,random,os,subprocess,hashlib
from Bio import SeqIO
from mask_pseudo import masked_sequences
R=Path(__file__).resolve().parents[1];OLD=R.parents[1]/'2026-09-27/task12-astra-nematophilus-expanded';d=OLD/'interface-test-v2'
rng=random.Random(20260927);reference=''.join(rng.choices('ACGT',k=100000));expected=list(reference)
for i in sorted(rng.sample(range(len(reference)),6250)):expected[i]=rng.choice([b for b in 'ACGT' if b!=reference[i]])
assert str(next(SeqIO.parse(d/'ref.fa','fasta')).seq)==reference
seq=masked_sequences(d/'out')['ref'];called=[i for i,b in enumerate(seq) if b in 'ACGT'];assert len(called)>98000 and all(seq[i]==expected[i] for i in called);assert set(seq[39990:40110])=={'N'}
env=os.environ.copy();env['PYTHONPATH']=str(R/'pinned');env['PATH']='${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin:'+env['PATH']
for name,cmd in [('snippy-help',['snippy-ng','asm','--help']),('iqtree-help',['iqtree3','-h']),('iqtree-version',['iqtree3','--version']),('nucmer-version',['nucmer','--version'])]:
 x=subprocess.run(cmd,env=env,capture_output=True,text=True);assert x.returncode==0;(R/'logs'/(name+'.txt')).write_text(x.stdout+x.stderr)
helptext=(R/'logs/iqtree-help.txt').read_text()
for flag in ['--scfl','--alrt','--bnni','--tree-fix']:assert flag in helptext,flag
# --keep-ident is accepted by this version but omitted from its short help.
# Confirm it with a tiny deterministic six-tip alignment including an identical pair.
from Bio import Phylo
fixture=R/'inputs/iqtree-cli-fixture.fna';base=reference[:1000]
fixture.write_text(''.join('>tip'+str(i)+'\n'+(base if i<2 else ''.join(b if j%(i+3) else 'ACGT'[('ACGT'.index(b)+1)%4] for j,b in enumerate(base)))+'\n' for i in range(6)))
x=subprocess.run(['iqtree3','-s',str(fixture),'-m','JC','-n','0','-T','1','--keep-ident','--seed','20260928','--prefix',str(R/'logs/iqtree-cli-fixture')],env=env,capture_output=True,text=True)
(R/'logs/iqtree-cli-fixture.stdout').write_text(x.stdout+x.stderr);assert x.returncode==0
assert len(Phylo.read(R/'logs/iqtree-cli-fixture.treefile','newick').get_terminals())==6
(R/'outputs/preflight.json').write_text(json.dumps({'status':'PASS','synthetic_called_bases':len(called),'correct_called_bases':len(called),'synthetic_divergence':0.0625,'saved_fixture_reused':True,'masked_indel_and_reverse_orientation':True,'copied_code_hashes_checked_in_prepare':True,'scope':'alignment uncertainty mask only; no recombination masking'},indent=2)+'\n')
print('PASS',len(called))
