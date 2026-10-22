"""Rerun only the AU test interrupted by Slurm; preserve its original outputs."""
from pathlib import Path
import json,os,re,shutil,subprocess,time,hashlib
r=Path(__file__).resolve().parents[1]
t=r/'local/GCA_026344075.1/tests';dest=t/'topology-masked-recovery';dest.mkdir(exist_ok=False)
expected=['RY-ml','scf-20260927','scf-20260928','118C5-original','118C5-alternate','gubbins','masked-ml','masked-scf','primary-boundary_with_nematophilus','primary-boundary_with_bugandensis','primary-boundary_outside_both','primary-AU','masked-boundary_with_nematophilus','masked-boundary_with_bugandensis','masked-boundary_outside_both']
for name in expected:
 rec=json.loads((t/'logs'/(name+'.command.json')).read_text());assert rec['exit_code']==0,name
for name in ['candidate-trees.nwk','candidate-order.json']:shutil.copy2(t/'topology-masked'/name,dest/name)
model=re.search(r'Best-fit model according to BIC: (\S+)',(t.parent/'tree/primary.iqtree').read_text()).group(1)
argv=['${USER_HOME}/sw/miniforge3/envs/snippy-ng/bin/iqtree3','-s',str(t/'nonrecombinant.core.fna'),'-z',str(dest/'candidate-trees.nwk'),'-n','0','-m',model,'-T','4','-zb','10000','-au','--seed','20260927','--keep-ident','--prefix',str(dest/'AU')]
record={'argv':argv,'cwd':str(dest),'job':os.getenv('SLURM_JOB_ID'),'reason':'28999982 reached its wall-time limit during final masked AU test; all preceding outputs reused','input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [t/'nonrecombinant.core.fna',dest/'candidate-trees.nwk']}}
(dest/'command-start.json').write_text(json.dumps(record,indent=2)+'\n')
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
b=time.monotonic()
with (dest/'stdout').open('w') as so,(dest/'stderr').open('w') as se:
 code=subprocess.run(['/usr/bin/time','-v','-o',str(dest/'resource.time')]+argv,cwd=dest,env=env,stdout=so,stderr=se).returncode
record.update(exit_code=code,seconds=time.monotonic()-b);(dest/'command.json').write_text(json.dumps(record,indent=2)+'\n');assert code==0
assert 'p-AU' in (dest/'AU.iqtree').read_text()
validation={'status':'PASS','primary_tips':38,'model':model,'scf_primary_seeds':[20260927,20260928],'scf_quartets_per_branch':100,'hypotheses':['boundary_with_nematophilus','boundary_with_bugandensis','boundary_outside_both'],'AU_resamples':10000,'masked_and_unmasked_tested':True,'118C5_individual_sensitivities':True,'recombination':'5 iterations maximum; inspect convergence log','species_rank_test':False,'recovered_final_AU':str(dest),'preceding_successful_commands':expected}
assert not (t/'validation.json').exists()
(t/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
print('PASS final AU recovery; collection/review still required')
