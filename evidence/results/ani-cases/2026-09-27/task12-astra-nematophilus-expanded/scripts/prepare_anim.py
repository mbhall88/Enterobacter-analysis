import csv,json,subprocess
from pathlib import Path
r=Path(__file__).resolve().parents[1];rs=list(csv.DictReader((r/'inputs/panel-all.tsv').open(),delimiter='\t'));g={x['accession']:x for x in rs}
new=[x['accession'] for x in rs if x['role'] not in ['existing_focal','shared_neighbour_clade','hold_poor_completeness_frameshifts']]
assert len(new)==8
pairs=[(q,t,'additional_to_type') for q in new for t in ['GCA_019046905.1','GCA_026344075.1']]+[('GCA_034107645.1','GCA_050773425.1','118C5_assembly_sensitivity'),('GCA_030510045.1','GCA_043966805.1','boundary_to_boundary')]
rows=[]
for q,t,reason in pairs:
 a,b=sorted([q,t]);folder=r/'anim/fasta'/(a+'__'+b);folder.mkdir(parents=True)
 for acc in [a,b]:(folder/(acc+'.fna')).symlink_to(g[acc]['source'])
 rows.append(dict(a=a,b=b,query=q,reference=t,selection_reasons=reason,reused=False,source='new'))
with (r/'anim/selected-pairs.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
(r/'anim/pairs').mkdir()
ledger=[];after='28999934'
for i,x in enumerate(rows,1):
 pair=x['a']+'__'+x['b'];cmd=f"cd ${PROJECT_ROOT} && .work/envs/task12-pyani/bin/python scripts/run_group_b_anim.py --run {r} --workers 1 --pair {pair}"
 job=subprocess.check_output(['python',str(r/'scripts/submit.py'),f't12-nem-ani-{i:02d}',cmd,'--cpus','1','--mem','1G','--time','10m','--after',after],text=True).strip();ledger.append(dict(pair=pair,job_id=job));after=job
 (r/'logs/anim-jobs.json').write_text(json.dumps(ledger,indent=2)+'\n');print(pair,job,flush=True)
