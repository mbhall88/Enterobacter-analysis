#!/usr/bin/env python3
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--submit',action='store_true');a=p.parse_args()
r=Path('results/ani-cases/2026-10-01/task12-astra-genus-tree')
jobs=[('fasttree',4,'8G','2h','bash '+str(r.resolve()/'run-fasttree.sh'))]+[(x+'-iqtree',16,'32G','6h','bash '+str(r.resolve()/'run-iqtree.sh')+' '+x) for x in ['primary-full80','canonical-masked']]
for label,cpu,mem,time,cmd in jobs:
 args=['ssubmit','--json']+([] if a.submit else ['--dry-run'])+['genus_'+label,cmd,'--mem',mem,'--time',time,'--output',str(r/'logs'/f'{label}-%j.out'),'--error',str(r/'logs'/f'{label}-%j.err'),'--',f'--cpus-per-task={cpu}']
 out=subprocess.run(args,capture_output=True,text=True)
 status='submit' if a.submit else 'plan'
 (r/'logs'/f'{label}-{status}.stderr').write_text(out.stderr)
 data=json.loads(out.stdout);(r/f'{label}-{status}.json').write_text(json.dumps({'argv':args,'response':data},indent=2)+'\n')
 assert out.returncode==0 and data['ok'] and data['schema_version']==1
 print(label,data.get('submission',data['plan']['job']))
