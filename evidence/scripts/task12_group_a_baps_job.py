#!/usr/bin/env python3
"""One authorized fastbaps sensitivity per SLURM array element, with exact command ledger."""
import csv,json,os,subprocess,sys,time
from pathlib import Path
r=Path(sys.argv[1]).resolve(); i=int(os.environ.get('SLURM_ARRAY_TASK_ID',sys.argv[2] if len(sys.argv)>2 else '0'))
with (r/'inputs/datasets.tsv').open() as f: row=list(csv.DictReader(f,delimiter='\t'))[i]
name=row['dataset']+'--'+row['prior'];out=r/'fastbaps'/name
cmd=[str(Path('.work/envs/task12-populations/bin/Rscript').resolve()),'scripts/task12_group_a_fastbaps.R',row['alignment'],row['prior'],str(out)]
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
start=time.time(); rec={'argv':cmd,'cwd':os.getcwd(),'array_index':i,'job':os.environ.get('SLURM_JOB_ID'),'started_epoch':start}
p=r/'logs'/(name+'.command.json');p.write_text(json.dumps(rec,indent=2))
with (r/'logs'/(name+'.log')).open('w') as f: result=subprocess.run(['/usr/bin/time','-v',*cmd],stdout=f,stderr=subprocess.STDOUT)
rec.update(exit_status=result.returncode,elapsed_seconds=time.time()-start);p.write_text(json.dumps(rec,indent=2));sys.exit(result.returncode)
