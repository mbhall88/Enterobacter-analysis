#!/usr/bin/env python3
"""One ssubmit invocation per ANIm pair and scope; bounded concurrent scheduler jobs."""
import argparse,json,subprocess,time
from pathlib import Path
import task12_scope_s1_s2 as s
R=s.RUN
p=argparse.ArgumentParser();p.add_argument('--submit',action='store_true');a=p.parse_args()
for i,j in enumerate(s.rows(R/'s2-anim-jobs.tsv')):
 if j['reused']=='True':continue
 name=f"scope2-{i:03d}-{j['mode']}";folder=R/'scheduler'/name;folder.mkdir(exist_ok=True)
 cmd=['ssubmit','--json',name,'python '+str(s.ROOT/'scripts/task12_scope_s2_sensitivity.py')+' run --pair '+j['pair']+' --mode '+j['mode'],'--mem','2G','--time','10m','--output',str(folder/'job.%j.out'),'--error',str(folder/'job.%j.err'),'--','--cpus-per-task=1']
 def invoke(mode):
  args=cmd.copy()
  if mode:args.insert(1,mode)
  r=subprocess.run(args,capture_output=True,text=True);key={'--dry-run':'plan','--test-only':'test','':'submit'}[mode]
  (folder/(key+'.stdout')).write_text(r.stdout);(folder/(key+'.stderr')).write_text(r.stderr)
  (folder/(key+'.command.json')).write_text(json.dumps(args,indent=2));x=json.loads(r.stdout)
  assert r.returncode==0 and x['schema_version']==1 and x['ok'] and x['operation']==key,(args,r.stderr)
  return x
 invoke('--dry-run');invoke('--test-only')
 if a.submit:
  receipts=list((R/'scheduler').glob('scope2-*/submit.stdout'))
  while True:
   ids=[json.loads(f.read_text())['submission']['job_id'] for f in receipts]
   if not ids:break
   r=subprocess.run(['squeue','-h','-j',','.join(ids),'-o','%i'],capture_output=True,text=True);assert r.returncode==0
   if len(r.stdout.strip().splitlines())<12:break
   time.sleep(15)
  x=invoke('');print(j['pair'],j['mode'],x['submission']['job_id'],flush=True)
