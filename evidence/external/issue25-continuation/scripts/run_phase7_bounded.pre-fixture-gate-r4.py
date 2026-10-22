"""Step7 full saved-evidence replay, bounded within the frozen8CPU16GiB4h job."""
import argparse,pathlib,json,subprocess,hashlib,os,time,datetime,tarfile
from phase7_storage_budget import own_usage,parse_fileset_quota
root=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def write(p,x):
 with p.open('x') as h:json.dump(x,h,indent=2);h.write('\n')
ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
bp=root/'inputs/group-phase7-binding.json';b=json.loads(bp.read_text())
for p,want in b['source_sha256'].items():assert sha(p)==want,(p,'source changed')
assert json.loads((root/'checks/final-anim-summary.json').read_text())['gate_pass']
assert int(os.environ['SLURM_CPUS_PER_TASK'])>=8
out=root/'outputs/v2-r3';out.mkdir(parents=True,exist_ok=True)
if not args.pilot:
 p=json.loads((out/'task-0000.json').read_text());assert p['gate_pass'] and sha(out/'task-0000.tar.gz')==p['archive_sha256']
active={};remaining=[];manifest=json.loads((root/'inputs/final-v2-context-manifest.json').read_text());reused=[]
for task in range(1 if not args.pilot else 0,756 if not args.pilot else 1):
 archive=out/f'task-{task:04}.tar.gz';receipt=out/f'task-{task:04}.json'
 if archive.exists() or receipt.exists():
  assert not args.pilot and archive.exists() and receipt.exists() and not(out/f'.task-{task:04}.part').exists(),'Do not rerun successful/partial task'
  r=json.loads(receipt.read_text());assert r['gate_pass'] and r['failure'] is None and sha(archive)==r['archive_sha256'];prefix=f'task-{task:04}'
  with tarfile.open(archive) as t:
   hashes=json.load(t.extractfile(prefix+'/artifact-hashes.json'))
   for name,want in hashes.items():assert hashlib.sha256(t.extractfile(prefix+'/'+name).read()).hexdigest()==want
   stage=json.load(t.extractfile(prefix+'/stage.json'));assert stage['gate_pass'] and stage['completed_records']==stage['selected_records']==len(manifest[task*20:(task+1)*20]);assert stage['source_sha256']==b['source_sha256'][str(root/'scripts/run_v2.py')] if str(root/'scripts/run_v2.py') in b['source_sha256'] else stage['source_sha256']==next(v for k,v in b['source_sha256'].items() if k.endswith('/scripts/run_v2.py'));assert stage['manifest_sha256']==sha(root/'inputs/final-v2-context-manifest.json');assert {x['record_id'] for x in stage['results']}=={x['record_id'] for x in manifest[task*20:(task+1)*20]};assert all(x['model_pass'] and x['format_pass'] for x in stage['results'])
  reused.append({'task':task,'archive_sha256':r['archive_sha256'],'completed_records':r['completed_records']});continue
 remaining.append(task)
write(root/'checks/phase7-r4-success-reuse.json',{'gate_pass':True,'reused':reused,'task0_successful_pilot_reused':not args.pilot,'source_sha256':sha(__file__),'binding_sha256':sha(bp),'no_successful_case_replayed':True})
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');failure=None
try:
 while remaining or active:
  for task,proc in list(active.items()):
   if proc.poll() is None:continue
   assert proc.returncode==0,('Saved evidence task failed',task,proc.returncode)
   receipt=json.loads((out/f'task-{task:04}.json').read_text());assert receipt['gate_pass'] and sha(out/f'task-{task:04}.tar.gz')==receipt['archive_sha256'];del active[task]
   print(json.dumps({'task':task,'gate_pass':True,'completed_records':receipt['completed_records']}),flush=True)
  assert not (out/'gate-stop.json').exists(),'R2/P1/format gate stopped'
  while remaining and len(active)<4:
   task=remaining[0];assert not (out/f'task-{task:04}.tar.gz').exists() and not (out/f'.task-{task:04}.part').exists(),'No completed/partial task rerun'
   u=own_usage(root);raw=subprocess.run(['/usr/lpp/mmfs/bin/mmlsquota','-j','S0256','scratch'],capture_output=True,text=True,check=True).stdout;q=parse_fileset_quota(raw);reserve=(4*256+768)*1024**2;inode_need=10000+4*1000
   proof={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task':task,'own_usage':u,'quota':q,'reserved_peak_bytes':reserve,'cap_bytes':b['cap_bytes'],'admitted':u['allocated_bytes']+reserve<=b['cap_bytes'] and q['block_hard_bytes']-q['block_used_bytes']-q['block_in_doubt_bytes']>=reserve and q['files_hard']-q['files_used']-q['files_in_doubt']>=inode_need and 22581+u['file_and_directory_entries']+4000<50000};write(out/f'admission-task-{task:04}.json',proof);assert proof['admitted'],'Storage admission refused before further replay'
   probe=out/f'.write-probe-{task:04}';data=b'P4b step7 storage admission\n'*8192
   with probe.open('xb') as h:h.write(data);h.flush();os.fsync(h.fileno())
   assert probe.read_bytes()==data;probe.unlink()
   remaining.pop(0);active[task]=subprocess.Popen(['python3',str(root/'scripts/run_v2.py'),'--task',str(task)],env=env)
  time.sleep(1)
except Exception as e:
 failure=str(e)
 # Already admitted tasks finish with their own full archive/failure receipts; no new task begins.
 for task,proc in active.items():proc.wait()
finally:
 write(root/'checks'/('phase7-pilot-driver-r4.json' if args.pilot else 'phase7-full-driver-r4.json'),{'gate_pass':failure is None,'failure':failure,'remaining_tasks':remaining,'binding_sha256':sha(bp),'source_sha256':sha(__file__),'job_id':os.environ['SLURM_JOB_ID'],'new_measurements':False})
if failure:raise RuntimeError(failure)
