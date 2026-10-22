"""Step8 frozen genomovar actual R1/R5 replay, bounded within the frozen8CPU16GiB4h job."""
import argparse,pathlib,json,subprocess,hashlib,os,time,datetime
from phase7_storage_budget import own_usage,parse_fileset_quota
root=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def write(p,x):
 with p.open('x') as h:json.dump(x,h,indent=2);h.write('\n')
ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');args=ap.parse_args()
bp=root/'inputs/group-genomovar-binding.json';b=json.loads(bp.read_text())
for p,want in b['source_sha256'].items():assert sha(p)==want,(p,'source changed')
assert json.loads((root/'checks/step7-summary.json').read_text())['gate_pass']
assert int(os.environ['SLURM_CPUS_PER_TASK'])>=8
out=root/'outputs/genomovar-r3';out.mkdir(parents=True,exist_ok=True)
if not args.pilot:
 p=json.loads((out/'task-0000.json').read_text());assert p['gate_pass'] and sha(out/'task-0000.tar.gz')==p['archive_sha256']
active={};remaining=list(range(1 if not args.pilot else 0,43 if not args.pilot else 1));env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');failure=None
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
   remaining.pop(0);active[task]=subprocess.Popen(['python3',str(root/'scripts/run_genomovar_sets.py'),'--task',str(task)],env=env)
  time.sleep(1)
except Exception as e:
 failure=str(e)
 # Already admitted tasks finish with their own full archive/failure receipts; no new task begins.
 for task,proc in active.items():proc.wait()
finally:
 write(root/'checks'/('genomovar-pilot-driver-r3.json' if args.pilot else 'genomovar-full-driver-r3.json'),{'gate_pass':failure is None,'failure':failure,'remaining_tasks':remaining,'binding_sha256':sha(bp),'source_sha256':sha(__file__),'job_id':os.environ['SLURM_JOB_ID'],'new_measurements':'Only276 fixed fresh R1 derivative inputs;565 actual source R1 and173 actual R5 reused'})
if failure:raise RuntimeError(failure)
