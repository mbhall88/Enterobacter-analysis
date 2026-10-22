"""Storage-only entry: verify accepted binding, then exec unchanged frozen wrapper in local/tmp."""
import hashlib,json,os,pathlib,subprocess,sys,tempfile
root=pathlib.Path(__file__).resolve().parents[1];base=root.parent
bfile=base/'inputs/storage-binding.json';b=json.loads(bfile.read_text());receipt=json.loads((base/'inputs/coordinator-acceptance.json').read_text())
assert receipt['accepted'] is True and hashlib.sha256(bfile.read_bytes()).hexdigest()==receipt['binding_sha256']
assert len(sys.argv)==4 and pathlib.Path(sys.argv[1]).resolve()==pathlib.Path(b['frozen_config']).resolve()
for key,pathkey in [('config_sha256','frozen_config'),('wrapper_sha256','frozen_wrapper')]:assert hashlib.sha256(pathlib.Path(b[pathkey]).read_bytes()).hexdigest()==b[key]
for binary,want in b['candidate_binaries'].items():assert hashlib.sha256(pathlib.Path(binary).read_bytes()).hexdigest()==want
adapter=json.loads(pathlib.Path(b['frozen_config']).read_text())['anim_pair_argv'][1]
assert hashlib.sha256(pathlib.Path(adapter).read_bytes()).hexdigest()==b['adapter_sha256']
assert pathlib.Path(sys.argv[3]).resolve().is_relative_to(root.resolve())
fs=subprocess.run(['findmnt','-n','-T','/tmp','-o','FSTYPE'],capture_output=True,text=True,check=True).stdout.strip()
assert fs and fs not in ['gpfs','nfs','nfs4','lustre','cifs','fuse.sshfs'],fs
st=os.statvfs('/tmp');assert st.f_bavail*st.f_frsize>=512*1024**2
with tempfile.TemporaryDirectory(prefix='issue25-local-proof-',dir='/tmp') as d:
 x=pathlib.Path(d)/'probe';data=b'P4b local tmp readback\n'*2048
 with x.open('xb') as h:h.write(data);h.flush();os.fsync(h.fileno())
 assert x.read_bytes()==data
print(json.dumps({'storage_only_launcher':True,'node_tmp_root':'/tmp','filesystem':fs,'available_bytes':st.f_bavail*st.f_frsize,'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'task_index':os.environ.get('SLURM_ARRAY_TASK_ID')}),flush=True)
os.environ['TMPDIR']='/tmp'
os.execv('/usr/bin/bash',['bash',b['frozen_wrapper'],*sys.argv[1:]])
