"""Full per-task archival publication; keep failed native work and partial archive."""
import pathlib,hashlib,json,tarfile,shutil,os,socket,sys

def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()

def pack(out,work,task):
 out=pathlib.Path(out);work=pathlib.Path(work);prefix=f'task-{task:04d}';final=out/(prefix+'.tar.gz');part=out/('.'+prefix+'.part');assert not final.exists() and not part.exists()
 hashes={str(q.relative_to(work)):sha(q) for q in work.rglob('*') if q.is_file()};(work/'artifact-hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
 try:
  with tarfile.open(part,'w:gz') as t:t.add(work,arcname=prefix)
  with tarfile.open(part) as t:
   for name,want in hashes.items():assert hashlib.sha256(t.extractfile(prefix+'/'+name).read()).hexdigest()==want
  expected=sha(part);os.link(part,final);assert sha(final)==expected;part.unlink();shutil.rmtree(work);return {'archive':str(final),'archive_sha256':expected,'member_files':len(hashes)}
 except Exception:
  recovery=out/('.recovery-'+prefix+'.'+os.environ.get('SLURM_JOB_ID','local'))
  try:
   shutil.copytree(work,recovery)
   for name,want in hashes.items():assert sha(recovery/name)==want
   where=str(recovery)
  except Exception:where=socket.gethostname()+':'+str(work)
  print('P4b task archive publication failed; native work and partial archive retained at '+where,file=sys.stderr)
  raise
