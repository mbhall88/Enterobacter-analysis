import pathlib,tempfile,json,hashlib,os,tarfile,unittest.mock
import stage_archive as a
root=pathlib.Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='p4b-archive-fixtures-',dir='/tmp') as temp:
 base=pathlib.Path(temp);out=base/'out';out.mkdir();work=base/'success';work.mkdir();(work/'stdout.json').write_text('{"synthetic":true}\n');(work/'stderr').write_text('')
 result=a.pack(out,work,0);assert not work.exists();assert a.sha(pathlib.Path(result['archive']))==result['archive_sha256']
 failed=base/'failed';failed.mkdir();(failed/'native.txt').write_text('preserve synthetic native output\n')
 original=tarfile.open
 def broken(name,mode='r',*args,**kwargs):
  if mode=='w:gz':raise OSError('synthetic archive creation failure')
  return original(name,mode,*args,**kwargs)
 with unittest.mock.patch('stage_archive.tarfile.open',side_effect=broken):
  try:a.pack(out,failed,1)
  except OSError:pass
  else:raise AssertionError('archive failure was accepted')
 assert failed.exists() and (out/('.recovery-task-0001.'+os.environ.get('SLURM_JOB_ID','local'))/'native.txt').read_bytes()==(failed/'native.txt').read_bytes()
 collision=base/'collision';collision.mkdir();(collision/'native.txt').write_text('collision evidence\n')
 with unittest.mock.patch('stage_archive.os.link',side_effect=FileExistsError('synthetic no-clobber collision')):
  try:a.pack(out,collision,2)
  except FileExistsError:pass
  else:raise AssertionError('publication collision accepted')
 assert collision.exists() and (out/'.task-0002.part').exists() and not (out/'task-0002.tar.gz').exists()
 proof={'gate_pass':True,'full_successful_archive_hash_and_readback':True,'archive_creation_failure_retains_native_work':True,'publication_collision_retains_native_work_and_partial':True,'native_or_scientific_queries':0,'source_sha256':a.sha(root/'scripts/stage_archive.py')};(root/'checks/group-task-archive-tests.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
