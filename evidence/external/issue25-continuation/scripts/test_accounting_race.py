"""Narrow accounting-race fixtures; quota/I/Oerrors must not be hidden."""
import pathlib,tempfile,unittest.mock,json,hashlib,datetime,errno,shutil
import phase7_storage_budget as b
root=pathlib.Path(__file__).resolve().parents[1];folder=pathlib.Path(tempfile.mkdtemp(prefix='p4b-accounting-test-',dir='/tmp'));(folder/'stable').write_bytes(b'x'*8192);ghost=folder/'.archive.part';ghost.write_bytes(b'y'*8192);original=pathlib.Path.lstat
try:
 def vanished(p,*a,**kw):
  if p==ghost:raise FileNotFoundError(errno.ENOENT,'Concurrent publication removed own temporary archive',str(p))
  return original(p,*a,**kw)
 with unittest.mock.patch.object(pathlib.Path,'lstat',vanished):u=b.own_usage(folder)
 assert u['file_and_directory_entries']==1 and u['allocated_bytes']>=8192
 errors=[]
 for code in [errno.EDQUOT,errno.EIO,errno.EACCES]:
  def fatal(p,*a,**kw):
   if p==ghost:raise OSError(code,'Fixture fatal storage error')
   return original(p,*a,**kw)
  try:
   with unittest.mock.patch.object(pathlib.Path,'lstat',fatal):b.own_usage(folder)
  except OSError as e:assert e.errno==code;errors.append(code)
  else:raise AssertionError('Fatal storage error swallowed')
 out=root/'checks/accounting-race-fixtures.json';assert not out.exists();out.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':True,'disappearing_temp_name_safe':True,'fatal_errnos_preserved':errors,'source_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'storage_helper_sha256':hashlib.sha256((root/'scripts/phase7_storage_budget.py').read_bytes()).hexdigest(),'biological_work':False},indent=2)+'\n');print('Accounting race fixtures PASS')
finally:shutil.rmtree(folder)
