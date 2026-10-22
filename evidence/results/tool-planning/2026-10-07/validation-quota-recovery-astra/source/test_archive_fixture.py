import tempfile,json,hashlib,csv,sys,importlib.util,tarfile
from pathlib import Path
source=Path('/tmp/enterobacter-validation-quota-recovery-2026-10-07/archive_completed.py')
spec=importlib.util.spec_from_file_location('archive_fixture_target',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
base=Path(tempfile.mkdtemp(prefix='enterobacter-quota-fixture-'));run=base/'run';local=base/'local'
(run/'inputs').mkdir(parents=True);(run/'outputs').mkdir();local.mkdir()
for filename in ['execution-pins.tsv','driver-pins.tsv']:(run/'inputs'/filename).write_text('path\tsha256\n')
(run/'inputs/execution-manifest.tsv').write_text('coded_id\nSYNTHETIC-GOOD\nSYNTHETIC-FAILED\n')
d=run/'queries/SYNTHETIC-GOOD/attempt01';(d/'analyze/nested').mkdir(parents=True)
(d/'finished.json').write_text(json.dumps(dict(exit_status=0)))
files={}
for i in range(130):
 p=d/'analyze/nested'/f'{i}.txt';p.write_text('synthetic bytes '+str(i));files[str(p.relative_to(d/'analyze'))]=m.sha(p)
for name in m.KEEP-{'run.json'}:
 p=d/'analyze'/name;p.write_text('retained synthetic '+name);files[name]=m.sha(p)
(d/'analyze/run.json').write_text(json.dumps(dict(artifacts_sha256=files)))
bad=run/'queries/SYNTHETIC-FAILED/attempt01';bad.mkdir(parents=True);(bad/'finished.json').write_text(json.dumps(dict(exit_status=2)));(bad/'partial').write_text('do not touch')
m.RUN=run;m.LOCAL=local;m.verify_object=lambda:None;m.audit_query=lambda row,root:dict(coded_id=row['coded_id'],exit_status=0)
m.main()
assert (bad/'partial').read_text()=='do not touch'
assert all((d/'analyze'/name).is_file() for name in m.KEEP)
assert not (d/'analyze/nested').exists()
archive=run/'outputs/completed-native-ancillary.tar.gz'
with tarfile.open(archive,'r:gz') as tar:
 members=tar.getmembers();assert len(members)==130
 for member in members:
  rel=member.name.split('/analyze/')[1]
  assert hashlib.sha256(tar.extractfile(member).read()).hexdigest()==files[rel]
assert m.sha(archive)==m.sha(local/'completed-native-ancillary.tar.gz')
assert (run/'outputs/archive-verified-checkpoint.json').is_file()
record=dict(status='pass',fixture_root=str(base),archive_members=130,kept_core_files=7,failed_attempt_untouched=True,source_sha256=m.sha(source),real_data_modified=False)
Path('/tmp/enterobacter-validation-quota-recovery-2026-10-07/synthetic-archive-test.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
