"""Authenticate complete frozenstep7 classification and18 native integration gates."""
import pathlib,json,hashlib,tarfile,datetime,collections
root=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
v2=json.loads((root/'checks/full-v2-summary.json').read_text());assert v2['gate_pass'] and v2['contexts']==15109
manifest=json.loads((root/'inputs/integrated-context-manifest.json').read_text());assert len(manifest)==18;records=[];pins={};out=root/'outputs/integrated'
for i,row in enumerate(manifest):
 prefix=f'task-{i:04}';archive=out/(prefix+'.tar.gz');receipt=json.loads((out/(prefix+'.json')).read_text());assert receipt['gate_pass'] and receipt['failure'] is None and sha(archive)==receipt['archive_sha256'];pins[str(archive)]=receipt['archive_sha256']
 with tarfile.open(archive) as t:
  hashes=json.load(t.extractfile(prefix+'/artifact-hashes.json'))
  for name,want in hashes.items():assert hashlib.sha256(t.extractfile(prefix+'/'+name).read()).hexdigest()==want
  stage=json.load(t.extractfile(prefix+'/stage.json'));assert stage['gate_pass'] and stage['completed_records']==1 and stage['query_sha256']==row['query_sha256'];r=stage['result'];assert r['record_id']==row['record_id'] and r['P1_pass'];assert r['split_pass'] if i<15 else r['custom_unvalidated'] and r['explicit_default_decision_parity'];records.append(r)
assert sum(len(r['integer_checks']) for r in records[:15])==88
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':True,'saved_evidence_contexts':15109,'native_integrations':15,'postseal_RS7_integrations':11,'EN5_split_integer_orientations':88,'explicit_default_custom_cases':3,'CU2_default_decisions_equal':True,'scoring_executed':False,'source_sha256':sha(__file__),'full_v2_summary_sha256':sha(root/'checks/full-v2-summary.json'),'integration_manifest_sha256':sha(root/'inputs/integrated-context-manifest.json'),'archive_sha256':pins,'results':records};f=root/'checks/step7-summary.json';assert not f.exists();f.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k not in ['archive_sha256','results']}))
