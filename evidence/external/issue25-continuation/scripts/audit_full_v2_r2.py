"""Read back all actual R2/P1/format task evidence and exact record denominators."""
import pathlib,json,tarfile,hashlib,datetime,collections
from lower_rank_checks import assert_no_genomovar
ROOT=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
manifest=json.loads((ROOT/'inputs/final-v2-context-manifest.json').read_text());assert len(manifest)==15109;out=ROOT/'outputs/v2-r3';results=[];index={};pins={};scorer_rows=[]
for i in range((len(manifest)+19)//20):
 prefix=f'task-{i:04d}';archive=out/(prefix+'.tar.gz');receipt=json.loads((out/(prefix+'.json')).read_text());assert receipt['gate_pass'] and receipt['failure'] is None and sha(archive)==receipt['archive_sha256'];pins[str(archive)]=receipt['archive_sha256']
 with tarfile.open(archive) as t:
  hashes=json.load(t.extractfile(prefix+'/artifact-hashes.json'))
  for name,want in hashes.items():assert hashlib.sha256(t.extractfile(prefix+'/'+name).read()).hexdigest()==want
  stage=json.load(t.extractfile(prefix+'/stage.json'));assert stage['gate_pass'] and stage['selected_records']==stage['completed_records'];results.extend(stage['results'])
  for result in stage['results']:
   rid=result['record_id'];case=prefix+'/'+rid.replace('|','--');assert rid not in index;v2=json.load(t.extractfile(case+'/v2.stdout.json'));v1=json.load(t.extractfile(case+'/v1.stdout.json'));comparison=json.load(t.extractfile(case+'/p1-comparison.json'));formats=json.load(t.extractfile(case+'/format-check.json'));assert comparison['exact_match'] and not comparison['differences'] and all(formats.values());assert v2['input']['query_sha256']==result['query_sha256']
   if not v2['trigger']:assert v2['species_id']==v1['assignment']['species_id'] and v2['subspecies_id']==v1['assignment']['subspecies_id'],'RG2 retained call mismatch'
   if v2['call_basis']=='ANIm_assisted':assert v2['subspecies_id'] is None;assert_no_genomovar(v2)
   index[rid]={'archive':str(archive),'archive_sha256':receipt['archive_sha256'],'v2_member':case+'/v2.stdout.json','v1_member':case+'/v1.stdout.json','p1_member':case+'/p1-comparison.json','scorer_member':case+'/scorer-records.json','query_sha256':result['query_sha256'],'dataset':result['dataset'],'setting':result['setting']}
   records=json.load(t.extractfile(case+'/scorer-records.json'))
   for row in records:assert row['scoring_deferred_to_step9'] is True;scorer_rows.append(row)
assert len(results)==len(index)==len(manifest)==15109
assert {r['record_id'] for r in results}=={r['record_id'] for r in manifest};assert all(r['model_pass'] and r['format_pass'] for r in results)
(ROOT/'inputs/full-v2-result-index.json').write_text(json.dumps(index,separators=(',',':'))+'\n');(ROOT/'inputs/rescue-scorer-inputs.json').write_text(json.dumps(scorer_rows,separators=(',',':'))+'\n');summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':True,'contexts':15109,'datasets':dict(collections.Counter(r['dataset'] for r in results)),'model_and_format_exact':15109,'retained_v1_calls_and_rescue_lower_rank_gates':True,'scoring_contexts':dict(collections.Counter(r['record']['dataset'] for r in scorer_rows)),'scoring_executed':False,'input_archive_sha256':pins,'v2_result_index_sha256':sha(ROOT/'inputs/full-v2-result-index.json'),'scorer_input_sha256':sha(ROOT/'inputs/rescue-scorer-inputs.json'),'source_sha256':sha(pathlib.Path(__file__))};(ROOT/'checks/full-v2-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='input_archive_sha256'}))
