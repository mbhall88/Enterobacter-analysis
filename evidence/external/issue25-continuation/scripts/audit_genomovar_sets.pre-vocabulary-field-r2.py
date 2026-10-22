"""Authenticate every fixed actual-genomovar record before any biological scoring."""
import pathlib,json,tarfile,hashlib,collections,datetime
root=pathlib.Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
assert json.loads((root/'checks/step7-summary.json').read_text())['gate_pass'];allrows=[];index={};pins={}
for name,manifestfile,size in [('genomovar','genomovar-context-plan.json',20),('gv3','gv3-native-context-plan.json',50)]:
 manifest=json.loads((root/'inputs'/manifestfile).read_text());ids={r['record_id'] for r in manifest};assert len(ids)==len(manifest)
 for i in range((len(manifest)+size-1)//size):
  prefix=f'task-{i:04}';out=root/'outputs'/name;archive=out/(prefix+'.tar.gz');receipt=json.loads((out/(prefix+'.json')).read_text());assert receipt['gate_pass'] and receipt['failure'] is None and sha(archive)==receipt['archive_sha256'];pins[str(archive)]=receipt['archive_sha256']
  with tarfile.open(archive) as t:
   hashes=json.load(t.extractfile(prefix+'/artifact-hashes.json'))
   for member,want in hashes.items():assert hashlib.sha256(t.extractfile(prefix+'/'+member).read()).hexdigest()==want
   stage=json.load(t.extractfile(prefix+'/stage.json'));expected=manifest[i*size:(i+1)*size];assert stage['gate_pass'] and stage['selected_records']==stage['completed_records']==len(expected);assert {r['record_id'] for r in stage['results']}=={r['record_id'] for r in expected}
   for result in stage['results']:
    rid=result['record_id'];case=prefix+'/'+rid;comparison=json.load(t.extractfile(case+'/full-p3-comparison.json'));assert not comparison['model_differences'];bridge=comparison['oracle_package_bridge'];assert bridge['actual_parent_search_bindings_checked_before_alias'] and bridge['all_25_173_definition_identities_checked'];rec=json.load(t.extractfile(case+'/scorer-record.json'));assert rec==comparison['record'] and rec['record_id']==rid;assert type(rec['diagnostic']) is bool and type(rec['self_member']) is bool and type(rec['singleton_source']) is bool;assert rid not in index
    if rec['diagnostic']:assert rec['gv']['assigned'] is None
    assigned=rec['gv'].get('diagnostic_would_assign') if rec['diagnostic'] else rec['gv']['assigned']
    if assigned:assert rec['parent']['status']=='assigned' and (rec['parent']['species_id'],rec['parent']['subspecies_id'])==('ECTX-SP025','ECTX-SS007') and rec['parent']['call_basis']!='ANIm_assisted' and assigned not in rec['removed_groups']
    source={'archive':str(archive),'archive_sha256':receipt['archive_sha256'],'scorer_member':case+'/scorer-record.json','comparison_member':case+'/full-p3-comparison.json','query_sha256':result['query_sha256']};allrows.append({'record':rec,'source':source,'scoring_deferred_to_step9':True});index[rid]=source
 expectedcounts={'GV1b':173,'GV2':173,'GV6':173,'GV4':322} if name=='genomovar' else {'GV3-public':749,'GV3-study':165};actual=collections.Counter(x['record']['dataset'] for x in allrows if x['record']['record_id'] in ids);assert actual==expectedcounts
assert len(allrows)==1755 and sum(r['record']['self_member'] for r in allrows)==165
for name,data in [('genomovar-scorer-inputs.json',allrows),('genomovar-result-index.json',index)]:
 f=root/'inputs'/name;assert not f.exists();f.write_text(json.dumps(data,separators=(',',':'))+'\n')
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':True,'records':1755,'datasets':dict(collections.Counter(r['record']['dataset'] for r in allrows)),'self_members':165,'full_P3_actual_parent_search_and_cell_replays':1755,'biological_scoring_executed':False,'input_archive_sha256':pins,'source_sha256':sha(__file__),'scorer_input_sha256':sha(root/'inputs/genomovar-scorer-inputs.json')};f=root/'checks/genomovar-sets-summary.json';assert not f.exists();f.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='input_archive_sha256'}))
