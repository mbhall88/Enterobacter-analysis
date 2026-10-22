"""Synthetic fixtures through the production rescue_record adapter; no biological scoring."""
import copy,json,datetime
import evaluation_common as c
fixtures=json.loads((c.P4A/'fixtures/scoring-cases.json').read_text())
base=next(copy.deepcopy(x['record']) for x in fixtures['rescue'] if x['expected']['category']=='rescued_correct')
base['truth'].pop('structural_species',None);base['truth']['ani_species_passing']=[]
v1={};v2={}
for f in c.FIELDS:
 if f['record']!='rescue' or f['supplied_by']!='terobacter':continue
 value=c.get(base,f['scorer_field'])
 if f['rule'].startswith('map '):
  policy=f['rule'].split()[1];value=[next(reason for (p,field,reason),contract in c.VOCAB.items() if p==policy and field=='reasons' and contract==r) for r in value]
 c.put(v1 if f['terobacter_output']=='v1-report' else v2,f['terobacter_path'],value)
results=[]
for exposure in [False,True]:
 for correct in [False,True]:
  truth=copy.deepcopy(base['truth']);truth['anim_informed']=exposure;truth['species_id']=base['v2']['species_id'] if correct else 'TEST_WRONG_SPECIES'
  r=c.rescue_record(v1,v2,base['manifest'],truth,'SYNTHETIC-BOOL-TRANSFER','synthetic','as_run');assert type(r['truth']['anim_informed']) is bool and r['truth']['anim_informed'] is exposure
  category,violations=c.SCORER.score_rescue(r);expected=('rescued_correct_structural' if exposure else 'rescued_correct') if correct else 'rescued_wrong';assert category==expected and not violations;results.append({'anim_informed':exposure,'correct':correct,'category':category,'through_production_adapter':True})
truth=copy.deepcopy(base['truth']);truth['anim_informed']='False'
try:c.rescue_record(v1,v2,base['manifest'],truth,'SYNTHETIC-BOOL-STRING','synthetic','as_run')
except AssertionError:rejected=True
else:raise AssertionError('Production adapter accepted string False')
out=c.CONT/'checks/production-truth-transfer-fixtures.json';assert not out.exists();out.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':True,'results':results,'string_false_rejected':rejected,'source_sha256':c.sha(__file__),'production_adapter_sha256':c.sha(c.CONT/'scripts/evaluation_common.py'),'fixture_sha256':c.sha(c.P4A/'fixtures/scoring-cases.json'),'biological_scoring':False},indent=2)+'\n');print('Production truth transfer fixtures PASS')
