import json,pathlib,copy,datetime
import evaluation_common as c
case_root=c.ROOT/'results/tool-planning/2026-10-11/terobacter-issue17-opus-release-candidate/outputs/p1-replay-inputs';cases=[]
for path in sorted(case_root.iterdir()):
 if not path.is_dir():continue
 case=json.loads((path/'case.json').read_text());ev=json.loads((path/'evidence.json').read_text());inventory=json.loads((path/'inventory.json').read_text())['references'];records=list(map(json.loads,(path/'records.jsonl').read_text().splitlines())) if (path/'records.jsonl').exists() else [];config=c.p1_input(ev,inventory,records);config['warnings']=case['input'].get('warnings',False);config['scope']=case['input'].get('scope','not_assessed');actual=c.P1.evaluate(config);expected=c.P1.evaluate(case['input']);assert actual==expected,(case['id'],actual,expected);cases.append(case['id'])
assert len(cases)==50
fixtures=json.loads((c.P4A/'fixtures/scoring-cases.json').read_text());counts={}
for key,function in [('rescue',c.SCORER.score_rescue),('genomovar',c.SCORER.score_genomovar)]:
 for row in fixtures[key]:
  category,violations=function(row['record']);assert category==row['expected']['category'] and violations==row['expected']['violations'],row['id']
 counts[key]=len(fixtures[key])
for row in fixtures['rescue_raises']:
 try:c.SCORER.score_rescue(row['record'])
 except ValueError:pass
 else:raise AssertionError('Frozen refusal fixture accepted')
counts['rescue_raises']=len(fixtures['rescue_raises'])
for row in fixtures['summary']:
 result=c.SCORER.summarise_rescue(row['rows'])
 for key,value in row['expected'].items():assert (list(result[key]) if isinstance(result[key],tuple) else result[key])==value
counts['summary']=len(fixtures['summary'])
for row in fixtures['clusters']:
 result=c.SCORER.by_cluster(row['rows'],'cluster');result['ci95']=list(result['ci95']);assert result==row['expected']
counts['clusters']=len(fixtures['clusters'])
for row in fixtures['intervals']:
 result=c.SCORER.clopper_pearson(row['k'],row['n']);assert (list(result) if result is not None else None)==row['ci']
counts['intervals']=len(fixtures['intervals']);assert sum(counts.values())==84
# Actual typed transfer: explicit false cannot be treated as truthy text; true agreement structural,
# wrong rescue remains safety under either exposure value.
base=next(copy.deepcopy(x['record']) for x in fixtures['rescue'] if x['expected']['category']=='rescued_correct');base['truth'].pop('structural_species',None);base['truth']['ani_species_passing']=[]
exposure=[]
for informed in [False,True]:
 for correct in [True,False]:
  r=copy.deepcopy(base);r['truth']['anim_informed']=informed;r['truth']['species_id']=r['v2']['species_id'] if correct else 'TEST_WRONG_SPECIES';category,violations=c.SCORER.score_rescue(r);expected=('rescued_correct_structural' if informed else 'rescued_correct') if correct else 'rescued_wrong';assert category==expected and not violations;exposure.append(dict(anim_informed=informed,correct=correct,category=category))
truth=json.loads((c.CONT/'inputs/rs7-scoring-truth.json').read_text());assert len(truth)==11 and all(type(r['anim_informed']) is bool for r in truth) and sum(r['anim_informed'] for r in truth)==6
(c.CONT/'checks/adapter-tests.json').write_text(json.dumps(dict(gate_pass=True,P1_exact_integer_metric_adapter_cases=50,scorer_fixtures=counts,typed_exposure_tests=exposure,sealed_exposure_true=6,sealed_exposure_false=5,script_sha256=c.sha(__file__),adapter_sha256=c.sha(c.CONT/'scripts/evaluation_common.py'),fixture_sha256=c.sha(c.P4A/'fixtures/scoring-cases.json'),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2)+'\n');print('Adapter conformance PASS',len(cases),counts,exposure)
