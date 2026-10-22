"""Frozen model/scorer adapters; authenticate input and map the accepted vocabularies."""
import csv,hashlib,importlib.util,pathlib
from fractions import Fraction
CONT=pathlib.Path(__file__).resolve().parents[1];RUN=CONT.parents[1];ROOT=pathlib.Path('${PROJECT_ROOT}');P4A=ROOT/'results/tool-planning/2026-10-08/p4a-evaluation-design-claude'
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def load(name,p):
 assert sha(p)==SOURCE_PINS[str(p)], 'Frozen source hash changed';spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
SOURCE_PINS={str(ROOT/'results/tool-planning/2026-10-08/issue20-calling-contract-astra/source/contract_model.py'):'76ea8b6195afc92307fcdcaa7b56be4f5531d03891e01942418e963b2edac438',str(ROOT/'results/tool-planning/2026-10-08/p3a-genomovar-contract-astra-1609/source/spec_model.py'):'17b9b60eebc7ab3a22e7c6ea6f6b56c3c8c70ea89bcfaab80cf3bb49296fcd0d',str(P4A/'scripts/score_changed_behaviour.py'):'a3bab5d8824cf91d201b640850a7f24d8e2867d694f378524046fef18a4377d3'}
for filename,expected in [('reason-vocabulary.tsv','b988ab4f895741cd57dcaac765707d38fd7187c49f3401deb0b28e66f67431d1'),('scorer-fields.tsv','2f738d2cba149cf08d1c1c7b16966e73708f2279518185b9ec73f3936d5c611e')]:assert sha(RUN/'build/source/docs/contracts'/filename)==expected
P1=load('p4b_frozen_p1',ROOT/'results/tool-planning/2026-10-08/issue20-calling-contract-astra/source/contract_model.py');P3=load('p4b_frozen_p3',ROOT/'results/tool-planning/2026-10-08/p3a-genomovar-contract-astra-1609/source/spec_model.py');SCORER=load('p4b_frozen_scorer',P4A/'scripts/score_changed_behaviour.py')
VOCAB={};FIELDS=[]
with (RUN/'build/source/docs/contracts/reason-vocabulary.tsv').open() as h:
 for r in csv.DictReader(h,delimiter='\t'):
  key=(r['policy'],r['field'],r['terobacter_reason']);assert key not in VOCAB;VOCAB[key]=r['contract_reason']
with (RUN/'build/source/docs/contracts/scorer-fields.tsv').open() as h:FIELDS=list(csv.DictReader(h,delimiter='\t'))
def translate(policy,values):return list(dict.fromkeys(VOCAB[policy,'reasons',v] for v in values))
def get(x,path):
 for field in path.split('.'):x=x[field]
 return x
def put(x,path,value):
 *parents,last=path.split('.')
 for field in parents:x=x.setdefault(field,{})
 x[last]=value

def p1_input(evidence,inventory,records,mask=(),custom=None):
 rows={r['reference_id']:r for r in evidence['ranked_reference_evidence'] if r['reference_id'] not in mask};refs=[]
 for r in inventory:
  if r['reference_id'] in mask:continue
  row=rows.get(r['reference_id']);refs.append(dict(id=r['reference_id'],species=r.get('species_id') or None,subspecies=r.get('subspecies_id') or None,lineage=r.get('lineage_id') or None,ani=row['ani_percent'] if row else None,aq=row['aligned_fraction_query_percent'] if row else None,ar=row['aligned_fraction_reference_percent'] if row else None))
 assert set(rows)<={r['id'] for r in refs};anim={}
 for r in records:
  rid=r['reference']['id'];assert rid not in mask and rid not in anim;assert r['query']['sha256']==evidence['provenance']['query_sha256'];directions=[]
  for orientation in ['QR','RQ']:
   d=r['orientations'].get(orientation)
   if d is None:continue
   item=dict(id=orientation,status={'aligned':'complete','no_alignment':'no_alignment'}.get(d['status'],'error'))
   if item['status']=='complete':
    m=d['metrics'];assert type(m['D']) is int and type(m['E']) is int and m['D']>0;item.update(ani=str(Fraction(100*(m['D']-2*m['E']),m['D'])),aq=str(Fraction(100*m['query_covered'],m['query_length'])),ar=str(Fraction(100*m['reference_covered'],m['reference_length'])))
   directions.append(item)
  anim[rid]=dict(directions=directions)
 config=dict(references=refs,anim=anim,warnings=bool(evidence.get('external_problem_warnings',[])),scope=evidence.get('scope_assessment',{}).get('outcome','not_assessed'))
 if custom:config['custom']=custom
 return config

def compare_p1(output,expected):
 actual=dict(species=output['species_id'],subspecies=output['subspecies_id'],basis=output['call_basis'],reasons=sorted(translate('conservative-assembly-anim-v2',output['reasons'])),trigger=output['trigger'],mode=output['policy_mode'],display=[r['reference_id'] for r in output['top_references']],challenges=sorted(output['challenges']),competitor_evaluation=output['competitor_evaluation'],validation_status=output['validation_status'])
 return {k:dict(expected=v,actual=actual.get(k)) for k,v in expected.items() if actual.get(k)!=v}

def rescue_record(v1,v2,trigger,truth,record_id,dataset,setting):
 record=dict(record_id=record_id,dataset=dataset,setting=setting,truth=truth,manifest=dict(state=trigger['state'],parent_id=trigger['parent_id'] or None));sources={'v1-report':v1,'v2-report':v2}
 for f in FIELDS:
  if f['record']!='rescue' or f['supplied_by']!='terobacter':continue
  value=get(sources[f['terobacter_output']],f['terobacter_path'])
  if f['rule'].startswith('map '):value=translate(f['rule'].split()[1],value)
  else:assert f['rule']=='copy'
  put(record,f['scorer_field'],value)
 if 'anim_informed' in truth:assert type(truth['anim_informed']) is bool
 return record
