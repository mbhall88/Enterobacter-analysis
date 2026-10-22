"""Frozenstep9 actual scoring;every planned record retained;stop before rates on specification failure."""
import pathlib,json,datetime,collections,traceback
import evaluation_common as c
root=c.CONT
assert json.loads((root/'checks/step7-summary.json').read_text())['gate_pass']
assert json.loads((root/'checks/genomovar-sets-summary.json').read_text())['gate_pass']
assert json.loads((root/'checks/saved-gv-perturbation-summary.json').read_text())['gate_pass']
assert json.loads((root/'checks/adapter-tests.json').read_text())['gate_pass']
assert json.loads((root/'checks/production-truth-transfer-fixtures.json').read_text())['gate_pass']
rescue=json.loads((root/'inputs/rescue-scorer-inputs.json').read_text());genomovar=json.loads((root/'inputs/genomovar-scorer-inputs.json').read_text())+json.loads((root/'inputs/saved-gv-scorer-inputs.json').read_text());blocks=json.loads((root/'inputs/layer2-summary-block-contexts.json').read_text());manifest={r['record_id']:r for r in json.loads((root/'inputs/final-v2-context-manifest.json').read_text())};assert len(rescue)==15220;out=root/'outputs/scores';out.mkdir();scored_rescue=[];scored_gv=[];seen=set();failure=None
try:
 for row in rescue:
  r=row['record'];key=(r['dataset'],r['setting'],r['record_id']);assert key not in seen;seen.add(key);context=row['bookkeeping'];truth=r['truth']
  if 'anim_informed' in truth:assert type(truth['anim_informed']) is bool
  category,violations=c.SCORER.score_rescue(r);entry={'record_id':r['record_id'],'dataset':r['dataset'],'setting':r['setting'],'manifest_state':r['manifest']['state'],'manifest_parent':r['manifest']['parent_id'],'category':category,'violations':violations,'truth':truth,'bookkeeping':context,'source_v1_v2_record':r};scored_rescue.append(entry);assert category!='spec_violation' and not violations,'Frozenrescue scorer specification violation'
  if r['dataset']=='layer2':
   b=blocks[r['record_id']];assert b['query_sha256']==manifest[r['record_id']]['query_sha256'];entry['block']=b['clone_block'];entry['omitted_species']=truth.get('species_id');entry['applicable']=context['omission_applicable'];assert type(entry['applicable']) is bool
  else:entry['block']=truth.get('block') or truth.get('component_id');entry['applicable']=True
  assert entry['block'],'No recorded sensitivity block'
 for row in genomovar:
  r=row['record'];assert type(r['diagnostic']) is bool and type(r['self_member']) is bool;category,violations=c.SCORER.score_genomovar(r);scored_gv.append({'record_id':r['record_id'],'dataset':r['dataset'],'category':category,'violations':violations,'record':r,'source':row['source']});assert category!='spec_violation' and not violations,'Frozengenomovar scorer specification violation'
except Exception as e:failure={'message':str(e),'traceback':traceback.format_exc()}
finally:
 (out/'rescue-records.json').write_text(json.dumps(scored_rescue,separators=(',',':'))+'\n');(out/'genomovar-records.json').write_text(json.dumps(scored_gv,separators=(',',':'))+'\n');receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':failure is None,'failure':failure,'planned_rescue_records':len(rescue),'scored_rescue_records':len(scored_rescue),'planned_genomovar_records':len(genomovar),'scored_genomovar_records':len(scored_gv),'source_sha256':c.sha(__file__),'scorer_source_sha256':c.SOURCE_PINS[str(c.P4A/'scripts/score_changed_behaviour.py')],'rates_published':False};(out/'scoring-gate.json').write_text(json.dumps(receipt,indent=2)+'\n')
if failure:raise RuntimeError('Scoringstopped;partialrecordsandfailurepreserved;no rates published')
groups=collections.defaultdict(list)
for r in scored_rescue:groups[r['dataset'],r['setting'],r['applicable']].append(r)
summary=[]
for key,rows in sorted(groups.items()):
 if key[0]=='layer2' and not key[2]:
  summary.append({'dataset':key[0],'setting':key[1],'applicable':False,'records':len(rows),'descriptive_category_counts':dict(collections.Counter(r['category'] for r in rows)),'excluded_from_arm_denominators':True});continue
 x={'dataset':key[0],'setting':key[1],'applicable':key[2],'records':len(rows),'record_summary':c.SCORER.summarise_rescue(rows),'clone_or_component_block_summary':c.SCORER.summarise_rescue(c.SCORER.thin_by_block(rows))}
 if key[0]=='layer2':x['species_level']=c.SCORER.by_cluster(rows,'omitted_species');by=collections.defaultdict(list)
 if key[0]=='layer2':
  for row in rows:by[row['omitted_species']].append(row)
  x['per_species']=[{'species_id':species,'records':len(items),'summary':c.SCORER.summarise_rescue(items)} for species,items in sorted(by.items(),key=lambda x:str(x[0]))]
 summary.append(x)
(out/'rescue-summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'genomovar-category-counts.json').write_text(json.dumps({dataset:dict(collections.Counter(x['category'] for x in scored_gv if x['dataset']==dataset)) for dataset in sorted({r['dataset'] for r in scored_gv})},indent=2)+'\n');print(json.dumps(receipt))

gv_safety=[]
for dataset in ['GV4','GV6']:
 rows=[r for r in scored_gv if r['dataset']==dataset];classified=[]
 for r in rows:
  rec=r['record'];assert rec['source_group'];classified.append({'record_id':r['record_id'],'block':rec['source_group'],'wrong':r['category'].split(':')[-1] in ['gv_other_group','gv_other_group_under_omission']})
 k=sum(r['wrong'] for r in classified);thin=c.SCORER.thin_by_block(classified);kb=sum(r['wrong'] for r in thin);gv_safety.append({'dataset':dataset,'records':len(rows),'wrong_group':k,'wrong_group_rate':f'{k}/{len(rows)}','ci95':c.SCORER.clopper_pearson(k,len(rows)),'group_sensitivity_records':len(thin),'group_sensitivity_wrong':kb,'group_sensitivity_ci95':c.SCORER.clopper_pearson(kb,len(thin)),'claim_class':'Safety;definition-relative counts,not independent accuracy'})
(out/'genomovar-safety-summary.json').write_text(json.dumps(gv_safety,indent=2)+'\n');receipt['rates_published']=True;(out/'scoring-gate.json').write_text(json.dumps(receipt,indent=2)+'\n')
