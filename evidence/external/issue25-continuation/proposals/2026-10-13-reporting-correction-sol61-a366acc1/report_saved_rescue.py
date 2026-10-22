"""R25-1: read-only attribution of authenticated, already-scored native outputs."""
import argparse, collections, datetime, hashlib, json, os, pathlib, tarfile, traceback
from fractions import Fraction as F
ORIGINS=('unreturned','returned_below_band','displayed_top_five','returned_outside_top_five')
G=pathlib.Path('${EVALUATION_ROOT}/issue25-storage-audit-sol61-9c546e16/continuation-staged')
R=pathlib.Path('${PROJECT_ROOT}/results/tool-planning/2026-10-11/issue25-frozen-evaluation-sol61-ab252002')
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def rawsha(b):return hashlib.sha256(b).hexdigest()
def write(p,x):
 with pathlib.Path(p).open('x') as h:json.dump(x,h,indent=2);h.write('\n')
def parent(r):return ('species',r['species_id']) if r.get('species_id') else ('lineage',r.get('lineage_id') or r['reference_id'])
def sk_gate(r,t):return F(str(r['ani_percent']))>=t and min(F(str(r['aligned_fraction_query_percent'])),F(str(r['aligned_fraction_reference_percent'])))>=80
def orientation(d,t,strict=False,challenge=False):
 if d is None or d['status']!='aligned':return False
 m=d['metrics'];assert all(type(m[k]) is int for k in ['D','E','query_covered','query_length','reference_covered','reference_length']);assert m['D']>0 and m['query_length']>0 and m['reference_length']>0
 ani=F(100*(m['D']-2*m['E']),m['D']);ar=F(100*m['reference_covered'],m['reference_length']);aq=F(100*m['query_covered'],m['query_length'])
 return (ani>t if strict else ani>=t) and ar>=80 and (challenge or aq>=80)
def finished(pair):return pair is not None and all(pair['orientations'].get(o,{}).get('status') in ['aligned','no_alignment'] for o in ['QR','RQ'])
def strict(pair,t):return pair is not None and all(orientation(pair['orientations'].get(o),t,strict=True) for o in ['QR','RQ'])
def origin(rid,hits,display,lower):
 hit=hits.get(rid);flags={'returned_by_skani':hit is not None,'displayed_top_five_including_ties':rid in display,'below_trigger_lower_band':hit is not None and F(str(hit['ani_percent']))<lower}
 category='unreturned' if hit is None else 'returned_below_band' if flags['below_trigger_lower_band'] else 'displayed_top_five' if rid in display else 'returned_outside_top_five'
 return {'origin':category,**flags}
def derive(evidence,inventory,mask,records,t,lower,conflict):
 active={r['reference_id']:r for r in inventory if r['reference_id'] not in mask};assert len(active)==len(inventory)-len(mask);assert set(mask)<={r['reference_id'] for r in inventory}
 allhits={r['reference_id']:r for r in evidence['ranked_reference_evidence']};assert len(allhits)==len(evidence['ranked_reference_evidence']);assert set(allhits)<={r['reference_id'] for r in inventory}
 hits={k:v for k,v in allhits.items() if k in active};ordered=sorted(hits.values(),key=lambda r:(-F(str(r['ani_percent'])),r['reference_id']));nearest=ordered[0] if ordered else None;band=bool(nearest and lower<=F(str(nearest['ani_percent']))<t)
 display={r['reference_id'] for r in ordered if len(ordered)<=5 or F(str(r['ani_percent']))>=F(str(ordered[4]['ani_percent']))}
 pairs={r['reference']['id']:r for r in records};assert len(pairs)==len(records) and set(pairs)<=set(active)
 candidates=[r['reference_id'] for r in ordered if parent(active[r['reference_id']])==parent(active[nearest['reference_id']]) and active[r['reference_id']].get('species_id') and lower<=F(str(r['ani_percent']))<t and min(F(str(r['aligned_fraction_query_percent'])),F(str(r['aligned_fraction_reference_percent'])))>=80] if nearest else []
 passed=[rid for rid in candidates if strict(pairs.get(rid),t)];challenges={};diagnostics={}
 if band:
  for rid,r in active.items():
   if parent(r)==parent(active[nearest['reference_id']]):continue
   h=hits.get(rid);sk=bool(h and F(str(h['ani_percent']))>=conflict and F(str(h['aligned_fraction_reference_percent']))>=80);an=[o for o in ['QR','RQ'] if orientation(pairs.get(rid,{}).get('orientations',{}).get(o),conflict,challenge=True)]
   if sk or an:challenges[rid]={'skani_challenges':sk,'anim_challenging_orientations':an,**origin(rid,hits,display,lower)}
 for rid,h in hits.items():
  if rid not in pairs:continue
  p=pairs[rid];sk=sk_gate(h,t);gates={o:orientation(p['orientations'].get(o),t) if p['orientations'].get(o,{}).get('status') in ['aligned','no_alignment'] else None for o in ['QR','RQ']};disagreement=any(v is not None and v!=sk for v in gates.values());expected=band and rid in candidates and strict(p,t)
  diagnostics[rid]={'skani_meets_species_gate':sk,'anim_meets_species_gate':gates,'gate_disagreement':disagreement,'expected_ANI_method_crossing':disagreement and expected,'vetoes':disagreement and not expected}
 return dict(active=active,hits=hits,display=display,band=band,candidates=candidates,passed=passed,challenges=challenges,diagnostics=diagnostics,pairs=pairs)
def check_native(v,d):
 assert type(v['trigger']) is bool and v['trigger']==d['band'];assert {r['reference_id'] for r in v['top_references']}==d['display']
 w={r['reference_id']:r for r in v['witness_assessment']};assert len(w)==len(v['witness_assessment']) and set(w)==set(d['candidates'])
 for rid,r in w.items():assert type(r['strict_pass']) is bool and r['strict_pass']==(rid in d['passed']);assert r['finished']==finished(d['pairs'].get(rid))
 assert set(v['challenges'])==set(d['challenges']) and len(v['challenges'])==len(d['challenges']);challenge={r['reference_id']:r for r in v['challenge_evidence']};assert set(challenge)==set(d['challenges']) and len(challenge)==len(v['challenge_evidence'])
 for rid,r in challenge.items():
  assert r['skani_challenges']==d['challenges'][rid]['skani_challenges'] and r['anim_challenging_orientations']==d['challenges'][rid]['anim_challenging_orientations']
  kind,concept=parent(d['active'][rid]);assert r['parent']['kind']==kind and r['parent']['concept_id']==concept
 diagnostics={r['reference_id']:r for r in v['method_diagnostics']};assert len(diagnostics)==len(v['method_diagnostics']) and set(diagnostics)==set(d['diagnostics'])
 for rid,r in diagnostics.items():
  for k,want in d['diagnostics'][rid].items():assert r[k]==want,(rid,k,r[k],want)
 assert bool(any(r['vetoes'] for r in diagnostics.values()))==('method_gate_disagreement' in v['reasons'])
 assert {rid for r in v['panel_summary'] for rid in r['challenging_reference_ids']}==set(d['challenges'])
 return diagnostics,challenge
class Reader:
 def __init__(self,path,expected,pins):
  self.path=pathlib.Path(path);assert sha(self.path)==expected;self.expected=expected;pins[str(self.path)]=expected;self.data={};hash_manifests=[]
  suffixes=('/artifact-hashes.json','/v2.stdout.json','/v2.command.json','/production-evidence.json','/requested-pairs.jsonl','/p1-comparison.json','/scorer-records.json')
  with tarfile.open(self.path,mode='r|gz') as t:
   for m in t:
    if m.isfile() and m.name.endswith(suffixes):
     assert m.name not in self.data;self.data[m.name]=t.extractfile(m).read()
     if m.name.endswith('/artifact-hashes.json'):hash_manifests.append(m.name)
  assert len(hash_manifests)==1;self.prefix=hash_manifests[0].rsplit('/',1)[0];self.hashes=json.loads(self.data[hash_manifests[0]]);self.member_pins={}
 def raw(self,name):
  assert name.startswith(self.prefix+'/');relative=name[len(self.prefix)+1:];want=self.hashes[relative];raw=self.data[name];assert rawsha(raw)==want;self.member_pins[name]=want;return raw
 def json(self,name):return json.loads(self.raw(name))
 def close(self):assert sha(self.path)==self.expected;self.data.clear()

def execute(packet):
 binding=json.load(open(packet/'binding.json'));assert sha(__file__)==binding['source_sha256'][str(packet/'report_saved_rescue.py')]
 for p,w in binding['source_sha256'].items():assert sha(p)==w,('pin',p)
 assert not (packet/'outputs').exists();out=packet/'outputs';out.mkdir();pins={};memberpins={};contextrows=[];challenge_rows=[];method_rows=[];veto_contexts=[];seen=set();failure=None
 manifest={r['record_id']:r for r in json.load(open(G/'inputs/final-v2-context-manifest.json'))};index=json.load(open(G/'inputs/full-v2-result-index.json'));assert set(index)==set(manifest) and len(manifest)==15109
 inventory=json.load(open(R/'inputs/original-primary-inventory.json'))['references'];assert len(inventory)==140
 scored=json.load(open(G/'outputs/scores/rescue-records.json'));by_id=collections.defaultdict(list)
 for r in scored:
  key=(r['dataset'],r['setting'],r['record_id']);assert key not in seen;seen.add(key);by_id[r['record_id']].append(r)
 assert len(scored)==15220 and set(by_id)==set(manifest)
 catalogue={}
 for raw in open(G/'inputs/final-anim-catalogue.jsonl'):
  r=json.loads(raw);key=(r['query']['sha256'],r['reference']['id']);assert key not in catalogue;catalogue[key]=(rawsha(raw.rstrip('\n').encode()),r)
 archives=collections.defaultdict(list)
 for rid,x in index.items():archives[x['archive']].append((rid,x))
 try:
  for number,(path,items) in enumerate(sorted(archives.items())):
   wants={x['archive_sha256'] for rid,x in items};assert len(wants)==1;reader=Reader(path,wants.pop(),pins)
   for rid,x in items:
    m=manifest[rid];assert x['query_sha256']==m['query_sha256'] and x['dataset']==m['dataset'] and x['setting']==m['setting'];case=x['v2_member'].rsplit('/',1)[0];v=reader.json(x['v2_member']);e_raw=reader.raw(case+'/production-evidence.json');assert rawsha(e_raw)==m['production_evidence_sha256'];e=json.loads(e_raw)
    assert v['input']['query_sha256']==e['provenance']['query_sha256']==m['query_sha256'];assert v['input']['evidence_sha256']==rawsha(e_raw);assert set(v['input']['mask']['excluded_reference_ids'])==set(m['mask']);assert v['input']['verification']=='input_and_artifact_hashes_checked; reference_package_verified';assert v['reference_package']['package']['manifest_sha256']==binding['package_manifest_sha256'];assert v['policy_id']=='conservative-assembly-anim-v2'
    assert reader.json(case+'/v2.command.json')['exit_status']==0;assert reader.json(x['p1_member'])['exact_match'] is True
    pairs_raw=reader.raw(case+'/requested-pairs.jsonl');pairs=[json.loads(line) for line in pairs_raw.splitlines()];scope={r['reference_id']:r['reference_sha256'] for r in m['anim_scope']};assert len(scope)==len(m['anim_scope']);assert len(pairs)==len(scope);assert set(scope)=={r['reference']['id'] for r in pairs}==set(v['input']['anim']['bound_pairs'])
    for r in pairs:
     assert r['query']['sha256']==m['query_sha256'] and r['reference']['sha256']==scope[r['reference']['id']];retained=catalogue[m['query_sha256'],r['reference']['id']][1];assert r==retained
    if pairs:assert v['input']['anim_record_files'][0]['sha256']==rawsha(pairs_raw) and v['input']['anim']['method_profile']==binding['anim_profile']
    t=F(v['thresholds']['species_ani']);lower=F(v['thresholds']['trigger_lower_ani']);conflict=F(v['thresholds']['conflict_ani']);assert (t,lower,conflict)==(F(96),F('94.5'),F(96));d=derive(e,inventory,set(m['mask']),pairs,t,lower,conflict);diagnostics,challenges=check_native(v,d)
    assert len(d['active'])==v['input']['mask']['assessed_panel_reference_count'];assert v['counts']['required']==len(scope) and v['counts']['failed']==0
    actual_scorer=reader.json(x['scorer_member']);assert len(actual_scorer)==len(by_id[rid])
    actual_by_key={(r['record']['dataset'],r['record']['setting'],r['record']['record_id']):r for r in actual_scorer};assert len(actual_by_key)==len(actual_scorer)
    provenance={'archive':path,'archive_sha256':x['archive_sha256'],'v2_member':x['v2_member'],'v2_sha256':reader.member_pins[x['v2_member']],'primary_evidence_member':case+'/production-evidence.json','primary_evidence_sha256':rawsha(e_raw),'original_primary_evidence_sha256':m['production_evidence_sha256'],'original_source_archive':m['source_archive'],'original_source_archive_sha256':m['source_archive_sha256'],'scope_pairs_member':case+'/requested-pairs.jsonl','scope_pairs_sha256':rawsha(pairs_raw),'mask_reference_ids':m['mask'],'reference_inventory_sha256':binding['source_sha256'][str(R/'inputs/original-primary-inventory.json')],'resolution':'actual final result index; per-member authenticated archived copy of original full primary evidence, irrespective of deleted loose paths or prior failed aggregate stages'}
    method_vetos=[r for r in diagnostics.values() if r['vetoes']];crossings=[r for r in diagnostics.values() if r['expected_ANI_method_crossing']];incomplete=bool(v['evaluation_status']=='incomplete' or v['rescue_status']=='incomplete' or any('incomplete' in x for x in v['reasons']))
    for s in by_id[rid]:
     key=(s['dataset'],s['setting'],rid);assert actual_by_key[key]['record']==s['source_v1_v2_record'] and actual_by_key[key]['bookkeeping']==s['bookkeeping'];assert s['source_v1_v2_record']['v2']['trigger']==v['trigger'];base={'record_id':rid,'coded_id':m['coded_id'],'query_sha256':m['query_sha256'],'dataset':s['dataset'],'setting':s['setting'],'category':s['category'],'applicable':s['applicable'],'block':s['block'],'bookkeeping':s['bookkeeping'],'source':provenance};counts={'triggered':int(v['trigger']),'manifest_witness_candidate_references':int(m['initial_trigger']['witness_candidates']),'reported_witness_candidate_references':len(d['candidates']),'triggered_witness_candidate_references':len(d['candidates']) if v['trigger'] else 0,'passing_witness_references':len(d['passed']),'contexts_with_passing_witness':int(bool(d['passed'])),'foreign_challenge_references':len(challenges),'contexts_with_foreign_challenge':int(bool(challenges)),'method_veto_references':len(method_vetos),'contexts_with_method_veto':int(bool(method_vetos)),'expected_crossing_references':len(crossings),'contexts_with_expected_crossing':int(bool(crossings)),'incomplete':int(incomplete),'rescued':int(v['call_basis']=='ANIm_assisted'),'complete_panel':int(v['competitor_evaluation']=='complete')};contextrows.append({**base,**counts,'native_reasons':v['reasons'],'candidate_reference_ids':d['candidates'],'passing_witness_reference_ids':d['passed'],'challenge_reference_ids':sorted(challenges),'retained_method_veto_reference_ids':[r['reference_id'] for r in method_vetos],'thresholds':v['thresholds'],'member_compatibility_not_assessed_by_this_report':True})
     for ref,entry in challenges.items():challenge_rows.append({**base,**d['challenges'][ref],'reference_id':ref,'actual_challenge_evidence':entry,'complete_panel':v['competitor_evaluation']=='complete','skani_primary_evidence':d['hits'].get(ref),'reciprocal_metrics':{o:d['pairs'][ref]['orientations'][o] for o in ['QR','RQ']} if ref in d['pairs'] else None,'trigger_lower_ani':str(lower)})
     for ref,entry in diagnostics.items():
      if entry['expected_ANI_method_crossing'] or entry['vetoes']:method_rows.append({**base,'reference_id':ref,'native_diagnostic':entry,'independent_diagnostic':d['diagnostics'][ref],'reciprocal_metrics':{o:d['pairs'][ref]['orientations'][o] for o in ['QR','RQ']},'native_reasons':v['reasons'],'competitor_evaluation':v['competitor_evaluation'],'thresholds':v['thresholds']})
     if method_vetos:veto_contexts.append({**base,'retained_method_veto_reference_ids':[r['reference_id'] for r in method_vetos],'native_reasons':v['reasons'],'rescued':v['call_basis']=='ANIm_assisted'})
   memberpins[path]=reader.member_pins;reader.close()
   if number%50==0:print(json.dumps({'archives_read':number+1,'contexts':len(contextrows)}),flush=True)
   assert sum(p.stat().st_size for p in packet.rglob('*') if p.is_file())<binding['output_cap_bytes']
  assert len(contextrows)==15220;groups=collections.defaultdict(list)
  for r in contextrows:groups[r['dataset'],r['setting'],r['applicable']].append(r)
  summaries=[];fields=list(counts)
  for (dataset,setting,applicable),rows in sorted(groups.items()):
   keys={r['record_id'] for r in rows};cr=[r for r in challenge_rows if r['dataset']==dataset and r['setting']==setting and r['record_id'] in keys];mr=[r for r in method_rows if r['dataset']==dataset and r['setting']==setting and r['record_id'] in keys];summaries.append({'dataset':dataset,'setting':setting,'applicable':applicable,'records':len(rows),'distinct_query_sha256':len({r['query_sha256'] for r in rows}),'counts':{k:sum(r[k] for r in rows) for k in fields},'scorer_categories':dict(collections.Counter(r['category'] for r in rows)),'complete_panel_challenge_origins':{o:sum(r['complete_panel'] and r['origin']==o for r in cr) for o in ORIGINS},'other_panel_challenge_origins_descriptive':{o:sum(not r['complete_panel'] and r['origin']==o for r in cr) for o in ORIGINS},'actual_reference_diagnostics':{'expected_ANI_method_crossing':sum(r['native_diagnostic']['expected_ANI_method_crossing'] for r in mr),'retained_method_disagreement_veto':sum(r['native_diagnostic']['vetoes'] for r in mr)},'denominator_scope':'record-setting observations; applicable=false excluded from omission arm rates; query/settings and source overlap retained, no independent-accuracy claim'})
  for name,rows in [('RS1-context-counts.jsonl',contextrows),('RS8-challenge-origin-catalogue.jsonl',challenge_rows),('RS9-method-crossing-veto-catalogue.jsonl',method_rows),('RS9-retained-veto-contexts.jsonl',veto_contexts)]:
   with (out/name).open('x') as h:
    for r in rows:h.write(json.dumps(r,separators=(',',':'))+'\n')
  write(out/'dataset-setting-summary.json',summaries);write(out/'source-archive-member-pins.json',memberpins)
  complete=[r for r in contextrows if r['complete_panel']];report={'scored_record_settings':len(contextrows),'unique_classification_contexts':len(manifest),'distinct_query_sha256':len({r['query_sha256'] for r in contextrows}),'archive_count':len(pins),'counts':{k:sum(r[k] for r in contextrows) for k in fields},'complete_panel_contexts':len(complete),'complete_panel_reference_denominator':sum(len(manifest[r['record_id']]['anim_scope']) for r in complete),'complete_panel_challenge_rows':sum(r['complete_panel'] for r in challenge_rows),'complete_panel_challenge_origins':{o:sum(r['complete_panel'] and r['origin']==o for r in challenge_rows) for o in ORIGINS},'other_panel_challenge_origins_descriptive':{o:sum(not r['complete_panel'] and r['origin']==o for r in challenge_rows) for o in ORIGINS},'expected_crossing_rows':sum(r['native_diagnostic']['expected_ANI_method_crossing'] for r in method_rows),'retained_method_veto_rows':sum(r['native_diagnostic']['vetoes'] for r in method_rows),'retained_method_veto_contexts':len(veto_contexts),'actual_counts_not_rates':True,'overlapping_observations_not_independent':True};write(out/'report-summary.json',report)
  for p,w in binding['source_sha256'].items():assert sha(p)==w,('post_read_pin',p)
  size=sum(p.stat().st_size for p in packet.rglob('*') if p.is_file());assert size<binding['output_cap_bytes'];gate={'gate_pass':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_measurements':False,'new_classifications_or_scores':False,'source_binding_sha256':sha(packet/'binding.json'),'job_id':os.environ.get('SLURM_JOB_ID'),'allocated_cpus':os.environ.get('SLURM_CPUS_PER_TASK'),'source_archive_sha256':pins,'output_sha256':{str(p):sha(p) for p in out.iterdir()},'retained_bytes':size,'full_primary_evidence_independent_attribution':True,'all_scored_record_settings_accounted':len(contextrows)==15220};write(packet/'reporting-gate.json',gate);print(json.dumps(report),flush=True)
 except Exception as e:
  write(packet/'reporting-failure.json',{'gate_pass':False,'failure':str(e),'traceback':traceback.format_exc(),'completed_contexts':len(contextrows),'source_archive_sha256':pins});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--packet',type=pathlib.Path,required=True);a=p.parse_args();execute(a.packet)
