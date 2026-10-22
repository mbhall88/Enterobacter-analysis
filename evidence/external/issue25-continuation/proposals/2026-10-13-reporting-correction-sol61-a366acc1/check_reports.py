"""Separate readback of report coverage, zero categories, overlap and cross-catalogue totals."""
import collections,json,pathlib,hashlib,sys,datetime
D=pathlib.Path(sys.argv[1]);G=D.parents[1]
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def rows(name):return [json.loads(line) for line in open(D/'outputs'/name)]
checks=[]
def check(name,predicate):assert predicate,name;checks.append(name)
gate=json.load(open(D/'reporting-gate.json'));b=json.load(open(D/'binding.json'));check('report extraction pass',gate['gate_pass'] and gate['new_measurements'] is False and gate['new_classifications_or_scores'] is False);check('binding same',gate['source_binding_sha256']==sha(D/'binding.json'));check('allocated >=2 CPUs',int(gate['allocated_cpus'])>=2)
for p,w in gate['output_sha256'].items():check('output '+pathlib.Path(p).name,sha(p)==w)
for p,w in b['source_sha256'].items():check('source preserved '+pathlib.Path(p).name,sha(p)==w)
c=rows('RS1-context-counts.jsonl');ch=rows('RS8-challenge-origin-catalogue.jsonl');m=rows('RS9-method-crossing-veto-catalogue.jsonl');v=rows('RS9-retained-veto-contexts.jsonl');s=json.load(open(G/'outputs/scores/rescue-records.json'));summaries=json.load(open(D/'outputs/dataset-setting-summary.json'));tot=json.load(open(D/'outputs/report-summary.json'))
key=lambda r:(r['dataset'],r['setting'],r['record_id']);ck={key(r):r for r in c};sk={key(r):r for r in s};check('exact 15220 scored identities',len(ck)==len(c)==len(sk)==15220 and set(ck)==set(sk));check('exact15109 original classification contexts',len({r['record_id'] for r in c})==15109)
for k,r in ck.items():
 original=sk[k];assert all(r[f]==original[f] for f in ['dataset','setting','record_id','category','applicable','block','bookkeeping']);assert r['passing_witness_references']==len(r['passing_witness_reference_ids']);assert r['foreign_challenge_references']==len(r['challenge_reference_ids']);assert r['method_veto_references']==len(r['retained_method_veto_reference_ids']);assert r['contexts_with_passing_witness']==bool(r['passing_witness_reference_ids']);assert r['contexts_with_foreign_challenge']==bool(r['challenge_reference_ids']);assert r['contexts_with_method_veto']==bool(r['retained_method_veto_reference_ids']);assert r['complete_panel'] in [0,1]
checks.append('all individual identities/bookkeeping/count vectors')
for catalogue in [ch,m]:
 identities=[(*key(r),r['reference_id']) for r in catalogue];check('no duplicate reference/context identities',len(set(identities))==len(identities));check('all catalogue rows match context sources',all(key(r) in ck and r['source']==ck[key(r)]['source'] for r in catalogue))
origins=['unreturned','returned_below_band','displayed_top_five','returned_outside_top_five'];check('explicit all4 origin zeros/values',list(tot['complete_panel_challenge_origins'])==origins and list(tot['other_panel_challenge_origins_descriptive'])==origins)
for r in ch:
 expected='unreturned' if not r['returned_by_skani'] else 'returned_below_band' if r['below_trigger_lower_band'] else 'displayed_top_five' if r['displayed_top_five_including_ties'] else 'returned_outside_top_five';assert r['origin']==expected;assert r['reference_id'] in ck[key(r)]['challenge_reference_ids']
checks.append('challenge attribution precedence/denominator')
check('all reference challenge counts reconstructed',sum(r['foreign_challenge_references'] for r in c)==len(ch));check('all method veto reference counts reconstructed',sum(r['method_veto_references'] for r in c)==sum(r['native_diagnostic']['vetoes'] for r in m));check('all crossings reconstructed',sum(r['expected_crossing_references'] for r in c)==sum(r['native_diagnostic']['expected_ANI_method_crossing'] for r in m));check('exact veto context catalogue',{key(r) for r in v}=={key(r) for r in c if r['contexts_with_method_veto']} and len(v)==len({key(r) for r in v}));check('retained veto never rescue',all(not r['rescued'] and 'method_gate_disagreement' in r['native_reasons'] for r in v))
for r in m:assert r['native_diagnostic']['expected_ANI_method_crossing']!=r['native_diagnostic']['vetoes'];assert all(r['native_diagnostic'][k]==want for k,want in r['independent_diagnostic'].items())
checks.append('crossing exception/veto mutually exclusive')
fields=list(tot['counts']);group=collections.defaultdict(list)
for r in c:group[r['dataset'],r['setting'],r['applicable']].append(r)
check('all dataset setting applicability summaries',set(group)=={(r['dataset'],r['setting'],r['applicable']) for r in summaries})
for r in summaries:
 x=group[r['dataset'],r['setting'],r['applicable']];assert r['records']==len(x);assert r['counts']=={f:sum(z[f] for z in x) for f in fields};assert set(r['complete_panel_challenge_origins'])==set(origins);assert set(r['other_panel_challenge_origins_descriptive'])==set(origins)
checks.append('all 20 grouped totals and explicit zeros');check('aggregate reconstructed',tot['counts']=={f:sum(r[f] for r in c) for f in fields});check('complete panel reference denominator',tot['complete_panel_reference_denominator']==sum(140-len(r['source']['mask_reference_ids']) for r in c if r['complete_panel']))
overlap=collections.defaultdict(list)
for r in c:overlap[r['query_sha256']].append({'record_id':r['record_id'],'coded_id':r['coded_id'],'dataset':r['dataset'],'setting':r['setting'],'applicable':r['applicable'],'block':r['block']})
physical_witness={(r['query_sha256'],ref) for r in c for ref in r['passing_witness_reference_ids']};physical_challenge={(r['query_sha256'],ref) for r in c for ref in r['challenge_reference_ids']};physical_method={(r['query_sha256'],r['reference_id']) for r in m if r['native_diagnostic']['vetoes']};physical_crossing={(r['query_sha256'],r['reference_id']) for r in m if r['native_diagnostic']['expected_ANI_method_crossing']}
with (D/'outputs/overlap-bookkeeping.json').open('x') as h:json.dump({'distinct_query_sha256':len(overlap),'multiple_record_setting_queries':sum(len(r)>1 for r in overlap.values()),'query_contexts':overlap,'unique_passing_witness_pairs_descriptive':len(physical_witness),'unique_challenge_pairs_descriptive':len(physical_challenge),'unique_veto_pairs_descriptive':len(physical_method),'unique_expected_crossing_pairs_descriptive':len(physical_crossing),'not_accuracy_denominators':True},h,separators=(',',':'));h.write('\n')
check('within 512MiB retained cap',sum(p.stat().st_size for p in D.rglob('*') if p.is_file())<b['output_cap_bytes']);receipt={'gate_pass':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'check_count':len(checks),'reporting_gate_sha256':sha(D/'reporting-gate.json'),'source_sha256':sha(__file__),'overlap_sha256':sha(D/'outputs/overlap-bookkeeping.json'),'counts':tot};
with (D/'readback.json').open('x') as h:json.dump(receipt,h,indent=2);h.write('\n')
print(json.dumps({'gate_pass':True,'check_count':len(checks),'counts':tot,'unique_pairs':{'passing_witness':len(physical_witness),'foreign_challenge':len(physical_challenge),'retained_veto':len(physical_method),'expected_crossing':len(physical_crossing)}}))
