"""Join only fully audited ANIm attempts; full mode requires all60 waves."""
import pathlib,json,hashlib,datetime,argparse
ROOT=pathlib.Path(__file__).resolve().parents[1];RUN=pathlib.Path('${PROJECT_ROOT}/results/tool-planning/2026-10-11/issue25-frozen-evaluation-sol61-ab252002');OLD=RUN/'continuations/2026-10-12-p4b-resume-sol61-013515cc';DURABLE=pathlib.Path('${EVALUATION_ROOT}/issue25-owned-quota-preservation-sol61-8b9385c9');PROFILE='reciprocal-mummer4-4.0.1-preserve-n-mum-1-span-v1'
ap=argparse.ArgumentParser();ap.add_argument('--verify-prefix-only',action='store_true');args=ap.parse_args()
def sha(p):
 with pathlib.Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def integers(r):return {'status':r['status'],'query':r['query']['sha256'],'reference':r['reference']['sha256'],'orientations':{k:{'status':v['status'],'metrics':v.get('metrics')} for k,v in r['orientations'].items()}}
sourcepins={};proofs=[];records={};aliases={};duplicates=[]
def ingest(path,label,proof=None):
 sourcepins[str(path)]=sha(path)
 if proof:
  gate=json.loads(proof.read_text());assert gate['gate_pass'] and gate['failure'] is None;sourcepins[str(proof)]=sha(proof);proofs.append({'record_file':str(path),'audit_file':str(proof),'counts':gate['counts'],'archive_sha256':gate['archive_sha256']})
 rows=[json.loads(line) for line in path.open()];assert not proof or len(rows)==gate['counts']['complete_pairs']
 for r in rows:
  assert r['status']=='complete' and r['method']['profile']==PROFILE and r['method']['amendment']=='P1-ANIM-MUMMER4-r1' and r['method']['implementation_sha256']=='84a7658063d4c023fb69e591ec25d5917bf03821d230c21b9a36f83d5e3042e7'
  assert all(r['orientations'][o]['status'] in ['aligned','no_alignment'] for o in ['QR','RQ'])
  key=hashlib.sha256('|'.join([PROFILE,r['query']['sha256'],r['reference']['sha256']]).encode()).hexdigest();identity=(r['query']['sha256'],r['reference']['id']);assert identity not in aliases or aliases[identity]==key;aliases[identity]=key
  if key in records:assert integers(records[key])==integers(r),'Successful attempts have differing integers';duplicates.append({'key':key,'last_successful_source':label})
  records[key]=r
 return len(rows)
resource=OLD/'checks/retained-anim-resource-binding.json';gate=json.loads(resource.read_text());assert gate['gate_pass'];sourcepins[str(resource)]=sha(resource)
assert ingest(OLD/'inputs/retained-anim-catalogue.jsonl','accepted_retained4381')==4381
assert ingest(DURABLE/'anim/records/plan-panel-postseal-wave-0001.jsonl','successful1000_superseding29656250',DURABLE/'checks/anim-integers-plan-panel-postseal-wave-0001.json')==1000
assert ingest(ROOT/'anim/records/plan-published399-wave-0001.jsonl','recovered399',ROOT/'checks/anim-integers-plan-published399-wave-0001.json')==399
assert ingest(RUN/'anim/records/plan-rs7-wave-0001.jsonl','postseal_RS7_requested33',RUN/'checks/anim-integers-plan-rs7-wave-0001.json')==33
assert ingest(RUN/'anim/records/plan-rs7-recheck-wave-0001.jsonl','postseal_EN4only1',RUN/'checks/anim-integers-plan-rs7-recheck-wave-0001.json')==1
if args.verify_prefix_only:
 assert len(records)==5814 and not duplicates
 receipt={'prefix_binding_pass':True,'retained_unique_pairs':len(records),'normal_candidate_pairs':5780,'postseal_RS7_candidate_pairs':33,'EN4_only_extra_pairs':1,'preseal_keys_ingested':False,'source_sha256':sourcepins,'no_new_measurements_or_classifications':True};(ROOT/'checks/final-anim-prefix-preparation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='source_sha256'}));raise SystemExit(0)
plan=ROOT/'anim/plan-panel-group-r1';fresh=0
for i in range(1,61):
 name=f'wave-{i:04d}';summary_path=plan/'collected'/f'{name}.summary.json';summary=json.loads(summary_path.read_text());n=100 if i<60 else 77;assert summary['planned']==summary['pairs']==n and summary['failed']==summary['missing']==0;sourcepins[str(summary_path)]=sha(summary_path)
 fresh+=ingest(ROOT/'anim/records'/f'{plan.name}-{name}.jsonl',f'new_group_wave_{i}',ROOT/'checks'/f'anim-integers-{plan.name}-{name}.json')
assert fresh==5977 and len(records)==11791 and not duplicates
manifest_path=ROOT/'inputs/final-v2-context-manifest.json';manifest=json.loads(manifest_path.read_text());sourcepins[str(manifest_path)]=sha(manifest_path);requested=set()
for c in manifest:
 for r in c['anim_scope']:
  key=aliases[c['query_sha256'],r['reference_id']];record=records[key];assert record['reference']['sha256']==r['reference_sha256'] and r['reference_id'] not in c['mask'];requested.add(key)
assert len(requested)==11790
out=ROOT/'inputs/final-anim-catalogue.jsonl';assert not out.exists();part=out.with_suffix('.jsonl.part');assert not part.exists()
with part.open('x') as h:
 for key in sorted(records):h.write(json.dumps(records[key],separators=(',',':'))+'\n')
part.rename(out);summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gate_pass':True,'normal_candidate_pairs':11757,'postseal_RS7_candidate_pairs':33,'EN4_only_extra_pairs':1,'unique_retained_pairs':len(records),'candidate_required_keys':len(requested),'fresh_group_pairs':fresh,'duplicate_success_attempts':duplicates,'all_requested_pairs_available_exact_scope':True,'preseal_RS7_keys_ingested':False,'input_sha256':sourcepins,'successful_audit_receipts':proofs,'final_catalogue_sha256':sha(out),'final_context_manifest_sha256':sha(manifest_path),'source_sha256':sha(pathlib.Path(__file__))};(ROOT/'checks/final-anim-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ['input_sha256','successful_audit_receipts']}))
