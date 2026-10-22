"""Full saved-evidence v2 and exact P1/format/scorer gates, no new measurements."""
import argparse,datetime,json,os,pathlib,shutil,subprocess,tarfile,time,traceback,csv,io
import evaluation_common as common
cont=common.CONT;run=common.RUN;sha=common.sha
ap=argparse.ArgumentParser();ap.add_argument('--task',type=int,required=True);a=ap.parse_args();binding=json.loads((cont/'inputs/binding.json').read_text());manifest=json.loads((cont/'inputs/final-v2-context-manifest.json').read_text());selected=manifest[a.task*20:(a.task+1)*20];assert selected;catalogue={(r['query']['sha256'],r['reference']['id']):r for r in map(json.loads,(cont/'inputs/final-anim-catalogue.jsonl').read_text().splitlines())};inventory=json.loads((run/'inputs/original-primary-inventory.json').read_text())['references'];truth_native={r['coded_id']:r for r in json.loads((cont/'inputs/native-scoring-contexts.json').read_text())};truth_l2=json.loads((cont/'inputs/layer2-scoring-contexts.json').read_text());out=cont/'outputs/v2';out.mkdir(exist_ok=True);work=out/f'work/task-{a.task:04}';work.mkdir(parents=True,exist_ok=False);results=[];commands=[];inputs={};failure=None;cache={}
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def command(argv,stem,allowed=(0,)):
 receipt=dict(argv=argv,started=datetime.datetime.now(datetime.timezone.utc).isoformat());start=time.monotonic()
 with pathlib.Path(str(stem)+'.stdout.json').open('x') as stdout,pathlib.Path(str(stem)+'.stderr').open('x') as stderr:done=subprocess.run(['/usr/bin/time','-v','-o',str(stem)+'.resources.txt',*argv],stdout=stdout,stderr=stderr)
 receipt.update(exit_status=done.returncode,seconds=time.monotonic()-start);write(pathlib.Path(str(stem)+'.command.json'),receipt);commands.append(receipt);assert done.returncode in allowed, 'R2 command failure';return json.loads(pathlib.Path(str(stem)+'.stdout.json').read_text())
def archived(r):
 path=r['source_archive']
 if path not in cache:
  assert sha(path)==r['source_archive_sha256'];inputs[path]=r['source_archive_sha256'];data={}
  with tarfile.open(path) as t:
   member=next(m.name for m in t.getmembers() if m.name.endswith('/artifact-hashes.json'));prefix=member.rsplit('/',1)[0];hashes=json.load(t.extractfile(member))
   for name,digest in hashes.items():
    raw=t.extractfile(prefix+'/'+name).read();assert __import__('hashlib').sha256(raw).hexdigest()==digest;data[prefix+'/'+name]=raw
  cache[path]=data
 return cache[path]
def formats(v2,reports):
 assert json.loads((reports/'report.json').read_text())==v2;rows=list(csv.DictReader((reports/'calls.tsv').open(),delimiter='\t'));assert len(rows)==1;r=rows[0]
 for field,value in dict(call_status=v2['status'],species_id=v2['species_id'] or '',subspecies_id=v2['subspecies_id'] or '',call_basis=v2['call_basis'] or '',policy_mode=v2['policy_mode'],validation_status=v2['validation_status'],competitor_evaluation=v2['competitor_evaluation']).items():assert r[field]==value,(field,r[field],value)
 assert set(filter(None,r['reasons'].split(';')))==set(v2['reasons']);md=(reports/'report.md').read_text();final=md.split('### Final call\n',1)[1].split('### ',1)[0];assert '**Call status:** '+v2['status'] in final
 for rank,id in [('Species',v2['species_id']),('Subspecies',v2['subspecies_id'])]:
  if id:assert '- **'+rank+':** assigned' in final and '('+id+')' in final
  else:assert '- **'+rank+':** unresolved (no assignment)' in final
 for reason in v2['reasons']:assert reason in final
 return dict(json_equal=True,tsv_equal=True,markdown_final_call_equal=True)
try:
 assert int(os.environ['SLURM_CPUS_PER_TASK'])>=2;assert json.loads((cont/'checks/final-anim-summary.json').read_text())['gate_pass']
 for p,d in binding['binaries'].items():assert sha(p)==d
 for r in selected:
  if (out/'gate-stop.json').exists():raise RuntimeError('Another full R2/model/scorer/format gate stopped; pending retained')
  case=work/r['record_id'].replace('|','--');case.mkdir();source=archived(r);old_v1=json.loads(source[r['production_plan_member'].rsplit('/',1)[0]+'/v1.stdout.json']);evidence=case/'production-evidence.json'
  if r['dataset']=='native':data=source[r['production_evidence_member']]
  else:data=pathlib.Path(r['production_evidence_path']).read_bytes()
  assert __import__('hashlib').sha256(data).hexdigest()==r['production_evidence_sha256'];evidence.write_bytes(data);evi=json.loads(data);assert evi['provenance']['query_sha256']==r['query_sha256'];assert evi['provenance']['reference_package_sha256']==binding['reference_package']['manifest_sha256'];mask=case/'reference-mask.txt';mask.write_text(''.join(x+'\n' for x in r['mask']));records=[]
  for request in r['anim_scope']:
   record=catalogue[r['query_sha256'],request['reference_id']];assert record['reference']['sha256']==request['reference_sha256'] and request['reference_id'] not in r['mask'];records.append(record)
  assert len(records)==len(r['anim_scope']);recordfile=case/'requested-pairs.jsonl';recordfile.write_text(''.join(json.dumps(x,separators=(',',':'))+'\n' for x in records));base=['--reference-package',binding['reference_package']['path'],'--expected-sha256',binding['reference_package']['manifest_sha256']];extras=['--exclude-references',str(mask)] if r['mask'] else [];dev=str(run/'build/bin/terobacter-dev');v1=command([dev,'saved-evidence',str(evidence),'--evidence-sha256',sha(evidence),*base,*extras],case/'v1')
  canonical=lambda x:dict(status=x['assignment']['status'],species_id=x['assignment']['species_id'],subspecies_id=x['assignment']['subspecies_id'],reasons=sorted(x['blocking_reasons']))
  assert canonical(v1)==canonical(old_v1),'Standalone v1 differs from original completed regression'
  if records:extras+=['--anim-records',str(recordfile),'--anim-method-profile',binding['anim_profile']]
  reports=case/'v2-reports';v2=command([dev,'saved-evidence',str(evidence),'--evidence-sha256',sha(evidence),'--policy','conservative-assembly-anim-v2',*base,*extras,'--report-dir',str(reports)],case/'v2');config=common.p1_input(evi,inventory,records,r['mask']);expected=common.P1.evaluate(config);diff=common.compare_p1(v2,expected);write(case/'p1-comparison.json',dict(exact_match=not diff,differences=diff,expected=expected,input=config));assert not diff,'Full frozen P1 model mismatch'
  assert v2['input']['query_sha256']==r['query_sha256'] and v2['input']['reference_package_sha256']==binding['reference_package']['manifest_sha256'];assert set(v2['input']['mask']['excluded_reference_ids'])==set(r['mask'])
  expected_counts=dict(required=len(records),attempted=len(records),completed=len(records),failed=0,missing_orientation=0,not_run=0,supplied_outside_plan=0,no_alignment=sum(any(x['orientations'][o]['status']=='no_alignment' for o in ['QR','RQ']) for x in records))
  assert all(v2['counts'][k]==v for k,v in expected_counts.items()), 'Pair attempt accounting differs from bound reciprocal records'
  panel_challenges={x for parent in v2['panel_summary'] for x in parent['challenging_reference_ids']};assert panel_challenges==set(v2['challenges']), 'Whole-panel summary omits known challenge'
  format_gate=formats(v2,reports)
  scorer_rows=[];contexts=truth_native[r['coded_id']]['scoring_contexts'] if r['dataset']=='native' else [truth_l2[r['record_id']]]
  for context in contexts:
   record=common.rescue_record(v1,v2,r['initial_trigger'],context['truth'],r['record_id'],context['dataset'],r['setting']);category,violations=common.SCORER.score_rescue(record);scorer_rows.append(dict(category=category,violations=violations,record=record,bookkeeping=context));assert not violations,'Frozen rescue scorer specification violation'
  write(case/'scorer-records.json',scorer_rows);write(case/'format-check.json',format_gate);results.append(dict(record_id=r['record_id'],coded_id=r['coded_id'],query_sha256=r['query_sha256'],dataset=r['dataset'],setting=r['setting'],manifest_state=r['initial_trigger']['state'],rescued=v2['call_basis']=='ANIm_assisted',required_pairs=len(r['anim_scope']),model_pass=True,format_pass=True,scores=[dict(dataset=x['record']['dataset'],category=x['category']) for x in scorer_rows]))
except Exception as e:
 failure=f'{type(e).__name__}: {e}';write(work/'traceback.json',dict(traceback=traceback.format_exc()))
 try:
  with (out/'gate-stop.json').open('x') as h:json.dump(dict(task=a.task,failure=failure,completed=len(results),planned=len(selected)),h,indent=2)
 except FileExistsError:pass
finally:
 write(work/'stage.json',dict(task=a.task,gate_pass=failure is None,failure=failure,selected_records=len(selected),completed_records=len(results),pending_records=len(selected)-len(results),results=results,commands=commands,input_sha256=inputs,binding_sha256=sha(cont/'inputs/binding.json'),manifest_sha256=sha(cont/'inputs/final-v2-context-manifest.json'),source_sha256=sha(__file__),adapter_sha256=sha(cont/'scripts/evaluation_common.py'),job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),new_measurements=False));hashes={str(p.relative_to(work)):sha(p) for p in work.rglob('*') if p.is_file()};write(work/'artifact-hashes.json',hashes);archive=out/f'task-{a.task:04}.tar.gz';part=out/f'.task-{a.task:04}.part';assert not archive.exists() and not part.exists()
 with tarfile.open(part,'w:gz') as t:t.add(work,arcname=f'task-{a.task:04}')
 with tarfile.open(part) as t:
  for name,digest in hashes.items():assert __import__('hashlib').sha256(t.extractfile(f'task-{a.task:04}/'+name).read()).hexdigest()==digest
 part.rename(archive);shutil.rmtree(work);write(out/f'task-{a.task:04}.json',dict(gate_pass=failure is None,failure=failure,selected_records=len(selected),completed_records=len(results),archive_sha256=sha(archive)))
print(json.dumps(dict(task=a.task,completed=len(results),failure=failure)))
if failure:raise SystemExit(1)
