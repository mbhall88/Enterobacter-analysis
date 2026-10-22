"""Frozen EN5/CU2 native integration exception; actual native/split/P1 comparisons."""
import pathlib,json,argparse,datetime,subprocess,tempfile,time,tarfile,hashlib,traceback,os,signal
import evaluation_common as c
from stage_archive import pack
root=c.CONT;run=c.RUN;sha=c.sha
ap=argparse.ArgumentParser();ap.add_argument('--task',type=int,required=True);a=ap.parse_args();manifest=json.loads((root/'inputs/integrated-context-manifest.json').read_text());row=manifest[a.task];parameters=json.loads((root/'inputs/integrated-argv-binding.json').read_text());assert json.loads((root/'checks/full-v2-summary.json').read_text())['gate_pass'];catalogue={(x['query']['sha256'],x['reference']['id']):x for x in map(json.loads,(root/'inputs/final-anim-catalogue.jsonl').read_text().splitlines())};inventory=json.loads((run/'inputs/original-primary-inventory.json').read_text())['references'];binding=json.loads((root/'inputs/binding.json').read_text());out=root/'outputs/integrated';out.mkdir(parents=True,exist_ok=True);assert not(out/f'task-{a.task:04}.tar.gz').exists() and not(out/f'.task-{a.task:04}.part').exists();work=pathlib.Path(tempfile.mkdtemp(prefix=f'p4b-integrated-{a.task:04}-',dir='/tmp'));failure=None;result={};argv=[];started=datetime.datetime.now(datetime.timezone.utc).isoformat()
def interrupted(sig,frame):raise RuntimeError('Owned integrated task interrupted '+str(sig))
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
try:
 assert int(os.environ['SLURM_CPUS_PER_TASK'])>=16;assert __import__('shutil').disk_usage('/tmp').free>=1024**3, 'Insufficient node-local temporary space';probe=work/'local-write-probe';data=b'P4b native private temporary directory\n'*8192
 with probe.open('xb') as h:h.write(data);h.flush();os.fsync(h.fileno())
 assert probe.read_bytes()==data;probe.unlink()
 for path,want in binding['binaries'].items():assert sha(path)==want
 assert sha(row['query_fasta'])==row['query_sha256'];assert sha(row['source_archive'])==row['source_archive_sha256'];assert sha(root/'inputs/final-anim-catalogue.jsonl')==json.loads((root/'checks/final-anim-summary.json').read_text())['final_catalogue_sha256']
 with tarfile.open(row['source_archive']) as t:raw=t.extractfile(row['production_evidence_member']).read();assert hashlib.sha256(raw).hexdigest()==row['production_evidence_sha256'];evidence=json.loads(raw)
 assert evidence['provenance']['query_sha256']==row['query_sha256']
 argv=parameters['immutable_argv_prefix']+['--anim-output',str(work/'native-pairs'),'--anim-cache',str(work/'native-cache'),'--report-dir',str(work/'reports')]
 if row['custom_thresholds']:argv+=['--species-ani',row['custom_thresholds']['species_ani'],'--subspecies-ani',row['custom_thresholds']['subspecies_ani']]
 argv+=[row['query_fasta']];start=time.monotonic()
 with(work/'native.stdout.json').open('x') as stdout,(work/'native.stderr').open('x') as stderr:done=subprocess.run(['/usr/bin/time','-v','-o',str(work/'native.resources.txt'),*argv],stdout=stdout,stderr=stderr)
 (work/'native.command.json').write_text(json.dumps({'argv':argv,'exit_status':done.returncode,'wall_seconds':time.monotonic()-start},indent=2)+'\n');assert done.returncode==0,'Native integrated route failed';native=json.loads((work/'native.stdout.json').read_text());computed=[];integer_checks=[]
 for path in sorted((work/'native-pairs').rglob('result.json')):
  record=json.loads(path.read_text())
  if record.get('schema')!='terobacter-anim-pair-v1':continue
  assert record['status']=='complete' and record['query']['sha256']==row['query_sha256'];computed.append(record);reference=record['reference']['id'];split=catalogue[row['query_sha256'],reference]
  for key in ['profile','amendment','upstream_commit','build_sha256','parser_sha256','source_sha256','query_threads']:assert record['method'][key]==split['method'][key]
  assert record['method']['implementation_sha256']==binding['binaries'][str(run/'build/bin/tero')]
  for orientation in ['QR','RQ']:
   actual=record['orientations'][orientation];expected=split['orientations'][orientation];assert actual['status']==expected['status'] and actual.get('metrics')==expected.get('metrics'),'Integrated/split integer mismatch';integer_checks.append({'reference_id':reference,'orientation':orientation,'exact_match':True})
 if row['custom_thresholds'] is None:assert {x['reference']['id'] for x in computed}=={x['reference_id'] for x in row['anim_scope']},'Integrated requested-set mismatch'
 config=c.p1_input(evidence,inventory,computed,custom=row['custom_thresholds']);expected=c.P1.evaluate(config);differences=c.compare_p1(native,expected);assert not differences,'Integrated full P1 mismatch';assert json.loads((work/'reports/report.json').read_text())==native
 if row['custom_thresholds']:
  assert native['policy_mode']=='custom' and native['validation_status']=='custom_unvalidated'
 index=json.loads((root/'inputs/full-v2-result-index.json').read_text());source=index[row['coded_id'] if row['custom_thresholds'] else row['record_id']];assert sha(source['archive'])==source['archive_sha256']
 with tarfile.open(source['archive']) as t:split_output=json.load(t.extractfile(source['v2_member']))
 fields=['status','species_id','subspecies_id','call_basis','reasons','trigger','competitor_evaluation','challenges','counts','top_references','panel_summary']
 if row['custom_thresholds'] is None:fields+=['policy_mode','validation_status']
 assert {k:native[k] for k in fields}=={k:split_output[k] for k in fields},'Native/split or explicit-default decision mismatch'
 result={'record_id':row['record_id'],'coded_id':row['coded_id'],'kind':row['integration_kind'],'query_sha256':row['query_sha256'],'P1_pass':True,'split_pass':row['custom_thresholds'] is None,'integer_checks':integer_checks,'native_pairs':len(computed),'custom_unvalidated':bool(row['custom_thresholds']),'explicit_default_decision_parity':bool(row['custom_thresholds']),'native_output_member':f'task-{a.task:04}/native.stdout.json','scoring_deferred_to_step9':True};(work/'integration-comparison.json').write_text(json.dumps({'result':result,'model_expected':expected,'model_differences':differences},indent=2)+'\n')
except Exception as e:
 failure=f'{type(e).__name__}: {e}';(work/'traceback.json').write_text(json.dumps({'traceback':traceback.format_exc()},indent=2)+'\n')
finally:
 (work/'stage.json').write_text(json.dumps({'task':a.task,'gate_pass':failure is None,'failure':failure,'selected_records':1,'completed_records':int(failure is None),'result':result,'argv':argv,'started':started,'ended':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(pathlib.Path(__file__)),'binding_sha256':sha(root/'inputs/integrated-argv-binding.json'),'manifest_sha256':sha(root/'inputs/integrated-context-manifest.json'),'query_sha256':row['query_sha256'],'job_id':os.environ.get('SLURM_JOB_ID'),'allocated_cpus':os.environ.get('SLURM_CPUS_PER_TASK'),'frozen_native_integration_exception':True},indent=2)+'\n');publication=pack(out,work,a.task);(out/f'task-{a.task:04}.json').write_text(json.dumps({'gate_pass':failure is None,'failure':failure,'archive_sha256':publication['archive_sha256']},indent=2)+'\n')
print(json.dumps({'task':a.task,'result':result,'failure':failure}))
if failure:raise SystemExit(1)
