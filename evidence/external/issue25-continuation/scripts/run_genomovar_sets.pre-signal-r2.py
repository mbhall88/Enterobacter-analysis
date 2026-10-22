"""Frozen GV1b/GV2/GV6/GV4; actual native R1/R5 and full P3a replay,scoring deferred."""
import pathlib,json,argparse,datetime,subprocess,tempfile,time,tarfile,hashlib,traceback,os
import evaluation_common as c
import p3_adapter
from stage_archive import pack
root=c.CONT;run=c.RUN;sha=c.sha
ap=argparse.ArgumentParser();ap.add_argument('--task',type=int,required=True);a=ap.parse_args();binding=json.loads((root/'inputs/binding.json').read_text());manifest=json.loads((root/'inputs/genomovar-context-plan.json').read_text());selected=manifest[a.task*20:(a.task+1)*20];assert selected;out=root/'outputs/genomovar';out.mkdir(parents=True,exist_ok=True);assert not(out/f'task-{a.task:04}.tar.gz').exists() and not(out/f'.task-{a.task:04}.part').exists(),'Refuse completed or partial task rerun';work=pathlib.Path(tempfile.mkdtemp(prefix=f'p4b-gv-task-{a.task:04}-',dir='/tmp'));results=[];commands=[];pins={};cache={};failure=None

def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def command(argv,stem):
 receipt={'argv':argv,'started':datetime.datetime.now(datetime.timezone.utc).isoformat()};start=time.monotonic()
 with pathlib.Path(str(stem)+'.stdout.json').open('x') as stdout,pathlib.Path(str(stem)+'.stderr').open('x') as stderr:r=subprocess.run(['/usr/bin/time','-v','-o',str(stem)+'.resources.txt',*argv],stdout=stdout,stderr=stderr)
 receipt.update(exit_status=r.returncode,wall_seconds=time.monotonic()-start);write(pathlib.Path(str(stem)+'.command.json'),receipt);commands.append(receipt);assert r.returncode==0,f'Actual CLI failure {stem.name}: {r.returncode}';return json.loads(pathlib.Path(str(stem)+'.stdout.json').read_text())
def archived(source,stem):
 path=source['archive']
 if path not in cache:
  assert sha(path)==source['archive_sha256'];pins[path]=source['archive_sha256'];data={}
  with tarfile.open(path) as t:
   m=next(m.name for m in t if m.name.endswith('/artifact-hashes.json'));prefix=m.rsplit('/',1)[0];h=json.load(t.extractfile(m))
   for name,want in h.items():
    b=t.extractfile(prefix+'/'+name).read();assert hashlib.sha256(b).hexdigest()==want;data[prefix+'/'+name]=b
  cache[path]=data
 data=cache[path][source['report_member']];pathlib.Path(str(stem)+'.stdout.json').write_bytes(data);write(pathlib.Path(str(stem)+'.reuse.json'),{'source':source,'source_report_sha256':hashlib.sha256(data).hexdigest(),'actual_completed_output_reused':True,'no_new_measurements':True});return json.loads(data)
try:
 assert int(os.environ['SLURM_CPUS_PER_TASK'])>=2
 assert json.loads((root/'checks/step7-summary.json').read_text())['gate_pass'],'Frozenstep7 including integrated/custom must pass before GVsets'
 for path,want in binding['binaries'].items():assert sha(path)==want
 for row in selected:
  if(out/'gate-stop.json').exists():raise RuntimeError('Another genomovar task stopped the gate')
  case=work/row['record_id'];case.mkdir();assert sha(row['query_fasta'])==row['query_sha256'];pins[row['query_fasta']]=row['query_sha256'];base=['--reference-package',binding['reference_package']['path'],'--expected-sha256',binding['reference_package']['manifest_sha256']]
  report=archived(row['prior_R1'],case/'r1') if row['prior_R1'] else command([str(run/'build/bin/tero'),*base,'--policy','conservative-assembly-v1','--threads','2','--report-dir',str(case/'r1-reports'),row['query_fasta']],case/'r1')
  assert report['execution']['status']=='complete' and report['input']['query_sha256']==row['query_sha256'] and report['input']['reference_package_sha256']==binding['reference_package']['manifest_sha256']
  search=report['genomovar']['search']
  # A not-run search is kept as actual parent-veto evidence,never replaced by invented rows.
  if search['status']!='completed':
   write(case/'capability-gap.json',{'record_id':row['record_id'],'dataset':row['dataset'],'kind':row.get('kind'),'actual_search':search,'actual_parent':report['genomovar'].get('parent'),'deployed_genomovar':report['genomovar'].get('genomovar_id'),'R5_requires_completed_native_search':True,'no_fabricated_or_relabelled_search':True});raise RuntimeError('R5 capability gap: actual native R1 genotype search not completed;coordinator disposition required')
  assert search['query_sha256']==row['query_sha256'] and search['reference_package_sha256']==binding['reference_package']['manifest_sha256'];mask=case/'primary-mask.txt';mask.write_text(''.join(x+'\n' for x in row['primary_mask']))
  if row.get('actual_R5_reuse'):r5=archived(row['actual_R5_reuse'],case/'r5')
  else:
   argv=[str(run/'build/bin/terobacter-dev'),'genomovar-diagnostic','--report',str(case/'r1.stdout.json'),'--expected-report-sha256',sha(case/'r1.stdout.json'),*base,'--setting',row['setting'],'--exclude-references',str(mask),'--exclude-references-sha256',sha(mask)]
   for x in row['masked_members']:argv+=['--mask-member',x]
   for x in row['removed_groups']:argv+=['--remove-group',x]
   r5=command(argv,case/'r5')
  assert r5['genomovar_id'] is None and r5['setting']==row['setting'];assert set(r5['parent_recomputation']['excluded_reference_ids'])==set(row['primary_mask']);section=r5['diagnostic'];assert section['parent']['query_sha256']==row['query_sha256'];assert set(section['diagnostic_actions']['masked_members'])==set(row['masked_members']);assert set(section['diagnostic_actions']['removed_groups'])==set(row['removed_groups'])
  evaluated=p3_adapter.evaluate(section,search,row['record_id'],row['dataset'],source_group=row['source_group'],singleton_source=row['singleton_source']);write(case/'full-p3-comparison.json',evaluated);assert not evaluated['model_differences'],'Full P3 contract replay mismatch';write(case/'scorer-record.json',evaluated['record']);results.append({'record_id':row['record_id'],'dataset':row['dataset'],'kind':row.get('kind'),'query_sha256':row['query_sha256'],'model_pass':True,'actual_parent':section['parent'],'masked_members':row['masked_members'],'removed_groups':row['removed_groups'],'deployment_genomovar_id':r5['genomovar_id'],'diagnostic_would_assign':r5['diagnostic_would_assign_genomovar_id'],'scoring_deferred_to_step9':True,'prior_R1_reused':row['prior_R1'] is not None,'actual_R5_reused':bool(row.get('actual_R5_reuse'))})
except Exception as e:
 failure=f'{type(e).__name__}: {e}';write(work/'traceback.json',{'traceback':traceback.format_exc()})
 try:
  with(out/'gate-stop.json').open('x') as h:json.dump({'task':a.task,'failure':failure,'completed':len(results),'selected':len(selected)},h,indent=2)
 except FileExistsError:pass
finally:
 write(work/'stage.json',{'task':a.task,'gate_pass':failure is None,'failure':failure,'selected_records':len(selected),'completed_records':len(results),'pending_records':len(selected)-len(results),'results':results,'commands':commands,'input_sha256':pins,'source_sha256':sha(pathlib.Path(__file__)),'binding_sha256':sha(root/'inputs/binding.json'),'genomovar_manifest_sha256':sha(root/'inputs/genomovar-context-plan.json'),'P3_adapter_sha256':sha(root/'scripts/p3_adapter.py'),'job_id':os.environ.get('SLURM_JOB_ID'),'allocated_cpus':os.environ.get('SLURM_CPUS_PER_TASK'),'biological_scoring':False});publication=pack(out,work,a.task);write(out/f'task-{a.task:04}.json',{'gate_pass':failure is None,'failure':failure,'selected_records':len(selected),'completed_records':len(results),'archive_sha256':publication['archive_sha256']})
print(json.dumps({'task':a.task,'completed':len(results),'failure':failure}))
if failure:raise SystemExit(1)
