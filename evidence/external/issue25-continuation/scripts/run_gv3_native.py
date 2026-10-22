"""GV3 fixed native contexts; authenticate reused actual R1 and RG2 parent,full P3 replay."""
import pathlib,json,tarfile,hashlib,argparse,tempfile,traceback,os,signal
import evaluation_common as c
import p3_adapter
from stage_archive import pack
root=c.CONT;sha=c.sha
ap=argparse.ArgumentParser();ap.add_argument('--task',type=int,required=True);a=ap.parse_args();manifest=json.loads((root/'inputs/gv3-native-context-plan.json').read_text());selected=manifest[a.task*50:(a.task+1)*50];assert selected;out=root/'outputs/gv3';out.mkdir(parents=True,exist_ok=True);assert not(out/f'task-{a.task:04}.tar.gz').exists() and not(out/f'.task-{a.task:04}.part').exists();work=pathlib.Path(tempfile.mkdtemp(prefix='p4b-gv3-',dir='/tmp'));results=[];cache={};pins={};failure=None

def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def archived(source,member):
 path=source['archive'];assert sha(path)==source['archive_sha256'];pins[path]=source['archive_sha256']
 if path not in cache:
  data={}
  with tarfile.open(path) as t:
   m=next(x.name for x in t.getmembers() if x.name.endswith('/artifact-hashes.json'));prefix=m.rsplit('/',1)[0];hashes=json.load(t.extractfile(m))
   for name,want in hashes.items():
    raw=t.extractfile(prefix+'/'+name).read();assert hashlib.sha256(raw).hexdigest()==want;data[prefix+'/'+name]=raw
  cache[path]=data
 return json.loads(cache[path][member])
def interrupt(sig,frame):raise RuntimeError('Owned GV3 task interrupted '+str(sig))
signal.signal(signal.SIGTERM,interrupt);signal.signal(signal.SIGINT,interrupt)
try:
 assert json.loads((root/'checks/step7-summary.json').read_text())['gate_pass'],'All frozenstep7 including integrated/custom cases must pass first'
 v2index=json.loads((root/'inputs/full-v2-result-index.json').read_text());assert len(manifest)==914
 for row in selected:
  case=work/row['record_id'];case.mkdir();report=archived(row['actual_R1'],row['actual_R1']['report_member']);idx=v2index[row['coded_id']];v2=archived(idx,idx['v2_member']);assert idx['dataset']=='native' and idx['setting']=='as_run';assert v2['species_id']==report['assignment']['species_id'] and v2['subspecies_id']==report['assignment']['subspecies_id'] and not v2['trigger'];assert v2['call_basis']!='ANIm_assisted'
  section=report['genomovar'];search=section['search'];assert report['execution']['status']=='complete' and report['input']['query_sha256']==row['query_sha256']==section['parent']['query_sha256']==search['query_sha256'];assert search['status']=='completed' and search['diagnostic_exclusion'] is False;assert (section['parent']['species_id'],section['parent']['subspecies_id'])==('ECTX-SP025','ECTX-SS007')
  evaluated=p3_adapter.evaluate(section,search,row['record_id'],row['dataset'],source_group=row['source_group'],self_member=row['self_member'],singleton_source=row['singleton_source'],parent_basis=v2['call_basis']);write(case/'full-p3-comparison.json',evaluated);write(case/'actual-genomovar-section.json',section);write(case/'source-bindings.json',{'actual_R1':row['actual_R1'],'actual_v2':idx,'parent_equivalence_from_RG2':True,'source_accession':row['source_accession']});assert not evaluated['model_differences'],'GV3 full frozen P3 mismatch';write(case/'scorer-record.json',evaluated['record']);results.append({'record_id':row['record_id'],'dataset':row['dataset'],'self_member':row['self_member'],'query_sha256':row['query_sha256'],'model_pass':True,'scoring_deferred_to_step9':True})
except Exception as e:
 failure=f'{type(e).__name__}: {e}';write(work/'traceback.json',{'traceback':traceback.format_exc()})
finally:
 write(work/'stage.json',{'task':a.task,'gate_pass':failure is None,'failure':failure,'selected_records':len(selected),'completed_records':len(results),'pending_records':len(selected)-len(results),'results':results,'input_sha256':pins,'source_sha256':sha(__file__),'adapter_sha256':sha(root/'scripts/p3_adapter.py'),'manifest_sha256':sha(root/'inputs/gv3-native-context-plan.json'),'job_id':os.environ.get('SLURM_JOB_ID'),'new_measurements':False,'scoring_executed':False});published=pack(out,work,a.task);write(out/f'task-{a.task:04}.json',{'gate_pass':failure is None,'failure':failure,'completed_records':len(results),'archive_sha256':published['archive_sha256']})
if failure:raise RuntimeError(failure)
