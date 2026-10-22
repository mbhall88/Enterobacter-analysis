"""Retained Linux fixed60 baseline: exact native rows,44reciprocal pairs,and R1/R2 calls."""
import pathlib,json,csv,tarfile,hashlib,tempfile,datetime,platform,collections
import evaluation_common as c
from stage_archive import pack
root=c.CONT;run=c.RUN;sha=c.sha
assert platform.machine()=='x86_64';assert json.loads((root/'checks/step7-summary.json').read_text())['gate_pass'];manifest=json.loads((root/'inputs/linux-platform-manifest.json').read_text());assert len(manifest)==60
frozen=c.P4A/'outputs/fixed-sets-attempt01/platform-panel.tsv';expected=list(csv.DictReader(frozen.open(),delimiter='\t'));assert {(r['coded_id'],r['layer'],r['state'],r['reason']) for r in manifest}=={(r['coded_id'],r['layer'],r['state'],r['reason']) for r in expected};assert len({r['coded_id'] for r in manifest})==60
engineplan=list(csv.DictReader((run/'inputs/engine-parity-manifest-r7.tsv').open(),delimiter='\t'));locations={r['coded_id']:(i//20,r) for i,r in enumerate(engineplan) if r['dataset']=='native'};enginegate=json.loads((run/'checks/engine-parity-summary.json').read_text());assert enginegate['gate_pass'] and enginegate['failed_comparisons']==0;index=json.loads((root/'inputs/full-v2-result-index.json').read_text());contexts={r['record_id']:r for r in json.loads((root/'inputs/final-v2-context-manifest.json').read_text())};catalogue={(r['query']['sha256'],r['reference']['id']):r for r in map(json.loads,(root/'inputs/final-anim-catalogue.jsonl').read_text().splitlines())};work=pathlib.Path(tempfile.mkdtemp(prefix='p4b-linux60-',dir='/tmp'));pins={};results=[];pairkeys=set();out=root/'outputs/platform-linux';out.mkdir(exist_ok=True)
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def readarchive(path,want,members):
 path=pathlib.Path(path);assert sha(path)==want;pins[str(path)]=want
 with tarfile.open(path) as t:
  hm=next(x.name for x in t.getmembers() if x.name.endswith('/artifact-hashes.json'));prefix=hm.rsplit('/',1)[0];hashes=json.load(t.extractfile(hm));data={}
  for name,digest in hashes.items():
   raw=t.extractfile(prefix+'/'+name).read();assert hashlib.sha256(raw).hexdigest()==digest
   if prefix+'/'+name in members:data[prefix+'/'+name]=raw
  assert set(data)==set(members);return data
for row in manifest:
 code=row['coded_id'];case=work/code;case.mkdir();assert sha(pathlib.Path(row['query_fasta']))==row['query_sha256'];task,eng=locations[code];assert eng['query_sha256']==row['query_sha256'];prefix=f'task-{task:04d}/'+code;ep=run/'outputs/engine-parity-recovery01'/f'task-{task:04}.json'
 if not ep.exists() or not json.loads(ep.read_text())['gate_pass']:ep=run/'outputs/engine-parity'/f'task-{task:04}.json'
 er=json.loads(ep.read_text());assert er['gate_pass'];archive=ep.with_suffix('.tar.gz');assert str(archive) in enginegate['input_sha256'] or str(ep) in enginegate['input_sha256'];wanted=[]
 for component in ['primary_whole','genomovar_chromosome']:wanted.extend(prefix+'/'+component+x for x in ['.embedded.stdout','.comparison.json','.embedded.command.json'])
 data=readarchive(archive,er['archive_sha256'],wanted)
 for name,raw in data.items():(case/name.rsplit('/',1)[-1]).write_bytes(raw)
 for component in ['primary_whole','genomovar_chromosome']:assert json.loads(data[prefix+'/'+component+'.comparison.json'])['exact_match']
 data=readarchive(row['source_R1_archive'],row['source_R1_archive_sha256'],[row['source_R1_report_member']]);v1=json.loads(data[row['source_R1_report_member']]);assert v1['input']['query_sha256']==row['query_sha256'];write(case/'R1.json',v1);idx=index[code];data=readarchive(idx['archive'],idx['archive_sha256'],[idx['v2_member']]);v2=json.loads(data[idx['v2_member']]);assert v2['input']['query_sha256']==row['query_sha256'];write(case/'R2.json',v2);pairs=[]
 if row['state']=='triggered_with_witness_candidate':
  scope=contexts[code]['anim_scope']
  for request in scope:
   r=catalogue[row['query_sha256'],request['reference_id']];assert r['reference']['sha256']==request['reference_sha256'];pairs.append(r);pairkeys.add((row['query_sha256'],request['reference_id']))
 write(case/'ANIm.json',pairs);write(case/'source-bindings.json',{'fixed_manifest_row':row,'engine_archive':str(archive),'engine_archive_sha256':er['archive_sha256'],'R2_source':idx});results.append({'coded_id':code,'query_sha256':row['query_sha256'],'state':row['state'],'engine_components':2,'anim_pairs':len(pairs),'v1_v2_retained':True})
assert len(pairkeys)==44;write(work/'stage.json',{'gate_pass':True,'records':60,'engine_components':120,'anim_pairs':44,'anim_orientations':88,'states':dict(collections.Counter(r['state'] for r in results)),'results':results,'input_sha256':pins,'source_sha256':sha(__file__),'frozen_manifest_sha256':sha(frozen),'kernel':platform.platform(),'architecture':platform.machine(),'source_scope':'Actual checksum-verified Linux native measurement and R1/R2 outputs reused on exact original inputs; no biological reruns','cross_platform_acceptance':False,'macOS_issue':'mbhall88/Enterobacter#30'});publication=pack(out,work,0);receipt={'gate_pass':True,'Linux_baseline_records':60,'primary_and_member_query_components':120,'ANIm_pairs':44,'ANIm_orientations':88,'archive_sha256':publication['archive_sha256'],'source_sha256':sha(__file__),'cross_platform_acceptance':False,'Mac_execution_pending_issue30':True};f=root/'checks/linux-platform-summary.json';assert not f.exists();write(f,receipt);print(json.dumps(receipt))
