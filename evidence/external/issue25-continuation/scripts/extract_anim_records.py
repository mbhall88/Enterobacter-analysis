"""Authenticate collected pairs; audit integer metrics against unchanged pyANI 0.3.0 parser."""
import argparse,ast,collections,hashlib,io,json,os,pathlib,subprocess,tarfile,tempfile,typing
import intervaltree
run=pathlib.Path('${PROJECT_ROOT}/results/tool-planning/2026-10-11/issue25-frozen-evaluation-sol61-ab252002');packet=pathlib.Path(__file__).resolve().parents[1];root=pathlib.Path('${PROJECT_ROOT}');ap=argparse.ArgumentParser();ap.add_argument('plan');ap.add_argument('wave');args=ap.parse_args();plan=pathlib.Path(args.plan);wave=args.wave;collected=plan/'collected';profile='reciprocal-mummer4-4.0.1-preserve-n-mum-1-span-v1';snapshot=root/'results/tool-planning/2026-10-06/pilot-coverage-probes-astra/source/installed-pyani-anim.py';snapshot_sha='92f95c578dfe5931a2b8b44b4b3cc803b0d1161e3c125428729d66a3b850fceb'
def sha(p):
 h=hashlib.sha256()
 with pathlib.Path(p).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
assert sha(snapshot)==snapshot_sha;node=next(n for n in ast.parse(snapshot.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='parse_delta');namespace=dict(Path=pathlib.Path,Tuple=typing.Tuple,defaultdict=collections.defaultdict,intervaltree=intervaltree);exec(compile(ast.Module(body=[node],type_ignores=[]),str(snapshot),'exec'),namespace);parse=namespace['parse_delta'];summary=json.loads((collected/f'{wave}.summary.json').read_text());archive=collected/f'{wave}.tar';assert sha(archive)==summary['tar_sha256'];index={r['pair_key']:r for r in map(json.loads,(collected/f'{wave}.index.jsonl').read_text().splitlines())};records=packet/'anim/records';records.mkdir(exist_ok=True);target=records/f'{plan.name}-{wave}.jsonl';part=records/f'.{plan.name}-{wave}.jsonl.part';assert not target.exists() and not part.exists();checks=[];keys=[];counts=collections.Counter();failure=None;binding=json.loads((run/'inputs/p4b-bindings.json').read_text());dependencies={str(p):sha(p) for package in [pathlib.Path(intervaltree.__file__).parent] for p in package.glob('*.py')}
seqkit=pathlib.Path('${USER_HOME}/sw/bin/seqkit').resolve();seqkit_sha=sha(seqkit);length_cache={}
def fasta_lengths(data,digest):
 assert hashlib.sha256(data).hexdigest()==digest
 if digest not in length_cache:
  result=subprocess.run([str(seqkit),'fx2tab','--name','--only-id','--length','--threads','1'],input=data,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
  rows=[line.split('\t') for line in result.stdout.decode().splitlines()]
  assert all(len(x)==2 for x in rows) and len({x[0] for x in rows})==len(rows)
  length_cache[digest]={x[0]:int(x[1]) for x in rows}
 return length_cache[digest]
def oracle(data,ref_lengths,query_lengths):
 blocks=[];refid=queryid=None
 for line in data.decode().splitlines():
  fields=line.split()
  if not fields:continue
  if fields[0].startswith('>'):
   refid=fields[0][1:];queryid=fields[1];assert int(fields[2])==ref_lengths[refid] and int(fields[3])==query_lengths[queryid]
  elif len(fields)==7:
   a,b,c,d,e=map(int,fields[:5]);assert refid and queryid;assert 1<=min(a,b)<=max(a,b)<=ref_lengths[refid] and 1<=min(c,d)<=max(c,d)<=query_lengths[queryid] and e>=0;blocks.append((a,b,c,d,e))
 total=sum(abs(b-a)+abs(d-c)+2 for a,b,c,d,e in blocks);errors=sum(e for a,b,c,d,e in blocks)
 with tempfile.TemporaryDirectory(prefix='p4b-oracle-',dir=os.environ.get('TMPDIR','/tmp')) as work:
  path=pathlib.Path(work)/'native.filter';path.write_bytes(data)
  if total:
   ar,aq,identity,nerrors=parse(path);assert nerrors==errors and abs(identity-(total-2*errors)/total)<1e-12
  else:
   try:parse(path)
   except ZeroDivisionError:pass
   else:raise AssertionError('Reference parser unexpectedly accepts empty alignment')
   ar=aq=0
 return dict(D=total,E=errors,blocks=len(blocks),reference_covered=ar,query_covered=aq,reference_length=sum(ref_lengths.values()),query_length=sum(query_lengths.values()))
try:
 with part.open('x') as output,tarfile.open(archive) as outer:
  for member in outer:
   if member.name.startswith('partial/'):continue
   assert member.isfile() and member.name.endswith('.tar.zst');key=member.name[:-8];entry=index[key];data=outer.extractfile(member).read();assert hashlib.sha256(data).hexdigest()==entry['archive_sha256'];raw=subprocess.run(['zstd','-dc'],input=data,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True).stdout
   with tarfile.open(fileobj=io.BytesIO(raw)) as inner:
    names={m.name:m for m in inner.getmembers() if m.isfile()};task=json.load(inner.extractfile('pair/task.json'));assert task['pair_key']==key and task['exit_status']==entry['exit_status'];assert hashlib.sha256('|'.join([profile,task['query_sha256'],task['reference_sha256']]).encode()).hexdigest()==key
    result_name='pair/out/result.json'
    if result_name not in names:
     assert task['exit_status']!=0;counts['failed_without_result']+=1;continue
    result=json.load(inner.extractfile(result_name));assert result['method']['profile']==profile and result['method']['amendment']=='P1-ANIM-MUMMER4-r1';assert result['query']['sha256']==task['query_sha256'] and result['reference']['sha256']==task['reference_sha256'];assert result['method']['implementation_sha256']==binding['binaries'][str(run/'build/bin/terobacter-dev')]
    for which in ['query','reference']:
     lengths=fasta_lengths(inner.extractfile(f'pair/out/inputs/{which}.fa').read(),result[which]['sha256']);assert lengths==result[which]['record_lengths'] and sum(lengths.values())==result[which]['total_bases']
    for orientation in ['QR','RQ']:
     observed=result['orientations'][orientation];receipt=observed.get('receipt');metrics=observed.get('metrics')
     if metrics is None:continue
     assert receipt is not None and receipt['exit_code']==0;prefix=f'pair/out/{orientation}/'
     for name,digest in receipt['artifacts'].items():assert hashlib.sha256(inner.extractfile(prefix+name).read()).hexdigest()==digest
     ref_lengths=result['reference']['record_lengths'] if orientation=='QR' else result['query']['record_lengths'];query_lengths=result['query']['record_lengths'] if orientation=='QR' else result['reference']['record_lengths'];expected=oracle(inner.extractfile(prefix+'native.filter').read(),ref_lengths,query_lengths)
     if orientation=='RQ':
      expected['reference_covered'],expected['query_covered']=expected['query_covered'],expected['reference_covered'];expected['reference_length'],expected['query_length']=expected['query_length'],expected['reference_length']
     differences={k:dict(expected=v,actual=metrics.get(k)) for k,v in expected.items() if metrics.get(k)!=v};checks.append(dict(pair_key=key,orientation=orientation,exact_match=not differences,differences=differences,no_alignment=expected['D']==0));assert not differences,(key,orientation,differences);counts['audited_orientations']+=1
    if task['exit_status']==0:
     assert task['command_status']==0 and result['status']=='complete';assert all(result['orientations'][x]['status'] in ['aligned','no_alignment'] for x in ['QR','RQ']);assert json.loads(inner.extractfile('pair/stdout.txt').read())==result;output.write(json.dumps(result,separators=(',',':'))+'\n');keys.append(key);counts['complete_pairs']+=1
    else:counts['failed_with_result']+=1
 assert len(keys)==len(set(keys));assert sha(snapshot)==snapshot_sha
 for p,digest in dependencies.items():assert sha(p)==digest
 assert sha(seqkit)==seqkit_sha
 part.rename(target);(records/f'{plan.name}-{wave}.completed-keys.txt').write_text(''.join(k+'\n' for k in keys))
except Exception as error:failure=f'{type(error).__name__}: {error}';raise
finally:
 (packet/'checks'/f'anim-integers-{plan.name}-{wave}.json').write_text(json.dumps(dict(gate_pass=failure is None,failure=failure,planned=summary['planned'],archived_pairs=summary['pairs'],collector_failed=summary['failed'],collector_missing=summary['missing'],counts=dict(counts),checks=checks,archive_sha256=summary['tar_sha256'],parser_snapshot_sha256=snapshot_sha,parser_loaded='Exact AST function node from unchanged pyANI0.3.0 snapshot; unrelated plotting imports omitted. Exact D/E sums are independently computed from the same block headers; pyANI checks E, identity and per-record union coverage. Empty-filter division-by-zero agrees with P1 no-alignment semantics.',dependency_sha256=dependencies,sequence_parser_argv=[str(seqkit),'fx2tab','--name','--only-id','--length','--threads','1'],sequence_parser_sha256=seqkit_sha,independent_input_lengths_checked=True,source_sha256=sha(__file__),job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK')),indent=2)+'\n')
print(json.dumps(dict(counts=dict(counts),failure=failure)))
