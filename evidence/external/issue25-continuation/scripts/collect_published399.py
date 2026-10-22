import hashlib,json,pathlib,shutil,os,subprocess
root=pathlib.Path(__file__).resolve().parents[1]
base=root.parent
accept=json.loads((base/'inputs/coordinator-acceptance.json').read_text())
b=base/'inputs/storage-binding.json'
assert accept['accepted'] is True and hashlib.sha256(b.read_bytes()).hexdigest()==accept['binding_sha256']
assert json.loads((root/'checks/pilot-compute.json').read_text())['comparison_executed'] is False
inventory=json.loads((root/'inputs/published399-source-inventory.json').read_text());assert len(inventory['published'])==399
plan=root/'anim/plan-published399';dest=plan/'pairs/wave-0001';copied=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(1048576),b''):h.update(z)
 return h.hexdigest()
for e in inventory['published']:
 source=pathlib.Path(e['path']);target=dest/source.name
 assert sha(source)==e['sha256'] and source.stat().st_size==e['bytes']
 with source.open('rb') as src,target.open('xb') as dst:shutil.copyfileobj(src,dst,1048576);dst.flush();os.fsync(dst.fileno())
 assert sha(target)==e['sha256'] and sha(source)==e['sha256']
 copied.append({'source':str(source),'copy':str(target),'sha256':e['sha256'],'bytes':e['bytes']})
(root/'checks/published399-copy-readback.json').write_text(json.dumps({'copied':copied,'originals_unchanged':True,'biological_reruns':0},indent=2)+'\n')
settings=json.loads(b.read_text());frozen=pathlib.Path(settings['frozen_collector']);assert sha(frozen)==settings['collector_sha256']
argv=['python3',str(frozen),str(plan/'waves/wave-0001.tsv'),str(dest),str(plan/'collected')]
(root/'receipts/published399-frozen-collector.argv.json').write_text(json.dumps(argv,indent=2)+'\n')
subprocess.run(argv,check=True)
summary=json.loads((plan/'collected/wave-0001.summary.json').read_text());assert summary['planned']==summary['pairs']==399 and summary['failed']==summary['missing']==0
print(json.dumps({'published399_collection_pass':True,'summary':summary,'native_comparisons':0}))
