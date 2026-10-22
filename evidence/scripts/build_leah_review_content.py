#!/usr/bin/env python3
"""Build the private review's explicitly curated content and bundled public-genome figures."""
import json, hashlib, shutil, argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'review-site');p.add_argument('--record-dir',type=Path,default=ROOT/'results/project-review/2026-09-30/task16-astra-review-site');a=p.parse_args()
 data=json.loads((ROOT/'config/leah-review-content.json').read_text()); manifest=[]
 for g in data['groups']:
  for f in g['figures']:
   source=ROOT/f.pop('source'); assert source.is_file(),source
   name=g['id']+'-'+source.name; dest=a.output/'public/evidence'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest)
   f['url']='/evidence/'+name
   f['sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
   manifest.append({'group':g['id'],'source':str(source.relative_to(ROOT)),'asset':f['url'],'sha256':f['sha256']})
   pdf=source.with_suffix('.pdf')
   if pdf.exists():
    target=dest.with_suffix('.pdf');shutil.copyfile(pdf,target);f['pdf']='/evidence/'+target.name
    manifest.append({'group':g['id'],'source':str(pdf.relative_to(ROOT)),'asset':f['pdf'],'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()})
  for src in g.pop('localSources',[]):
   s=ROOT/src;assert s.is_file(),s
   manifest.append({'group':g['id'],'source':src,'sha256':hashlib.sha256(s.read_bytes()).hexdigest()})
 canonical=json.dumps(data,sort_keys=True,ensure_ascii=False).encode();data['revision']=hashlib.sha256(canonical).hexdigest()[:16]
 (a.output/'data').mkdir(exist_ok=True)
 (a.output/'data/review.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
 run=a.record_dir;run.mkdir(parents=True,exist_ok=True)
 (run/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({'groups':len(data['groups']),'questions':sum(len(g['questions']) for g in data['groups']),'figures':sum(len(g['figures']) for g in data['groups']),'revision':data['revision'],'source_records':len(manifest)}))
if __name__=='__main__':main()
