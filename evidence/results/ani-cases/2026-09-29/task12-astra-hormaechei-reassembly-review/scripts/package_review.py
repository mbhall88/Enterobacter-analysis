"""Export only reviewed public-genome evidence, without raw reads or issue history."""
from pathlib import Path
import shutil,json,hashlib,re
R=Path(__file__).resolve().parents[1];P=R/'package';P.mkdir(exist_ok=True)
files=['README.md','manuscript-contribution.md','commands.sh','validation-final.json','inputs/reads-sha256.txt','logs/versions.txt',
'outputs/span-subreads-allchains.json','outputs/span-hifi-allchains.json','outputs/span-controls.json','outputs/bam-inventory.json','outputs/power-check.tsv','outputs/repeat-collapsed-bq20.json','outputs/homopolymer-flank-check.json','outputs/construction.json','outputs/tag-check.json','outputs/repeat-focus-sites.json','outputs/published-47-alleles.tsv','outputs/independent-draft-bridge.json','outputs/terminal-overlap.paf',
'type-census/README.md','type-census/ncbi-type-material-census.tsv','type-census/atcc-assembly-metadata.json','type-census/atcc-replicons.tsv','type-census/ncbi-replicons.tsv',
'code-audit/verdict.md','code-audit/hifiasm-graph-check.json','code-audit/check_graph.py','code-audit/fetch-record.json']
files += [str(p.relative_to(R)) for p in sorted((R/'scripts').glob('*')) if p.suffix in ['.py','.sh'] and p.name not in ['publish.py','package_review.py','handover.py']]
manifest=[]
for f in files:
 src=R/f;dst=P/f;dst.parent.mkdir(parents=True,exist_ok=True);data=src.read_bytes()
 if f=='code-audit/verdict.md':data=data.replace(b'[sources](sources/)',b'source snapshots (retained in the local run)')
 dst.write_bytes(data);manifest.append(dict(path=f,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
for md in P.rglob('*.md'):
 for target in re.findall(r'\]\(([^)]+)\)',md.read_text()):
  if target.startswith(('https:','http:','#')):continue
  assert (md.parent/target.split('#')[0]).exists(),(md,target)
(P/'export-manifest.json').write_text(json.dumps(dict(scope='Reviewed public type-genome evidence only; no raw reads, private clinical metadata, issue history, or source checkout',files=manifest),indent=2)+'\n')
print('Packaged',len(files),'reviewed files; relative Markdown links pass')
