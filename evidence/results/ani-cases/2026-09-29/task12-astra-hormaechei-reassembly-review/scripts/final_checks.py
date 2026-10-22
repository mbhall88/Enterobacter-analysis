from pathlib import Path
import json,csv,subprocess,pysam,hashlib
from read_audit import columns
R=Path(__file__).resolve().parents[1]
truth=json.loads((R/'outputs/power-truth.json').read_text());rows=[]
for lab in ['final-sim','collapsed-sim']:
 obs={x['pos1']:x for x in map(json.loads,(R/'outputs'/f'repeat-{lab}-bq20.jsonl').open())}
 for v in truth['variants']:
  pos=v['repeat_pos1']+(4165525 if lab=='final-sim' else 0);x=obs.get(pos,{});alt=v['alternate'];depth=x.get('depth',0);n=x.get('indels' if alt.startswith('-') else 'bases',{}).get(alt,0)
  rows.append(dict(reference=lab,repeat_pos1=v['repeat_pos1'],alternate=alt,depth=depth,alt_count=n,frequency=n/depth if depth else None,detected=depth>=20 and .2<=n/depth<=.8))
with (R/'outputs/power-check.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
out=[]
for label,c,s in [('collapsed','R_once',32145),('final','chromosome_U1',4165525+32145),('final','unit_U2',32145-20888),('final','chromosome_U1',4181845),('final','unit_U2',584600)]:
 b=pysam.AlignmentFile(str(R/'work'/f'{label}.bam'));f=pysam.FastaFile(str(R/'work'/f'{label}.fa'));v=list(columns(b,f,c,s,s+1,20))[0];v.update(label=label,context=f.fetch(c,s-20,s+21));out.append(v)
(R/'outputs/repeat-focus-sites.json').write_text(json.dumps(out,indent=2)+'\n')
f=pysam.FastaFile(str(R/'work/published.fa'));s=f.fetch('CP077309.1')
(R/'work/overlap-head.fa').write_text('>head\n'+s[:40000]+'\n');(R/'work/overlap-tail.fa').write_text('>tail\n'+s[-40000:]+'\n')
cmd=['${USER_HOME}/sw/miniforge3/envs/autocycler/bin/minimap2','-t','1','-cx','asm5',str(R/'work/overlap-tail.fa'),str(R/'work/overlap-head.fa')]
with (R/'outputs/terminal-overlap.paf').open('w') as out,(R/'logs/terminal-overlap.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
(R/'logs/terminal-overlap-command.json').write_text(json.dumps(cmd)+'\n')
p=(R/'outputs/terminal-overlap.paf').read_text().split('\t');assert [int(p[i]) for i in [2,3,7,8,9,10]]==[0,17418,22571,40000,17401,17436]
print('Power table, focus sites and independent terminal-overlap check complete')
