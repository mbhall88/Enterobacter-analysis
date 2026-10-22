"""Final computational checks. Human figure/scientific review remains separate."""
from pathlib import Path
import csv, json, hashlib
from Bio import SeqIO
r=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
rows=list(csv.DictReader((r/'inputs/panel-qc.tsv').open(),delimiter='\t'))
assert len(rows)==41 and len({x['accession'] for x in rows})==41
for row in rows:assert sha(row['source'])==row['sha256'],row['accession']
checks=['outputs/ani-phase-validation.json','outputs/tree-validation.json',
        'interface-test-v2/validation-masked-reproduced.json']
checks += ['figures/local-trees-'+c+'-validation.json' for c in ['primary','masked','RY']]
for p in checks:assert json.loads((r/p).read_text())['status']=='PASS',p
for ref in ['GCA_019046905.1','GCA_026344075.1']:
    o=r/'local'/ref
    v=json.loads((o/'validation.json').read_text());assert v['strict_core_bases']>1000000
    assert json.loads((o/'tests/validation.json').read_text())['status']=='PASS'
    for folder in [o/'logs',o/'tests/logs']:
        for p in folder.glob('*.command.json'):
            rec=json.loads(p.read_text())
            if 'exit_code' in rec:assert rec['exit_code']==0,str(p)
    counts={len(rec.seq) for rec in SeqIO.parse(o/'primary.core.fna','fasta')}
    assert counts=={v['strict_core_bases']}
paths=list((r/'scripts').glob('*.py'))+list((r/'scripts').glob('*.sh'))
paths+=list((r/'inputs').glob('*.tsv'))+list((r/'outputs').glob('*'))
paths=[p for p in paths if p.is_file() and p.name!='final-validation.json']
manifest={str(p.relative_to(r)):sha(p) for p in sorted(paths)}
(r/'inputs/final-code-and-output-hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
(r/'outputs/final-validation.json').write_text(json.dumps(dict(
    status='PASS',input_hashes=41,checks=checks,dddh='not_run',
    human_review='pending: plots, convergence, model assumptions and interpretation',
    task_complete=False),indent=2)+'\n')
print('PASS computational checks; human review is still required')
