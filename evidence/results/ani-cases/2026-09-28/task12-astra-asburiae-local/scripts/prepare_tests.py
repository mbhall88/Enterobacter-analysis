"""Freeze clade hypotheses and count BAC120 information before local inference completes."""
from pathlib import Path
from io import StringIO
import csv,json,hashlib
import numpy as np
from Bio import SeqIO,Phylo
R=Path(__file__).resolve().parents[1];rows=list(csv.DictReader((R/'inputs/panel.tsv').open(),delimiter='\t'));groups={}
for x in rows:groups.setdefault(x['group'],set()).add(x['accession'])
ids={x['accession'] for x in rows};assert len(ids)==93
A,M,D,P,Q=[groups[k] for k in ['asburiae_type_core','shared_middle','dykesii_type_core','peripheral_three','peripheral_six']]
hypotheses={'middle_with_dykesii':M|D,'middle_with_asburiae':M|A,'P3_with_asburiae':P|A,'P6_with_asburiae_and_P3':Q|A|P,'P6_with_middle_and_dykesii':Q|M|D}
(R/'constraints').mkdir();ledger=[]
for label,ingroup in hypotheses.items():
 text='('+','.join(sorted(ids-ingroup))+',('+','.join(sorted(ingroup))+'));\n';tree=Phylo.read(StringIO(text),'newick');assert {x.name for x in tree.get_terminals()}==ids;assert {x.name for x in tree.common_ancestor(sorted(ingroup)).get_terminals()}==ingroup
 path=R/'constraints'/(label+'.nwk');path.write_text(text);ledger.append(dict(hypothesis=label,clade_members=sorted(ingroup),constraint=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),tips=93))
(R/'inputs/hypotheses.json').write_text(json.dumps(ledger,indent=2)+'\n')
source=Path('${PROJECT_ROOT}/results/full-cohort/2026-09-22/task11-opus-full/marker/filtered/alignment.faa');seqs={x.id:str(x.seq) for x in SeqIO.parse(source,'fasta') if x.id in ids};assert len(seqs)==93
out=[]
for label,ss in [('all93',ids),('focal90',A|M|D|P|Q),('A9',A),('M69',M),('D3',D),('P3',P),('P6',Q)]:
 arr=np.array([list(seqs[q]) for q in sorted(ss)]);counts=np.array([(arr==b).sum(axis=0) for b in 'ACDEFGHIKLMNPQRSTVWY']);states=(counts>0).sum(axis=0);pi=(counts>=2).sum(axis=0)>=2
 out.append(dict(panel=label,taxa=len(ss),columns=arr.shape[1],variable_columns=int((states>=2).sum()),parsimony_informative_columns=int(pi.sum()),constant_observed_columns=int((states==1).sum()),no_standard_amino_acid_columns=int((states==0).sum()),distinct_sequences=len({seqs[q] for q in ss}),scope='existing_full_cohort_filtered_BAC120_amino_acid_alignment'))
with (R/'outputs/bac120-information.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
with (R/'inputs/bac120-matched93.faa').open('w') as f:
 for q in sorted(ids):f.write('>'+q+'\n'+seqs[q]+'\n')
(R/'outputs/test-design-validation.json').write_text(json.dumps(dict(status='PASS',hypotheses=5,taxa=93,marker_input=str(source),marker_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),blocks={k:len(v) for k,v in groups.items()}),indent=2)+'\n')
print(json.dumps(out,indent=2))
