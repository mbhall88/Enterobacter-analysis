from io import StringIO
from pathlib import Path
from Bio import Phylo
from result_helpers import splits,key,au_table,tips
x=Phylo.read(StringIO('((A:1,B:1)99:1,(C:1,D:1)88:1,E:1);'),'newick');s=splits(x);assert key({'A','B'},tips(x)) in s;assert key({'A','C'},tips(x)) not in s
before=set(s);x.root_with_outgroup('E');assert set(splits(x))==before
sample='USER TREES\n  1 -10622448.06 3.3136e-05   0.502 +  0.521 +  0.829 +       0.5 +        1 + \n  2 -10627622.79  5174.7       0 -      0 -      0 -         0 - 2.99e-05 - \n';a=au_table(sample);assert len(a)==2 and a[1]['p_AU']==1 and a[2]['p_AU']==2.99e-5
old=Path('${PROJECT_ROOT}/results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded/local/GCA_019046905.1/tests/topology-primary/AU.iqtree');assert len(au_table(old.read_text()))==4
print('PASS: unrooted split membership, reroot invariance, absent clade, synthetic and existing AU parsing')
