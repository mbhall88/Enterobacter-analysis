import sys, json
from pathlib import Path
sys.path.insert(0, 'scripts')
from summarise_concordance import read_alignment
from codon_partitions import gene_ranges
R = Path('results/conflict-resolution/2026-09-23/task17-opus55-resolve')
N = {'GCA_001750725.1': 'lud', 'GCA_000025565.1': 'clo', 'GCA_002939185.1': 'sich', 'GCA_001521715.1': 'asb',
     'GCA_019048125.1': 'Kaer', 'GCA_019048245.1': 'horm', 'GCA_000211415.1': 'mori', 'GCA_003594935.2': 'huax'}
out = {}
for panel, quartets in [('C1-ludwigii', [('GCA_002939185.1', 'GCA_019048125.1'), ('GCA_001521715.1', 'GCA_019048125.1')]),
                        ('C12x-extended-sampling', [('GCA_001521715.1', o) for o in ('GCA_019048125.1', 'GCA_019048245.1', 'GCA_000211415.1', 'GCA_003594935.2')])]:
    L = R / panel / 'core-genome/ledger-panel'
    aln = read_alignment(L / 'alignment.fna')
    genes = [(g, a, b) for g, a, b in gene_ranges(L / 'concatenation-genes.nex') if (b - a + 1) % 3 == 0]
    for x, og in quartets:
        q = ['GCA_001750725.1', 'GCA_000025565.1', x, og]
        s = [aln[t] for t in q]
        counts = {p: {'lud+clo': 0, 'lud+X (clo outside)': 0, 'clo+X (lud outside)': 0} for p in (1, 2, 3)}
        for _, a, b in genes:
            for i in range(a - 1, b):
                c = [seq[i].upper() for seq in s]
                if any(ch not in 'ACGT' for ch in c) or len(set(c)) != 2:
                    continue
                pos = (i - (a - 1)) % 3 + 1
                if c[0] == c[1] and c[2] == c[3]: counts[pos]['lud+clo'] += 1
                elif c[0] == c[2] and c[1] == c[3]: counts[pos]['lud+X (clo outside)'] += 1
                elif c[0] == c[3] and c[1] == c[2]: counts[pos]['clo+X (lud outside)'] += 1
        key = f"{panel}: lud, clo, {N[x]} | {N[og]}"
        out[key] = {}
        for p, d in counts.items():
            tot = sum(d.values())
            out[key][f'pos{p}'] = {'informative': tot, **{k: round(v / tot, 3) for k, v in d.items()}}
        print(key, json.dumps(out[key]))
json.dump(out, open(R / 'C12x-extended-sampling/core-genome/c1-site-patterns-by-codon-position.json', 'w'), indent=1)
