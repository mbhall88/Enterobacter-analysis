"""Display-only aggregation of frozen classification comparisons. No new labels."""
from collections import Counter, defaultdict
import csv, json
from prepare import ROOT, RUN, rows, sha

KRAKEN=ROOT/'results/standard-kraken/2026-10-05/task14-astra-scientific-review/tables'
SCORING=ROOT/'results/tool-planning/2026-10-07/validation-runs-astra/outputs/scoring'
METRICS=ROOT/'results/standard-kraken/2026-10-08/issue26-kraken-brief-codex-061402/tables/metrics.tsv'
SOURCES=[KRAKEN/'taxon-project-pairs.tsv',KRAKEN/'observed-taxonomy-audit.tsv',METRICS]+[SCORING/x/'scored-genomes.tsv' for x in ['layer1','layer3','expert']]

def write_tsv(path, data):
    assert data
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)

def kraken_data():
    taxonomy={r['taxid']:r for r in rows(KRAKEN/'observed-taxonomy-audit.tsv')}
    comparisons=[r for r in rows(KRAKEN/'taxon-project-pairs.tsv') if r['input_mode']=='native' and r['confidence']=='0']
    cells=Counter(); transformation=[]
    for r in comparisons:
        tax=taxonomy[r['raw_taxid']]
        # Collapse strains only to their recorded subspecies/species ancestor.
        # No inferred synonym mapping or new taxonomic interpretation is applied.
        name=tax['subspecies_name'] or tax['species_name'] or r['raw_name']
        project=(r['project_species']+' subsp. '+r['project_subspecies']) if r['project_subspecies'] else (r['project_species'] or 'Species unassigned')
        cells[project,name]+=int(r['bases'])
        transformation.append({**r,'display_kraken_name':name,'display_study_name':project})
    totals=Counter()
    for (project,name),n in cells.items():totals[project]+=n
    data=[{'study_assignment':p,'kraken_name':k,'bases':n,'study_group_bases':totals[p],
           'percent_study_group_bases':100*n/totals[p]}
          for (p,k),n in sorted(cells.items())]
    assert sum(totals.values())==5382995503
    return data,transformation

def pterobacter_data():
    study=rows(SCORING/'layer1/scored-genomes.tsv');public=rows(SCORING/'layer3/scored-genomes.tsv')
    print('Study species outcomes',dict(Counter(r['species_outcome'] for r in study)))
    print('Study label sources',dict(Counter((r['label_status'],r['species_outcome']) for r in study)))
    print('Study subspecies',dict(Counter(r['subspecies_outcome'] for r in study)))
    print('Public species outcomes',dict(Counter(r['species_outcome'] for r in public)))
    return study,public

if __name__=='__main__':
    out=RUN/'outputs/figures/source-data';out.mkdir(parents=True,exist_ok=True)
    data,transformation=kraken_data()
    write_tsv(out/'figure3-and-S13-kraken-taxonomy.tsv',data)
    write_tsv(out/'figure3-and-S13-taxonomy-transformation.tsv',transformation)
    print('Kraken rows',len({r['study_assignment'] for r in data}),'columns',len({r['kraken_name'] for r in data}))
    print('Study assignments',sorted({r['study_assignment'] for r in data}))
    pterobacter_data()
    (RUN/'inputs/v6-data-sources.json').write_text(json.dumps([{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in SOURCES],indent=2)+'\n')
