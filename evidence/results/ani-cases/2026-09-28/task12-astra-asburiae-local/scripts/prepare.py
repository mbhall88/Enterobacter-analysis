from pathlib import Path
import csv,json,hashlib,shutil,sys
R=Path(__file__).resolve().parents[1];ROOT=Path('${PROJECT_ROOT}');OLD=ROOT/'results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded'
sys.path.insert(0,str(ROOT/'scripts'));from gtdb_context import parse_newick
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rs):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
source=ROOT/'results/ani-cases/2026-09-26/task12-astra-asburiae-dykesii-review/outputs/genomes.tsv';review=rows(source);assert len(review)==90
g={x['accession']:x for x in rows(ROOT/'results/ani-cases/2026-09-23/task12-astra-evidence/inputs/resolved/genomes.tsv')};focal={x['accession'] for x in review};groups={x['accession']:x['display_block'] for x in review};extra={'GCA_025502525.1':'outside_boundary','GCA_019968765.1':'outside_boundary_partner','GCA_034008295.1':'chinensis_context'}
panel=[]
for q in sorted(focal|set(extra)):
 x=g[q];p=Path(x['source']).resolve();assert sha(p)==x['sha256'];panel.append(dict(accession=q,source=str(p),sha256=x['sha256'],group=groups.get(q,extra.get(q)),focal=q in focal))
write(R/'inputs/panel.tsv',panel);write(R/'inputs/previous-genome-evidence.tsv',review)
# Selection justified by the existing two-parent BAC120 context; other methods retained for reference.
TREES=ROOT/'results/interactive-trees/2026-09-25/task16-astra-treeprofiler/data';contexts={}
for m in ['bac120_fasttree','verticall_fastme','mashtree_consensus']:
 p=TREES/m/'tree.nwk';t=parse_newick(p.read_text());ix={t.label[k]:k for k in t.tips()};n=t.mrca([ix[q] for q in focal]);assert {t.label[k] for k in t.descendant_tips(n)}==focal
 pp=t.parent[t.parent[n]];contexts[m]={'two_parent_tips':len(t.descendant_tips(pp)),'source':str(p),'sha256':sha(p)}
 if m=='bac120_fasttree':assert {t.label[k] for k in t.descendant_tips(pp)}==set(g['accession'] for g in panel)
# Immutable copy of the exact previously exercised assembly-calling implementation.
shutil.copytree(OLD/'pinned',R/'pinned',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
for rel,h in json.loads((OLD/'pinned/source-hashes.json').read_text()).items():assert sha(R/'pinned'/rel)==h
shutil.copyfile(OLD/'scripts/mask_pseudo.py',R/'scripts/mask_pseudo.py')
shutil.copyfile(OLD/'scripts/submit.py',R/'scripts/submit.py')
(R/'inputs/provenance.json').write_text(json.dumps({'panel':93,'focal':90,'contexts':contexts,'review_source':str(source),'review_sha256':sha(source),'reused_pinned_implementation':str(OLD/'pinned'),'mask_function_sha256':sha(R/'scripts/mask_pseudo.py'),'reused_fixture_validation':str(OLD/'interface-test-v2/validation-masked.json'),'naming':'open; user leaning not final','dddh':'not_run','gubbins':'not_run'},indent=2)+'\n')
print('93 FASTA hashes and pinned source hashes passed')
