import csv,json,hashlib,sys,shutil,subprocess
from pathlib import Path
ROOT=Path.cwd();r=ROOT/'results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded'
sys.path.insert(0,str(ROOT/'scripts'));from gtdb_context import parse_newick
read=lambda p:list(csv.DictReader(open(p),delimiter='\t'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,rs):
 with open(p,'w') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
old=read(ROOT/'results/ani-cases/2026-09-26/task12-astra-bugandensis-nematophilus-review/outputs/genomes.tsv');focal={x['accession'] for x in old};assert len(focal)==28
src=ROOT/'results/ani-cases/2026-09-23/task12-astra-evidence/inputs/resolved/genomes.tsv';available={x['accession']:x for x in read(src)}
treedir=ROOT/'results/interactive-trees/2026-09-25/task16-astra-treeprofiler/data';contexts={}
for m in ['bac120_fasttree','verticall_fastme','mashtree_consensus']:
 t=parse_newick((treedir/m/'tree.nwk').read_text());ix={t.label[n]:n for n in t.tips()};n=t.mrca([ix[q] for q in focal]);parent=t.parent[n];contexts[m]=set(t.label[k] for k in t.descendant_tips(parent))-focal
assert len(set.union(*contexts.values()))==4 and all(len(v)==4 for v in contexts.values())
rows=[];labels={x['accession']:x['display_block'] for x in old}
for q in sorted(focal|set.union(*contexts.values())):
 x=available[q];assert sha(x['source'])==x['sha256'];rows.append(dict(accession=q,source=str(Path(x['source']).resolve()),sha256=x['sha256'],group=labels.get(q,'context'),role='existing_focal' if q in focal else 'shared_neighbour_clade',primary=True,reason='frozen_review_focal' if q in focal else 'same_immediate_enclosing_clade_in_all_three_trees'))
meta={x['genbank_accession']:x for x in read(ROOT/'results/ani-cases/2026-09-27/task12-astra-nematophilus-followup/candidates.tsv')}
for x in read(r/'sequences/type-reference-sequence-scopes.tsv'):
 q=x['accession_version'];assert sha(x['chromosome_scope_path'])==x['chromosome_scope_sha256'];primary=q not in ['GCA_901485125.1','GCA_050773425.1'];group='boundary' if q=='GCA_030510045.1' else 'nematophilus_candidate'
 reason='additional_candidate'
 if q=='GCA_901485125.1':reason='QC_hold_64pct_completeness_frameshift_warnings'
 if q=='GCA_050773425.1':reason='118C5_alternate_assembly_identity_unresolved_sensitivity_only'
 if q=='GCA_034107645.1':reason='118C5_original_epidemiology_submission_primary_sensitivity_unit_not_identity_adjudication'
 rows.append(dict(accession=q,source=str(Path(x['chromosome_scope_path']).resolve()),sha256=x['chromosome_scope_sha256'],group=group,role=meta[q]['recommendation'],primary=primary,reason=reason))
assert len(rows)==41 and sum(x['primary'] for x in rows)==39
(r/'fasta').mkdir()
for x in rows:(r/'fasta'/(x['accession']+'.fna')).symlink_to(x['source'])
write(r/'inputs/panel-all.tsv',rows)
write(r/'inputs/panel-primary.tsv',[x for x in rows if x['primary']])
(r/'inputs/all.list').write_text(''.join(str(r/'fasta'/(x['accession']+'.fna'))+'\n' for x in rows))
(r/'inputs/new.list').write_text(''.join(str(r/'fasta'/(q+'.fna'))+'\n' for q in meta if q not in focal))
(r/'inputs/marker.batch.tsv').write_text(''.join(str(r/'fasta'/(q+'.fna'))+'\t'+q+'\n' for q in meta if q not in focal))
(r/'inputs/provenance.json').write_text(json.dumps({'cohort_sha256':sha('manifests/full.tsv'),'types_sha256':sha('manifests/types.tsv'),'resolved_inputs_sha256':sha(src),'contexts':{k:sorted(v) for k,v in contexts.items()},'scope':'NCBI chromosome-assigned sequences; draft scope can include unrecognized plasmids','dddh':'excluded_by_user'},indent=2)+'\n')
cmd=['seqkit','stats','--all','--tabular','--threads','1']+[x['source'] for x in rows]
(r/'inputs/fasta-stats.tsv').write_text(subprocess.check_output(cmd,text=True));print('41 inputs; 39 primary; four shared neighbouring genomes',contexts)
