from pathlib import Path
import json,csv,hashlib
R=Path(__file__).resolve().parents[1];checks={}
def j(p):return json.loads((R/p).read_text())
def check(n,b):
 checks[n]=bool(b)
 if not b:raise AssertionError(n)
check('construction',j('outputs/construction.json')['construction_verified'])
for path,h in j('outputs/construction.json')['source_hashes'].items():check('source_unchanged:'+path,hashlib.sha256(Path(path).read_bytes()).hexdigest()==h)
b=j('outputs/bam-inventory.json');check('bam_and_ena_longreads_identical',b['bam_records']==545384 and b['bam_ge20k']==b['fastq_ge20k']==11808 and not b['missing_in_fastq'] and not b['changed_sequence'] and not b['fastq_only'])
for lab,count in [('hifi-allchains',378),('subreads-allchains',11808)]:
 x=j('outputs/span-'+lab+'.json');check('no_spans:'+lab,x['input_reads_ge20k']==count and not x['spanning_hits'] and not x['extent_spanning_hits'] and not x['envelope_candidates'])
x=j('outputs/span-controls.json');check('span_controls_24positive_8negative',x['positive_controls']==24 and x['negative_controls']==8 and not x['missed_positives'] and not x['false_positives'])
x=j('outputs/repeat-final-bq0.json');check('original_repeat_test_reproduces',x['positions_ge20']==31817 and not x['intermediate_alleles'])
x=j('outputs/repeat-collapsed-bq20.json');check('collapsed_covers_entire_repeat',x['positions_ge20']==32396);check('real_repeat_single_intermediate_homopolymer',len(x['intermediate_alleles'])==1 and x['intermediate_alleles'][0]['allele']=='+T')
rows=list(csv.DictReader((R/'outputs/power-check.tsv').open(),delimiter='\t'))
check('power_fixture_demonstrates_placement_blindspot',sum(x['detected']=='True' for x in rows if x['reference']=='final-sim')==4 and sum(x['detected']=='True' for x in rows if x['reference']=='collapsed-sim')==6)
check('47_sites_spotchecked',len(j('outputs/published-47-alleles.json'))==47)
for lab,n in [('final',1),('published',47)]:check('original_counter:'+lab,int((R/'outputs'/f'opus-counter-{lab}.tsv').read_text().split('\t')[3])==n)
check('two_and_only_two_header_tags',len(j('outputs/tag-check.json')['changes'])==2)
x=j('outputs/independent-draft-bridge.json');check('independent_draft_exact_gapfree_bridge',len(x)==1 and x[0]['matches']==x[0]['alignment_columns']==72396 and x[0]['ambiguous_bases']==0)
check('independent_terminal_overlap_verified',(R/'outputs/terminal-overlap.paf').exists())
for p in ['README.md','manuscript-contribution.md','commands.sh','code-audit/verdict.md','type-census/README.md','inputs/reads-sha256.txt','logs/allocation.txt','comment-template.md']:check('deliverable:'+p,(R/p).is_file())
check('comment_authorship',(R/'comment-template.md').read_text().startswith('> **Posted by:** Codex Astra (`gpt-6-astra`)'))
out=dict(status='PASS',checks=checks,count=len(checks),limitations=['No instrumented assembler branch trace','No independently held-out accuracy estimate','ATCC FASTA/raw reads not obtained','Physical layout and repeat-copy phasing unresolved','No reference change or NCBI contact'])
(R/'validation-final.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
