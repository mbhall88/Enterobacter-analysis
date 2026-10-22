from pathlib import Path
import json
R=Path(__file__).resolve().parents[1];url=json.loads((R/'logs/comment-post.json').read_text())['html_url']
brief=f'**[Hormaechei type reassembly independently reviewed, 29 Sep]({url})**, Codex Astra (`gpt-6-astra`): broader spanning search negative (378 HiFi/11,808 subreads ≥20kb), BAM/ENA long-read equality checked,31 validation checks pass. Original47-versus1 disagreement positions reproduce; repeat identity/accuracy/topology claims qualified after a placement-power test. Final FASTA leaves4,355bp of one repeat copy unrepresented. Independent ATCC ONT/Illumina assembly found (reported circles4,787,102/111,400bp; entitlement required, raw reads unavailable);2016 public draft exactly spans one proposed fused connection. Recommend retain frozen baseline with all-replicon sensitivity, reassembly as corroboration; Michael’s reference/scope decision pending. No reference/cohort change or NCBI contact; full Task12 open.'
p=Path('docs/analysis-execution-plan.md');s=p.read_text();lines=s.splitlines()
for i,line in enumerate(lines):
 if line.startswith('| 12 |'):
  line=line.replace('| 12 |','| 12 | **Independent hormaechei type reassembly review complete29Sep; reference/scope decision pending, Task12 open (Codex Astra, `gpt-6-astra`).** ',1);lines[i]=line;break
else:raise AssertionError('Task12 row absent')
s='\n'.join(lines)+'\n';head='## 12 — ANI/Galah cases, nomenclature and corrected numbers\n';assert head in s;s=s.replace(head,head+'\n'+brief+'\n',1);p.write_text(s)
p=Path('docs/project-hub.md');s=p.read_text();head='# Enterobacter project hub\n';assert s.startswith(head);p.write_text(s.replace(head,head+'\n'+brief+'\n',1))
p=Path('docs/group-review-pack-plan.md');s=p.read_text();head='## Current readiness — 29 September\n';assert head in s;s=s.replace(head,head+'\n'+brief+'\n',1)
old='Michael has chosen to attempt read-based reassembly first; revisit reference scope after that validation, with NCBI reporting deferred.'
new=f'Michael’s reassembly and the [independent review]({url}) are complete. Present the unresolved32-kb repeat layout, omitted4,355bp repeat copy sequence, placement-power caveat and independent ATCC/2016-draft evidence. Reference/scope choice remains pending; NCBI reporting remains deferred. Deposited linear metadata do not establish physical linearity or circularity.'
assert old in s;s=s.replace(old,new,1);p.write_text(s)
p=Path('docs/task-status/12-cases-naming-denominators.md')
with p.open('a') as f:f.write('\n## 29 September — independent hormaechei reassembly review completed\n\n'+brief+'\n\nRun: `'+str(R.relative_to(Path.cwd()))+'`. The original Opus task/output and read tree were preserved. Commands, input hashes, mappings/pileups, source snapshots, synthetic fixtures,31-check validation, public-only export and the single reviewed issue comment are retained. The initially slow pileup readers were stopped and replaced by equivalent cached decoding; SRA-header normalization was corrected before drawing a BAM/FASTQ conclusion. No jobs remain from this review.\n\nLimitations: no instrumented assembler execution trace; ATCC sequence/raw reads unavailable; no held-out accuracy estimate; repeat-copy phasing/physical topology unresolved. No proof that a prospective reference change leaves ANI/AF unchanged. Next unblocked work: Michael’s explicit reference/scope decision and Task16 integration of this reviewed evidence into the hormaechei panel. Optional ATCC sequence comparison requires an existing permitted access route. Overall Task12 remains open.\n')
print(url)
