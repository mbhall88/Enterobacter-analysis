from pathlib import Path
import hashlib,json,sys
R=Path(__file__).resolve().parents[1];D=Path('data/reads/hormaechei-type-SAMN16357575');O=Path('results/ani-cases/2026-09-29/task12-opus55-hormaechei-reassembly')
def read(p):
 out={};h=None
 for line in p.open():
  if line.startswith('>'):h=line[1:].split()[0];out[h]=''
  else:out[h]+=line.strip().upper()
 return out
def write(p,s):
 with p.open('w') as f:
  for n,seq in s.items():f.write('>'+n+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
paths={'default':D/'autocycler-hifi/autocycler_out/consensus_assembly.fasta','tagged':D/'autocycler-hifi-tagged/autocycler_out/consensus_assembly.fasta','unpolished':D/'autocycler-hifi-tagged/hormaechei_type_hifi_linear.fasta','final':D/'autocycler-hifi-tagged/hormaechei_type_hifi_linear_polished.fasta','published':Path('MoniqueC_ECC_Study/all_assemblies/All_Enterobacter_NCBI/GCA_019048245.1_ASM1904824v1_genomic.fna'),'polypolish':D/'autocycler-hifi-tagged/polish/05_polypolish_polished.fasta','pypolca':D/'autocycler-hifi-tagged/polish/06_pypolca_polished.fasta'}
ss={k:read(p) for k,p in paths.items()}
for k in ['default','final','published','tagged','unpolished']:write(R/'work'/f'{k}.fa',ss[k])
s=ss['default']['1'];repeat=s[484988:517384];assert len(repeat)==32396 and repeat==s[4682910:4715306]
fl={'U2L':s[464988:484988],'U1R':s[517384:537384],'U1L':s[4662910:4682910],'U2R':s[4715306:4735306]}
models={a+'-R-'+b:fl[a]+repeat+fl[b] for a in ['U1L','U2L'] for b in ['U1R','U2R']}
write(R/'work/junction-models.fa',models)
# Separated unique segments plus one full repeat: remove alignment competition between repeat copies.
f=ss['final'];rep=f['chromosome_U1'][-32396:];u2=f['unit_U2'];assert u2[:11508]==rep[20888:] and u2[-16533:]==rep[:16533]
collapsed={'U1_unique':f['chromosome_U1'][:-32396],'U2_unique':u2[11508:-16533],'R_once':rep,'circle_111kb':f['circle_111kb']}
write(R/'work/collapsed.fa',collapsed)
# Verify exact construction without executing upstream script or writing there.
t=ss['tagged'];up=ss['unpolished'];assert up['chromosome_U1']==t['1'][3101764:]+t['1'][:3101764]
assert up['unit_U2']==t['2'][:584600]+'CG'+t['2'][584600:]
assert up['circle_111kb']==t['3']==f['circle_111kb']
x=f['circle_111kb'];y=ss['published']['CP077310.1'];rc=y.translate(str.maketrans('ACGT','TGCA'))[::-1];circle_rotation=(x+x).find(y);circle_rc_rotation=(x+x).find(rc);assert len(x)==len(y) and max(circle_rotation,circle_rc_rotation)>=0
assert f==ss['pypolca']==ss['polypolish'];assert up['unit_U2']==f['unit_U2']
u,v=up['chromosome_U1'],f['chromosome_U1'];i=next(i for i,(a,b) in enumerate(zip(u,v)) if a!=b);assert u[:i]+u[i+1:]==v
record=dict(source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},construction_verified=True,circle_rotation=circle_rotation,circle_rc_rotation=circle_rc_rotation,rotation_cut0=3101764,CG_insertion_before0=584600,polishing_deletion_index0=i,polishing_deleted_base=u[i],rep_len=len(rep),missing_repeat_copy_bases=4355,final_lengths={k:len(v) for k,v in f.items()},collapsed_lengths={k:len(v) for k,v in collapsed.items()})
(R/'outputs/construction.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
