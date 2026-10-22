"""Mask ambiguous alignment overlaps and entire indel neighbourhoods on fixed coordinates."""
import gzip
from pathlib import Path
import numpy as np
from Bio import SeqIO

def masked_sequences(folder, padding=10):
 folder=Path(folder);seqs={x.id:np.frombuffer(str(x.seq).upper().encode(),dtype='S1').copy() for x in SeqIO.parse(folder/'snippy.pseudo.fna','fasta')}
 masks={k:np.zeros(len(v),dtype=bool) for k,v in seqs.items()}
 with gzip.open(folder/'snippy.vcf.gz','rt') as f:
  for l in f:
   if l.startswith('#'):continue
   fields=l.rstrip().split('\t');chrom,pos,ref,alt=fields[0],int(fields[1])-1,fields[3],fields[4]
   if len(ref)!=len(alt) or alt.startswith('<'):
    end=pos+len(ref)
    if alt.startswith('<'):
     info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x);end=int(info.get('END',end))
    masks[chrom][max(0,pos-padding):min(len(masks[chrom]),end+padding)]=True
 # One-to-one alignment can retain overlapping local blocks; do not fill these from reference.
 depth={k:np.zeros(len(v)+1,dtype=np.int16) for k,v in seqs.items()}
 for l in (folder/'snippy.coords.tsv').open():
  t=l.rstrip().split('\t')
  if len(t)<13:continue
  k=t[-2];start=min(int(t[0]),int(t[1]))-1;end=max(int(t[0]),int(t[1]));depth[k][start]+=1;depth[k][end]-=1
  masks[k][max(0,start-padding):min(len(masks[k]),start+padding)]=True
  masks[k][max(0,end-padding):min(len(masks[k]),end+padding)]=True
 for k in seqs:
  masks[k]|=np.cumsum(depth[k])[:-1]!=1
  masks[k]|=~np.isin(seqs[k],np.array(list('ACGT'),dtype='S1'))
  seqs[k][masks[k]]=b'N'
 return {k:v.tobytes().decode() for k,v in seqs.items()}
