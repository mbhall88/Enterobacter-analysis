#!/usr/bin/env python3
"""Source-linked S1 diagnostic coverage plot."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import task12_scope_s1_s2 as s
rs=[r for r in s.rows(s.RUN/'s1-records.tsv') if r['status']=='screened']
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(8.2,6.2))
colors={'chromosome_like':'#D55E00','plasmid_like':'#0072B2','shared_or_mixed':'#009E73','unresolved_low_or_intermediate_coverage':'#999999'}
for k,c in colors.items():
 rr=[r for r in rs if r['classification']==k]
 ax.scatter([float(r['chromosome_pct']) for r in rr],[float(r['plasmid_pct']) for r in rr],s=16 if k!='chromosome_like' else 40,c=c,alpha=.5 if k!='chromosome_like' else .9,label=k.replace('_',' ')+f' (n={len(rr)})',edgecolors='none')
shared=[r for r in rs if r['classification']=='chromosome_like' and int(r['self_chromosome_aligned_bp'])/int(r['bp'])>=.9]
ax.scatter([float(r['chromosome_pct']) for r in shared],[float(r['plasmid_pct']) for r in shared],s=90,facecolors='none',edgecolors='black',label='Flag also >=90% covered by own chromosome (n=4)')
ax.axvline(50,color='#D55E00',ls='--',lw=.8);ax.axhline(20,color='#D55E00',ls='--',lw=.8)
ax.set(xlim=(-2,102),ylim=(-2,102),xlabel='Record bases aligning to neighbour chromosome (%)',ylabel='Record bases aligning to neighbour plasmid(s) (%)',title='S1: label-specific interval unions at >=95% identity')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.16),frameon=False,fontsize=8,ncol=1)
fig.text(.5,.015,'1,367 screened records; 1,308 below 50 kb and 164 without eligible neighbour not screened.\nMatches may reflect mobile/shared sequence; they do not prove physical linkage.',ha='center',fontsize=8)
fig.subplots_adjust(bottom=.32,left=.12,right=.97,top=.92)
for ext in ['png','pdf']:fig.savefig(s.RUN/('s1-coverage-diagnostics.'+ext),dpi=200)
