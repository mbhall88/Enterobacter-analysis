#!/usr/bin/env python3
"""Freeze an INTERNAL development reserve and balanced marker-discovery panel.
Not an untouched external test or a production reference selector. All paths
come from explicit manifests. Config supplies thresholds; no taxonomy changes.
"""
import argparse, collections, csv, hashlib, json, math
from pathlib import Path

def read(p):
 with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rows,fields=None):
 with open(p,'w') as h:
  w=csv.DictWriter(h,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def digest(p):
 with open(p,'rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
class Union:
 def __init__(self,ids):self.p={x:x for x in ids}
 def find(self,x):
  while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
  return x
 def join(self,a,b):
  a,b=self.find(a),self.find(b)
  if a!=b:self.p[max(a,b)]=min(a,b)
def seeded(seed,x):return hashlib.sha256(f'{seed}|{x}'.encode()).hexdigest()
def main(config,out):
 cfg=json.loads(Path(config).read_text());out=Path(out);out.mkdir(exist_ok=False)
 inputs=cfg['inputs'];hashes=[dict(path=p,sha256=digest(p)) for p in [config,__file__,*inputs.values()]]
 recovery={r['accession']:r for r in read(inputs['recovery'])};ids=set(recovery)
 assignments={r['accession']:r for r in read(inputs['assignments'])};study={r['accession_version']:r for r in read(inputs['study'])}
 meta={}
 for key in ['type_metadata','study_metadata']:
  for line in open(inputs[key]):
   r=json.loads(line);meta[r['accession']]=r
 assert ids<=set(meta) and ids<=set(assignments)
 assert not ids & set(cfg['excluded'])
 nodes={}
 for a in sorted(ids):
  m=meta[a];ass=assignments[a];info=m.get('assembly_info',{});sample=info.get('biosample',{});qc=m.get('checkm_info',{})
  comp=qc.get('completeness');cont=qc.get('contamination');qpass=comp is not None and cont is not None and comp>=cfg['min_completeness'] and cont<=cfg['max_contamination']
  species=ass['project_species'];sub=ass['project_subspecies'];named=species.startswith('Enterobacter ') and species!='Enterobacter sp.' and not ass['usage_restrictions']
  group=species+(' subsp. '+sub if sub else '') if named else 'context:'+ (ass['review_lineage'] or species or a)
  nodes[a]=dict(accession=a,group=group,project_species=species,project_subspecies=sub,named_target=int(named),type_anchor=int(recovery[a]['type_anchor']),biosample=sample.get('accession',''),bioproject=info.get('bioproject_accession',''),strain=sample.get('strain',''),entity=study.get(a,{}).get('selected_inference_entity_id',''),completeness=comp if comp is not None else '',contamination=cont if cont is not None else '',qc_pass=int(qpass),quality_flags=ass['quality_flags'],assembly_level=info.get('assembly_level',''),whole_records=m.get('assembly_stats',{}).get('number_of_contigs',''),source_role=ass['record_role'])
 # Split grouping is conservative: edges protect against leakage, not prove isolate identity.
 clone=Union(ids);edges=[]
 def join_keys(field,union,reason):
  seen={}
  for a in sorted(ids):
   value=nodes[a].get(field,'')
   if not value:continue
   if value in seen:union.join(a,seen[value]);edges.append(dict(a=a,b=seen[value],reason=reason,value=value))
   else:seen[value]=a
 join_keys('biosample',clone,'same_BioSample');join_keys('entity',clone,'adjudicated_inference_entity')
 whole_seen={}
 for a in sorted(ids):
  value=recovery[a]['recommended_whole_sha256']
  if value in whole_seen:clone.join(a,whole_seen[value]);edges.append(dict(a=a,b=whole_seen[value],reason='same_whole_file_hash_not_identity_claim',value=value))
  else:whole_seen[value]=a
 pairs={};pair_count=0
 with open(inputs['pairs']) as h:
  for r in csv.DictReader(h,delimiter='\t'):
   a,b=r['a'],r['b']
   if a not in ids or b not in ids or a==b:continue
   key=tuple(sorted((a,b)));assert key not in pairs,('duplicate pair',key)
   val=float(r['ani']) if r['status']=='measured' and r['ani'] else None
   pairs[key]=val;pair_count+=1
   if val is not None and val>=cfg['near_clone_ani'] and min(float(r['af_a']),float(r['af_b']))>=cfg['near_clone_min_af']:
    clone.join(a,b);edges.append(dict(a=a,b=b,reason='chromosome_ANI_similarity_block_not_identity',value=f"{val};AF={r['af_a']},{r['af_b']}"))
 assert pair_count==len(ids)*(len(ids)-1)//2,('incomplete pair ledger',pair_count,len(ids))
 blocks=Union(ids)
 for a in ids:blocks.join(a,clone.find(a))
 join_keys('bioproject',blocks,'same_BioProject_conservative_split_block')
 # One reproducible 20% hash allocation per linked block; types always development.
 groups=collections.defaultdict(list)
 for a,n in nodes.items():
  n['clone_block']=clone.find(a);n['split_block']=blocks.find(a)
  if n['named_target']:groups[n['group']].append(a)
 forced={blocks.find(a) for a,n in nodes.items() if n['type_anchor']}
 for group,aa in groups.items():
  eligible=[a for a in aa if nodes[a]['qc_pass'] or nodes[a]['type_anchor']]
  if eligible and not any(blocks.find(a) in forced for a in eligible):
   forced.add(min({blocks.find(a) for a in eligible},key=lambda b:seeded(cfg['seed'],b)))
 reserve={b for b in {blocks.find(a) for a in ids} if b not in forced and int(seeded(cfg['seed'],b)[:16],16)/16**16<cfg['reserve_fraction']}
 for a,n in nodes.items():n['partition']='internal_reserve' if blocks.find(a) in reserve else 'development'
 selected={a for a,n in nodes.items() if n['type_anchor'] and n['named_target']}
 reasons={a:'historical_type_anchor;QC_exception_explicit' if not nodes[a]['qc_pass'] else 'historical_type_anchor' for a in selected}
 def distance(a,b):
  if a==b:return 0.0
  ani=pairs[tuple(sorted((a,b)))];return 30.0 if ani is None else 100-ani
 for group,aa in sorted(groups.items()):
  chosen=sorted(set(aa)&selected)
  candidates={a for a in aa if nodes[a]['partition']=='development' and nodes[a]['qc_pass']}
  def qc_key(a):return (-float(nodes[a]['completeness']),float(nodes[a]['contamination']),int(nodes[a]['whole_records'] or 99999),a)
  # First supplemental anchor is the within-stratum medoid if no historical type exists.
  if not chosen and candidates:
   first=min(candidates,key=lambda a:(sum(distance(a,b) for b in candidates)/len(candidates),qc_key(a)));chosen.append(first);selected.add(first);reasons[first]='development_medoid_no_historical_type'
  while len(chosen)<cfg['max_per_stratum']:
   used={clone.find(a) for a in chosen};available=[a for a in candidates-selected if clone.find(a) not in used]
   if not available:break
   a=min(available,key=lambda a:(-min(distance(a,b) for b in chosen),qc_key(a)));chosen.append(a);selected.add(a);reasons[a]='farthest_first_development_diversity'
 for a,n in nodes.items():
  n['in_discovery']=int(a in selected);n['selection_reason']=reasons.get(a,'not_selected_for_discovery')
  n['evaluation_role']='internal_concordance_only' if n['partition']=='internal_reserve' else 'development_only'
  if not n['named_target']:n['evaluation_role']+=';context_or_unnamed_challenge_not_marker_training'
  if not n['qc_pass']:n['evaluation_role']+=';QC_missing_or_below_development_gate'
  n['whole_fasta']=recovery[a]['recommended_whole_path'];n['whole_sha256']=recovery[a]['recommended_whole_sha256']
 boundary=[]
 for e in edges:
  if nodes[e['a']]['partition']!=nodes[e['b']]['partition']:boundary.append(e)
 assert not boundary
 assert all(nodes[a]['partition']=='development' for a in selected)
 assert {a for a,n in nodes.items() if n['type_anchor'] and n['named_target']}<=selected
 # Never expose an internal reserve as independent accuracy or train on it.
 write(out/'split-manifest.tsv',list(nodes.values()));write(out/'discovery-panel.tsv',[nodes[a] for a in sorted(selected)]);write(out/'leakage-edges.tsv',edges)
 counts=[]
 for group,aa in sorted(groups.items()):
  counts.append(dict(group=group,total=len(aa),split_blocks=len({blocks.find(a) for a in aa}),clone_blocks=len({clone.find(a) for a in aa}),development=sum(nodes[a]['partition']=='development' for a in aa),internal_reserve=sum(nodes[a]['partition']=='internal_reserve' for a in aa),discovery=sum(a in selected for a in aa),types=sum(nodes[a]['type_anchor'] for a in aa),QC_pass=sum(nodes[a]['qc_pass'] for a in aa)))
 write(out/'group-accounting.tsv',counts);write(out/'input-hashes.tsv',hashes)
 summary=dict(genomes=len(ids),discovery=len(selected),types_in_discovery=sum(nodes[a]['type_anchor'] for a in selected),strata=len(groups),partitions=dict(collections.Counter(n['partition'] for n in nodes.values())),split_blocks=len({blocks.find(a) for a in ids}),clone_blocks=len({clone.find(a) for a in ids}),max_split_block=max(collections.Counter(blocks.find(a) for a in ids).values()),cross_partition_edges=len(boundary),near_clone_ani=cfg['near_clone_ani'],near_clone_min_af=cfg['near_clone_min_af'],pair_ledger_scope='historical_chromosome;not_whole_assembly_distance',pair_ledger_complete=True,unseen_external_test='not_available_not_claimed',missing_group_challenge='future_leave_one_project_species_out_refit;do_not_reuse_full_marker_discovery_as_independent_test',stages='split_and_discovery_frozen_before_new_orthology')
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args();main(a.config,a.output)
