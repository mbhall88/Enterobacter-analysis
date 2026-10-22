#!/usr/bin/env python3
"""Rebuild all 23 numbered TSVs from saved upstream evidence; standard library only."""
import argparse,csv,gzip,hashlib,io,json,re
from collections import Counter
from pathlib import Path
import dataset_flow as flow

def sha(raw):return hashlib.sha256(raw).hexdigest()
def table(path,values,fields=None):
    # Match the historical csv writer's quoting, order and CRLF line ending.
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(values[0]),delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(values)
class Input:
    def __init__(self,item,base):self.item=item;self.base=base
    def read_bytes(self):
        raw=(self.base/self.item['file']).read_bytes()
        if sha(raw)!=self.item['export_sha256']:raise ValueError('Compressed input checksum mismatch: '+self.item['key'])
        value=gzip.decompress(raw)
        if sha(value)!=self.item['decoded_export_sha256']:raise ValueError('Decoded input checksum mismatch: '+self.item['key'])
        return value
    def read_text(self):return self.read_bytes().decode('utf-8')
def rows(p):return list(csv.DictReader(io.StringIO(p.read_text(),newline=''),delimiter='\t'))
class TablePaths:
    def __init__(self,out,index):self.out=out;self.names={x['original_source_table']:x['file'] for x in index}
    def __truediv__(self,name):return self.out/self.names[name]

def execute(base,out):
    index=json.loads((base/'expected-tables.json').read_text())
    inventory=json.loads((base/'input-inventory.json').read_text())
    inputs={r['key']:Input(r,base) for r in inventory}
    if len(inputs)!=len(inventory):raise ValueError('Duplicate input keys')
    # Validate the entire allowlist before writing results.
    for p in inputs.values():p.read_bytes()
    out.mkdir(exist_ok=False)
    tab=TablePaths(out,index)
    assignments=rows(inputs['assignments']);marker={r['tip_id']:r for r in rows(inputs['marker'])}
    fields=['accession','in_frozen_study','record_role','assignment_status','project_species','project_subspecies','review_lineage','historical_ncbi','current_ncbi','historical_kraken_chromosome','historical_kraken_whole','formal_type_names','nomenclature_and_boundary_flags','decision_authority','decision_or_retention_basis','usage_restrictions','quality_flags','genomovar_id']
    # Original order includes current_ncbi_date immediately after current_ncbi.
    fields.insert(fields.index('current_ncbi')+1,'current_ncbi_date')
    for a in assignments:
        m=marker.get(a['accession'])
        a.update({'frozen_unique_bac120':m['number_unique_genes'] if m else 'not_in_1090_preparation','frozen_missing_bac120':m['number_missing_genes'] if m else 'not_in_1090_preparation','supplementary_low_marker_flag':str(int(m['number_missing_genes'])>=5).lower() if m else 'not_assessed','marker_scope':'Task11 1090-input preparation; repaired inherited/extracted chromosome and draft-anchor inputs; September 2026'})
    fields+=['frozen_unique_bac120','frozen_missing_bac120','supplementary_low_marker_flag','marker_scope']
    table(tab/'S01-accessions.tsv',assignments,fields)
    type_ids=set(re.findall(r'GCA_\d+\.\d+',inputs['type_projection'].read_text()));by={r['accession']:r for r in assignments};display=[]
    for acc in sorted(type_ids):
        a=by[acc];species=a['project_species'];sub=a['project_subspecies'];current=species+(' subsp. '+sub if sub else '')
        if not current:current=a['historical_ncbi']
        if not current:current={'GCA_000982825.1':'Cronobacter sakazakii (outgroup)','GCA_019048125.1':'Klebsiella aerogenes (outgroup)'}[acc]
        display.append({'accession':acc,'current_project_label':current,'historical_type_names':a['formal_type_names'],'historical_ncbi':a['historical_ncbi'],'flags':a['nomenclature_and_boundary_flags'],'source_scope':'35-tip projection of frozen 1089 inference; current labels from final additive register'})
    assert len(display)==35
    table(tab/'S02-type-display-map.tsv',[{'type_display_id':f'T{i+1:02}',**r} for i,r in enumerate(display)])
    # Literal count routine extracted from the original dataset_flow.py, with only documentary receipt writing removed.
    flow.PATHS={k:inputs['assignments' if k=='accessions' else k] for k in ['accessions','initial_transitions','reference_panel','development_split','genomovar_members','native_execution','public_pool','confirmation','linked_hold','public_reserve','public_truth','public_held_truth','native_preparation']}
    flow.rows=rows;flow.Counter=Counter;flow.json=json;flow.csv=csv;flow.tab=tab;flow.dataset_flow()
    # Original manuscript copied these saved measurement/register products verbatim.
    for key,p in inputs.items():
        if key.startswith('saved_'):(out/f'Table{key[6:]}.tsv').write_bytes(p.read_bytes())
    # Original code_catalogues projection, preserving nested diagnostics and original full source digest.
    keys=['record_id','coded_id','query_sha256','dataset','setting','category','applicable','triggered','manifest_witness_candidate_references','reported_witness_candidate_references','triggered_witness_candidate_references','passing_witness_references','contexts_with_passing_witness','foreign_challenge_references','contexts_with_foreign_challenge','method_veto_references','contexts_with_method_veto','expected_crossing_references','contexts_with_expected_crossing','incomplete','rescued','complete_panel','reference_id','competitor_evaluation']
    nested=['counts','complete_panel_challenge_origins','actual_reference_diagnostics','native_reasons','candidate_reference_ids','passing_witness_reference_ids','challenge_reference_ids','retained_method_veto_reference_ids','thresholds','member_compatibility_not_assessed_by_this_report','native_diagnostic','independent_diagnostic']
    for n in [12,13,14]:
        p=inputs[f'jsonl_S{n:02}'];values=[]
        for line in p.read_text().splitlines():
            r=json.loads(line);view={k:r.get(k,'') for k in keys};view['source_record_sha256']=p.item['source_sha256']
            for k in nested:
                if k in r:view[k]=json.dumps(r[k],sort_keys=True)
            if 'independent_diagnostic' in r:
                for k in ['expected_ANI_method_crossing','vetoes','gate_disagreement']:view[k]=r['independent_diagnostic'].get(k,'')
            if 'reciprocal_metrics' in r:view['reciprocal_metrics']=json.dumps({orientation:{k:v for k,v in metric.items() if k!='receipt'} for orientation,metric in r['reciprocal_metrics'].items()},sort_keys=True)
            values.append(view)
        if not values:(out/f'TableS{n:02}.tsv').write_text('coded_id\tdataset\tsetting\tcategory\tforeign_challenge_references\n')
        else:
            fields=list(dict.fromkeys(k for x in values for k in x));table(out/f'TableS{n:02}.tsv',[{k:x.get(k,'') for k in fields} for x in values])
    cases=json.loads(inputs['safety'].read_text())
    table(out/'TableS15.tsv',[{'coded_id':c['record_id'],'dataset':c['dataset'],'source_group':c['source_group'],'omission_match':c['diagnostic_would_assign'],'deployment_genomovar':'null','full_reference_outcome':'withheld_overlap (two compatible groups)','P3_conformance':'passed','scientific_disposition':'approved conditional matching; missing-reference limitation; thresholds unchanged'} for c in cases])
    for n in [17,18,19]:
        records=rows(inputs[f'cost_S{n:02}']);fields=[k for k in records[0] if k not in ['source_archive','resource_member']]
        table(out/f'TableS{n:02}.tsv',[{k:r[k] for k in fields} for r in records])
    verification=[]
    for x in index:
        actual=sha((out/x['file']).read_bytes());verification.append({'table':x['table'],'file':x['file'],'expected_sha256':x['sha256'],'actual_sha256':actual,'exact':actual==x['sha256']})
    report={'status':'pass' if all(x['exact'] for x in verification) else 'fail','tables':verification,'inputs_verified':len(inputs),'scope':'Reconstruction of 23 manuscript table transformations from saved measurements/aggregates; no new upstream scientific analysis.','software':'Python standard library only'}
    (out/'reconstruction-checks.json').write_text(json.dumps(report,indent=2)+'\n')
    if report['status']!='pass':raise ValueError('Table mismatch: '+', '.join(x['table'] for x in verification if not x['exact']))
    print(json.dumps({'status':'pass','exact_tables':len(verification),'inputs_verified':len(inputs)},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--output',type=Path,required=True,help='New output directory; existing directories are refused');a=p.parse_args();execute(a.inputs,a.output)
