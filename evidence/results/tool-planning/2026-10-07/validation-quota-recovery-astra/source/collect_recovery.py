"""Merge verified original outputs and infrastructure recovery; retain both attempt histories."""
from pathlib import Path
from collections import Counter
import csv,json,sys,tarfile,hashlib
ORIGINAL=Path('${PROJECT_ROOT}/results/tool-planning/2026-10-07/validation-runs-astra')
RECOVERY=ORIGINAL.with_name('validation-quota-recovery-astra')
sys.path.insert(0,str(ORIGINAL/'source'))
from common import read,write,sha,json_write,verify_object,command
from collect import audit_query
csv.field_size_limit(1_000_000_000)


def verify_archive():
    record=json.loads((ORIGINAL/'outputs/ancillary-archive.json').read_text())
    archive=ORIGINAL/'outputs/completed-native-ancillary.tar.gz'
    inventory=ORIGINAL/'outputs/ancillary-inventory.tsv'
    if sha(archive)!=record['archive_sha256'] or sha(inventory)!=record['inventory_sha256']:
        raise ValueError('archival storage changed')
    expected={r['member']:r for r in read(inventory)};seen=set()
    if len(expected)!=record['files_archived']:raise ValueError('archive count differs')
    with tarfile.open(archive,'r|gz') as handle:
        for member in handle:
            if not member.isfile() or member.name not in expected or member.name in seen:raise ValueError('archive member mismatch')
            if hashlib.file_digest(handle.extractfile(member),'sha256').hexdigest()!=expected[member.name]['sha256']:
                raise ValueError('archived member digest mismatch')
            seen.add(member.name)
    if seen!=set(expected):raise ValueError('archive missing member')
    checkpoint=json.loads((ORIGINAL/'outputs/archive-verified-checkpoint.json').read_text())
    if checkpoint['inventory_sha256']!=sha(inventory):raise ValueError('checkpoint inventory differs')
    ledger=ORIGINAL/'outputs/prearchive-successful-ledger.tsv'
    if sha(ledger)!=checkpoint['ledger_sha256']:raise ValueError('prearchive audit ledger changed')
    if checkpoint['archive_sha256']!=record['archive_sha256']:raise ValueError('archive/checkpoint mismatch')
    return expected,{r['genome_id']:r for r in read(ledger)}


def archived_original(row,previous,packed):
    """Reverify original successful evidence across loose files plus verified archive."""
    native=ORIGINAL/'queries'/row['coded_id']/'attempt01/analyze'
    finished=native.parent/'finished.json'
    if sha(finished)!=previous['finished_sha256'] or sha(native/'run.json')!=previous['native_run_sha256']:
        raise ValueError('original audit record changed')
    record=json.loads(finished.read_text());run=json.loads((native/'run.json').read_text())
    if (record['query_sha256']!=row['query_sha256'] or record['genome_id']!=row['genome_id']
        or record['coded_id']!=row['coded_id'] or record['exit_status']!=0
        or record['query_unchanged'] is not True or record['infrastructure_error']
        or record['argv']!=command(Path(row['fasta_path']),native)
        or run['positive_calling_policy_enabled'] is not False or run['status']!='complete'):
        raise ValueError('original outcome binding differs')
    loose={str(p.relative_to(native)) for p in native.rglob('*') if p.is_file() and p!=native/'run.json'}
    if loose & set(packed) or loose | set(packed)!=set(run['artifacts_sha256']):
        raise ValueError('loose/archive artifact inventory differs')
    for name,digest in run['artifacts_sha256'].items():
        actual=packed[name]['sha256'] if name in packed else sha(native/name)
        if actual!=digest:raise ValueError('original artifact digest changed')
    if sha(previous['report_tsv'])!=previous['report_sha256']:raise ValueError('report changed')
    return dict(previous)


def main():
    verify_object()
    for pin in read(RECOVERY/'inputs/recovery-driver-pins.tsv')+read(RECOVERY/'inputs/execution-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('recovery input changed')
    allrows=[r for r in read(ORIGINAL/'inputs/execution-manifest.tsv') if r['layer']=='layer3']
    todo=read(RECOVERY/'inputs/execution-manifest.tsv');ids={r['genome_id'] for r in todo}
    if len(ids)!=len(todo):raise ValueError('duplicate recovery identities')
    for index in range((len(todo)+19)//20):
        worker=json.loads((RECOVERY/'logs'/f'layer3-worker-{index:04d}/finished.json').read_text())
        if worker['infrastructure_error'] or worker['incomplete_coded_ids'] or worker['unattempted_coded_ids']:
            raise ValueError('recovery worker incomplete')
    archived,prior=verify_archive()
    if ids & set(prior):raise ValueError('previously successful input rerun')
    if ids!={r['genome_id'] for r in allrows}-set(prior):raise ValueError('recovery set does not exactly cover missing public results')
    by_code={}
    for name,item in archived.items():
        parts=Path(name).parts
        if len(parts)<5 or parts[0]!='queries' or parts[2:4]!=('attempt01','analyze') or '..' in parts:
            raise ValueError('unexpected archive member location')
        by_code.setdefault(parts[1],{})['/'.join(parts[4:])]=item
    reason_rows=read(RECOVERY/'inputs/recovery-dispositions.tsv')
    reasons={r['genome_id']:r for r in reason_rows}
    if len(reasons)!=len(reason_rows) or set(reasons)!=ids:raise ValueError('disposition identities/counts differ')
    ledger=[]
    for row in allrows:
        ident=row['genome_id']
        if ident in ids:
            item=audit_query(row,RECOVERY)
            item.update(attempt_source='recovery',previous_attempt_state=reasons[ident]['previous_state'])
        else:
            item=archived_original(row,prior[ident],by_code[row['coded_id']])
            item.update(attempt_source='original',previous_attempt_state='successful')
        ledger.append(item)
    if len(ledger)!=4857:raise ValueError('public denominator changed')
    destination=ORIGINAL/'outputs/layer3-collection';destination.mkdir(exist_ok=False)
    write(destination/'attempt-ledger.tsv',ledger)
    fields=['genome_id','report_tsv','report_sha256']
    successful=[r for r in ledger if int(r['exit_status'])==0]
    write(destination/'results.tsv',[{k:r[k] for k in fields} for r in successful],fields)
    summary=dict(layer='layer3',genomes=len(ledger),successful=len(successful),failed=len(ledger)-len(successful),
        all_accounted=True,near_clones=dict(Counter(r['near_clone'] for r in ledger)),
        native_seconds=sum(float(r['seconds']) for r in ledger),production_assignments_enabled=False,
        infrastructure_recovery=dict(inputs=len(todo),previous_states=dict(Counter(r['previous_state'] for r in reasons.values())),
            successful_inputs_rerun=0,earlier_attempts_retained=True),archive_sha256=sha(ORIGINAL/'outputs/completed-native-ancillary.tar.gz'),
        collector_sha256=sha(__file__))
    json_write(destination/'summary.json',summary)
    json_write(RECOVERY/'outputs/merged-collection.json',summary)
    print(json.dumps(summary))

if __name__=='__main__':main()
