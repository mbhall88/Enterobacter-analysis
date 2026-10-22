"""Pack ONLY this run's verified successful ancillary native outputs; retain all bytes."""
import sys,os,json,tarfile,hashlib,csv,shutil,time,datetime
from pathlib import Path
ROOT=Path('${PROJECT_ROOT}')
RUN=ROOT/'results/tool-planning/2026-10-07/validation-runs-astra'
LOCAL=Path('/tmp/enterobacter-validation-quota-recovery-2026-10-07')
sys.path.insert(0,str(RUN/'source'))
from common import read,sha,verify_object
from collect import audit_query
KEEP={'report.tsv','report.json','evidence.json','run.json','experimental-eligibility.json','mlst.json','genomovar-compatibility.json'}

def fsync_directory(path):
    fd=os.open(path,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)

def main():
    started=time.monotonic();verify_object();paths={};ledger=[]
    for pin in read(RUN/'inputs/execution-pins.tsv')+read(RUN/'inputs/driver-pins.tsv'):
        if sha(pin['path'])!=pin['sha256']:raise ValueError('execution input/driver pin changed')
    for i,row in enumerate(read(RUN/'inputs/execution-manifest.tsv')):
        directory=RUN/'queries'/row['coded_id']/'attempt01'
        if not (directory/'finished.json').exists():continue
        r=json.loads((directory/'finished.json').read_text())
        if r.get('exit_status')!=0:continue
        item=audit_query(row,RUN);ledger.append(item)
        native=directory/'analyze';manifest=json.loads((native/'run.json').read_text())['artifacts_sha256']
        for rel,digest in manifest.items():
            path=native/rel
            if rel in KEEP:continue
            if path.is_symlink():raise ValueError('unexpected symlink; preserve and inspect before packing')
            member=str(path.relative_to(RUN))
            if '..' in Path(member).parts or not path.resolve().is_relative_to(RUN.resolve()):
                raise ValueError('archive source outside owned run')
            paths[member]=dict(path=path,sha256=digest,bytes=path.stat().st_size)
        if len(ledger)%500==0:print('audited',len(ledger),'successful genomes',flush=True)
    with (LOCAL/'ancillary-inventory.tsv').open('x') as handle:
        writer=csv.DictWriter(handle,fieldnames=['member','sha256','bytes'],delimiter='\t',lineterminator='\n');writer.writeheader()
        for name,r in sorted(paths.items()):writer.writerow(dict(member=name,sha256=r['sha256'],bytes=r['bytes']))
    with (LOCAL/'successful-ledger.tsv').open('x') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(ledger[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(ledger)
    archive=LOCAL/'completed-native-ancillary.tar.gz'
    with tarfile.open(archive,'x:gz',compresslevel=1,dereference=True) as handle:
        for name,r in sorted(paths.items()):handle.add(r['path'],arcname=name,recursive=False)
    verified=set()
    with tarfile.open(archive,'r|gz') as handle:
        for member in handle:
            if not member.isfile() or member.name not in paths or member.name in verified:raise ValueError('archive inventory differs')
            reader=handle.extractfile(member)
            digest=hashlib.file_digest(reader,'sha256').hexdigest()
            if digest!=paths[member.name]['sha256']:raise ValueError('archive byte mismatch')
            verified.add(member.name)
    if verified!=set(paths):raise ValueError('archive missing files')
    archive_sha=sha(archive)
    destination=RUN/'outputs/completed-native-ancillary.tar.gz'
    checkpoint=dict(stage='fully_verified_before_any_unlink',time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        archive_path=str(archive),archive_sha256=archive_sha,files=len(paths),successful_genomes=len(ledger),
        inventory_sha256=sha(LOCAL/'ancillary-inventory.tsv'),ledger_sha256=sha(LOCAL/'successful-ledger.tsv'),
        source_sha256=sha(__file__),planned_destination=str(destination))
    with archive.open('rb') as handle:os.fsync(handle.fileno())
    for name in ['ancillary-inventory.tsv','successful-ledger.tsv']:
        with (LOCAL/name).open('rb') as handle:os.fsync(handle.fileno())
    with (LOCAL/'archive-verified-checkpoint.json').open('x') as handle:
        json.dump(checkpoint,handle,indent=2);handle.write('\n');handle.flush();os.fsync(handle.fileno())
    fsync_directory(LOCAL)
    print('archive verified',len(verified),'files',archive.stat().st_size,'compressed bytes',flush=True)
    # Free only enough entries to create a durable scratch copy. Every removed
    # byte is already held in the fully hash-verified local archive.
    ordered=sorted(paths)
    for name in ordered[:100]:
        p=paths[name]['path']
        if sha(p)!=paths[name]['sha256']:raise ValueError('source changed after packing')
        p.unlink()
    with destination.open('xb') as out,archive.open('rb') as inp:
        shutil.copyfileobj(inp,out,1024*1024);out.flush();os.fsync(out.fileno())
    if sha(destination)!=archive_sha:raise ValueError('scratch archive copy differs; stop, local archive retained')
    for source,target in [('ancillary-inventory.tsv','ancillary-inventory.tsv'),
                          ('successful-ledger.tsv','prearchive-successful-ledger.tsv'),
                          ('archive-verified-checkpoint.json','archive-verified-checkpoint.json')]:
        with (LOCAL/source).open('rb') as inp,(RUN/'outputs'/target).open('xb') as out:
            shutil.copyfileobj(inp,out);out.flush();os.fsync(out.fileno())
        if sha(LOCAL/source)!=sha(RUN/'outputs'/target):raise ValueError('archive provenance copy mismatch')
    fsync_directory(RUN/'outputs')
    # With the scratch archive verified, compact the remaining loose copies.
    for name in ordered[100:]:
        p=paths[name]['path']
        if sha(p)!=paths[name]['sha256']:raise ValueError('source changed after packing')
        p.unlink()
    for row in ledger:
        d=RUN/'queries'/row['coded_id']/'attempt01/analyze'
        for parent,dirs,files in os.walk(d,topdown=False):
            for child in dirs:
                candidate=Path(parent)/child
                if not any(candidate.iterdir()):candidate.rmdir()
    record=dict(status='complete',time=datetime.datetime.now(datetime.timezone.utc).isoformat(),successful_genomes=len(ledger),
        files_archived=len(paths),uncompressed_bytes=sum(r['bytes'] for r in paths.values()),
        compressed_bytes=destination.stat().st_size,archive_sha256=archive_sha,inventory_sha256=sha(RUN/'outputs/ancillary-inventory.tsv'),
        kept_native_files=sorted(KEEP),seconds=time.monotonic()-started,source_sha256=sha(__file__),rdm_touched=False,
        note='No outputs discarded. All archived member bytes checked against native manifests; source checked again before unlink. Restore to original paths by tar -xzf archive -C runroot.')
    for p in [LOCAL/'archive-result.json',RUN/'outputs/ancillary-archive.json']:
        with p.open('x') as out:json.dump(record,out,indent=2);out.write('\n')
    print(json.dumps(record),flush=True)

if __name__=='__main__':main()
