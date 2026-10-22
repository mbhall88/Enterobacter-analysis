"""One attempt per frozen input. Infrastructure reruns require separate authorisation records."""
import argparse
import datetime
import json
import os
import socket
import subprocess
import time
from pathlib import Path

from common import ROOT, RUN, batch, command, json_write, read, sha, verify_object


def execute(row, root, runner=subprocess.run):
    output = root / 'queries' / row['coded_id'] / 'attempt01'
    output.mkdir(parents=True, exist_ok=False)  # A duplicate job cannot overwrite/rerun an input.
    fasta = Path(row['fasta_path'])
    argv = command(fasta, output / 'analyze')
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    initial = dict(coded_id=row['coded_id'], genome_id=row['genome_id'], query_sha256=row['query_sha256'],
                   argv=argv, cwd=str(ROOT), started=started, hostname=socket.gethostname(),
                   slurm_job_id=os.environ.get('SLURM_JOB_ID'), slurm_array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'),
                   source_verification='worker start/end plus preflight', production_assignments_enabled=False)
    json_write(output / 'started.json', initial)
    env = dict(os.environ, PYTHONPATH=str(ROOT / 'src'), PYTHONDONTWRITEBYTECODE='1',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    timer = time.monotonic()
    result = dict(initial, exit_status=None, process_started=False, query_unchanged=None,
                  report_tsv='', report_sha256='', infrastructure_error='')
    try:
        if sha(fasta) != row['query_sha256']:
            raise ValueError(f"query hash changed: {row['coded_id']}")
        with (output / 'stdout.txt').open('x') as out, (output / 'stderr.txt').open('x') as err:
            completed = runner(['/usr/bin/time', '-v', '-o', str(output / 'resources.txt'), *argv],
                               cwd=ROOT, env=env, stdout=out, stderr=err, check=False)
        result.update(process_started=True, exit_status=completed.returncode)
        result['query_unchanged'] = sha(fasta) == row['query_sha256']
        report = output / 'analyze/report.tsv'
        if report.is_file():
            result.update(report_tsv=str(report), report_sha256=sha(report))
        if not result['query_unchanged']:
            raise ValueError('query changed during execution; outputs quarantined')
        if completed.returncode == 0 and not report.is_file():
            raise ValueError('zero exit without report; review infrastructure before proceeding')
    except Exception as error:
        result['infrastructure_error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        result.update(finished=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      seconds=time.monotonic() - timer)
        json_write(output / 'finished.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--layer', choices=['layer1', 'layer3'], required=True)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--chunk-size', type=int, default=20)
    args = p.parse_args()
    log = RUN / 'logs' / f'{args.layer}-worker-{args.index:04d}'
    log.mkdir(exist_ok=False)
    check = verify_object()
    pins = read(RUN / 'inputs/execution-pins.tsv')
    # Input metadata and driver pins are cheap; no repeated hashing of tree tables here.
    for record in pins:
        if sha(record['path']) != record['sha256']:
            raise ValueError('frozen execution input changed: ' + record['path'])
    for record in read(RUN / 'inputs/driver-pins.tsv'):
        if sha(record['path']) != record['sha256']:
            raise ValueError('execution driver changed before batch')
    rows = [r for r in read(RUN / 'inputs/execution-manifest.tsv') if r['layer'] == args.layer]
    selected = batch(rows, args.index, args.chunk_size)
    json_write(log / 'started.json', dict(layer=args.layer, index=args.index,
               coded_ids=[r['coded_id'] for r in selected], verification=check))
    results, error_text = [], ''
    try:
        for row in selected:
            results.append(execute(row, RUN))
        verify_object()
    except Exception as error:
        error_text = f'{type(error).__name__}: {error}'
        raise
    finally:
        json_write(log / 'finished.json', dict(layer=args.layer, index=args.index, completed=len(results),
                   successful=sum(r['exit_status'] == 0 for r in results),
                   failed=sum(r['exit_status'] != 0 for r in results), infrastructure_error=error_text,
                   incomplete_coded_ids=[r['coded_id'] for r in selected[len(results):]],
                   unattempted_coded_ids=[r['coded_id'] for r in selected
                     if not (RUN / 'queries' / r['coded_id'] / 'attempt01/started.json').exists()]))
    print(json.dumps(dict(layer=args.layer, batch=args.index, attempted=len(results),
                         failed=sum(r['exit_status'] != 0 for r in results))))


if __name__ == '__main__':
    main()
