"""Synthetic execution tests; never call a biological program on validation inputs."""
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pytest

from common import batch, command, sha
from prepare import coded
from run_batch import execute


@pytest.mark.parametrize('count', [946, 4857])
def test_array_partitions_exactly_once(count):
    rows = list(range(count))
    chunks = [batch(rows, i, 20) for i in range((count + 19) // 20)]
    assert [r for chunk in chunks for r in chunk] == rows
    assert all(1 <= len(chunk) <= 20 for chunk in chunks)
    with pytest.raises(ValueError):
        batch(rows, len(chunks), 20)


def test_frozen_command_contract():
    argv = command(Path('/synthetic/query.fa'), Path('/synthetic/output'))
    assert argv[1:4] == ['-m', 'enterobacter_classifier', 'analyze']
    flags = [x for x in argv if x.startswith('--')]
    assert set(flags) == {'--reference-package', '--expected-package-sha256', '--experimental-rank-rule',
                         '--experimental-eligibility-config', '--expected-eligibility-config-sha256',
                         '--experimental-genomovars', '--mlst-config', '--expected-mlst-config-sha256',
                         '--skani-binary', '--threads', '--query', '--output'}
    assert len(flags) == len(set(flags))
    assert argv[argv.index('--threads') + 1] == '2'
    assert not any('checkm2' in x or 'quality-evidence' in x for x in argv)


def test_coded_ids_are_deterministic_and_hide_accessions():
    code = coded('synthetic-accession', 'l1')
    assert code == coded('synthetic-accession', 'l1')
    assert code != coded('different-accession', 'l1')
    assert re.fullmatch(r'VAL-L1-[A-F0-9]{12}', code)


def row_for(tmp_path):
    query = tmp_path / 'fixture.fa'
    query.write_text('>synthetic\nACGT\n')
    return dict(genome_id='synthetic-only', coded_id='TEST-001', fasta_path=str(query), query_sha256=sha(query))


def test_success_preserves_exact_command_and_refuses_duplicate(tmp_path):
    row = row_for(tmp_path)
    calls = []
    def fake(argv, **kwargs):
        calls.append(argv)
        output = Path(argv[argv.index('--output') + 1])
        output.mkdir()
        (output / 'report.tsv').write_text('synthetic\nonly\n')
        assert kwargs['env']['PYTHONDONTWRITEBYTECODE'] == '1'
        return SimpleNamespace(returncode=0)
    result = execute(row, tmp_path, fake)
    assert result['exit_status'] == 0 and result['query_unchanged']
    assert result['report_sha256'] == sha(result['report_tsv'])
    assert result['argv'] == calls[0][4:]
    with pytest.raises(FileExistsError):
        execute(row, tmp_path, fake)
    assert len(calls) == 1


def test_failed_run_is_recorded_without_retry(tmp_path):
    row = row_for(tmp_path)
    result = execute(row, tmp_path, lambda *a, **k: SimpleNamespace(returncode=2))
    assert result['exit_status'] == 2 and result['report_tsv'] == ''
    assert json.loads((tmp_path / 'queries/TEST-001/attempt01/finished.json').read_text()) == result


def test_input_hash_mismatch_stops_before_analysis(tmp_path):
    row = row_for(tmp_path)
    Path(row['fasta_path']).write_text('>changed\nTTTT\n')
    def forbidden(*a, **k):
        pytest.fail('must not run a changed input')
    with pytest.raises(ValueError, match='query hash changed'):
        execute(row, tmp_path, forbidden)
    record = json.loads((tmp_path / 'queries/TEST-001/attempt01/finished.json').read_text())
    assert not record['process_started'] and record['exit_status'] is None
    assert 'query hash changed' in record['infrastructure_error']


def test_launch_exception_records_outcome_and_preserves_duplicate_guard(tmp_path):
    row = row_for(tmp_path)
    def unavailable(*a, **k):
        raise FileNotFoundError('synthetic missing executable')
    with pytest.raises(FileNotFoundError):
        execute(row, tmp_path, unavailable)
    record = json.loads((tmp_path / 'queries/TEST-001/attempt01/finished.json').read_text())
    assert record['exit_status'] is None and record['seconds'] >= 0
    assert not record['process_started'] and 'FileNotFoundError' in record['infrastructure_error']
    with pytest.raises(FileExistsError):
        execute(row, tmp_path, unavailable)


def test_input_mutation_during_execution_is_quarantined(tmp_path):
    row = row_for(tmp_path)
    def fake(*a, **k):
        Path(row['fasta_path']).write_text('changed during run')
        return SimpleNamespace(returncode=2)
    with pytest.raises(ValueError, match='query changed during execution'):
        execute(row, tmp_path, fake)
    record = json.loads((tmp_path / 'queries/TEST-001/attempt01/finished.json').read_text())
    assert record['query_unchanged'] is False


def test_missing_report_is_not_silent_success(tmp_path):
    with pytest.raises(ValueError, match='zero exit without report'):
        execute(row_for(tmp_path), tmp_path, lambda *a, **k: SimpleNamespace(returncode=0))


def test_post_run_missing_input_retains_actual_child_exit(tmp_path):
    row = row_for(tmp_path)
    def fake(*args, **kwargs):
        Path(row['fasta_path']).unlink()
        return SimpleNamespace(returncode=2)
    with pytest.raises(FileNotFoundError):
        execute(row, tmp_path, fake)
    record = json.loads((tmp_path / 'queries/TEST-001/attempt01/finished.json').read_text())
    assert record['exit_status'] == 2 and record['process_started']
    assert 'FileNotFoundError' in record['infrastructure_error']
