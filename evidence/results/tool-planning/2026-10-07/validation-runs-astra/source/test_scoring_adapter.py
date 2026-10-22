"""Synthetic-only manifest selection: never lose failed genomes from frozen truth."""
import pytest
from score_layers import select_results
from score_layers import logged_process
import json
import score_layers
from types import SimpleNamespace


def test_failures_omitted_from_results_but_truth_untouched():
    truth=[dict(genome_id='one'),dict(genome_id='two')]
    ledger=[dict(genome_id='one',exit_status='0',report_tsv='synthetic-report',report_sha256='a'*64),
            dict(genome_id='two',exit_status='2',report_tsv='partial-report',report_sha256='b'*64)]
    assert [r['genome_id'] for r in select_results(truth,ledger)]==['one']
    assert len(truth)==2


def test_missing_ledger_and_duplicate_ledger_fail():
    truth=[dict(genome_id='one')]
    with pytest.raises(KeyError):select_results(truth,[])
    row=dict(genome_id='one',exit_status='0',report_tsv='file',report_sha256='a'*64)
    with pytest.raises(ValueError,match='duplicate'):select_results(truth,[row,row])


def test_subprocess_launch_failure_has_durable_record(tmp_path,monkeypatch):
    def unavailable(*args,**kwargs):raise FileNotFoundError('synthetic')
    monkeypatch.setattr(score_layers.subprocess,'run',unavailable)
    with pytest.raises(FileNotFoundError):logged_process(['synthetic'],tmp_path)
    result=json.loads((tmp_path/'finished.json').read_text())
    assert result['exit_status'] is None and result['seconds']>=0
    assert 'FileNotFoundError' in result['infrastructure_error']


def test_original_names_split_combined_sampling_group():
    truth=[dict(genome_id='a',sample_group='unnamed'),dict(genome_id='b',sample_group='unnamed')]
    ledger=[dict(genome_id=x,exit_status='2') for x in ['a','b']]
    metadata={'a':dict(ncbi_name_normalised='Enterobacter sp.'),
              'b':dict(ncbi_name_normalised='Enterobacter cloacae complex sp.')}
    scorer=SimpleNamespace(score=lambda t,r:t,summarise=lambda rows:len(rows))
    result=score_layers.expert_name_breakdown(truth,ledger,metadata,scorer)
    assert result=={'Enterobacter sp.':1,'Enterobacter cloacae complex sp.':1}
    assert all(t['sample_group']=='unnamed' for t in truth)
