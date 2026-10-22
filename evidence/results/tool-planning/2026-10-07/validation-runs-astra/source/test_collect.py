"""Synthetic collector checks; no validation sequences are analysed."""
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from common import sha, json_write, write
from run_batch import execute
from collect import audit_query


def fixture(tmp_path):
    query=tmp_path/'query.fa'; query.write_text('>synthetic\nACGT\n')
    row=dict(genome_id='synthetic',coded_id='TEST-001',layer='layer1',fasta_path=str(query),query_sha256=sha(query))
    def fake(argv, **kwargs):
        out=Path(argv[argv.index('--output')+1]); out.mkdir()
        write(out/'report.tsv',[dict(query_sha256=row['query_sha256'],assigned_species_id='',assigned_subspecies_id='')])
        json_write(out/'evidence.json',dict(ranked_reference_evidence=[dict(ani_percent=99.9,aligned_fraction_query_percent=90,aligned_fraction_reference_percent=90)]))
        json_write(out/'run.json',dict(status='complete',positive_calling_policy_enabled=False,child_peak_rss_kib=123,
            artifacts_sha256={p.name:sha(p) for p in out.iterdir()}))
        return SimpleNamespace(returncode=0)
    execute(row,tmp_path,fake)
    return row


def test_valid_report_and_inclusive_near_clone(tmp_path):
    row=fixture(tmp_path)
    item=audit_query(row,tmp_path)
    assert item['near_clone']=='True' and item['exit_status']==0


def test_artifact_change_refused(tmp_path):
    row=fixture(tmp_path)
    (tmp_path/'queries/TEST-001/attempt01/analyze/evidence.json').write_text('{}')
    with pytest.raises(ValueError,match='native artifact changed'):
        audit_query(row,tmp_path)


def test_missing_native_manifest_member_refused(tmp_path):
    row=fixture(tmp_path)
    p=tmp_path/'queries/TEST-001/attempt01/analyze/run.json'
    data=json.loads(p.read_text());data['artifacts_sha256'].pop('evidence.json');p.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='native artifact manifest incomplete'):
        audit_query(row,tmp_path)


def test_unattempted_genome_not_silently_failed(tmp_path):
    row=fixture(tmp_path);row['coded_id']='TEST-002'
    with pytest.raises(FileNotFoundError):audit_query(row,tmp_path)


def test_nonzero_partial_report_not_a_result(tmp_path):
    row=fixture(tmp_path)
    p=tmp_path/'queries/TEST-001/attempt01/finished.json'
    data=json.loads(p.read_text());data['exit_status']=2;p.write_text(json.dumps(data))
    result=audit_query(row,tmp_path)
    assert result['exit_status']==2 and result['report_tsv']=='' and result['report_sha256']==''
