"""Synthetic-only saved-ANI adapter and subset tests."""
import json
import pytest
from common import PINS, load_module
from layer2 import CAL, masked_context, replay_row
from enterobacter_classifier.reporting import ReportingContext
from enterobacter_classifier.eligibility import primary_inventory


def synthetic_context():
    c=object.__new__(ReportingContext);c.package_sha256=PINS['package'][1]
    c._concepts={name:dict(concept_id=name,rank='species',display_name=name,nomenclature_status='synthetic',identifier_status='synthetic') for name in ['S1','S2']}
    c._references={}
    c._source_rows={}
    for rid,sp,block in [('r1','S1','b1'),('r2','S1','b1'),('r3','S2','b2')]:
        c._references[rid]=dict(reference_id=rid,accession=rid,species_id=sp,subspecies_id='',type_anchor=rid!='r2',search_label='/'+rid,source=dict(sha256='a'*64))
        c._source_rows[rid]=dict(clone_block=block,split_block=block)
    return c


def test_subset_recomputes_representation_without_mutating_context():
    context=synthetic_context();subset=masked_context(context,{'r2','r3'})
    inventory=primary_inventory(subset)
    assert set(context._references)=={'r1','r2','r3'}
    assert set(subset._references)=={'r2','r3'}
    s=inventory['representation']['S1']
    assert s['reference_ids']==['r2'] and s['reference_assemblies']==1 and s['type_anchor_assemblies']==0
    assert s['recorded_blocks']['clone_block']['distinct_recorded']==1


def test_saved_pairs_adapter_and_unmeasured_competitor():
    context=synthetic_context()
    refs=[dict(accession=r,whole_fasta='/'+r,whole_sha256='a'*64) for r in context._references]
    query=dict(whole_fasta='/q',whole_sha256='b'*64)
    def hit(r,ani):return dict(Query_file='/q',Ref_file='/'+r,ANI=str(ani),Align_fraction_ref='90',Align_fraction_query='90')
    rank=json.loads(PINS['rank'][0].read_text());config=json.loads(PINS['config'][0].read_text())
    report,result,_=replay_row([hit('r1',99),hit('r3',90)],query,refs,context,rank,config)
    assert report['hypothetical_species_id']=='S1' and result['production_enabled'] is False
    subset=masked_context(context,{'r1','r2'})
    report,result,_=replay_row([hit('r1',99)],query,refs[:2],subset,rank,config)
    assert not report['hypothetical_species_id']
    assert 'species_competitor_unmeasured' in report['experimental_eligibility_reasons']


def test_saved_reference_hash_mismatch_refused():
    context=synthetic_context()
    refs=[dict(accession=r,whole_fasta='/'+r,whole_sha256='c'*64) for r in context._references]
    with pytest.raises(ValueError,match='saved reference bytes differ'):
        replay_row([],dict(whole_fasta='/q',whole_sha256='b'*64),refs,context,{}, {})


def test_original_exclusion_contract_includes_aliases_and_blocks():
    old=load_module(CAL/'source/analyze.py','calibration_test_helpers')
    q=dict(accession='q',biosample='shared',entity='',clone_block='c1',split_block='s1',named_target='1',project_species='S1',project_subspecies='sub1',group='S1/sub1')
    r=dict(q,accession='r')
    assert old.exclude(q,r,'self_alias','represented')=='self_or_declared_alias'
    r.update(biosample='different',clone_block='c2')
    assert old.exclude(q,r,'self_alias','represented')==''
    assert old.exclude(q,r,'split_block','represented')=='shared_split_block'
    assert old.exclude(q,r,'self_alias','omit_reviewed_species')=='entire_reviewed_species'
    assert old.exclude(q,r,'self_alias','omit_reviewed_subspecies')=='entire_reviewed_subspecies'
