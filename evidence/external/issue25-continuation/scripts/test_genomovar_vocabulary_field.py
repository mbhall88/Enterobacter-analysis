import pathlib,json,datetime
import evaluation_common as c
import p3_adapter as p
rows=[(key,value) for key,value in c.VOCAB.items() if key[0]=='p3a-genomovar-contract-r1' and key[1]=='genomovar.reasons'];assert rows
for (policy,field,reason),mapped in rows:assert p.translate_genomovar_reasons(policy,[reason])==[mapped]
try:p.translate_genomovar_reasons('p3a-genomovar-contract-r1',['UNKNOWN_SYNTHETIC_REASON'])
except KeyError:pass
else:raise AssertionError('Unmapped reason silently accepted')
out=c.CONT/'checks/genomovar-vocabulary-field-fixtures.json';assert not out.exists();out.write_text(json.dumps({'gate_pass':True,'mapped_reasons':len(rows),'unknown_reason_rejected':True,'exact_policy_and_field':True,'source_sha256':c.sha(__file__),'adapter_sha256':c.sha(c.CONT/'scripts/p3_adapter.py'),'vocabulary_sha256':c.sha(c.RUN/'build/source/docs/contracts/reason-vocabulary.tsv'),'biological_scoring':False,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n');print('Genomovar vocabulary field fixtures PASS',len(rows))
