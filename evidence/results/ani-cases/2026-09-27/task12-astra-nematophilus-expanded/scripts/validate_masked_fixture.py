"""Audit saved deterministic synthetic alignment after production masking.

The raw-caller checks and their failures remain in validate_local_interface.py.
This test does not rerun or overwrite the alignment used to diagnose masking.
"""
from pathlib import Path
import json, random
from mask_pseudo import masked_sequences
from Bio import SeqIO

r = Path(__file__).resolve().parents[1]
d = r / 'interface-test-v2'
rng = random.Random(20260927)
reference = ''.join(rng.choices('ACGT', k=100000))
expected = list(reference)
for i in sorted(rng.sample(range(len(reference)), 6250)):
    expected[i] = rng.choice([b for b in 'ACGT' if b != reference[i]])
assert str(next(SeqIO.parse(d/'ref.fa', 'fasta')).seq) == reference
query = ''.join(expected[:40000] + expected[40100:])
query = query.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
assert str(next(SeqIO.parse(d/'query.fa', 'fasta')).seq) == query
observed = masked_sequences(d/'out')['ref']
assert len(observed) == len(reference)
called = [i for i, b in enumerate(observed) if b in 'ACGT']
assert len(called) > 98000
assert all(observed[i] == expected[i] for i in called)
assert set(observed[39990:40110]) == {'N'}
result = dict(status='PASS', callable=len(called), correct=len(called),
              reverse_orientation=True, indel_neighbourhood_masked=True,
              divergence=0.0625, saved_alignment_reused=True)
(d/'validation-masked-reproduced.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result))
