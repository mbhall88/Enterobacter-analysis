"""Fail closed on record selection and preservation, including a changed base."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from task12_focused_scope import additive
class AdditiveScopeTests(unittest.TestCase):
 def test_only_selected_added(self):
  original={'c':('c chromosome','ACGT')};whole={**original,'p':('p plasmid','TTTT'),'q':('q plasmid','GG')}
  result=additive(original,whole,'p');self.assertEqual(set(result),{'c','p'});self.assertEqual(result['c'],original['c']);self.assertEqual(len(original),1)
 def test_changed_retained_sequence_rejected(self):
  with self.assertRaises(AssertionError):additive({'c':('c','ACGT')},{'c':('c','ACGA'),'p':('p','TT')},'p')
 def test_already_retained_rejected(self):
  with self.assertRaises(AssertionError):additive({'c':('c','ACGT')},{'c':('c','ACGT')},'c')
 def test_absent_record_rejected(self):
  with self.assertRaises(AssertionError):additive({'c':('c','ACGT')},{'c':('c','ACGT')},'p')
if __name__=='__main__':unittest.main()

from task12_focused_scope_validate import covered, delta
import tempfile
class IndependentDeltaTests(unittest.TestCase):
 def test_overlap_reverse_adjacent(self):
  self.assertEqual(covered([(5,1),(3,7),(10,8),(12,12)]),11)
 def test_contigs_separate_identity_and_union(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'x.delta';p.write_text('ref query\nNUCMER\n>r1 q1 20 20\n1 10 11 1 5 5 0\n0\n5 10 5 10 0 0 0\n0\n>r2 q2 10 10\n1 6 6 1 0 0 0\n0\n')
   r,q,ani,e=delta(p);self.assertEqual((r,q,e),(16,17,5));self.assertAlmostEqual(ani,100*(1-10/45))
