import io,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from Bio import Phylo
from task12_focused_scope_trees import splits
class SplitTests(unittest.TestCase):
 def test_reroot_invariant(self):
  t=Phylo.read(io.StringIO('((a:1,b:1):1,(c:1,(d:1,e:1):1):1);'),'newick');a=splits(t);t.root_with_outgroup('e');self.assertEqual(splits(t),a)
 def test_zero_edge_not_resolved(self):
  t=Phylo.read(io.StringIO('((a:1,b:1):0,(c:1,(d:1,e:1):0):0);'),'newick');self.assertEqual(splits(t,True),set());self.assertEqual(len(splits(t)),2)
 def test_different_quartet(self):
  a=Phylo.read(io.StringIO('((a:1,b:1):1,(c:1,d:1):1);'),'newick');b=Phylo.read(io.StringIO('((a:1,c:1):1,(b:1,d:1):1);'),'newick');self.assertEqual(len(splits(a)^splits(b)),2)
if __name__=='__main__':unittest.main()
