"""Read unrooted splits and IQ-TREE AU tables without treating a root as evidence."""
def tips(tree):return {x.name for x in tree.get_terminals()}
def key(side,all_tips):return min(tuple(sorted(side)),tuple(sorted(all_tips-set(side))))
def splits(tree):
 all_tips=tips(tree);out={}
 for n in tree.get_nonterminals():
  s=tips(n)
  if min(len(s),len(all_tips-s))>=2:out[key(s,all_tips)]=n
 return out

def au_table(text):
 section=text.split('USER TREES',1)[1];out={}
 for line in section.splitlines():
  t=line.split()
  if len(t)==13 and t[0].isdigit() and t[-1] in ['+','-']:
   out[int(t[0])]=dict(logL=float(t[1]),deltaL=float(t[2]),p_AU=float(t[-2]),raw_row=line.strip())
 if not out:raise ValueError('No AU rows parsed')
 return out
