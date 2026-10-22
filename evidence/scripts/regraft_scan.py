#!/usr/bin/env python3
"""Every placement of one tip on an otherwise fixed tree, labelled by the hypotheses it satisfies.

For C3 the constrained ML searches wander: on the 53-genome core panel the search constrained to
the resolution the unconstrained tree already has ended 360 log-likelihood units below it, so a
constrained optimum is not a reliable representative of its hypothesis when the differences are
small.  The *E. hormaechei* subsp. *hoffmannii* type is a single-genome lineage, so the trio
question is where that one genome attaches.  This writes the unconstrained ML tree with the tip
pruned and regrafted onto every edge (the original placement included), as topologies without
branch lengths for IQ-TREE to optimise, and a table of which named hypotheses each carries.  IQ-TREE
then scores all of them on the same alignment and model (``--trees``), so the best placement under
each hypothesis differs from the others only in where that tip sits.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import Tree, parse_newick  # noqa: E402
from summarise_topology_tests import satisfies  # noqa: E402


def to_newick(tree: Tree, node: int = 0) -> str:
    if not tree.children[node]:
        return tree.label[node]
    return "(" + ",".join(to_newick(tree, c) for c in tree.children[node]) + ")"


def pruned(tree: Tree, tip_label: str) -> Tree:
    t = Tree(parent=list(tree.parent), children=[list(c) for c in tree.children],
             label=list(tree.label), blen=list(tree.blen))
    tip = next(n for n in t.tips() if t.label[n] == tip_label)
    parent = t.parent[tip]
    t.children[parent].remove(tip)
    t.parent[tip] = -1
    if len(t.children[parent]) == 1 and parent != 0:
        (only,) = t.children[parent]
        grand = t.parent[parent]
        t.children[grand] = [only if c == parent else c for c in t.children[grand]]
        t.parent[only] = grand
        t.children[parent] = []
        t.parent[parent] = -1
    elif len(t.children[0]) == 1:  # the root lost a child and became unary
        (only,) = t.children[0]
        t.children[0] = list(t.children[only])
        for c in t.children[0]:
            t.parent[c] = 0
        t.children[only] = []
        t.parent[only] = -1
    return t


def reachable(t: Tree) -> list[int]:
    out, stack = [], [0]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(t.children[n])
    return out


def regrafts(base: Tree, tip_label: str) -> list[str]:
    trees = []
    for edge in reachable(base):
        if edge == 0:
            continue
        t = Tree(parent=list(base.parent), children=[list(c) for c in base.children],
                 label=list(base.label), blen=list(base.blen))
        new = len(t.parent)
        leaf = new + 1
        t.parent += [t.parent[edge], new]
        t.children += [[edge, leaf], []]
        t.label += ["", tip_label]
        t.blen += [None, None]
        up = t.parent[edge]
        t.children[up] = [new if c == edge else c for c in t.children[up]]
        t.parent[edge] = new
        trees.append(to_newick(t) + ";")
    return trees


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--tip", required=True)
    parser.add_argument("--hypotheses", type=Path, required=True)
    parser.add_argument("--scope", choices=["subset", "panel"], required=True)
    parser.add_argument("--outgroups", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True, help="output prefix")
    args = parser.parse_args(argv)
    tree = parse_newick(args.tree.read_text().strip())
    hyps = json.loads(args.hypotheses.read_text())
    base = pruned(tree, args.tip)
    trees = regrafts(base, args.tip)
    original = to_newick(tree) + ";"
    rows = []
    for i, nw in enumerate(trees, 1):
        rows.append({"tree": i, "carries": [n for n, h in hyps.items()
                                            if satisfies(h, args.scope, nw, args.outgroups)]})
    tips = {len(parse_newick(nw).tips()) for nw in trees}
    if tips != {len(tree.tips())}:
        raise SystemExit(f"regrafted trees have {tips} tips, the source {len(tree.tips())}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{args.output}.trees").write_text("\n".join(trees) + "\n")
    summary = {"source_tree": str(args.tree), "tip": args.tip, "placements": len(trees),
               "per_hypothesis": {n: sum(1 for r in rows if n in r["carries"]) for n in hyps},
               "source_carries": [n for n, h in hyps.items()
                                  if satisfies(h, args.scope, original, args.outgroups)],
               "rows": rows}
    Path(f"{args.output}.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
