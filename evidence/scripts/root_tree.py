#!/usr/bin/env python3
"""Root an unrooted ML tree on reviewed outgroups, and measure whether the root is stable.

Task 06 selected two outgroups from *different genera* -- *Klebsiella aerogenes* and *Cronobacter
sakazakii* -- precisely so that rooting stability can be **shown** rather than assumed.  This tool
therefore does three rootings, not one: both outgroups together, and each outgroup alone.  If the
two outgroups are not sister on the inferred tree, the three roots land on different edges and that
disagreement is reported as a result rather than hidden by picking one.

Rooting never changes the unrooted topology; it only chooses an edge.  So the comparison that
matters is *which edge the root sits on* and whether the ingroup clade set is identical under each
rooting.  Both are reported.

Support values in a Newick internal label annotate the **branch** below that node, not the node, so
this reroots on an edge-labelled graph: when the root is placed on an edge, that edge's support is
carried to both halves and the duplication is recorded.  Rerooting on a node-labelled array would
silently shift support values onto the wrong branches.

The Newick parser, formatter and tree container are Task 05's (``scripts/gtdb_context.py``); this
adds only the rerooting and the stability comparison.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import Tree, format_newick, parse_newick  # noqa: E402


class RootingError(RuntimeError):
    pass


# --------------------------------------------------------------------------------------
# Undirected edge-labelled view
# --------------------------------------------------------------------------------------


def to_graph(tree: Tree):
    """Return (adjacency, tip_name_of_node).

    ``adjacency[u]`` is a list of ``(v, length, support)``.  ``support`` is the label that
    annotated the edge in the input -- for an edge from a parent to an internal child that is the
    child's Newick label, which is where IQ-TREE writes ``aLRT/UFboot``.
    """
    adjacency: dict[int, list[tuple[int, float, str]]] = {i: [] for i in range(tree.size)}
    tip_name: dict[int, str] = {}
    for node in range(tree.size):
        if not tree.children[node]:
            tip_name[node] = tree.label[node]
    for node in range(tree.size):
        parent = tree.parent[node]
        if parent == -1:
            continue
        length = tree.blen[node] or 0.0
        support = tree.label[node] if tree.children[node] else ""
        adjacency[parent].append((node, length, support))
        adjacency[node].append((parent, length, support))
    return adjacency, tip_name


def tip_partitions(adjacency, tip_name):
    """For each undirected edge (u, v), the set of tip names on v's side when u is removed."""
    tips = set(tip_name.values())
    memo: dict[tuple[int, int], frozenset] = {}

    def side(u: int, v: int) -> frozenset:
        """Tips reachable from v without going back through u."""
        key = (u, v)
        if key in memo:
            return memo[key]
        # Iterative DFS; the ML trees here are deep enough that recursion is a real risk.
        found = set()
        stack = [(u, v)]
        seen = {(u, v)}
        while stack:
            prev, cur = stack.pop()
            if cur in tip_name:
                found.add(tip_name[cur])
            for nxt, _, _ in adjacency[cur]:
                if nxt == prev:
                    continue
                if (cur, nxt) not in seen:
                    seen.add((cur, nxt))
                    stack.append((cur, nxt))
        memo[key] = frozenset(found)
        return memo[key]

    return side, tips


def find_root_edge(adjacency, tip_name, outgroup: set[str]):
    """Find the unique edge whose removal separates exactly ``outgroup`` from the rest.

    Returns (u, v, length, support) with v on the outgroup side, or raises if the outgroup is not
    a clade of the unrooted tree.
    """
    side, tips = tip_partitions(adjacency, tip_name)
    missing = outgroup - tips
    if missing:
        raise RootingError(f"outgroup tips absent from the tree: {sorted(missing)}")
    ingroup = tips - outgroup
    if not ingroup:
        raise RootingError("the outgroup is every tip in the tree")
    for u in adjacency:
        for v, length, support in adjacency[u]:
            if side(u, v) == outgroup:
                return u, v, length, support
    raise RootingError(
        "the outgroup is not monophyletic on this tree: no single edge separates "
        f"{sorted(outgroup)} from the remaining {len(ingroup)} tips"
    )


def reroot(adjacency, tip_name, u: int, v: int, length: float, support: str,
           position: float = 0.5) -> Tree:
    """Place a new root on the edge (u, v) and rebuild a rooted :class:`Tree`."""
    if not 0.0 <= position <= 1.0:
        raise RootingError("position must be between 0 and 1")
    new_root = max(adjacency) + 1
    local = {k: [e for e in val] for k, val in adjacency.items()}
    local[new_root] = []
    local[u] = [e for e in local[u] if e[0] != v]
    local[v] = [e for e in local[v] if e[0] != u]
    # The edge's support annotates the bipartition, which both halves still express.
    local[new_root].append((u, length * (1.0 - position), support))
    local[u].append((new_root, length * (1.0 - position), support))
    local[new_root].append((v, length * position, support))
    local[v].append((new_root, length * position, support))

    parent: list[int] = []
    children: list[list[int]] = []
    label: list[str] = []
    blen: list[float | None] = []
    index: dict[int, int] = {}

    def add(node: int, parent_index: int, length_: float | None, support_: str) -> int:
        idx = len(parent)
        index[node] = idx
        parent.append(parent_index)
        children.append([])
        label.append(tip_name.get(node, support_))
        blen.append(length_)
        if parent_index != -1:
            children[parent_index].append(idx)
        return idx

    root_index = add(new_root, -1, None, "")
    stack = [(new_root, root_index)]
    while stack:
        node, idx = stack.pop()
        for nxt, length_, support_ in local[node]:
            if nxt in index and index[nxt] < idx:
                continue
            if nxt in index:
                continue
            child_idx = add(nxt, idx, length_, support_)
            stack.append((nxt, child_idx))
    return Tree(parent=parent, children=children, label=label, blen=blen)


def clade_set(tree: Tree, restrict_to: set[str] | None = None) -> set[frozenset]:
    """Every non-trivial clade of a rooted tree, as a set of tip-name frozensets."""
    tips_below: dict[int, set[str]] = {}
    clades: set[frozenset] = set()
    for node in tree.postorder():
        if not tree.children[node]:
            name = tree.label[node]
            tips_below[node] = {name} if (restrict_to is None or name in restrict_to) else set()
        else:
            merged: set[str] = set()
            for child in tree.children[node]:
                merged |= tips_below[child]
            tips_below[node] = merged
        if len(tips_below[node]) > 1:
            clades.add(frozenset(tips_below[node]))
    return clades


# --------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------


def root_and_compare(newick: str, outgroups: dict[str, list[str]]) -> dict:
    """``outgroups`` maps a rooting name to its tip ids; all are attempted independently."""
    tree = parse_newick(newick)
    adjacency, tip_name = to_graph(tree)
    all_tips = set(tip_name.values())

    results = {}
    rooted_trees = {}
    for name, tips in outgroups.items():
        entry: dict = {"outgroup_tips": sorted(tips)}
        try:
            u, v, length, support = find_root_edge(adjacency, tip_name, set(tips))
        except RootingError as exc:
            entry["rooted"] = False
            entry["reason"] = str(exc)
            results[name] = entry
            continue
        rooted = reroot(adjacency, tip_name, u, v, length, support)
        side, _ = tip_partitions(adjacency, tip_name)
        entry.update({
            "rooted": True,
            "root_edge_length": length,
            "root_edge_support": support,
            "outgroup_side_tip_count": len(side(u, v)),
            "ingroup_side_tip_count": len(all_tips) - len(side(u, v)),
        })
        results[name] = entry
        rooted_trees[name] = rooted

    # Stability: do the rootings agree on the ingroup clade set?
    reference = None
    for name in outgroups:
        if name in rooted_trees:
            reference = name
            break
    stability: dict = {"comparable_rootings": sorted(rooted_trees)}
    if reference and len(rooted_trees) > 1:
        union_outgroup = set().union(*(set(t) for t in outgroups.values()))
        ingroup_only = all_tips - union_outgroup
        ref_clades = clade_set(rooted_trees[reference], ingroup_only)
        for name, rooted in rooted_trees.items():
            if name == reference:
                continue
            other = clade_set(rooted, ingroup_only)
            stability[f"{reference}_vs_{name}"] = {
                "ingroup_clades_identical": ref_clades == other,
                "clades_only_in_reference": len(ref_clades - other),
                "clades_only_in_other": len(other - ref_clades),
            }
    return {"rootings": results, "stability": stability, "trees": rooted_trees,
            "tip_count": len(all_tips)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", required=True, type=Path, help="unrooted Newick tree")
    parser.add_argument("--outgroup", action="append", required=True, metavar="NAME=TIP[,TIP...]",
                        help="a named rooting; repeat for each rooting to compare")
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args(argv)

    outgroups: dict[str, list[str]] = {}
    for spec in args.outgroup:
        if "=" not in spec:
            parser.error(f"--outgroup needs NAME=TIP[,TIP...]; got {spec!r}")
        name, _, tips = spec.partition("=")
        outgroups[name] = [t for t in tips.split(",") if t]
        if not outgroups[name]:
            parser.error(f"--outgroup {name} lists no tips")

    result = root_and_compare(args.tree.read_text(), outgroups)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, rooted in result.pop("trees").items():
        path = Path(f"{args.output_prefix}.{name}.rooted.nwk")
        path.write_text(format_newick(rooted) + "\n")
        result["rootings"][name]["tree"] = str(path)
    report = Path(f"{args.output_prefix}.rooting.json")
    report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
