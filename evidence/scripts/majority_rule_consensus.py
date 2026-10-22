#!/usr/bin/env python3
"""Build a majority-rule consensus tree from replicate trees, and report what it contains.

**SET-19.** Task 10's "seed-stable" Mashtree comparator was not a consensus: it kept the seed-42
topology and deleted splits found in fewer than 19 of 20 seeds, so it could only ever *remove*.
Recomputing from the replicates showed the cost — 642 splits appeared in at least one replicate but
not in the seed-42 tree and were therefore unrecoverable by that filter, 25 of them in a majority of
the seeds and 3 in 19 of 20. A majority-rule consensus contains those and privileges no draw. This
script builds it.

What a majority-rule consensus is, and what it is not
----------------------------------------------------
It contains every split present in **more than half** the replicates. Such splits are guaranteed
pairwise compatible — two splits each in more than half the replicates must co-occur in at least one
— so they always form a tree. That guarantee is *verified here rather than assumed*, because a bug
in split extraction would otherwise produce a silently wrong tree rather than an error.

The frequencies written onto the internal branches are **replicate frequencies**. For a Mash seed
sweep they are **sketch stability** — how reproducible a split is when mash re-draws its
10,000-hash sketch — and SET-18 rules 2 and 6 forbid putting them in the same column as, or on the
same scale as, the marker tree's bootstrap or SH-like support. Nothing here calls them support.

Branch lengths are the **mean over the replicates that contain the split**, so the consensus keeps
units in which SET-18 rule 10's per-method near-zero threshold still means something. A split's
length is averaged only over the replicates that have it, never over all replicates, because
averaging in a zero for a replicate that lacks the split would invent a shorter branch.

Outputs
-------
``consensus.nwk``           the consensus, internal labels = replicate frequency
``consensus-counts.nwk``    the same topology, internal labels = "k/n"
``split-frequencies.tsv``   every split seen in any replicate, with its count and whether it is in
                            the consensus — including the ones a filter-the-reference approach
                            could never recover
``consensus-summary.json``  counts, the compatibility verification, and the replicate ledger
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import Tree, format_newick, parse_newick  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical(mask: int, full: int) -> int:
    """The side not containing tip 0: a genuine canonical form for any tip set."""
    return mask if not (mask & 1) else (full ^ mask)


def compatible(x: int, y: int, full: int) -> bool:
    """Compatible when one of the four quadrant intersections is empty."""
    xc, yc = full ^ x, full ^ y
    return not (x & y) or not (x & yc) or not (xc & y) or not (xc & yc)


def tip_names(tree: Tree) -> list[str]:
    return [tree.label[node] for node in tree.tips()]


def splits_of(tree: Tree, index: dict[str, int]) -> dict[int, float | None]:
    """Non-trivial unrooted splits, keyed by canonical bitmask, valued by branch length."""
    n = len(index)
    full = (1 << n) - 1
    mask = [0] * tree.size
    for node in tree.postorder():
        if not tree.children[node]:
            name = tree.label[node]
            if name not in index:
                raise SystemExit(f"tip {name!r} is not in the reference tip set")
            mask[node] = 1 << index[name]
        else:
            acc = 0
            for child in tree.children[node]:
                acc |= mask[child]
            mask[node] = acc
    found: dict[int, float | None] = {}
    for node in range(1, tree.size):
        if not tree.children[node]:
            continue
        size = bin(mask[node]).count("1")
        if size < 2 or size > n - 2:
            continue
        key = canonical(mask[node], full)
        # A tree can carry the same split twice only across a zero-length duplication; keep the
        # first length seen so a repeat cannot double-count.
        found.setdefault(key, tree.blen[node])
    return found


def terminal_lengths(tree: Tree) -> dict[str, float | None]:
    return {tree.label[node]: tree.blen[node] for node in tree.tips()}


def build_from_splits(order: list[str], masks: list[int], labels: dict[int, str],
                      lengths: dict[int, float | None],
                      tip_lengths: dict[str, float | None]) -> Tree:
    """Assemble a tree from a pairwise-compatible split set, largest clade first.

    Every tip starts as a child of the root, so the tree begins as a star. Splits are inserted in
    decreasing size, which guarantees that when a split arrives, the node containing all of its
    members already exists: its clade is formed by moving that node's children which lie wholly
    inside the split under one new node. If those children's union is not exactly the split, the
    split set was not compatible, and that is an error rather than a rounding of the topology.
    """
    n = len(order)
    full = (1 << n) - 1
    parent = [-1]
    children: list[list[int]] = [[]]
    label = [""]
    blen: list[float | None] = [None]
    node_mask = [full]

    def add(parent_index: int, mask: int, text: str, length: float | None) -> int:
        index = len(parent)
        parent.append(parent_index)
        children.append([])
        label.append(text)
        blen.append(length)
        node_mask.append(mask)
        children[parent_index].append(index)
        return index

    for position, name in enumerate(order):
        add(0, 1 << position, name, tip_lengths.get(name))

    for mask in sorted(masks, key=lambda m: -bin(m).count("1")):
        anchor, descending = 0, True
        while descending:
            descending = False
            for child in children[anchor]:
                if node_mask[child] & mask == mask:     # this child contains the whole split
                    anchor, descending = child, True
                    break
        movers = [child for child in children[anchor]
                  if node_mask[child] and node_mask[child] & mask == node_mask[child]]
        union = 0
        for child in movers:
            union |= node_mask[child]
        if union != mask:
            raise SystemExit(
                f"a split of {bin(mask).count('1')} tips does not sit under one common parent; "
                f"the split set is not pairwise compatible")
        if len(movers) == len(children[anchor]):
            # The split is already this node: nothing to insert, and inserting would create a
            # unifurcation. This happens only for a duplicate mask, which is filtered upstream.
            continue
        new_node = add(anchor, mask, labels.get(mask, ""), lengths.get(mask))
        for child in movers:
            children[anchor].remove(child)
            children[new_node].append(child)
            parent[child] = new_node
    return Tree(parent=parent, children=children, label=label, blen=blen)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicate", type=Path, action="append", required=True, default=[],
                        help="one replicate tree; repeat for each")
    parser.add_argument("--expected-tips", type=Path,
                        help="one tip id per line; every replicate must carry exactly this set")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="a split enters the consensus when its frequency EXCEEDS this "
                             "(default 0.5, the majority rule; compatibility is only guaranteed "
                             "at 0.5 or above and is verified either way)")
    parser.add_argument("--statistic-name", default="replicate_frequency",
                        help="what the internal labels mean; for a Mash seed sweep this is sketch "
                             "stability and is never branch support")
    parser.add_argument("--decimals", type=int, default=6)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    if len(args.replicate) < 2:
        raise SystemExit("a consensus needs at least two replicates")

    trees, ledger = [], []
    for path in args.replicate:
        tree = parse_newick(path.read_text())
        trees.append(tree)
        ledger.append({"file": str(path), "sha256": sha256_file(path),
                       "tips": len(tree.tips())})

    order = sorted(tip_names(trees[0]))
    if args.expected_tips:
        expected = sorted({line.strip() for line in args.expected_tips.read_text().splitlines()
                           if line.strip()})
        if order != expected:
            missing = sorted(set(expected) - set(order))
            extra = sorted(set(order) - set(expected))
            raise SystemExit(
                f"replicate tip set differs from {args.expected_tips}: "
                f"{len(missing)} missing {missing[:5]}, {len(extra)} unexpected {extra[:5]}")
    index = {name: position for position, name in enumerate(order)}
    for path, tree in zip(args.replicate, trees):
        if sorted(tip_names(tree)) != order:
            raise SystemExit(f"{path}: tip set differs from the first replicate")

    n = len(order)
    full = (1 << n) - 1
    counts: Counter = Counter()
    lengths: dict[int, list[float]] = defaultdict(list)
    tip_lengths: dict[str, list[float]] = defaultdict(list)
    for tree in trees:
        for mask, length in splits_of(tree, index).items():
            counts[mask] += 1
            if length is not None:
                lengths[mask].append(length)
        for name, length in terminal_lengths(tree).items():
            if length is not None:
                tip_lengths[name].append(length)

    replicates = len(trees)
    chosen = [mask for mask, count in counts.items() if count / replicates > args.threshold]

    # Verified, not assumed: a bug in split extraction would otherwise yield a silently wrong tree.
    incompatible = []
    for i, x in enumerate(chosen):
        for y in chosen[i + 1:]:
            if not compatible(x, y, full):
                incompatible.append((x, y))
    if incompatible:
        raise SystemExit(
            f"{len(incompatible)} incompatible pairs among {len(chosen)} majority splits. "
            f"Splits above half the replicates must be pairwise compatible, so this is a bug in "
            f"split extraction or a replicate with a different tip set — not a consensus.")

    mean = lambda values: sum(values) / len(values) if values else None  # noqa: E731
    labels = {mask: f"{counts[mask] / replicates:.{args.decimals}g}" for mask in chosen}
    count_labels = {mask: f"{counts[mask]}/{replicates}" for mask in chosen}
    branch_lengths = {mask: mean(lengths[mask]) for mask in chosen}
    mean_tip_lengths = {name: mean(tip_lengths[name]) for name in order}

    consensus = build_from_splits(order, chosen, labels, branch_lengths, mean_tip_lengths)
    counted = build_from_splits(order, chosen, count_labels, branch_lengths, mean_tip_lengths)

    # The consensus must contain exactly the splits chosen; a construction error would show here.
    rebuilt = set(splits_of(consensus, index))
    if rebuilt != set(chosen):
        raise SystemExit(
            f"the assembled consensus carries {len(rebuilt)} splits, {len(chosen)} were chosen "
            f"({len(set(chosen) - rebuilt)} lost, {len(rebuilt - set(chosen))} invented)")

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_name(prefix.name + ".nwk").write_text(
        format_newick(consensus, decimals=args.decimals) + "\n")
    prefix.with_name(prefix.name + "-counts.nwk").write_text(
        format_newick(counted, decimals=args.decimals) + "\n")

    rows = []
    for mask, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        size = bin(mask).count("1")
        rows.append({
            "split_rank": len(rows) + 1,
            "replicate_count": count,
            "replicate_frequency": round(count / replicates, args.decimals),
            "smaller_side_tips": min(size, n - size),
            "in_consensus": mask in set(chosen),
            "mean_branch_length": "" if not lengths[mask] else f"{mean(lengths[mask]):.10g}",
            "members": ";".join(order[position] for position in range(n) if mask & (1 << position)),
        })
    with (prefix.with_name(prefix.name + "-split-frequencies.tsv")).open(
            "w", encoding="utf-8", newline="") as handle:
        fields = ["split_rank", "replicate_count", "replicate_frequency", "smaller_side_tips",
                  "in_consensus", "mean_branch_length", "members"]
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[field]) for field in fields) + "\n")

    distribution = Counter(count for count in counts.values())
    summary = {
        "replicates": replicates,
        "tips": n,
        "threshold_exceeded": args.threshold,
        "statistic_name": args.statistic_name,
        "statistic_is_branch_support": False,
        "distinct_splits_across_replicates": len(counts),
        "splits_in_consensus": len(chosen),
        "splits_in_every_replicate": sum(1 for c in counts.values() if c == replicates),
        "splits_in_exactly_one_replicate": distribution.get(1, 0),
        "maximum_possible_splits": n - 3,
        "consensus_resolution_fraction": round(len(chosen) / (n - 3), 6) if n > 3 else None,
        "incompatible_majority_pairs": 0,
        "compatibility_verified": True,
        "reconstruction_verified": True,
        "replicate_count_distribution": dict(sorted(distribution.items())),
        "replicates_ledger": ledger,
        "outputs": {
            "consensus": str(prefix.with_name(prefix.name + ".nwk")),
            "consensus_counts": str(prefix.with_name(prefix.name + "-counts.nwk")),
            "split_frequencies": str(prefix.with_name(prefix.name + "-split-frequencies.tsv")),
        },
    }
    prefix.with_name(prefix.name + "-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items()
                      if key != "replicates_ledger"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
