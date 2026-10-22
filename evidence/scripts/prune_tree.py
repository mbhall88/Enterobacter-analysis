#!/usr/bin/env python3
"""Prune a tree to a stated set of tips, preserving patristic distances.

Task 06 declared that the 13 `declared_out_of_genus_context` tips are drawn in the tree but are
excluded from the rooting argument, which is interpreted on the remaining 375. Rooting therefore
has to be tested on that pruned tree, not on all 388 -- with the 13 present, the complement of the
two outgroups is not a clade and the question cannot even be asked.

Pruning is Task 05's tested implementation (`scripts/gtdb_context.py`), which sums the branch
lengths of suppressed nodes into the surviving branch and returns collapsed internal labels rather
than dropping them silently.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import format_newick, parse_newick, prune_to_tips  # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--keep", required=True, type=Path,
                        help="file of tip ids to keep, one per line")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--decimals", type=int, default=6,
                        help="decimal places for branch lengths in the output (default 6, kept "
                             "so earlier command ledgers reproduce byte-for-byte). Distance trees "
                             "need more: Verticall/FastME branches go down to ~1e-08, and the "
                             "script refuses rather than round a non-zero branch to zero")
    parser.add_argument("--allow-missing", action="store_true",
                        help="keep the intersection of --keep and the tree's tips instead of "
                             "failing. Needed for a sensitivity tree that deliberately omits some "
                             "of the scope; the omitted tips are still reported.")
    args = parser.parse_args(argv)

    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")
    tree = parse_newick(args.tree.read_text())
    keep_names = {l.strip() for l in args.keep.read_text().splitlines() if l.strip()}
    by_name = {tree.label[i]: i for i in tree.tips()}
    missing = keep_names - set(by_name)
    if missing and not args.allow_missing:
        raise SystemExit(f"{len(missing)} requested tips are not in the tree: "
                         f"{sorted(missing)[:10]}")
    keep_names -= missing
    if len(keep_names) < 3:
        raise SystemExit(f"only {len(keep_names)} of the requested tips are in the tree")
    result = prune_to_tips(tree, {by_name[n] for n in keep_names})
    # Rounding a branch to zero changes the tree's meaning, not just its precision: every
    # comparison in this project collapses branches at or below a near-zero threshold, and a
    # rounded-away branch becomes a collapsed one. So it is an error, not a formatting choice.
    zeroed = [length for length in result.tree.blen
              if length is not None and length != 0 and round(abs(length), args.decimals) == 0]
    if zeroed:
        raise SystemExit(
            f"{len(zeroed)} non-zero branch lengths would round to zero at --decimals "
            f"{args.decimals} (smallest {min(abs(x) for x in zeroed):.3e}); raise --decimals")
    args.output.write_text(format_newick(result.tree, decimals=args.decimals) + "\n")
    report = {
        "input_tips": len(by_name),
        "kept_tips": len(keep_names),
        "dropped_tips": sorted(set(by_name) - keep_names),
        "requested_tips_absent_from_this_tree": sorted(missing),
        "internal_labels_collapsed": sum(len(v) for v in result.collapsed_labels.values()),
        "output": str(args.output),
    }
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
