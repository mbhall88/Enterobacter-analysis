#!/usr/bin/env python3
"""Read one run_topology_tests.sh directory: which resolution each candidate tree carries, and
whether the data reject it.

The AU, SH, KH, weighted KH/SH and ELW columns come from IQ-TREE's ``au.iqtree`` table as
written.  Each candidate is labelled by the hypotheses it satisfies, found by tip set: a clade
hypothesis is satisfied when the tree has its bipartition (ignoring free tips); a quartet
hypothesis when the tree, restricted to the two paired types, the third type and the outgroups,
puts the pair together.  A candidate that satisfies none of the named hypotheses is labelled so.

A p-value here is a topology-test statistic: rejection at 0.05 means the alignment rejects that
resolution under the stated model.  It is not a support value (SET-18 rule 2) and never shares a
column with one.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from summarise_concordance import find_branch, normalise, tree_splits  # noqa: E402

COLUMNS = ["logL", "deltaL", "bp-RELL", "p-KH", "p-SH", "p-WKH", "p-WSH", "c-ELW", "p-AU"]


def parse_au(path: Path) -> list[dict]:
    """IQ-TREE's USER TREES table, columns named from its own header row.

    Without ``--test-weight`` IQ-TREE omits p-WKH and p-WSH, so the layout is read, not assumed.
    """
    text = path.read_text()
    block = text[text.index("Tree      logL"):]
    lines = block.splitlines()
    names = lines[0].split()[1:]
    rows = []
    for line in lines[2:]:
        if not line.strip() or not re.match(r"\s*\d+\s", line):
            if rows:
                break
            continue
        fields = line.split()
        number, values = fields[0], [f for f in fields[1:] if f not in "+-"]
        rows.append({"tree": int(number), **dict(zip(names, values))})
    return rows


def satisfies(hyp: dict, scope: str, newick: str, outgroups: list[str]) -> bool:
    tips, splits = tree_splits(newick)
    bip = hyp[scope]
    free = set(hyp.get("free_in_constraint", []))
    if hyp["kind"] == "quartet":
        pair = set(bip["side_a"])
        single = set(bip["side_b_named"])
        outs = set(outgroups) & tips
        if not outs:
            raise SystemExit(f"none of the outgroups {outgroups} is a tip; a quartet hypothesis "
                             f"cannot be read without one")
        universe = frozenset(pair | single | outs)
        target = normalise(frozenset(pair), universe)
        return any(normalise(side & universe, universe) == target for side, _, _ in splits
                   if 2 <= len(side & universe) <= len(universe) - 2)
    return find_branch(set(bip["side_a"]), free, splits, tips) is not None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--hypotheses", type=Path, required=True)
    parser.add_argument("--scope", choices=["subset", "panel"], required=True)
    parser.add_argument("--outgroups", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    hypotheses = json.loads(args.hypotheses.read_text())
    candidates = [line.split("\t") for line in
                  (args.topology / "candidates.tsv").read_text().splitlines() if line.strip()]
    au = parse_au(args.topology / "au.iqtree")
    if len(au) != len(candidates):
        raise SystemExit(f"{len(au)} rows in au.iqtree, {len(candidates)} candidates")
    model = ""
    for line in (args.topology / "environment.txt").read_text().splitlines():
        if line.startswith("model:"):
            model = line.split(":", 1)[1].strip()
    rows = []
    for (label, path), result in zip(candidates, au):
        newick = Path(path).read_text().strip()
        carried = [n for n, h in hypotheses.items() if satisfies(h, args.scope, newick, args.outgroups)]
        rows.append({"candidate": label, "constrained_to": "" if label == "ml" else label,
                     "carries_hypotheses": ";".join(carried) or "none_of_the_named",
                     **{k: result.get(k, "") for k in COLUMNS},
                     "rejected_by_AU_at_0.05": float(result["p-AU"]) < 0.05,
                     "rejected_by_SH_at_0.05": float(result["p-SH"]) < 0.05,
                     "model": model, "scope": args.scope})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print("\t".join(str(row[k]) for k in ("candidate", "carries_hypotheses", "logL", "deltaL",
                                               "p-KH", "p-SH", "p-AU", "rejected_by_AU_at_0.05")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
