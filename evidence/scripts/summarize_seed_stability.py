#!/usr/bin/env python3
"""Measure how much of a topology survives re-seeding, and annotate a tree with the answer.

**This does not produce branch support, and its output must never be reported as support.**
Mashtree's sketching is seeded: `mash` draws a 10,000-hash sketch of each genome and the draw
depends on the seed. Re-running with a different seed re-draws every sketch, so the distances and
therefore the tree move. The frequency with which a split reappears across independent seeds is a
**sketch-stability** statistic. It says how reproducible that split is under the method's own
randomness. It says nothing about how much evidence supports it, it is not a bootstrap, it is not
comparable with an SH-like local support or an ultrafast bootstrap, and a tree annotated by this
tool must never be tabulated on the same scale as one carrying either.

Task 09 established that the seed matters -- one alternate seed moved 22.3% of splits, while a
same-seed replicate was byte-identical -- and recorded that one draw bounds the magnitude without
estimating a distribution. This estimates the distribution, which is what makes a cross-method
difference interpretable: a difference smaller than a method's own seed-induced spread is not
evidence about method.

Two things are reported and kept separate:

* **Per-split reproducibility** -- the fraction of replicate trees containing each split of the
  reference tree, written both as a table and as internal labels on a copy of the reference tree.
* **Replicate-to-replicate distance** -- the distribution of Robinson-Foulds distances among the
  replicates themselves, which is the baseline any between-method distance has to clear.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_phylogenies import (  # noqa: E402
    ComparisonError, SupportConvention, canonical, extract_splits, tip_names)
from gtdb_context import format_newick, parse_newick  # noqa: E402

NO_SUPPORT = SupportConvention(
    kind="none", statistic="none; splits are read topologically for this measurement")


def load(path: Path):
    tree = parse_newick(path.read_text())
    return tree, tip_names(tree)


def splits_of(tree, index: dict[str, int]) -> set[int]:
    # near_zero below zero retains every branch: this measurement is purely topological, and
    # collapsing here would confound branch length with reproducibility.
    return set(extract_splits(tree, index, NO_SUPPORT, -1.0))


def summarize(reference: Path, replicates: list[Path], thresholds: list[float]) -> dict:
    ref_tree, ref_tips = load(reference)
    order = sorted(ref_tips)
    index = {tip: i for i, tip in enumerate(order)}
    n = len(order)

    replicate_splits = []
    for path in replicates:
        tree, tips = load(path)
        if sorted(tips) != order:
            raise ComparisonError(
                f"{path} does not carry the reference tip set: "
                f"{len(set(order) - set(tips))} missing, {len(set(tips) - set(order))} extra")
        replicate_splits.append(splits_of(tree, index))
    if len(replicate_splits) < 2:
        raise ComparisonError("at least two replicate trees are needed to estimate a spread")

    ref_splits = splits_of(ref_tree, index)
    frequency = {
        key: sum(1 for s in replicate_splits if key in s) / len(replicate_splits)
        for key in ref_splits
    }

    union: set[int] = set()
    for s in replicate_splits:
        union |= s
    replicate_frequency = {
        key: sum(1 for s in replicate_splits if key in s) / len(replicate_splits)
        for key in union
    }

    distances = []
    for i in range(len(replicate_splits)):
        for j in range(i + 1, len(replicate_splits)):
            a, b = replicate_splits[i], replicate_splits[j]
            raw = len(a - b) + len(b - a)
            distances.append({
                "raw_robinson_foulds": raw,
                "normalised_robinson_foulds": raw / (len(a) + len(b)) if (a or b) else 0.0,
                "percent_of_a_splits_absent_from_b": 100 * len(a - b) / len(a) if a else None,
            })
    percents = [d["percent_of_a_splits_absent_from_b"] for d in distances
                if d["percent_of_a_splits_absent_from_b"] is not None]
    raws = [d["raw_robinson_foulds"] for d in distances]

    def spread(values):
        ordered = sorted(values)
        return {"n": len(ordered), "min": ordered[0],
                "median": round(statistics.median(ordered), 4),
                "mean": round(statistics.fmean(ordered), 4), "max": ordered[-1]}

    values = sorted(frequency.values())
    return {
        "statistic": "seed_stability",
        "what_this_is_not": (
            "This is not branch support. It is the fraction of independent sketch seeds whose "
            "tree contains the split. It must never be compared with, merged into or tabulated "
            "on the same scale as an SH-like local support or an ultrafast bootstrap."),
        "reference_tree": str(reference),
        "replicate_trees": [str(p) for p in replicates],
        "replicates": len(replicate_splits),
        "tips": n,
        "reference_non_trivial_splits": len(ref_splits),
        "splits_appearing_in_at_least_one_replicate": len(union),
        "splits_in_a_replicate_but_not_the_reference": len(union - ref_splits),
        "reference_split_reproducibility": {
            "min": values[0] if values else None,
            "median": round(statistics.median(values), 4) if values else None,
            "mean": round(statistics.fmean(values), 4) if values else None,
            "max": values[-1] if values else None,
            "at_or_above": {
                f"{t:g}": sum(1 for v in values if v >= t) for t in thresholds},
            "percent_at_or_above": {
                f"{t:g}": round(100 * sum(1 for v in values if v >= t) / len(values), 1)
                for t in thresholds} if values else {},
        },
        "replicate_to_replicate_distance": {
            "pairs_compared": len(distances),
            "raw_robinson_foulds": spread(raws) if raws else None,
            "percent_of_splits_moved": spread(percents) if percents else None,
            "note": ("The baseline any between-method topological distance has to clear before "
                     "it can be attributed to method rather than to this method's own seed."),
        },
        "_frequency": frequency,
        "_replicate_frequency": replicate_frequency,
        "_order": order,
        "_index": index,
    }


def annotate(reference: Path, frequency: dict[int, float], index: dict[str, int],
             decimals: int) -> str:
    """A copy of the reference tree with each internal label replaced by its seed frequency."""
    tree = parse_newick(reference.read_text())
    n = len(index)
    full = (1 << n) - 1
    below: dict[int, int] = {}
    for node in tree.postorder():
        if not tree.children[node]:
            bit = index.get(tree.label[node])
            below[node] = 0 if bit is None else 1 << bit
        else:
            mask = 0
            for child in tree.children[node]:
                mask |= below[child]
            below[node] = mask
    for node in range(tree.size):
        if node == 0 or not tree.children[node]:
            continue
        key = canonical(below[node], full)
        value = frequency.get(key)
        tree.label[node] = "" if value is None else f"{value:.{decimals}f}"
    return format_newick(tree) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, type=Path,
                        help="the tree whose splits are scored; normally the primary run")
    parser.add_argument("--replicate", action="append", required=True, type=Path,
                        help="a replicate tree differing from the reference only in its seed; "
                             "repeat for each replicate")
    parser.add_argument("--threshold", action="append", type=float, default=None,
                        help="report the count of splits at or above this frequency; repeatable")
    parser.add_argument("--decimals", type=int, default=3)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args(argv)

    thresholds = args.threshold or [0.5, 0.7, 0.95, 1.0]
    try:
        result = summarize(args.reference, args.replicate, thresholds)
    except ComparisonError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    frequency = result.pop("_frequency")
    replicate_frequency = result.pop("_replicate_frequency")
    order = result.pop("_order")
    index = result.pop("_index")

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{args.output_prefix}.annotated.nwk").write_text(
        annotate(args.reference, frequency, index, args.decimals))
    rows = []
    n = len(order)
    for key, value in sorted(frequency.items(), key=lambda kv: (kv[1], -bin(kv[0]).count("1"))):
        size = bin(key).count("1")
        side = [order[i] for i in range(n) if key >> i & 1]
        if len(side) > n - len(side):
            side = [order[i] for i in range(n) if not (key >> i & 1)]
        rows.append({"smaller_side_tips": min(size, n - size),
                     "seed_reproducibility": round(value, 6),
                     "in_reference": True,
                     "example_tips": ";".join(sorted(side)[:4])})
    for key, value in sorted(replicate_frequency.items(), key=lambda kv: -kv[1]):
        if key in frequency:
            continue
        size = bin(key).count("1")
        side = [order[i] for i in range(n) if key >> i & 1]
        if len(side) > n - len(side):
            side = [order[i] for i in range(n) if not (key >> i & 1)]
        rows.append({"smaller_side_tips": min(size, n - size),
                     "seed_reproducibility": round(value, 6),
                     "in_reference": False,
                     "example_tips": ";".join(sorted(side)[:4])})
    import csv
    with Path(f"{args.output_prefix}.splits.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    Path(f"{args.output_prefix}.stability.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
