#!/usr/bin/env python3
"""Attach a within-method baseline to every cross-method conflict figure, or refuse to print it.

SET-18 rule 3: *quote every conflict figure with its within-method baseline in the same statistic.*
A conflict percentage on its own is not interpretable — a difference smaller than a method's own
spread under one changed parameter is not evidence about method at all. This turns that rule from
something a writer has to remember into something the table cannot violate: a cross-method figure
with no baseline for one of its two sides is **omitted and listed as unreportable**, not printed
bare.

How the baseline for a side is chosen
-------------------------------------
Baseline comparison sets are the ones whose name begins with ``baseline-``. Each compares one tree
against another form of *the same method* with one thing changed: the program, the model, the
distance mode, the sketch seed. For a given tree, its baseline is the **largest** such figure
available — the most conservative reading, because a cross-method difference has to clear the
method's own worst within-method spread before it means anything.

The statistic is the same one on both sides: ``percent_a_retained_conflicting_in_b``, the share of a
tree's retained splits that are incompatible with a split the other tree also retains. Percentages
computed under different collapsing rules are never mixed: rows are matched on
``collapsing_rule`` as well as on the trees.

Nothing here merges two support conventions, ranks methods, or converts one method's units into
another's. It only refuses to print a number without the number it has to be read against.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_pair(pair: str) -> tuple[str, str] | None:
    if "_vs_" not in pair:
        return None
    left, _, right = pair.partition("_vs_")
    return left, right


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary", type=Path, required=True,
                        help="summary.tsv from scripts/compare_phylogenies.py")
    parser.add_argument("--config", type=Path,
                        help="the comparison configuration, for each tree's method description")
    parser.add_argument("--baseline-prefix", default="baseline-",
                        help="comparison sets whose name starts with this are baselines")
    parser.add_argument("--set", action="append", default=[],
                        help="restrict the report to these comparison sets; repeatable")
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    rows = read_tsv(args.summary)
    methods: dict[str, str] = {}
    classes: dict[str, str] = {}
    if args.config:
        config = json.loads(args.config.read_text())
        for name, spec in config.get("trees", {}).items():
            methods[name] = spec.get("method", "")
            classes[name] = spec.get("inference_class", "")

    # tree -> collapsing rule -> list of (baseline set, other tree, percent)
    baselines: dict[tuple[str, str], list[tuple[str, str, float]]] = defaultdict(list)
    for row in rows:
        if not row["comparison_set"].startswith(args.baseline_prefix):
            continue
        pair = split_pair(row["pair"])
        if pair is None:
            continue
        a, b = pair
        rule = row["collapsing_rule"]
        for tree, other, column in ((a, b, "percent_a_retained_conflicting_in_b"),
                                    (b, a, "percent_b_retained_conflicting_in_a")):
            value = number(row.get(column, ""))
            if value is not None:
                baselines[(tree, rule)].append((row["comparison_set"], other, value))

    reportable, unreportable = [], []
    for row in rows:
        if row["comparison_set"].startswith(args.baseline_prefix):
            continue
        if args.set and row["comparison_set"] not in args.set:
            continue
        pair = split_pair(row["pair"])
        if pair is None:
            continue
        a, b = pair
        rule = row["collapsing_rule"]
        for tree, other, column, retained in (
                (a, b, "percent_a_retained_conflicting_in_b", "retained_splits_a"),
                (b, a, "percent_b_retained_conflicting_in_a", "retained_splits_b")):
            value = number(row.get(column, ""))
            entry = {
                "comparison_set": row["comparison_set"],
                "collapsing_rule": rule,
                "tips": row.get("tips", ""),
                "tree": tree,
                "against": other,
                "tree_method": methods.get(tree, ""),
                "tree_inference_class": classes.get(tree, ""),
                "against_inference_class": classes.get(other, ""),
                "retained_splits": row.get(retained, ""),
                "percent_conflicting": "" if value is None else f"{value:g}",
            }
            available = baselines.get((tree, rule), [])
            if value is None or not available:
                entry["reason_unreportable"] = (
                    "no conflict figure" if value is None
                    else f"no within-method baseline for {tree} under collapsing rule {rule!r}")
                unreportable.append(entry)
                continue
            worst = max(available, key=lambda item: item[2])
            entry.update({
                "baseline_percent": f"{worst[2]:g}",
                "baseline_set": worst[0],
                "baseline_against": worst[1],
                "baselines_available": len(available),
                "all_baselines": "; ".join(f"{name}={value:g}" for name, _, value in
                                           sorted(available, key=lambda item: -item[2])),
                "exceeds_own_baseline": str(value > worst[2]),
                "ratio_to_baseline": f"{value / worst[2]:.2f}" if worst[2] else "",
                "same_inference_class": str(classes.get(tree, "") == classes.get(other, "")
                                            and bool(classes.get(tree, ""))),
            })
            reportable.append(entry)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fields = ["comparison_set", "collapsing_rule", "tips", "tree", "against", "tree_method",
              "tree_inference_class", "against_inference_class", "same_inference_class",
              "retained_splits", "percent_conflicting", "baseline_percent", "baseline_set",
              "baseline_against", "baselines_available", "all_baselines", "exceeds_own_baseline",
              "ratio_to_baseline"]
    with prefix.with_name(prefix.name + ".tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for entry in reportable:
            writer.writerow({field: entry.get(field, "") for field in fields})
    with prefix.with_name(prefix.name + "-unreportable.tsv").open(
            "w", newline="", encoding="utf-8") as handle:
        fields2 = ["comparison_set", "collapsing_rule", "tree", "against", "percent_conflicting",
                   "reason_unreportable"]
        writer = csv.DictWriter(handle, fieldnames=fields2, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for entry in unreportable:
            writer.writerow({field: entry.get(field, "") for field in fields2})

    summary = {
        "rule": ("SET-18 rule 3: no cross-method conflict figure is printed without the "
                 "within-method baseline it must be read against, in the same statistic and under "
                 "the same collapsing rule."),
        "statistic": "percent of a tree's retained splits incompatible with a split the other "
                     "tree also retains",
        "reportable_figures": len(reportable),
        "unreportable_figures": len(unreportable),
        "trees_with_a_baseline": sorted({tree for tree, _ in baselines}),
        "figures_exceeding_their_own_baseline": sum(
            1 for entry in reportable if entry["exceeds_own_baseline"] == "True"),
        "same_inference_class_figures": sum(
            1 for entry in reportable if entry["same_inference_class"] == "True"),
        "outputs": {"reportable": str(prefix.with_name(prefix.name + ".tsv")),
                    "unreportable": str(prefix.with_name(prefix.name + "-unreportable.tsv"))},
    }
    prefix.with_name(prefix.name + ".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
