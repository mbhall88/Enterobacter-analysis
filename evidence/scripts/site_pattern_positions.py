#!/usr/bin/env python3
"""Line 6 of Task 17: do the genes carrying a minority history sit together on the chromosome?

For one quartet of named genomes (the fourth is the outgroup), every core gene's informative
columns are counted into the three pairings, exactly as ``summarise_concordance.site_patterns``
counts them.  Two comparisons are then made along a complete reference chromosome, using the gene
positions ``build_locus_ledger.py core --positional-reference`` wrote:

* **Runs.**  Genes whose columns favour pairing A over pairing B (and the reverse), in chromosome
  order, and the number of runs of same-state genes against 10,000 shuffles.  Fewer runs than
  the shuffles give is clustering: the signature of recombination or introgression tracts, as
  opposed to incomplete lineage sorting, which scatters discordant genes independently.
* **Neighbour correlation.**  The per-gene log-ratio of A to B columns, and its correlation with
  the next gene's (lag-1, Pearson), against the same shuffles.

Also reported: how much of the genome carries each history -- the share of informative columns,
and of deciding genes, on each side -- and a windowed profile (50 genes) for plotting.  Nothing
here is a support value.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from summarise_concordance import NT, read_alignment, read_loci  # noqa: E402


def per_gene_patterns(aln: dict[str, str], quartet: list[str], loci) -> list[tuple[str, list[int]]]:
    a, b, c, d = (aln[t] for t in quartet)
    out = []
    for locus, start, end in loci:
        counts = [0, 0, 0]
        for i in range(start - 1, end):
            x = (a[i].upper(), b[i].upper(), c[i].upper(), d[i].upper())
            if not all(ch in NT for ch in x) or len(set(x)) != 2:
                continue
            if x[0] == x[1] and x[2] == x[3]:
                counts[0] += 1
            elif x[0] == x[2] and x[1] == x[3]:
                counts[1] += 1
            elif x[0] == x[3] and x[1] == x[2]:
                counts[2] += 1
        out.append((locus, counts))
    return out


def count_runs(seq: list[str]) -> int:
    return 1 + sum(1 for x, y in zip(seq, seq[1:]) if x != y) if seq else 0


def lag1(values: list[float]) -> float:
    if len(values) < 3:
        return 0.0
    x, y = values[:-1], values[1:]
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sx = math.sqrt(sum((v - mx) ** 2 for v in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    return sum((p - mx) * (q - my) for p, q in zip(x, y)) / (sx * sy) if sx and sy else 0.0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True, help="concatenation-genes.nex")
    parser.add_argument("--positions", type=Path, required=True, help="gene-positions.tsv")
    parser.add_argument("--quartet", required=True, help="four comma-separated tips")
    parser.add_argument("--compare", nargs=2, type=int, default=[0, 1],
                        help="indexes (0,1,2) of the two pairings to compare: 0 = q1q2|q3q4, "
                             "1 = q1q3|q2q4, 2 = q1q4|q2q3")
    parser.add_argument("--label", required=True)
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument("--shuffles", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.output / args.label
    if stem.with_suffix(".json").exists():
        raise SystemExit(f"refusing to overwrite {stem}.json")
    quartet = args.quartet.split(",")
    aln = read_alignment(args.alignment)
    loci = read_loci(args.genes)
    with args.positions.open(newline="") as handle:
        pos = {r["locus"]: int(r["genome_position"]) for r in csv.DictReader(handle, delimiter="\t")
               if r["placed"] == "True"}
    patterns = per_gene_patterns(aln, quartet, loci)
    labels = [f"{quartet[0]},{quartet[1]}|{quartet[2]},{quartet[3]}",
              f"{quartet[0]},{quartet[2]}|{quartet[1]},{quartet[3]}",
              f"{quartet[0]},{quartet[3]}|{quartet[1]},{quartet[2]}"]
    ia, ib = args.compare
    placed = sorted((pos[l], l, c) for l, c in patterns if l in pos)
    rows, states, ratios = [], [], []
    for position, locus, counts in placed:
        a, b = counts[ia], counts[ib]
        state = "A" if a > b else "B" if b > a else "tie"
        ratio = math.log((a + 0.5) / (b + 0.5))
        rows.append({"genome_position": position, "locus": locus, "pattern_0": counts[0],
                     "pattern_1": counts[1], "pattern_2": counts[2], "A": a, "B": b,
                     "state": state, "log_ratio_A_B": round(ratio, 4)})
        ratios.append(ratio)
        if state != "tie":
            states.append(state)
    rng = random.Random(args.seed)
    observed_runs = count_runs(states)
    observed_r = lag1(ratios)
    runs_null, r_null = [], []
    work_s, work_r = list(states), list(ratios)
    for _ in range(args.shuffles):
        rng.shuffle(work_s)
        rng.shuffle(work_r)
        runs_null.append(count_runs(work_s))
        r_null.append(lag1(work_r))
    totals = [sum(c[k] for _, c in patterns) for k in range(3)]
    window = []
    for i in range(0, len(rows), args.window):
        block = rows[i:i + args.window]
        a = sum(r["A"] for r in block)
        b = sum(r["B"] for r in block)
        window.append({"first_position": block[0]["genome_position"],
                       "last_position": block[-1]["genome_position"], "genes": len(block),
                       "A": a, "B": b, "share_A_of_A_plus_B": round(a / (a + b), 4) if a + b else ""})
    with stem.with_suffix(".genes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with stem.with_suffix(".windows.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(window[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(window)
    result = {
        "quartet": quartet, "pairing_A": labels[ia], "pairing_B": labels[ib],
        "genes": len(patterns), "genes_placed": len(placed),
        "informative_columns": dict(zip(labels, totals)),
        "share_of_columns": {l: round(t / sum(totals), 4) for l, t in zip(labels, totals)},
        "genes_favouring_A": states.count("A"), "genes_favouring_B": states.count("B"),
        "genes_tied": sum(1 for r in rows if r["state"] == "tie"),
        "runs": {"observed": observed_runs, "null_median": statistics.median(runs_null),
                 "null_2.5%": sorted(runs_null)[int(0.025 * len(runs_null))],
                 "p_clustered": round((sum(1 for v in runs_null if v <= observed_runs) + 1)
                                      / (args.shuffles + 1), 5)},
        "lag1_correlation": {"observed": round(observed_r, 4),
                             "null_97.5%": round(sorted(r_null)[int(0.975 * len(r_null))], 4),
                             "p_positive": round((sum(1 for v in r_null if v >= observed_r) + 1)
                                                 / (args.shuffles + 1), 5)},
        "shuffles": args.shuffles, "seed": args.seed, "window_genes": args.window,
    }
    stem.with_suffix(".json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
