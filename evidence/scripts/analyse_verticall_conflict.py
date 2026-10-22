#!/usr/bin/env python3
"""Line 1 of Task 17: what Verticall's own record says about one named conflict.

Verticall's distance for a pair is built from the windows it paints *vertical*, so its horizontal
windows are already excluded -- and its tree still agrees with Mashtree against the marker tree.
"The distance methods are pulled by gene transfer" therefore cannot be assumed; it has to be shown
from the record.  Three things are computed, all from Task 11's 1,089-genome table (plus the
second *E. sichuanensis* assembly's own run), none from a tree:

1. **Lineage-pair summaries.**  For every pair of named lineages, the distribution over all genome
   pairs of the vertical-window distance FastME used (``median_vertical_window_distance``), the
   all-window ``mean_distance``, ``r/m``, and the horizontal and aligned fractions.

2. **Four-point tests.**  For each quartet of lineages and every combination of one genome per
   lineage, the pairing whose two within-pair distances sum smallest is the one the distances
   support locally -- the condition neighbour joining and minimum evolution act on.  Scored on the
   Jukes-Cantor-corrected vertical distance FastME actually read, on the uncorrected vertical
   distance, and on the all-window distance, so the effect of the correction and of the vertical
   filter are both measured rather than assumed.  The margin is reported as the gap between the best and second
   pairing, relative to the best sum.

3. **Where the horizontal regions are.**  For a complete reference genome R and each lineage Y,
   R's horizontal regions against every member of Y, placed on R's chromosome.  Reported: the
   fraction of R horizontal against Y; the part of R horizontal against at least half of Y (a
   lineage-level signal, not one genome's); whether those regions are clustered, by the index of
   dispersion over 50-kb bins against 1,000 random placements of the same intervals; and how much
   of the lineage-level horizontal sequence two lineages share with R in the same place.

No statistic here is a support value, and none is put on a scale with one.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import random
import statistics
import sys
from pathlib import Path

import numpy as np

csv.field_size_limit(sys.maxsize)

DIST_V = "median_vertical_window_distance"
DIST_ALL = "mean_distance"
# What Task 11's FastME actually read: ``verticall matrix`` applies a Jukes-Cantor correction to
# the chosen column unless --no_jukes_cantor is given (verticall/matrix.py:34), so the matrix cell
# for a pair is JC(median_vertical_window_distance), averaged over both directions afterwards.
DIST_V_JC = "jc:median_vertical_window_distance"


def jukes_cantor(p: float) -> float:
    return -0.75 * math.log(1 - 4.0 * p / 3.0)


def pct(text: str) -> float:
    return float(text.rstrip("%")) / 100 if text.endswith("%") else float(text)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Pairs:
    """Directed Verticall rows, with symmetric distances averaged over both directions."""

    def __init__(self):
        self.rows: dict[tuple[str, str], dict[str, str]] = {}

    def add(self, row: dict[str, str]) -> None:
        if row.get("result_level", "primary") != "primary":
            return
        key = (row["assembly_a"], row["assembly_b"])
        if key in self.rows:
            raise SystemExit(f"duplicate primary row {key}")
        self.rows[key] = row

    def value(self, a: str, b: str, column: str) -> float | None:
        corrected = column.startswith("jc:")
        column = column[3:] if corrected else column
        values = [pct(self.rows[k][column]) for k in ((a, b), (b, a))
                  if k in self.rows and self.rows[k][column] not in ("", "-")]
        if corrected:
            values = [jukes_cantor(v) for v in values]
        return sum(values) / len(values) if values else None


def load_pairs(paths: list[Path]) -> Pairs:
    pairs = Pairs()
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pairs.add(row)
    return pairs


def resolve_lineages(spec: dict, subset_dir: Path) -> dict[str, list[str]]:
    cores = json.loads((subset_dir / "lineage-cores.json").read_text())
    units = read_tsv(subset_dir / "units.tsv")
    out: dict[str, list[str]] = {}
    for name, entry in spec["lineages"].items():
        if "tips" in entry:
            out[name] = list(entry["tips"])
        elif "core_of" in entry:
            key = ",".join(sorted(entry["core_of"]))
            if key not in cores:
                raise SystemExit(f"no lineage core for {key} in {subset_dir}")
            out[name] = cores[key]["core"]
        elif "core_of_shared_unit" in entry:
            tip = entry["core_of_shared_unit"]
            owners = [u for u in units if tip in u["members"].split(";")]
            if len(owners) != 1:
                raise SystemExit(f"{tip} is in {len(owners)} units")
            out[name] = owners[0]["members"].split(";")
        else:
            raise SystemExit(f"lineage {name} has no recognised form")
    for tip in spec.get("each_tip_as_lineage", []):
        out.setdefault(tip, [tip])
    return out


def summarise(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    values = sorted(values)
    q = statistics.quantiles(values, n=4) if len(values) >= 2 else [values[0]] * 3
    return {"n": len(values), "min": values[0], "q1": q[0], "median": statistics.median(values),
            "q3": q[2], "max": values[-1]}


def lineage_pair_rows(lineages: dict[str, list[str]], pairs: Pairs) -> list[dict]:
    rows = []
    for x, y in itertools.combinations(sorted(lineages), 2):
        dv, da, rm, hz, af = [], [], [], [], []
        for a in lineages[x]:
            for b in lineages[y]:
                if a == b:
                    continue
                for key in ((a, b), (b, a)):
                    row = pairs.rows.get(key)
                    if row is None:
                        continue
                    dv.append(pct(row[DIST_V]))
                    da.append(pct(row[DIST_ALL]))
                    rm.append(pct(row["r/m"]))
                    hz.append(pct(row["alignments_horizontal_fraction"]))
                    af.append(pct(row["aligned_fraction"]))
        out = {"lineage_a": x, "lineage_b": y}
        for label, values in (("vertical_window_distance", dv), ("all_window_distance", da),
                              ("r_over_m", rm), ("horizontal_fraction_of_alignment", hz),
                              ("aligned_fraction", af)):
            s = summarise(values)
            out[f"{label}_n"] = s["n"]
            for k in ("median", "q1", "q3"):
                out[f"{label}_{k}"] = f"{s[k]:.6f}" if s["n"] else ""
        rows.append(out)
    return rows


def parse_topology(text: str) -> tuple[frozenset, frozenset]:
    left, right = text.split("|")
    return frozenset(left.split(",")), frozenset(right.split(","))


def four_point(quartet: dict, fourth: str, lineages: dict[str, list[str]], pairs: Pairs,
               column: str, rng: random.Random, cap: int = 200000) -> dict:
    names = [fourth if n == "X" else n for n in quartet["lineages"]]
    topologies = {}
    for label, text in quartet["topologies"].items():
        left, right = parse_topology(text)
        topologies[label] = (frozenset(fourth if n == "X" else n for n in left),
                             frozenset(fourth if n == "X" else n for n in right))
    members = [lineages[n] for n in names]
    combos = list(itertools.product(*members))
    if len(combos) > cap:
        combos = rng.sample(combos, cap)
    wins = {label: 0 for label in topologies}
    gaps: dict[str, list[float]] = {label: [] for label in topologies}
    missing = ties = 0
    for combo in combos:
        pick = dict(zip(names, combo))
        if len(set(combo)) < 4:
            continue
        sums = {}
        for label, (left, right) in topologies.items():
            l1, l2 = sorted(left)
            r1, r2 = sorted(right)
            d1 = pairs.value(pick[l1], pick[l2], column)
            d2 = pairs.value(pick[r1], pick[r2], column)
            if d1 is None or d2 is None:
                sums = None
                break
            sums[label] = d1 + d2
        if sums is None:
            missing += 1
            continue
        ordered = sorted(sums.items(), key=lambda kv: kv[1])
        if ordered[0][1] == ordered[1][1]:
            ties += 1
            continue
        wins[ordered[0][0]] += 1
        gaps[ordered[0][0]].append((ordered[1][1] - ordered[0][1]) / ordered[0][1])
    scored = sum(wins.values())
    return {"quartet": quartet["label"], "fourth": fourth, "distance_column": column,
            "combinations": len(combos), "scored": scored, "missing": missing, "ties": ties,
            **{f"fraction_{k}": (round(v / scored, 4) if scored else None) for k, v in wins.items()},
            **{f"median_relative_gap_{k}": (round(statistics.median(g), 5) if g else None)
               for k, g in gaps.items()}}


# --------------------------------------------------------------------------------------------
# Positional analysis
# --------------------------------------------------------------------------------------------

def contig_offsets(fasta: Path) -> tuple[dict[str, int], int]:
    """Offsets of each record when the records are laid end to end in file order."""
    offsets, total, name = {}, 0, None
    lengths: dict[str, int] = {}
    order = []
    with fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                order.append(name)
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())
    for name in order:
        offsets[name] = total
        total += lengths[name]
    return offsets, total


def parse_regions(text: str, offsets: dict[str, int]) -> list[tuple[int, int]]:
    out = []
    if not text:
        return out
    for item in text.split(","):
        contig, span = item.rsplit(":", 1)
        start, end = span.split("-")
        base = offsets[contig]
        out.append((base + int(start), base + int(end)))
    return out


def coverage_counts(intervals_per_member: list[list[tuple[int, int]]], length: int):
    """How many members call each position horizontal (a numpy difference array)."""
    diff = np.zeros(length + 1, dtype=np.int32)
    for intervals in intervals_per_member:
        for start, end in merge(intervals):
            diff[start] += 1
            diff[min(end, length)] -= 1
    return np.cumsum(diff[:-1])


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(s, e) for s, e in out]


def runs(mask) -> list[tuple[int, int]]:
    padded = np.concatenate(([False], np.asarray(mask, dtype=bool), [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return [(int(edges[i]), int(edges[i + 1])) for i in range(0, len(edges), 2)]


def intersect_bp(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    i = j = total = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def dispersion(intervals: list[tuple[int, int]], length: int, bin_size: int) -> float:
    bins = [0] * (length // bin_size + 1)
    for start, end in intervals:
        pos = start
        while pos < end:
            b = pos // bin_size
            stop = min(end, (b + 1) * bin_size)
            bins[b] += stop - pos
            pos = stop
    mean = statistics.fmean(bins)
    return statistics.pvariance(bins) / mean if mean else 0.0


def clustering_test(intervals: list[tuple[int, int]], length: int, rng: random.Random,
                    bin_size: int = 50000, permutations: int = 1000) -> dict:
    observed = dispersion(intervals, length, bin_size)
    null = []
    for _ in range(permutations):
        placed = []
        for start, end in intervals:
            size = end - start
            new = rng.randrange(0, max(length - size, 1))
            placed.append((new, new + size))
        null.append(dispersion(placed, length, bin_size))
    exceed = sum(1 for v in null if v >= observed)
    return {"index_of_dispersion": round(observed, 1),
            "null_median": round(statistics.median(null), 1) if null else None,
            "p_permutation": round((exceed + 1) / (permutations + 1), 4),
            "bin_size": bin_size, "permutations": permutations}


def positional(reference: str, fasta: Path, lineages: dict[str, list[str]], names: list[str],
               regions: dict[tuple[str, str], dict[str, str]], rng: random.Random) -> dict:
    offsets, length = contig_offsets(fasta)
    per_lineage = {}
    consistent_masks = {}
    for name in names:
        members = [m for m in lineages[name] if m != reference and (reference, m) in regions]
        if not members:
            continue
        horizontal = [parse_regions(regions[(reference, m)]["assembly_a_horizontal_regions"],
                                    offsets) for m in members]
        fractions = [sum(e - s for s, e in merge(h)) / length for h in horizontal]
        counts = coverage_counts(horizontal, length)
        need = (len(members) + 1) // 2
        consistent = runs(counts >= need)
        consistent_masks[name] = consistent
        consistent_bp = sum(e - s for s, e in consistent)
        entry = {"members_compared": len(members),
                 "horizontal_fraction_of_reference": summarise(fractions),
                 "consistent_horizontal_bp": consistent_bp,
                 "consistent_horizontal_fraction": round(consistent_bp / length, 5),
                 "consistent_regions": len(consistent)}
        if consistent:
            entry["clustering"] = clustering_test(consistent, length, rng)
            entry["largest_consistent_regions"] = [
                {"start": s, "end": e, "bp": e - s}
                for s, e in sorted(consistent, key=lambda r: r[0] - r[1])[:5]]
        per_lineage[name] = entry
    overlaps = []
    for x, y in itertools.combinations(sorted(consistent_masks), 2):
        a, b = consistent_masks[x], consistent_masks[y]
        both = intersect_bp(a, b)
        either = sum(e - s for s, e in a) + sum(e - s for s, e in b) - both
        overlaps.append({"lineage_a": x, "lineage_b": y, "shared_bp": both,
                         "jaccard": round(both / either, 4) if either else None})
    return {"reference": reference, "reference_length": length,
            "reference_records": len(offsets), "per_lineage": per_lineage,
            "consistent_region_overlap": overlaps}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--conflict", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, action="append", required=True,
                        help="Verticall tables without regions (repeatable)")
    parser.add_argument("--regions", type=Path, action="append", required=True,
                        help="Verticall tables with regions (repeatable)")
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text())
    spec = config["conflicts"][args.conflict]
    subset_dir = args.run_dir / spec["directory"] / "subset"
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    line1 = spec["line1"]
    lineages = resolve_lineages(line1, subset_dir)
    members = {r["tip_id"] for r in read_tsv(subset_dir / "subset.tsv")}
    for name, tips in lineages.items():
        stray = [t for t in tips if t not in members]
        if stray:
            raise SystemExit(f"lineage {name} has tips outside the subset: {stray[:3]}")

    pairs = load_pairs(args.summary)
    rng = random.Random(args.seed)
    lp = lineage_pair_rows(lineages, pairs)
    with (args.output / "lineage-pairs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(lp[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(lp)

    fp_rows = []
    for quartet in line1["quartets"]:
        for fourth in quartet["fourth"]:
            for column in (DIST_V_JC, DIST_V, DIST_ALL):
                fp_rows.append(four_point(quartet, fourth, lineages, pairs, column, rng))
    fields = sorted({k for r in fp_rows for k in r}, key=lambda k: (not k.startswith(("quartet", "fourth", "distance", "comb", "scored", "missing", "ties")), k))
    with (args.output / "four-point.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(fp_rows)

    regions: dict[tuple[str, str], dict[str, str]] = {}
    for path in args.regions:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("result_level", "primary") == "primary":
                    regions[(row["assembly_a"], row["assembly_b"])] = row
    subset_rows = {r["tip_id"]: r for r in read_tsv(subset_dir / "subset.tsv")}
    root = Path(".")
    pos = []
    for ref in line1["positional_references"]:
        fasta = root / subset_rows[ref]["derivative_fasta_path"]
        pos.append(positional(ref, fasta, lineages, line1["positional_lineages"], regions, rng))
    (args.output / "positional.json").write_text(json.dumps(pos, indent=1) + "\n")
    record = {"conflict": args.conflict, "lineages": {k: len(v) for k, v in lineages.items()},
              "summary_tables": [str(p) for p in args.summary],
              "region_tables": [str(p) for p in args.regions], "seed": args.seed,
              "directed_rows_loaded": len(pairs.rows), "region_rows_loaded": len(regions)}
    (args.output / "record.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
