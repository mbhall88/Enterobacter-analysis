#!/usr/bin/env python3
"""Is the disagreement the data or the method?  Rebuild distance trees on one focused subset.

Both whole-genome methods agree against the marker tree, and Verticall's distance already
excludes horizontal windows (line 1).  A remaining explanation is the method itself: neighbour
joining and minimum evolution read a distance matrix, and a matrix of distances that saturate --
because only the alignable part of two genomes is compared, or because rate variation among sites
is ignored -- is not additive, which pairs long branches.  That explanation makes a prediction
that could fail: the distance topology should appear from uncorrected distances and disappear as
the correction approaches the model the ML tree uses, *on the same data*.

FastME is run exactly as Tasks 08 and 11 ran it (``-m B -n B -s``: TaxAdd_BalME start, balanced
NNI, SPR), on:

* Verticall's table restricted to the subset -- uncorrected vertical-window distance; the
  Jukes-Cantor-corrected one ``verticall matrix`` writes by default (what Task 11's FastME read);
  and gamma-corrected Jukes-Cantor at each ``--gamma-alpha``;
* an alignment (the core-gene concatenation) -- uncorrected p-distance with pairwise deletion,
  JC, and gamma-JC at each alpha;
* any IQ-TREE ``.mldist`` given, which is the ML distance under that run's fitted model.

Each tree is read for the hypotheses it carries, by tip set.  Nothing here is a support value.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import format_newick, parse_newick, prune_to_tips  # noqa: E402
from summarise_topology_tests import satisfies  # noqa: E402

csv.field_size_limit(sys.maxsize)


def jc(p: float) -> float:
    x = 1 - 4.0 * p / 3.0
    return float("inf") if x <= 0 else -0.75 * math.log(x)


def gamma_jc(p: float, alpha: float) -> float:
    x = 1 - 4.0 * p / 3.0
    return float("inf") if x <= 0 else 0.75 * alpha * (x ** (-1.0 / alpha) - 1)


def pct(text: str) -> float:
    return float(text.rstrip("%")) / 100 if text.endswith("%") else float(text)


def verticall_matrices(tables: list[Path], ids: list[str], alphas: list[float]) -> dict:
    wanted = set(ids)
    directed: dict[tuple[str, str], float] = {}
    for path in tables:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("result_level", "primary") != "primary":
                    continue
                a, b = row["assembly_a"], row["assembly_b"]
                if a in wanted and b in wanted:
                    directed[(a, b)] = pct(row["median_vertical_window_distance"])
    transforms = {"verticall-p": lambda p: p, "verticall-jc": jc}
    for alpha in alphas:
        transforms[f"verticall-gammajc-a{alpha:g}"] = lambda p, a=alpha: gamma_jc(p, a)
    out = {}
    for label, f in transforms.items():
        n = len(ids)
        m = np.zeros((n, n))
        for i, a in enumerate(ids):
            for j, b in enumerate(ids):
                if i < j:
                    vals = [f(directed[k]) for k in ((a, b), (b, a)) if k in directed]
                    if not vals:
                        raise SystemExit(f"no Verticall distance for {a} {b}")
                    m[i, j] = m[j, i] = sum(vals) / len(vals)
        out[label] = m
    return out


def alignment_matrices(path: Path, ids: list[str], alphas: list[float]) -> dict:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip().upper())
    missing = [t for t in ids if t not in seqs]
    if missing:
        raise SystemExit(f"alignment lacks {missing[:3]}")
    codes = {c: i for i, c in enumerate("ACGT")}
    lut = np.full(256, 255, dtype=np.uint8)
    for c, i in codes.items():
        lut[ord(c)] = i
    arr = np.stack([lut[np.frombuffer("".join(seqs[t]).encode(), dtype=np.uint8)] for t in ids])
    present = arr != 255
    # Variable columns carry every mismatch; constant columns only add to the denominators.
    variable = np.zeros(arr.shape[1], dtype=bool)
    for state in range(4):
        has = ((arr == state) & present).any(axis=0)
        variable |= has & ((arr != state) & present).any(axis=0)
    var = arr[:, variable]
    var_present = present[:, variable]
    n = len(ids)
    shared = np.zeros((n, n))
    pres = present.astype(np.float32)
    for start in range(0, arr.shape[1], 500000):
        block = pres[:, start:start + 500000]
        shared += block @ block.T
    mism = np.zeros((n, n))
    for i in range(n):
        both = var_present[i] & var_present
        mism[i] = (both & (var[i] != var)).sum(axis=1)
    p = np.divide(mism, shared, out=np.zeros_like(mism), where=shared > 0)
    out = {"alignment-p": p,
           "alignment-jc": np.vectorize(jc)(p)}
    for alpha in alphas:
        out[f"alignment-gammajc-a{alpha:g}"] = np.vectorize(lambda x, a=alpha: gamma_jc(x, a))(p)
    info = {"columns": int(arr.shape[1]), "variable_columns": int(variable.sum()),
            "min_shared_columns": int(shared[np.triu_indices(n, 1)].min())}
    return out, info


def read_mldist(path: Path, ids: list[str]) -> np.ndarray:
    lines = path.read_text().split("\n")
    n = int(lines[0].split()[0])
    rows = {}
    for line in lines[1:n + 1]:
        parts = line.split()
        rows[parts[0]] = [float(x) for x in parts[1:]]
    order = [l.split()[0] for l in lines[1:n + 1]]
    index = {t: i for i, t in enumerate(order)}
    m = np.zeros((len(ids), len(ids)))
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            m[i, j] = rows[a][index[b]]
    return m


def read_phylip(path: Path, ids: list[str]) -> np.ndarray:
    lines = [l for l in path.read_text().split("\n") if l.strip()]
    n = int(lines[0].split()[0])
    names = [l.split()[0] for l in lines[1:n + 1]]
    index = {t: i for i, t in enumerate(names)}
    missing = [t for t in ids if t not in index]
    if missing:
        raise SystemExit(f"{path} lacks {missing[:3]}")
    full = [lines[1 + index[t]].split()[1:] for t in ids]
    return np.array([[float(row[index[b]]) for b in ids] for row in full])


def restrict(newick: str, keep: set[str]) -> str:
    tree = parse_newick(newick)
    nodes = {tree.label[n] for n in tree.tips()}
    missing = keep - nodes
    if missing:
        raise SystemExit(f"tree lacks {sorted(missing)[:3]}")
    result = prune_to_tips(tree, {n for n in tree.tips() if tree.label[n] in keep})
    return format_newick(result.tree, decimals=12)


def write_phylip(path: Path, ids: list[str], m: np.ndarray) -> None:
    if not np.isfinite(m).all():
        raise SystemExit(f"{path.name}: saturated distance (infinite after correction)")
    # JC(0) is -0.0, which prints as "-0.000000000" and FastME refuses it as a non-number.
    m = np.where(m == 0, 0.0, m)
    with path.open("w") as handle:
        handle.write(f"{len(ids)}\n")
        for i, t in enumerate(ids):
            handle.write(t + "\t" + "\t".join(f"{v:.9f}" for v in m[i]) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ids", type=Path, required=True)
    parser.add_argument("--verticall", type=Path, action="append", default=[])
    parser.add_argument("--alignment", type=Path)
    parser.add_argument("--mldist", action="append", default=[], help="label=path")
    parser.add_argument("--phylip", action="append", default=[],
                        help="label=path of an existing PHYLIP matrix, restricted to --ids")
    parser.add_argument("--evaluate-on", type=Path,
                        help="ids to prune each tree to before reading hypotheses (default --ids)")
    parser.add_argument("--gamma-alpha", type=float, action="append", default=[])
    parser.add_argument("--fastme", type=Path, required=True)
    parser.add_argument("--hypotheses", type=Path, required=True)
    parser.add_argument("--scope", choices=["subset", "panel"], required=True)
    parser.add_argument("--outgroups", nargs="+", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    (args.output / "matrices").mkdir(parents=True, exist_ok=True)
    (args.output / "trees").mkdir(exist_ok=True)
    ids = [l.split()[0] for l in args.ids.read_text().splitlines() if l.strip()]
    matrices: dict[str, np.ndarray] = {}
    record: dict = {"ids": str(args.ids), "tips": len(ids), "gamma_alphas": args.gamma_alpha}
    if args.verticall:
        matrices.update(verticall_matrices(args.verticall, ids, args.gamma_alpha))
    if args.alignment:
        mats, info = alignment_matrices(args.alignment, ids, args.gamma_alpha)
        matrices.update(mats)
        record["alignment"] = {"path": str(args.alignment), **info}
    for spec in args.mldist:
        label, path = spec.split("=", 1)
        matrices[f"mldist-{label}"] = read_mldist(Path(path), ids)
    for spec in args.phylip:
        label, path = spec.split("=", 1)
        matrices[label] = read_phylip(Path(path), ids)
    evaluate = set(ids)
    if args.evaluate_on:
        evaluate = {l.split()[0] for l in args.evaluate_on.read_text().splitlines() if l.strip()}
        record["evaluated_on"] = {"ids": str(args.evaluate_on), "tips": len(evaluate)}
    hypotheses = json.loads(args.hypotheses.read_text())
    version = subprocess.run([str(args.fastme), "--version"], capture_output=True, text=True)
    record["fastme"] = (version.stdout or version.stderr).strip().splitlines()[0:1]
    rows = []
    for label, m in matrices.items():
        phylip = args.output / "matrices" / f"{label}.phylip"
        tree = args.output / "trees" / f"{label}.nwk"
        write_phylip(phylip, ids, m)
        result = subprocess.run([str(args.fastme), "-i", str(phylip), "-o", str(tree),
                                 "-I", str(args.output / "trees" / f"{label}.info"),
                                 "-m", "B", "-n", "B", "-s", "-T", str(args.threads)],
                                capture_output=True, text=True)
        if result.returncode:
            raise SystemExit(f"fastme failed on {label}: {result.stderr}")
        newick = tree.read_text().strip()
        if evaluate != set(ids):
            newick = restrict(newick, evaluate)
        carried = [n for n, h in hypotheses.items() if satisfies(h, args.scope, newick, args.outgroups)]
        upper = m[np.triu_indices(len(ids), 1)]
        rows.append({"matrix": label, "carries_hypotheses": ";".join(carried) or "none_of_the_named",
                     "max_distance": round(float(upper.max()), 6),
                     "median_distance": round(float(np.median(upper)), 6)})
    with (args.output / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    record["results"] = rows
    (args.output / "record.json").write_text(json.dumps(record, indent=1) + "\n")
    for row in rows:
        print(f"{row['matrix']:28s} {row['carries_hypotheses']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
