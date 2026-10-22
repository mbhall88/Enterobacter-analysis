#!/usr/bin/env python3
"""Regression test for coordinator review R1: G+C from exact integer base counts.

Both G+C producers (type_genome_gc.sh and circumscription_stats.measure) must
return the value computed from integer A/C/G/T counts summed over contigs. The
synthetic cases are chosen so that rebuilding counts from seqkit's rounded
per-contig percentages (-B) gives a different two-decimal G+C; the test checks
that the rounded route really differs, so it would catch a return to it. The
real-genome case is the oharae type (GCA_001729705.1), whose exact counts were
also derived independently by the coordinator (55.58 mol%, not 55.59).

Usage: test_gc_exact.py [--scripts DIR] [--fasta-paths FILE]
  --scripts      directory holding the two producers (default: this directory);
                 point it at superseded/2026-10-08-pre-review/scripts to show
                 the earlier code fails.
Exit status 0 only if every case passes.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent

CASES = {  # name: list of contig sequences
    # one contig, C=1 G=1 T=4: exact 33.33; rounded route 33.34
    "single_contig": ["CGTTTT"],
    # two contigs: exact 41.21; rounded route 41.20
    "two_contigs": ["A" * 34 + "C" * 15 + "G" * 29 + "T" * 32, "A" * 36 + "C" * 15 + "G" * 23 + "T" * 15],
    # soft-masked bases count; N and IUPAC codes are excluded from the denominator
    "masked_and_ambiguous": ["acgtNNNNACGGGCRYT", "ggccaattNNNN" + "G" * 7],
}
OHARAE = ("GCA_001729705.1", (1049958, 1302384, 1323445, 1048529), "55.58")


def exact(seqs: list[str]) -> tuple[tuple[int, int, int, int], str]:
    s = "".join(seqs).upper()
    a, c, g, t = (s.count(x) for x in "ACGT")
    return (a, c, g, t), f"{100 * (g + c) / (a + c + g + t):.2f}"


def rounded_route(fa: Path) -> str:
    out = subprocess.run(["seqkit", "fx2tab", "-n", "-l", "-B", "A", "-B", "C", "-B", "G", "-B", "T", str(fa)],
                         check=True, capture_output=True, text=True).stdout
    tot = [0.0] * 4
    for line in out.splitlines():
        f = line.split("\t")
        n = float(f[1])
        for i in range(4):
            tot[i] += float(f[2 + i]) * n / 100
    return f"{100 * (tot[1] + tot[2]) / sum(tot):.2f}"


def load_measure(scripts: Path):
    spec = importlib.util.spec_from_file_location("circ", scripts / "circumscription_stats.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.measure


def run_type_gc(scripts: Path, acc: str, fa: Path, tmp: Path) -> dict:
    out = tmp / f"{acc}.tsv"
    subprocess.run(["bash", str(scripts / "type_genome_gc.sh"), str(out), f"{acc}={fa}"], check=True,
                   capture_output=True, text=True)
    return next(csv.DictReader(open(out), delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scripts", type=Path, default=HERE)
    ap.add_argument("--fasta-paths", type=Path, default=HERE.parent / "inputs/type-genome-fasta-paths.txt")
    a = ap.parse_args()
    measure = load_measure(a.scripts)
    results = []

    def check(name: str, ok: bool, detail: str):
        results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}\t{name}\t{detail}")

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for name, seqs in CASES.items():
            fa = tmp / f"{name}.fa"
            fa.write_text("".join(f">{name}_{i}\n{s}\n" for i, s in enumerate(seqs, 1)))
            counts, gc = exact(seqs)
            rr = rounded_route(fa)
            if name != "masked_and_ambiguous":
                check(f"{name}:rounded_route_differs", rr != gc, f"exact={gc} rounded_route={rr}")
            row = run_type_gc(a.scripts, name, fa, tmp)
            got = tuple(int(float(row[k])) for k in "ACGT")
            check(f"{name}:type_genome_gc.sh", row["gc_mol_percent_acgt"] == gc and got == counts,
                  f"gc={row['gc_mol_percent_acgt']} expected={gc} counts={got} expected={counts}")
            m = measure((name, str(fa), ""))
            check(f"{name}:circumscription_stats.measure", f"{m[5]:.2f}" == gc,
                  f"gc={m[5]:.4f} expected={gc}")
        paths = dict(line.strip().split("=", 1) for line in open(a.fasta_paths) if "=" in line)
        acc, counts, gc = OHARAE
        fa = Path(paths[acc])
        row = run_type_gc(a.scripts, acc, fa, tmp)
        got = tuple(int(float(row[k])) for k in "ACGT")
        check(f"{acc}:type_genome_gc.sh", row["gc_mol_percent_acgt"] == gc and got == counts,
              f"gc={row['gc_mol_percent_acgt']} expected={gc} counts={got} expected={counts}")
        m = measure((acc, str(fa), ""))
        check(f"{acc}:circumscription_stats.measure", f"{m[5]:.2f}" == gc, f"gc={m[5]:.4f} expected={gc}")
        check(f"{acc}:rounded_route_differs", rounded_route(fa) != gc, f"rounded_route={rounded_route(fa)} exact={gc}")
    print(f"{sum(results)}/{len(results)} passed (scripts: {a.scripts})")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
