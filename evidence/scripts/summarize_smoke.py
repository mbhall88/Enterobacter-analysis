#!/usr/bin/env python3
"""Validate eight-genome engineering outputs and report their limitations."""
import collections
import csv
import gzip
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MARKER = ROOT / "results/marker-smoke/2026-09-11"
VERTICALL = ROOT / "results/verticall-smoke/2026-09-11"


def fasta(path):
    seqs = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line.strip()[1:]
                assert name not in seqs
                seqs[name] = ""
            else:
                seqs[name] += line.strip()
    assert len({len(s) for s in seqs.values()}) == 1
    return seqs


def alignment_summary(seqs):
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    variable = informative = empty = 0
    for column in zip(*seqs.values()):
        counts = collections.Counter(c.upper() for c in column if c.upper() in aa)
        variable += len(counts) > 1
        informative += sum(n >= 2 for n in counts.values()) >= 2
        empty += not counts
    return {"tips": len(seqs), "columns": len(next(iter(seqs.values()))),
            "variable_columns": variable, "parsimony_informative_columns": informative,
            "columns_without_standard_amino_acids": empty,
            "standard_amino_acid_percent": {k: 100 * sum(c.upper() in aa for c in s) / len(s)
                                           for k, s in seqs.items()}}


def main():
    summary = {"purpose": "Engineering test only; not a frozen common-input scientific comparison."}
    canonical = fasta(MARKER / "align-canonical/align/gtdbtk.bac120.user_msa.fasta.gz")
    full = fasta(MARKER / "align-full/align/gtdbtk.bac120.user_msa.fasta.gz")
    summary["canonical"] = alignment_summary(canonical)
    summary["full"] = alignment_summary(full)
    assert set(canonical) == set(full)
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    # Prespecified simple engineering sensitivity: at least 80% standard-AA occupancy.
    # This is not a full assessment of alignment quality or marker orthology.
    keep = [i for i, col in enumerate(zip(*full.values()))
            if sum(c.upper() in aa for c in col) / len(col) >= 0.8]
    filtered = {k: "".join(s[i] for i in keep) for k, s in full.items()}
    (MARKER / "full-occupancy80.faa").write_text("".join(f">{k}\n{s}\n" for k, s in filtered.items()))
    summary["full_occupancy80"] = alignment_summary(filtered)
    with (VERTICALL / "pairwise.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pairs = collections.Counter((r["assembly_a"], r["assembly_b"]) for r in rows)
    expected = {(a, b) for a in canonical for b in canonical if a != b}
    assert set(pairs) == expected, "Missing or unexpected directed pairs"
    fields = [k for k in rows[0] if not k.endswith("_regions")]
    with (VERTICALL / "pairwise-qc.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    lines = (VERTICALL / "distances.phylip").read_text().splitlines()
    n = int(lines[0])
    names = [line.split()[0] for line in lines[1:] if line.strip()]
    matrix = [[float(x) for x in line.split()[1:]] for line in lines[1:] if line.strip()]
    assert n == len(names) == len(set(names)) == len(canonical)
    assert set(names) == set(canonical)
    assert all(len(row) == n for row in matrix)
    assert all(math.isfinite(x) and x >= 0 for row in matrix for x in row)
    assert all(matrix[i][i] == 0 for i in range(n))
    asymmetry = max(abs(matrix[i][j] - matrix[j][i]) for i in range(n) for j in range(n))
    assert asymmetry < 1e-10
    summary["verticall"] = {"directed_pairs": len(pairs), "rows": len(rows),
        "pairs_with_multiple_results": sum(v > 1 for v in pairs.values()),
        "rows_with_multiple_detected_mass_peaks": sum("," in r["mass_peaks"] for r in rows),
        "result_levels": dict(collections.Counter(r["result_level"] for r in rows)),
        "aligned_fraction_range": [min(float(r["aligned_fraction"]) for r in rows),
                                   max(float(r["aligned_fraction"]) for r in rows)],
        "assembly_a_vertical_percent_range": [min(float(r["assembly_a_vertical_fraction"].rstrip("%")) for r in rows),
                                              max(float(r["assembly_a_vertical_fraction"].rstrip("%")) for r in rows)],
        "matrix_passed_finite_nonnegative_zero_diagonal_symmetry_and_tip_checks": True,
        "matrix_max_asymmetry": asymmetry}
    (MARKER / "smoke-qc.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
