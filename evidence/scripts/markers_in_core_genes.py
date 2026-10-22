#!/usr/bin/env python3
"""Where do the BAC120 markers sit among the core genes' competing histories?

If the marker tree and genome-scale data disagree because the markers happen to sit in the part of
the genome that carries a minority history, the markers' core genes should be enriched for that
history.  This joins, on one complete reference chromosome, GTDB-Tk's marker-gene calls (its own
prodigal coordinates, from Task 11's ``identify`` run) to the core genes placed by
``build_locus_ledger.py core`` and classified by ``site_pattern_positions.py``, by coordinate
overlap, and compares the markers' states with every other deciding core gene's.

The comparison is Fisher's exact test (two-sided, computed longhand), and the markers are also
read directly: the summed site-pattern columns of their core genes.  Only the markers Task 11's
occupancy filter retained are counted.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


def gtdbtk_markers(marker_dir: Path, retained: set[str]) -> list[dict]:
    coords = {}
    for line in next(marker_dir.glob("*_protein.faa")).read_text().splitlines():
        if line.startswith(">"):
            parts = [p.strip() for p in line[1:].split("#")]
            gid = parts[0]
            contig = gid.rsplit("_", 1)[0]
            coords[gid] = (contig, int(parts[1]), int(parts[2]))
    hits = []
    for name in ("pfam", "tigrfam"):
        path = next(marker_dir.glob(f"*_{name}_tophit.tsv"))
        for row in csv.DictReader(path.open(), delimiter="\t"):
            for hit in row["Top hits (Family id,e-value,bitscore)"].split(";"):
                family = hit.split(",")[0]
                key = family.replace(".", "_")
                if key in retained:
                    contig, start, end = coords[row["Gene Id"]]
                    hits.append({"marker": key, "gene": row["Gene Id"], "contig": contig,
                                 "start": start, "end": end})
    return hits


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p for the 2x2 table [[a, b], [c, d]]."""
    n = a + b + c + d
    row1, col1 = a + b, a + c

    def prob(x):
        return math.exp(math.lgamma(row1 + 1) + math.lgamma(n - row1 + 1) + math.lgamma(col1 + 1)
                        + math.lgamma(n - col1 + 1) - math.lgamma(n + 1) - math.lgamma(x + 1)
                        - math.lgamma(row1 - x + 1) - math.lgamma(col1 - x + 1)
                        - math.lgamma(n - row1 - col1 + x + 1))
    observed = prob(a)
    lo, hi = max(0, row1 + col1 - n), min(row1, col1)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= observed * (1 + 1e-9)))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--marker-dir", type=Path, required=True,
                        help="GTDB-Tk identify intermediate_results/marker_genes/<genome>")
    parser.add_argument("--partitions", type=Path, required=True,
                        help="Task 11 filtered partitions.nex (the retained markers)")
    parser.add_argument("--positions", type=Path, required=True, help="gene-positions.tsv")
    parser.add_argument("--states", type=Path, required=True, help="site_pattern_positions *.genes.tsv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    retained = set(re.findall(r"charset\s+(\S+)\s*=", args.partitions.read_text()))
    markers = gtdbtk_markers(args.marker_dir, retained)
    core = [r for r in csv.DictReader(args.positions.open(), delimiter="\t") if r["placed"] == "True"]
    states = {r["locus"]: r for r in csv.DictReader(args.states.open(), delimiter="\t")}
    marker_loci = {}
    for m in markers:
        overlapping = [r["locus"] for r in core if r["contig"] == m["contig"]
                       and int(r["start"]) <= m["end"] and int(r["end"]) >= m["start"]]
        m["core_gene"] = ";".join(overlapping)
        m["state"] = states[overlapping[0]]["state"] if len(overlapping) == 1 and overlapping[0] in states else ""
        for locus in overlapping:
            marker_loci[locus] = m["marker"]
    in_markers = [states[l] for l in marker_loci if l in states]
    others = [r for l, r in states.items() if l not in marker_loci]
    count = lambda rows, s: sum(1 for r in rows if r["state"] == s)
    a, b = count(in_markers, "A"), count(in_markers, "B")
    c, d = count(others, "A"), count(others, "B")
    ma = sum(int(r["A"]) for r in in_markers)
    mb = sum(int(r["B"]) for r in in_markers)
    result = {"markers_retained": len(retained), "marker_hits_on_reference": len(markers),
              "markers_matched_to_one_core_gene": sum(1 for m in markers if m["state"]),
              "marker_core_genes": {"A": a, "B": b, "tie": count(in_markers, "tie")},
              "other_core_genes": {"A": c, "B": d, "tie": count(others, "tie")},
              "share_A_of_deciding": {"markers": round(a / (a + b), 4) if a + b else None,
                                      "other_core_genes": round(c / (c + d), 4) if c + d else None},
              "fisher_two_sided_p": fisher_two_sided(a, b, c, d),
              "marker_core_gene_columns": {"A": ma, "B": mb,
                                           "share_A": round(ma / (ma + mb), 4) if ma + mb else None}}
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "markers-in-core-genes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(markers[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(markers)
    (args.output / "markers-in-core-genes.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
