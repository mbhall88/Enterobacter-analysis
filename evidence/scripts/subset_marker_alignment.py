#!/usr/bin/env python3
"""Cut Task 11's occupancy-filtered BAC120 alignment down to one frozen focused subset.

The filtered alignment (SET-16) holds 1,090 sequences because it predates SET-20.  A subset is cut
from it by id, never by position, and the output is verified against the subset's ids file.

A genome outside the cohort (the second *E. sichuanensis* assembly) cannot be cut from it.  Its
untrimmed GTDB-Tk alignment has the same 41,084 columns as Task 11's, because GTDB-Tk aligns every
genome to the same marker HMMs; Task 11's ``column-map.tsv`` records which of those columns the
occupancy filter kept.  The projection is checked before it is trusted: every cohort tip in the
subset is projected from Task 11's untrimmed alignment through the same map and must reproduce its
filtered row exactly.

Per-marker site counts are written for the subset and for any ``--panel`` ids, which is what the
line-4 marker ledger reads.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, list[str]] = {}
    order: list[str] = []
    with opener(path, "rt") as handle:
        name = None
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise SystemExit(f"{path}: duplicate sequence {name}")
                records[name] = []
                order.append(name)
            elif name is not None:
                records[name].append(line.strip())
    return {n: "".join(records[n]) for n in order}


def read_ids(path: Path) -> list[str]:
    return [line.split()[0] for line in path.read_text().splitlines() if line.strip()]


def read_partitions(path: Path) -> list[tuple[str, int, int]]:
    parts = []
    for match in re.finditer(r"charset\s+(\S+)\s*=\s*(\d+)-(\d+)\s*;", path.read_text()):
        parts.append((match.group(1), int(match.group(2)), int(match.group(3))))
    if not parts:
        raise SystemExit(f"{path}: no charsets")
    return parts


def column_stats(rows: list[str], start: int, end: int) -> dict:
    variable = informative = 0
    for col in range(start - 1, end):
        counts: dict[str, int] = {}
        for row in rows:
            residue = row[col].upper()
            if residue in STANDARD:
                counts[residue] = counts.get(residue, 0) + 1
        if len(counts) > 1:
            variable += 1
            if sum(1 for c in counts.values() if c >= 2) >= 2:
                informative += 1
    return {"columns": end - start + 1, "variable": variable, "parsimony_informative": informative}


def occupancy(row: str) -> float:
    return sum(1 for r in row if r.upper() in STANDARD) / len(row)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--filtered", type=Path, required=True)
    parser.add_argument("--untrimmed", type=Path, required=True,
                        help="Task 11's untrimmed user_msa, used to check the projection")
    parser.add_argument("--column-map", type=Path, required=True)
    parser.add_argument("--partitions", type=Path, required=True)
    parser.add_argument("--ids", type=Path, required=True)
    parser.add_argument("--panel", type=Path, help="a second id set for per-marker counts")
    parser.add_argument("--extra-untrimmed", type=Path, action="append", default=[],
                        help="untrimmed GTDB-Tk alignment of genomes outside the cohort")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)

    ids = read_ids(args.ids)
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate ids")
    filtered = read_fasta(args.filtered)
    untrimmed = read_fasta(args.untrimmed)
    width = {len(s) for s in filtered.values()}
    if len(width) != 1:
        raise SystemExit("filtered alignment is ragged")
    (fwidth,) = width

    with args.column_map.open(newline="") as handle:
        cmap = list(csv.DictReader(handle, delimiter="\t"))
    kept = [int(r["original_column"]) - 1 for r in cmap if r["decision"] == "retained"]
    if len(kept) != fwidth:
        raise SystemExit(f"column map keeps {len(kept)} columns, alignment has {fwidth}")
    uwidth = len(cmap)

    def project(row: str) -> str:
        if len(row) != uwidth:
            raise SystemExit(f"untrimmed row has {len(row)} columns, map expects {uwidth}")
        return "".join(row[i] for i in kept)

    extras: dict[str, str] = {}
    for path in args.extra_untrimmed:
        for name, row in read_fasta(path).items():
            if name in filtered:
                raise SystemExit(f"{name} is already a cohort tip; do not project it")
            extras[name] = project(row)

    checked = 0
    for tip in ids:
        if tip in filtered:
            if tip not in untrimmed:
                raise SystemExit(f"{tip} missing from the untrimmed alignment")
            if project(untrimmed[tip]) != filtered[tip]:
                raise SystemExit(f"projection check failed for {tip}")
            checked += 1

    out_rows = {}
    for tip in ids:
        if tip in filtered:
            out_rows[tip] = filtered[tip]
        elif tip in extras:
            out_rows[tip] = extras[tip]
        else:
            raise SystemExit(f"{tip} is in neither the filtered alignment nor an extra alignment")
    with (args.output / "alignment.faa").open("w") as handle:
        for tip in ids:
            handle.write(f">{tip}\n{out_rows[tip]}\n")
    (args.output / "partitions.nex").write_text(args.partitions.read_text())

    parts = read_partitions(args.partitions)
    if parts[-1][2] != fwidth:
        raise SystemExit("partitions do not end at the alignment width")
    panel = read_ids(args.panel) if args.panel else []
    stats_rows = []
    for name, start, end in parts:
        row = {"marker": name, "start": start, "end": end}
        for label, tips in (("subset", ids), ("panel", panel)):
            if not tips:
                continue
            seqs = [out_rows[t] for t in tips]
            s = column_stats(seqs, start, end)
            present = sum(1 for q in seqs if any(c.upper() in STANDARD for c in q[start - 1:end]))
            row.update({f"{label}_{k}": v for k, v in s.items() if k != "columns"})
            row[f"{label}_tips_with_residues"] = present
        row["columns"] = end - start + 1
        stats_rows.append(row)
    with (args.output / "per-marker-sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stats_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(stats_rows)

    whole = column_stats(list(out_rows.values()), 1, fwidth)
    digest = hashlib.sha256((args.output / "alignment.faa").read_bytes()).hexdigest()
    record = {
        "filtered": str(args.filtered), "untrimmed": str(args.untrimmed),
        "ids": str(args.ids), "tips": len(ids), "columns": fwidth, **whole,
        "extra_tips": sorted(t for t in ids if t in extras), "projection_checked_on_cohort_tips": checked,
        "cohort_alignment_tips": len(filtered),
        "min_tip_occupancy": min((occupancy(r), t) for t, r in out_rows.items()),
        "alignment_sha256": digest,
    }
    (args.output / "subset-record.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
