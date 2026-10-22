#!/usr/bin/env python3
"""Pull the Verticall rows a focused analysis needs out of a full pairwise table, in one pass.

Task 11's 1,089-genome table is 16 GB because every directed row carries the coordinates of both
assemblies' vertical, horizontal and unaligned regions.  Task 17 needs two views of it:

* ``--summary-out``: every directed row whose two assemblies are both in ``--ids``, with the six
  region columns dropped.  This is the r/m, fraction and distance table for whole subsets.
* ``--regions-out``: the full row, regions included, when assembly_a is in ``--region-ids`` and
  assembly_b in ``--ids``.  This is what the positional analysis reads.

Only ``result_level == primary`` rows are kept -- the row Verticall's matrix uses under
``--multi first`` -- and secondary-peak rows are counted, not silently merged.

The expected directed pairs are accounted for: every ordered pair of distinct ``--ids`` members
that are present in ``--universe`` must appear exactly once, and a missing or duplicated pair is an
error, not a gap to be filled.  A pair absent from a table is not a distance of zero.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

REGION_COLUMNS = ("assembly_a_vertical_regions", "assembly_a_horizontal_regions",
                  "assembly_a_unaligned_regions", "assembly_b_vertical_regions",
                  "assembly_b_horizontal_regions", "assembly_b_unaligned_regions")


def read_ids(path: Path) -> set[str]:
    return {line.split()[0] for line in path.read_text().splitlines() if line.strip()}


def extract(pairwise: Path, ids: set[str], region_ids: set[str], universe: set[str],
            summary_out: Path, regions_out: Path) -> dict:
    wanted = ids & universe
    expected = len(wanted) * (len(wanted) - 1)
    seen: set[tuple[str, str]] = set()
    duplicates = secondary = 0
    rows_read = region_rows = 0
    started = time.monotonic()
    with pairwise.open(encoding="utf-8") as handle, \
            summary_out.open("w", encoding="utf-8") as summary, \
            regions_out.open("w", encoding="utf-8") as regions:
        header = handle.readline().rstrip("\n").split("\t")
        keep = [i for i, name in enumerate(header) if name not in REGION_COLUMNS]
        level = header.index("result_level")
        if len(keep) != len(header) - len(REGION_COLUMNS):
            raise SystemExit(f"{pairwise}: expected all six region columns in the header")
        summary.write("\t".join(header[i] for i in keep) + "\n")
        regions.write("\t".join(header) + "\n")
        for line in handle:
            rows_read += 1
            a, b, _ = line.split("\t", 2)
            if b not in wanted or a not in wanted:
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[level] != "primary":
                secondary += 1
                continue
            if (a, b) in seen:
                duplicates += 1
                continue
            seen.add((a, b))
            summary.write("\t".join(parts[i] for i in keep) + "\n")
            if a in region_ids:
                regions.write(line)
                region_rows += 1
    missing = expected - len(seen)
    return {"pairwise": str(pairwise), "rows_read": rows_read, "ids": len(ids),
            "ids_in_universe": len(wanted), "ids_not_in_universe": sorted(ids - universe),
            "expected_directed_pairs": expected, "found_directed_pairs": len(seen),
            "missing_directed_pairs": missing, "duplicate_primary_rows": duplicates,
            "non_primary_rows_skipped": secondary,
            "region_rows": region_rows, "seconds": round(time.monotonic() - started, 1)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pairwise", type=Path, required=True)
    parser.add_argument("--ids", type=Path, required=True, help="assemblies whose pairs to keep")
    parser.add_argument("--region-ids", type=Path, required=True,
                        help="assembly_a values whose rows keep their region columns")
    parser.add_argument("--universe", type=Path, required=True,
                        help="ids the table was built on; --ids outside it are reported, not expected")
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--regions-out", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--expected-sha256", help="refuse a table whose hash differs")
    args = parser.parse_args(argv)
    for out in (args.summary_out, args.regions_out, args.record):
        if out.exists():
            raise SystemExit(f"refusing to overwrite {out}")
    if args.expected_sha256:
        got = sha256(args.pairwise)
        if got != args.expected_sha256:
            raise SystemExit(f"{args.pairwise}: sha256 {got} != expected {args.expected_sha256}")
    record = extract(args.pairwise, read_ids(args.ids), read_ids(args.region_ids),
                     read_ids(args.universe), args.summary_out, args.regions_out)
    record["pairwise_sha256_checked"] = bool(args.expected_sha256)
    args.record.write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record, indent=1))
    if record["missing_directed_pairs"] or record["duplicate_primary_rows"]:
        print("error: pair accounting failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
