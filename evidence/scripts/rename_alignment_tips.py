#!/usr/bin/env python3
"""Rename alignment tips through an explicit two-column map, with no partial renames.

GTDB-Tk refuses genome ids that collide with its own reference genome ids, so Task 07 submits the
frozen panel under prefixed ids and maps them back here.  A tree whose tips do not match the frozen
panel exactly fails this task's completion criteria, so this refuses anything ambiguous: an
unmapped name, a map entry that matches nothing, or two names mapping to one output name.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_map(path: Path, from_col: str, to_col: str) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path} is empty")
    for col in (from_col, to_col):
        if col not in rows[0]:
            raise ValueError(f"{path} has no column {col!r}")
    mapping: dict[str, str] = {}
    for row in rows:
        src, dst = row[from_col], row[to_col]
        if src in mapping:
            raise ValueError(f"duplicate {from_col} in {path}: {src}")
        mapping[src] = dst
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(f"{to_col} values in {path} are not unique")
    return mapping


def rename_fasta(in_path: Path, out_path: Path, mapping: dict[str, str]) -> dict:
    seen: set[str] = set()
    written = 0
    with open_maybe_gzip(in_path) as src, open_maybe_gzip(out_path, "wt") as dst:
        for line in src:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if name not in mapping:
                    raise ValueError(
                        f"{in_path} contains {name!r}, which the map does not cover; refusing a "
                        "partial rename"
                    )
                if name in seen:
                    raise ValueError(f"duplicate sequence {name!r} in {in_path}")
                seen.add(name)
                dst.write(f">{mapping[name]}\n")
                written += 1
            else:
                dst.write(line)
    unused = sorted(set(mapping) - seen)
    return {"sequences_renamed": written, "map_entries_unused": unused}


def rename_newick(in_path: Path, out_path: Path, mapping: dict[str, str]) -> dict:
    """Replace whole tip labels in a Newick string.

    Labels are replaced longest-first so that no id that is a prefix of another is corrupted, and
    each replacement is counted so an unused map entry is visible.
    """
    text = in_path.read_text()
    counts = {}
    for src in sorted(mapping, key=len, reverse=True):
        dst = mapping[src]
        n = text.count(src)
        if n:
            text = text.replace(src, dst)
        counts[src] = n
    out_path.write_text(text)
    return {
        "labels_replaced": sum(1 for n in counts.values() if n),
        "map_entries_unused": sorted(k for k, n in counts.items() if n == 0),
        "labels_seen_more_than_once": sorted(k for k, n in counts.items() if n > 1),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--map", required=True, type=Path)
    parser.add_argument("--from-column", default="gtdbtk_genome_id")
    parser.add_argument("--to-column", default="tip_id")
    parser.add_argument("--format", choices=["fasta", "newick"], default="fasta")
    parser.add_argument("--allow-unused-map-entries", action="store_true",
                        help="permit map rows that match nothing in the input")
    args = parser.parse_args(argv)

    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")
    mapping = read_map(args.map, args.from_column, args.to_column)
    if args.format == "fasta":
        report = rename_fasta(args.input, args.output, mapping)
    else:
        report = rename_newick(args.input, args.output, mapping)
    if report["map_entries_unused"] and not args.allow_unused_map_entries:
        args.output.unlink()
        raise SystemExit(
            f"{len(report['map_entries_unused'])} map entries matched nothing in {args.input}: "
            f"{report['map_entries_unused'][:10]}"
        )
    print(f"{args.input} -> {args.output}: {report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
