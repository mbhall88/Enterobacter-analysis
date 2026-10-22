#!/usr/bin/env python3
"""Check that splitting ``verticall pairwise`` with ``--part`` changes nothing but the schedule.

The whole full-cohort design rests on this. Verticall 0.4.2 has no resume, so a 2,840 CPU-hour run
is only safe if it can be cut into independently re-runnable parts — and that is only sound if the
parts are a *partition*: every directed pair analysed exactly once, by exactly one part, with the
same result as the undivided run. This script refuses to assume any of that.

Three properties, all required (Michael, DEC-009 Decision 1, 22 September 2026)
------------------------------------------------------------------------------
**1. The split is a true partition.** Every directed pair of the whole run appears in exactly one
part: none twice, none missing, and no part inventing a pair the whole run does not have.

**2. Identical directed-pair counts.** Not only the set of pairs but the number of *rows* per pair:
Verticall writes one row per distance peak, so a pair with a secondary peak has two rows, and a
comparison that collapsed them would hide a difference.

**3. An identical distance matrix.** Compared value by value after aligning names, and separately as
a canonical hash over sorted (name, name, value) triples. Task 09 found that Mashtree's matrix row
order is not reproducible, so no tool's row order is assumed here: a matrix is a mapping from name
pairs to values, and that is what is compared.

A difference in any of the three is a failure. Nothing is tolerated as "close enough": these runs
analyse the same sequence with the same parameters, so the only correct answer is *the same answer*.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

HEADER_FIRST_TWO = ("assembly_a", "assembly_b")


def read_pairwise(path: Path) -> tuple[list[str], Counter, int, int]:
    """Return (header fields, directed-pair row counts, data rows, header lines seen)."""
    header: list[str] = []
    counts: Counter = Counter()
    rows = header_lines = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if tuple(parts[:2]) == HEADER_FIRST_TWO:
                header_lines += 1
                if not header:
                    header = parts
                continue
            counts[(parts[0], parts[1])] += 1
            rows += 1
    return header, counts, rows, header_lines


def read_phylip(path: Path) -> tuple[list[str], dict[tuple[str, str], float]]:
    """Read a PHYLIP distance matrix as a name->name mapping, not as an ordered grid."""
    lines = [line.rstrip("\n") for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise SystemExit(f"{path}: empty matrix")
    try:
        declared = int(lines[0].split()[0])
    except (ValueError, IndexError):
        raise SystemExit(f"{path}: first line is not a taxon count")
    body = lines[1:]
    if len(body) != declared:
        raise SystemExit(f"{path}: declares {declared} taxa, carries {len(body)} rows")
    names = [line.split()[0] for line in body]
    if len(set(names)) != len(names):
        raise SystemExit(f"{path}: repeats a taxon name")
    values: dict[tuple[str, str], float] = {}
    for line in body:
        parts = line.split()
        row_name, cells = parts[0], parts[1:]
        if len(cells) != declared:
            raise SystemExit(f"{path}: row {row_name} has {len(cells)} cells, expected {declared}")
        for column_name, cell in zip(names, cells):
            values[(row_name, column_name)] = float(cell)
    return names, values


def matrix_hash(values: dict[tuple[str, str], float]) -> str:
    """Order-independent hash: sort the triples, then hash. Row order must not matter."""
    digest = hashlib.sha256()
    for (a, b), value in sorted(values.items()):
        digest.update(f"{a}\t{b}\t{value!r}\n".encode())
    return digest.hexdigest()


def compare_matrices(whole: Path, parted: Path, tolerance: float) -> dict:
    whole_names, whole_values = read_phylip(whole)
    part_names, part_values = read_phylip(parted)
    differences = []
    for key in sorted(set(whole_values) | set(part_values)):
        left, right = whole_values.get(key), part_values.get(key)
        if left is None or right is None:
            differences.append({"pair": list(key), "whole": left, "parted": right,
                                "reason": "present in only one matrix"})
        elif not math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance):
            differences.append({"pair": list(key), "whole": left, "parted": right,
                                "reason": "values differ"})
    return {
        "whole_taxa": len(whole_names),
        "parted_taxa": len(part_names),
        "taxon_sets_identical": sorted(whole_names) == sorted(part_names),
        "row_order_identical": whole_names == part_names,
        "cells_compared": len(set(whole_values) | set(part_values)),
        "differing_cells": len(differences),
        "first_differences": differences[:20],
        "whole_matrix_sha256": matrix_hash(whole_values),
        "parted_matrix_sha256": matrix_hash(part_values),
        "identical": not differences,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--whole", type=Path, required=True,
                        help="the --part 1/1 pairwise TSV")
    parser.add_argument("--part", type=Path, action="append", default=[], required=True,
                        help="one partitioned pairwise TSV; repeat once per part")
    parser.add_argument("--whole-matrix", type=Path,
                        help="PHYLIP matrix built from --whole")
    parser.add_argument("--parted-matrix", type=Path,
                        help="PHYLIP matrix built from the merged parts")
    parser.add_argument("--expected-ids", type=Path,
                        help="one tip id per line; every id must appear in the whole run")
    parser.add_argument("--tolerance", type=float, default=0.0,
                        help="absolute tolerance on a matrix cell (default: 0.0, exact)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    whole_header, whole_counts, whole_rows, whole_headers = read_pairwise(args.whole)
    if whole_headers != 1:
        raise SystemExit(f"{args.whole}: {whole_headers} header lines, expected exactly 1")

    per_part = []
    combined: Counter = Counter()
    seen_in: dict[tuple[str, str], list[int]] = {}
    for index, path in enumerate(args.part, start=1):
        header, counts, rows, headers = read_pairwise(path)
        if header and whole_header and header != whole_header:
            raise SystemExit(f"{path}: header differs from {args.whole}")
        expected_headers = 1 if index == 1 else 0
        for pair in counts:
            seen_in.setdefault(pair, []).append(index)
        combined.update(counts)
        per_part.append({
            "part": index, "file": str(path), "data_rows": rows, "distinct_pairs": len(counts),
            "header_lines": headers, "header_lines_expected": expected_headers,
            "header_as_expected": headers == expected_headers,
        })

    duplicated = sorted(pair for pair, parts in seen_in.items() if len(parts) > 1)
    missing = sorted(set(whole_counts) - set(combined))
    extra = sorted(set(combined) - set(whole_counts))
    count_mismatch = sorted(pair for pair in set(whole_counts) & set(combined)
                            if whole_counts[pair] != combined[pair])

    ids_result = None
    if args.expected_ids:
        expected = [line.strip() for line in args.expected_ids.read_text().splitlines()
                    if line.strip()]
        observed = {name for pair in whole_counts for name in pair}
        n = len(expected)
        ids_result = {
            "expected_ids": n,
            "observed_names": len(observed),
            "missing_from_run": sorted(set(expected) - observed),
            "unexpected_in_run": sorted(observed - set(expected)),
            "expected_directed_pairs": n * (n - 1),
            "observed_directed_pairs": len(whole_counts),
            "directed_pairs_complete": len(whole_counts) == n * (n - 1),
        }

    matrix_result = None
    if args.whole_matrix and args.parted_matrix:
        matrix_result = compare_matrices(args.whole_matrix, args.parted_matrix, args.tolerance)

    checks = {
        "every_pair_appears_in_exactly_one_part": not duplicated,
        "no_pair_of_the_whole_run_is_missing": not missing,
        "no_part_invents_a_pair": not extra,
        "directed_pair_row_counts_identical": not count_mismatch,
        "total_rows_identical": whole_rows == sum(combined.values()),
        "headers_as_expected": all(entry["header_as_expected"] for entry in per_part),
    }
    if ids_result:
        checks["every_expected_id_is_present"] = not ids_result["missing_from_run"]
        checks["no_unexpected_id_is_present"] = not ids_result["unexpected_in_run"]
        checks["the_whole_run_is_complete"] = ids_result["directed_pairs_complete"]
    if matrix_result:
        checks["distance_matrices_identical"] = matrix_result["identical"]
        checks["matrix_taxon_sets_identical"] = matrix_result["taxon_sets_identical"]

    record = {
        "whole": {"file": str(args.whole), "data_rows": whole_rows,
                  "distinct_directed_pairs": len(whole_counts), "header_lines": whole_headers},
        "parts": per_part,
        "parts_total": len(args.part),
        "combined": {"data_rows": sum(combined.values()),
                     "distinct_directed_pairs": len(combined)},
        "partition": {
            "pairs_in_more_than_one_part": duplicated[:50],
            "pairs_in_more_than_one_part_count": len(duplicated),
            "pairs_missing_from_every_part": missing[:50],
            "pairs_missing_from_every_part_count": len(missing),
            "pairs_no_part_should_have": extra[:50],
            "pairs_no_part_should_have_count": len(extra),
            "pairs_with_different_row_counts": [
                {"pair": list(pair), "whole": whole_counts[pair], "parted": combined[pair]}
                for pair in count_mismatch[:50]],
            "pairs_with_different_row_counts_count": len(count_mismatch),
        },
        "tip_coverage": ids_result,
        "matrix": matrix_result,
        "checks": checks,
        "passed": all(checks.values()),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"checks": checks, "passed": record["passed"]}, indent=2))
    if not record["passed"]:
        failed = [name for name, ok in checks.items() if not ok]
        print(f"FAILED: {', '.join(failed)}. See {args.output}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
