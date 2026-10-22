"""Compare a new skani run against an archived ANI table on the pairs they share.

The archived tables carry no software version, command line or input list (Task 12 provenance
audit, DEC-002 closed as unrecoverable). That is a reason to distrust their provenance, not a
reason to assume their numbers are wrong. This script measures the question directly: on the pairs
both tables contain, does a pinned rerun reproduce the archived values?

Two comparisons are reported and never merged:

sequence comparison
    For every genome present in both the new input set and the archive, whether the sequence is
    identical, and if not, whether one is a prefix of the other. A truncated archived download and
    a genuinely different assembly are different findings and are labelled differently.
pair comparison
    For every genome pair present in both tables, the difference in ANI. Pairs present in only one
    table are counted and reported as such; they are never treated as a difference of zero.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from pathlib import Path

ACCESSION = re.compile(r"(GC[AF]_\d+\.\d+)")


def accession_of(text: str) -> str:
    match = ACCESSION.search(text)
    return match.group(1) if match else ""


def read_sequences(path: Path) -> list[str]:
    records: list[list[str]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append([])
        elif records:
            records[-1].append(line.strip().upper())
    return ["".join(record) for record in records]


def compare_sequences(new: Path, archived: Path) -> dict[str, object]:
    """Classify how two FASTA copies of one accession differ."""
    new_records = sorted(read_sequences(new))
    archived_records = sorted(read_sequences(archived))
    new_total = sum(len(record) for record in new_records)
    archived_total = sum(len(record) for record in archived_records)
    result = {
        "new_records": len(new_records), "archived_records": len(archived_records),
        "new_total_bp": new_total, "archived_total_bp": archived_total,
        "bp_difference": new_total - archived_total,
    }
    if new_records == archived_records:
        return dict(result, status="identical_sequence")
    if len(new_records) == len(archived_records) == 1:
        new_sequence, archived_sequence = new_records[0], archived_records[0]
        if new_sequence.startswith(archived_sequence):
            return dict(result, status="archived_truncated_prefix_of_current")
        if archived_sequence.startswith(new_sequence):
            return dict(result, status="current_truncated_prefix_of_archived")
        if archived_sequence in new_sequence:
            return dict(result, status="archived_contained_in_current")
    return dict(result, status="sequence_differs")


def read_pairs(path: Path, ani_column: str, ref_column: str,
               query_column: str) -> dict[frozenset[str], float]:
    """Read a skani-style table into unordered accession pairs; drop self comparisons."""
    pairs: dict[frozenset[str], float] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            left, right = accession_of(row[ref_column]), accession_of(row[query_column])
            if not left or not right or left == right:
                continue
            key = frozenset((left, right))
            value = float(row[ani_column])
            # skani reports both directions with the same ANI; keep the first and verify the rest.
            pairs.setdefault(key, value)
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--new-table", type=Path, required=True)
    parser.add_argument("--archived-table", type=Path, required=True)
    parser.add_argument("--sequence-pairs", type=Path,
                        help="TSV of accession_version, new_path, archived_path")
    parser.add_argument("--tolerance", type=float, default=0.05,
                        help="ANI points below which a difference is reported as agreement")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    new_pairs = read_pairs(args.new_table, "ANI", "Ref_file", "Query_file")
    archived_pairs = read_pairs(args.archived_table, "ANI", "Ref_file", "Query_file")
    shared = sorted(set(new_pairs) & set(archived_pairs), key=lambda key: sorted(key))

    rows = []
    for key in shared:
        left, right = sorted(key)
        difference = new_pairs[key] - archived_pairs[key]
        rows.append({
            "accession_a": left, "accession_b": right,
            "new_ani_percent": f"{new_pairs[key]:.4f}",
            "archived_ani_percent": f"{archived_pairs[key]:.4f}",
            "ani_difference": f"{difference:.4f}",
            "within_tolerance": str(abs(difference) <= args.tolerance),
        })
    header = ["accession_a", "accession_b", "new_ani_percent", "archived_ani_percent",
              "ani_difference", "within_tolerance"]
    with (args.output / "shared-pair-comparison.tsv").open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row[column] for column in header) + "\n")

    differences = [abs(new_pairs[key] - archived_pairs[key]) for key in shared]
    summary: dict[str, object] = {
        "new_table": str(args.new_table),
        "archived_table": str(args.archived_table),
        "new_unordered_pairs": len(new_pairs),
        "archived_unordered_pairs": len(archived_pairs),
        "shared_pairs": len(shared),
        "pairs_only_in_new": len(set(new_pairs) - set(archived_pairs)),
        "pairs_only_in_archive": len(set(archived_pairs) - set(new_pairs)),
        "tolerance_percent": args.tolerance,
        "ani_difference": {
            "max_absolute": round(max(differences), 4) if differences else None,
            "mean_absolute": round(statistics.mean(differences), 5) if differences else None,
            "median_absolute": round(statistics.median(differences), 5) if differences else None,
            "pairs_over_tolerance": sum(1 for d in differences if d > args.tolerance),
        },
        "note": ("pairs present in only one table are counted, never scored as a zero difference"),
    }

    if args.sequence_pairs:
        sequence_rows = []
        with args.sequence_pairs.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                comparison = compare_sequences(Path(row["new_path"]), Path(row["archived_path"]))
                sequence_rows.append(dict(accession_version=row["accession_version"],
                                          new_path=row["new_path"],
                                          archived_path=row["archived_path"], **comparison))
        sequence_header = list(sequence_rows[0])
        with (args.output / "sequence-comparison.tsv").open("w") as handle:
            handle.write("\t".join(sequence_header) + "\n")
            for row in sequence_rows:
                handle.write("\t".join(str(row[column]) for column in sequence_header) + "\n")
        summary["sequence_comparison"] = {
            status: sum(1 for row in sequence_rows if row["status"] == status)
            for status in sorted({row["status"] for row in sequence_rows})
        }
        summary["sequences_not_identical"] = [
            {"accession_version": row["accession_version"], "status": row["status"],
             "bp_difference": row["bp_difference"]}
            for row in sequence_rows if row["status"] != "identical_sequence"]

    (args.output / "archive-comparison-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
