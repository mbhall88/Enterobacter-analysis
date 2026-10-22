"""Measure whether an archived all-versus-all ANI table still covers a study cohort.

`compare_ani_to_archive.py` asks whether the archived *values* are reproducible. This script asks
the prior question: does the archived table cover the cohort at all, and do its inputs still exist?
A table can be numerically perfect and still be unusable, either because its coverage is
systematically missing the genomes that matter or because the files it was computed from are gone.

Three things are reported separately, because they have different remedies:

genome coverage
    Which cohort accessions appear in the table, and which table accessions are outside the cohort.
pair coverage
    How many of the cohort's unordered pairs the table actually contains, and per genome how many
    of its pairs are absent. Absent pairs are counted, never read as low identity: a screening
    threshold and a genuine distant comparison look identical in a table that simply omits the row.
input availability
    Whether the input file each table row names can still be found. The table records paths; if
    that directory no longer exists, the archived numbers cannot be tied to any sequence held now.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import re
from collections import Counter
from pathlib import Path

ACCESSION = re.compile(r"(GC[AF]_\d+\.\d+)")


def accession_of(text: str) -> str:
    match = ACCESSION.search(text)
    return match.group(1) if match else ""


def read_archive(path: Path, ref_column: str = "Ref_file",
                 query_column: str = "Query_file") -> tuple[set[str], set[frozenset[str]],
                                                             Counter, int]:
    """Return the table's accessions, its unordered non-self pairs, its path prefixes and row count."""
    accessions: set[str] = set()
    pairs: set[frozenset[str]] = set()
    prefixes: Counter = Counter()
    rows = 0
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            left, right = accession_of(row[ref_column]), accession_of(row[query_column])
            for field in (row[ref_column], row[query_column]):
                prefixes[field.rsplit("/", 1)[0] if "/" in field else "<no directory>"] += 1
            if not left or not right:
                continue
            accessions.update((left, right))
            if left != right:
                pairs.add(frozenset((left, right)))
    return accessions, pairs, prefixes, rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--archived-table", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True,
                        help="study manifest TSV with an accession_version column")
    parser.add_argument("--cohort-file-columns", default="resolved_chromosome_path,chromosome_path",
                        help="comma-separated manifest columns, first non-empty wins")
    parser.add_argument("--input-root", type=Path,
                        help="directory the archived table's relative paths were written against; "
                             "its existence is reported, not assumed")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    columns = [column.strip() for column in args.cohort_file_columns.split(",")]
    cohort_rows = list(csv.DictReader(args.cohort.open(), delimiter="\t"))
    cohort = [row["accession_version"] for row in cohort_rows]
    cohort_files = {
        row["accession_version"]: next((row[column] for column in columns if row.get(column)), "")
        for row in cohort_rows}
    cohort_set = set(cohort)

    accessions, pairs, prefixes, rows = read_archive(args.archived_table)
    in_cohort_pairs = {pair for pair in pairs if pair <= cohort_set}
    expected = len(cohort) * (len(cohort) - 1) // 2

    missing_by_genome: Counter = Counter()
    for left, right in itertools.combinations(sorted(cohort_set), 2):
        if frozenset((left, right)) not in in_cohort_pairs:
            missing_by_genome[left] += 1
            missing_by_genome[right] += 1

    per_genome = []
    for accession in sorted(cohort_set):
        missing = missing_by_genome.get(accession, 0)
        per_genome.append({
            "accession_version": accession,
            "in_archived_table": str(accession in accessions),
            "pairs_expected": len(cohort) - 1,
            "pairs_missing_from_archive": missing,
            "fraction_missing": f"{missing / (len(cohort) - 1):.4f}",
            "all_pairs_missing": str(missing == len(cohort) - 1),
            "cohort_input_file": cohort_files[accession],
        })
    per_genome.sort(key=lambda row: (-row["pairs_missing_from_archive"], row["accession_version"]))
    header = list(per_genome[0])
    with (args.output / "archive-pair-coverage-by-genome.tsv").open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in per_genome:
            handle.write("\t".join(str(row[column]) for column in header) + "\n")

    referenced_directories = {}
    for prefix in prefixes:
        if prefix == "<no directory>":
            continue
        candidate = (args.input_root / prefix) if args.input_root else Path(prefix)
        referenced_directories[prefix] = {
            "rows_referencing": prefixes[prefix],
            "resolved_path": str(candidate),
            "exists_now": candidate.is_dir(),
        }

    # Can each cohort accession's archived input still be found in a directory the table names?
    findable = set()
    for accession in cohort_set:
        for directory in referenced_directories.values():
            root = Path(directory["resolved_path"])
            if root.is_dir() and any(root.glob(f"{accession}_*")):
                findable.add(accession)
                break

    summary = {
        "archived_table": str(args.archived_table),
        "archived_rows": rows,
        "archived_distinct_accessions": len(accessions),
        "archived_unordered_non_self_pairs": len(pairs),
        "cohort_manifest": str(args.cohort),
        "cohort_size": len(cohort),
        "cohort_accessions_absent_from_table": sorted(cohort_set - accessions),
        "table_accessions_outside_cohort": sorted(accessions - cohort_set),
        "cohort_pairs_expected": expected,
        "cohort_pairs_present_in_archive": len(in_cohort_pairs),
        "cohort_pairs_missing_from_archive": expected - len(in_cohort_pairs),
        "cohort_pairs_missing_fraction": round((expected - len(in_cohort_pairs)) / expected, 6),
        "genomes_with_every_pair_missing": [row["accession_version"] for row in per_genome
                                            if row["all_pairs_missing"] == "True"],
        "genomes_with_any_pair_missing": sum(1 for row in per_genome
                                             if row["pairs_missing_from_archive"] > 0),
        "referenced_input_directories": referenced_directories,
        "cohort_inputs_findable_in_referenced_directories": len(findable),
        "cohort_inputs_not_findable_in_referenced_directories": sorted(cohort_set - findable),
        "note": ("a pair absent from the table is absent, not a low ANI; the archived table's own "
                 "screening threshold and a genuinely distant comparison are indistinguishable "
                 "once the row is omitted"),
    }

    (args.output / "archive-pair-coverage-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
