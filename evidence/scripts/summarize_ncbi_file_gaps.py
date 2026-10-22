#!/usr/bin/env python3
"""Join Task 01 file-gap questions to current NCBI component records."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


RULES = {
    "chromosome_extraction_absent",
    "no_local_sequence_for_roster_accession",
    "whole_assembly_file_absent",
}


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--questions", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--assembly-reports", required=True, type=Path)
    parser.add_argument("--sequence-reports", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    table = args.output / "file-gap-status.tsv"
    summary_path = args.output / "file-gap-summary.json"
    if table.exists() or summary_path.exists():
        parser.error("refusing to overwrite existing file-gap summary")

    reports = {r["accession"]: r for r in map(json.loads, args.assembly_reports.open())}
    sequences = defaultdict(list)
    for record in map(json.loads, args.sequence_reports.open()):
        sequences[record["assembly_accession"]].append(record)
    local = defaultdict(lambda: defaultdict(list))
    with args.audit.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            local[row["accession_version"]][row["role"]].append(row["path"])
    question_ids = defaultdict(list)
    with args.questions.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["rule"] not in RULES:
                continue
            for accession in row["entity_ids"].split(";"):
                question_ids[accession].append(row["question_id"])

    rows = []
    for accession in sorted(reports):
        report = reports[accession]
        info = report["assembly_info"]
        chromosomes = [
            sequence for sequence in sequences[accession]
            if sequence.get("assigned_molecule_location_type") == "Chromosome"
        ]
        has_whole = bool(local[accession]["public_whole"])
        has_chromosome = bool(local[accession]["study_chromosome_candidate"])
        if has_whole and not has_chromosome:
            action = "extract_named_chromosome_then_qc"
        elif not has_whole and not has_chromosome:
            action = "retrieve_complete_assembly_then_qc"
        elif not has_whole and has_chromosome:
            action = "chromosome_ready;retrieve_whole_only_if_needed"
        else:
            action = "no_file_gap"
        rows.append({
            "accession": accession,
            "question_ids": ";".join(sorted(question_ids[accession])),
            "assembly_status": info.get("assembly_status", ""),
            "current_accession": report.get("current_accession", ""),
            "assembly_level": info.get("assembly_level", ""),
            "component_count": len(sequences[accession]),
            "chromosome_component_count": len(chromosomes),
            "chromosome_genbank_accession": ";".join(x.get("genbank_accession", "") for x in chromosomes),
            "chromosome_refseq_accession": ";".join(x.get("refseq_accession", "") for x in chromosomes),
            "local_whole_assembly": str(has_whole).lower(),
            "local_chromosome_extraction": str(has_chromosome).lower(),
            "local_whole_path": ";".join(local[accession]["public_whole"]),
            "local_chromosome_path": ";".join(local[accession]["study_chromosome_candidate"]),
            "ncbi_genome_notes": ";".join(info.get("genome_notes", [])),
            "recommended_action": action,
        })
    write_tsv(table, rows)
    counts = defaultdict(int)
    for row in rows:
        counts[row["recommended_action"]] += 1
    summary = {
        "accessions": len(rows),
        "all_current": all(r["assembly_status"] == "current" and r["current_accession"] == r["accession"] for r in rows),
        "all_complete_genome": all(r["assembly_level"] == "Complete Genome" for r in rows),
        "all_have_one_named_chromosome": all(r["chromosome_component_count"] == 1 for r in rows),
        "actions": dict(sorted(counts.items())),
        "accessions_with_ncbi_genome_notes": sum(bool(r["ncbi_genome_notes"]) for r in rows),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    expected_actions = {
        "chromosome_ready;retrieve_whole_only_if_needed": 16,
        "extract_named_chromosome_then_qc": 10,
        "retrieve_complete_assembly_then_qc": 1,
    }
    if not summary["all_current"] or not summary["all_complete_genome"] or not summary["all_have_one_named_chromosome"] or summary["actions"] != expected_actions:
        raise SystemExit(f"unexpected file-gap summary: {summary}")


if __name__ == "__main__":
    main()
