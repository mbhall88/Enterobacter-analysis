#!/usr/bin/env python3
"""Summarize NCBI evidence for Task 01 public provenance questions."""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


MISSING = {"", "missing", "not collected", "not applicable", "na", "none"}


def norm(value):
    if value is None or str(value).strip().casefold() in MISSING:
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).casefold())


def relation(left, right):
    if not norm(left) and not norm(right):
        return "both_missing"
    if not norm(left) or not norm(right):
        return "one_missing"
    if str(left) == str(right):
        return "exact_match"
    if norm(left) == norm(right):
        return "normalized_match"
    return "different"


def attrs(report):
    sample = report.get("assembly_info", {}).get("biosample", {})
    return {item["name"]: item.get("value", "") for item in sample.get("attributes", [])}


def flatten(report):
    info = report.get("assembly_info", {})
    sample = info.get("biosample", {})
    sample_attrs = attrs(report)
    organism = report.get("organism", {})
    stats = report.get("assembly_stats", {})
    checkm = report.get("checkm_info", {})
    ani = report.get("average_nucleotide_identity", {})
    return {
        "accession": report.get("accession", ""),
        "current_accession": report.get("current_accession", ""),
        "assembly_status": info.get("assembly_status", ""),
        "paired_accession": report.get("paired_accession", ""),
        "organism": organism.get("organism_name", ""),
        "strain": organism.get("infraspecific_names", {}).get("strain", ""),
        "biosample": sample.get("accession", ""),
        "bioproject": info.get("bioproject_accession", ""),
        "submitter": info.get("submitter", ""),
        "sample_owner": sample.get("owner", {}).get("name", ""),
        "sample_submission_date": sample.get("submission_date", ""),
        "collection_date": sample_attrs.get("collection_date", ""),
        "collection_year": sample_attrs.get("collection_date", "")[:4],
        "collected_by": sample_attrs.get("collected_by", ""),
        "geo_loc_name": sample_attrs.get("geo_loc_name", ""),
        "lat_lon": sample_attrs.get("lat_lon", ""),
        "host": sample_attrs.get("host", ""),
        "host_disease": sample_attrs.get("host_disease", ""),
        "isolation_source": sample_attrs.get("isolation_source", ""),
        "assembly_release_date": info.get("release_date", ""),
        "assembly_method": info.get("assembly_method", ""),
        "sequencing_tech": info.get("sequencing_tech", ""),
        "component_count": stats.get("number_of_component_sequences", ""),
        "total_sequence_length": stats.get("total_sequence_length", ""),
        "contig_n50": stats.get("contig_n50", ""),
        "checkm_completeness": checkm.get("completeness", ""),
        "checkm_contamination": checkm.get("contamination", ""),
        "taxonomy_check_status": ani.get("taxonomy_check_status", ""),
        "genome_notes": ";".join(info.get("genome_notes", [])) or "none",
    }


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--questions", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--assembly-reports", required=True, type=Path)
    parser.add_argument("--sequence-reports", required=True, type=Path)
    parser.add_argument("--routine-reports", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    outputs = [
        args.output / "accession-metadata.tsv",
        args.output / "pair-field-comparisons.tsv",
        args.output / "routine-pair-verification.tsv",
        args.output / "automated-summary.json",
    ]
    if any(path.exists() for path in outputs):
        parser.error("refusing to overwrite existing summary outputs")

    reports = {r["accession"]: r for r in map(json.loads, args.assembly_reports.open())}
    flat = {key: flatten(value) for key, value in reports.items()}
    sequences = defaultdict(list)
    for record in map(json.loads, args.sequence_reports.open()):
        sequences[record["assembly_accession"]].append(record)
    routine = {r["accession"]: r for r in map(json.loads, args.routine_reports.open())}
    with args.questions.open(newline="") as handle:
        questions = list(csv.DictReader(handle, delimiter="\t"))
    with args.pairs.open(newline="") as handle:
        pairs = list(csv.DictReader(handle, delimiter="\t"))

    metadata_fields = list(next(iter(flat.values())))
    write_tsv(outputs[0], [flat[key] for key in sorted(flat)], metadata_fields)

    compare_fields = [
        "organism", "strain", "collection_year", "collection_date", "collected_by",
        "geo_loc_name", "lat_lon", "host", "host_disease", "isolation_source",
        "biosample", "bioproject", "submitter", "sample_owner", "assembly_method",
        "sequencing_tech", "component_count", "total_sequence_length", "contig_n50",
        "checkm_completeness", "checkm_contamination", "taxonomy_check_status",
    ]
    comparison_rows = []
    for pair in pairs:
        left, right = pair["left_accession"], pair["right_accession"]
        for field in compare_fields:
            comparison_rows.append({
                "pair_id": pair["pair_id"], "question_ids": pair["question_ids"],
                "left_accession": left, "right_accession": right, "field": field,
                "left_value": flat[left][field], "right_value": flat[right][field],
                "automated_relation": relation(flat[left][field], flat[right][field]),
            })
        left_lengths = sorted(int(x["length"]) for x in sequences[left])
        right_lengths = sorted(int(x["length"]) for x in sequences[right])
        left_accessions = sorted(x.get("genbank_accession", "") for x in sequences[left])
        right_accessions = sorted(x.get("genbank_accession", "") for x in sequences[right])
        for field, lv, rv in [
            ("component_length_multiset", ";".join(map(str, left_lengths)), ";".join(map(str, right_lengths))),
            ("component_genbank_accessions", ";".join(left_accessions), ";".join(right_accessions)),
        ]:
            comparison_rows.append({
                "pair_id": pair["pair_id"], "question_ids": pair["question_ids"],
                "left_accession": left, "right_accession": right, "field": field,
                "left_value": lv, "right_value": rv, "automated_relation": relation(lv, rv),
            })
    comparison_fields = list(comparison_rows[0])
    write_tsv(outputs[1], comparison_rows, comparison_fields)

    routine_rows = []
    for question in questions:
        if question["rule"] != "multiple_accession_representations_of_one_assembly":
            continue
        gca, expected_gcf = question["entity_ids"].split(";")
        report = routine[gca]
        observed = report.get("paired_accession", "")
        status = report.get("assembly_info", {}).get("assembly_status", "")
        current = report.get("current_accession", "")
        routine_rows.append({
            "question_id": question["question_id"], "gca_accession": gca,
            "expected_gcf_accession": expected_gcf, "observed_paired_accession": observed,
            "assembly_status": status, "current_accession": current,
            "verified": str(observed == expected_gcf and status == "current" and current == gca).lower(),
        })
    write_tsv(outputs[2], routine_rows, list(routine_rows[0]))

    summary = {
        "assembly_reports": len(reports),
        "ambiguous_pairs": len(pairs),
        "pair_comparison_rows": len(comparison_rows),
        "routine_pair_questions": len(routine_rows),
        "routine_pair_failures": sum(row["verified"] != "true" for row in routine_rows),
    }
    outputs[3].write_text(json.dumps(summary, indent=2) + "\n")
    if summary != {
        "assembly_reports": 9, "ambiguous_pairs": 4, "pair_comparison_rows": 96,
        "routine_pair_questions": 52, "routine_pair_failures": 0,
    }:
        raise SystemExit(f"unexpected summary: {summary}")


if __name__ == "__main__":
    main()
