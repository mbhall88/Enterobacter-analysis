#!/usr/bin/env python3
"""Resolve explicit whole-assembly/chromosome file gaps without editing source data."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path):
    header = None
    sequence = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper()
                header, sequence = line[1:], []
            elif line.strip():
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                sequence.append("".join(line.split()))
    if header is None:
        raise ValueError(f"{path}: no FASTA records")
    yield header, "".join(sequence).upper()


def accession_token(header):
    return header.split()[0]


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def write_fasta(path, header, sequence):
    with path.open("w") as handle:
        handle.write(f">{header}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--sequence-reports", required=True, type=Path)
    parser.add_argument("--download-package", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-total", type=int, default=27)
    parser.add_argument("--expected-existing", type=int, default=16)
    parser.add_argument("--expected-derived", type=int, default=11)
    args = parser.parse_args()
    manifest_path = args.output / "resolved-file-manifest.tsv"
    chromosome_dir = args.output / "chromosomes"
    if manifest_path.exists() or chromosome_dir.exists():
        parser.error("refusing to overwrite resolved outputs")

    with args.status.open(newline="") as handle:
        status_rows = list(csv.DictReader(handle, delimiter="\t"))
    sequence_reports = defaultdict(list)
    for report in map(json.loads, args.sequence_reports.open()):
        sequence_reports[report["assembly_accession"]].append(report)
    assembly_report_path = args.download_package / "ncbi_dataset/data/assembly_data_report.jsonl"
    downloaded_reports = {r["accession"]: r for r in map(json.loads, assembly_report_path.open())}
    chromosome_dir.mkdir(parents=True)

    output_rows = []
    for row in status_rows:
        accession = row["accession"]
        chromosomes = [r for r in sequence_reports[accession] if r.get("assigned_molecule_location_type") == "Chromosome"]
        if len(chromosomes) != 1:
            raise ValueError(f"{accession}: expected one NCBI chromosome, found {len(chromosomes)}")
        component = chromosomes[0]
        expected_accession = component["genbank_accession"]
        expected_length = int(component["length"])

        downloaded_whole = ""
        if accession in downloaded_reports:
            report = downloaded_reports[accession]
            info = report["assemblyInfo"]
            downloaded_whole = args.download_package / "ncbi_dataset/data" / accession / f"{accession}_{info['assemblyName']}_genomic.fna"
            if not downloaded_whole.is_file():
                raise FileNotFoundError(downloaded_whole)
        if row["local_whole_assembly"] == "true":
            whole_path = Path(row["local_whole_path"])
            whole_origin = "student_archive"
        else:
            whole_path = downloaded_whole
            whole_origin = "ncbi_download_2026-09-15"
        if not whole_path or not whole_path.is_file():
            raise FileNotFoundError(f"{accession}: no whole assembly")

        source_matches = [record for record in fasta_records(whole_path) if accession_token(record[0]) == expected_accession]
        if len(source_matches) != 1:
            raise ValueError(f"{accession}: expected one {expected_accession} record in {whole_path}, found {len(source_matches)}")
        source_header, source_sequence = source_matches[0]
        if len(source_sequence) != expected_length:
            raise ValueError(f"{accession}: NCBI length {expected_length}, observed {len(source_sequence)}")

        if row["local_chromosome_extraction"] == "true":
            chromosome_path = Path(row["local_chromosome_path"])
            existing = list(fasta_records(chromosome_path))
            if len(existing) != 1:
                raise ValueError(f"{accession}: local chromosome has {len(existing)} FASTA records")
            chromosome_header, chromosome_sequence = existing[0]
            if chromosome_sequence != source_sequence:
                raise ValueError(f"{accession}: local chromosome differs from downloaded named chromosome")
            chromosome_origin = "student_archive_existing"
        else:
            chromosome_path = chromosome_dir / f"{accession}_{expected_accession}_chromosome.fna"
            write_fasta(chromosome_path, source_header, source_sequence)
            chromosome_header, chromosome_sequence = source_header, source_sequence
            chromosome_origin = "derived_exact_accession"

        output_rows.append({
            "accession": accession,
            "question_ids": row["question_ids"],
            "whole_path": str(whole_path.resolve()),
            "whole_origin": whole_origin,
            "whole_sha256": sha256(whole_path),
            "chromosome_path": str(chromosome_path.resolve()),
            "chromosome_origin": chromosome_origin,
            "chromosome_file_sha256": sha256(chromosome_path),
            "chromosome_sequence_sha256": sequence_sha256(chromosome_sequence),
            "chromosome_accession": expected_accession,
            "chromosome_length": len(chromosome_sequence),
            "chromosome_header": chromosome_header,
            "ncbi_genome_notes": row["ncbi_genome_notes"],
            "validation": "pass",
        })

    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    expected = {
        "student_archive_existing": args.expected_existing,
        "derived_exact_accession": args.expected_derived,
    }
    observed = {kind: sum(row["chromosome_origin"] == kind for row in output_rows) for kind in expected}
    if len(output_rows) != args.expected_total or observed != expected:
        raise SystemExit(f"unexpected resolution count: rows={len(output_rows)}, origins={observed}")


if __name__ == "__main__":
    main()
