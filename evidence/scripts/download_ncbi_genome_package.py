#!/usr/bin/env python3
"""Download, validate and unpack an NCBI Datasets genome package."""

import argparse
import hashlib
import json
import subprocess
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", required=True, type=Path)
    parser.add_argument("--accessions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if not args.datasets.is_file():
        parser.error(f"datasets executable is absent: {args.datasets}")
    requested = [line.strip() for line in args.accessions.read_text().splitlines() if line.strip()]
    if not requested or len(requested) != len(set(requested)):
        parser.error("accessions must be a non-empty unique list")
    args.output.mkdir(parents=True, exist_ok=True)
    archive = args.output / "genomes.zip"
    unpacked = args.output / "package"
    managed = [archive, unpacked, args.output / "download.stdout", args.output / "download.stderr", args.output / "download-record.json"]
    if any(path.exists() for path in managed):
        parser.error("refusing to overwrite an existing download record")

    version = subprocess.run([str(args.datasets), "version"], text=True, capture_output=True, check=True)
    command = [
        str(args.datasets), "download", "genome", "accession",
        "--inputfile", str(args.accessions), "--include", "genome,seq-report",
        "--assembly-source", "GenBank", "--filename", str(archive), "--no-progressbar",
    ]
    started = datetime.now(timezone.utc)
    before = time.monotonic()
    result = subprocess.run(command, text=True, capture_output=True)
    elapsed = time.monotonic() - before
    (args.output / "download.stdout").write_text(result.stdout)
    (args.output / "download.stderr").write_text(result.stderr)
    if result.returncode:
        raise SystemExit(result.returncode)

    with zipfile.ZipFile(archive) as package:
        bad = package.testzip()
        if bad:
            raise SystemExit(f"zip checksum failed for {bad}")
        root = unpacked.resolve()
        for member in package.infolist():
            destination = (unpacked / member.filename).resolve()
            if root != destination and root not in destination.parents:
                raise SystemExit(f"unsafe zip member: {member.filename}")
        report_path = "ncbi_dataset/data/assembly_data_report.jsonl"
        reports = [json.loads(line) for line in package.read(report_path).decode().splitlines() if line.strip()]
        returned = [report["accession"] for report in reports]
        package.extractall(unpacked)

    record = {
        "retrieved_at_utc": started.isoformat(),
        "elapsed_seconds": round(elapsed, 3),
        "datasets_version": version.stdout.strip(),
        "command": command,
        "exit_status": result.returncode,
        "requested_accessions": requested,
        "returned_accessions": returned,
        "missing_accessions": sorted(set(requested) - set(returned)),
        "unexpected_accessions": sorted(set(returned) - set(requested)),
        "accessions_sha256": sha256(args.accessions),
        "archive_sha256": sha256(archive),
        "archive_bytes": archive.stat().st_size,
    }
    (args.output / "download-record.json").write_text(json.dumps(record, indent=2) + "\n")
    if set(returned) != set(requested):
        raise SystemExit("download package did not contain exactly the requested accessions")


if __name__ == "__main__":
    main()
