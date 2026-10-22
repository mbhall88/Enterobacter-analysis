#!/usr/bin/env python3
"""Fetch NCBI Datasets assembly reports with a reproducible command record."""

import argparse
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", required=True, type=Path)
    parser.add_argument("--accessions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", choices=("genome", "sequence"), default="genome")
    args = parser.parse_args()

    if not args.datasets.is_file():
        parser.error(f"datasets executable is absent: {args.datasets}")
    accessions = [line.strip() for line in args.accessions.read_text().splitlines() if line.strip()]
    if not accessions:
        parser.error("accession file is empty")
    if len(accessions) != len(set(accessions)):
        parser.error("accession file contains duplicates")

    args.output.mkdir(parents=True, exist_ok=True)
    report_name = "assembly-reports.jsonl" if args.report == "genome" else "sequence-reports.jsonl"
    managed = [report_name, "datasets.stderr", "fetch-record.json"]
    if any((args.output / name).exists() for name in managed):
        parser.error("refusing to overwrite an existing fetch record")

    version = subprocess.run(
        [str(args.datasets), "version"], text=True, capture_output=True, check=True
    )
    command = [
        str(args.datasets),
        "summary",
        "genome",
        "accession",
        *accessions,
        "--report",
        args.report,
        "--as-json-lines",
    ]
    started = datetime.now(timezone.utc)
    before = time.monotonic()
    result = subprocess.run(command, text=True, capture_output=True)
    elapsed = time.monotonic() - before
    (args.output / report_name).write_text(result.stdout)
    (args.output / "datasets.stderr").write_text(result.stderr)

    parsed = []
    parse_error = None
    try:
        parsed = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
    except json.JSONDecodeError as exc:
        parse_error = str(exc)
    accession_field = "accession" if args.report == "genome" else "assembly_accession"
    returned = [report.get(accession_field) for report in parsed]
    returned_unique = sorted(set(returned))
    record = {
        "retrieved_at_utc": started.isoformat(),
        "elapsed_seconds": round(elapsed, 3),
        "datasets_version": version.stdout.strip(),
        "report": args.report,
        "command": command,
        "exit_status": result.returncode,
        "requested_accessions": accessions,
        "returned_accessions": returned_unique,
        "missing_accessions": sorted(set(accessions) - set(returned_unique)),
        "unexpected_accessions": sorted(set(returned_unique) - set(accessions)),
        "parse_error": parse_error,
        "input_sha256": {str(args.accessions): sha256(args.accessions)},
    }
    (args.output / "fetch-record.json").write_text(json.dumps(record, indent=2) + "\n")

    if result.returncode or parse_error or set(returned_unique) != set(accessions):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
