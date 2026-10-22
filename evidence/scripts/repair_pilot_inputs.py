#!/usr/bin/env python3
"""Produce one ``verticall repair`` derivative per frozen pilot tip, and record both hashes.

``verticall pairwise`` requires unambiguous input, so the Verticall arm of the comparison cannot use
the original FASTAs.  If only that arm used repaired sequence, a topology difference between the
three methods could be the repair rather than the method.  So every method gets the same derivative
set, and this script is where the source-to-derivative mapping is created and recorded.

What repair actually does (Verticall 0.4.2, ``verticall/repair.py``) is worth stating plainly,
because it is not a no-op: each record is upper-cased and then split on every run of non-ACGT
characters, which are **deleted**, and duplicate contig names get ``_1``/``_2`` suffixes.  A single
ambiguous base inside a closed chromosome therefore becomes an extra contig boundary and one lost
base.  This script measures that per genome - bases removed, records gained - so the cost is a
number in a table rather than an assurance.

Safety properties
-----------------
* Sources are opened read-only and a derivative is never written over its source; the script aborts
  if any derivative path collides with any source path.
* Tip ids and derivative paths must both be unique, and every frozen row must get exactly one
  derivative.
* The derivative directory is created fresh per run and must be empty, so no minimap2 index built
  from earlier contents can survive next to different sequence.  ``repair`` itself builds no index;
  the directory this script fills is the one ``verticall pairwise`` will index, which is why it is
  per-run and writable.
* An existing complete record is reused only with ``--resume``, and then only when the recorded
  source hash still matches.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import time
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_stats(path: Path) -> tuple[int, int, int]:
    """Return (records, total bases, ambiguous bases) for a FASTA."""
    records = ambiguous = total = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                seq = line.strip().upper()
                total += len(seq)
                ambiguous += sum(1 for base in seq if base not in "ACGT")
    return records, total, ambiguous


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", "--pilot", dest="manifest", type=Path, required=True,
                        help="the frozen manifest to repair: manifests/pilot.tsv or "
                             "manifests/full.tsv. --pilot is kept as an alias so the Task 08/09 "
                             "command ledgers still run, but nothing here is pilot-specific.")
    parser.add_argument("--verticall", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="run record directory")
    parser.add_argument("--limit", type=int, help="repair only the first N rows (benchmark mode)")
    parser.add_argument("--resume", action="store_true",
                        help="reuse derivatives already recorded against an unchanged source hash")
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    if not rows:
        raise SystemExit(f"{args.manifest}: no rows")
    tips = [row["tip_id"] for row in rows]
    if len(set(tips)) != len(tips):
        raise SystemExit("pilot.tsv has duplicate tip ids; refusing to repair")
    derivatives = [Path(row["derivative_fasta_path"]) for row in rows]
    if len(set(derivatives)) != len(derivatives):
        raise SystemExit("pilot.tsv has duplicate derivative paths; refusing to repair")
    sources = {Path(row["source_fasta_path"]).resolve() for row in rows}
    clashes = sorted(str(path) for path in derivatives if path.resolve() in sources)
    if clashes:
        raise SystemExit(f"derivative paths collide with source FASTAs: {clashes}")

    derivative_dir = derivatives[0].parent
    derivative_dir.mkdir(parents=True, exist_ok=True)
    args.output.mkdir(parents=True, exist_ok=True)
    record_path = args.output / "source-derivative-map.tsv"
    previous = {}
    if args.resume and record_path.is_file():
        previous = {row["tip_id"]: row for row in read_tsv(record_path)}
    elif not args.resume:
        existing = [path for path in derivative_dir.iterdir()] if derivative_dir.is_dir() else []
        if existing and not args.limit:
            raise SystemExit(
                f"{derivative_dir} is not empty ({len(existing)} entries). Use a fresh directory so "
                f"no index or derivative from different content is reused, or pass --resume.")

    version = subprocess.run([str(args.verticall), "--version"], text=True, capture_output=True)
    logs = args.output / "repair-logs"
    logs.mkdir(exist_ok=True)

    selected = rows[:args.limit] if args.limit else rows
    mapped = []
    started = time.monotonic()
    for row in selected:
        tip = row["tip_id"]
        source = Path(row["source_fasta_path"])
        derivative = Path(row["derivative_fasta_path"])
        source_hash = row["source_fasta_sha256"]
        if sha256(source) != source_hash:
            raise SystemExit(f"{tip}: {source} no longer matches the hash frozen in {args.manifest}")
        cached = previous.get(tip)
        if (args.resume and cached and cached["source_fasta_sha256"] == source_hash
                and derivative.is_file() and sha256(derivative) == cached["derivative_sha256"]):
            mapped.append(cached)
            continue
        before_records, before_bases, ambiguous = fasta_stats(source)
        begin = time.monotonic()
        result = subprocess.run(
            [str(args.verticall), "repair", "-i", str(source), "-o", str(derivative)],
            text=True, capture_output=True)
        elapsed = time.monotonic() - begin
        (logs / f"{tip}.repair.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise SystemExit(f"{tip}: verticall repair exited {result.returncode}\n{result.stderr}")
        after_records, after_bases, after_ambiguous = fasta_stats(derivative)
        if after_ambiguous:
            raise SystemExit(f"{tip}: {after_ambiguous} ambiguous bases survived repair")
        mapped.append({
            "tip_id": tip,
            "source_fasta_path": str(source),
            "source_fasta_sha256": source_hash,
            "source_records": str(before_records),
            "source_bp": str(before_bases),
            "source_ambiguous_bp": str(ambiguous),
            "derivative_fasta_path": str(derivative),
            "derivative_sha256": sha256(derivative),
            "derivative_records": str(after_records),
            "derivative_bp": str(after_bases),
            "records_gained": str(after_records - before_records),
            "bases_removed": str(before_bases - after_bases),
            "derivative_ambiguous_bp": "0",
            "sequence_preserved": str(after_bases == before_bases and after_records == before_records),
            "file_identical_to_source": str(sha256(derivative) == source_hash),
            "repair_seconds": f"{elapsed:.3f}",
        })
    total_elapsed = time.monotonic() - started

    fields = ["tip_id", "source_fasta_path", "source_fasta_sha256", "source_records", "source_bp",
              "source_ambiguous_bp", "derivative_fasta_path", "derivative_sha256",
              "derivative_records", "derivative_bp", "derivative_ambiguous_bp", "records_gained",
              "bases_removed", "sequence_preserved", "file_identical_to_source", "repair_seconds"]
    with record_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(mapped)

    changed = [row for row in mapped if row["sequence_preserved"] == "False"]
    record = {
        "verticall_version": version.stdout.strip() or version.stderr.strip(),
        "manifest": str(args.manifest),
        "manifest_sha256": sha256(args.manifest),
        "rows_in_pilot": len(rows),
        "derivatives_written": len(mapped),
        "complete": len(mapped) == len(rows),
        "derivative_dir": str(derivative_dir),
        "elapsed_seconds": round(total_elapsed, 3),
        "seconds_per_genome": round(total_elapsed / max(len(mapped), 1), 3),
        "sequence_preserved": len(mapped) - len(changed),
        "sequence_changed_by_repair": len(changed),
        "file_rewritten_but_sequence_preserved": sum(
            1 for row in mapped
            if row["sequence_preserved"] == "True" and row["file_identical_to_source"] == "False"),
        "total_bases_removed": sum(int(row["bases_removed"]) for row in mapped),
        "total_records_gained": sum(int(row["records_gained"]) for row in mapped),
        "genomes_with_ambiguous_bases": sum(1 for row in mapped if int(row["source_ambiguous_bp"])),
        "max_bases_removed": max((int(row["bases_removed"]) for row in mapped), default=0),
        "max_records_gained": max((int(row["records_gained"]) for row in mapped), default=0),
    }
    (args.output / "repair-record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
