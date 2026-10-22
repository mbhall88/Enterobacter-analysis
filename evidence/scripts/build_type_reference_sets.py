"""Build replicon-scoped FASTA sets for the Task 02 type anchors.

Two scopes are produced for every assembly, so that a type-versus-study comparison can be run at a
declared replicon scope instead of silently mixing whole assemblies with chromosome-only study
files:

``whole``
    The NCBI Datasets ``*_genomic.fna`` exactly as downloaded.
``chromosome``
    Every record whose sequence report entry has ``assignedMoleculeLocationType == "Chromosome"``.

The ``role`` field matters and is never collapsed into the location type. In a complete genome the
chromosome record has ``role == "assembled-molecule"``. In a draft assembly the same location type
is carried by ``unplaced-scaffold`` contigs: NCBI is saying "not plasmid", not "this is one
assembled chromosome". Both are written, and the per-assembly record reports the role composition
and contig count so that fragmentation is visible to the reader rather than hidden inside a scope
name.

Output FASTA headers are reduced to the bare GenBank accession, which is the convention the
archived ``*_parsed.fna`` study chromosome files already use.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

CHROMOSOME = "Chromosome"


@dataclass
class RepliconSummary:
    """What the sequence report says about one assembly's records."""

    accession: str
    chromosome_ids: list[str] = field(default_factory=list)
    chromosome_bp: int = 0
    chromosome_roles: dict[str, int] = field(default_factory=dict)
    other_bp: int = 0
    other_location_types: dict[str, int] = field(default_factory=dict)
    total_records: int = 0

    @property
    def scope_status(self) -> str:
        if not self.chromosome_ids:
            return "no_chromosome_assigned_record"
        if self.chromosome_roles == {"assembled-molecule": len(self.chromosome_ids)}:
            return "assembled_chromosome"
        if "assembled-molecule" in self.chromosome_roles:
            return "mixed_assembled_and_unplaced"
        return "chromosome_assigned_draft_contigs"


def read_sequence_report(path: Path) -> RepliconSummary:
    accession = ""
    summary = RepliconSummary(accession="")
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        accession = accession or record["assemblyAccession"]
        summary.total_records += 1
        length = int(record["length"])
        location = record.get("assignedMoleculeLocationType") or "unassigned"
        if location == CHROMOSOME:
            summary.chromosome_ids.append(record["genbankAccession"])
            summary.chromosome_bp += length
            role = record.get("role") or "unknown"
            summary.chromosome_roles[role] = summary.chromosome_roles.get(role, 0) + 1
        else:
            summary.other_bp += length
            summary.other_location_types[location] = (
                summary.other_location_types.get(location, 0) + 1
            )
    summary.accession = accession
    return summary


def read_fasta(path: Path) -> list[tuple[str, list[str]]]:
    records: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line[1:].split()[0], []))
        elif records:
            records[-1][1].append(line.strip())
    return records


def write_chromosome_fasta(source: Path, wanted: list[str], destination: Path) -> int:
    """Write the wanted records in sequence-report order; fail loudly on a missing record."""
    available = dict(read_fasta(source))
    missing = [identifier for identifier in wanted if identifier not in available]
    if missing:
        raise SystemExit(f"{source}: sequence report lists records absent from the FASTA: {missing}")
    written = 0
    with destination.open("w") as handle:
        for identifier in wanted:
            handle.write(f">{identifier}\n")
            for chunk in available[identifier]:
                handle.write(f"{chunk}\n")
            written += 1
    return written


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-dir", type=Path, required=True,
                        help="ncbi_dataset/data directory from a datasets download")
    parser.add_argument("--accessions", type=Path, required=True,
                        help="one accession.version per line; every one must be present")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    accessions = [line.strip() for line in args.accessions.read_text().splitlines() if line.strip()]
    chromosome_dir = args.output / "chromosome"
    chromosome_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for accession in accessions:
        directory = args.datasets_dir / accession
        if not directory.is_dir():
            raise SystemExit(f"{accession}: no downloaded directory at {directory}")
        fasta_candidates = sorted(directory.glob("*_genomic.fna"))
        if len(fasta_candidates) != 1:
            raise SystemExit(f"{accession}: expected one *_genomic.fna, found {fasta_candidates}")
        whole = fasta_candidates[0]
        summary = read_sequence_report(directory / "sequence_report.jsonl")
        if summary.accession != accession:
            raise SystemExit(
                f"{accession}: sequence report reports {summary.accession}; refusing to continue")

        chromosome_path = ""
        chromosome_file_sha = ""
        if summary.chromosome_ids:
            destination = chromosome_dir / f"{accession}_chromosome_scope.fna"
            write_chromosome_fasta(whole, summary.chromosome_ids, destination)
            chromosome_path = str(destination)
            chromosome_file_sha = sha256(destination)

        rows.append({
            "accession_version": accession,
            "whole_path": str(whole),
            "whole_sha256": sha256(whole),
            "whole_records": str(summary.total_records),
            "whole_bp": str(summary.chromosome_bp + summary.other_bp),
            "chromosome_scope_path": chromosome_path,
            "chromosome_scope_sha256": chromosome_file_sha,
            "chromosome_records": str(len(summary.chromosome_ids)),
            "chromosome_bp": str(summary.chromosome_bp),
            "chromosome_roles": ";".join(
                f"{role}={count}" for role, count in sorted(summary.chromosome_roles.items())),
            "non_chromosome_records": str(summary.total_records - len(summary.chromosome_ids)),
            "non_chromosome_bp": str(summary.other_bp),
            "non_chromosome_location_types": ";".join(
                f"{name}={count}" for name, count in sorted(summary.other_location_types.items())),
            "chromosome_scope_status": summary.scope_status,
        })

    table = args.output / "type-reference-sequence-scopes.tsv"
    header = list(rows[0])
    with table.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row[column] for column in header) + "\n")
    print(f"wrote {len(rows)} assemblies to {table}")


if __name__ == "__main__":
    main()
