"""Strict assembly structure, without inferring biological quality or taxonomy.

Lengths describe FASTA records (which may be scaffolds), including ambiguity
characters. The file hash covers supplied bytes, including gzip compression.
Completeness and contamination require separately validated tool evidence.
"""
from __future__ import annotations

from collections import Counter
import gzip
from pathlib import Path
import re

from .evidence import sha256

SCHEMA_VERSION = "1.0.0-dev-assembly-quality"
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")


def assess_assembly(path: Path) -> dict:
    """Return structural metrics for a plain or ``.gz`` DNA FASTA.

    Reject empty files/records, duplicate or empty IDs, sequence before a
    header, and characters outside the IUPAC DNA alphabet. The first
    decompressed byte must be ``>``. Empty sequence lines are ignored;
    sequence whitespace other than LF/CRLF line endings is invalid. Bare CR
    line endings are rejected because the unchanged ANI input uses LF records. Lowercase DNA
    is accepted. No records are filtered or rewritten.
    A syntactically valid assembly always has unknown biological quality.
    """
    path = Path(path).resolve(strict=True)
    initial_hash = sha256(path)
    opener = gzip.open if path.suffix == ".gz" else open
    identifiers: set[str] = set()
    lengths: list[int] = []
    records: list[dict] = []
    counts: Counter[str] = Counter()
    identifier: str | None = None
    length = 0
    record_acgt = longest_run = trailing_run = 0
    try:
        with opener(path, "rt", encoding="ascii", newline="") as handle:
            for number, raw in enumerate(handle, 1):
                if raw.endswith("\r"):
                    raise ValueError("bare CR line ending is not supported; use LF or CRLF")
                if number == 1 and not raw.startswith(">"):
                    raise ValueError("assembly FASTA must start with a header at its first byte")
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if identifier is not None:
                        if not length:
                            raise ValueError(f"empty FASTA record: {identifier}")
                        lengths.append(length)
                        records.append({"record_id": identifier, "length_bp": length,
                                        "acgt_bp": record_acgt,
                                        "max_unambiguous_acgt_run_bp": longest_run})
                    fields = line[1:].split()
                    if not fields or fields[0] in identifiers:
                        raise ValueError(f"empty or duplicate FASTA ID at line {number}")
                    identifier = fields[0]
                    identifiers.add(identifier)
                    length = 0
                    record_acgt = longest_run = trailing_run = 0
                else:
                    sequence = line.upper()
                    if identifier is None or not set(sequence) <= DNA_ALPHABET:
                        raise ValueError(f"invalid nucleotide assembly FASTA at line {number}")
                    counts.update(sequence)
                    length += len(sequence)
                    runs = re.split("[^ACGT]+", sequence)
                    record_acgt += sum(map(len, runs))
                    longest_run = max(longest_run, trailing_run + len(runs[0]),
                                      max(map(len, runs)))
                    trailing_run = trailing_run + len(sequence) if len(runs) == 1 else len(runs[-1])
    except (UnicodeError, gzip.BadGzipFile, EOFError) as exc:
        raise ValueError("assembly is not a readable ASCII nucleotide FASTA") from exc
    if identifier is None or not length:
        raise ValueError("assembly FASTA is empty or ends with an empty record")
    lengths.append(length)
    records.append({"record_id": identifier, "length_bp": length, "acgt_bp": record_acgt,
                    "max_unambiguous_acgt_run_bp": longest_run})
    if sha256(path) != initial_hash:
        raise ValueError("assembly changed while its structure was being assessed")
    total = sum(lengths)
    cumulative = 0
    n50 = l50 = 0
    for l50, size in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += size
        if cumulative * 2 >= total:
            n50 = size
            break
    acgt = sum(counts[base] for base in "ACGT")
    gc = counts["G"] + counts["C"]
    n = counts["N"]
    ambiguous = total - acgt
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "structurally_valid",
        "query_path": str(path),
        "query_sha256": initial_hash,
        "structural_metrics": {
            "record_count": len(lengths), "total_bp": total,
            "min_record_bp": min(lengths), "max_record_bp": max(lengths),
            "record_n50_bp": n50, "record_l50": l50,
            "max_unambiguous_acgt_run_bp": max(r["max_unambiguous_acgt_run_bp"] for r in records),
            "gc_bp": gc, "acgt_bp": acgt, "n_bp": n,
            "other_ambiguity_bp": ambiguous - n, "ambiguous_bp": ambiguous,
            "gc_percent_of_total_bp": 100 * gc / total,
            "gc_percent_of_acgt_bp": 100 * gc / acgt if acgt else None,
            "n_percent_of_total_bp": 100 * n / total,
            "ambiguous_percent_of_total_bp": 100 * ambiguous / total,
        },
        "record_metrics": records,
        "biological_quality": {
            "status": "unknown", "completeness_percent": None,
            "contamination_percent": None, "evidence": [],
        },
    }
