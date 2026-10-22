"""Extend the paper Verticall table with the 100 expert genomes (validation-design.md section 7.2).

stage-map  writes the derivative map and ID list for scripts/run_verticall_cohort.sh --stage stage:
           the 1,089 cohort tips from Task 11's repair map (which carries derivative hashes), checked
           against the 1,089-tip freeze, plus the new genomes' repair map, whose source hashes must
           equal the query manifest's (the same FASTA bytes skani and analyze read).
skip-list  writes the paper table's directed pairs as a two-column TSV for
           `verticall pairwise --existing_tsv`, so no part re-reads the 16 GB table.
merge      concatenates the paper table and the new parts into one table under one header, refusing
           a part count other than expected, a repeated paper pair, a row without a new genome or a
           stray header. Write it to RUN_DIR/<label>-pairwise.tsv, the table --stage downstream reads.

Written and tested before download; hashed at freeze.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
T11 = ROOT / "results/full-cohort/2026-09-22/task11-opus-full"
REPAIR_MAP = T11 / "repair/source-derivative-map.tsv"
FREEZE_MAP = T11 / "freeze-1089/source-derivative-map.tsv"
PAPER_TABLE = T11 / "verticall-1089/full-1089-pairwise.tsv"
PAPER_TABLE_SHA256 = T11 / "verticall-1089/full-1089-pairwise.sha256"
HEADER_START = "assembly_a\tassembly_b"
csv.field_size_limit(1_000_000_000)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def expected_sha(sha_file):
    return Path(sha_file).read_text().split()[0]


def stage_map(repair_map, freeze_map, new_map, query_manifest, output_map, output_ids):
    """Cohort rows from the repair map, restricted to and checked against the freeze; new rows appended."""
    freeze = {r["tip_id"]: r for r in read_tsv(freeze_map)}
    repair = {r["tip_id"]: r for r in read_tsv(repair_map)}
    new = read_tsv(new_map)
    query = {r["accession"]: r["query_sha256"] for r in read_tsv(query_manifest)}
    if {r["tip_id"] for r in new} != set(query):
        raise ValueError("new genomes differ from the query manifest")
    if any(r["source_fasta_sha256"] != query[r["tip_id"]] for r in new):
        raise ValueError("a new genome's repair source is not the query manifest's FASTA")
    if set(freeze) - set(repair):
        raise ValueError("freeze tips missing from the repair map")
    rows = []
    for tip in sorted(freeze):
        r = repair[tip]
        if r["derivative_fasta_path"] != freeze[tip]["derivative_fasta_path"]:
            raise ValueError(f"derivative path differs between repair map and freeze for {tip}")
        rows.append(dict(tip_id=tip, derivative_fasta_path=r["derivative_fasta_path"], derivative_sha256=r["derivative_sha256"],
                         origin="paper_cohort"))
    clash = {r["tip_id"] for r in new} & set(freeze)
    if clash:
        raise ValueError(f"new genomes clash with cohort tips: {sorted(clash)}")
    if len({r["tip_id"] for r in new}) != len(new) or any(not r["derivative_sha256"] for r in new):
        raise ValueError("new map has duplicate tips or missing hashes")
    rows += [dict(tip_id=r["tip_id"], derivative_fasta_path=r["derivative_fasta_path"], derivative_sha256=r["derivative_sha256"],
                  origin="expert_sample") for r in sorted(new, key=lambda x: x["tip_id"])]
    with Path(output_map).open("x", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    Path(output_ids).open("x").write("".join(r["tip_id"] + "\n" for r in rows))
    return dict(cohort=len(freeze), new=len(new), total=len(rows))


def skip_list(table, table_sha_file, output, expected_pairs=None):
    if sha256(table) != expected_sha(table_sha_file):
        raise ValueError("paper table hash mismatch")
    pairs = set()
    with Path(table).open() as handle:
        header = handle.readline()
        if not header.startswith(HEADER_START):
            raise ValueError("paper table does not start with the Verticall header")
        for line in handle:
            a, b = line.split("\t", 2)[:2]
            pairs.add((a, b))
    if expected_pairs is not None and len(pairs) != expected_pairs:
        raise ValueError(f"{len(pairs)} directed pairs, expected {expected_pairs}")
    with Path(output).open("x") as handle:
        handle.write(HEADER_START + "\n")
        for a, b in sorted(pairs):
            handle.write(f"{a}\t{b}\n")
    return dict(directed_pairs=len(pairs), skip_list_sha256=sha256(output), paper_table_sha256=expected_sha(table_sha_file))


def merge(table, table_sha_file, parts, new_ids, output, expected_parts=None):
    """Paper table, then every part in order; part 1 alone carries a header, equal to the paper header."""
    if expected_parts is not None and len(parts) != expected_parts:
        raise ValueError(f"{len(parts)} parts given, expected {expected_parts}")
    if sha256(table) != expected_sha(table_sha_file):
        raise ValueError("paper table hash mismatch")
    new = {line.strip() for line in Path(new_ids).read_text().splitlines() if line.strip()}
    paper_pairs, paper_rows = set(), 0
    ledger = []
    with Path(output).open("x") as out, Path(table).open() as handle:
        header = handle.readline()
        if not header.startswith(HEADER_START):
            raise ValueError("paper table does not start with the Verticall header")
        out.write(header)
        for line in handle:
            a, b = line.split("\t", 2)[:2]
            if a in new or b in new:
                raise ValueError("the paper table already names a new genome")
            paper_pairs.add((a, b))
            paper_rows += 1
            out.write(line)
        new_rows = 0
        for i, part in enumerate(parts, start=1):
            rows = 0
            with Path(part).open() as handle:
                first = handle.readline()
                if i == 1:
                    if first != header:
                        raise ValueError("part 1 header differs from the paper table header")
                elif first.startswith(HEADER_START):
                    raise ValueError(f"part {i} carries a header; only part 1 may")
                lines = handle if i == 1 else _chain(first, handle)
                for line in lines:
                    if line.startswith(HEADER_START):
                        raise ValueError(f"part {i} carries a second header")
                    a, b = line.split("\t", 2)[:2]
                    if (a, b) in paper_pairs:
                        raise ValueError(f"part {i} repeats a paper pair {a} {b}")
                    if a not in new and b not in new:
                        raise ValueError(f"part {i} row names no new genome: {a} {b}")
                    out.write(line)
                    rows += 1
            ledger.append(dict(part=i, file=str(part), data_rows=rows, sha256=sha256(part)))
            new_rows += rows
    with Path(output).open() as handle:
        total = sum(1 for _ in handle) - 1
    if total != paper_rows + new_rows:
        raise ValueError("merged row count does not equal paper plus parts")
    return dict(paper_rows=paper_rows, new_rows=new_rows, merged_rows=total, merged_sha256=sha256(output), parts=ledger)


def _chain(first, handle):
    if first:
        yield first
    yield from handle


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("stage-map")
    s.add_argument("--new-map", type=Path, required=True, help="repair map of the new genomes (repair_pilot_inputs.py output)")
    s.add_argument("--query-manifest", type=Path, required=True,
                   help="accession, query_sha256 (the label_transfer.py skani manifest)")
    s.add_argument("--output-map", type=Path, required=True)
    s.add_argument("--output-ids", type=Path, required=True)
    k = sub.add_parser("skip-list")
    k.add_argument("--output", type=Path, required=True)
    k.add_argument("--expected-pairs", type=int, default=1_184_832)
    m = sub.add_parser("merge")
    m.add_argument("--parts", type=Path, nargs="+", required=True, help="part_1_of_N ... part_N_of_N, in order")
    m.add_argument("--expected-parts", type=int, default=100)
    m.add_argument("--new-ids", type=Path, required=True)
    m.add_argument("--output", type=Path, required=True)
    m.add_argument("--record", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "stage-map":
        result = stage_map(REPAIR_MAP, FREEZE_MAP, args.new_map, args.query_manifest, args.output_map, args.output_ids)
    elif args.command == "skip-list":
        result = skip_list(PAPER_TABLE, PAPER_TABLE_SHA256, args.output, args.expected_pairs)
    else:
        result = merge(PAPER_TABLE, PAPER_TABLE_SHA256, args.parts, args.new_ids, args.output, args.expected_parts)
        Path(args.record).open("x").write(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "parts"}, indent=2))


if __name__ == "__main__":
    main()
