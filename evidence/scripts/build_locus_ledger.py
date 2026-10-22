#!/usr/bin/env python3
"""Decide, and record, which loci carry a gene tree on a panel -- markers or core genes.

A gene tree from a locus with no informative site on the panel is a star, and concordance
factors computed over stars say nothing.  So each locus is included or excluded by a stated rule,
and the ledger keeps every locus with its reason, including the excluded ones.

``markers``: the BAC120 panel alignment and its marker partitions.  A marker carries a gene tree
if it has at least ``--min-informative`` parsimony-informative columns on the panel and residues in
at least ``--min-taxa-fraction`` of the panel.

``core``: Panaroo's per-gene core alignments (``aligned_gene_sequences/``) for the ids given.
Genes are taken in Panaroo's core-alignment order.  Every gene enters the concatenation, with gaps
for an absent isolate; a gene carries a gene tree if it also passes the informative-site and
taxon-fraction rules.  An isolate with several sequences in one gene is
handled by ``merge_fragments``: non-overlapping pieces of a broken gene are merged, overlapping
copies mask that isolate for that gene, and a gene is excluded only when more than
``--max-masked-fraction`` of isolates are masked.  With ``--positional-reference`` each gene is placed on
that complete genome's chromosome through Panaroo's presence/absence table and the reference's
pyrodigal GFF, which is what the positional test of line 6 reads.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

AA = set("ACDEFGHIKLMNPQRSTVWY")
NT = set("ACGT")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, list[str]]] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                out.append((line[1:].split()[0], []))
            elif out:
                out[-1][1].append(line.strip())
    return [(n, "".join(s)) for n, s in out]


def site_counts(rows: list[str], alphabet: set[str]) -> tuple[int, int]:
    variable = informative = 0
    if not rows:
        return 0, 0
    for col in range(len(rows[0])):
        counts: dict[str, int] = {}
        for row in rows:
            c = row[col].upper()
            if c in alphabet:
                counts[c] = counts.get(c, 0) + 1
        if len(counts) > 1:
            variable += 1
            if sum(1 for v in counts.values() if v >= 2) >= 2:
                informative += 1
    return variable, informative


def write_nexus(path: Path, charsets: list[tuple[str, int, int]]) -> None:
    lines = ["#nexus", "begin sets;"]
    lines += [f"  charset {name} = {start}-{end};" for name, start, end in charsets]
    lines.append("end;")
    path.write_text("\n".join(lines) + "\n")


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def decide(taxa_present: int, n_taxa: int, informative: int, args) -> str:
    if informative < args.min_informative:
        return f"excluded:informative_sites<{args.min_informative}"
    if taxa_present < args.min_taxa_fraction * n_taxa:
        return f"excluded:taxa_present<{args.min_taxa_fraction:g}_of_panel"
    return "included"


def run_markers(args) -> dict:
    records = read_fasta(args.alignment)
    ids = [n for n, _ in records]
    rows = [s for _, s in records]
    parts = [(m.group(1), int(m.group(2)), int(m.group(3))) for m in
             re.finditer(r"charset\s+(\S+)\s*=\s*(\d+)-(\d+)\s*;", args.partitions.read_text())]
    ledger, loci = [], []
    for name, start, end in parts:
        block = [r[start - 1:end] for r in rows]
        present = sum(1 for b in block if any(c.upper() in AA for c in b))
        variable, informative = site_counts(block, AA)
        decision = decide(present, len(ids), informative, args)
        ledger.append({"locus": name, "start": start, "end": end, "columns": end - start + 1,
                       "taxa_present": present, "variable": variable,
                       "parsimony_informative": informative, "decision": decision})
        if decision == "included":
            loci.append((name, start, end))
    write_nexus(args.output / "loci.nex", loci)
    return {"mode": "markers", "alignment": str(args.alignment), "taxa": len(ids),
            "loci_total": len(parts), "loci_included": len(loci), "ledger": ledger}


def merge_fragments(copies: list[str]) -> str | None:
    """One row from several aligned pieces of one gene in one isolate, or None if they conflict.

    Panaroo clusters the pieces of a gene broken by a frameshift or an internal stop into one
    gene; aligned, they occupy different columns.  Pieces that never both call a base in the same
    column are merged column by column.  Pieces that overlap are two copies, not two fragments,
    and the isolate is masked for that gene rather than guessed at.
    """
    width = len(copies[0])
    out = []
    for col in range(width):
        called = {c[col] for c in copies if c[col] not in "-Nn?"}
        called_positions = sum(1 for c in copies if c[col] not in "-Nn?")
        if called_positions > 1:
            return None
        out.append(called.pop() if called else "-")
    return "".join(out)


def core_gene_order(embl: Path) -> list[str]:
    return re.findall(r"/label=(\S+?)\.aln\b", embl.read_text())


def run_core(args) -> dict:
    ids = [l.split()[0] for l in args.ids.read_text().splitlines() if l.strip()]
    wanted = set(ids)
    genes = core_gene_order(args.embl)
    if not genes:
        raise SystemExit(f"{args.embl}: no core genes listed")
    concat = {t: [] for t in ids}
    ledger, loci = [], []
    position = 0
    for gene in genes:
        path = args.gene_dir / f"{gene}.aln.fas"
        if not path.is_file():
            raise SystemExit(f"missing per-gene alignment {path}")
        per_isolate: dict[str, list[str]] = {}
        for name, seq in read_fasta(path):
            isolate = name.split(";")[0]
            if isolate in wanted:
                per_isolate.setdefault(isolate, []).append(seq)
        width = len(next(iter(per_isolate.values()))[0]) if per_isolate else 0
        merged = masked = 0
        for isolate, copies in list(per_isolate.items()):
            if len(copies) == 1:
                continue
            combined = merge_fragments(copies)
            if combined is None:
                del per_isolate[isolate]
                masked += 1
            else:
                per_isolate[isolate] = [combined]
                merged += 1
        entry = {"locus": gene, "columns": width, "taxa_present": len(per_isolate),
                 "isolates_fragments_merged": merged, "isolates_masked_conflicting_copies": masked}
        if masked > args.max_masked_fraction * len(ids):
            entry.update({"start": "", "end": "", "variable": "", "parsimony_informative": "",
                          "decision": "excluded:conflicting_copies_in_many_isolates",
                          "in_concatenation": False})
            ledger.append(entry)
            continue
        if width == 0:
            entry.update({"start": "", "end": "", "variable": "", "parsimony_informative": "",
                          "decision": "excluded:absent_from_panel", "in_concatenation": False})
            ledger.append(entry)
            continue
        block = [per_isolate[t][0] if t in per_isolate else "-" * width for t in ids]
        variable, informative = site_counts([per_isolate[t][0] for t in per_isolate], NT)
        start, end = position + 1, position + width
        position = end
        for t, row in zip(ids, block):
            concat[t].append(row)
        decision = decide(len(per_isolate), len(ids), informative, args)
        entry.update({"start": start, "end": end, "variable": variable,
                      "parsimony_informative": informative, "decision": decision,
                      "in_concatenation": True})
        ledger.append(entry)
        if decision == "included":
            loci.append((gene, start, end))
    with (args.output / "alignment.fna").open("w") as handle:
        for t in ids:
            handle.write(f">{t}\n{''.join(concat[t])}\n")
    write_nexus(args.output / "loci.nex", loci)
    in_concat = [(e["locus"], e["start"], e["end"]) for e in ledger if e["in_concatenation"]]
    write_nexus(args.output / "concatenation-genes.nex", in_concat)
    record = {"mode": "core", "ids": str(args.ids), "taxa": len(ids),
              "core_genes_listed": len(genes), "genes_in_concatenation": len(in_concat),
              "loci_included": len(loci), "concatenated_columns": position, "ledger": ledger}
    if args.positional_reference:
        record["positions"] = positions(args, genes)
    return record


def positions(args, genes: list[str]) -> dict:
    ref = args.positional_reference
    with args.presence.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        col = header.index(ref)
        annot = {row[0]: row[col] for row in reader}
    coords = {}
    for line in args.reference_gff.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) == 9 and parts[2] == "CDS":
            gid = re.search(r"ID=([^;]+)", parts[8]).group(1)
            coords[gid] = (parts[0], int(parts[3]), int(parts[4]), parts[6])
    offsets, total, name = {}, 0, None
    lengths: dict[str, int] = {}
    order = []
    for line in args.reference_fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            order.append(name)
            lengths[name] = 0
        else:
            lengths[name] += len(line.strip())
    for n in order:
        offsets[n] = total
        total += lengths[n]
    rows = []
    for gene in genes:
        ids = [i for i in annot.get(gene, "").split(";") if i]
        if len(ids) != 1 or ids[0] not in coords:
            rows.append({"locus": gene, "reference_gene": ";".join(ids), "contig": "",
                         "start": "", "end": "", "strand": "", "genome_position": "",
                         "placed": False})
            continue
        contig, start, end, strand = coords[ids[0]]
        rows.append({"locus": gene, "reference_gene": ids[0], "contig": contig, "start": start,
                     "end": end, "strand": strand,
                     "genome_position": offsets[contig] + (start + end) // 2, "placed": True})
    write_tsv(args.output / "gene-positions.tsv", rows)
    return {"reference": ref, "reference_length": total,
            "placed": sum(1 for r in rows if r["placed"]), "genes": len(rows)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    for mode in ("markers", "core"):
        p = sub.add_parser(mode)
        p.add_argument("--output", type=Path, required=True)
        p.add_argument("--min-informative", type=int, default=1)
        p.add_argument("--min-taxa-fraction", type=float, default=0.8)
        if mode == "markers":
            p.add_argument("--alignment", type=Path, required=True)
            p.add_argument("--partitions", type=Path, required=True)
        else:
            p.add_argument("--ids", type=Path, required=True)
            p.add_argument("--max-masked-fraction", type=float, default=0.05)
            p.add_argument("--embl", type=Path, required=True)
            p.add_argument("--gene-dir", type=Path, required=True)
            p.add_argument("--presence", type=Path)
            p.add_argument("--positional-reference")
            p.add_argument("--reference-gff", type=Path)
            p.add_argument("--reference-fasta", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    record = run_markers(args) if args.mode == "markers" else run_core(args)
    ledger = record.pop("ledger")
    write_tsv(args.output / "ledger.tsv", ledger)
    record["rule"] = {"min_informative": args.min_informative,
                      "min_taxa_fraction": args.min_taxa_fraction}
    for name in ("loci.nex", "alignment.fna"):
        path = args.output / name
        if path.exists():
            record[f"{name}_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (args.output / "ledger-record.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
