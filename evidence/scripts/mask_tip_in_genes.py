#!/usr/bin/env python3
"""Mask one genome's sequence in named genes of a concatenation, for a falsification test.

C3: the *E. hormaechei* subsp. *hoffmannii* type genome has a terminal branch ~19x the trio
clade's median and a Verticall r/m ~7.6x that between the other two trio types, the signature of
imported divergent segments.  If its basal placement in the concatenated ML tree comes from those
segments, masking its sequence in the core genes that overlap its Verticall horizontal regions
should remove the preference; if the preference survives, it is not the imports.  This writes the
concatenation with that genome's columns set to gaps in exactly those genes, and nothing else
changed, plus a record of what was masked.  The gene list comes from Verticall's regions, fixed
before the masked data are analysed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True, help="concatenation-genes.nex")
    parser.add_argument("--mask-genes", type=Path, required=True, help="one locus per line")
    parser.add_argument("--tip", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")
    ranges = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in
              re.finditer(r"charset\s+(\S+)\s*=\s*(\d+)-(\d+)\s*;", args.genes.read_text())}
    wanted = [l.strip() for l in args.mask_genes.read_text().splitlines() if l.strip()]
    missing = [g for g in wanted if g not in ranges]
    names, seqs, name = [], {}, None
    for line in args.alignment.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            names.append(name)
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    if args.tip not in seqs:
        raise SystemExit(f"{args.tip} not in the alignment")
    row = list("".join(seqs[args.tip]))
    masked_columns = 0
    for gene in wanted:
        if gene in ranges:
            start, end = ranges[gene]
            for i in range(start - 1, end):
                if row[i] != "-":
                    masked_columns += 1
                row[i] = "-"
    with args.output.open("w") as handle:
        for n in names:
            handle.write(f">{n}\n{''.join(row) if n == args.tip else ''.join(seqs[n])}\n")
    record = {"alignment": str(args.alignment), "tip": args.tip, "genes_requested": len(wanted),
              "genes_masked": len(wanted) - len(missing), "genes_not_in_concatenation": missing,
              "called_columns_masked": masked_columns,
              "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest()}
    args.output.with_suffix(".json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k != "genes_not_in_concatenation"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
