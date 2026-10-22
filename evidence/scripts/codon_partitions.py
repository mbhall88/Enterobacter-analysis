#!/usr/bin/env python3
"""Write codon-position partition schemes for a core-gene concatenation.

Panaroo's core alignments are codon alignments (protein-aligned, back-translated), so within each
gene the columns keep their reading frame.  A single substitution model over a genus-level
concatenation treats saturated third positions like conserved second ones, which is the kind of
misspecification that pulls long branches together.  So the core-genome topology tests use a
partition per codon position, built gene by gene (``start-end\\3``) so the phase survives the
concatenation, and a sensitivity that drops third positions entirely.

Schemes:

``codon3``  first, second and third positions, each with its own model; genes whose aligned length
            is not a multiple of three (Panaroo re-aligned them as DNA) form a fourth partition.
``pos12``   first and second positions only, of the in-frame genes.  Columns in no charset are not
            analysed by IQ-TREE, so third positions and the out-of-frame genes are dropped.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def gene_ranges(nexus: Path) -> list[tuple[str, int, int]]:
    return [(m.group(1), int(m.group(2)), int(m.group(3))) for m in
            re.finditer(r"charset\s+(\S+)\s*=\s*(\d+)-(\d+)\s*;", nexus.read_text())]


def build(genes: list[tuple[str, int, int]], scheme: str, model: str) -> tuple[str, dict]:
    in_frame = [(g, a, b) for g, a, b in genes if (b - a + 1) % 3 == 0]
    off_frame = [(g, a, b) for g, a, b in genes if (b - a + 1) % 3]
    sets: dict[str, list[str]] = {"pos1": [], "pos2": [], "pos3": []}
    for _, a, b in in_frame:
        sets["pos1"].append(f"{a}-{b}\\3")
        sets["pos2"].append(f"{a + 1}-{b}\\3")
        sets["pos3"].append(f"{a + 2}-{b}\\3")
    if scheme == "pos12":
        sets.pop("pos3")
    elif off_frame:
        sets["noncodon"] = [f"{a}-{b}" for _, a, b in off_frame]
    lines = ["#nexus", "begin sets;"]
    for name, ranges in sets.items():
        lines.append(f"  charset {name} = {' '.join(ranges)};")
    lines.append("  charpartition scheme = " +
                 ", ".join(f"{model}: {name}" for name in sets) + ";")
    lines.append("end;")
    info = {"scheme": scheme, "model_per_partition": model, "genes": len(genes),
            "genes_in_frame": len(in_frame), "genes_out_of_frame": len(off_frame),
            "out_of_frame_genes": [g for g, _, _ in off_frame],
            "partitions": list(sets),
            "columns_per_partition": {
                "pos1": sum((b - a + 1) // 3 for _, a, b in in_frame),
                "pos2": sum((b - a + 1) // 3 for _, a, b in in_frame),
                **({"pos3": sum((b - a + 1) // 3 for _, a, b in in_frame)} if "pos3" in sets else {}),
                **({"noncodon": sum(b - a + 1 for _, a, b in off_frame)} if "noncodon" in sets else {})}}
    return "\n".join(lines) + "\n", info


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", type=Path, required=True, help="concatenation-genes.nex")
    parser.add_argument("--scheme", choices=["codon3", "pos12"], required=True)
    parser.add_argument("--model", default="GTR+F+R4")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")
    text, info = build(gene_ranges(args.genes), args.scheme, args.model)
    args.output.write_text(text)
    args.output.with_suffix(".json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps({k: v for k, v in info.items() if k != "out_of_frame_genes"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
