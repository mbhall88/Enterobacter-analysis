#!/usr/bin/env python3
"""Measure how much of a panel an alignment can tell apart, and which type anchors it cannot.

SET-16 made the untrimmed, occupancy-filtered alignment primary on *resolution* grounds, not on
support: at pilot scale GTDB's canonical mask left 183 distinct sequences among 388 genomes and made
two pairs of differently named nomenclatural type anchors byte-identical, so it could not address
the study's question even in principle. That ratio is a measurement, not a property of the mask, so
it has to be re-measured whenever the panel changes — which is what this script is for.

What it reports, and why each number is separate
------------------------------------------------
``distinct_sequences``
    How many different sequences the alignment actually contains. Two genomes with the same
    sequence are one tip's worth of evidence however they are labelled.
``not_uniquely_resolved``
    Genomes sharing their sequence with at least one other genome, as a count and a fraction. This
    is the number to quote for "how much of the panel the alignment cannot separate" — it is *not*
    the same as ``distinct_sequences`` and the two are easy to conflate.
``redundant``
    ``genomes - distinct_sequences``: how many genomes carry no sequence of their own. This is the
    quantity Task 07 published as "not uniquely resolved" (52.8% masked, 21.4% untrimmed at pilot
    scale), and it is always the smaller of the two — a group of *k* identical genomes contributes
    *k* to ``not_uniquely_resolved`` but only *k − 1* to ``redundant``. Both are reported under
    their own names so a full-scale figure is never quietly compared with a pilot figure computed
    the other way.
``largest_identical_group``
    The worst case, which a percentage hides.
``identical_type_anchor_pairs``
    Pairs of **differently named** nomenclatural type anchors that are byte-identical in this
    alignment. This is the publication-critical one: if two type anchors are identical, no tree
    built on this alignment can separate the names they anchor, whatever its support values say.

Ambiguity handling is explicit: sequences are compared after upper-casing, and gap and ambiguity
characters are compared as themselves rather than being treated as wildcards. Two sequences that
differ only in where their gaps fall are therefore *different*, which is the conservative direction
— it can only make the alignment look more resolving than it is, never less.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: dict[str, list[str]] = {}
    name = None
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise SystemExit(f"{path}: duplicate sequence name {name!r}")
                records[name] = []
            elif name is not None:
                records[name].append(line.strip())
    if not records:
        raise SystemExit(f"{path}: no sequences")
    return {name: "".join(parts).upper() for name, parts in records.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", type=Path, action="append", required=True, default=[],
                        help="LABEL=PATH; repeat for each alignment form to compare")
    parser.add_argument("--panel", type=Path,
                        help="frozen panel manifest, for formal type-anchor names")
    parser.add_argument("--id-map", type=Path,
                        help="TSV with gtdbtk_genome_id and tip_id, if the alignment is prefixed")
    parser.add_argument("--expected-tips", type=Path,
                        help="one tip id per line; every alignment must carry exactly this set")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    id_map: dict[str, str] = {}
    if args.id_map:
        with args.id_map.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                id_map[row["gtdbtk_genome_id"]] = row["tip_id"]

    formal: dict[str, list[str]] = {}
    if args.panel:
        with args.panel.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("formal_names"):
                    formal[row["tip_id"]] = [name for name in row["formal_names"].split(";") if name]

    expected = None
    if args.expected_tips:
        expected = sorted({line.strip() for line in args.expected_tips.read_text().splitlines()
                           if line.strip()})

    report = {}
    for spec in args.alignment:
        label, _, path_text = str(spec).partition("=")
        if not path_text:
            raise SystemExit(f"--alignment needs LABEL=PATH, got {spec!r}")
        path = Path(path_text)
        records = read_fasta(path)
        renamed = {id_map.get(name, name): sequence for name, sequence in records.items()}
        if len(renamed) != len(records):
            raise SystemExit(f"{path}: the id map collapses two sequences onto one tip id")
        if expected is not None and sorted(renamed) != expected:
            missing = sorted(set(expected) - set(renamed))
            extra = sorted(set(renamed) - set(expected))
            raise SystemExit(f"{path}: tip set differs from {args.expected_tips}: "
                             f"{len(missing)} missing {missing[:5]}, "
                             f"{len(extra)} unexpected {extra[:5]}")

        widths = {len(sequence) for sequence in renamed.values()}
        if len(widths) != 1:
            raise SystemExit(f"{path}: sequences are not all the same length: {sorted(widths)[:5]}")

        groups: dict[str, list[str]] = defaultdict(list)
        for tip, sequence in renamed.items():
            groups[sequence].append(tip)
        sizes = Counter(len(members) for members in groups.values())
        shared = sum(len(members) for members in groups.values() if len(members) > 1)
        multi = sorted((sorted(members) for members in groups.values() if len(members) > 1),
                       key=lambda members: (-len(members), members[0]))

        anchor_collisions = []
        for members in multi:
            named = [(tip, formal[tip]) for tip in members if tip in formal]
            names = {name for _, entries in named for name in entries}
            if len(named) > 1 and len(names) > 1:
                anchor_collisions.append({
                    "tips": [tip for tip, _ in named],
                    "formal_names": sorted(names),
                    "group_size": len(members),
                })

        report[label] = {
            "path": str(path),
            "genomes": len(renamed),
            "alignment_columns": widths.pop(),
            "distinct_sequences": len(groups),
            "not_uniquely_resolved_genomes": shared,
            "not_uniquely_resolved_fraction": round(shared / len(renamed), 6),
            "redundant_genomes": len(renamed) - len(groups),
            "redundant_fraction": round((len(renamed) - len(groups)) / len(renamed), 6),
            "largest_identical_group": max(sizes) if sizes else 0,
            "identical_groups": len(multi),
            "group_size_distribution": dict(sorted(sizes.items())),
            "identical_type_anchor_pairs": anchor_collisions,
            "identical_type_anchor_pair_count": len(anchor_collisions),
            "largest_identical_groups": [
                {"size": len(members), "tips": members[:12],
                 "truncated": len(members) > 12} for members in multi[:5]],
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    brief = {label: {key: value[key] for key in
                     ("genomes", "alignment_columns", "distinct_sequences",
                      "not_uniquely_resolved_genomes", "not_uniquely_resolved_fraction",
                      "redundant_genomes", "redundant_fraction",
                      "largest_identical_group", "identical_type_anchor_pair_count")}
             for label, value in report.items()}
    print(json.dumps(brief, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
