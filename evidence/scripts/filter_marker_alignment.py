#!/usr/bin/env python3
"""Occupancy-filter a GTDB-Tk concatenated marker alignment and remap its partitions.

This generalizes the eight-genome occupancy analysis in ``scripts/summarize_smoke.py`` into a
reusable interface for Task 07.

What it does
------------
The GTDB-Tk ``--skip_trimming`` alignment is a concatenation of the BAC120 marker alignments laid
out in ``gtdbtk.bac120.marker_info.tsv`` order, each marker occupying exactly ``Length (bp)``
columns.  That invariant is *checked*, not assumed: if the marker lengths do not sum to the
alignment width the run fails rather than producing a silently wrong partition file.

Columns are kept or dropped on **occupancy alone** -- the fraction of tips carrying a standard
amino acid.  Invariant (constant) columns are retained: no variability criterion is applied
anywhere, and no SNP-only ascertainment correction is appropriate for the output.  Tips are never
dropped; low-occupancy tips are reported so they can be judged before inference, because the
downstream completion criterion is that the tip set still matches the frozen panel exactly.

Because every marker occupies a contiguous block of the original alignment and filtering preserves
column order, each marker's retained columns remain contiguous in the filtered alignment.  A marker
that retains no columns is removed from the partition file entirely and recorded as removed.  The
partitions are verified to cover the retained columns exactly: no gap, no overlap, no column
outside a partition.

Outputs
-------
``alignment.faa``                filtered FASTA, tips in input order
``column-map.tsv``               every original column: marker, occupancy, kept/dropped, new index
``partitions.nex``               remapped marker partitions (Nexus sets block, for IQ-TREE ``-p``)
``per-tip-occupancy.tsv``        per tip, standard-AA occupancy before and after filtering
``per-marker-occupancy.tsv``     per marker, columns kept, tip occupancy, removal decision
``filter-summary.json``          site counts before/after, decisions, parameters, validation
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import sys
from pathlib import Path

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: Path) -> dict[str, str]:
    """Read an aligned FASTA.  Rejects duplicate names and ragged alignments."""
    seqs: dict[str, list[str]] = {}
    name = None
    with open_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0] if line[1:].split() else line[1:]
                if name in seqs:
                    raise ValueError(f"duplicate sequence name in {path}: {name}")
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError(f"sequence data before first header in {path}")
                seqs[name].append(line.strip())
    if not seqs:
        raise ValueError(f"no sequences in {path}")
    joined = {k: "".join(v) for k, v in seqs.items()}
    widths = {len(s) for s in joined.values()}
    if len(widths) != 1:
        raise ValueError(f"ragged alignment in {path}: widths {sorted(widths)}")
    return joined


def read_marker_info(path: Path) -> list[tuple[str, str, int]]:
    """Return [(marker_id, marker_name, length), ...] in file order."""
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Marker ID", "Name", "Length (bp)"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        markers = []
        seen = set()
        for row in reader:
            marker_id = row["Marker ID"].strip()
            if not marker_id:
                continue
            if marker_id in seen:
                raise ValueError(f"duplicate marker id in {path}: {marker_id}")
            seen.add(marker_id)
            markers.append((marker_id, row["Name"].strip(), int(row["Length (bp)"])))
    if not markers:
        raise ValueError(f"no markers in {path}")
    return markers


def column_stats(columns) -> dict:
    """Variable / parsimony-informative / constant / empty counts over an iterable of columns.

    A column is *variable* when it shows more than one distinct standard amino acid, and
    *parsimony-informative* when at least two distinct standard amino acids each appear at least
    twice.  Columns with no standard amino acid at all are counted separately and are neither
    variable nor constant.
    """
    variable = informative = constant = empty = 0
    total = 0
    for column in columns:
        total += 1
        counts = collections.Counter(c for c in column if c in STANDARD_AA)
        if not counts:
            empty += 1
            continue
        if len(counts) > 1:
            variable += 1
            if sum(n >= 2 for n in counts.values()) >= 2:
                informative += 1
        else:
            constant += 1
    return {
        "columns": total,
        "variable_columns": variable,
        "parsimony_informative_columns": informative,
        "constant_columns": constant,
        "columns_without_standard_amino_acids": empty,
    }


def sanitize_charset_name(marker_id: str) -> str:
    """Nexus charset names must not carry characters IQ-TREE's parser treats as syntax."""
    out = "".join(c if (c.isalnum() or c == "_") else "_" for c in marker_id)
    if out and out[0].isdigit():
        out = "m_" + out
    return out


def filter_alignment(
    seqs: dict[str, str],
    markers: list[tuple[str, str, int]],
    min_column_occupancy: float,
    min_marker_tip_occupancy: float = 0.0,
) -> dict:
    """Core computation.  Returns everything the writers need; touches no files."""
    tips = list(seqs)
    n_tips = len(tips)
    width = len(seqs[tips[0]])

    total_marker_length = sum(length for _, _, length in markers)
    if total_marker_length != width:
        raise ValueError(
            "marker lengths do not sum to the alignment width -- this alignment and this "
            f"marker_info file do not describe the same concatenation: {total_marker_length} "
            f"marker columns vs {width} alignment columns. If the alignment was masked/trimmed, "
            "the per-marker boundaries no longer apply and partitions cannot be derived."
        )

    # Uppercase once; GTDB writes uppercase, but a hand-made input must not change the counts.
    upper = {t: seqs[t].upper() for t in tips}
    columns = [tuple(upper[t][i] for t in tips) for i in range(width)]

    occupancy = [sum(1 for c in col if c in STANDARD_AA) / n_tips for col in columns]

    # Assign each original column to its marker.
    col_marker: list[int] = []
    for idx, (_, _, length) in enumerate(markers):
        col_marker.extend([idx] * length)

    # Stage 1: whole-marker removal on per-marker tip occupancy (off by default).
    marker_tip_occupancy = []
    for idx, (_, _, length) in enumerate(markers):
        start = sum(m[2] for m in markers[:idx])
        block = columns[start:start + length]
        if length == 0:
            marker_tip_occupancy.append(0.0)
            continue
        present = sum(
            1 for t_i in range(n_tips)
            if any(block[c][t_i] in STANDARD_AA for c in range(length))
        )
        marker_tip_occupancy.append(present / n_tips)

    marker_removed_by_rule = [occ < min_marker_tip_occupancy for occ in marker_tip_occupancy]

    # Stage 2: per-column occupancy.
    keep = [
        (occupancy[i] >= min_column_occupancy) and not marker_removed_by_rule[col_marker[i]]
        for i in range(width)
    ]

    retained_indices = [i for i in range(width) if keep[i]]
    new_index = {}
    for new, orig in enumerate(retained_indices, start=1):
        new_index[orig] = new

    filtered = {t: "".join(upper[t][i] for i in retained_indices) for t in tips}

    # Per-marker accounting and contiguous remapped ranges.
    marker_rows = []
    offset = 0
    for idx, (marker_id, name, length) in enumerate(markers):
        block = list(range(offset, offset + length))
        offset += length
        kept_cols = [i for i in block if keep[i]]
        row = {
            "marker_id": marker_id,
            "marker_name": name,
            "charset_name": sanitize_charset_name(marker_id),
            "original_length": length,
            "original_start": block[0] + 1 if block else 0,
            "original_end": block[-1] + 1 if block else 0,
            "columns_retained": len(kept_cols),
            "columns_dropped": length - len(kept_cols),
            "tip_occupancy": round(marker_tip_occupancy[idx], 6),
            "removed_entirely": len(kept_cols) == 0,
            "removal_reason": "",
            "filtered_start": new_index[kept_cols[0]] if kept_cols else 0,
            "filtered_end": new_index[kept_cols[-1]] if kept_cols else 0,
        }
        if not kept_cols:
            row["removal_reason"] = (
                f"marker_tip_occupancy<{min_marker_tip_occupancy}"
                if marker_removed_by_rule[idx]
                else f"no_column_reached_occupancy>={min_column_occupancy}"
            )
        marker_rows.append(row)

    # Validation: every retained marker range is contiguous and the ranges tile 1..len exactly.
    retained_markers = [r for r in marker_rows if not r["removed_entirely"]]
    for row in retained_markers:
        span = row["filtered_end"] - row["filtered_start"] + 1
        if span != row["columns_retained"]:
            raise AssertionError(
                f"marker {row['marker_id']} is not contiguous in the filtered alignment: "
                f"span {span} != {row['columns_retained']} retained columns"
            )
    covered = 0
    expected_next = 1
    for row in sorted(retained_markers, key=lambda r: r["filtered_start"]):
        if row["filtered_start"] != expected_next:
            raise AssertionError(
                f"partition gap or overlap before marker {row['marker_id']}: expected filtered "
                f"column {expected_next}, got {row['filtered_start']}"
            )
        expected_next = row["filtered_end"] + 1
        covered += row["columns_retained"]
    if covered != len(retained_indices):
        raise AssertionError(
            f"partitions cover {covered} columns but {len(retained_indices)} were retained"
        )

    # Per-tip occupancy, before and after.
    tip_rows = []
    for t in tips:
        before = sum(1 for c in upper[t] if c in STANDARD_AA)
        after = sum(1 for c in filtered[t] if c in STANDARD_AA)
        tip_rows.append({
            "tip_id": t,
            "original_columns": width,
            "original_standard_aa": before,
            "original_occupancy": round(before / width, 6) if width else 0.0,
            "filtered_columns": len(retained_indices),
            "filtered_standard_aa": after,
            "filtered_occupancy": round(after / len(retained_indices), 6) if retained_indices else 0.0,
            "markers_with_any_residue": sum(
                1 for r, (_, _, length) in zip(marker_rows, markers)
                if length and any(
                    upper[t][i] in STANDARD_AA
                    for i in range(r["original_start"] - 1, r["original_end"])
                )
            ),
        })

    non_standard = collections.Counter(
        c for t in tips for c in upper[t] if c not in STANDARD_AA
    )

    return {
        "tips": tips,
        "filtered": filtered,
        "marker_rows": marker_rows,
        "tip_rows": tip_rows,
        "occupancy": occupancy,
        "keep": keep,
        "col_marker": col_marker,
        "new_index": new_index,
        "retained_indices": retained_indices,
        "columns": columns,
        "before_stats": column_stats(columns),
        "after_stats": column_stats(columns[i] for i in retained_indices),
        "non_standard_character_counts": dict(sorted(non_standard.items())),
        "markers": markers,
    }


def write_outputs(result: dict, out_dir: Path, params: dict) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    tips = result["tips"]
    markers = result["markers"]

    with (out_dir / "alignment.faa").open("w") as handle:
        for t in tips:
            handle.write(f">{t}\n{result['filtered'][t]}\n")

    with (out_dir / "column-map.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "original_column", "marker_id", "marker_column", "tips_with_standard_aa",
            "occupancy", "decision", "filtered_column",
        ])
        marker_offsets = {}
        offset = 0
        for idx, (_, _, length) in enumerate(markers):
            marker_offsets[idx] = offset
            offset += length
        n_tips = len(tips)
        for i, occ in enumerate(result["occupancy"]):
            m = result["col_marker"][i]
            kept = result["keep"][i]
            writer.writerow([
                i + 1,
                markers[m][0],
                i - marker_offsets[m] + 1,
                int(round(occ * n_tips)),
                f"{occ:.6f}",
                "retained" if kept else "dropped",
                result["new_index"].get(i, 0),
            ])

    with (out_dir / "per-marker-occupancy.tsv").open("w", newline="") as handle:
        fields = list(result["marker_rows"][0]) if result["marker_rows"] else []
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result["marker_rows"])

    with (out_dir / "per-tip-occupancy.tsv").open("w", newline="") as handle:
        fields = list(result["tip_rows"][0]) if result["tip_rows"] else []
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result["tip_rows"])

    retained_markers = [r for r in result["marker_rows"] if not r["removed_entirely"]]
    with (out_dir / "partitions.nex").open("w") as handle:
        handle.write("#nexus\n")
        handle.write("begin sets;\n")
        for row in retained_markers:
            handle.write(
                f"  charset {row['charset_name']} = "
                f"{row['filtered_start']}-{row['filtered_end']};\n"
            )
        handle.write("end;\n")

    summary = {
        "parameters": params,
        "tips": len(tips),
        "markers_in_marker_info": len(markers),
        "markers_retained": len(retained_markers),
        "markers_removed_entirely": [
            {"marker_id": r["marker_id"], "reason": r["removal_reason"]}
            for r in result["marker_rows"] if r["removed_entirely"]
        ],
        "before": result["before_stats"],
        "after": result["after_stats"],
        "columns_retained": len(result["retained_indices"]),
        "columns_dropped": result["before_stats"]["columns"] - len(result["retained_indices"]),
        "non_standard_character_counts": result["non_standard_character_counts"],
        "validation": {
            "marker_lengths_sum_to_alignment_width": True,
            "every_retained_marker_is_contiguous_in_the_filtered_alignment": True,
            "partitions_cover_retained_columns_exactly": True,
            "tips_unchanged": True,
            "invariant_columns_preserved": (
                "no variability criterion is applied; constant columns are dropped only if their "
                "occupancy is below the threshold"
            ),
        },
        "tip_occupancy_range": [
            min(r["filtered_occupancy"] for r in result["tip_rows"]),
            max(r["filtered_occupancy"] for r in result["tip_rows"]),
        ],
        "lowest_occupancy_tips": sorted(
            ({"tip_id": r["tip_id"], "filtered_occupancy": r["filtered_occupancy"],
              "markers_with_any_residue": r["markers_with_any_residue"]}
             for r in result["tip_rows"]),
            key=lambda r: r["filtered_occupancy"],
        )[:15],
    }
    (out_dir / "filter-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True, type=Path,
                        help="concatenated marker alignment (FASTA, optionally gzipped)")
    parser.add_argument("--marker-info", required=True, type=Path,
                        help="gtdbtk.bac120.marker_info.tsv describing the concatenation")
    parser.add_argument("--min-column-occupancy", type=float, default=0.8,
                        help="retain a column when this fraction of tips carries a standard amino "
                             "acid (default 0.8; a starting sensitivity setting, not an optimum)")
    parser.add_argument("--min-marker-tip-occupancy", type=float, default=0.0,
                        help="remove a marker entirely when fewer than this fraction of tips carry "
                             "any residue for it (default 0.0, i.e. off)")
    parser.add_argument("--expected-tips", type=Path, default=None,
                        help="file of expected tip ids, one per line; the run fails if the "
                             "alignment's tip set differs")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    if not 0.0 <= args.min_column_occupancy <= 1.0:
        parser.error("--min-column-occupancy must be between 0 and 1")
    if not 0.0 <= args.min_marker_tip_occupancy <= 1.0:
        parser.error("--min-marker-tip-occupancy must be between 0 and 1")
    if args.output.exists() and any(args.output.iterdir()):
        parser.error(f"refusing to write into non-empty {args.output}")

    seqs = read_fasta(args.alignment)
    markers = read_marker_info(args.marker_info)

    if args.expected_tips:
        expected = {l.strip() for l in args.expected_tips.read_text().splitlines() if l.strip()}
        got = set(seqs)
        if expected != got:
            raise SystemExit(
                f"tip set mismatch: {len(expected - got)} expected tips absent from the alignment "
                f"({sorted(expected - got)[:10]}), {len(got - expected)} unexpected "
                f"({sorted(got - expected)[:10]})"
            )

    result = filter_alignment(
        seqs, markers, args.min_column_occupancy, args.min_marker_tip_occupancy,
    )
    summary = write_outputs(result, args.output, {
        "alignment": str(args.alignment),
        "marker_info": str(args.marker_info),
        "min_column_occupancy": args.min_column_occupancy,
        "min_marker_tip_occupancy": args.min_marker_tip_occupancy,
    })
    json.dump(summary, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
