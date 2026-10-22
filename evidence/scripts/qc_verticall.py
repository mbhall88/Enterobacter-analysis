#!/usr/bin/env python3
"""QC for a Verticall pairwise run and the PHYLIP matrices derived from it.

This generalises the eight-genome checks in ``scripts/summarize_smoke.py`` to an
arbitrary frozen panel (Task 08 runs it on the 388-tip Task 06 panel, where
n(n-1) = 150,156 directed results span 75,078 unordered pairs).

What it is for
--------------
``verticall pairwise`` silently *omits* a directed pair whose alignment produced
no sliding windows (``pairwise.py``: ``if window_count == 0: return
all_log_text, []``), and ``verticall matrix`` then writes an **empty field** for
that cell and continues with exit status 0 (``matrix.py: save_matrix``, which
only emits a ``warning(...)``).  Worse, ``matrix.py: make_symmetrical`` copies a
single available direction into both cells, so an incomplete pairwise run can
produce a matrix that looks complete and perfectly symmetric.  Matrix symmetry
is therefore *not* evidence of complete directed-pair coverage, and this script
checks the two things separately.

Failure versus review
---------------------
Hard failures (non-zero exit) are **absent or unusable output**:

* a directed pair with no row at all (a zero-result pair);
* a row whose chosen distance is empty or non-finite;
* a matrix with a missing cell, a non-finite or negative value, a non-zero
  diagonal, the wrong dimension, or the wrong tip set;
* a matrix that is asymmetric when it was expected to be symmetric.

A **low aligned fraction is never a pair failure.** On a genus-wide,
type-anchored panel with two out-of-genus outgroups, a low aligned fraction is
the expected biological signal for a distant pair, not a defect. Low-AF pairs
are reported in review bands (``--af-warn``/``--af-review``) and counted per
assembly; no row is ever dropped, rescored or filtered by this script.

Nothing here imputes a value. A missing distance stays missing and is reported.

Usage
-----
    python3 scripts/qc_verticall.py \
        --pairwise  results/.../pairwise.tsv \
        --expected-ids manifests/pilot.ids.txt \
        --matrix primary=results/.../median-first.phylip \
        --matrix low=results/.../median-low.phylip \
        --matrix asymmetric=results/.../asymmetric.phylip \
        --expect-asymmetric asymmetric \
        --output results/.../qc

Exit status: 0 all hard checks pass, 1 at least one hard check failed,
2 usage/IO error.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import math
import statistics
import sys
from pathlib import Path

# Columns emitted by Verticall 0.4.2 (verticall/pairwise.py: get_table_header).
REQUIRED_COLUMNS = (
    "assembly_a", "assembly_b", "alignment_count", "n50_alignment_length",
    "aligned_fraction", "mean_distance", "window_size", "window_count",
    "mean_window_distance", "median_window_distance", "mass_peaks", "result_level",
    "peak_window_distance", "peak_mass", "alignments_vertical_fraction",
    "alignments_horizontal_fraction", "mean_vertical_window_distance",
    "median_vertical_window_distance", "mean_vertical_distance", "r/m",
    "assembly_a_vertical_fraction", "assembly_a_horizontal_fraction",
    "assembly_a_unaligned_fraction", "assembly_b_vertical_fraction",
    "assembly_b_horizontal_fraction", "assembly_b_unaligned_fraction",
)

# Columns written as a percent string with a trailing '%' (2 dp) rather than a fraction.
PERCENT_COLUMNS = (
    "alignments_vertical_fraction", "alignments_horizontal_fraction",
    "assembly_a_vertical_fraction", "assembly_a_horizontal_fraction",
    "assembly_a_unaligned_fraction", "assembly_b_vertical_fraction",
    "assembly_b_horizontal_fraction", "assembly_b_unaligned_fraction",
)

# verticall/paint.py: get_r_over_m returns these strings instead of a number.
RM_SENTINELS = ("undef", "inf")

# verticall/matrix.py: jukes_cantor() returns exactly 25.0 for any d >= 0.75.
# A cell at this value is a saturation sentinel, not a measured distance.
JUKES_CANTOR_CAP = 25.0

RESULT_LEVELS = ("primary", "secondary")


class QCError(Exception):
    """A usage or input error that prevents the QC from running at all."""


# --------------------------------------------------------------------------- parsing


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="", encoding="utf-8")


def parse_float(value):
    """Return a finite float, or None for an empty/unparseable/non-finite field."""
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def parse_percent(value):
    """Parse Verticall's '12.34%' percent strings. Returns percent units, not a fraction."""
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    if text.endswith("%"):
        text = text[:-1]
    return parse_float(text)


def parse_r_over_m(value):
    """Return (numeric_or_None, sentinel_or_None) for the r/m column."""
    text = (value or "").strip()
    if text in RM_SENTINELS:
        return None, text
    return parse_float(text), None


def count_mass_peaks(value):
    """Number of detected peaks in the comma-delimited mass_peaks field."""
    text = (value or "").strip()
    if not text:
        return 0
    return len([part for part in text.split(",") if part.strip()])


def read_expected_ids(path: Path):
    ids, seen = [], set()
    with open_text(path) as handle:
        for line in handle:
            tip = line.strip()
            if not tip or tip.startswith("#"):
                continue
            if tip in seen:
                raise QCError(f"duplicate tip id {tip!r} in {path}")
            seen.add(tip)
            ids.append(tip)
    if not ids:
        raise QCError(f"no tip ids found in {path}")
    return ids


def load_pairwise(path: Path):
    """Read pairwise.tsv, dropping the six painted-region columns.

    Those columns are about 98% of the file - 16.5 kB of the 16.8 kB an
    eight-genome row occupied - so at 388 tips the file is roughly 2.5 GB while
    the columns this QC actually reads are about 48 MB. Nothing here interprets
    a region string, so they are never loaded. ``scripts/summarize_smoke.py``
    made the same exclusion when it wrote its reduced table.
    """
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise QCError(f"{path} is empty")
        missing = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise QCError(f"{path} is missing expected Verticall columns: {', '.join(missing)}")
        kept = [c for c in reader.fieldnames if not c.endswith("_regions")]
        rows = [{column: record.get(column, "") for column in kept} for record in reader]
    if not rows:
        raise QCError(f"{path} has a header but no result rows")
    return rows, kept


# --------------------------------------------------------------------------- statistics


def percentile(sorted_values, fraction):
    """Linear-interpolated percentile of an already-sorted, non-empty list."""
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def summarise(values):
    """Descriptive summary of a list of floats; None entries are counted, not imputed."""
    present = sorted(v for v in values if v is not None)
    absent = sum(1 for v in values if v is None)
    if not present:
        return {"n": 0, "n_missing": absent}
    return {
        "n": len(present),
        "n_missing": absent,
        "min": present[0],
        "p05": percentile(present, 0.05),
        "q1": percentile(present, 0.25),
        "median": percentile(present, 0.50),
        "q3": percentile(present, 0.75),
        "p95": percentile(present, 0.95),
        "max": present[-1],
        "mean": statistics.fmean(present),
        "n_zero": sum(1 for v in present if v == 0.0),
    }


def histogram(values, edges):
    """Counts in [edges[i], edges[i+1]) with a final closed bin. Returns a list of dicts."""
    present = [v for v in values if v is not None]
    bins = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i == len(edges) - 2:
            count = sum(1 for v in present if lo <= v <= hi)
        else:
            count = sum(1 for v in present if lo <= v < hi)
        bins.append({"from": lo, "to": hi, "count": count})
    return bins


# --------------------------------------------------------------------------- PHYLIP


def load_phylip(path: Path):
    """Parse a Verticall-written PHYLIP matrix.

    ``verticall matrix`` writes ``n`` then one tab-delimited ``name\\tv1...vn``
    line per sample, using an **empty field** where a distance is missing. The
    empty field is preserved here as ``None`` and is never replaced with zero.

    Returns (declared_n, names, rows, problems) where ``rows`` is a list of
    lists of float-or-None and ``problems`` lists structural parse errors.
    """
    problems = []
    with open_text(path) as handle:
        lines = handle.read().splitlines()
    lines = [line for line in lines if line.strip()]
    if not lines:
        raise QCError(f"{path} is empty")
    try:
        declared = int(lines[0].strip().split()[0])
    except (ValueError, IndexError):
        raise QCError(f"{path}: first line {lines[0]!r} is not a sample count")

    names, rows = [], []
    for line_number, line in enumerate(lines[1:], start=2):
        if "\t" in line:
            fields = line.split("\t")
        else:
            # Tolerate a whitespace-delimited matrix from another writer. An empty
            # cell is not representable in that form, so a short row is an error.
            fields = line.split()
        name, values = fields[0], fields[1:]
        parsed = []
        for column, raw in enumerate(values, start=1):
            text = raw.strip()
            if text == "":
                parsed.append(None)
                continue
            try:
                number = float(text)
            except ValueError:
                problems.append(
                    f"line {line_number} column {column}: {raw!r} is not a number")
                parsed.append(None)
                continue
            if not math.isfinite(number):
                problems.append(
                    f"line {line_number} column {column}: {text} is not finite")
                parsed.append(None)
                continue
            parsed.append(number)
        names.append(name)
        rows.append(parsed)
    return declared, names, rows, problems


def check_matrix(label, path, expected_ids, expect_symmetric, symmetry_tolerance):
    """Validate one PHYLIP matrix. Returns a report dict with ``status``."""
    report = {
        "label": label,
        "path": str(path),
        "expect_symmetric": expect_symmetric,
        "failures": [],
        "warnings": [],
    }
    if not path.exists():
        report["failures"].append("matrix file does not exist")
        report["status"] = "fail"
        return report

    declared, names, rows, problems = load_phylip(path)
    report["declared_sample_count"] = declared
    report["row_count"] = len(names)
    report["parse_problems"] = problems
    if problems:
        report["failures"].append(
            f"{len(problems)} unparseable or non-finite cells (first: {problems[0]})")

    expected_n = len(expected_ids)
    if declared != len(names):
        report["failures"].append(
            f"header declares {declared} samples but the file has {len(names)} rows")
    if len(names) != expected_n:
        report["failures"].append(
            f"matrix has {len(names)} tips, expected {expected_n}")

    duplicate_names = sorted(n for n, c in collections.Counter(names).items() if c > 1)
    if duplicate_names:
        report["failures"].append(
            f"{len(duplicate_names)} duplicated tip names (first: {duplicate_names[0]})")
    report["duplicate_tip_names"] = duplicate_names

    expected_set, observed_set = set(expected_ids), set(names)
    missing_tips = sorted(expected_set - observed_set)
    extra_tips = sorted(observed_set - expected_set)
    report["tips_missing_from_matrix"] = missing_tips
    report["tips_not_in_expected_panel"] = extra_tips
    if missing_tips:
        report["failures"].append(
            f"{len(missing_tips)} expected tips absent (first: {missing_tips[0]})")
    if extra_tips:
        report["failures"].append(
            f"{len(extra_tips)} tips are not in the expected panel (first: {extra_tips[0]})")

    ragged = [names[i] for i, row in enumerate(rows) if len(row) != len(names)]
    report["rows_with_wrong_width"] = ragged
    if ragged:
        report["failures"].append(
            f"{len(ragged)} rows are not {len(names)} cells wide (first: {ragged[0]})")

    missing_cells, negative_cells, capped_cells = [], [], []
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            if j >= len(names):
                break
            if value is None:
                missing_cells.append([names[i], names[j]])
            elif value < 0:
                negative_cells.append([names[i], names[j], value])
            elif value == JUKES_CANTOR_CAP:
                capped_cells.append([names[i], names[j]])

    report["missing_cell_count"] = len(missing_cells)
    report["missing_cells_sample"] = missing_cells[:20]
    report["negative_cell_count"] = len(negative_cells)
    report["negative_cells_sample"] = negative_cells[:20]
    report["jukes_cantor_capped_cell_count"] = len(capped_cells)
    report["jukes_cantor_capped_cells_sample"] = capped_cells[:20]

    if missing_cells:
        # verticall matrix writes an empty field and exits 0. This is the gate.
        report["failures"].append(
            f"{len(missing_cells)} missing distances: the matrix is incomplete and is "
            "rejected. Missing distances are not zero and are not imputed")
    if negative_cells:
        report["failures"].append(f"{len(negative_cells)} negative distances")
    if capped_cells:
        report["warnings"].append(
            f"{len(capped_cells)} cells are exactly {JUKES_CANTOR_CAP}, the Jukes-Cantor "
            "saturation cap applied to any uncorrected distance >= 0.75; these are "
            "sentinels, not measured distances")

    diagonal_problems = []
    for i in range(min(len(rows), len(names))):
        if i < len(rows[i]):
            value = rows[i][i]
            if value is None or value != 0.0:
                diagonal_problems.append([names[i], value])
    report["nonzero_diagonal_count"] = len(diagonal_problems)
    report["nonzero_diagonal_sample"] = diagonal_problems[:20]
    if diagonal_problems:
        report["failures"].append(
            f"{len(diagonal_problems)} diagonal cells are not exactly zero")

    max_asymmetry, asymmetry_pair = 0.0, None
    compared = 0
    for i in range(len(rows)):
        for j in range(i + 1, len(names)):
            if j >= len(rows[i]) or i >= len(rows[j]):
                continue
            a, b = rows[i][j], rows[j][i]
            if a is None or b is None:
                continue
            compared += 1
            gap = abs(a - b)
            if gap > max_asymmetry:
                max_asymmetry, asymmetry_pair = gap, [names[i], names[j], a, b]
    report["symmetry_comparisons"] = compared
    report["max_asymmetry"] = max_asymmetry
    report["max_asymmetry_pair"] = asymmetry_pair

    if expect_symmetric:
        if max_asymmetry > symmetry_tolerance:
            report["failures"].append(
                f"matrix is asymmetric (max |d(a,b)-d(b,a)| = {max_asymmetry:.3e} > "
                f"{symmetry_tolerance:.3e}) but was declared symmetric")
    else:
        if compared and max_asymmetry <= symmetry_tolerance:
            report["warnings"].append(
                "matrix was declared asymmetric (--asymmetrical) but every reciprocal "
                "pair is identical; check that the intended flag was used")

    finite_values = [v for i, row in enumerate(rows) for j, v in enumerate(row)
                     if v is not None and i != j and j < len(names)]
    report["offdiagonal_distance_summary"] = summarise(finite_values)
    report["status"] = "fail" if report["failures"] else "pass"
    return report


# --------------------------------------------------------------------------- pairwise QC


def index_pairwise(rows, expected_ids):
    """Group rows by directed pair and classify row multiplicity.

    Distinguishes three things the plan asks be kept apart:

    * a **zero-result** pair - no row at all, because ``window_count == 0``;
    * **secondary accepted results** - extra rows with ``result_level ==
      'secondary'`` alongside the one ``primary`` row, which is normal output
      for a pair with a close second peak;
    * **duplicate** rows - more than one ``primary`` row for the same directed
      pair, which Verticall cannot emit and which means the TSV was
      concatenated or ``--part`` ranges overlapped.

    Multiple *detected* peaks (``mass_peaks``) is a fourth, separate quantity: a
    pair can detect several peaks and still accept exactly one result.
    """
    expected_set = set(expected_ids)
    by_pair = collections.defaultdict(list)
    self_pairs, unknown_tip_rows = [], []
    for row_number, row in enumerate(rows, start=2):
        a, b = row["assembly_a"].strip(), row["assembly_b"].strip()
        if a == b:
            self_pairs.append([row_number, a])
            continue
        if a not in expected_set or b not in expected_set:
            unknown_tip_rows.append([row_number, a, b])
            continue
        by_pair[(a, b)].append(row)
    return by_pair, self_pairs, unknown_tip_rows


def classify_pair_rows(pair_rows, distance_column, af_warn, af_review):
    """Per-directed-pair facts used by both the review table and the summary."""
    levels = collections.Counter((r["result_level"] or "").strip() for r in pair_rows)
    primary_rows = [r for r in pair_rows if (r["result_level"] or "").strip() == "primary"]
    secondary_rows = [r for r in pair_rows if (r["result_level"] or "").strip() == "secondary"]
    other_rows = [r for r in pair_rows
                  if (r["result_level"] or "").strip() not in RESULT_LEVELS]

    # verticall matrix --multi first takes the first row for the pair, and
    # pairwise.py sorts results most-to-least massive, so 'first' is the primary
    # result whenever the TSV has not been reordered. Record both so a reordered
    # file is visible rather than silently changing which distance is used.
    first_row = pair_rows[0] if pair_rows else None
    chosen = primary_rows[0] if primary_rows else first_row

    detected_peaks = max((count_mass_peaks(r["mass_peaks"]) for r in pair_rows), default=0)
    af = parse_float(chosen["aligned_fraction"]) if chosen is not None else None
    distance = parse_float(chosen[distance_column]) if chosen is not None else None

    flags = []
    if not pair_rows:
        flags.append("zero_result")
    if len(primary_rows) > 1:
        flags.append("duplicate_primary_rows")
    if secondary_rows:
        flags.append("secondary_accepted_results")
    if other_rows:
        flags.append("unrecognised_result_level")
    if detected_peaks > 1 and len(pair_rows) == 1:
        flags.append("multiple_detected_peaks_one_accepted_result")
    if distance is None:
        flags.append("distance_missing_or_nonfinite")
    if af is not None:
        if af < af_review:
            flags.append("aligned_fraction_review")
        elif af < af_warn:
            flags.append("aligned_fraction_warn")
    return {
        "rows": len(pair_rows),
        "primary_rows": len(primary_rows),
        "secondary_rows": len(secondary_rows),
        "other_level_rows": len(other_rows),
        "result_levels": dict(levels),
        "detected_peaks": detected_peaks,
        "first_row_is_primary": bool(primary_rows) and first_row is primary_rows[0],
        "aligned_fraction": af,
        "distance": distance,
        "peak_mass": parse_float(chosen["peak_mass"]) if chosen is not None else None,
        "flags": flags,
        "chosen_row": chosen,
    }


def relative_difference(a, b):
    """Symmetric relative difference; None when both directions are ~zero."""
    if a is None or b is None:
        return None
    scale = max(abs(a), abs(b))
    if scale == 0.0:
        return 0.0
    return abs(a - b) / scale


# --------------------------------------------------------------------------- assemblies


def blank_observation():
    return {
        "vertical_percent": [], "horizontal_percent": [], "unaligned_percent": [],
        "aligned_fraction_as_query": [], "primary_peak_mass": [],
        "rows_as_query": 0, "rows_as_target": 0,
    }


def assembly_observations(rows, expected_ids):
    """Collect per-assembly painted fractions, coverage and primary peak mass."""
    obs = {tip: blank_observation() for tip in expected_ids}
    for row in rows:
        a, b = row["assembly_a"].strip(), row["assembly_b"].strip()
        level = (row["result_level"] or "").strip()
        for tip, side in ((a, "a"), (b, "b")):
            if tip not in obs:
                continue
            obs[tip]["vertical_percent"].append(
                parse_percent(row[f"assembly_{side}_vertical_fraction"]))
            obs[tip]["horizontal_percent"].append(
                parse_percent(row[f"assembly_{side}_horizontal_fraction"]))
            obs[tip]["unaligned_percent"].append(
                parse_percent(row[f"assembly_{side}_unaligned_fraction"]))
            if level == "primary":
                obs[tip]["primary_peak_mass"].append(parse_float(row["peak_mass"]))
        if a in obs:
            obs[a]["rows_as_query"] += 1
            # aligned_fraction is the query coverage of assembly_a, so it is a
            # directional property of a, not a property of the unordered pair.
            obs[a]["aligned_fraction_as_query"].append(parse_float(row["aligned_fraction"]))
        if b in obs:
            obs[b]["rows_as_target"] += 1
    return obs


# --------------------------------------------------------------------------- output


def write_tsv(path: Path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value, digits=9):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


PAIR_FIELDS = [
    "assembly_a", "assembly_b", "status", "flags",
    "rows_ab", "rows_ba", "primary_rows_ab", "primary_rows_ba",
    "secondary_rows_ab", "secondary_rows_ba",
    "detected_peaks_ab", "detected_peaks_ba",
    "aligned_fraction_ab", "aligned_fraction_ba", "aligned_fraction_abs_difference",
    "distance_ab", "distance_ba",
    "distance_abs_reciprocal_difference", "distance_rel_reciprocal_difference",
    "distance_mean_before_jukes_cantor",
    "peak_mass_ab", "peak_mass_ba",
]

ASSEMBLY_FIELDS = [
    "assembly_id", "status", "flags",
    "directed_rows_as_query", "directed_rows_as_target",
    "expected_directed_pairs_each_direction",
    "zero_result_pairs", "pairs_with_secondary_results", "pairs_with_duplicate_rows",
    "pairs_aligned_fraction_below_warn", "pairs_aligned_fraction_below_review",
    "aligned_fraction_as_query_min", "aligned_fraction_as_query_median",
    "aligned_fraction_as_query_max",
    "vertical_percent_min", "vertical_percent_median", "vertical_percent_max",
    "horizontal_percent_median", "horizontal_percent_max",
    "unaligned_percent_min", "unaligned_percent_median", "unaligned_percent_max",
    "primary_peak_mass_min", "primary_peak_mass_median",
]


def run_qc(args):
    expected_ids = read_expected_ids(args.expected_ids)
    n = len(expected_ids)
    expected_directed = n * (n - 1)
    expected_unordered = n * (n - 1) // 2

    labels = [label for label, _ in args.matrix]
    repeated = sorted(label for label, count in collections.Counter(labels).items() if count > 1)
    if repeated:
        raise QCError(f"--matrix labels must be unique; repeated: {', '.join(repeated)}")
    unknown_asymmetric = sorted(set(args.expect_asymmetric) - set(labels))
    if unknown_asymmetric:
        raise QCError("--expect-asymmetric names a label with no --matrix: "
                      + ", ".join(unknown_asymmetric))

    rows, fieldnames = load_pairwise(args.pairwise)
    if args.distance_column not in fieldnames:
        raise QCError(f"distance column {args.distance_column!r} is not in {args.pairwise}")

    by_pair, self_pairs, unknown_tip_rows = index_pairwise(rows, expected_ids)
    classified = {
        pair: classify_pair_rows(pair_rows, args.distance_column, args.af_warn, args.af_review)
        for pair, pair_rows in by_pair.items()
    }

    expected_pairs = [(a, b) for a in expected_ids for b in expected_ids if a != b]
    missing_directed = [pair for pair in expected_pairs if pair not in classified]

    # Exact duplicate row content is a different error from two primary rows.
    row_signatures = collections.Counter(
        tuple(row.get(column, "") for column in fieldnames) for row in rows)
    exact_duplicate_rows = sum(count - 1 for count in row_signatures.values() if count > 1)

    failures, warnings = [], []
    if missing_directed:
        failures.append(
            f"{len(missing_directed)} of {expected_directed} directed pairs have no result "
            "row (zero-result pairs). A pair with no row is not a distance of zero")
    if self_pairs:
        failures.append(f"{len(self_pairs)} rows compare an assembly with itself")
    if unknown_tip_rows:
        failures.append(
            f"{len(unknown_tip_rows)} rows name a tip that is not in the expected panel")
    if exact_duplicate_rows:
        failures.append(
            f"{exact_duplicate_rows} byte-identical duplicate rows: the TSV looks concatenated")

    duplicate_pairs = [p for p, c in classified.items() if c["primary_rows"] > 1]
    if duplicate_pairs:
        failures.append(
            f"{len(duplicate_pairs)} directed pairs have more than one primary row; "
            "Verticall emits exactly one, so --part ranges probably overlapped")

    bad_distance_pairs = [p for p, c in classified.items()
                          if "distance_missing_or_nonfinite" in c["flags"]]
    if bad_distance_pairs:
        failures.append(
            f"{len(bad_distance_pairs)} directed pairs have an empty or non-finite "
            f"{args.distance_column}")

    unrecognised = [p for p, c in classified.items() if c["other_level_rows"]]
    if unrecognised:
        failures.append(
            f"{len(unrecognised)} directed pairs carry a result_level that is neither "
            "primary nor secondary")

    reordered = [p for p, c in classified.items()
                 if c["primary_rows"] == 1 and not c["first_row_is_primary"]]
    if reordered:
        warnings.append(
            f"{len(reordered)} directed pairs do not list their primary result first; "
            "'verticall matrix --multi first' would take a secondary distance for them")

    secondary_pairs = [p for p, c in classified.items() if c["secondary_rows"]]
    multipeak_single_result = [p for p, c in classified.items()
                               if "multiple_detected_peaks_one_accepted_result" in c["flags"]]
    low_af_warn = [p for p, c in classified.items() if "aligned_fraction_warn" in c["flags"]]
    low_af_review = [p for p, c in classified.items() if "aligned_fraction_review" in c["flags"]]
    if low_af_warn or low_af_review:
        warnings.append(
            f"{len(low_af_warn) + len(low_af_review)} directed pairs are below the "
            f"aligned-fraction review bands (<{args.af_warn} / <{args.af_review}). "
            "Low aligned fraction is a biological review category on a genus-wide "
            "panel, not a pair failure; no row was dropped")

    # ---- reciprocal differences, computed before any averaging -------------
    # Rows are streamed to disk rather than accumulated: a run that died partway
    # through can flag every one of the 75,078 unordered pairs at 388 tips.
    args.output.mkdir(parents=True, exist_ok=True)
    reciprocal_distance_gaps, reciprocal_af_gaps = [], []
    one_directional = []
    pair_review_rows = 0
    pair_handle = (args.output / "pair-review.tsv").open("w", newline="", encoding="utf-8")
    pair_writer = csv.DictWriter(pair_handle, fieldnames=PAIR_FIELDS, delimiter="\t",
                                 lineterminator="\n", extrasaction="ignore")
    pair_writer.writeheader()
    for i, a in enumerate(expected_ids):
        for b in expected_ids[i + 1:]:
            ab, ba = classified.get((a, b)), classified.get((b, a))
            flags = set()
            for direction, result in (("ab", ab), ("ba", ba)):
                if result is None:
                    flags.add(f"zero_result_{direction}")
                else:
                    flags.update(f"{flag}_{direction}" for flag in result["flags"])
            d_ab = ab["distance"] if ab else None
            d_ba = ba["distance"] if ba else None
            af_ab = ab["aligned_fraction"] if ab else None
            af_ba = ba["aligned_fraction"] if ba else None

            gap = abs(d_ab - d_ba) if d_ab is not None and d_ba is not None else None
            rel = relative_difference(d_ab, d_ba)
            if gap is not None:
                reciprocal_distance_gaps.append(gap)
            if af_ab is not None and af_ba is not None:
                reciprocal_af_gaps.append(abs(af_ab - af_ba))
            if (ab is None) != (ba is None):
                one_directional.append([a, b])
                flags.add("one_direction_only")
            if gap is not None and gap > args.reciprocal_warn:
                flags.add("reciprocal_distance_disagreement")

            mean_before_jc = None
            if d_ab is not None and d_ba is not None:
                mean_before_jc = (d_ab + d_ba) / 2.0
            elif d_ab is not None or d_ba is not None:
                # matrix.py: make_symmetrical copies the single available
                # direction into both cells rather than leaving it missing.
                mean_before_jc = d_ab if d_ab is not None else d_ba

            hard = {f for f in flags if f.startswith("zero_result")
                    or f.startswith("duplicate_primary_rows")
                    or f.startswith("distance_missing_or_nonfinite")
                    or f.startswith("unrecognised_result_level")
                    or f == "one_direction_only"}
            status = "fail" if hard else ("review" if flags else "ok")
            if status == "ok" and not args.write_all_pairs:
                continue
            pair_review_rows += 1
            pair_writer.writerow({
                "assembly_a": a, "assembly_b": b, "status": status,
                "flags": ";".join(sorted(flags)),
                "rows_ab": ab["rows"] if ab else 0, "rows_ba": ba["rows"] if ba else 0,
                "primary_rows_ab": ab["primary_rows"] if ab else 0,
                "primary_rows_ba": ba["primary_rows"] if ba else 0,
                "secondary_rows_ab": ab["secondary_rows"] if ab else 0,
                "secondary_rows_ba": ba["secondary_rows"] if ba else 0,
                "detected_peaks_ab": ab["detected_peaks"] if ab else 0,
                "detected_peaks_ba": ba["detected_peaks"] if ba else 0,
                "aligned_fraction_ab": fmt(af_ab), "aligned_fraction_ba": fmt(af_ba),
                "aligned_fraction_abs_difference": fmt(
                    abs(af_ab - af_ba) if af_ab is not None and af_ba is not None else None),
                "distance_ab": fmt(d_ab), "distance_ba": fmt(d_ba),
                "distance_abs_reciprocal_difference": fmt(gap),
                "distance_rel_reciprocal_difference": fmt(rel),
                "distance_mean_before_jukes_cantor": fmt(mean_before_jc),
                "peak_mass_ab": fmt(ab["peak_mass"] if ab else None),
                "peak_mass_ba": fmt(ba["peak_mass"] if ba else None),
            })
    pair_handle.close()
    if one_directional:
        failures.append(
            f"{len(one_directional)} unordered pairs were measured in one direction only. "
            "'verticall matrix' would copy that single direction into both cells and the "
            "matrix would look complete")

    # ---- per-assembly ------------------------------------------------------
    obs = assembly_observations(rows, expected_ids)
    zero_by_assembly = collections.Counter()
    for a, b in missing_directed:
        zero_by_assembly[a] += 1
        zero_by_assembly[b] += 1
    secondary_by_assembly = collections.Counter()
    for a, b in secondary_pairs:
        secondary_by_assembly[a] += 1
        secondary_by_assembly[b] += 1
    duplicate_by_assembly = collections.Counter()
    for a, b in duplicate_pairs:
        duplicate_by_assembly[a] += 1
        duplicate_by_assembly[b] += 1
    lowaf_warn_by_assembly = collections.Counter(a for a, _ in low_af_warn)
    lowaf_review_by_assembly = collections.Counter(a for a, _ in low_af_review)

    assembly_table = []
    for tip in expected_ids:
        record = obs[tip]
        af_stats = summarise(record["aligned_fraction_as_query"])
        vertical = summarise(record["vertical_percent"])
        horizontal = summarise(record["horizontal_percent"])
        unaligned = summarise(record["unaligned_percent"])
        peak_mass = summarise(record["primary_peak_mass"])
        flags = []
        if zero_by_assembly[tip]:
            flags.append("has_zero_result_pairs")
        if duplicate_by_assembly[tip]:
            flags.append("has_duplicate_rows")
        if record["rows_as_query"] == 0 and record["rows_as_target"] == 0:
            flags.append("absent_from_pairwise_output")
        if lowaf_review_by_assembly[tip]:
            flags.append("low_aligned_fraction_review")
        elif lowaf_warn_by_assembly[tip]:
            flags.append("low_aligned_fraction_warn")
        if unaligned.get("median") is not None and \
                unaligned["median"] > args.unaligned_review_percent:
            flags.append("high_median_unaligned_fraction")
        if peak_mass.get("median") is not None and peak_mass["median"] < args.peak_mass_review:
            flags.append("low_median_primary_peak_mass")
        status = "fail" if ("has_zero_result_pairs" in flags
                            or "has_duplicate_rows" in flags
                            or "absent_from_pairwise_output" in flags) else (
            "review" if flags else "ok")
        assembly_table.append({
            "assembly_id": tip, "status": status, "flags": ";".join(flags),
            "directed_rows_as_query": record["rows_as_query"],
            "directed_rows_as_target": record["rows_as_target"],
            "expected_directed_pairs_each_direction": n - 1,
            "zero_result_pairs": zero_by_assembly[tip],
            "pairs_with_secondary_results": secondary_by_assembly[tip],
            "pairs_with_duplicate_rows": duplicate_by_assembly[tip],
            "pairs_aligned_fraction_below_warn": lowaf_warn_by_assembly[tip],
            "pairs_aligned_fraction_below_review": lowaf_review_by_assembly[tip],
            "aligned_fraction_as_query_min": fmt(af_stats.get("min")),
            "aligned_fraction_as_query_median": fmt(af_stats.get("median")),
            "aligned_fraction_as_query_max": fmt(af_stats.get("max")),
            "vertical_percent_min": fmt(vertical.get("min"), 2),
            "vertical_percent_median": fmt(vertical.get("median"), 2),
            "vertical_percent_max": fmt(vertical.get("max"), 2),
            "horizontal_percent_median": fmt(horizontal.get("median"), 2),
            "horizontal_percent_max": fmt(horizontal.get("max"), 2),
            "unaligned_percent_min": fmt(unaligned.get("min"), 2),
            "unaligned_percent_median": fmt(unaligned.get("median"), 2),
            "unaligned_percent_max": fmt(unaligned.get("max"), 2),
            "primary_peak_mass_min": fmt(peak_mass.get("min")),
            "primary_peak_mass_median": fmt(peak_mass.get("median")),
        })

    # ---- distributions -----------------------------------------------------
    af_values = [parse_float(r["aligned_fraction"]) for r in rows]
    distance_values = [parse_float(r[args.distance_column]) for r in rows]
    peak_mass_values = [parse_float(r["peak_mass"]) for r in rows]
    rm_numeric, rm_sentinels = [], collections.Counter()
    for row in rows:
        number, sentinel = parse_r_over_m(row["r/m"])
        if sentinel:
            rm_sentinels[sentinel] += 1
        else:
            rm_numeric.append(number)

    # ---- matrices ----------------------------------------------------------
    matrix_reports = []
    for label, path in args.matrix:
        report = check_matrix(label, path, expected_ids,
                              label not in args.expect_asymmetric, args.symmetry_tolerance)
        matrix_reports.append(report)
        if report["status"] == "fail":
            failures.append(f"matrix '{label}' failed: {report['failures'][0]}")
        for warning in report["warnings"]:
            warnings.append(f"matrix '{label}': {warning}")

    summary = {
        "tool": "scripts/qc_verticall.py",
        "inputs": {
            "pairwise_tsv": str(args.pairwise),
            "expected_ids": str(args.expected_ids),
            "distance_column": args.distance_column,
            "matrices": {label: str(path) for label, path in args.matrix},
            "matrices_expected_asymmetric": sorted(args.expect_asymmetric),
        },
        "thresholds": {
            "aligned_fraction_warn": args.af_warn,
            "aligned_fraction_review": args.af_review,
            "reciprocal_distance_warn": args.reciprocal_warn,
            "unaligned_review_percent": args.unaligned_review_percent,
            "peak_mass_review": args.peak_mass_review,
            "matrix_symmetry_tolerance": args.symmetry_tolerance,
            "note": "These are reporting bands, not filters. No row is dropped, "
                    "rescored or imputed by this script.",
        },
        "pair_accounting": {
            "tips": n,
            "expected_directed_pairs": expected_directed,
            "expected_unordered_pairs": expected_unordered,
            "observed_directed_pairs": len(classified),
            "result_rows": len(rows),
            "zero_result_directed_pairs": len(missing_directed),
            "directed_pairs_with_secondary_results": len(secondary_pairs),
            "directed_pairs_with_duplicate_primary_rows": len(duplicate_pairs),
            "exact_duplicate_rows": exact_duplicate_rows,
            "self_comparison_rows": len(self_pairs),
            "rows_naming_a_tip_outside_the_panel": len(unknown_tip_rows),
            "unordered_pairs_measured_in_one_direction_only": len(one_directional),
            "directed_pairs_not_listing_primary_first": len(reordered),
            "result_level_counts": dict(collections.Counter(
                (r["result_level"] or "").strip() for r in rows)),
            "detected_peak_counts": dict(collections.Counter(
                count_mass_peaks(r["mass_peaks"]) for r in rows)),
            "directed_pairs_multiple_detected_peaks_one_accepted_result":
                len(multipeak_single_result),
        },
        "distributions": {
            "aligned_fraction": summarise(af_values),
            "aligned_fraction_histogram": histogram(
                af_values, [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]),
            args.distance_column: summarise(distance_values),
            "peak_mass": summarise(peak_mass_values),
            "r_over_m_numeric": summarise(rm_numeric),
            "r_over_m_sentinels": dict(rm_sentinels),
            "reciprocal_distance_absolute_difference": summarise(reciprocal_distance_gaps),
            "reciprocal_aligned_fraction_absolute_difference":
                summarise(reciprocal_af_gaps),
            "reciprocal_note": "aligned_fraction is the query coverage of assembly_a, so a "
                               "reciprocal difference is expected and reflects genome size and "
                               "coverage, not an error. A reciprocal difference in the distance "
                               "is measured before matrix averaging.",
        },
        "aligned_fraction_review": {
            "directed_pairs_below_warn": len(low_af_warn),
            "directed_pairs_below_review": len(low_af_review),
            "policy": "Low aligned fraction is a biological warning/review category, never "
                      "an automatic pair failure. Absent or non-finite output is a failure.",
        },
        "matrices": matrix_reports,
        "failures": failures,
        "warnings": warnings,
    }
    summary["status"] = "fail" if failures else "pass"
    summary["outputs"] = {
        "pair_review_rows": pair_review_rows,
        "pair_review_scope": "every unordered pair" if args.write_all_pairs
                             else "flagged unordered pairs only",
        "assembly_review_rows": len(assembly_table),
    }

    (args.output / "qc-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    write_tsv(args.output / "assembly-review.tsv", ASSEMBLY_FIELDS, assembly_table)
    if missing_directed:
        write_tsv(args.output / "zero-result-pairs.tsv",
                  ["assembly_a", "assembly_b"],
                  [{"assembly_a": a, "assembly_b": b} for a, b in sorted(missing_directed)])
    return summary


def parse_matrix_argument(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"--matrix expects LABEL=PATH, got {value!r}")
    label, path = value.split("=", 1)
    label = label.strip()
    if not label or not path.strip():
        raise argparse.ArgumentTypeError(f"--matrix expects LABEL=PATH, got {value!r}")
    return label, Path(path.strip())


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pairwise", type=Path, required=True,
                        help="pairwise.tsv written by 'verticall pairwise'")
    parser.add_argument("--expected-ids", type=Path, required=True,
                        help="one expected tip id per line (manifests/pilot.ids.txt)")
    parser.add_argument("--output", type=Path, required=True,
                        help="output directory for the summary and review tables")
    parser.add_argument("--matrix", type=parse_matrix_argument, action="append", default=[],
                        metavar="LABEL=PATH",
                        help="a PHYLIP matrix to validate; repeatable")
    parser.add_argument("--expect-asymmetric", action="append", default=[], metavar="LABEL",
                        help="matrix label produced with 'verticall matrix --asymmetrical'")
    parser.add_argument("--distance-column", default="median_vertical_window_distance",
                        help="pairwise column matching the matrix --distance_type "
                             "(default: %(default)s)")
    parser.add_argument("--af-warn", type=float, default=0.5,
                        help="aligned fractions below this are a review band (default: %(default)s)")
    parser.add_argument("--af-review", type=float, default=0.1,
                        help="aligned fractions below this are the lower review band "
                             "(default: %(default)s)")
    parser.add_argument("--reciprocal-warn", type=float, default=0.005,
                        help="flag an unordered pair whose two directed distances differ by "
                             "more than this before averaging (default: %(default)s)")
    parser.add_argument("--unaligned-review-percent", type=float, default=50.0,
                        help="flag an assembly whose median unaligned fraction exceeds this "
                             "percent (default: %(default)s)")
    parser.add_argument("--peak-mass-review", type=float, default=0.5,
                        help="flag an assembly whose median primary peak mass is below this "
                             "(default: %(default)s)")
    parser.add_argument("--symmetry-tolerance", type=float, default=1e-9,
                        help="maximum allowed |d(a,b)-d(b,a)| in a symmetric matrix "
                             "(default: %(default)s)")
    parser.add_argument("--write-all-pairs", action="store_true",
                        help="write every unordered pair to pair-review.tsv, not only flagged "
                             "ones (75,078 rows at 388 tips)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.af_review > args.af_warn:
        print("Error: --af-review must not exceed --af-warn", file=sys.stderr)
        return 2
    try:
        summary = run_qc(args)
    except QCError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    except OSError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    accounting = summary["pair_accounting"]
    print(f"status: {summary['status']}")
    print(f"tips: {accounting['tips']}  "
          f"expected directed pairs: {accounting['expected_directed_pairs']}  "
          f"expected unordered pairs: {accounting['expected_unordered_pairs']}")
    print(f"observed directed pairs: {accounting['observed_directed_pairs']}  "
          f"rows: {accounting['result_rows']}  "
          f"zero-result pairs: {accounting['zero_result_directed_pairs']}")
    for failure in summary["failures"]:
        print(f"FAIL: {failure}")
    for warning in summary["warnings"]:
        print(f"warn: {warning}")
    return 1 if summary["status"] == "fail" else 0


if __name__ == "__main__":
    sys.exit(main())
