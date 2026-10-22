#!/usr/bin/env python3
"""Report the BAC120 marker distribution across a frozen phylogeny panel.

Task 06 froze a panel in ``manifests/pilot.tsv`` and handed on a bounded follow-up: most tips had
no BAC120 marker QC.  Task 07's own ``gtdbtk identify`` step produces those counts as a
by-product, so this joins the identify output back onto the frozen panel rather than running a
separate screen.

It reports.  It never drops a tip: the downstream completion criterion is that the inferred tip set
matches the frozen panel exactly, so a low-marker tip is *flagged before inference* and its
disposition is stated, not silently removed.  A marker count is not a completeness estimate.

The identify genome ids may carry a prefix (GTDB-Tk refuses ids that collide with its own reference
genomes), so ``--id-map`` maps them back to the frozen tip ids.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_id_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    rows = read_tsv(path)
    required = {"gtdbtk_genome_id", "tip_id"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} needs columns {sorted(required)}")
    mapping = {r["gtdbtk_genome_id"]: r["tip_id"] for r in rows}
    if len(mapping) != len(rows):
        raise ValueError(f"duplicate gtdbtk_genome_id in {path}")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(f"duplicate tip_id in {path}")
    return mapping


def parse_prior_count(status: str) -> int | None:
    """``screened:118_unique_bac120`` -> 118; ``not_screened`` -> None."""
    if not status or not status.startswith("screened:"):
        return None
    token = status.split(":", 1)[1].split("_", 1)[0]
    return int(token)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", required=True, type=Path,
                        help="frozen panel manifest (manifests/pilot.tsv); read-only")
    parser.add_argument("--markers-summary", required=True, type=Path,
                        help="gtdbtk.bac120.markers_summary.tsv from this run's identify")
    parser.add_argument("--translation-tables", type=Path, default=None,
                        help="gtdbtk.translation_table_summary.tsv (headerless, two columns)")
    parser.add_argument("--failed-genomes", type=Path, default=None,
                        help="gtdbtk.failed_genomes.tsv from this run's identify")
    parser.add_argument("--id-map", type=Path, default=None,
                        help="TSV with gtdbtk_genome_id and tip_id columns")
    parser.add_argument("--followup", type=Path, default=None,
                        help="Task 06 tables/marker-qc-followup.tsv, to confirm it is fully absorbed")
    parser.add_argument("--flag-below", type=int, default=100,
                        help="flag a tip with fewer than this many unique markers (default 100)")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)

    panel = read_tsv(args.panel)
    if not panel:
        raise SystemExit(f"{args.panel} is empty")
    panel_by_tip = {r["tip_id"]: r for r in panel}
    if len(panel_by_tip) != len(panel):
        raise SystemExit(f"duplicate tip_id in {args.panel}")

    id_map = read_id_map(args.id_map)
    markers = read_tsv(args.markers_summary)
    counts: dict[str, dict] = {}
    for row in markers:
        raw = row["name"]
        tip = id_map.get(raw, raw)
        if tip in counts:
            raise SystemExit(f"two identify rows map to tip {tip}")
        counts[tip] = row

    tables = {}
    if args.translation_tables and args.translation_tables.exists():
        for line in args.translation_tables.read_text().splitlines():
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise SystemExit(f"malformed translation-table row: {line!r}")
            tables[id_map.get(fields[0], fields[0])] = fields[1]

    failed = []
    if args.failed_genomes and args.failed_genomes.exists():
        for line in args.failed_genomes.read_text().splitlines():
            if line.strip():
                failed.append(line.split("\t")[0])

    missing = sorted(set(panel_by_tip) - set(counts))
    extra = sorted(set(counts) - set(panel_by_tip))

    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for tip, prow in panel_by_tip.items():
        mrow = counts.get(tip)
        prior = parse_prior_count(prow.get("marker_qc_status", ""))
        unique = int(mrow["number_unique_genes"]) if mrow else None
        rows.append({
            "tip_id": tip,
            "panel_scope": prow.get("panel_scope", ""),
            "panel_roles": prow.get("panel_roles", ""),
            "ingroup_status": prow.get("ingroup_status", ""),
            "formal_names": prow.get("formal_names", ""),
            "prior_marker_qc_status": prow.get("marker_qc_status", ""),
            "prior_unique_bac120": "" if prior is None else prior,
            "number_unique_genes": "" if unique is None else unique,
            "number_multiple_genes": mrow["number_multiple_genes"] if mrow else "",
            "number_multiple_unique_genes": mrow["number_multiple_unique_genes"] if mrow else "",
            "number_missing_genes": mrow["number_missing_genes"] if mrow else "",
            "translation_table": tables.get(tip, ""),
            "identify_failed": tip in failed,
            "newly_screened": prior is None and unique is not None,
            "agrees_with_prior": (
                "" if prior is None or unique is None else str(prior == unique)
            ),
            "flag": (
                "identify_failed" if tip in failed
                else "no_identify_result" if unique is None
                else f"below_{args.flag_below}_unique_markers" if unique < args.flag_below
                else ""
            ),
        })
    rows.sort(key=lambda r: (r["number_unique_genes"] == "", r["number_unique_genes"], r["tip_id"]))

    with (args.output / "panel-marker-qc.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    screened = [r for r in rows if r["number_unique_genes"] != ""]
    distribution = collections.Counter(int(r["number_unique_genes"]) for r in screened)
    with (args.output / "marker-count-distribution.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["unique_bac120_markers", "tips"])
        for k in sorted(distribution):
            writer.writerow([k, distribution[k]])

    followup_absorbed = None
    if args.followup and args.followup.exists():
        wanted = {r["tip_id"] for r in read_tsv(args.followup)}
        resolved = {r["tip_id"] for r in rows
                    if r["tip_id"] in wanted and r["number_unique_genes"] != ""}
        followup_absorbed = {
            "handed_over_by_task_06": len(wanted),
            "now_screened": len(resolved),
            "still_unscreened": sorted(wanted - resolved),
        }

    disagreements = [r for r in rows if r["agrees_with_prior"] == "False"]
    flagged = [r for r in rows if r["flag"]]
    summary = {
        "panel_tips": len(panel_by_tip),
        "tips_with_identify_result": len(screened),
        "tips_without_identify_result": missing,
        "identify_rows_not_in_panel": extra,
        "identify_failed_genomes": failed,
        "newly_screened": sum(1 for r in rows if r["newly_screened"]),
        "previously_screened": sum(1 for r in rows if r["prior_unique_bac120"] != ""),
        "unique_marker_range": (
            [min(distribution), max(distribution)] if distribution else []
        ),
        "median_unique_markers": (
            sorted(int(r["number_unique_genes"]) for r in screened)[len(screened) // 2]
            if screened else None
        ),
        "distribution": {str(k): distribution[k] for k in sorted(distribution)},
        "flag_threshold": args.flag_below,
        "flagged_tips": [
            {"tip_id": r["tip_id"], "unique": r["number_unique_genes"], "flag": r["flag"],
             "formal_names": r["formal_names"], "ingroup_status": r["ingroup_status"]}
            for r in flagged
        ],
        "prior_screen_disagreements": [
            {"tip_id": r["tip_id"], "prior": r["prior_unique_bac120"],
             "now": r["number_unique_genes"]} for r in disagreements
        ],
        "translation_tables": dict(collections.Counter(
            r["translation_table"] for r in rows if r["translation_table"])),
        "task_06_followup": followup_absorbed,
        "note": (
            "A unique-marker count is not a completeness estimate and no tip is dropped here. "
            "Counts come from this run's own gtdbtk identify, so the Task 06 follow-up is "
            "absorbed rather than screened separately."
        ),
    }
    (args.output / "panel-marker-qc-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "distribution"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
