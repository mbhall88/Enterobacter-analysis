#!/usr/bin/env python3
"""Build the two Task 06 builder inputs that change when the panel goes from pilot to full cohort.

``scripts/freeze_phylogeny_inputs.py`` is the audited place where a tip set is decided, and Task 11
must not replace it.  Two of its inputs are panel-specific, though, and this script produces them so
that the change from a 388-tip pilot to a full cohort is a *documented derivation* rather than a
hand-edited table:

``--out/full-study-membership.tsv``
    The ``--study-draft`` input.  For the pilot that was an 80-genome deterministic sampling draft.
    Under **SET-17** there is no representative cap, so full coverage means *every* study candidate
    is a member and the selection rule is one sentence rather than a sampling design.

``--out/full-dispositions.tsv``
    The ``--dispositions`` input: the excludes, the restricted-context tip and the declared
    out-of-genus tips.  The settled exclusions and restrictions are **carried forward verbatim**
    from the pilot's dispositions, because they are settled decisions (SET-01, SET-02, SET-03,
    SET-13) and re-deriving a settled decision is how wording drifts.  The **out-of-genus set is
    recomputed**, because it is a property of the panel and the pilot's count is a pilot count.

The out-of-genus rule, recomputed
---------------------------------
A tip is ``retain_declared_out_of_genus_context`` when either:

**(a)** GTDB R232 places it in a genus other than the ingroup genus.  ``g__Enterobacter_D`` counts
as *other*: it is a GTDB genus split sister to *Pseudescherichia*, not the genus under study.  This
is a counting-and-rooting declaration, never a nomenclatural act, and the archived names are
preserved unchanged (SET-13).

**(b)** it is absent from GTDB R232 and its best type-anchor ANI is below **every** genome GTDB does
place in the ingroup genus.  The boundary is therefore read off the data rather than set by hand:
genomes deposited after the release must not be treated as in-genus merely because one release
happened not to contain them.  The size of the gap between the two classes is reported, so a reader
can see whether the rule is sharp on this cohort or merely happens to separate it.

Every rule carries a count, every declared tip carries its numbers, and the full ledger of the
evidence for all candidates is written out, not only for the ones the rule moves.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

CARRIED_DISPOSITIONS = ("exclude", "retain_restricted_context")
OUT_OF_GENUS = "retain_declared_out_of_genus_context"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def classify_ingroup(study_accessions, gtdb, ani, ingroup_genus, excluded):
    """Return (ledger rows, out-of-genus accessions, boundary record).

    ``excluded`` accessions still get a ledger row — their evidence is why they were excluded in
    several cases — but they cannot influence the in-genus minimum, because a genome that is not a
    tip must not move a boundary that decides the status of tips.
    """
    ledger = []
    in_genus_best = []
    for accession in sorted(study_accessions):
        placement = gtdb.get(accession, {})
        measures = ani.get(accession, {})
        genus = placement.get("gtdb_genus", "")
        best = number(measures.get("best_ani_percent", ""))
        if genus == ingroup_genus and accession not in excluded and best is not None:
            in_genus_best.append((best, accession))
    if not in_genus_best:
        raise SystemExit("no genome is placed in the ingroup genus; refusing to derive a boundary")
    boundary_value, boundary_accession = min(in_genus_best)

    out_of_genus: list[str] = []
    for accession in sorted(study_accessions):
        placement = gtdb.get(accession, {})
        measures = ani.get(accession, {})
        genus = placement.get("gtdb_genus", "")
        best = number(measures.get("best_ani_percent", ""))
        af = number(measures.get("best_af_min_percent", ""))
        if genus and genus != ingroup_genus:
            rule, declared = "gtdb_genus_outside_the_ingroup_genus", True
        elif genus == ingroup_genus:
            rule, declared = "gtdb_places_it_in_the_ingroup_genus", False
        elif best is None:
            raise SystemExit(
                f"{accession}: absent from GTDB and carries no best type-anchor ANI, so neither "
                f"rule can be applied. Refusing to guess its genus.")
        elif best < boundary_value:
            rule, declared = "absent_from_gtdb_and_below_every_in_genus_best_type_ani", True
        else:
            rule, declared = "absent_from_gtdb_but_at_or_above_the_in_genus_minimum", False
        if declared and accession not in excluded:
            out_of_genus.append(accession)
        ledger.append({
            "accession_version": accession,
            "gtdb_genus": genus,
            "gtdb_species": placement.get("gtdb_species", ""),
            "gtdb_evidence_level": placement.get("evidence_level", ""),
            "best_type_accession_version": measures.get("best_type_accession_version", ""),
            "best_type_taxon": measures.get("best_type_taxon", ""),
            "best_ani_percent": measures.get("best_ani_percent", ""),
            "best_af_min_percent": measures.get("best_af_min_percent", ""),
            "rule_applied": rule,
            "declared_out_of_genus": str(declared),
            "excluded_from_panel": str(accession in excluded),
            "margin_to_in_genus_minimum_ani": (
                "" if best is None else f"{best - boundary_value:+.2f}"),
            "af_min_percent_note": "" if af is None else f"{af:.2f}",
        })

    declared_best = [number(row["best_ani_percent"]) for row in ledger
                     if row["declared_out_of_genus"] == "True"
                     and row["excluded_from_panel"] == "False"
                     and number(row["best_ani_percent"]) is not None]
    boundary = {
        "ingroup_genus": ingroup_genus,
        "in_genus_minimum_best_type_ani_percent": boundary_value,
        "in_genus_minimum_set_by": boundary_accession,
        "declared_out_of_genus_maximum_best_type_ani_percent": (
            max(declared_best) if declared_best else None),
        "gap_percent_points": (
            round(boundary_value - max(declared_best), 4) if declared_best else None),
        "rule_is_sharp_on_this_cohort": bool(declared_best) and max(declared_best) < boundary_value,
    }
    return ledger, out_of_genus, boundary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study-manifest", type=Path, required=True,
                        help="manifests/study-candidates.tsv")
    parser.add_argument("--gtdb-membership", type=Path, required=True,
                        help="Task 05 gtdb-membership.tsv")
    parser.add_argument("--type-ani", type=Path, required=True,
                        help="Task 04 type-assignment-summary.tsv")
    parser.add_argument("--carry-dispositions", type=Path, required=True,
                        help="the pilot dispositions whose settled excludes/restrictions carry over")
    parser.add_argument("--ingroup-genus", default="Enterobacter")
    parser.add_argument("--membership-reason", required=True,
                        help="the one sentence that justifies full coverage, with its record")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    study = read_tsv(args.study_manifest)
    accessions = [row["accession_version"] for row in study]
    if len(set(accessions)) != len(accessions):
        raise SystemExit("the study manifest repeats an accession; refusing to build a cohort")
    gtdb = {row["accession_version"]: row for row in read_tsv(args.gtdb_membership)
            if row.get("role") == "study"}
    ani = {row["study_accession_version"]: row for row in read_tsv(args.type_ani)}

    carried = [row for row in read_tsv(args.carry_dispositions)
               if row["disposition"] in CARRIED_DISPOSITIONS]
    carried_present = [row for row in carried if row["accession_version"] in set(accessions)]
    carried_absent = [row["accession_version"] for row in carried
                      if row["accession_version"] not in set(accessions)]
    excluded = {row["accession_version"] for row in carried_present
                if row["disposition"] == "exclude"}

    ledger, out_of_genus, boundary = classify_ingroup(
        set(accessions), gtdb, ani, args.ingroup_genus, excluded)

    membership = [{"accession_version": row["accession_version"],
                   "genome": row["accession_version"],
                   "selection_reason": args.membership_reason,
                   "revised_species": row.get("legacy_revised_species", ""),
                   "revised_subspecies": row.get("legacy_revised_subspecies", "")}
                  for row in study]

    dispositions = list(carried_present)
    carried_ids = {row["accession_version"] for row in carried_present}
    pilot_out = {row["accession_version"] for row in read_tsv(args.carry_dispositions)
                 if row["disposition"] == OUT_OF_GENUS}
    for accession in out_of_genus:
        if accession in carried_ids:
            raise SystemExit(f"{accession}: both carried forward and recomputed as out-of-genus")
        entry = next(row for row in ledger if row["accession_version"] == accession)
        if entry["rule_applied"] == "gtdb_genus_outside_the_ingroup_genus":
            reason = (
                f"GTDB R232 places it in g__{entry['gtdb_genus']} "
                f"(s__{entry['gtdb_species']}, {entry['gtdb_evidence_level']}), outside "
                f"g__{args.ingroup_genus}. Task 04 agrees independently: its best match among the "
                f"type anchors is {entry['best_type_taxon']} at {entry['best_ani_percent']}% ANI "
                f"and {entry['af_min_percent_note']}% minimum aligned fraction, "
                f"{entry['margin_to_in_genus_minimum_ani']} points against the lowest best-type ANI "
                f"of any genome GTDB does place in the genus. Retained as a tip because removing it "
                f"would delete the misidentification evidence the study exists to report, but "
                f"declared out of genus: not counted in any {args.ingroup_genus} ingroup total and "
                f"not used to judge the root. Archived names unchanged; this is a "
                f"counting-and-rooting declaration, not a nomenclatural act.")
        else:
            reason = (
                f"Absent from GTDB R232, so no release-based genus is available, but its best match "
                f"among the type anchors is {entry['best_type_taxon']} at "
                f"{entry['best_ani_percent']}% ANI and {entry['af_min_percent_note']}% minimum "
                f"aligned fraction — {entry['margin_to_in_genus_minimum_ani']} points against the "
                f"lowest best-type ANI of any genome GTDB does place in g__{args.ingroup_genus}, so "
                f"it is below every in-genus genome. Given the same treatment as the genomes GTDB "
                f"could place, so the out-of-genus set is defined by evidence rather than by what "
                f"one release happened to contain.")
        dispositions.append({
            "accession_version": accession,
            "disposition": OUT_OF_GENUS,
            "disposition_reason": reason,
            "decision_record": (
                "SET-12 item 4; SET-13 (declared out-of-genus context); recomputed at full cohort "
                "scale by scripts/build_full_cohort_inputs.py from the Task 05 GTDB R232 membership "
                "table and the Task 04 type-anchored ANI summary"),
        })
    dispositions.sort(key=lambda row: (row["disposition"], row["accession_version"]))

    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "full-study-membership.tsv", membership,
              ["accession_version", "genome", "selection_reason", "revised_species",
               "revised_subspecies"])
    write_tsv(args.output / "full-dispositions.tsv", dispositions,
              ["accession_version", "disposition", "disposition_reason", "decision_record"])
    write_tsv(args.output / "ingroup-scope-ledger.tsv", ledger, list(ledger[0]))

    recomputed = set(out_of_genus)
    record = {
        "ingroup_genus": args.ingroup_genus,
        "study_candidates": len(study),
        "membership_rule": args.membership_reason,
        "carried_dispositions_applied": len(carried_present),
        "carried_dispositions_not_in_this_manifest": sorted(carried_absent),
        "excluded": sorted(excluded),
        "declared_out_of_genus": sorted(recomputed),
        "counts_by_rule": dict(sorted(Counter(row["rule_applied"] for row in ledger).items())),
        "counts_by_disposition": dict(sorted(
            Counter(row["disposition"] for row in dispositions).items())),
        "boundary": boundary,
        "comparison_with_the_pilot": {
            "pilot_declared_out_of_genus": sorted(pilot_out),
            "recomputed_declared_out_of_genus": sorted(recomputed),
            "in_both": sorted(pilot_out & recomputed),
            "new_at_full_scale": sorted(recomputed - pilot_out),
            "in_the_pilot_but_not_recomputed": sorted(pilot_out - recomputed),
        },
    }
    (args.output / "full-cohort-inputs.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
