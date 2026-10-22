#!/usr/bin/env python3
"""Freeze one common tip set for BAC120 ML, Verticall/FastME and matched Mashtree.

The three methods must receive the *same* sequences under the *same* tip labels, or a topology
difference between them cannot be attributed to the method.  This script is the single place where
that set is decided, and it is deliberately unforgiving: nothing enters the panel without a recorded
reason and nothing leaves it without one either.

What it merges
--------------
``--study-draft``      the deterministic 80-genome study sampling draft (Task 06's earlier stage)
``--additions``        the tiered ANI boundary additions from Task 04
``--types``            the approved nomenclatural type anchors from Task 02
``--outgroups``        the outgroup assemblies selected under SET-12
``--dispositions``     the explicit keep/exclude decisions this freeze applies, each with a reason

Four things this script refuses to do silently
----------------------------------------------
**An addition without a reason.**  Every accession entering the panel carries at least one non-empty
reason string from the manifest that proposed it.  A blank reason is an error, not an empty cell.

**An exclusion without a reason.**  Exclusions come only from ``--dispositions``; each needs a
non-empty reason and a decision record (``SET-01``, ``SET-03``, a run path, an issue comment).  An
exclusion naming an accession that is not in the union is an error too, because that is how a typo
in an accession quietly fails to exclude anything.

**Two tips for one sequence.**  Task 01 established that four identical-assembly pairs are one
inference unit each.  Any panel accession that the exclusions manifest records as an alternate
submission is folded onto its selected inference entity, and after that the multiset of sequences in
every source FASTA is hashed and collisions are fatal.  Identical sequence must not become two tips
carrying different names and voting twice.

**Two tips for one accession, or one tip for two names.**  A species type is also the type of its
autonymic subspecies, so one accession can anchor two names (*E. cloacae* and *E. cloacae* subsp.
*cloacae*; *E. hormaechei* and *E. hormaechei* subsp. *hormaechei*).  That is one tip carrying two
formal names, never two tips.  Formal nomenclatural names and operational/database labels are kept
in separate columns and never merged into one "species" field.

Sequence scope and source precedence
------------------------------------
Every tip is chromosome scope (DEC-007/SET-10).  Where an accession is both a study genome and an
approved type anchor, the verified NCBI-derived chromosome record wins over the archived study copy:
Task 04 found one archived type copy truncated, and a type anchor must be the verified record.  The
chosen source and the reason for it are recorded per tip.

Outputs
-------
``manifests/<prefix>.tsv``       one row per tip, with roles, names, reasons and hashes
``manifests/<prefix>.list``      one derivative FASTA path per line
``manifests/<prefix>.batch.tsv`` ``<derivative path>\\t<tip id>``, the GTDB-Tk batchfile form
``manifests/<prefix>.ids.txt``   one tip id per line
run directory                    the same tables plus ``freeze.json`` with every input hash

``--prefix`` defaults to ``pilot``.  Task 11 passes ``--prefix full`` together with a study draft
listing every study candidate rather than the 80-genome sampling draft, so the full cohort is
frozen by this same audited code path rather than by a second implementation.  Only names change.

``pilot.list`` and ``pilot.batch.tsv`` name the *derivative* files that
``scripts/repair_pilot_inputs.py`` writes, because those are what the three methods consume.  The
source path and hash stay in ``pilot.tsv`` so the mapping is never lost.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ACCESSION_RE = re.compile(r"GC[AF]_\d{9}\.\d+")

TIER_ORDER = ["tier_1_ani_cannot_decide", "tier_2_boundary_context", "tier_3_group_representation"]

# Dispositions this freeze understands.  ``exclude`` removes a tip; the others keep it and say how it
# may be read.  ``declared_out_of_genus_context`` tips are tips, but they are not Enterobacter for
# the purpose of any ingroup count and they cannot be used to judge the root.
DISPOSITIONS = {
    "retain_ingroup": "ingroup",
    "retain_declared_out_of_genus_context": "declared_out_of_genus_context",
    "retain_restricted_context": "restricted_context",
    "exclude": None,
}


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_multiset(path: Path) -> tuple[str, int, int]:
    """Order-independent hash of a FASTA's sequence content, ignoring headers and case.

    Two files that differ only in record order, line wrapping, header text or case describe the
    same sequence set and must collide here.  Returns (hash, record count, total bases).
    """
    records: list[str] = []
    current: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    records.append("".join(current))
                    current = []
            else:
                current.append(line.strip())
    if current:
        records.append("".join(current))
    if not records:
        raise SystemExit(f"{path}: no sequence records")
    normalised = sorted(record.upper() for record in records)
    digest = hashlib.sha256("\n".join(normalised).encode()).hexdigest()
    return digest, len(records), sum(len(record) for record in normalised)


class Tip:
    """One inference unit.  Roles, names and reasons accumulate; the accession never changes."""

    __slots__ = ("accession", "roles", "reasons", "tiers", "in_draft", "formal", "operational",
                 "folded_from")

    def __init__(self, accession: str) -> None:
        self.accession = accession
        self.roles: set[str] = set()
        self.reasons: set[str] = set()
        self.tiers: set[str] = set()
        self.in_draft = False
        self.formal: list[tuple[str, str]] = []       # (taxon name, rank)
        self.operational: set[str] = set()
        self.folded_from: set[str] = set()

    def add(self, role: str, reason: str, source: str) -> None:
        if not reason.strip():
            raise SystemExit(
                f"{self.accession}: {source} proposes this genome with an empty reason. "
                f"Every addition needs a reason; refusing to freeze an unexplained tip.")
        self.roles.add(role)
        self.reasons.add(reason.strip())

    def scope(self) -> str:
        """The smallest panel scope that already contains this tip."""
        if self.in_draft or "type_anchor" in self.roles or "outgroup" in self.roles:
            return "core"
        for tier in TIER_ORDER:
            if tier in self.tiers:
                return tier
        return "core"


def collect(args) -> tuple[dict[str, Tip], dict[str, str]]:
    """Build the candidate union.  Returns the tips and the alternate->selected fold map."""
    tips: dict[str, Tip] = {}

    def tip(accession: str) -> Tip:
        return tips.setdefault(accession, Tip(accession))

    for row in read_tsv(args.study_draft):
        entry = tip(row["accession_version"])
        entry.in_draft = True
        entry.add("study_genome", row.get("selection_reason", ""), "the study draft")
        for column in ("revised_species", "revised_subspecies"):
            if row.get(column):
                entry.operational.add(row[column])

    for row in read_tsv(args.additions):
        entry = tip(row["accession_version"])
        entry.tiers.add(row["priority_tier"])
        role = "type_anchor" if row["entity_kind"] == "type_reference" else "boundary_addition"
        entry.add(role, row.get("addition_reason", ""), f"the {row['priority_tier']} additions")
        if row.get("taxon_or_label"):
            entry.operational.add(row["taxon_or_label"])

    for row in read_tsv(args.types):
        if row["reference_role"] != "nomenclatural_type_anchor":
            continue
        accession = row["selected_accession_version"].strip()
        if not accession:
            continue
        entry = tip(accession)
        entry.add("type_anchor",
                  row.get("selection_basis", "") or row.get("selection_status", ""),
                  "the type anchors")
        # One accession can anchor two names: a species type is also its autonymic subspecies type.
        entry.formal.append((row["taxon_name"], row["rank"]))

    for row in read_tsv(args.outgroups):
        if not row["selection_status"].startswith("selected_outgroup"):
            continue
        entry = tip(row["accession_version"])
        entry.add("outgroup", row.get("selection_basis", ""), "the selected outgroups")
        entry.formal.append((row["taxon_name"], row["rank"]))

    # Fold duplicate submissions onto the inference unit Task 01 selected.
    fold: dict[str, str] = {}
    for row in read_tsv(args.exclusions):
        if row.get("entity_role") != "alternate_submission_record":
            continue
        alternate, selected = row["entity_id"], row.get("selected_inference_entity_id", "")
        if alternate in tips:
            if not selected:
                raise SystemExit(
                    f"{alternate}: recorded as an alternate submission with no selected inference "
                    f"entity. Refusing to guess which of an identical pair is the tip.")
            fold[alternate] = selected
            source = tips.pop(alternate)
            target = tips.setdefault(selected, Tip(selected))
            target.roles |= source.roles
            target.reasons |= source.reasons
            target.tiers |= source.tiers
            target.in_draft = target.in_draft or source.in_draft
            target.operational |= source.operational
            target.formal.extend(source.formal)
            target.folded_from.add(alternate)
            target.reasons.add(
                f"identical_submission_folded:{alternate}_is_one_inference_unit_with_{selected}")
    return tips, fold


def resolve_sources(tips: dict[str, Tip], args) -> dict[str, dict[str, str]]:
    """Choose one chromosome-scope FASTA per tip and say why that source won."""
    study_index = {row["accession_version"]: row for row in read_tsv(args.study_index)}
    type_scopes = {row["accession_version"]: row for row in read_tsv(args.type_scopes)}
    outgroup_scopes = {row["accession_version"]: row for row in read_tsv(args.outgroup_scopes)}

    sources: dict[str, dict[str, str]] = {}
    missing = []
    for accession, entry in tips.items():
        if "outgroup" in entry.roles and accession in outgroup_scopes:
            row = outgroup_scopes[accession]
            source, reason = row["chromosome_scope_path"], (
                "ncbi_chromosome_scope:selected_outgroup_verified_record")
        elif "type_anchor" in entry.roles and accession in type_scopes:
            row = type_scopes[accession]
            source, reason = row["chromosome_scope_path"], (
                "ncbi_chromosome_scope:type_anchor_uses_the_verified_ncbi_record_not_the_archive")
        elif accession in study_index:
            row = study_index[accession]
            source, reason = row["fasta_path"], (
                f"study_chromosome_index:{row.get('path_origin', 'archived_or_resolved')}")
        else:
            missing.append(accession)
            continue
        if not source:
            missing.append(accession)
            continue
        sources[accession] = {"path": source, "reason": reason,
                              "usage_restrictions": row.get("usage_restrictions", "")}
    if missing:
        raise SystemExit(
            f"no chromosome-scope sequence for {len(missing)} panel accessions: {sorted(missing)}")
    return sources


def apply_dispositions(tips: dict[str, Tip], args) -> tuple[dict[str, dict[str, str]], list[dict]]:
    rows = read_tsv(args.dispositions)
    by_accession: dict[str, dict[str, str]] = {}
    unknown = []
    for row in rows:
        accession = row["accession_version"]
        disposition = row["disposition"]
        if disposition not in DISPOSITIONS:
            raise SystemExit(f"{accession}: unknown disposition {disposition!r}; "
                             f"expected one of {sorted(DISPOSITIONS)}")
        if not row.get("disposition_reason", "").strip():
            raise SystemExit(f"{accession}: disposition {disposition} has no reason. "
                             f"Every exclusion and every declared status needs one.")
        if not row.get("decision_record", "").strip():
            raise SystemExit(f"{accession}: disposition {disposition} names no decision record.")
        if accession in by_accession:
            raise SystemExit(f"{accession}: listed twice in the dispositions")
        if accession not in tips:
            unknown.append(accession)
        by_accession[accession] = row
    if unknown:
        raise SystemExit(
            f"dispositions name {len(unknown)} accessions that are not in the candidate union: "
            f"{sorted(unknown)}. A disposition that matches nothing excludes nothing.")

    excluded = []
    for accession, row in by_accession.items():
        if row["disposition"] == "exclude":
            entry = tips.pop(accession)
            excluded.append({
                "accession_version": accession,
                "panel_roles": ";".join(sorted(entry.roles)),
                "priority_tiers": ";".join(sorted(entry.tiers)),
                "in_study_draft": str(entry.in_draft),
                "exclusion_reason": row["disposition_reason"],
                "decision_record": row["decision_record"],
            })
    return by_accession, sorted(excluded, key=lambda row: row["accession_version"])


def build_rows(tips, sources, dispositions, study_manifest, gtdb, marker_qc, derivative_dir, tiers):
    rows = []
    for accession in sorted(tips):
        entry = tips[accession]
        scope = entry.scope()
        if scope != "core" and scope not in tiers:
            continue
        source = sources[accession]
        path = Path(source["path"])
        multiset, records, bases = sequence_multiset(path)
        study = study_manifest.get(accession, {})
        placement = gtdb.get(accession, {})
        disposition = dispositions.get(accession, {})
        operational = set(entry.operational)
        for column in ("historical_ncbi", "historical_kraken_chromosome", "legacy_revised_species",
                       "legacy_revised_subspecies"):
            if study.get(column):
                operational.add(study[column])
        ingroup = DISPOSITIONS.get(disposition.get("disposition", "retain_ingroup"), "ingroup")
        if "outgroup" in entry.roles:
            ingroup = "outgroup"
        restrictions = [value for value in
                        (source.get("usage_restrictions", ""), study.get("usage_restrictions", ""))
                        if value]
        rows.append({
            "tip_id": accession,
            "accession_version": accession,
            "panel_scope": scope,
            "priority_tiers": ";".join(sorted(entry.tiers)),
            "panel_roles": ";".join(sorted(entry.roles)),
            "ingroup_status": ingroup,
            "in_study_draft": str(entry.in_draft),
            "formal_names": ";".join(name for name, _ in sorted(entry.formal)),
            "formal_name_ranks": ";".join(f"{name}={rank}" for name, rank in sorted(entry.formal)),
            "operational_labels": ";".join(sorted(operational)),
            "gtdb_species": placement.get("gtdb_species", ""),
            "gtdb_genus": placement.get("gtdb_genus", ""),
            "gtdb_evidence_level": placement.get("evidence_level", ""),
            "addition_reasons": " | ".join(sorted(entry.reasons)),
            "disposition": disposition.get("disposition", "retain_ingroup"),
            "disposition_reason": disposition.get("disposition_reason", "default_ingroup_member"),
            "decision_record": disposition.get("decision_record", ""),
            "folded_identical_submissions": ";".join(sorted(entry.folded_from)),
            "sequence_scope": "chromosome_only",
            "sequence_source": source["reason"],
            "source_fasta_path": str(path),
            "source_fasta_sha256": sha256(path),
            "source_sequence_multiset_sha256": multiset,
            "source_records": str(records),
            "source_bp": str(bases),
            "derivative_fasta_path": str(Path(derivative_dir) / f"{accession}.fna"),
            "marker_qc_status": marker_qc.get(accession, "not_screened"),
            "usage_restrictions": ";".join(
                dict.fromkeys(part for part in ";".join(restrictions).split(";") if part)),
        })

    seen_tip: dict[str, str] = {}
    seen_sequence: dict[str, str] = {}
    for row in rows:
        if row["tip_id"] in seen_tip:
            raise SystemExit(f"duplicate tip id {row['tip_id']}")
        seen_tip[row["tip_id"]] = row["accession_version"]
        digest = row["source_sequence_multiset_sha256"]
        if digest in seen_sequence:
            raise SystemExit(
                f"identical sequence in two tips: {seen_sequence[digest]} and {row['tip_id']} share "
                f"multiset hash {digest}. Identical sequence is one inference unit; resolve it in "
                f"the exclusions manifest or the dispositions before freezing.")
        seen_sequence[digest] = row["tip_id"]
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study-draft", type=Path, required=True)
    parser.add_argument("--additions", type=Path, required=True)
    parser.add_argument("--types", type=Path, required=True)
    parser.add_argument("--outgroups", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--dispositions", type=Path, required=True)
    parser.add_argument("--study-manifest", type=Path, required=True)
    parser.add_argument("--study-index", type=Path, required=True)
    parser.add_argument("--type-scopes", type=Path, required=True)
    parser.add_argument("--outgroup-scopes", type=Path, required=True)
    parser.add_argument("--gtdb-membership", type=Path)
    parser.add_argument("--marker-qc", type=Path, action="append", default=[],
                        help="BAC120 markers_summary.tsv already produced for some panel members")
    parser.add_argument("--tiers", default=",".join(TIER_ORDER),
                        help="comma-separated priority tiers to freeze; core members are always in")
    parser.add_argument("--derivative-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest-dir", type=Path, default=Path("manifests"))
    parser.add_argument("--prefix", default="pilot",
                        help="manifest basename stem: 'pilot' writes pilot.tsv/.list/.batch.tsv/"
                             ".ids.txt, 'full' writes full.*. The freeze logic is identical; only "
                             "the names change, so one audited builder serves both panels.")
    parser.add_argument("--write-manifests", action="store_true")
    args = parser.parse_args()

    if "/" in args.prefix or not args.prefix.strip():
        parser.error("--prefix is a filename stem, not a path")

    tiers = {tier.strip() for tier in args.tiers.split(",") if tier.strip()}
    unknown = tiers - set(TIER_ORDER)
    if unknown:
        parser.error(f"unknown tiers: {sorted(unknown)}")

    tips, fold = collect(args)
    candidate_count = len(tips)
    dispositions, excluded = apply_dispositions(tips, args)
    sources = resolve_sources(tips, args)

    study_manifest = {row["entity_id"]: row for row in read_tsv(args.study_manifest)}
    gtdb = {}
    if args.gtdb_membership:
        gtdb = {row["accession_version"]: row for row in read_tsv(args.gtdb_membership)
                if row.get("role") == "study"}
    marker_qc: dict[str, str] = {}
    for path in args.marker_qc:
        for row in read_tsv(path):
            match = ACCESSION_RE.search(row["name"])
            if not match:
                raise SystemExit(f"{path}: marker row {row['name']!r} has no versioned accession")
            marker_qc[match.group(0)] = f"screened:{row['number_unique_genes']}_unique_bac120"

    rows = build_rows(tips, sources, dispositions, study_manifest, gtdb, marker_qc,
                      args.derivative_dir, tiers)

    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.prefix
    fields = list(rows[0])
    write_tsv(args.output / f"{stem}.tsv", rows, fields)
    (args.output / f"{stem}.list").write_text(
        "\n".join(row["derivative_fasta_path"] for row in rows) + "\n")
    (args.output / f"{stem}.batch.tsv").write_text(
        "\n".join(f"{row['derivative_fasta_path']}\t{row['tip_id']}" for row in rows) + "\n")
    (args.output / f"{stem}.ids.txt").write_text("\n".join(row["tip_id"] for row in rows) + "\n")
    write_tsv(args.output / f"{stem}-exclusions.tsv", excluded,
              ["accession_version", "panel_roles", "priority_tiers", "in_study_draft",
               "exclusion_reason", "decision_record"])
    write_tsv(args.output / "source-derivative-map.tsv", rows,
              ["tip_id", "source_fasta_path", "source_fasta_sha256",
               "source_sequence_multiset_sha256", "derivative_fasta_path"])

    by_scope = defaultdict(int)
    by_status = defaultdict(int)
    by_role = defaultdict(int)
    for row in rows:
        by_scope[row["panel_scope"]] += 1
        by_status[row["ingroup_status"]] += 1
        for role in row["panel_roles"].split(";"):
            by_role[role] += 1
    inputs = [args.study_draft, args.additions, args.types, args.outgroups, args.exclusions,
              args.dispositions, args.study_manifest, args.study_index, args.type_scopes,
              args.outgroup_scopes, Path(__file__)]
    if args.gtdb_membership:
        inputs.append(args.gtdb_membership)
    record = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_prefix": args.prefix,
        "frozen_tiers": sorted(tiers),
        "candidate_union": candidate_count,
        "excluded": len(excluded),
        "tips": len(rows),
        "tips_by_panel_scope": dict(sorted(by_scope.items())),
        "tips_by_ingroup_status": dict(sorted(by_status.items())),
        "tips_by_role": dict(sorted(by_role.items())),
        "folded_identical_submissions": fold,
        "sequence_scope": "chromosome_only",
        "marker_qc_screened": sum(row["marker_qc_status"] != "not_screened" for row in rows),
        "marker_qc_not_screened": sum(row["marker_qc_status"] == "not_screened" for row in rows),
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "manifests_written": bool(args.write_manifests),
    }
    if args.write_manifests:
        args.manifest_dir.mkdir(parents=True, exist_ok=True)
        names = tuple(f"{stem}{suffix}" for suffix in
                      (".tsv", ".list", ".batch.tsv", ".ids.txt"))
        for name in names:
            (args.manifest_dir / name).write_bytes((args.output / name).read_bytes())
        record["manifest_sha256"] = {
            str(args.manifest_dir / name): sha256(args.manifest_dir / name) for name in names}
    (args.output / "freeze.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
