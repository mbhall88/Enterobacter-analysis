"""Select type-anchored ANI boundary cases and the pilot genomes they require.

The input is a skani ``dist`` table of study chromosomes against the Task 02 type anchors. The
output is an auditable chain: an expected-pair ledger, a per-genome type assignment table, an
enumerated boundary-case table, a threshold-sensitivity grid, and the pilot additions that follow
from them.

Three things this module refuses to do, because each one has already produced a wrong number in
this project's archived analyses:

* It never treats an absent pair as ANI zero. A pair that skani did not report is carried through
  the ledger as ``not_measured`` and is excluded from every denominator it would otherwise distort.
* It never compares an aligned fraction against a threshold without naming the unit. skani reports
  both aligned fractions as percentages, and every column this module writes is suffixed
  ``_percent`` so a later reader cannot silently mix a fraction with a percentage.
* It never reports a single threshold outcome as the answer. The 95% ANI and 80% aligned-fraction
  rules are the project's legacy operational baselines and are under review (DEC-007), so the
  sensitivity grid reports how the counts move, including the strict ``>`` versus inclusive ``>=``
  boundary that the Task 12 audit found changes 40 historical pairs.

Selection is driven by the evidence, not by a desired tree: every rule is applied to every study
genome, and the reason each genome enters the pilot is written next to it.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import difflib
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# The competing groups the execution plan makes mandatory anchors for Task 04. Membership is by
# species binomial, so a subspecies type joins its species' group.
DISPUTED_GROUPS = {
    "hormaechei_complex": [
        "Enterobacter hormaechei", "Enterobacter hoffmannii", "Enterobacter intestinihominis",
        "Enterobacter quasihormaechei",
    ],
    "cloacae_dissolvens": ["Enterobacter cloacae"],
    "roggenkampii_quasiroggenkampii": [
        "Enterobacter roggenkampii", "Enterobacter quasiroggenkampii",
        "Enterobacter pseudoroggenkampii",
    ],
    "asburiae_dykesii": ["Enterobacter asburiae", "Enterobacter dykesii"],
}


def species_binomial(taxon_name: str) -> str:
    """'Enterobacter hormaechei subsp. hoffmannii' -> 'Enterobacter hormaechei'."""
    match = re.match(r"^(\S+\s+\S+)", taxon_name.strip())
    return match.group(1) if match else taxon_name.strip()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], header: list[str]) -> None:
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if row.get(c) is None else str(row.get(c, "")) for c in header) + "\n")


# ---------------------------------------------------------------------------------------------
# Legacy label resolution
#
# The archived study labels are working notes, not nomenclature. They contain orthographic
# variants ("Enterobacter intestinhominis" for E. intestinihominis, "Enterobacter Kobei"),
# operational placeholders ("Novel Species 3"), GTDB placeholders ("Enterobacter sp001888805"),
# uncertainty markers ("subsp?") and compound labels naming two epithets at once. Comparing those
# strings directly against a type anchor's name manufactures discordance that is really spelling,
# and hides discordance that is really taxonomy. Every label is therefore resolved by an ordered
# set of named rules, and the rule that fired is written next to the result.
# ---------------------------------------------------------------------------------------------

PLACEHOLDER_OPERATIONAL = re.compile(r"^novel\s+species\s+\d+$", re.IGNORECASE)
PLACEHOLDER_GENOMIC = re.compile(r"^\S+\s+sp\d+$", re.IGNORECASE)
UNCERTAINTY_SUFFIX = re.compile(r"\s+(subsp\.?\??\d*|sub\d+|subsp\.?\s*\d+)$", re.IGNORECASE)


@dataclass
class ResolvedLabel:
    raw: str
    resolved_name: str
    rule: str
    similarity: float = 0.0

    @property
    def species(self) -> str:
        return species_binomial(self.resolved_name) if self.resolved_name else ""

    @property
    def is_nomenclatural(self) -> bool:
        return bool(self.resolved_name)


class LabelResolver:
    """Resolve an archived working label onto a reviewed nomenclatural name, or refuse to."""

    def __init__(self, nomenclatural_names: set[str], minimum_similarity: float = 0.95):
        self.names = {name.lower(): name for name in nomenclatural_names}
        self.epithets = {
            name.lower().split()[1] for name in nomenclatural_names if len(name.split()) > 1}
        self.minimum_similarity = minimum_similarity

    def resolve(self, raw: str) -> ResolvedLabel:
        label = " ".join((raw or "").split())
        if not label:
            return ResolvedLabel(raw=raw, resolved_name="", rule="absent")
        lowered = label.lower()
        if lowered in self.names:
            return ResolvedLabel(raw=label, resolved_name=self.names[lowered],
                                 rule="exact_case_insensitive", similarity=1.0)
        if PLACEHOLDER_OPERATIONAL.match(label):
            return ResolvedLabel(raw=label, resolved_name="", rule="operational_placeholder")
        if PLACEHOLDER_GENOMIC.match(label):
            return ResolvedLabel(raw=label, resolved_name="", rule="genomic_placeholder")
        named_epithets = [token for token in lowered.split()[1:] if token in self.epithets]
        if len(named_epithets) > 1:
            return ResolvedLabel(raw=label, resolved_name="", rule="compound_label_two_epithets")
        stripped = UNCERTAINTY_SUFFIX.sub("", label).strip()
        if stripped.lower() != lowered and stripped.lower() in self.names:
            return ResolvedLabel(raw=label, resolved_name=self.names[stripped.lower()],
                                 rule="uncertainty_suffix_removed", similarity=1.0)
        candidates = difflib.get_close_matches(
            stripped.lower(), list(self.names), n=2, cutoff=self.minimum_similarity)
        if len(candidates) == 1:
            similarity = difflib.SequenceMatcher(None, stripped.lower(), candidates[0]).ratio()
            return ResolvedLabel(raw=label, resolved_name=self.names[candidates[0]],
                                 rule="orthographic_variant", similarity=round(similarity, 4))
        if len(candidates) > 1:
            return ResolvedLabel(raw=label, resolved_name="", rule="ambiguous_near_match")
        return ResolvedLabel(raw=label, resolved_name="", rule="unresolved_legacy_label")


def load_nomenclatural_names(crosswalk: Path) -> set[str]:
    return {row["exact_name"].strip() for row in read_tsv(crosswalk)
            if row["rank"] in ("species", "subspecies") and row["exact_name"].strip()}


@dataclass
class TypeAnchor:
    accession: str
    taxon_names: list[str] = field(default_factory=list)
    ranks: list[str] = field(default_factory=list)
    chromosome_path: str = ""
    whole_path: str = ""
    scope_status: str = ""
    chromosome_records: int = 0

    @property
    def species(self) -> str:
        return species_binomial(self.taxon_names[0]) if self.taxon_names else ""

    @property
    def label(self) -> str:
        return " / ".join(self.taxon_names)

    @property
    def disputed_group(self) -> str:
        for group, members in DISPUTED_GROUPS.items():
            if self.species in members:
                return group
        return ""


@dataclass
class Measurement:
    ani: float
    af_type_percent: float
    af_study_percent: float
    type_contigs: int
    study_contigs: int

    @property
    def af_min_percent(self) -> float:
        return min(self.af_type_percent, self.af_study_percent)


def load_type_anchors(types_manifest: Path, scope_table: Path) -> dict[str, TypeAnchor]:
    anchors: dict[str, TypeAnchor] = {}
    for row in read_tsv(types_manifest):
        accession = row["selected_accession_version"].strip()
        if not accession:
            continue
        anchor = anchors.setdefault(accession, TypeAnchor(accession=accession))
        anchor.taxon_names.append(row["taxon_name"])
        anchor.ranks.append(row["rank"])
    for row in read_tsv(scope_table):
        anchor = anchors.get(row["accession_version"])
        if anchor is None:
            raise SystemExit(f"{row['accession_version']} has sequence scopes but no types.tsv row")
        anchor.chromosome_path = row["chromosome_scope_path"]
        anchor.whole_path = row["whole_path"]
        anchor.scope_status = row["chromosome_scope_status"]
        anchor.chromosome_records = int(row["chromosome_records"])
    missing = [a for a, anchor in anchors.items() if not anchor.chromosome_path]
    if missing:
        raise SystemExit(f"type anchors without a chromosome-scope FASTA: {missing}")
    # Sort names so the species row precedes its subspecies row, for a stable label.
    for anchor in anchors.values():
        anchor.taxon_names.sort(key=len)
    return anchors


def load_measurements(table: Path, study_by_path: dict[str, str],
                      type_by_path: dict[str, str]) -> dict[tuple[str, str], Measurement]:
    """Read a skani dist table keyed by (study accession, type accession).

    skani reports ANI and both aligned fractions as percentages in the range 0-100.
    """
    measurements: dict[tuple[str, str], Measurement] = {}
    for row in read_tsv(table):
        study = study_by_path.get(row["Query_file"]) or study_by_path.get(row["Ref_file"])
        type_accession = type_by_path.get(row["Ref_file"]) or type_by_path.get(row["Query_file"])
        if study is None or type_accession is None:
            raise SystemExit(f"unrecognised file pair in {table}: "
                             f"{row['Query_file']} / {row['Ref_file']}")
        study_is_query = row["Query_file"] in study_by_path
        af_study = float(row["Align_fraction_query"] if study_is_query else row["Align_fraction_ref"])
        af_type = float(row["Align_fraction_ref"] if study_is_query else row["Align_fraction_query"])
        study_contigs = int(row["Num_query_contigs"] if study_is_query else row["Num_ref_contigs"])
        type_contigs = int(row["Num_ref_contigs"] if study_is_query else row["Num_query_contigs"])
        key = (study, type_accession)
        if key in measurements:
            raise SystemExit(f"{table}: duplicate measurement for {key}")
        measurements[key] = Measurement(
            ani=float(row["ANI"]), af_type_percent=af_type, af_study_percent=af_study,
            type_contigs=type_contigs, study_contigs=study_contigs)
    return measurements


def passes(measurement: Measurement, ani_threshold: float, af_threshold: float,
           inclusive: bool) -> bool:
    """Both-genome aligned-fraction rule, with the >= / > boundary made explicit."""
    if inclusive:
        return (measurement.ani >= ani_threshold
                and measurement.af_type_percent >= af_threshold
                and measurement.af_study_percent >= af_threshold)
    return (measurement.ani > ani_threshold
            and measurement.af_type_percent > af_threshold
            and measurement.af_study_percent > af_threshold)


def rank_types(study: str, anchors: dict[str, TypeAnchor],
               measurements: dict[tuple[str, str], Measurement]) -> list[tuple[str, Measurement]]:
    scored = [(accession, measurements[(study, accession)])
              for accession in anchors if (study, accession) in measurements]
    scored.sort(key=lambda item: (-item[1].ani, item[0]))
    return scored



# ---------------------------------------------------------------------------------------------
# Pilot tiers
#
# Every flagged genome is reported, but not every flag is a reason to spend a pilot tip. A tier
# says what the tree is being asked to do:
#
# tier_1_ani_cannot_decide
#     ANI has returned an answer that ANI itself cannot arbitrate: no anchor above threshold, a
#     value on the threshold, two named species inseparable within the tightest margin, a call
#     that the aligned-fraction filter or the replicon scope flips, or a case the coordinator has
#     already marked for review. A tree is the only remaining arbiter.
# tier_2_boundary_context
#     The case sits near a rule boundary without being decided by it. These genomes make the
#     boundary visible rather than resolving it.
# tier_3_group_representation
#     A disputed-group member or a label conflict whose ANI position is otherwise unambiguous.
#     Including all of them would mean most of the cohort, so representation is by a stated
#     deterministic rule: within each (best type species, resolved label) stratum, the genomes
#     with the smallest margin to a competing species, with the accession as the tiebreak. The
#     rule is fixed before the topology is seen and does not consult a preferred grouping.
# ---------------------------------------------------------------------------------------------

TIER_1_FLAGS = {
    "no_type_at_or_above_ani_threshold",
    "best_type_ani_near_species_threshold",
    "aligned_fraction_filter_changes_assignment",
    "strict_versus_inclusive_boundary_changes_outcome",
    "replicon_scope_changes_best_type",
    "mandatory_legacy_label_review",
    "best_pair_missing_in_reciprocal_direction",
    "reciprocal_direction_ani_disagreement",
    "no_type_measured",
}
TIER_2_FLAGS = {
    "nearest_competing_species_within_margin",
    "aligned_fraction_near_filter_boundary",
    "replicon_scope_changes_rule_outcome",
}


def assign_tier(flags: list[str], margin: float | None, tier_one_margin: float) -> str:
    if any(flag in TIER_1_FLAGS for flag in flags):
        return "tier_1_ani_cannot_decide"
    if ("nearest_competing_species_within_margin" in flags
            and margin is not None and margin <= tier_one_margin):
        return "tier_1_ani_cannot_decide"
    if any(flag in TIER_2_FLAGS or flag.startswith("legacy_label_not_nomenclatural")
           for flag in flags):
        return "tier_2_boundary_context"
    return "tier_3_group_representation"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ani", type=Path, required=True,
                        help="skani dist table, study chromosomes versus type chromosome scope")
    parser.add_argument("--reciprocal", type=Path,
                        help="the same comparison run in the opposite query/reference direction")
    parser.add_argument("--alternate-scope", type=Path,
                        help="study chromosomes versus whole type assemblies")
    parser.add_argument("--study-index", type=Path, required=True)
    parser.add_argument("--type-scopes", type=Path, required=True)
    parser.add_argument("--types", type=Path, required=True)
    parser.add_argument("--name-crosswalk", type=Path, required=True,
                        help="manifests/name-alias-crosswalk.tsv: the reviewed nomenclatural names")
    parser.add_argument("--study-manifest", type=Path, required=True)
    parser.add_argument("--review-cases", type=Path, required=True,
                        help="one study accession per line: the unresolved legacy label cases")
    parser.add_argument("--ani-threshold", type=float, default=95.0)
    parser.add_argument("--af-threshold", type=float, default=80.0)
    parser.add_argument("--ani-window", type=float, default=1.0,
                        help="half-width in ANI points of the near-threshold band")
    parser.add_argument("--af-window", type=float, default=5.0,
                        help="half-width in percentage points of the near-filter band")
    parser.add_argument("--tier-one-margin", type=float, default=0.1,
                        help="ANI points within which two named species are treated as "
                             "inseparable by ANI alone")
    parser.add_argument("--representatives-per-stratum", type=int, default=2,
                        help="tier 3 genomes kept per (best type species, resolved label) stratum")
    parser.add_argument("--margin", type=float, default=0.5,
                        help="ANI points within which a competing species counts as ambiguous")
    parser.add_argument("--direction-tolerance", type=float, default=0.1,
                        help="ANI points of query/reference asymmetry treated as agreement")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    anchors = load_type_anchors(args.types, args.type_scopes)
    resolver = LabelResolver(load_nomenclatural_names(args.name_crosswalk))
    study_rows = read_tsv(args.study_index)
    study_by_path = {row["fasta_path"]: row["accession_version"] for row in study_rows}
    study_order = [row["accession_version"] for row in study_rows]
    study_meta = {row["accession_version"]: row for row in study_rows}
    type_chromosome_by_path = {a.chromosome_path: a.accession for a in anchors.values()}
    type_whole_by_path = {a.whole_path: a.accession for a in anchors.values()}

    labels = {row["accession_version"]: row for row in read_tsv(args.study_manifest)}
    review_cases = {line.strip() for line in args.review_cases.read_text().splitlines() if line.strip()}
    unknown_review = review_cases - set(study_order)

    measurements = load_measurements(args.ani, study_by_path, type_chromosome_by_path)
    reciprocal = (load_measurements(args.reciprocal, study_by_path, type_chromosome_by_path)
                  if args.reciprocal else {})
    alternate = (load_measurements(args.alternate_scope, study_by_path, type_whole_by_path)
                 if args.alternate_scope else {})

    # ---- expected-pair ledger -------------------------------------------------------------
    ledger = []
    for study, type_accession in itertools.product(study_order, sorted(anchors)):
        measurement = measurements.get((study, type_accession))
        ledger.append({
            "study_accession_version": study,
            "type_accession_version": type_accession,
            "type_taxon": anchors[type_accession].label,
            "measurement_status": "measured" if measurement else "not_measured",
            "ani_percent": f"{measurement.ani:.4f}" if measurement else "",
            "af_type_percent": f"{measurement.af_type_percent:.4f}" if measurement else "",
            "af_study_percent": f"{measurement.af_study_percent:.4f}" if measurement else "",
            "reciprocal_status": ("measured" if (study, type_accession) in reciprocal
                                  else ("not_measured" if reciprocal else "not_run")),
            "note": "" if measurement else "absent from skani output; not an ANI of zero",
        })
    write_tsv(args.output / "expected-pair-ledger.tsv", ledger, list(ledger[0]))

    # ---- per-genome assignment ------------------------------------------------------------
    assignments = []
    flags_by_study: dict[str, list[str]] = defaultdict(list)
    endpoints_by_study: dict[str, set[str]] = defaultdict(set)

    for study in study_order:
        ranked = rank_types(study, anchors, measurements)
        row: dict[str, object] = {
            "study_accession_version": study,
            "entity_role": study_meta[study]["entity_role"],
            "usage_restrictions": study_meta[study]["usage_restrictions"],
            "types_measured": len(ranked),
            "types_expected": len(anchors),
        }
        label_row = labels.get(study, {})
        legacy = resolver.resolve(label_row.get("legacy_revised_species", ""))
        ncbi = resolver.resolve(label_row.get("historical_ncbi", ""))
        row["legacy_revised_species"] = legacy.raw
        row["legacy_label_resolution_rule"] = legacy.rule
        row["legacy_label_resolved_name"] = legacy.resolved_name
        row["legacy_label_match_similarity"] = legacy.similarity or ""
        row["legacy_revised_subspecies"] = (label_row.get("legacy_revised_subspecies") or "").strip()
        row["historical_ncbi_label"] = ncbi.raw
        row["historical_ncbi_resolution_rule"] = ncbi.rule
        row["historical_ncbi_resolved_name"] = ncbi.resolved_name

        if not ranked:
            row.update({"best_type_accession_version": "", "assignment_status": "no_type_measured"})
            assignments.append(row)
            flags_by_study[study].append("no_type_measured")
            continue

        best_accession, best = ranked[0]
        best_anchor = anchors[best_accession]
        row.update({
            "best_type_accession_version": best_accession,
            "best_type_taxon": best_anchor.label,
            "best_type_species": best_anchor.species,
            "best_type_scope_status": best_anchor.scope_status,
            "best_type_chromosome_records": best_anchor.chromosome_records,
            "best_ani_percent": f"{best.ani:.4f}",
            "best_af_type_percent": f"{best.af_type_percent:.4f}",
            "best_af_study_percent": f"{best.af_study_percent:.4f}",
            "best_af_min_percent": f"{best.af_min_percent:.4f}",
        })

        competing = next(((a, m) for a, m in ranked[1:]
                          if anchors[a].species != best_anchor.species), None)
        if competing:
            competing_accession, competing_measurement = competing
            row.update({
                "nearest_competing_type_accession_version": competing_accession,
                "nearest_competing_type_taxon": anchors[competing_accession].label,
                "nearest_competing_species": anchors[competing_accession].species,
                "nearest_competing_ani_percent": f"{competing_measurement.ani:.4f}",
                "ani_margin_to_competing_species": f"{best.ani - competing_measurement.ani:.4f}",
            })

        second = ranked[1] if len(ranked) > 1 else None
        if second:
            row["second_best_type_accession_version"] = second[0]
            row["second_best_ani_percent"] = f"{second[1].ani:.4f}"

        inclusive_pass = passes(best, args.ani_threshold, args.af_threshold, inclusive=True)
        strict_pass = passes(best, args.ani_threshold, args.af_threshold, inclusive=False)
        row["best_pair_passes_inclusive_rule"] = str(inclusive_pass)
        row["best_pair_passes_strict_rule"] = str(strict_pass)
        row["assignment_status"] = (
            "type_supported" if inclusive_pass
            else ("ani_supported_af_below_filter"
                  if best.ani >= args.ani_threshold else "no_type_at_or_above_ani_threshold"))

        # ---- boundary rules ----
        flags = []
        if abs(best.ani - args.ani_threshold) <= args.ani_window:
            flags.append("best_type_ani_near_species_threshold")
        if best.ani < args.ani_threshold:
            flags.append("no_type_at_or_above_ani_threshold")
        if competing and (best.ani - competing[1].ani) <= args.margin:
            flags.append("nearest_competing_species_within_margin")
        if best.ani >= args.ani_threshold and not inclusive_pass:
            flags.append("aligned_fraction_filter_changes_assignment")
        if abs(best.af_min_percent - args.af_threshold) <= args.af_window:
            flags.append("aligned_fraction_near_filter_boundary")
        if inclusive_pass != strict_pass:
            flags.append("strict_versus_inclusive_boundary_changes_outcome")
        if best_anchor.disputed_group:
            flags.append(f"disputed_group:{best_anchor.disputed_group}")
        if study in review_cases:
            flags.append("mandatory_legacy_label_review")
        # A label only contradicts a type anchor once it has been resolved to a nomenclatural
        # name. An unresolved label is its own finding, not evidence of a taxonomic conflict.
        if legacy.rule not in ("absent", "exact_case_insensitive", "orthographic_variant",
                               "uncertainty_suffix_removed"):
            flags.append(f"legacy_label_not_nomenclatural:{legacy.rule}")
        elif legacy.is_nomenclatural and best.ani >= args.ani_threshold:
            if legacy.species != best_anchor.species:
                flags.append("legacy_label_discordant_with_best_type")
        if (ncbi.is_nomenclatural and best.ani >= args.ani_threshold
                and ncbi.species.startswith("Enterobacter ")
                and ncbi.species != best_anchor.species):
            flags.append("ncbi_label_discordant_with_best_type")

        if reciprocal:
            back = reciprocal.get((study, best_accession))
            if back is None:
                flags.append("best_pair_missing_in_reciprocal_direction")
            elif abs(back.ani - best.ani) > args.direction_tolerance:
                flags.append("reciprocal_direction_ani_disagreement")
                row["reciprocal_best_ani_percent"] = f"{back.ani:.4f}"
                row["reciprocal_ani_difference"] = f"{back.ani - best.ani:.4f}"

        if alternate:
            alternate_ranked = sorted(
                ((a, alternate[(study, a)]) for a in anchors if (study, a) in alternate),
                key=lambda item: (-item[1].ani, item[0]))
            if alternate_ranked:
                alternate_best, alternate_measurement = alternate_ranked[0]
                row["whole_scope_best_type_accession_version"] = alternate_best
                row["whole_scope_best_ani_percent"] = f"{alternate_measurement.ani:.4f}"
                row["whole_scope_best_af_study_percent"] = (
                    f"{alternate_measurement.af_study_percent:.4f}")
                if alternate_best != best_accession:
                    flags.append("replicon_scope_changes_best_type")
                elif (passes(alternate_measurement, args.ani_threshold, args.af_threshold,
                             inclusive=True) != inclusive_pass):
                    flags.append("replicon_scope_changes_rule_outcome")

        row["boundary_flags"] = ";".join(flags)
        assignments.append(row)
        flags_by_study[study] = flags
        if flags:
            endpoints_by_study[study].add(best_accession)
            if competing:
                endpoints_by_study[study].add(competing[0])

    header = list(dict.fromkeys(key for row in assignments for key in row))
    write_tsv(args.output / "type-assignment-summary.tsv", assignments, header)

    # ---- boundary cases -------------------------------------------------------------------
    assignment_by_study = {row["study_accession_version"]: row for row in assignments}
    cases = []
    for study in study_order:
        flags = flags_by_study.get(study) or []
        if not flags:
            continue
        row = assignment_by_study[study]
        cases.append({
            "study_accession_version": study,
            "priority_tier": assign_tier(
                flags,
                float(row["ani_margin_to_competing_species"])
                if row.get("ani_margin_to_competing_species") else None,
                args.tier_one_margin),
            "boundary_flags": ";".join(flags),
            "flag_count": len(flags),
            "assignment_status": row.get("assignment_status", ""),
            "best_type_accession_version": row.get("best_type_accession_version", ""),
            "best_type_taxon": row.get("best_type_taxon", ""),
            "best_ani_percent": row.get("best_ani_percent", ""),
            "best_af_min_percent": row.get("best_af_min_percent", ""),
            "nearest_competing_type_accession_version":
                row.get("nearest_competing_type_accession_version", ""),
            "nearest_competing_type_taxon": row.get("nearest_competing_type_taxon", ""),
            "nearest_competing_ani_percent": row.get("nearest_competing_ani_percent", ""),
            "ani_margin_to_competing_species": row.get("ani_margin_to_competing_species", ""),
            "legacy_revised_species": row.get("legacy_revised_species", ""),
            "historical_ncbi_label": row.get("historical_ncbi_label", ""),
            "entity_role": row.get("entity_role", ""),
            "usage_restrictions": row.get("usage_restrictions", ""),
        })
    cases.sort(key=lambda row: (-row["flag_count"], row["study_accession_version"]))
    write_tsv(args.output / "boundary-cases.tsv", cases, list(cases[0]) if cases else
              ["study_accession_version", "boundary_flags"])

    # ---- estimator behaviour: direction and replicon scope ---------------------------------
    # skani's two query/reference directions are the same estimator run twice, so they are a
    # diagnostic of estimator and coverage behaviour, not two independent observations.
    direction_rows = []
    for (study, type_accession), forward in sorted(measurements.items()):
        back = reciprocal.get((study, type_accession))
        direction_rows.append({
            "study_accession_version": study,
            "type_accession_version": type_accession,
            "ani_study_as_query_percent": f"{forward.ani:.4f}",
            "ani_type_as_query_percent": f"{back.ani:.4f}" if back else "",
            "ani_difference": f"{back.ani - forward.ani:.4f}" if back else "",
            "af_study_study_as_query_percent": f"{forward.af_study_percent:.4f}",
            "af_study_type_as_query_percent": f"{back.af_study_percent:.4f}" if back else "",
            "status": "measured_both_directions" if back else "measured_one_direction_only",
        })
    if reciprocal:
        write_tsv(args.output / "direction-comparison.tsv", direction_rows, list(direction_rows[0]))

    scope_rows = []
    for (study, type_accession), chromosome in sorted(measurements.items()):
        whole = alternate.get((study, type_accession))
        scope_rows.append({
            "study_accession_version": study,
            "type_accession_version": type_accession,
            "type_chromosome_scope_status": anchors[type_accession].scope_status,
            "ani_type_chromosome_scope_percent": f"{chromosome.ani:.4f}",
            "ani_type_whole_scope_percent": f"{whole.ani:.4f}" if whole else "",
            "ani_difference": f"{whole.ani - chromosome.ani:.4f}" if whole else "",
            "af_type_chromosome_scope_percent": f"{chromosome.af_type_percent:.4f}",
            "af_type_whole_scope_percent": f"{whole.af_type_percent:.4f}" if whole else "",
            "af_type_difference": f"{whole.af_type_percent - chromosome.af_type_percent:.4f}"
                                  if whole else "",
        })
    if alternate:
        write_tsv(args.output / "replicon-scope-comparison.tsv", scope_rows, list(scope_rows[0]))

    # ---- threshold sensitivity ------------------------------------------------------------
    sensitivity = []
    for ani_threshold in (94.0, 95.0, 95.5, 96.0, 97.0):
        for af_threshold in (0.0, 50.0, 80.0, 90.0):
            for inclusive in (True, False):
                supported = sum(
                    1 for study in study_order
                    for ranked in [rank_types(study, anchors, measurements)]
                    if ranked and passes(ranked[0][1], ani_threshold, af_threshold, inclusive))
                pairs = sum(1 for measurement in measurements.values()
                            if passes(measurement, ani_threshold, af_threshold, inclusive))
                sensitivity.append({
                    "ani_threshold_percent": ani_threshold,
                    "af_threshold_percent": af_threshold,
                    "boundary_operator": ">=" if inclusive else ">",
                    "study_genomes_with_supported_best_type": supported,
                    "study_genomes_without": len(study_order) - supported,
                    "study_type_pairs_passing": pairs,
                })
    write_tsv(args.output / "threshold-sensitivity.tsv", sensitivity, list(sensitivity[0]))

    # ---- pilot additions ------------------------------------------------------------------
    tier_by_study: dict[str, str] = {}
    margin_by_study: dict[str, float | None] = {}
    for study in study_order:
        flags = flags_by_study.get(study) or []
        if not flags:
            continue
        raw_margin = assignment_by_study[study].get("ani_margin_to_competing_species")
        margin = float(raw_margin) if raw_margin else None
        margin_by_study[study] = margin
        tier_by_study[study] = assign_tier(flags, margin, args.tier_one_margin)

    # Deterministic tier-3 representation: smallest margin first, accession as the tiebreak.
    strata: dict[tuple[str, str], list[str]] = defaultdict(list)
    for study, tier in tier_by_study.items():
        if tier != "tier_3_group_representation":
            continue
        row = assignment_by_study[study]
        key = (row.get("best_type_species", ""),
               row.get("legacy_label_resolved_name") or row.get("legacy_revised_species", ""))
        strata[key].append(study)
    tier_three_selected: set[str] = set()
    for key, members in strata.items():
        members.sort(key=lambda s: (margin_by_study.get(s) if margin_by_study.get(s) is not None
                                    else float("inf"), s))
        tier_three_selected.update(members[:args.representatives_per_stratum])

    additions = []
    for study in study_order:
        flags = flags_by_study.get(study) or []
        if not flags:
            continue
        tier = tier_by_study[study]
        if tier == "tier_3_group_representation" and study not in tier_three_selected:
            continue
        row = assignment_by_study[study]
        restrictions = study_meta[study]["usage_restrictions"]
        additions.append({
            "priority_tier": tier,
            "stratum_key": "|".join((row.get("best_type_species", ""),
                                     row.get("legacy_label_resolved_name")
                                     or row.get("legacy_revised_species", ""))),
            "entity_id": study,
            "entity_kind": "study_genome",
            "accession_version": study,
            "taxon_or_label": row.get("legacy_revised_species") or row.get("historical_ncbi_label", ""),
            "replicon_scope": "chromosome_only_archived_or_resolved",
            "fasta_path": study_meta[study]["fasta_path"],
            "addition_reason": ";".join(flags),
            "paired_endpoints": ";".join(sorted(endpoints_by_study[study])),
            "best_type_accession_version": row.get("best_type_accession_version", ""),
            "best_ani_percent": row.get("best_ani_percent", ""),
            "best_af_min_percent": row.get("best_af_min_percent", ""),
            "ani_margin_to_competing_species": row.get("ani_margin_to_competing_species", ""),
            "pilot_eligibility": ("restricted_pending_placement" if restrictions
                                  else "eligible_pilot_tip"),
            "eligibility_note": restrictions,
        })
    # Only genomes actually kept as pilot additions drag their type endpoints in with them.
    retained_studies = {row["entity_id"] for row in additions}
    implicated = sorted({accession for study, endpoints in endpoints_by_study.items()
                         if study in retained_studies for accession in endpoints})
    for accession in implicated:
        anchor = anchors[accession]
        additions.append({
            "priority_tier": "tier_1_ani_cannot_decide",
            "stratum_key": anchor.species,
            "entity_id": accession,
            "entity_kind": "type_reference",
            "accession_version": accession,
            "taxon_or_label": anchor.label,
            "replicon_scope": f"chromosome_scope:{anchor.scope_status}",
            "fasta_path": anchor.chromosome_path,
            "addition_reason": "boundary_pair_endpoint:type_anchor",
            "paired_endpoints": ";".join(
                sorted(study for study, endpoints in endpoints_by_study.items()
                       if accession in endpoints and study in retained_studies)),
            "best_type_accession_version": "",
            "best_ani_percent": "",
            "best_af_min_percent": "",
            "ani_margin_to_competing_species": "",
            "pilot_eligibility": "eligible_pilot_tip",
            "eligibility_note": (
                f"{anchor.chromosome_records} chromosome-assigned records; "
                "fragmentation affects aligned fraction"),
        })
    write_tsv(args.output / "pilot-boundary-additions.tsv", additions, list(additions[0]))

    summary = {
        "study_genomes": len(study_order),
        "type_anchors": len(anchors),
        "expected_pairs": len(study_order) * len(anchors),
        "measured_pairs": len(measurements),
        "unmeasured_pairs": len(study_order) * len(anchors) - len(measurements),
        "reciprocal_pairs": len(reciprocal),
        "alternate_scope_pairs": len(alternate),
        "rules": {
            "ani_threshold_percent": args.ani_threshold,
            "af_threshold_percent": args.af_threshold,
            "af_rule": "both genomes",
            "ani_window_percent": args.ani_window,
            "af_window_percent": args.af_window,
            "competing_species_margin_percent": args.margin,
            "direction_tolerance_percent": args.direction_tolerance,
            "status": "legacy operational baselines under review in DEC-007; not an approved policy",
        },
        "study_genomes_with_any_boundary_flag": sum(1 for f in flags_by_study.values() if f),
        "flag_counts": {
            flag: sum(1 for flags in flags_by_study.values() if flag in flags)
            for flag in sorted({flag for flags in flags_by_study.values() for flag in flags})
        },
        "tier_counts": {
            tier: sum(1 for value in tier_by_study.values() if value == tier)
            for tier in sorted(set(tier_by_study.values()))
        },
        "tier_three_strata": len(strata),
        "tier_three_represented": len(tier_three_selected),
        "pilot_additions": {
            "study_genomes": sum(1 for row in additions if row["entity_kind"] == "study_genome"),
            "type_references": sum(1 for row in additions if row["entity_kind"] == "type_reference"),
            "restricted_pending_placement": sum(
                1 for row in additions if row["pilot_eligibility"] == "restricted_pending_placement"),
        },
        "direction_agreement": {
            "pairs_compared": sum(1 for r in direction_rows if r["ani_difference"]),
            "max_absolute_ani_difference": max(
                (abs(float(r["ani_difference"])) for r in direction_rows if r["ani_difference"]),
                default=None),
            "pairs_over_direction_tolerance": sum(
                1 for r in direction_rows
                if r["ani_difference"] and abs(float(r["ani_difference"])) > args.direction_tolerance),
        },
        "replicon_scope_effect": {
            "pairs_compared": sum(1 for r in scope_rows if r["ani_difference"]),
            "max_absolute_ani_difference": max(
                (abs(float(r["ani_difference"])) for r in scope_rows if r["ani_difference"]),
                default=None),
            "max_type_af_drop_from_including_plasmids": min(
                (float(r["af_type_difference"]) for r in scope_rows if r["af_type_difference"]),
                default=None),
        },
        "review_cases_supplied": len(review_cases),
        "review_cases_not_in_study_manifest": sorted(unknown_review),
    }
    (args.output / "selection-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
