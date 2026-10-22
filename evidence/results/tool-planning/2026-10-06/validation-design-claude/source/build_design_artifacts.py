"""Metadata-only artifacts for the Task 15 validation design.

Reads the frozen split manifest, the remaining-reserve metadata, the
exposure-audit dispositions and the reference package concept metadata.
Writes the Layer 1 manifest (every study genome that is not a package
reference), the NCBI name crosswalk, Layer 3 strata counts, and two seeded
100-genome expert-sample options from the effective reserve: stratified
across NCBI names (chosen by Michael, 7 Oct 2026), and simple random. Also
writes the labelling comparator manifest (all study genomes, whole assemblies)
and the cohort tree-tip manifest (chromosome-only sequences of the paper
trees) used by the expert-labelling procedure. No genome sequence is read.
Outputs go to a new directory.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "src"))
from enterobacter_classifier.reporting import ReportingContext  # noqa: E402

SPLIT = ROOT / "results/tool-planning/2026-10-04/balanced-marker-discovery-astra/design-v2/split-manifest.tsv"
POOL = ROOT / "results/tool-planning/2026-10-06/prospective-confirmation-astra/inputs/execution-freeze/design/remaining-reserve-metadata.tsv"
ASSIGNMENTS = ROOT / ("results/ani-cases/2026-10-03/task12-astra-cp102-approved-final-reconciliation/"
                      "accession-assignments.tsv")
# Saved chromosome-to-chromosome skani pairs among the 1,090 paper-tree tips (Task 12 evidence; DEC-007 scope).
PAIRS = ROOT / "results/ani-cases/2026-09-23/task12-astra-evidence/skani/pairs-complete.tsv"
DISPOSITIONS = ROOT / ("results/tool-planning/2026-10-06/validation-design-claude/outputs/exposure-audit-attempt03/"
                       "reserve-exposure-dispositions.tsv")
PACKAGE = ROOT / "results/tool-planning/2026-10-06/portable-reference-package-astra/relocated-package/portable.json"
PACKAGE_SHA256 = "12aa0705b774de8ea262869fdc4c85d513c3aa97d39c674f08e85d1c9055b88f"
SEED = "Enterobacter-validation-expert-20261006-v1"
SIZES = (100,)  # Michael, 7 Oct 2026: 100
UNNAMED_SHARE = 0.25  # of the slots left after the fixed core; rounded half up
# Layer 1 genomes that shaped the coverage and residual rules (ani-scope-contract-astra, issue 15), and C45's split-block mate.
RULE_DESIGN = {"GCA_900322725.1": "rule_design_focal (C45; replay control)",
               "GCA_030406505.2": "rule_design_focal",
               "GCA_900322715.1": "split_block_of_rule_design_genome"}
# Crosswalk rows whose kind depends on a decision still open.
OPEN_DECISION = {"Enterobacter hormaechei subsp. hoffmannii": "DEC-004 A/B species rank (Leah)",
                 "Enterobacter intestinihominis": "DEC-004 A/B species rank (Leah)"}
csv.field_size_limit(100_000_000)

# NCBI name (strain suffix removed) -> (kind, expected project species IDs, expected subspecies IDs, basis).
# Kinds: equivalent | combination_change | broader | unnamed | outside_genus.
CROSSWALK = {
    "Enterobacter asburiae": ("equivalent", "ECTX-SP002", "", "same name"),
    "Enterobacter bugandensis": ("equivalent", "ECTX-SP003", "", "same name; nematophilus kept separate (issue 2, 27 Sep)"),
    "Enterobacter cancerogenus": ("equivalent", "ECTX-SP004", "", "same name"),
    "Enterobacter chengduensis": ("equivalent", "ECTX-SP005", "", "same species; NCBI gives no subspecies"),
    "Enterobacter chuandaensis": ("equivalent", "ECTX-SP007", "", "same name"),
    "Enterobacter cloacae": ("equivalent", "ECTX-SP008", "", "same species; NCBI gives no subspecies"),
    "Enterobacter cloacae subsp. cloacae": ("equivalent", "ECTX-SP008", "ECTX-SS003", "same name"),
    "Enterobacter cloacae subsp. dissolvens": ("equivalent", "ECTX-SP008", "ECTX-SS004", "same name"),
    "Enterobacter cloacae complex sp.": ("unnamed", "", "", "complex-level label; no species expectation"),
    "Enterobacter hormaechei": ("broader", "ECTX-SP009;ECTX-SP025", "",
                                "NCBI sensu lato covers project hormaechei and xiangfangensis (Group B and hoffmannii group "
                                "separated from the hormaechei type lineage, issue 2, 24 Sep)"),
    "Enterobacter hormaechei subsp. hormaechei": ("equivalent", "ECTX-SP009", "", "type lineage of project hormaechei"),
    "Enterobacter hormaechei subsp. hoffmannii": ("combination_change", "ECTX-SP025", "ECTX-SS007",
                                                  "project xiangfangensis subsp. hoffmannii (final name open)"),
    "Enterobacter hormaechei subsp. oharae": ("combination_change", "ECTX-SP025", "ECTX-SS008", "project xiangfangensis subsp. oharae"),
    "Enterobacter hormaechei subsp. steigerwaltii": ("combination_change", "ECTX-SP025", "ECTX-SS009",
                                                     "project xiangfangensis subsp. steigerwaltii"),
    "Enterobacter hormaechei subsp. xiangfangensis": ("combination_change", "ECTX-SP025", "ECTX-SS010",
                                                      "project xiangfangensis subsp. xiangfangensis"),
    "Enterobacter huaxiensis": ("equivalent", "ECTX-SP010", "", "same name"),
    "Enterobacter intestinihominis": ("combination_change", "ECTX-SP025", "ECTX-SS007",
                                      "project hoffmannii/intestinihominis group (issue 2, 24 Sep; final name open)"),
    "Enterobacter kobei": ("equivalent", "ECTX-SP011", "", "same name"),
    "Enterobacter ludwigii": ("equivalent", "ECTX-SP012", "", "same name"),
    "Enterobacter mori": ("equivalent", "ECTX-SP013", "", "same name"),
    "Enterobacter oligotrophicus": ("equivalent", "ECTX-SP015", "", "same name"),
    "Enterobacter pasteurii": ("equivalent", "ECTX-SP016", "", "same name"),
    "Enterobacter pseudoroggenkampii": ("equivalent", "ECTX-SP017", "", "same name"),
    "Enterobacter quasihormaechei": ("equivalent", "ECTX-SP018", "", "same name"),
    "Enterobacter quasiroggenkampii": ("combination_change", "ECTX-SP019", "ECTX-SS005", "project roggenkampii subsp. quasiroggenkampii"),
    "Enterobacter roggenkampii": ("equivalent", "ECTX-SP019", "", "same species; NCBI gives no subspecies"),
    "Enterobacter sichuanensis": ("equivalent", "ECTX-SP021", "", "same name"),
    "Enterobacter soli": ("equivalent", "ECTX-SP022", "", "same name"),
    "Enterobacter sp.": ("unnamed", "", "", "genus-level label; no species expectation"),
    "Enterobacter vonholyi": ("equivalent", "ECTX-SP023", "", "same name"),
    "Enterobacter wuhouensis": ("equivalent", "ECTX-SP024", "", "same name"),
    "Lelliottia nimipressuralis": ("outside_genus", "", "", "no Lelliottia reference; no Enterobacter result expected"),
    "Lelliottia sp.": ("outside_genus", "", "", "no Lelliottia reference; no Enterobacter result expected"),
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def half_up(x):
    return math.floor(x + 0.5)


def order(*parts):
    return hashlib.sha256("|".join((SEED,) + parts).encode()).hexdigest()


def normalise(name):
    """Drop strain designations; keep species and subspecies epithets and the 'sp.' forms."""
    name = name.strip()
    for prefix in ("Enterobacter cloacae complex sp.", "Enterobacter sp.", "Lelliottia sp."):
        if name.startswith(prefix):
            return prefix
    m = re.match(r"^(\S+ \S+)( subsp\. \S+)?", name)
    return m.group(0)


def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open("x", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite output")
    context = ReportingContext(PACKAGE, expected_sha256=PACKAGE_SHA256)
    refs = context.evidence_references()
    reference_accessions = {r["reference_id"].split(":", 1)[1] for r in refs}
    concepts = {r["species_id"]: r["species_name"] for r in refs} | {r["subspecies_id"]: r["subspecies_name"] for r in refs if r["subspecies_id"]}
    for kind, sp, ss, _ in CROSSWALK.values():
        for cid in filter(None, (sp.split(";") + [ss])):
            assert cid in concepts, cid

    with SPLIT.open() as handle:
        split = list(csv.DictReader(handle, delimiter="\t"))
    assert reference_accessions <= {r["accession"] for r in split}
    reserve = [r for r in split if r["partition"] == "internal_reserve"]
    assert len(reserve) == 52 and not reference_accessions & {r["accession"] for r in reserve}
    dev_blocks = {(k, r[k]) for r in split if r["partition"] == "development" for k in ("clone_block", "split_block")}
    for r in reserve:
        assert ("clone_block", r["clone_block"]) not in dev_blocks and ("split_block", r["split_block"]) not in dev_blocks
    ref_blocks = {(k, r[k]) for r in split if r["accession"] in reference_accessions for k in ("clone_block", "split_block")}
    # Layer 1 (Michael, 7 Oct 2026): every study genome that is not a package reference. The 52 internal-reserve
    # genomes were never used to fit a threshold and are reported as their own subset.
    layer1 = []
    for r in sorted((r for r in split if r["accession"] not in reference_accessions), key=lambda x: x["accession"]):
        stratum = ("context_or_unnamed_challenge" if r["named_target"] != "1" else
                   "named_high_QC" if r["qc_pass"] == "1" else "named_below_development_QC_gate")
        layer1.append(dict(accession=r["accession"], partition=r["partition"], stratum=stratum,
                           project_species=r["project_species"], project_subspecies=r["project_subspecies"],
                           quality_flags=r["quality_flags"], completeness=r["completeness"], contamination=r["contamination"],
                           assembly_level=r["assembly_level"], clone_block=r["clone_block"], split_block=r["split_block"],
                           reference_in_clone_block=("clone_block", r["clone_block"]) in ref_blocks,
                           reference_in_split_block=("split_block", r["split_block"]) in ref_blocks,
                           rule_design_role=RULE_DESIGN.get(r["accession"], ""), whole_fasta=r["whole_fasta"],
                           whole_sha256=r["whole_sha256"]))
    assert len(layer1) == len(split) - len(reference_accessions)
    with ASSIGNMENTS.open() as handle:
        assignments = {r["accession"]: r for r in csv.DictReader(handle, delimiter="\t")}
    for r in layer1:
        r["label_status"] = assignments[r["accession"]]["assignment_status"]
    # Descriptive Layer 1 grouping: which named species pass 95/80 against the other study genomes in the saved
    # chromosome pairs. Where the only passing species is the label, a wrong species result is near-impossible.
    named_label = {r["accession"]: r["project_species"] for r in split if r["named_target"] == "1"}
    passing = defaultdict(set)
    with PAIRS.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["status"] != "measured" or r["a"] == r["b"] or r["passes_95_80"] != "True":
                continue
            for x, y in ((r["a"], r["b"]), (r["b"], r["a"])):
                if y in named_label:
                    passing[x].add(named_label[y])
    for r in layer1:
        species = passing.get(r["accession"], set())
        r["saved_pairs_species_passing"] = ";".join(sorted(species))
        r["saved_pairs_ani_class"] = ("" if r["stratum"] == "context_or_unnamed_challenge" else
                                      "none" if not species else "multiple" if len(species) > 1 else
                                      "sole_species_is_label" if species == {r["project_species"]} else "sole_species_differs")

    # Labelling comparators: every study genome (named or not), whole assembly, as in the tool's own search.
    comparators = [dict(accession=r["accession"], label_kind="named" if r["named_target"] == "1" else "unassigned_or_context",
                        project_species=r["project_species"] if r["named_target"] == "1" else "",
                        project_subspecies=r["project_subspecies"] if r["named_target"] == "1" else "",
                        is_package_reference=r["accession"] in reference_accessions,
                        whole_fasta=r["whole_fasta"], whole_sha256=r["whole_sha256"])
                   for r in sorted(split, key=lambda x: x["accession"])]
    # Cohort tree tips: the chromosome-only sequences of the paper BAC120 and Verticall trees (Michael, 7 Oct 2026).
    tree_tips = []
    # A tip is labelled when it is a named study genome (the comparators' definition) or a type anchor outside the
    # study; "Enterobacter sp." lineages, context, outgroup and excluded tips are unassigned.
    for a, r in sorted(assignments.items()):
        if r["frozen_tree_sequence_scope"] != "chromosome_only":
            continue
        named = a in named_label or (r["record_role"] == "type_anchor_nonstudy" and r["project_species"]
                                     and not r["project_species"].endswith(" sp."))
        if a in named_label:
            assert r["project_species"] == named_label[a], a
        tree_tips.append(dict(accession=a, record_role=r["record_role"], assignment_status=r["assignment_status"],
                              tip_label_kind="named" if named else "unassigned_or_context",
                              project_species=r["project_species"] if named else "",
                              project_subspecies=r["project_subspecies"] if named else "",
                              frozen_tree_fasta=r["frozen_tree_fasta"], frozen_tree_fasta_sha256=r["frozen_tree_fasta_sha256"]))
    assert len(tree_tips) == 1089
    assert {r["accession"] for r in layer1 if r["rule_design_role"]} == set(RULE_DESIGN)

    with POOL.open() as handle:
        pool = list(csv.DictReader(handle, delimiter="\t"))
    assert len(pool) == 4857
    with DISPOSITIONS.open() as handle:
        disposition = {r["assembly_unit"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert set(disposition) == {r["assembly_unit"] for r in pool}
    for r in pool:
        r["disposition"] = disposition[r["assembly_unit"]]["disposition"]
        # Michael, 7 Oct 2026: near-clones of references are good examples, so every effective unit is eligible.
        r["sample_eligible"] = r["disposition"] == "effective_reserve"
        r["ncbi_name_normalised"] = normalise(r["ncbi_names"])
        assert r["ncbi_name_normalised"] in CROSSWALK, r["ncbi_names"]
        r["crosswalk_kind"] = CROSSWALK[r["ncbi_name_normalised"]][0]
    effective = [r for r in pool if r["disposition"] == "effective_reserve"]
    by_name_all, by_name = defaultdict(list), defaultdict(list)
    for r in pool:
        by_name_all[r["ncbi_name_normalised"]].append(r)
    for r in pool:
        if r["sample_eligible"]:
            by_name[r["ncbi_name_normalised"]].append(r)
    crosswalk = []
    for name, (kind, sp, ss, basis) in sorted(CROSSWALK.items()):
        members = [m for m in by_name_all.get(name, []) if m["disposition"] == "effective_reserve"]
        crosswalk.append(dict(ncbi_name_normalised=name, kind=kind, expected_species_ids=sp,
                              expected_species_names=";".join(concepts[c] for c in sp.split(";") if c),
                              expected_subspecies_id=ss, expected_subspecies_name=concepts.get(ss, ""),
                              effective_assemblies=len(members), effective_components=len({m["component_id"] for m in members}),
                              moderate_or_unreported=sum(m["qc_stratum"] in ("moderate_quality_stress", "unreported_quality_stress")
                                                         for m in members),
                              held_assemblies=sum(m["disposition"] != "effective_reserve" for m in by_name_all.get(name, [])),
                              open_decision=OPEN_DECISION.get(name, ""), basis=basis))
    strata = Counter((r["disposition"], r["crosswalk_kind"], r["qc_stratum"], r["assembly_level"]) for r in pool)
    strata_rows = [dict(disposition=k[0], crosswalk_kind=k[1], qc_stratum=k[2], assembly_level=k[3], assemblies=v)
                   for k, v in sorted(strata.items())]

    # Stratified expert-sample option, drawn from effective-reserve units only.
    # At most one assembly per linked component within a name group. Linked components are provenance
    # groups (shared BioProject/BioSample/strain links) that often span several species, so a global
    # one-per-component rule would exclude rare names; instead, components not yet used anywhere in the
    # sample are preferred, and the component is recorded so results can be checked for clustering.
    # Fixed core: both outside-genus records and one per rare named group (<=5 eligible assemblies).
    # The remaining slots: UNNAMED_SHARE (rounded half up) to unnamed records, the rest to common named
    # groups, each at least one, weighted by sqrt(linked components); largest remainder rounding.
    # Within a group given two or more slots, ceil(k/3) come from moderate or unreported quality where
    # available; single-slot groups take the first eligible record in seeded order.
    low = ("moderate_quality_stress", "unreported_quality_stress")
    rare = sorted(n for n, m in by_name.items() if len(m) <= 5 and CROSSWALK[n][0] != "outside_genus")
    outside = sorted(n for n in by_name if CROSSWALK[n][0] == "outside_genus")
    common = sorted(n for n in by_name if n not in rare and n not in outside and CROSSWALK[n][0] != "unnamed")
    members = {n: by_name[n] for n in common + rare + outside}
    members["unnamed"] = [r for n in by_name if CROSSWALK[n][0] == "unnamed" for r in by_name[n]]

    def allocate(total, weights):
        raw = {g: total * w / sum(weights.values()) for g, w in weights.items()}
        quota = {g: max(1, int(raw[g])) for g in weights}
        # Remainders are measured against the quota, so a group raised to its minimum of one is not raised again first.
        for g in sorted(weights, key=lambda g: (-(raw[g] - quota[g]), g)):
            if sum(quota.values()) >= total:
                break
            quota[g] += 1
        assert sum(quota.values()) == total, "too few slots for one per common group"
        return quota

    options, summary = {}, {}
    for size in SIZES:
        used, used_components, picks = set(), set(), []

        def take(group, k):
            ranked = sorted(members[group], key=lambda r: order(group, r["assembly_unit"]))
            ranked = [r for r in ranked if r["component_id"] not in used_components] + \
                     [r for r in ranked if r["component_id"] in used_components]
            want_low = math.ceil(k / 3) if k >= 2 else 0
            got = []
            for candidates, cap in (([r for r in ranked if r["qc_stratum"] in low], want_low), (ranked, k)):
                for r in candidates:
                    if len(got) >= cap:
                        break
                    if (group, r["component_id"]) not in used and r not in got:
                        got.append(r)
                        used.add((group, r["component_id"]))
                        used_components.add(r["component_id"])
            picks.extend(dict(r, sample_group=group) for r in got)
            return k - len(got)

        shortfall = sum(take(n, 1) for n in outside + rare)
        remaining = size - len(picks)
        unnamed_k = half_up(remaining * UNNAMED_SHARE)
        quota = allocate(remaining - unnamed_k, {g: math.sqrt(len({r["component_id"] for r in members[g]})) for g in common})
        shortfall += take("unnamed", unnamed_k)
        for g in common:
            shortfall += take(g, quota[g])
        assert shortfall == 0, "a group ran out of unlinked assemblies"
        assert all(p["sample_eligible"] for p in picks)
        rows = [dict(sample_size_option=size, assembly_unit=p["assembly_unit"], sample_group=p["sample_group"],
                     ncbi_names=p["ncbi_names"], crosswalk_kind=p["crosswalk_kind"], qc_stratum=p["qc_stratum"],
                     assembly_level=p["assembly_level"], component_id=p["component_id"],
                     selection_key=order(p["sample_group"], p["assembly_unit"])) for p in picks]
        assert len(rows) == size and len({r["assembly_unit"] for r in rows}) == len(rows)
        assert len({(r["sample_group"], r["component_id"]) for r in rows}) == len(rows)
        options[f"stratified-n{size}"] = rows

    # Simple random option: the first 100 eligible units in seeded order, with no stratification or component rule.
    eligible = sorted((r for r in pool if r["sample_eligible"]), key=lambda r: order("simple_random", r["assembly_unit"]))
    options["random-n100"] = [dict(sample_size_option=100, assembly_unit=r["assembly_unit"], sample_group="simple_random",
                                   ncbi_names=r["ncbi_names"], crosswalk_kind=r["crosswalk_kind"], qc_stratum=r["qc_stratum"],
                                   assembly_level=r["assembly_level"], component_id=r["component_id"],
                                   selection_key=order("simple_random", r["assembly_unit"])) for r in eligible[:100]]
    for name, rows in options.items():
        summary[name] = dict(selected=len(rows), name_groups=len({normalise(r["ncbi_names"]) for r in rows}),
                             project_species_expected=len({c for r in rows for c in CROSSWALK[normalise(r["ncbi_names"])][1].split(";") if c}),
                             moderate_or_unreported=sum(r["qc_stratum"] in low for r in rows),
                             distinct_components=len({r["component_id"] for r in rows}),
                             max_per_component=max(Counter(r["component_id"] for r in rows).values()),
                             by_kind=dict(Counter(r["crosswalk_kind"] for r in rows)),
                             by_name=dict(sorted(Counter(normalise(r["ncbi_names"]) for r in rows).items())))

    args.output.mkdir()
    write_tsv(args.output / "layer1-manifest.tsv", layer1)
    write_tsv(args.output / "labelling-comparators.tsv", comparators)
    write_tsv(args.output / "cohort-tree-tips.tsv", tree_tips)
    write_tsv(args.output / "ncbi-name-crosswalk.tsv", crosswalk)
    write_tsv(args.output / "layer3-pool-strata.tsv", strata_rows)
    for name, rows in options.items():
        write_tsv(args.output / f"expert-sample-option-{name}.tsv", rows)
    result = dict(
        layer1={part: dict(genomes=len(rows), strata=dict(Counter(r["stratum"] for r in rows)),
                           named_species=len({r["project_species"] for r in rows if r["stratum"] != "context_or_unnamed_challenge"}),
                           named_with_subspecies=sum(bool(r["project_subspecies"]) for r in rows
                                                     if r["stratum"] != "context_or_unnamed_challenge"),
                           clone_blocks=len({r["clone_block"] for r in rows if r["stratum"] != "context_or_unnamed_challenge"}),
                           split_blocks=len({r["split_block"] for r in rows if r["stratum"] != "context_or_unnamed_challenge"}),
                           named_sharing_clone_block_with_reference=sum(r["reference_in_clone_block"] for r in rows
                                                                        if r["stratum"] != "context_or_unnamed_challenge"),
                           named_sharing_split_block_with_reference=sum(r["reference_in_split_block"] for r in rows
                                                                        if r["stratum"] != "context_or_unnamed_challenge"))
                for part, rows in (("all_nonreference", layer1),
                                   ("development_nonreference", [r for r in layer1 if r["partition"] == "development"]),
                                   ("internal_reserve", [r for r in layer1 if r["partition"] == "internal_reserve"]))}
               | dict(references_among_layer1=0, reserve_development_blocks_shared=0, rule_design_genomes=sorted(RULE_DESIGN),
                      named_label_status=dict(Counter(r["label_status"] for r in layer1
                                                      if r["stratum"] != "context_or_unnamed_challenge")),
                      saved_pairs_ani_class={part: dict(Counter(r["saved_pairs_ani_class"] for r in layer1 if r["saved_pairs_ani_class"]
                                                                and (part == "all" or r["partition"] == part)))
                                             for part in ("all", "development", "internal_reserve")}),
        labelling_comparators=dict(genomes=len(comparators), by_kind=dict(Counter(r["label_kind"] for r in comparators))),
        cohort_tree_tips=dict(tips=len(tree_tips), by_kind=dict(Counter(r["tip_label_kind"] for r in tree_tips)),
                              by_role=dict(Counter(r["record_role"] for r in tree_tips)),
                              named_species=len({r["project_species"] for r in tree_tips if r["tip_label_kind"] == "named"}),
                              unassigned_by_status=dict(Counter(f'{r["record_role"]}|{r["assignment_status"]}|{assignments[r["accession"]]["project_species"]}'
                                                                for r in tree_tips if r["tip_label_kind"] != "named"))),
        layer3=dict(reserve_assemblies=len(pool), reserve_components=len({r["component_id"] for r in pool}),
                    held_assemblies=len(pool) - len(effective),
                    effective_assemblies=len(effective), effective_components=len({r["component_id"] for r in effective}),
                    sample_eligible=sum(r["sample_eligible"] for r in pool),
                    effective_by_kind=dict(Counter(r["crosswalk_kind"] for r in effective)),
                    effective_by_qc_stratum=dict(Counter(r["qc_stratum"] for r in effective)),
                    effective_by_level=dict(Counter(r["assembly_level"] for r in effective)),
                    effective_name_forms=len({r["ncbi_name_normalised"] for r in effective}),
                    largest_effective_components=[c for c, _ in Counter(r["component_id"] for r in effective).most_common(3)],
                    largest_effective_component_sizes=[n for _, n in Counter(r["component_id"] for r in effective).most_common(3)],
                    rare_name_groups=rare, outside_genus=outside, common_name_groups=common),
        expert_options={str(k): v for k, v in summary.items()}, seed=SEED,
        inputs={str(SPLIT): sha(SPLIT), str(POOL): sha(POOL), str(DISPOSITIONS): sha(DISPOSITIONS), str(ASSIGNMENTS): sha(ASSIGNMENTS),
                str(PAIRS): sha(PAIRS), str(PACKAGE): PACKAGE_SHA256},
        source_sha256=sha(__file__),
        sequence_files_read=0)
    (args.output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
