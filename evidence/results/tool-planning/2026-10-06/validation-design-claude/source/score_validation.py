"""Frozen scorer for the Task 15 validation design (section 3 of validation-design.md).

Scores `analyze` report.tsv rows against a truth table and writes a per-genome
table plus summary counts with exact two-sided 95% Clopper-Pearson intervals.
Written and tested before any validation result exists; its SHA256 is posted
at freeze and it is not edited afterwards.

Truth table (TSV), one row per genome:
  genome_id, layer, truth_kind, label_source, truth_species_id,
  truth_subspecies_id, ani_species_passing, subspecies_structural,
  expected_species_ids, expected_subspecies_id, block, component_id,
  query_sha256, subspecies_label_source (expert stratum: automatic,
  michael or empty), plus any grouping columns named with --group.
query_sha256 is the SHA256 of the FASTA the genome was run from; a report
whose query_sha256 differs is refused.
block is the relatedness unit for the block-level sensitivity: the clone
block for study genomes, the linked component for public genomes. Every
accuracy and crosswalk row must have one.
--group takes a column, or several joined by commas for a crossed grouping
(for example crosswalk_kind,qc_stratum,assembly_level).
label_source is one of:
  project_label         project study label (Layer 1)
  auto_single_species   expert stratum: automatic, one species passing 95/80
  auto_tree_resolved    expert stratum: automatic, the trees chose among
                        several species passing 95/80
  michael               expert stratum: labelled by Michael
  none                  Layer 3 rows
Accuracy rows (named, outside_genus, no_reference) need a label source other
than none; crosswalk rows must use none; descriptive rows (unresolved,
context) may carry any.
ani_species_passing (expert stratum only) lists the species IDs with a
labelling comparator at skani ANI >= 95 and both aligned fractions >= 80.
A named label equal to the only passing species is "structural": every
package reference is a comparator and the tool's species passes the same
test, so a wrong species result is impossible by construction. Structural
rows count for coverage, withheld reasons and subspecies outcomes but are
excluded from the wrong-result rate. Layer 1 rows are never marked
structural; Layer 1 is reported by its saved-pairs grouping instead.
For expert rows, a tool species outside ani_species_passing is impossible
when the labelling and tool skani runs match, so it is counted separately
whatever the label source.
subspecies_structural (expert stratum) marks a subspecies label, whoever
made it, equal to the subspecies of the species' untied highest-ANI
comparator when that comparator is a package reference: the tool's
subspecies can then only be the label or absent. Those rows are counted
apart from the wrong-subspecies proportion. truth_subspecies_id
"undetermined" (subspecies_label_source michael only, which includes
Michael's subspecies on an automatic species) leaves the subspecies unscored
(outcome truth_undetermined); an empty truth_subspecies_id means the
genome has no subspecies, so a subspecies result is an overcall signal.
truth_kind is one of:
  named            project or expert label (species, optionally subspecies)
  outside_genus    expert label "outside Enterobacter"
  no_reference     a species without a package reference, or a novel lineage
  unresolved       expert label "unresolved"; descriptive only
  context          study context or unnamed record; descriptive only
  crosswalk        Layer 3: expected IDs come from the NCBI-name crosswalk
  crosswalk_none   Layer 3 unnamed or outside-genus NCBI names; distribution only
Results manifest (TSV): genome_id, report_tsv, report_sha256. A genome with no
manifest row, or whose run failed, is scored as run_failed and stays counted.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

SEED = "Enterobacter-validation-block-20261007-v1"
TRUTH_KINDS = {"named", "outside_genus", "no_reference", "unresolved", "context", "crosswalk", "crosswalk_none"}
TRUTH_LAYERS = ("named", "outside_genus", "no_reference")  # kinds entering accuracy denominators
BLOCKED = TRUTH_LAYERS + ("crosswalk",)  # kinds entering the block sensitivity
LABEL_SOURCES = {"project_label", "auto_single_species", "auto_tree_resolved", "michael", "none"}
EXPERT_SOURCES = ("auto_single_species", "auto_tree_resolved", "michael")
UNDETERMINED = "undetermined"  # expert truth_subspecies_id when Michael could not determine the subspecies
REPORT_COLUMNS = ("experimental_eligibility_status", "experimental_eligibility_reasons", "hypothetical_species_id",
                  "hypothetical_subspecies_id", "near_boundary_species_id", "query_sha256")
csv.field_size_limit(1_000_000_000)


def _split(value):
    return [v for v in (value or "").split(";") if v]


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_result(path, expected_sha256):
    """One analyze report.tsv with exactly one row, verified by hash."""
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected_sha256:
        raise ValueError(f"report hash mismatch: {path}")
    rows = read_tsv(path)
    if len(rows) != 1 or any(c not in rows[0] for c in REPORT_COLUMNS):
        raise ValueError(f"unexpected report shape: {path}")
    return rows[0]


def score(truth, result):
    """Per-genome outcome; result is a report.tsv row or None for a failed run."""
    kind, source = truth["truth_kind"], truth.get("label_source", "none")
    if kind not in TRUTH_KINDS:
        raise ValueError(f"unknown truth_kind {kind!r}")
    if source not in LABEL_SOURCES:
        raise ValueError(f"unknown label_source {source!r}")
    if kind in TRUTH_LAYERS and source == "none":
        raise ValueError("accuracy rows need a label_source")
    if kind in ("crosswalk", "crosswalk_none") and source != "none":
        raise ValueError("crosswalk rows must use label_source none")
    if truth.get("subspecies_structural", "") not in ("", "True", "False"):
        raise ValueError("subspecies_structural must be True, False or empty")
    passing = set(_split(truth.get("ani_species_passing")))
    if passing and source not in EXPERT_SOURCES:
        raise ValueError("ani_species_passing is for expert-stratum rows only")
    if source == "auto_single_species" and (kind != "named" or passing != {truth["truth_species_id"]}):
        raise ValueError("auto_single_species needs exactly the label's species passing")
    if source == "auto_tree_resolved" and (kind != "named" or len(passing) < 2 or truth["truth_species_id"] not in passing):
        raise ValueError("auto_tree_resolved needs several passing species including the label's")
    structural = kind == "named" and source in EXPERT_SOURCES and passing == {truth["truth_species_id"]}
    out = dict(genome_id=truth["genome_id"], layer=truth["layer"], truth_kind=kind, label_source=source,
               structural=structural, subspecies_structural=truth.get("subspecies_structural") == "True",
               result_outside_passing=False, species_outcome="",
               subspecies_outcome="", note_outcome="", result_species_id="", result_subspecies_id="",
               reasons="", note_species_id="")
    if result is None or result["experimental_eligibility_status"] == "failed":
        out.update(species_outcome="run_failed", reasons=result["experimental_eligibility_reasons"] if result else "")
        return out
    if not truth.get("query_sha256") or truth["query_sha256"] != result["query_sha256"]:
        raise ValueError(f"report query_sha256 does not match the truth row for {truth['genome_id']}")
    sp, ss = result["hypothetical_species_id"], result["hypothetical_subspecies_id"]
    note = result["near_boundary_species_id"]
    if ss and not sp:
        raise ValueError("subspecies result without species result")
    if note and sp:
        raise ValueError("near-boundary note alongside a species result")
    out.update(result_species_id=sp, result_subspecies_id=ss, reasons=result["experimental_eligibility_reasons"],
               note_species_id=note)
    if not sp and not out["reasons"]:
        raise ValueError("withheld record without a reason")
    out["result_outside_passing"] = bool(sp) and source in EXPERT_SOURCES and sp not in passing

    if kind == "named":
        truth_sp, truth_ss = truth["truth_species_id"], truth["truth_subspecies_id"]
        if not truth_sp:
            raise ValueError("named truth without a species")
        out["species_outcome"] = "withheld" if not sp else "correct" if sp == truth_sp else "wrong"
        if truth_ss == UNDETERMINED:
            if truth.get("subspecies_label_source") != "michael":
                raise ValueError("only Michael's subspecies labels can be undetermined")
            out["subspecies_outcome"] = "truth_undetermined"
        elif truth_ss:
            out["subspecies_outcome"] = ("withheld_at_species" if not sp else "wrong_species" if sp != truth_sp else
                                         "species_only" if not ss else "correct" if ss == truth_ss else "wrong_subspecies")
        else:
            out["subspecies_outcome"] = "overcall_signal" if ss else "not_applicable"
        if not sp:
            out["note_outcome"] = "absent" if not note else "names_truth_species" if note == truth_sp else "names_other_species"
    elif kind in ("outside_genus", "no_reference"):
        out["species_outcome"] = "wrong" if sp else "correct_withheld"
        out["subspecies_outcome"] = "not_applicable"
        out["note_outcome"] = "" if sp else "absent" if not note else "names_a_species"
    elif kind in ("unresolved", "context", "crosswalk_none"):
        out["species_outcome"] = "descriptive_result" if sp else "descriptive_withheld"
        out["note_outcome"] = "" if sp else "absent" if not note else "names_a_species"
    else:  # crosswalk
        expected = set(_split(truth["expected_species_ids"]))
        if not expected:
            raise ValueError("crosswalk truth without expected species")
        out["species_outcome"] = "withheld" if not sp else "consistent" if sp in expected else "discordant"
        exp_ss = truth["expected_subspecies_id"]
        if exp_ss and out["species_outcome"] == "consistent":
            out["subspecies_outcome"] = ("consistent_at_species_no_subspecies" if not ss else
                                         "consistent" if ss == exp_ss else "discordant")
        out["note_outcome"] = "" if sp else "absent" if not note else (
            "names_expected_species" if note in expected else "names_other_species")
    return out


def _binom_cdf(k, n, p):
    if p <= 0:
        return 1.0
    if p >= 1:
        return 1.0 if k >= n else 0.0
    return min(1.0, sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                                 + i * math.log(p) + (n - i) * math.log1p(-p)) for i in range(k + 1)))


def clopper_pearson(k, n, alpha=0.05):
    """Exact two-sided interval by bisection on the binomial CDF; None when n is 0."""
    if n == 0:
        return None

    def solve(f):
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) else (lo, mid)
        return (lo + hi) / 2
    lower = 0.0 if k == 0 else solve(lambda p: 1 - _binom_cdf(k - 1, n, p) < alpha / 2)
    upper = 1.0 if k == n else solve(lambda p: _binom_cdf(k, n, p) > alpha / 2)
    return [round(lower, 6), round(upper, 6)]


def summarise(rows):
    """Counts for one group of scored rows."""
    acc = [r for r in rows if r["truth_kind"] in TRUTH_LAYERS]
    evaluable = [r for r in acc if not r["structural"]]
    results = [r for r in evaluable if r["species_outcome"] in ("correct", "wrong")]
    wrong = sum(r["species_outcome"] == "wrong" for r in evaluable)
    failed = sum(r["species_outcome"] == "run_failed" for r in rows)
    named = [r for r in acc if r["truth_kind"] == "named"]
    named_results = [r for r in named if r["species_outcome"] in ("correct", "wrong")]
    withheld_reasons = Counter(x for r in acc if r["species_outcome"] in ("withheld", "correct_withheld")
                               for x in _split(r["reasons"]))
    sub_scored = [r for r in named if r["subspecies_outcome"] in ("correct", "wrong_subspecies")]
    sub_open = [r for r in sub_scored if not r["subspecies_structural"]]
    out = dict(genomes=len(rows), run_failed=failed, accuracy_denominator=len(acc),
               species_outcomes=dict(Counter(r["species_outcome"] for r in rows)),
               subspecies_outcomes=dict(Counter(r["subspecies_outcome"] for r in rows if r["subspecies_outcome"])),
               note_outcomes=dict(Counter(r["note_outcome"] for r in rows if r["note_outcome"])),
               withheld_reasons=dict(sorted(withheld_reasons.items())),
               label_sources=dict(Counter(r["label_source"] for r in acc)),
               species_results=len(results), wrong_results=wrong,
               wrong_among_results_ci95=clopper_pearson(wrong, len(results)),
               structural_rows=sum(bool(r["structural"]) for r in acc),
               wrong_result_basis="excludes structural rows (label is the only species passing 95/80)",
               structural_rows_wrong=sum(r["species_outcome"] == "wrong" for r in acc if r["structural"]),
               tool_species_outside_passing=sum(bool(r["result_outside_passing"]) for r in rows),
               subspecies_structural_outcomes=dict(Counter(r["subspecies_outcome"] for r in named
                                                           if r["subspecies_structural"] and r["subspecies_outcome"])),
               wrong_subspecies=[sum(r["subspecies_outcome"] == "wrong_subspecies" for r in sub_open), len(sub_open)],
               wrong_subspecies_ci95=clopper_pearson(sum(r["subspecies_outcome"] == "wrong_subspecies" for r in sub_open),
                                                     len(sub_open)),
               wrong_subspecies_basis="correct species with a subspecies result; excludes subspecies_structural rows",
               named_coverage=[len(named_results), len(named)],
               named_coverage_ci95=clopper_pearson(len(named_results), len(named)))
    layer3 = [r for r in rows if r["truth_kind"] == "crosswalk"]
    if layer3:
        out["crosswalk_species"] = dict(Counter(r["species_outcome"] for r in layer3))
        out["crosswalk_returned_species"] = dict(Counter(r["result_species_id"] or "withheld" for r in layer3))
    return out


def block_sensitivity(rows):
    """One accuracy or crosswalk genome per block (clone block or linked component), chosen by seeded order."""
    chosen = {}
    for r in sorted(rows, key=lambda r: hashlib.sha256(f"{SEED}|{r['genome_id']}".encode()).hexdigest()):
        if r["truth_kind"] in BLOCKED:
            if not r["block"]:
                raise ValueError(f"no block for {r['genome_id']}")
            chosen.setdefault(r["block"], r)
    return summarise(list(chosen.values()))


def component_weighted(rows):
    """Layer 3: each linked component carries weight 1, split among its scored members."""
    sizes = Counter(r["component_id"] for r in rows)
    weights = defaultdict(float)
    for r in rows:
        weights[r["species_outcome"]] += 1 / sizes[r["component_id"]]
    return {k: round(v, 6) for k, v in sorted(weights.items())}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True, help="genome_id, report_tsv, report_sha256")
    parser.add_argument("--group", action="append", default=[], help="extra truth column to group by (repeatable)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("refusing to overwrite output")
    truth = read_tsv(args.truth)
    header = set(truth[0]) if truth else set()
    unknown = sorted({c for g in args.group for c in g.split(",")} - header)
    if unknown:
        raise ValueError(f"--group columns not in the truth table: {unknown}")
    ids = [t["genome_id"] for t in truth]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate genome_id in truth")
    manifest = {m["genome_id"]: m for m in read_tsv(args.results)}
    if set(manifest) - set(ids):
        raise ValueError("results for genomes not in truth")
    for t in truth:
        for column in ("block", "component_id", "expected_species_ids", "expected_subspecies_id", "truth_species_id",
                       "truth_subspecies_id", "ani_species_passing", "subspecies_structural", "query_sha256",
                       "subspecies_label_source",
                       *(c for g in args.group for c in g.split(","))):
            t.setdefault(column, "")
        if t["truth_kind"] in BLOCKED and not t["block"]:
            raise ValueError(f"no block for {t['genome_id']}")
        t.setdefault("label_source", "none")
    scored = []
    for t in truth:
        m = manifest.get(t["genome_id"])
        result = load_result(m["report_tsv"], m["report_sha256"]) if m else None
        columns = list(dict.fromkeys(c for g in args.group for c in g.split(",")))
        scored.append(dict(score(t, result), block=t["block"], component_id=t["component_id"], **{c: t[c] for c in columns}))
    groups = defaultdict(list)
    for r in scored:
        groups[r["layer"]].append(r)
        for g in args.group:
            groups[f"{r['layer']}|" + ",".join(f"{c}={r[c]}" for c in g.split(","))].append(r)
    summary = dict(scorer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   truth_sha256=hashlib.sha256(args.truth.read_bytes()).hexdigest(),
                   results_manifest_sha256=hashlib.sha256(args.results.read_bytes()).hexdigest(),
                   interval="exact two-sided 95% Clopper-Pearson; genomes treated as independent",
                   groups={k: summarise(v) for k, v in sorted(groups.items())},
                   block_sensitivity={layer: block_sensitivity(v) for layer, v in sorted(groups.items()) if "|" not in layer},
                   component_weighted={key: component_weighted([r for r in v if r["truth_kind"] == "crosswalk"])
                                       for key, v in sorted(groups.items())
                                       if any(r["truth_kind"] == "crosswalk" for r in v)})
    args.output.mkdir(parents=True)
    with (args.output / "scored-genomes.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scored[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(scored)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
