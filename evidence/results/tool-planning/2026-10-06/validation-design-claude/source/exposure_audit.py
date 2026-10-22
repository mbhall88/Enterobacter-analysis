"""Metadata-only exposure audit of the 4,857-unit public reserve.

Finds reserve units already present in earlier project analyses and units
linked to earlier analysed genomes, then holds their whole linked components,
following the confirmation protocol's conservative link rule (accession
lineage, BioSample, BioProject, exact normalised strain label). Also flags
culture-collection strain designations, NCBI type-material records and
strain labels matching an alias of a package type strain (re-sequencings of a
reference). Reads text files only; no genome sequence is opened. Outputs go
to a new directory.
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "src"))
from enterobacter_classifier.reporting import ReportingContext  # noqa: E402

B = ROOT / "results/tool-planning/2026-10-06"
POOL = B / "prospective-confirmation-astra/inputs/execution-freeze/design/remaining-reserve-metadata.tsv"
SEEDS = B / "evaluation-freeze-design-astra/outputs/prior-accession-seeds-with-roles.tsv"
SUPPLEMENT = B / "confirmation-relatedness-astra/outputs/new-prior-exposure-supplement.tsv"
SELECTED = B / "evaluation-freeze-design-astra/reserve-protected-alternative/outputs/confirmation-manifest.tsv"
HOLD = B / "evaluation-freeze-design-astra/reserve-protected-alternative/outputs/confirmation-linked-hold.tsv"
NCBI = ROOT / "results/ani-cases/2026-10-02/task12-astra-sampling-audit/inputs/ncbi-enterobacter-current.jsonl"
TYPES = ROOT / "manifests/types.tsv"
PACKAGE = B / "portable-reference-package-astra/relocated-package/portable.json"
PACKAGE_SHA256 = "12aa0705b774de8ea262869fdc4c85d513c3aa97d39c674f08e85d1c9055b88f"
OVERLAY = {"GCA_003057745.2", "GCA_965138355.1"}  # dated confirmation-review correction
SCAN_ROOTS = ["results", "manifests", "config", "docs"]
# Runs that list reserve units by construction (this design and its review).
EXCLUDED_RUNS = ("results/tool-planning/2026-10-06/validation-design-claude/",
                 "results/tool-planning/2026-10-06/validation-design-review-claude/")
# Files that mention genomes only as metadata (no project sequence analysis or decision discussion).
# Anything not matched here counts as an exposure, so unclassified files stay conservative.
METADATA_FILES = [
    (r"(^|/)gtdb-r232-clusters\.(tsv|json)$", "GTDB species-cluster membership table"),
    (r"(^|/)(boundary-)?gtdb-cluster\.tsv$", "GTDB species-cluster membership table"),
    (r"(^|/)GTDB-clusters\.tsv$", "GTDB species-cluster membership table"),
    (r"/reference/qc_failed_r232\.tsv$", "GTDB release QC-failed table"),
    (r"/tables/focus-genus-clusters\.tsv$", "GTDB cluster summary"),
    (r"/database/library_report\.tsv$", "Kraken database library list"),
    (r"evaluation-freeze-design-astra/outputs/confirmation-(accessions\.txt|manifest\.tsv)$",
     "superseded initial confirmation draw; replaced before any sequence access"),
    (r"\.jsonl$", "NCBI metadata report"),
    (r"task12-astra-gtdb-naming-map/outputs/", "GTDB/NCBI name map"),
    (r"task12-astra-taxonomy-reconciliation/", "NCBI/GTDB taxonomy metadata reconciliation"),
    (r"task12-astra-sampling-audit/", "NCBI sampling metadata audit"),
    (r"ncbi-ANI-context\.tsv$", "NCBI-reported ANI metadata"),
    (r"(selection-validation\.(json|log)$|later-pool-lossless/|strain-delimiter-cases\.tsv$|semicolon-source-example\.json$)",
     "selection or metadata-parser validation"),
    (r"task02-genus-census-sol/", "type-material census metadata"),
    (r"^manifests/types\.tsv$", "type-material candidate list"),
    (r"^docs/task-status/15-validation-design", "this design's own review claim"),
]
TYPE_MATERIAL_FILES = r"(task02-genus-census-sol/|^manifests/types\.tsv$|task02-type-assembly-qc-astra/)"
LISTING_POOL_UNITS = 1000      # a file naming >= this many reserve units is a metadata listing
LISTING_ACCESSIONS = 1000      # a file naming >= this many accessions contributes no non-reserve exposures
CULTURE = re.compile(r"\b(ATCC|DSM|CIP|LMG|NCTC|JCM|CCUG|NBRC|KCTC|CECT|FDAARGOS|NCIMB|CGMCC|CCTCC|BCRC|ICMP|NRRL|IAM|IFO|NCDO|"
                     r"CCM|CCOS|GDMCC|CNCTC|NCCB|CFBP|NCPPB|HAMBI|NBIMCC|WDCM|HUT)(?:[ _:-]?(?:AB|BAA|B)?[ _:-]?)\d[\d.]*", re.I)
EXCLUDING_TYPE_LABELS = {"TYPE_MATERIAL", "TYPE_MATERIAL_NEOTYPE", "TYPE_MATERIAL_REFERENCE"}  # not heterotypic synonyms
NULL_LABELS = {"", "missing", "not provided", "not applicable", "unknown", "na", "n/a"}
csv.field_size_limit(100_000_000)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def alias_key(text):
    """Collection designations compare without spaces, punctuation or a trailing type-strain T."""
    key = re.sub(r"[\s_:.\-]", "", text.upper())
    return re.sub(r"(?<=\d)T$", "", key)


def type_strain_aliases():
    """Alias key -> package type anchors whose type strain carries it. Short or purely alphabetic aliases are skipped."""
    refs = ReportingContext(PACKAGE, expected_sha256=PACKAGE_SHA256).evidence_references()
    package = {r["reference_id"].split(":", 1)[1] for r in refs}
    aliases = defaultdict(set)
    for r in rows(TYPES):
        if r["selected_accession_version"] in package:
            for alias in r["type_strain_aliases"].split(";"):
                key = alias_key(alias)
                if len(key) >= 4 and re.search(r"\d", key):
                    aliases[key].add(f"{r['selected_accession_version']} {r['taxon_name']} ({alias.strip()})")
    return aliases


def core(accession):
    return accession.split("_", 1)[1].rsplit(".", 1)[0]


def snake(x):
    if isinstance(x, list):
        return [snake(v) for v in x]
    if isinstance(x, dict):
        return {re.sub(r"(?<!^)(?=[A-Z])", "_", k).lower(): snake(v) for k, v in x.items()}
    return x


def tokens(raw):
    """Same token rule as confirmation-results-review-astra/source/audit_exposure_overlay.py."""
    r = snake(raw)
    info = r["assembly_info"]
    bio = info.get("biosample", {})
    accessions = {r["accession"], r.get("paired_accession", ""), info.get("paired_assembly", {}).get("accession", "")} - {""}
    labels = [r.get("organism", {}).get("infraspecific_names", {}).get("strain", ""), bio.get("strain", "")] + \
             [v.get("value", "") for v in bio.get("attributes", []) if v.get("name") in {"strain", "isolate", "culture_collection"}]
    projects = {info.get("bioproject_accession", "")} - {""}
    stack = [info.get("bioproject_lineage", []), bio.get("bioprojects", [])]
    while stack:
        item = stack.pop()
        if isinstance(item, list):
            stack += item
        elif isinstance(item, dict):
            if re.fullmatch(r"PRJ(?:NA|EB|DB)[0-9]+", str(item.get("accession", ""))):
                projects.add(item["accession"])
            stack += list(item.values())
    return {"accession_lineage": {v.rsplit(".", 1)[0] for v in accessions}, "biosample": {bio.get("accession", "")} - {""},
            "bioproject": projects, "strain_label_conservative_hold": {v.strip().casefold() for v in labels} - NULL_LABELS}


def frozen(r):
    return {"accession_lineage": {x.rsplit(".", 1)[0] for x in r["accessions"].split(";") if x},
            "biosample": set(filter(None, r["biosample"].split(";"))),
            "bioproject": set(filter(None, r["all_known_bioprojects"].split(";"))),
            "strain_label_conservative_hold": set(json.loads(r["strain_labels_json"]))}


def scan():
    """Map each scanned file to the set of versioned GCA/GCF accessions it mentions."""
    globs = ["-g", "*.{tsv,csv,json,jsonl,txt,md,nwk,newick,treefile,tree,log,list}"]
    out = subprocess.run(["rg", "--no-ignore", "-o", "-N", "--with-filename", "--max-filesize", "400M", *globs,
                          r"GC[AF]_[0-9]{9}\.[0-9]+", *SCAN_ROOTS], cwd=ROOT, capture_output=True, text=True)
    assert out.returncode in (0, 1), out.stderr[:500]
    files = defaultdict(set)
    for line in out.stdout.splitlines():
        path, _, accession = line.rpartition(":")
        if not path.startswith(EXCLUDED_RUNS):
            files[path].add(accession)
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite output")
    pool = rows(POOL)
    assert len(pool) == 4857
    unit_by_core = {core(a): r["assembly_unit"] for r in pool for a in r["accessions"].split(";") if a}
    by_unit = {r["assembly_unit"]: r for r in pool}
    known = {core(r["accession"]) for r in rows(SEEDS)} | {core(r["accession"]) for r in rows(SUPPLEMENT)}
    for path in (SELECTED, HOLD):
        known |= {core(a) for r in rows(path) for a in r["accessions"].split(";") if a}
    index = defaultdict(set)
    for r in pool:
        for kind, values in frozen(r).items():
            for v in values:
                index[kind, v].add(r["assembly_unit"])

    aliases = type_strain_aliases()
    files = scan()
    direct = defaultdict(set)        # reserve unit -> files naming it
    others = defaultdict(set)        # earlier-analysed non-reserve accession -> files
    type_material = defaultdict(set)
    listings, metadata_only = [], []
    for path, accessions in sorted(files.items()):
        pool_units = {unit_by_core[core(a)] for a in accessions if core(a) in unit_by_core}
        if re.search(TYPE_MATERIAL_FILES, path):
            for u in pool_units:
                type_material[u].add(path)
        if len(pool_units) >= LISTING_POOL_UNITS:
            listings.append(dict(path=path, reserve_units=len(pool_units), accessions=len(accessions), reason="bulk listing"))
            continue
        rule = next((why for pattern, why in METADATA_FILES if re.search(pattern, path)), None)
        if rule:
            if pool_units:
                metadata_only.append(dict(path=path, reserve_units=len(pool_units), accessions=len(accessions), reason=rule,
                                          units=";".join(sorted(pool_units))))
            continue
        for u in pool_units:
            direct[u].add(path)
        if len(accessions) < LISTING_ACCESSIONS:
            for a in accessions:
                if core(a) not in unit_by_core and core(a) not in known:
                    others[a].add(path)

    metadata = {}
    with NCBI.open() as handle:
        for line in handle:
            if line.strip():
                record = json.loads(line)
                metadata[record["accession"]] = record
    type_label = {}
    for accession, record in metadata.items():
        label = snake(record).get("type_material", {}).get("type_label")
        if label and core(accession) in unit_by_core:
            type_label[unit_by_core[core(accession)]] = label
    links, unlinkable = [], []
    for accession, paths in sorted(others.items()):
        record = metadata.get(accession)
        if record is None:
            unlinkable.append(dict(accession=accession, files=len(paths), example=sorted(paths)[0]))
            continue
        for kind, values in tokens(record).items():
            for v in values:
                for target in sorted(index[kind, v]):
                    links.append(dict(earlier_accession=accession, kind=kind, value=v, reserve_unit=target,
                                      component_id=by_unit[target]["component_id"], example_file=sorted(paths)[0]))

    held_components = {by_unit[u]["component_id"] for u in direct} | {l["component_id"] for l in links} | \
                      {by_unit[u]["component_id"] for u in OVERLAY}
    out = []
    for r in pool:
        u = r["assembly_unit"]
        labels = json.loads(r["strain_labels_json"]) + [r["strain"]]
        culture = sorted({m.group(0).strip() for x in labels for m in CULTURE.finditer(x)})
        matches = sorted({a for key in {alias_key(c) for c in culture} | {alias_key(x) for x in labels} for a in aliases.get(key, ())})
        reasons = []
        flags = [f for f, on in (("type_material_file", u in type_material), ("culture_collection_strain", bool(culture)),
                                 ("package_type_strain_alias", bool(matches)), ("ncbi_type_material", u in type_label)) if on]
        resequenced_reference = bool(matches) or type_label.get(u) in EXCLUDING_TYPE_LABELS
        if u in OVERLAY:
            reasons.append("confirmation_review_overlay")
        if u in direct:
            reasons.append("named_in_earlier_analysis")
        if any(l["reserve_unit"] == u for l in links):
            reasons.append("linked_to_earlier_analysed_genome")
        if r["component_id"] in held_components and not reasons:
            reasons.append("component_of_exposed_unit")
        out.append(dict(assembly_unit=u, component_id=r["component_id"], ncbi_names=r["ncbi_names"], qc_stratum=r["qc_stratum"],
                        disposition="held_exposed" if reasons else "effective_reserve", exposure_reasons=";".join(reasons),
                        direct_files=";".join(sorted(direct.get(u, ())))[:2000],
                        culture_collection_flag=";".join(culture), package_type_strain_alias=";".join(matches),
                        ncbi_type_label=type_label.get(u, ""), provenance_flags=";".join(flags),
                        expert_eligible=not reasons and not resequenced_reference))
    args.output.mkdir()
    fields = list(out[0])
    with (args.output / "reserve-exposure-dispositions.tsv").open("x", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    for name, table in (("metadata-links.tsv", links), ("listing-files.tsv", listings), ("metadata-only-files.tsv", metadata_only),
                        ("unlinkable-earlier-accessions.tsv", unlinkable)):
        with (args.output / name).open("x", newline="") as handle:
            if table:
                w = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(table)
    effective = [r for r in out if r["disposition"] == "effective_reserve"]
    summary = dict(
        reserve_units=len(pool), files_scanned_with_accessions=len(files), listing_files=len(listings),
        directly_named_units=len(direct), earlier_nonreserve_accessions=len(others),
        earlier_accessions_without_metadata=len(unlinkable), metadata_link_rows=len(links),
        linked_units=len({l["reserve_unit"] for l in links}), held_components=len(held_components),
        held_units=len(out) - len(effective), effective_units=len(effective),
        effective_components=len({r["component_id"] for r in effective}),
        culture_collection_flagged_effective=sum(bool(r["culture_collection_flag"]) for r in effective),
        type_material_flagged_effective=sum(bool(r["ncbi_type_label"]) for r in effective),
        package_type_alias_units=sum(bool(r["package_type_strain_alias"]) for r in out),
        package_type_alias_effective=sum(bool(r["package_type_strain_alias"]) for r in effective),
        units_without_ncbi_metadata=sorted(set(by_unit) - {unit_by_core[core(a)] for a in metadata if core(a) in unit_by_core}),
        expert_eligible_units=sum(r["expert_eligible"] for r in out),
        metadata_only_files=len(metadata_only), metadata_rules=[why for _, why in METADATA_FILES],
        inputs={str(p): sha(p) for p in (POOL, SEEDS, SUPPLEMENT, SELECTED, HOLD, NCBI, TYPES)} | {str(PACKAGE): PACKAGE_SHA256},
        source_sha256=sha(__file__), expert_eligibility="effective reserve, not a probable re-sequencing of a package type strain "
        "(alias match or NCBI type material other than a heterotypic synonym)",
        thresholds=dict(listing_reserve_units=LISTING_POOL_UNITS, listing_accessions=LISTING_ACCESSIONS),
        excluded_runs=list(EXCLUDED_RUNS), sequence_files_opened=0,
        interpretation="conservative metadata correlation holds; not isolate identity or sequence-independence proof")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "inputs"}, indent=2))


if __name__ == "__main__":
    main()
