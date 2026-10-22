"""Label transfer for the Task 15 expert stratum (validation-design.md section 7).

Written and tested before any public validation genome is downloaded; its
SHA256 is posted at freeze and it is not edited afterwards.

prepare  reads, for the 100 sampled genomes, the labelling skani hits (one
         search per genome against the 1,086 comparators), the rebuilt BAC120
         and Verticall trees, the paper trees, the GTDB-Tk identify summaries
         and the Verticall QC gate. It writes the pre-seal table (automatic
         labels and the edge-case triggers), Michael's packets and the trees
         relabelled with coded IDs.
seal     merges Michael's labels for the edge cases, checks them against the
         pre-seal table and the package's species and subspecies names, and
         writes the sealed label table and the scorer's expert truth table.
loo      replays the rule leave-one-out on the paper trees with the saved
         chromosome pairs, as a check that it is operational (in-sample).

The trees are read with the project's tested Newick parser and rerooting
(scripts/gtdb_context.py, scripts/root_tree.py).
"""
import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
from gtdb_context import format_newick, parse_newick  # noqa: E402
from root_tree import find_root_edge, reroot, to_graph  # noqa: E402

ARTIFACTS = ROOT / "results/tool-planning/2026-10-06/validation-design-claude/outputs/design-artifacts-attempt08"
COMPARATORS = ARTIFACTS / "labelling-comparators.tsv"
TREE_TIPS = ARTIFACTS / "cohort-tree-tips.tsv"
SAMPLE = ARTIFACTS / "expert-sample-option-stratified-n100.tsv"
PACKAGE = ROOT / "results/tool-planning/2026-10-06/portable-reference-package-astra/relocated-package/portable.json"
PACKAGE_SHA256 = "12aa0705b774de8ea262869fdc4c85d513c3aa97d39c674f08e85d1c9055b88f"

PASS_ANI, PASS_AF = 95.0, 80.0            # DEC-007
NEAR_ANI, NEAR_AF = 99.9, 90.0            # confirmation protocol's near-clone convention
MIN_UNIQUE_MARKERS = 100                  # Task 11 low-marker flag (unique BAC120 markers, GTDB-Tk identify)
OUTGROUPS = {"klebsiella_aerogenes": "GCA_019048125.1", "cronobacter_sakazakii": "GCA_000982825.1"}
PRIMARY = "klebsiella_aerogenes"
METHODS = ("bac120", "verticall")
CODE_SEED = "Enterobacter-validation-coded-id-20261007-v1"
MIXED, UNASSIGNED, AMBIGUOUS = "mixed", "unassigned", "ambiguous"
DECISIONS = {"species", "species_without_reference", "unresolved", "outside_genus", "novel_lineage"}
NO_SUBSPECIES, UNDETERMINED = "none", "undetermined"   # Michael's subspecies entries besides a package epithet
NEAR_NAME_EDITS = 2                       # a species_without_reference name this close to a package species is refused
csv.field_size_limit(1_000_000_000)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    with Path(path).open("x", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def coded_id(accession):
    return "EXP-" + hashlib.sha256(f"{CODE_SEED}|{accession}".encode()).hexdigest()[:8].upper()


# ---------------------------------------------------------------------------------------------
# Labels


def load_comparators(path=None):
    """Comparators keyed by whole-assembly path (the skani Ref_file) and by accession."""
    rows = read_tsv(path or COMPARATORS)
    by_path = {r["whole_fasta"]: r for r in rows}
    if len(by_path) != len(rows):
        raise ValueError("duplicate comparator path")
    return by_path, {r["accession"]: r for r in rows}


def load_tip_labels(path=None):
    """Cohort tip -> (species, subspecies epithet) for labelled tips, None for unassigned tips."""
    return {r["accession"]: ((r["project_species"], r["project_subspecies"]) if r["tip_label_kind"] == "named" else None)
            for r in read_tsv(path or TREE_TIPS)}


def species_with_subspecies(labels):
    return {lab[0] for lab in labels.values() if lab and lab[1]}


# ---------------------------------------------------------------------------------------------
# ANI evidence


def _pct(value, what):
    v = float(value)
    if not 0.0 <= v <= 100.0:
        raise ValueError(f"{what} outside 0-100: {value}")
    return v


def ani_evidence(hits, comparators_by_path, query_fasta):
    """S, U, near-clone flag and per-species top comparators from one genome's skani search rows."""
    rows = []
    for h in hits:
        if h["Query_file"] != query_fasta:
            raise ValueError(f"skani query {h['Query_file']} is not the sampled FASTA {query_fasta}")
        if h["Ref_file"] not in comparators_by_path:
            raise ValueError(f"skani reference outside the comparator manifest: {h['Ref_file']}")
        c = comparators_by_path[h["Ref_file"]]
        ani, afr, afq = (_pct(h[k], k) for k in ("ANI", "Align_fraction_ref", "Align_fraction_query"))
        rows.append(dict(accession=c["accession"], named=c["label_kind"] == "named", species=c["project_species"],
                         subspecies=c["project_subspecies"], reference=c["is_package_reference"] == "True",
                         ani=ani, af_ref=afr, af_query=afq, passes=ani >= PASS_ANI and min(afr, afq) >= PASS_AF))
    if len({r["accession"] for r in rows}) != len(rows):
        raise ValueError("duplicate comparator in skani hits")
    rows.sort(key=lambda r: (-r["ani"], r["accession"]))
    return dict(hits=rows,
                S=sorted({r["species"] for r in rows if r["named"] and r["passes"]}),
                U=sorted(r["accession"] for r in rows if not r["named"] and r["passes"]),
                near_clone=any(r["reference"] and r["ani"] >= NEAR_ANI and min(r["af_ref"], r["af_query"]) >= NEAR_AF
                               for r in rows))


def top_within_species(evidence, species):
    """The tool's nearest-reference logic applied to comparators: highest ANI within the species, ties kept."""
    within = [r for r in evidence["hits"] if r["named"] and r["species"] == species]
    if not within:
        return dict(subspecies=None, tie=False, is_reference=False, accession="", ani=None)
    top = [r for r in within if r["ani"] == within[0]["ani"]]
    groups = {r["subspecies"] for r in top}
    return dict(subspecies=within[0]["subspecies"] if len(groups) == 1 else None, tie=len(groups) > 1,
                is_reference=len(groups) == 1 and any(r["reference"] for r in top),
                accession=";".join(r["accession"] for r in top), ani=within[0]["ani"])


# ---------------------------------------------------------------------------------------------
# Trees


def rooted(newick, outgroup):
    tree = parse_newick(newick)
    adjacency, tip_name = to_graph(tree)
    u, v, length, support = find_root_edge(adjacency, tip_name, {outgroup})
    return reroot(adjacency, tip_name, u, v, length, support)


def _cohort_below(tree, cohort):
    below = [frozenset()] * tree.size
    for node in tree.postorder():
        if not tree.children[node]:
            below[node] = frozenset([tree.label[node]]) & cohort
        else:
            below[node] = frozenset().union(*(below[c] for c in tree.children[node]))
    return below


def _reading(tips, labels, has_subspecies):
    if any(labels[t] is None for t in tips):
        return UNASSIGNED, ""
    species = {labels[t][0] for t in tips}
    if len(species) > 1:
        return MIXED, ""
    (sp,) = species
    if sp not in has_subspecies:
        return sp, ""
    subs = {labels[t][1] for t in tips}
    return sp, (subs.pop() if len(subs) == 1 and "" not in subs else AMBIGUOUS)


def read_clades(tree, labels, focal):
    """Clade reading for each focal tip: the cohort tips under its lowest ancestor that has any.

    labels covers the cohort tips; tips in neither labels nor focal (other validation genomes) are
    ignored. A focal tip is never part of its own clade, which is what the leave-one-out check needs.
    """
    cohort = frozenset(labels)
    below = _cohort_below(tree, cohort)
    has_sub = species_with_subspecies(labels)
    index = {tree.label[n]: n for n in tree.tips()}
    out = {}
    for name in focal:
        node, clade = index[name], frozenset()
        while not clade:
            node = tree.parent[node]
            if node == -1:
                raise ValueError(f"no cohort tip shares a clade with {name}")
            clade = below[node] - {name}
        species, subspecies = _reading(clade, labels, has_sub)
        out[name] = dict(species=species, subspecies=subspecies, clade_size=len(clade))
    return out


def label_clades(tree, labels):
    """Each species and subspecies label -> cohort tips of the smallest clade holding all its tips."""
    cohort = frozenset(labels)
    below = _cohort_below(tree, cohort)
    index = {tree.label[n]: n for n in tree.tips()}
    groups = defaultdict(list)
    for tip, lab in labels.items():
        if lab:
            groups[("species", lab[0])].append(index[tip])
            if lab[1]:
                groups[("subspecies", lab[0], lab[1])].append(index[tip])
    return {key: below[tree.mrca(nodes)] for key, nodes in groups.items()}


def changed_labels(paper_tree, rebuilt_tree, labels):
    paper, rebuilt = label_clades(paper_tree, labels), label_clades(rebuilt_tree, labels)
    if set(paper) != set(rebuilt):
        raise ValueError("label sets differ between the paper and rebuilt trees")
    return {key for key in paper if paper[key] != rebuilt[key]}


# ---------------------------------------------------------------------------------------------
# The rule (section 7.3)


def decide(ev):
    """Automatic labels and triggers for one genome.

    ev: S (species names), U (unassigned comparators passing), readings[method][rooting] =
    {species, subspecies}, markers (unique_markers, translation_table, aligned), verticall_ok,
    changed (set of label keys changed in either tree, any rooting), top (callable species ->
    top_within_species result).
    """
    S, primary = set(ev["S"]), {m: ev["readings"][m][PRIMARY] for m in METHODS}
    triggers = []
    if not S:
        triggers.append("S_empty")
    if ev["U"]:
        triggers.append("unassigned_comparator_passes")
    for m in METHODS:
        if primary[m]["species"] in (MIXED, UNASSIGNED):
            triggers.append(f"{m}_clade_{primary[m]['species']}")
    tree_species = {primary[m]["species"] for m in METHODS}
    named = tree_species - {MIXED, UNASSIGNED}
    if len(named) == 2:
        triggers.append("trees_disagree")
    T = named.pop() if len(tree_species) == 1 and named else None
    if T is not None and T not in S:
        triggers.append("T_not_in_S")
    mk = ev["markers"]
    if not mk["aligned"]:
        triggers.append("gtdbtk_not_identified_or_aligned")
    else:
        if mk["unique_markers"] < MIN_UNIQUE_MARKERS:
            triggers.append("low_markers")
        if mk["translation_table"] != "11":
            triggers.append("translation_table_not_11")
    if not ev["verticall_ok"]:
        triggers.append("verticall_pair_missing")
    if any(("species", s) in ev["changed"] for s in S | ({T} if T else set())):
        triggers.append("changed_species_clade")
    for m in METHODS:
        if any(r["species"] != primary[m]["species"] for r in ev["readings"][m].values()):
            triggers.append(f"{m}_rooting_sensitive")
    out = dict(T=T or "", triggers=triggers, species_auto=T is not None and not triggers,
               species_label="", species_label_source="", subspecies_label="", subspecies_label_source="",
               subspecies_reasons=[], top_comparator="", top_comparator_is_reference=False)
    if not out["species_auto"]:
        return out
    out.update(species_label=T, species_label_source="auto_single_species" if S == {T} else "auto_tree_resolved")
    top = ev["top"](T)
    out.update(top_comparator=top["accession"], top_comparator_is_reference=top["is_reference"])
    subs = {primary[m]["subspecies"] for m in METHODS}
    if subs == {""}:
        return out  # species without subspecies in the cohort
    reasons = []
    if len(subs) > 1 or AMBIGUOUS in subs:
        reasons.append("trees_subspecies_disagree_or_ambiguous")
    tree_sub = subs.pop() if len(subs) == 1 else None
    for m in METHODS:
        if any(r["subspecies"] != primary[m]["subspecies"] for r in ev["readings"][m].values()):
            reasons.append(f"{m}_subspecies_rooting_sensitive")
    if tree_sub and tree_sub != AMBIGUOUS and ("subspecies", T, tree_sub) in ev["changed"]:
        reasons.append("changed_subspecies_clade")
    if top["tie"]:
        reasons.append("top_comparator_subspecies_tie")
    elif not top["subspecies"]:
        reasons.append("top_comparator_without_subspecies")
    elif top["subspecies"] != tree_sub:
        reasons.append("top_comparator_subspecies_differs")
    out["subspecies_reasons"] = reasons
    if not reasons:
        out.update(subspecies_label=tree_sub, subspecies_label_source="automatic")
    return out


# ---------------------------------------------------------------------------------------------
# prepare


def _markers(summary_path, table_path, aligned_ids_path, prefix):
    unique = {r["name"].removeprefix(prefix): int(r["number_unique_genes"]) for r in read_tsv(summary_path)}
    with Path(table_path).open() as handle:
        tables = {line.split("\t")[0].removeprefix(prefix): line.rstrip("\n").split("\t")[1] for line in handle if line.strip()}
    aligned = {line.strip().removeprefix(prefix) for line in Path(aligned_ids_path).read_text().splitlines() if line.strip()}
    return unique, tables, aligned


def _verticall_ok(review_path, genomes, expected):
    """Complete: the QC gate's per-genome review shows no missing or duplicated pair and at least the
    expected rows each way. Secondary results add rows, so a complete genome can exceed the count."""
    rows = {r["assembly_id"]: r for r in read_tsv(review_path)}
    ok = {}
    for g in genomes:
        r = rows.get(g)
        if r and r["expected_directed_pairs_each_direction"] != str(expected):
            raise ValueError(f"Verticall review expects {r['expected_directed_pairs_each_direction']} pairs per genome, not {expected}")
        ok[g] = bool(r) and r["zero_result_pairs"] == "0" and r["pairs_with_duplicate_rows"] == "0" and all(
            int(r[f"directed_rows_as_{side}"]) >= expected for side in ("query", "target"))
    return ok


def _packet(accession, code, ev, decision, readings, markers, verticall_ok):
    lines = [f"# Expert labelling packet {code}", "",
             "The NCBI name and the classifier output are withheld. Trees: see `trees/` (validation genomes by coded ID).", "",
             f"**Species passing 95/80 (S):** {', '.join(ev['S']) or 'none'}",
             f"**Unassigned comparators passing (U):** {', '.join(ev['U']) or 'none'}",
             f"**Triggers:** {', '.join(decision['triggers'] + decision['subspecies_reasons']) or 'none'}",
             f"**Automatic species (if any):** {decision['species_label'] or '-'}",
             f"**BAC120 markers:** {markers}; **Verticall pairs complete:** {verticall_ok}", "",
             "| Tree | Rooting | Clade species | Clade subspecies | Cohort tips in clade |", "|---|---|---|---|---:|"]
    for m in METHODS:
        for rooting, r in readings[m].items():
            lines.append(f"| {m} | {rooting} | {r['species']} | {r['subspecies'] or '-'} | {r['clade_size']} |")
    lines += ["", "| Comparator | Species | Subspecies | Package reference | ANI | AF ref | AF query | 95/80 |",
              "|---|---|---|---|---:|---:|---:|---|"]
    for h in ev["hits"][:25]:
        lines.append(f"| {h['accession']} | {h['species'] or '(unassigned)'} | {h['subspecies'] or '-'} | "
                     f"{'yes' if h['reference'] else ''} | {h['ani']:.2f} | {h['af_ref']:.2f} | {h['af_query']:.2f} | "
                     f"{'pass' if h['passes'] else ''} |")
    return "\n".join(lines) + "\n"


def _relabel(tree, names):
    for n in tree.tips():
        tree.label[n] = names.get(tree.label[n], tree.label[n])
    return format_newick(tree, decimals=12)


def prepare(args):
    if args.output.exists():
        raise SystemExit("refusing to overwrite output")
    by_path, by_acc = load_comparators()
    labels = load_tip_labels()
    sample = {r["assembly_unit"]: r for r in read_tsv(SAMPLE)}
    manifest = {r["accession"]: r for r in read_tsv(args.skani_manifest)}
    if set(manifest) != set(sample):
        raise ValueError("skani manifest must cover exactly the sampled genomes")
    evidence = {}
    for acc, m in manifest.items():
        if sha256(m["query_fasta"]) != m["query_sha256"] or sha256(m["hits_tsv"]) != m["hits_sha256"]:
            raise ValueError(f"hash mismatch for {acc}")
        evidence[acc] = ani_evidence(read_tsv(m["hits_tsv"]), by_path, m["query_fasta"])
    newicks = {"bac120": args.bac120_tree.read_text(), "verticall": args.verticall_tree.read_text()}
    papers = {"bac120": args.paper_bac120_tree.read_text(), "verticall": args.paper_verticall_tree.read_text()}
    readings, changed = {m: {} for m in METHODS}, set()
    for m in METHODS:
        for rooting, outgroup in OUTGROUPS.items():
            tree, paper = rooted(newicks[m], outgroup), rooted(papers[m], outgroup)
            tips = {tree.label[n] for n in tree.tips()}
            if tips != set(labels) | set(sample) or {paper.label[n] for n in paper.tips()} != set(labels):
                raise ValueError(f"{m} tree tips do not match the cohort plus the sample")
            readings[m][rooting] = read_clades(tree, labels, sample)
            changed |= changed_labels(paper, tree, labels)
    unique, tables, aligned = _markers(args.marker_summary, args.translation_table, args.aligned_ids, args.gtdbtk_prefix)
    vok = _verticall_ok(args.verticall_assembly_review, sample, len(labels) + len(sample) - 1)
    rows, packets = [], {}
    for acc in sorted(sample):
        ev = evidence[acc]
        rd = {m: {r: readings[m][r][acc] for r in OUTGROUPS} for m in METHODS}
        mk = dict(unique_markers=unique.get(acc, 0), translation_table=tables.get(acc, ""), aligned=acc in aligned)
        d = decide(dict(S=ev["S"], U=ev["U"], readings=rd, markers=mk, verticall_ok=vok[acc], changed=changed,
                        top=lambda sp, ev=ev: top_within_species(ev, sp)))
        code = coded_id(acc)
        species_edge, sub_edge = not d["species_auto"], bool(d["subspecies_reasons"])
        rows.append(dict(coded_id=code, accession=acc, query_sha256=manifest[acc]["query_sha256"],
                         S=";".join(ev["S"]), U=";".join(ev["U"]), near_clone=ev["near_clone"],
                         **{f"{m}_{r}_{k}": rd[m][r][k] for m in METHODS for r in OUTGROUPS for k in ("species", "subspecies")},
                         unique_markers=mk["unique_markers"], translation_table=mk["translation_table"], aligned=mk["aligned"],
                         verticall_ok=vok[acc], T=d["T"], triggers=";".join(d["triggers"]),
                         subspecies_reasons=";".join(d["subspecies_reasons"]),
                         species_label=d["species_label"], species_label_source=d["species_label_source"],
                         subspecies_label=d["subspecies_label"], subspecies_label_source=d["subspecies_label_source"],
                         top_comparator=d["top_comparator"], top_comparator_is_reference=d["top_comparator_is_reference"],
                         top_by_species=json.dumps({sp: top_within_species(ev, sp)
                                                    for sp in sorted({h["species"] for h in ev["hits"] if h["named"]})},
                                                   sort_keys=True),
                         to_michael="species" if species_edge else "subspecies" if sub_edge else ""))
        if species_edge or sub_edge:
            packets[code] = _packet(acc, code, ev, d, rd, mk["unique_markers"], vok[acc])
    args.output.mkdir(parents=True)
    write_tsv(args.output / "preseal-labels.tsv", rows)
    (args.output / "packets").mkdir()
    for code, text in sorted(packets.items()):
        (args.output / "packets" / f"{code}.md").write_text(text)
    (args.output / "trees").mkdir()
    names = {acc: coded_id(acc) for acc in sample} | {t: f"{t}|{lab[0].split()[-1]}{'|' + lab[1] if lab[1] else ''}"
                                                      for t, lab in labels.items() if lab}
    for m in METHODS:
        (args.output / "trees" / f"{m}.{PRIMARY}.coded.nwk").write_text(_relabel(rooted(newicks[m], OUTGROUPS[PRIMARY]), names) + "\n")
    summary = dict(genomes=len(rows), to_michael=dict(Counter(r["to_michael"] or "none" for r in rows)),
                   species_label_source=dict(Counter(r["species_label_source"] or "edge" for r in rows)),
                   triggers=dict(Counter(t for r in rows for t in r["triggers"].split(";") if t)),
                   subspecies_reasons=dict(Counter(t for r in rows for t in r["subspecies_reasons"].split(";") if t)),
                   changed_labels=sorted("|".join(k) for k in changed),
                   inputs={str(p): sha256(p) for p in (COMPARATORS, TREE_TIPS, SAMPLE, args.skani_manifest, args.bac120_tree,
                                                       args.verticall_tree, args.paper_bac120_tree, args.paper_verticall_tree,
                                                       args.marker_summary, args.translation_table, args.aligned_ids,
                                                       args.verticall_assembly_review)},
                   source_sha256=sha256(__file__))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("genomes", "to_michael", "species_label_source", "triggers")}, indent=2))


# ---------------------------------------------------------------------------------------------
# seal


def concept_ids():
    from enterobacter_classifier.reporting import ReportingContext
    refs = ReportingContext(PACKAGE, expected_sha256=PACKAGE_SHA256).evidence_references()
    species = {r["species_name"]: r["species_id"] for r in refs}
    subspecies = {r["subspecies_name"]: r["subspecies_id"] for r in refs if r["subspecies_id"]}
    return species, subspecies


def _edit_distance(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, start=1):
        current = [i]
        for j, y in enumerate(b, start=1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (x != y)))
        previous = current
    return previous[-1]


def _name_without_reference(code, entry, species_ids):
    """A species the package lacks: whitespace normalised; refused if empty, or a package name ignoring case
    or within NEAR_NAME_EDITS edits of one, which is a package species mistyped."""
    name = " ".join(entry.split())
    near = sorted(p for p in species_ids if _edit_distance(name.casefold(), p.casefold()) <= NEAR_NAME_EDITS)
    if not name or near:
        raise ValueError(f"{code}: species_without_reference needs a species name the package lacks"
                         + (f"; {name!r} is within {NEAR_NAME_EDITS} edits of {near}" if near else ""))
    return name


def _michael_subspecies(code, species, entry, subspecies_ids):
    """Michael's subspecies entry for a package species: an epithet, 'none' or 'undetermined' when the
    species has subspecies in the package, otherwise blank."""
    if not any(name.startswith(f"{species} subsp. ") for name in subspecies_ids):
        if entry:
            raise ValueError(f"{code}: {species} has no subspecies in the package; leave subspecies blank")
        return ""
    if entry in (NO_SUBSPECIES, UNDETERMINED) or f"{species} subsp. {entry}" in subspecies_ids:
        return entry
    raise ValueError(f"{code}: unknown subspecies {entry!r} for {species}: use a package epithet, "
                     f"'{NO_SUBSPECIES}' or '{UNDETERMINED}'")


def subspecies_structural(species, subspecies, top_by_species):
    """Whoever labelled: the species' untied top comparator is a package reference of the labelled subspecies,
    so the tool's subspecies can only be the label or absent (section 7.4)."""
    if subspecies in ("", NO_SUBSPECIES, UNDETERMINED):
        return False
    top = top_by_species.get(species)
    return bool(top) and not top["tie"] and top["is_reference"] and top["subspecies"] == subspecies


def seal_rows(preseal, michael, sample, species_ids, subspecies_ids):
    """Merge Michael's labels into the pre-seal rows; return (sealed rows, scorer truth rows)."""
    by_code = {r["coded_id"]: r for r in preseal}
    if set(michael) - set(by_code):
        raise ValueError("labels for unknown coded IDs")
    needed = {c for c, r in by_code.items() if r["to_michael"]}
    if set(michael) != needed:
        raise ValueError(f"Michael's labels must cover exactly the edge cases: missing {sorted(needed - set(michael))}, "
                         f"extra {sorted(set(michael) - needed)}")
    sealed, truth = [], []
    for code, r in sorted(by_code.items()):
        row = dict(r, michael_decision="", michael_confidence="", michael_rationale="")
        species, subspecies = r["species_label"], r["subspecies_label"]
        sp_source, ss_source = r["species_label_source"], r["subspecies_label_source"]
        decision = "species"
        if code in michael:
            m = michael[code]
            decision = m["decision"]
            if decision not in DECISIONS or m["confidence"] not in ("high", "medium", "low") or not m["rationale"].strip():
                raise ValueError(f"incomplete or invalid label for {code}")
            if decision not in ("species", "species_without_reference") and (m["species"] or m["subspecies"]):
                raise ValueError(f"{code}: a {decision} decision takes no species or subspecies")
            if r["to_michael"] == "subspecies":
                if decision != "species" or m["species"] != species:
                    raise ValueError(f"{code}: a subspecies-only edge keeps the automatic species")
                subspecies, ss_source = _michael_subspecies(code, species, m["subspecies"], subspecies_ids), "michael"
            elif decision == "species":
                if not m["species"]:
                    raise ValueError(f"{code}: a species decision needs a species")
                if m["species"] not in species_ids:
                    raise ValueError(f"{code}: {m['species']!r} is not a package species; "
                                     "a species the package lacks is decision species_without_reference")
                species = m["species"]
                subspecies = _michael_subspecies(code, species, m["subspecies"], subspecies_ids)
                sp_source, ss_source = "michael", "michael" if subspecies else ""
            elif decision == "species_without_reference":
                if m["subspecies"]:
                    raise ValueError(f"{code}: species_without_reference takes no subspecies")
                species = _name_without_reference(code, m["species"], species_ids)
                subspecies, sp_source, ss_source = "", "michael", ""
            else:
                species, subspecies, sp_source, ss_source = "", "", "michael", ""
            row.update(michael_decision=decision, michael_confidence=m["confidence"], michael_rationale=m["rationale"])
        elif species not in species_ids:
            raise ValueError(f"{code}: automatic species {species!r} is not a package species")
        S = [s for s in r["S"].split(";") if s]
        if decision == "species":
            kind = "named"
        elif decision in ("species_without_reference", "novel_lineage"):
            kind = "no_reference"
        else:
            kind = decision
        sub_name = f"{species} subsp. {subspecies}" if subspecies not in ("", NO_SUBSPECIES, UNDETERMINED) else ""
        if kind == "named" and sub_name and sub_name not in subspecies_ids:
            raise ValueError(f"{code}: unknown subspecies {sub_name}")
        structural = kind == "named" and subspecies_structural(species, subspecies, json.loads(r["top_by_species"]))
        row.update(species_label=species, subspecies_label=subspecies, species_label_source=sp_source,
                   subspecies_label_source=ss_source, truth_kind=kind, subspecies_structural=structural)
        sealed.append(row)
        s = sample[r["accession"]]
        truth.append(dict(genome_id=r["accession"], layer="expert", truth_kind=kind,
                          label_source=sp_source if kind == "named" else "michael",
                          truth_species_id=species_ids[species] if kind == "named" else "",
                          truth_subspecies_id=(UNDETERMINED if subspecies == UNDETERMINED else subspecies_ids.get(sub_name, ""))
                          if kind == "named" else "",
                          ani_species_passing=";".join(sorted(species_ids[x] for x in S)),
                          subspecies_structural=str(structural),
                          subspecies_label_source=ss_source if kind == "named" else "",
                          expected_species_ids="", expected_subspecies_id="", block=s["component_id"],
                          component_id=s["component_id"], query_sha256=r["query_sha256"], sample_group=s["sample_group"],
                          qc_stratum=s["qc_stratum"], near_clone=r["near_clone"], coded_id=code))
    return sealed, truth


def seal(args):
    if args.output.exists():
        raise SystemExit("refusing to overwrite output")
    preseal = read_tsv(args.preseal)
    michael = {r["coded_id"]: r for r in read_tsv(args.michael)}
    sample = {r["assembly_unit"]: r for r in read_tsv(SAMPLE)}
    species_ids, subspecies_ids = concept_ids()
    sealed, truth = seal_rows(preseal, michael, sample, species_ids, subspecies_ids)
    args.output.mkdir(parents=True)
    write_tsv(args.output / "sealed-labels.tsv", sealed)
    write_tsv(args.output / "expert-truth.tsv", truth)
    record = dict(sealed_labels_sha256=sha256(args.output / "sealed-labels.tsv"),
                  expert_truth_sha256=sha256(args.output / "expert-truth.tsv"),
                  preseal_sha256=sha256(args.preseal), michael_labels_sha256=sha256(args.michael), source_sha256=sha256(__file__))
    (args.output / "seal.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2))


# ---------------------------------------------------------------------------------------------
# loo: in-sample leave-one-out on the paper trees


def loo(args):
    if args.output.exists():
        raise SystemExit("refusing to overwrite output")
    labels = load_tip_labels()
    named = {t: lab for t, lab in labels.items() if lab}
    passing, best = defaultdict(set), defaultdict(list)
    with args.pairs.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["status"] != "measured" or r["a"] == r["b"]:
                continue
            for x, y in ((r["a"], r["b"]), (r["b"], r["a"])):
                if x in named and y in named:
                    best[x].append(dict(accession=y, named=True, species=named[y][0], subspecies=named[y][1],
                                        reference=False, ani=float(r["ani"]), af_ref=0.0, af_query=0.0, passes=False))
                    if r["passes_95_80"] == "True":
                        passing[x].add(named[y][0])
    readings = {m: {} for m in METHODS}
    for m, path in (("bac120", args.paper_bac120_tree), ("verticall", args.paper_verticall_tree)):
        newick = path.read_text()
        for rooting, outgroup in OUTGROUPS.items():
            tree = rooted(newick, outgroup)
            readings[m][rooting] = read_clades(tree, labels, named)
    rows = []
    for t, (sp, ss) in sorted(named.items()):
        hits = sorted(best[t], key=lambda h: (-h["ani"], h["accession"]))
        ev = dict(hits=hits)
        d = decide(dict(S=sorted(passing[t]), U=[], changed=set(), verticall_ok=True,
                        markers=dict(unique_markers=120, translation_table="11", aligned=True),
                        readings={m: {r: readings[m][r][t] for r in OUTGROUPS} for m in METHODS},
                        top=lambda s, ev=ev: top_within_species(ev, s)))
        species_outcome = ("edge" if not d["species_auto"] else
                           f"{d['species_label_source']}_{'match' if d['species_label'] == sp else 'differs'}")
        sub_outcome = ("" if not d["species_auto"] or not ss else
                       "edge" if d["subspecies_reasons"] else "match" if d["subspecies_label"] == ss else "differs")
        rows.append(dict(tip=t, species=sp, subspecies=ss, S=";".join(sorted(passing[t])), T=d["T"],
                         triggers=";".join(d["triggers"]), species_outcome=species_outcome,
                         subspecies_reasons=";".join(d["subspecies_reasons"]), subspecies_outcome=sub_outcome))
    args.output.mkdir(parents=True)
    write_tsv(args.output / "loo.tsv", rows)
    summary = dict(tips=len(rows), species_outcomes=dict(Counter(r["species_outcome"] for r in rows)),
                   subspecies_outcomes=dict(Counter(r["subspecies_outcome"] for r in rows if r["subspecies_outcome"])),
                   triggers=dict(Counter(x for r in rows for x in r["triggers"].split(";") if x)),
                   subspecies_reasons=dict(Counter(x for r in rows for x in r["subspecies_reasons"].split(";") if x)),
                   note="in-sample and partly circular: the cohort labels drew on these trees and pairs",
                   inputs={str(p): sha256(p) for p in (TREE_TIPS, args.pairs, args.paper_bac120_tree, args.paper_verticall_tree)},
                   source_sha256=sha256(__file__))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("tips", "species_outcomes", "subspecies_outcomes", "triggers")}, indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--skani-manifest", type=Path, required=True,
                   help="accession, query_fasta, query_sha256, hits_tsv, hits_sha256 (query = the downloaded FASTA)")
    for name in ("bac120-tree", "verticall-tree", "paper-bac120-tree", "paper-verticall-tree", "marker-summary",
                 "translation-table", "aligned-ids", "verticall-assembly-review"):
        p.add_argument(f"--{name}", type=Path, required=True)
    p.add_argument("--gtdbtk-prefix", default="")
    p.add_argument("--output", type=Path, required=True)
    s = sub.add_parser("seal")
    s.add_argument("--preseal", type=Path, required=True)
    s.add_argument("--michael", type=Path, required=True,
                   help="coded_id, decision (species|species_without_reference|unresolved|outside_genus|novel_lineage), species, "
                        "subspecies (package epithet, 'none' or 'undetermined' for a species with subspecies; else blank), "
                        "confidence, rationale")
    s.add_argument("--output", type=Path, required=True)
    o = sub.add_parser("loo")
    o.add_argument("--pairs", type=Path, required=True)
    o.add_argument("--paper-bac120-tree", type=Path, required=True)
    o.add_argument("--paper-verticall-tree", type=Path, required=True)
    o.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    {"prepare": prepare, "seal": seal, "loo": loo}[args.command](args)


if __name__ == "__main__":
    main()
