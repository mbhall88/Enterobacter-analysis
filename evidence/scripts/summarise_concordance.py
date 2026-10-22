#!/usr/bin/env python3
"""Read the named branches' concordance, gene by gene and site by site, for one Task 17 panel.

Four new statistics, none of them a support value (SET-18 rule 2), and none ever tabulated on one
scale with FastTree's SH-like support, IQ-TREE's ultrafast bootstrap or a Mash seed frequency:

**gCF / sCF on the named branches.**  For every hypothesis in ``hypotheses.json`` and every
reference tree from ``run_concordance.sh``, the branch whose bipartition equals the hypothesis
(ignoring the hypothesis's free tips) is found by its tip set, not by its position, and IQ-TREE's
gCF, gDF1, gDF2, gDFP and likelihood sCF, sDF1, sDF2 are read for it.  gCF is read twice: on the
gene trees with branches under one substitution collapsed (``collapsed_*``, the one to quote) and
as written (``raw_*``).  sCF is reported per seed
and as a range, so its Monte-Carlo stability is visible.

**Gene support per hypothesis.**  Each gene tree is restricted to the hypothesis's taxa and read
in the three states SET-18 uses for trees: *present* (the gene tree has the bipartition),
*contradicted* (it has an incompatible split on a branch of at least one expected substitution
across the gene -- the one-unit-of-evidence rule of SET-18 rule 10, in the gene's own units),
*unresolved* (neither).  A gene with fewer than two taxa on either side is *not decisive*.

**Gene quartets.**  For the same named quartets, each gene tree restricted to the four genomes and
read as one of the three pairings (on an induced branch of at least one substitution) or as
unresolved -- the gene-level counterpart of the site patterns.

**Quartet site patterns.**  For named quartets of single genomes, the counts of alignment columns
supporting each of the three pairings (both states present twice, all four called), with a
gene-block bootstrap.  With the fourth genome as an outgroup, the two minority patterns are the
ABBA and BABA of Patterson's D: equal under incomplete lineage sorting alone, unequal under
introgression.

**Position (line 6).**  Genes that decide between two competing hypotheses, ordered along a
complete reference chromosome: the number of runs of same-state genes against its permutation
distribution.  Fewer runs than expected is clustering by position.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import random
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import format_newick, parse_newick  # noqa: E402

NT = set("ACGT")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader((l for l in handle if not l.startswith("#")), delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------------------------
# Splits of a tree
# --------------------------------------------------------------------------------------------

def tree_splits(newick: str) -> tuple[set[str], list[tuple[frozenset[str], float, str]]]:
    """Tips, and every internal edge as (tips below, length, label) -- rooted as written."""
    tree = parse_newick(newick)
    below: dict[int, frozenset[str]] = {}
    for node in tree.postorder():
        if not tree.children[node]:
            below[node] = frozenset([tree.label[node]])
        else:
            below[node] = frozenset().union(*(below[c] for c in tree.children[node]))
    tips = set(below[0])
    out = []
    for node in range(tree.size):
        if node == 0 or not tree.children[node]:
            continue
        out.append((below[node], tree.blen[node] or 0.0, tree.label[node]))
    return tips, out


def normalise(side: frozenset[str], tips: frozenset[str]) -> frozenset[str]:
    """One canonical side of a bipartition: the side not holding the smallest tip."""
    anchor = min(tips)
    return side if anchor not in side else tips - side


def compatible(a: frozenset[str], b: frozenset[str], tips: frozenset[str]) -> bool:
    a2, b2 = tips - a, tips - b
    return not (a & b) or not (a & b2) or not (a2 & b) or not (a2 & b2)


def find_branch(hyp_side: set[str], free: set[str], splits, tips: set[str]):
    """The edge whose bipartition, ignoring free tips, equals the hypothesis (exact match first)."""
    core = frozenset(tips - free)
    target = normalise(frozenset(hyp_side) & core, core)
    candidates = []
    for side, length, label in splits:
        reduced = normalise(side & core, core)
        if reduced == target:
            candidates.append((len(side - core), side, length, label))
    if not candidates:
        return None
    candidates.sort(key=lambda c: c[0])
    return candidates[0]


# --------------------------------------------------------------------------------------------
# IQ-TREE concordance output
# --------------------------------------------------------------------------------------------

def cf_table(stat: Path) -> dict[str, dict[str, str]]:
    return {r["ID"]: r for r in read_tsv(stat)}


def pair_branch(pair: set[str], excluded: set[str], splits, tips: set[str]):
    """For a quartet hypothesis: the edge uniting the pair, away from the third type and outgroups.

    The smallest bipartition side holding both paired types and none of ``excluded``; its edge is
    the named branch whose concordance factors bear on the pairing.  None if the tree does not
    separate the pair from the others.
    """
    universe = frozenset(tips)
    best = None
    for side, length, label in splits:
        for candidate in (side, universe - side):
            if pair <= candidate and not (excluded & candidate):
                if best is None or len(candidate) < len(best[1]):
                    best = (0, candidate, length, label)
    return best


def named_branch_rows(hypotheses: dict, scope: str, conc: Path, references: list[str],
                      seeds: list[str], outgroups: list[str] = ()) -> list[dict]:
    rows = []
    for ref in references:
        branch_tree = (conc / f"{ref}.gcf.cf.branch").read_text().strip()
        tips, splits = tree_splits(branch_tree)
        gcf = cf_table(conc / f"{ref}.gcf.cf.stat")
        gcfc_path = conc / f"{ref}.gcfc.cf.stat"
        gcfc = cf_table(gcfc_path) if gcfc_path.exists() else {}
        scf = {s: cf_table(conc / f"{ref}.scf-{s}.cf.stat") for s in seeds}
        scf_tips, scf_splits = tree_splits((conc / f"{ref}.scf-{seeds[0]}.cf.branch").read_text().strip())
        for name, hyp in hypotheses.items():
            side = set(hyp[scope]["side_a"])
            free = set(hyp.get("free_in_constraint", []))
            row = {"reference": ref, "hypothesis": name, "kind": hyp["kind"]}
            if hyp["kind"] == "quartet":
                found = pair_branch(side, set(hyp[scope]["side_b_named"]) | set(outgroups),
                                    splits, tips)
                found_s = pair_branch(side, set(hyp[scope]["side_b_named"]) | set(outgroups),
                                      scf_splits, scf_tips)
            else:
                found = find_branch(side, free, splits, tips)
                found_s = find_branch(side, free, scf_splits, scf_tips)
            if not found:
                row["branch_in_reference"] = False
                rows.append(row)
                continue
            _, _, length, branch_id = found
            sid = found_s[3] if found_s else None
            g = gcf.get(branch_id, {})
            gc = gcfc.get(branch_id, {})
            keys = ("gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN")
            row.update({"branch_in_reference": True, "branch_id": branch_id,
                        "branch_length": length,
                        **{f"collapsed_{k}": gc.get(k, "") for k in keys},
                        **{f"raw_{k}": g.get(k, "") for k in keys}})
            values = {k: [] for k in ("sCF", "sDF1", "sDF2", "sN")}
            for s in seeds:
                entry = scf[s].get(sid, {})
                for k in values:
                    if entry.get(k, "") != "":
                        values[k].append(float(entry[k]))
                row[f"sCF_seed_{s}"] = entry.get("sCF", "")
            for k, vs in values.items():
                if vs:
                    row[f"{k}_mean"] = round(statistics.fmean(vs), 2)
                    row[f"{k}_range"] = round(max(vs) - min(vs), 2)
            rows.append(row)
    return rows


# --------------------------------------------------------------------------------------------
# Gene support
# --------------------------------------------------------------------------------------------

def read_loci(nexus: Path) -> list[tuple[str, int, int]]:
    return [(m.group(1), int(m.group(2)), int(m.group(3))) for m in
            re.finditer(r"charset\s+(\S+)\s*=\s*(\d+)-(\d+)\s*;", nexus.read_text())]


def restricted_lengths(universe: frozenset[str], splits) -> dict[frozenset[str], float]:
    """Each bipartition of the tree induced on ``universe``, with its induced branch length.

    Several edges of the full tree can induce the same bipartition on a subset of its tips (a path
    whose side branches hold no universe tip); the induced edge's length is their sum.
    """
    out: dict[frozenset[str], float] = {}
    for side, length, _ in splits:
        reduced = side & universe
        if len(reduced) < 2 or len(universe - reduced) < 2:
            continue
        key = normalise(reduced, universe)
        out[key] = out.get(key, 0.0) + length
    return out


def classify(side_a: set[str], side_b: set[str], tips: set[str], splits, min_length: float) -> str:
    """present / contradicted / unresolved / not_decisive, both calls needing one substitution.

    A bipartition counts as *present* only on a branch of at least ``min_length``: gene trees from
    near-identical sequences resolve zero-length branches arbitrarily, and counting those would
    manufacture support for whichever resolution the search happened to write.
    """
    a, b = side_a & tips, side_b & tips
    if len(a) < 2 or len(b) < 2:
        return "not_decisive"
    universe = frozenset(a | b)
    target = normalise(frozenset(a), universe)
    lengths = restricted_lengths(universe, splits)
    if lengths.get(target, 0.0) >= min_length:
        return "present"
    for key, length in lengths.items():
        if length >= min_length and not compatible(key, target, universe):
            return "contradicted"
    return "unresolved"


def collapse_short(newick: str, min_length: float) -> str:
    """The gene tree with every internal branch shorter than ``min_length`` collapsed."""
    tree = parse_newick(newick)
    changed = True
    while changed:
        changed = False
        for node in range(1, tree.size):
            if tree.parent[node] == -1 or not tree.children[node]:
                continue
            if (tree.blen[node] or 0.0) < min_length:
                parent = tree.parent[node]
                for child in tree.children[node]:
                    tree.parent[child] = parent
                tree.children[parent] = [c for c in tree.children[parent] if c != node] + \
                    tree.children[node]
                tree.children[node] = []
                tree.parent[node] = -1
                changed = True
    return format_newick(tree, decimals=12)


def gene_support(hypotheses: dict, scope: str, loci: list[tuple[str, int, int]], trees: list[str],
                 outgroups: list[str]) -> list[dict]:
    if len(loci) != len(trees):
        raise SystemExit(f"{len(loci)} loci but {len(trees)} gene trees")
    rows = []
    for (locus, start, end), newick in zip(loci, trees):
        tips, splits = tree_splits(newick)
        min_length = 1.0 / (end - start + 1)
        row = {"locus": locus, "columns": end - start + 1, "tips": len(tips)}
        for name, hyp in hypotheses.items():
            free = set(hyp.get("free_in_constraint", []))
            bip = hyp[scope]
            if hyp["kind"] == "quartet":
                a = set(bip["side_a"])
                b = set(bip["side_b_named"]) | set(outgroups)
                if not set(outgroups) & tips:
                    row[name] = "not_decisive"
                    continue
            else:
                a, b = set(bip["side_a"]) - free, set(bip["side_b"]) - free
            row[name] = classify(a, b, tips, splits, min_length)
        rows.append(row)
    return rows


def gene_quartets(quartet: list[str], loci, trees: list[str]) -> dict:
    """Per gene, which pairing of four named genomes the gene tree carries.

    The gene tree is restricted to the four genomes; a pairing counts only on an induced branch of
    at least one expected substitution across the gene, otherwise the gene is *unresolved* for
    this quartet.  A gene missing any of the four is *not decisive*.  This reads each gene at the
    level of single named genomes, where lineage-level clades are broken by one misplaced
    representative.
    """
    a, b, c, d = quartet
    labels = {frozenset((a, b)): f"{a},{b}|{c},{d}", frozenset((a, c)): f"{a},{c}|{b},{d}",
              frozenset((a, d)): f"{a},{d}|{b},{c}"}
    counts = {v: 0 for v in labels.values()}
    counts.update({"unresolved": 0, "not_decisive": 0})
    per_gene = []
    universe = frozenset(quartet)
    for (locus, start, end), newick in zip(loci, trees):
        tips, splits = tree_splits(newick)
        if not universe <= tips:
            counts["not_decisive"] += 1
            per_gene.append((locus, "not_decisive"))
            continue
        lengths = restricted_lengths(universe, splits)
        state = "unresolved"
        for key, length in lengths.items():
            if length >= 1.0 / (end - start + 1):
                pair = key if a in key else universe - key
                state = labels[frozenset(pair)]
        counts[state] += 1
        per_gene.append((locus, state))
    decisive = sum(v for k, v in counts.items() if "|" in k)
    return {"quartet": quartet, "genes": len(loci), "resolved": decisive, **counts,
            "shares_of_resolved": {k: round(v / decisive, 4) for k, v in counts.items()
                                   if "|" in k and decisive},
            "per_gene": per_gene}


def support_summary(rows: list[dict], hypotheses: dict) -> list[dict]:
    out = []
    for name in hypotheses:
        states = [r[name] for r in rows]
        decisive = [s for s in states if s != "not_decisive"]
        entry = {"hypothesis": name, "genes": len(states), "decisive": len(decisive)}
        for s in ("present", "contradicted", "unresolved"):
            n = decisive.count(s)
            entry[s] = n
            entry[f"{s}_fraction_of_decisive"] = round(n / len(decisive), 4) if decisive else ""
        out.append(entry)
    return out


def exclusive_pairs(rows: list[dict], hypotheses: dict) -> list[dict]:
    """For each pair of hypotheses: genes presenting one and not the other."""
    out = []
    for h1, h2 in itertools.combinations(hypotheses, 2):
        only1 = sum(1 for r in rows if r[h1] == "present" and r[h2] != "present")
        only2 = sum(1 for r in rows if r[h2] == "present" and r[h1] != "present")
        both = sum(1 for r in rows if r[h1] == "present" and r[h2] == "present")
        out.append({"hypothesis_1": h1, "hypothesis_2": h2, "genes_present_1_only": only1,
                    "genes_present_2_only": only2, "genes_present_both": both,
                    "share_1_of_deciding": round(only1 / (only1 + only2), 4) if only1 + only2 else ""})
    return out


# --------------------------------------------------------------------------------------------
# Site patterns
# --------------------------------------------------------------------------------------------

def read_alignment(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name:
                out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def site_patterns(aln: dict[str, str], quartet: list[str], loci, alphabet: set[str],
                  rng: random.Random, reps: int = 1000) -> dict:
    """Counts of AB|CD, AC|BD, AD|BC columns per locus, totals, and a gene-block bootstrap."""
    a, b, c, d = (aln[t] for t in quartet)
    per_locus = []
    for _locus, start, end in loci:
        counts = [0, 0, 0]
        for i in range(start - 1, end):
            x = (a[i].upper(), b[i].upper(), c[i].upper(), d[i].upper())
            if not all(ch in alphabet for ch in x) or len(set(x)) != 2:
                continue
            if x[0] == x[1] and x[2] == x[3]:
                counts[0] += 1
            elif x[0] == x[2] and x[1] == x[3]:
                counts[1] += 1
            elif x[0] == x[3] and x[1] == x[2]:
                counts[2] += 1
        per_locus.append(counts)
    total = [sum(c[k] for c in per_locus) for k in range(3)]
    boot = []
    for _ in range(reps):
        sample = [per_locus[rng.randrange(len(per_locus))] for _ in per_locus]
        boot.append([sum(c[k] for c in sample) for k in range(3)])
    labels = [f"{quartet[0]},{quartet[1]}|{quartet[2]},{quartet[3]}",
              f"{quartet[0]},{quartet[2]}|{quartet[1]},{quartet[3]}",
              f"{quartet[0]},{quartet[3]}|{quartet[1]},{quartet[2]}"]
    out = {"quartet": quartet, "loci": len(per_locus), "informative_columns": sum(total)}
    shares = []
    for k in range(3):
        frac = [bt[k] / sum(bt) for bt in boot if sum(bt)]
        frac.sort()
        share = total[k] / sum(total) if sum(total) else None
        shares.append(share)
        out[labels[k]] = {"columns": total[k], "share": round(share, 4) if share is not None else None,
                          "bootstrap_95": [round(frac[int(0.025 * len(frac))], 4),
                                           round(frac[int(0.975 * len(frac)) - 1], 4)] if frac else None}
    # Patterson's D with the fourth genome as outgroup: minority patterns relative to the majority.
    order = sorted(range(3), key=lambda k: -total[k])
    m1, m2 = order[1], order[2]
    if total[m1] + total[m2]:
        dval = (total[m1] - total[m2]) / (total[m1] + total[m2])
        dboot = [(bt[m1] - bt[m2]) / (bt[m1] + bt[m2]) for bt in boot if bt[m1] + bt[m2]]
        sd = statistics.pstdev(dboot) if len(dboot) > 1 else 0.0
        out["majority_pattern"] = labels[order[0]]
        out["patterson_d"] = {"minority_a": labels[m1], "minority_b": labels[m2],
                              "D": round(dval, 4), "bootstrap_sd": round(sd, 4),
                              "z": round(dval / sd, 2) if sd else None}
    return out


# --------------------------------------------------------------------------------------------
# Position
# --------------------------------------------------------------------------------------------

def runs_test(states: list[str], rng: random.Random, reps: int = 10000) -> dict:
    def count_runs(seq):
        return 1 + sum(1 for x, y in zip(seq, seq[1:]) if x != y) if seq else 0
    observed = count_runs(states)
    null = []
    work = list(states)
    for _ in range(reps):
        rng.shuffle(work)
        null.append(count_runs(work))
    below = sum(1 for v in null if v <= observed)
    return {"genes": len(states), "runs": observed, "null_median_runs": statistics.median(null),
            "p_clustered": round((below + 1) / (reps + 1), 4)}


def positional_rows(support: list[dict], positions: Path, h1: str, h2: str,
                    rng: random.Random) -> dict:
    pos = {r["locus"]: int(r["genome_position"]) for r in read_tsv(positions)
           if r["placed"] == "True"}
    ordered = []
    for row in support:
        if row["locus"] not in pos:
            continue
        s1, s2 = row[h1] == "present", row[h2] == "present"
        if s1 != s2:
            ordered.append((pos[row["locus"]], h1 if s1 else h2))
    ordered.sort()
    states = [s for _, s in ordered]
    out = {"hypothesis_1": h1, "hypothesis_2": h2, "genes_placed_and_deciding": len(states),
           "genes_for_1": states.count(h1), "genes_for_2": states.count(h2)}
    if len(set(states)) == 2:
        out["runs_test"] = runs_test(states, rng)
    out["ordered"] = [{"genome_position": p, "supports": s} for p, s in ordered]
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hypotheses", type=Path, required=True)
    parser.add_argument("--scope", default="panel", choices=["panel", "subset"])
    parser.add_argument("--concordance", type=Path, required=True, help="run_concordance.sh dir")
    parser.add_argument("--references", nargs="+", required=True)
    parser.add_argument("--seeds", nargs="+", default=["20260923", "20260924", "20260925"])
    parser.add_argument("--loci", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--seqtype", choices=["AA", "DNA"], required=True)
    parser.add_argument("--outgroups", nargs="+", required=True)
    parser.add_argument("--quartet", action="append", default=[],
                        help="four comma-separated tips; the fourth is the outgroup for D")
    parser.add_argument("--positions", type=Path)
    parser.add_argument("--positional-pair", nargs=2, metavar=("H1", "H2"), action="append",
                        default=[])
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    hypotheses = json.loads(args.hypotheses.read_text())

    named = named_branch_rows(hypotheses, args.scope, args.concordance, args.references, args.seeds,
                              args.outgroups)
    write_tsv(args.output / "named-branches.tsv", named)

    loci = read_loci(args.loci)
    trees = [l.strip() for l in (args.concordance / "loci.treefile").read_text().splitlines()
             if l.strip()]
    aln = read_alignment(args.alignment)
    alphabet = NT if args.seqtype == "DNA" else AA
    for (locus, start, end), newick in zip(loci, trees):
        tips, _ = tree_splits(newick)
        with_data = {t for t, s in aln.items() if any(ch.upper() in alphabet for ch in s[start - 1:end])}
        if tips != with_data:
            raise SystemExit(f"gene tree for {locus} has {len(tips)} tips, the locus has data for "
                             f"{len(with_data)}: trees are not in partition order")
    support = gene_support(hypotheses, args.scope, loci, trees, args.outgroups)
    write_tsv(args.output / "gene-support.tsv", support)
    summary = {"support": support_summary(support, hypotheses),
               "exclusive_pairs": exclusive_pairs(support, hypotheses)}
    write_tsv(args.output / "gene-support-summary.tsv", summary["support"])
    write_tsv(args.output / "gene-support-pairs.tsv", summary["exclusive_pairs"])

    patterns = [site_patterns(aln, q.split(","), loci, alphabet, rng) for q in args.quartet]
    (args.output / "site-patterns.json").write_text(json.dumps(patterns, indent=1) + "\n")
    gq = [gene_quartets(q.split(","), loci, trees) for q in args.quartet]
    (args.output / "gene-quartets.json").write_text(json.dumps(gq, indent=1) + "\n")

    positional = []
    if args.positions:
        for h1, h2 in args.positional_pair:
            positional.append(positional_rows(support, args.positions, h1, h2, rng))
    (args.output / "positional.json").write_text(json.dumps(positional, indent=1) + "\n")
    record = {"hypotheses": list(hypotheses), "references": args.references, "loci": len(loci),
              "gene_trees_checked_against_partition_order": len(trees), "seed": args.seed}
    (args.output / "record.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({"support": summary["support"], "pairs": summary["exclusive_pairs"],
                      "patterns": [{k: v for k, v in p.items() if k != "quartet"} for p in patterns],
                      "gene_quartets": [{k: v for k, v in g.items() if k != "per_gene"} for g in gq],
                      "positional": [{k: v for k, v in p.items() if k != "ordered"} for p in positional]},
                     indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
