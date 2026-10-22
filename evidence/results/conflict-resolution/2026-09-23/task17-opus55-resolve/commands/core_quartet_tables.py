"""Task 17: the core-genome quartet site patterns and gene quartets quoted in README §7.2, §7.4, §7.5.

These were first computed interactively (the README records that); this script regenerates every
one of them into two tables with one fixed bootstrap seed, from the frozen panel concatenations
and the panel gene trees:

  tables/core-site-patterns.tsv   informative-column shares per pairing, 95% gene-block bootstrap
                                  interval (1,000 replicates, seed 20260923), Patterson's D and z
  tables/core-gene-quartets.tsv   gene trees restricted to the four genomes: resolved genes (induced
                                  branch >= 1 expected substitution) per pairing, and unresolved

Run from the project root:  python3 <run>/commands/core_quartet_tables.py
"""
import csv
import random
import sys
from pathlib import Path

sys.path.insert(0, "scripts")
from summarise_concordance import NT, gene_quartets, read_alignment, read_loci, site_patterns  # noqa: E402

RUN = Path("results/conflict-resolution/2026-09-23/task17-opus55-resolve")
SEED, REPS = 20260923, 1000
N = {
    "GCA_001750725.1": "lud", "GCA_000025565.1": "clo", "GCA_002939185.1": "sich",
    "GCA_025002605.1": "sich-alt", "GCA_034008295.1": "chin", "GCA_001521715.1": "asb",
    "GCA_008364625.2": "dyk", "GCA_001729765.1": "kob", "GCA_001729805.1": "rog",
    "GCA_019048125.1": "Kaer", "GCA_000982825.1": "Csak", "GCA_000211415.1": "mori",
    "GCA_003594935.2": "huax", "GCA_019048245.1": "horm", "GCA_028890245.1": "past",
    "GCA_019047785.1": "canc", "GCA_001654845.1": "soli", "GCA_001729745.1": "hof",
    "GCA_040096145.1": "int", "GCA_052922205.1": "not", "GCA_004331385.1": "qhorm",
    "GCA_001729705.1": "ohar",
}
ID = {v: k for k, v in N.items()}
C1_KAER = [["lud", "clo", x, "Kaer"] for x in ("sich", "asb", "kob")]
C3 = [["hof", "int", "not", o] for o in ("horm", "qhorm", "ohar")]
SITES = {
    ("C1-ludwigii", "ledger-panel/alignment.fna"): C1_KAER,
    ("C2-sichuanensis", "ledger-panel/alignment.fna"): [
        ["sich", "chin", "clo", "Kaer"], ["sich-alt", "chin", "clo", "Kaer"], ["sich", "chin", "clo", "lud"],
        ["sich", "clo", "rog", "lud"], ["sich", "chin", "asb", "clo"], ["sich", "chin", "asb", "lud"],
        ["sich", "chin", "dyk", "clo"], ["sich", "sich-alt", "chin", "clo"]],
    ("C3-trio", "ledger-panel/alignment.fna"): C3,
    ("C3-trio", "hof-masked/alignment.fna"): C3,
}
GENES = {
    ("C1-ludwigii", "concordance-panel/loci.treefile"): C1_KAER,
    ("C2-sichuanensis", "concordance-panel/loci.treefile"): [
        ["sich", "chin", "clo", "Kaer"], ["sich", "chin", "clo", "lud"], ["sich", "chin", "asb", "clo"],
        ["sich-alt", "chin", "asb", "clo"]],
    ("C3-trio", "concordance-panel/loci.treefile"): C3,
    ("C12x-extended-sampling", "gene-trees-panel/loci.treefile"):
        [["lud", "clo", x, o] for x in ("asb", "kob", "sich")
         for o in ("horm", "mori", "huax", "past", "canc", "soli", "Kaer")]
        + [["sich", "chin", x, o] for x in ("clo", "asb") for o in ("horm", "mori", "huax")],
}


def pairings(q):
    a, b, c, d = q
    return [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]


def key(p):
    return "|".join(",".join(ID[x] for x in side) for side in p)


def label(p, fourth):
    side = p[1] if fourth in p[0] else p[0]
    return f"{side[0]}+{side[1]}"


def main():
    (RUN / "tables").mkdir(exist_ok=True)
    rows = []
    for (conflict, rel), quartets in SITES.items():
        core = RUN / conflict / "core-genome"
        aln = read_alignment(core / rel)
        loci = read_loci(core / "ledger-panel/concatenation-genes.nex")
        for q in quartets:
            p = site_patterns(aln, [ID[x] for x in q], loci, NT, random.Random(SEED), reps=REPS)
            row = {"conflict": conflict, "alignment": rel, "quartet": ",".join(q[:3]), "fourth": q[3],
                   "informative_columns": p["informative_columns"]}
            for i, pr in enumerate(pairings(q), 1):
                v = p[key(pr)] if key(pr) in p else p["|".join(reversed(key(pr).split("|")))]
                row[f"pairing{i}"] = label(pr, q[3])
                row[f"share{i}"] = v["share"]
                row[f"ci95_{i}"] = f"{v['bootstrap_95'][0]}-{v['bootstrap_95'][1]}"
            row["patterson_D"] = p["patterson_d"]["D"]
            row["D_z"] = p["patterson_d"]["z"]
            rows.append(row)
            print(conflict, rel, q, row["share1"], row["share2"], row["share3"], flush=True)
    write(RUN / "tables/core-site-patterns.tsv", rows)

    rows = []
    for (conflict, rel), quartets in GENES.items():
        core = RUN / conflict / "core-genome"
        loci = read_loci(core / "ledger-panel/loci.nex")
        trees = [line.strip() for line in (core / rel).read_text().splitlines() if line.strip()]
        assert len(loci) == len(trees), (conflict, len(loci), len(trees))
        for q in quartets:
            g = gene_quartets([ID[x] for x in q], loci, trees)
            row = {"conflict": conflict, "gene_trees": rel, "quartet": ",".join(q[:3]), "fourth": q[3],
                   "genes": g["genes"], "resolved": g["resolved"], "unresolved": g["unresolved"]}
            for i, pr in enumerate(pairings(q), 1):
                k = key(pr)
                n = g.get(k, g.get("|".join(reversed(k.split("|"))), 0))
                row[f"pairing{i}"] = label(pr, q[3])
                row[f"genes{i}"] = n
                row[f"share{i}"] = round(n / g["resolved"], 3) if g["resolved"] else ""
            rows.append(row)
            print(conflict, q, row["resolved"], row["share1"], row["share2"], row["share3"], flush=True)
    write(RUN / "tables/core-gene-quartets.tsv", rows)


def write(path, rows):
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
