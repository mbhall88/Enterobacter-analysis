"""Task 17: do C1's competing gene-quartet histories form tracts along the ludwigii chromosome?

C12x panel, final gene trees (gene-trees-panel/loci.treefile, written 21:00 on 23 Sep when its
job finished).  A first version of this analysis read that file at 20:25, while the job was still
running, and is kept only as superseded/positional-gene-quartets.from-intermediate-gene-trees.json.

For each quartet {ludwigii, cloacae, X, fourth} with X a 288-group type and the fourth an
in-genus reference lineage, every gene tree is read as one of the three pairings (induced branch
>= 1 expected substitution) or unresolved; the genes deciding between two pairings are ordered
along the ludwigii type chromosome and their runs counted against 10,000 shuffles.

Run from the project root; writes C12x-extended-sampling/core-genome/positional-gene-quartets.json
"""
import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, "scripts")
from summarise_concordance import gene_quartets, read_loci, runs_test  # noqa: E402

CORE = Path("results/conflict-resolution/2026-09-23/task17-opus55-resolve/C12x-extended-sampling/core-genome")
LUD, CLO = "GCA_001750725.1", "GCA_000025565.1"
THIRDS = {"asb": "GCA_001521715.1", "kob": "GCA_001729765.1", "sich": "GCA_002939185.1"}
FOURTHS = {"horm": "GCA_019048245.1", "mori": "GCA_000211415.1", "huax": "GCA_003594935.2",
           "past": "GCA_028890245.1", "canc": "GCA_019047785.1", "soli": "GCA_001654845.1"}


def main():
    loci = read_loci(CORE / "ledger-panel/loci.nex")
    trees = [line.strip() for line in (CORE / "gene-trees-panel/loci.treefile").read_text().splitlines()
             if line.strip()]
    assert len(loci) == len(trees)
    pos = {r["locus"]: int(r["genome_position"])
           for r in csv.DictReader(open(CORE / "ledger-panel/gene-positions.tsv"), delimiter="\t")
           if r["placed"] == "True"}
    out = {"gene_trees": str(CORE / "gene-trees-panel/loci.treefile"), "reference_chromosome": "ludwigii type",
           "tests": {}}
    for tn, third in THIRDS.items():
        for on, o in FOURTHS.items():
            g = gene_quartets([LUD, CLO, third, o], loci, trees)
            states = {"lud+clo": f"{LUD},{CLO}|{third},{o}",
                      "clo-outside": f"{LUD},{third}|{CLO},{o}",
                      "lud-outside": f"{LUD},{o}|{CLO},{third}"}
            for a, b in (("lud+clo", "clo-outside"), ("clo-outside", "lud-outside"), ("lud+clo", "lud-outside")):
                sa, sb = states[a], states[b]
                seq = [s for _, s in sorted((pos[locus], s) for locus, s in g["per_gene"]
                                            if locus in pos and s in (sa, sb))]
                rt = runs_test(seq, random.Random(1), reps=10000)
                key = f"{tn}|{on} {a} vs {b}"
                out["tests"][key] = {**rt, "n_a": seq.count(sa), "n_b": seq.count(sb)}
                print(key, seq.count(sa), seq.count(sb), rt["runs"], rt["null_median_runs"], rt["p_clustered"],
                      flush=True)
    (CORE / "positional-gene-quartets.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
