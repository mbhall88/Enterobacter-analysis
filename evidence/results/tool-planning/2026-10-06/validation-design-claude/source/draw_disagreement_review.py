"""Draw Michael's Layer 3 disagreement review (validation-design.md section 6), fixed before any result.

Eligible: effective-reserve Layer 3 rows (truth_kind crosswalk) outside the expert sample whose species or
subspecies outcome is discordant. If at most CAP are eligible, all are reviewed. Otherwise, with
order(x) = sha256(SEED|x):
  1. name forms in order(form) each get one slot while slots remain;
  2. the remaining slots go by largest remainder in proportion to each form's discordant count minus
     the slot already given (ties broken by order(form)), never exceeding that count;
  3. within a form, genomes are taken in order(accession).
Inputs are the frozen scorer's scored-genomes.tsv and the Layer 3 truth table, which carries the
crosswalk's ncbi_name_normalised column.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

SEED = "Enterobacter-validation-disagreement-20261007-v1"
CAP = 20


def order(value):
    return hashlib.sha256(f"{SEED}|{value}".encode()).hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def eligible(scored, truth, expert):
    names = {t["genome_id"]: t["ncbi_name_normalised"] for t in truth}
    out = defaultdict(list)
    for r in scored:
        if r["truth_kind"] != "crosswalk" or r["genome_id"] in expert:
            continue
        if r["species_outcome"] == "discordant" or r["subspecies_outcome"] == "discordant":
            out[names[r["genome_id"]]].append(r["genome_id"])
    return out


def allocate(counts, cap=CAP):
    """Slots per name form; counts maps form -> number of eligible discordant genomes."""
    total = sum(counts.values())
    if total <= cap:
        return dict(counts)
    forms = sorted((f for f in counts if counts[f] > 0), key=order)
    slots = {f: 0 for f in counts}
    for f in forms[:cap]:
        slots[f] = 1
    left = cap - sum(slots.values())
    if left:
        rest = {f: counts[f] - slots[f] for f in forms}
        base = sum(rest.values())
        quota = {f: left * rest[f] / base for f in forms}
        for f in forms:
            slots[f] += int(quota[f])
        extra = left - sum(int(q) for q in quota.values())
        for f in sorted(forms, key=lambda f: (-(quota[f] - int(quota[f])), order(f)))[:extra]:
            slots[f] += 1
    assert sum(slots.values()) == cap and all(slots[f] <= counts[f] for f in counts)
    return slots


def draw(scored, truth, expert, cap=CAP):
    pool = eligible(scored, truth, expert)
    slots = allocate({f: len(g) for f, g in pool.items()}, cap)
    chosen = []
    for form in sorted(pool, key=order):
        for rank, genome in enumerate(sorted(pool[form], key=order)[:slots[form]], start=1):
            chosen.append(dict(genome_id=genome, ncbi_name_normalised=form, rank_within_form=rank))
    return chosen, dict(eligible=sum(len(g) for g in pool.values()), forms=len(pool),
                        slots={f: s for f, s in sorted(slots.items()) if s}, chosen=len(chosen), cap=cap, seed=SEED)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scored", type=Path, required=True, help="scorer scored-genomes.tsv for the effective Layer 3 table")
    parser.add_argument("--truth", type=Path, required=True, help="Layer 3 truth table with ncbi_name_normalised")
    parser.add_argument("--expert-sample", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("refusing to overwrite output")
    expert = {r["assembly_unit"] for r in read_tsv(args.expert_sample)}
    chosen, record = draw(read_tsv(args.scored), read_tsv(args.truth), expert)
    args.output.mkdir(parents=True)
    with (args.output / "disagreement-review.tsv").open("x", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=["genome_id", "ncbi_name_normalised", "rank_within_form"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(chosen)
    record["inputs"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.scored, args.truth, args.expert_sample)}
    record["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output / "draw.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
