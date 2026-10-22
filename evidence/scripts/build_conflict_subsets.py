#!/usr/bin/env python3
"""Freeze the Task 17 focused subsets, and verify tip sets against them.

Task 17 (SET-22) tests three named conflicts between the BAC120 marker tree and the two
whole-genome distance trees.  Every later line of evidence runs on a small, frozen subset, so the
one thing this script must not do is let a subset be chosen to favour a side.  Three rules follow.

**Competing clades are derived from the trees, never listed.**  For each focal type set F and each
method, the competing clade is the smallest clade of that method's rooted tree that holds F and at
least one other type anchor -- the clade that decides F's nearest named type under that method.
Every genome in every competing clade enters the subset, taken from each method (the union), so
each side's clade is present in full.

**Everything else that enters carries its reason.**  Neighbours and outgroups come from the
configuration with a written justification; an outgroup is refused unless it lies outside the MRCA
of the subset's ingroup in all three rooted trees.

**Single-gene analyses get a thinned panel chosen without any method's distances.**  The subset is
decomposed into lineage units: the maximal clades recovered with identical membership by all three
trees.  Each unit keeps its type anchors and is topped up by a seeded random draw.  Genomes no
shared clade contains are kept individually -- they are where the methods disagree.

Hypotheses are written as bipartitions of subset tips and as IQ-TREE ``-g`` constraint trees, for
the full subset and for the thinned panel.  Tips named in ``free_in_constraints`` (the second
*E. sichuanensis* assembly) are left out of every constraint so the data place them.

``verify`` checks that a tree, alignment or id list carries exactly a frozen tip set.

Standard library only; the Newick machinery is Task 05's (``scripts/gtdb_context.py``).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import Tree, parse_newick  # noqa: E402

METHODS = ("bac120_fasttree", "verticall_fastme", "mashtree_consensus")


class SubsetError(RuntimeError):
    """A subset cannot be built as configured."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------------------------
# Tree helpers
# --------------------------------------------------------------------------------------------

class CladeIndex:
    """Every internal clade of a rooted tree as a frozenset of tip labels, with its node label."""

    def __init__(self, tree: Tree):
        self.tree = tree
        self.node_of_tip = {tree.label[n]: n for n in tree.tips()}
        if len(self.node_of_tip) != len(tree.tips()):
            raise SubsetError("tree has duplicate tip labels")
        self.below: dict[int, frozenset[str]] = {}
        for node in tree.postorder():
            if not tree.children[node]:
                self.below[node] = frozenset([tree.label[node]])
            else:
                self.below[node] = frozenset().union(*(self.below[c] for c in tree.children[node]))
        self.clades: dict[frozenset[str], str] = {
            self.below[n]: tree.label[n] for n in self.below if tree.children[n]}

    @property
    def tips(self) -> frozenset[str]:
        return self.below[0]

    def mrca_clade(self, tips) -> tuple[int, frozenset[str]]:
        missing = [t for t in tips if t not in self.node_of_tip]
        if missing:
            raise SubsetError(f"tips absent from tree: {missing}")
        node = self.tree.mrca([self.node_of_tip[t] for t in tips])
        return node, self.below[node]

    def lineage_clade(self, named: list[str], types: set[str]) -> frozenset[str]:
        """Largest clade holding ``named`` and no other type anchor."""
        node, clade = self.mrca_clade(named)
        if (clade & types) - set(named):
            raise SubsetError(f"the MRCA of {named} already holds another type anchor")
        while self.tree.parent[node] != -1:
            parent = self.tree.parent[node]
            if (self.below[parent] & types) - set(named):
                break
            node, clade = parent, self.below[parent]
        return clade

    def competing_clade(self, focal: list[str], types: set[str]) -> tuple[frozenset[str], str]:
        """Smallest clade holding ``focal`` and at least one type anchor outside ``focal``."""
        node, clade = self.mrca_clade(focal)
        while True:
            if (clade & types) - set(focal):
                return clade, self.tree.label[node]
            parent = self.tree.parent[node]
            if parent == -1:
                raise SubsetError(f"no clade holds {focal} and another type anchor")
            node, clade = parent, self.below[parent]


def load_trees(config: dict, root: Path) -> dict[str, CladeIndex]:
    trees = {}
    for method in METHODS:
        path = root / config["trees"][method]
        trees[method] = CladeIndex(parse_newick(path.read_text()))
    tip_sets = {m: t.tips for m, t in trees.items()}
    first = tip_sets[METHODS[0]]
    for method, tips in tip_sets.items():
        if tips != first:
            raise SubsetError(f"{method} tips differ from {METHODS[0]}: "
                              f"{sorted(tips ^ first)[:5]}")
    return trees


# --------------------------------------------------------------------------------------------
# Subset construction
# --------------------------------------------------------------------------------------------

def lineage_units(ingroup: set[str], trees: dict[str, CladeIndex], types: set[str],
                  type_groups: list[set[str]]) -> list[frozenset[str]]:
    """Maximal clades shared exactly by all three trees, inside ``ingroup``, naming one lineage.

    A clade qualifies only if its type anchors all belong to one of ``type_groups`` (a focal set or
    a lineage a hypothesis names); every other type anchor counts as its own group.  Without that
    rule a clade all three methods share *around* the conflict -- the 173-tip clade holding the
    whole trio -- would become one unit and be thinned away, which is the structure under test.
    Genomes that no qualifying clade (of two or more tips) contains become singleton units, so the
    units partition ``ingroup``.
    """
    group_of = {}
    for i, group in enumerate(type_groups):
        for tip in group:
            group_of[tip] = i
    shared = set(trees[METHODS[0]].clades)
    for method in METHODS[1:]:
        shared &= set(trees[method].clades)

    def one_lineage(clade: frozenset[str]) -> bool:
        named = clade & types
        return len({group_of.get(t, t) for t in named}) <= 1

    inside = [c for c in shared if c <= ingroup and len(c) >= 2 and one_lineage(c)]
    maximal = [c for c in inside if not any(c < other for other in inside)]
    covered: set[str] = set()
    for clade in maximal:
        if covered & clade:  # pragma: no cover - maximal clades of one tree are disjoint
            raise SubsetError("maximal shared clades overlap")
        covered |= clade
    singles = [frozenset([t]) for t in sorted(ingroup - covered)]
    return sorted(maximal, key=lambda c: (-len(c), min(c))) + singles


def lineage_core(named: list[str], trees: dict[str, CladeIndex],
                 types: set[str]) -> tuple[frozenset[str], frozenset[str], dict[str, int]]:
    """The genomes every tree places in a named lineage, and those only some trees do."""
    per_method = {m: t.lineage_clade(named, types) for m, t in trees.items()}
    clades = list(per_method.values())
    core = frozenset.intersection(*clades)
    disputed = frozenset.union(*clades) - core
    return core, disputed, {m: len(c) for m, c in per_method.items()}


def build_conflict(cid: str, spec: dict, config: dict, trees: dict[str, CladeIndex],
                   manifest: dict[str, dict], types: set[str], root: Path) -> dict:
    tree_tips = trees[METHODS[0]].tips
    extras = config.get("extra_genomes", {})
    reasons: dict[str, list[str]] = {}

    def add(tip: str, reason: str) -> None:
        if tip not in tree_tips and tip not in extras:
            raise SubsetError(f"{cid}: {tip} is neither a tree tip nor a declared extra genome")
        reasons.setdefault(tip, [])
        if reason not in reasons[tip]:
            reasons[tip].append(reason)

    competing = []
    for focal in spec["focal_sets"]:
        for tip in focal["tips"]:
            if tip not in types:
                raise SubsetError(f"{cid}: focal tip {tip} is not a type anchor")
            add(tip, f"focal_type:{focal['label']}")
        for method, index in trees.items():
            clade, support = index.competing_clade(focal["tips"], types)
            others = sorted((clade & types) - set(focal["tips"]))
            competing.append({
                "focal_set": focal["label"], "method": method, "tips": len(clade),
                "node_label": support, "other_type_anchors": others,
                "members": sorted(clade)})
            for tip in clade:
                add(tip, f"competing_clade:{method}:{focal['label']}")

    for group in spec["neighbours"]:
        for tip in group["tips"]:
            add(tip, "neighbour")
    outgroups: list[str] = []
    for group in spec["outgroups"]:
        for tip in group["tips"]:
            add(tip, "outgroup")
            outgroups.append(tip)

    members = set(reasons)
    in_trees = members & tree_tips
    ingroup = in_trees - set(outgroups)
    outgroup_check = {}
    for method, index in trees.items():
        _, clade = index.mrca_clade(sorted(ingroup))
        inside = sorted(set(outgroups) & clade)
        outgroup_check[method] = {"ingroup_mrca_tips": len(clade), "outgroups_inside": inside}
        if inside:
            raise SubsetError(f"{cid}: outgroup(s) {inside} fall inside the ingroup MRCA "
                              f"({len(clade)} tips) in {method}")

    type_groups = [set(f["tips"]) for f in spec["focal_sets"]]
    for hyp in spec["hypotheses"].values():
        for key in ("clade_of_cores", "clade_of_everything_except_cores"):
            type_groups += [set(named) for named in hyp.get(key, [])]
    units = lineage_units(ingroup, trees, types, type_groups)
    unit_id = {}
    for i, unit in enumerate(units, 1):
        for tip in unit:
            unit_id[tip] = f"U{i:03d}"

    rng = random.Random(f"{config['thinning']['seed']}:{cid}")
    cap = int(config["thinning"]["max_per_unit"])
    panel: set[str] = set(outgroups) | (members - tree_tips)
    for group in spec["neighbours"]:
        panel |= set(group["tips"])
    unit_rows = []
    for i, unit in enumerate(units, 1):
        typed = sorted(unit & types)
        rest = sorted(unit - set(typed))
        keep = set(typed)
        if len(unit) <= cap:
            keep = set(unit)
        else:
            need = max(cap - len(keep), 0)
            keep |= set(rng.sample(rest, need))
        panel |= keep
        unit_rows.append({
            "unit_id": f"U{i:03d}", "tips_in_subset": len(unit),
            "type_anchors": ";".join(typed),
            "type_names": ";".join(manifest[t]["formal_names"] for t in typed),
            "panel_tips": ";".join(sorted(keep & panel)),
            "members": ";".join(sorted(unit))})

    free = set(spec.get("free_in_constraints", []))
    cores: dict[str, dict] = {}

    def core_of(named: list[str]) -> tuple[frozenset[str], frozenset[str]]:
        key = ",".join(sorted(named))
        if key not in cores:
            core, disputed, sizes = lineage_core(named, trees, types)
            outside = sorted(core - members)
            if outside:
                raise SubsetError(f"{cid}: the core of lineage {key} has {len(outside)} genomes "
                                  f"outside the subset, e.g. {outside[:3]}")
            cores[key] = {"types": sorted(named), "core_tips": len(core),
                          "disputed_tips_in_subset": sorted(disputed & members),
                          "disputed_tips_outside_subset": len(disputed - members),
                          "lineage_clade_tips_per_method": sizes, "core": sorted(core)}
        entry = cores[key]
        return frozenset(entry["core"]), frozenset(entry["disputed_tips_in_subset"])

    hypotheses = {}
    for name, hyp in spec["hypotheses"].items():
        loose: set[str] = set()
        if "clade_of_cores" in hyp:
            clade: set[str] = set()
            for named in hyp["clade_of_cores"]:
                core, disputed = core_of(named)
                clade |= core
                loose |= disputed
            kind = "clade"
        elif "clade_of_everything_except_cores" in hyp:
            excluded: set[str] = set()
            for named in hyp["clade_of_everything_except_cores"]:
                core, disputed = core_of(named)
                excluded |= core
                loose |= disputed
            clade = ingroup - excluded - loose
            kind = "clade"
        elif "quartet_of_types" in hyp:
            pair, _single = hyp["quartet_of_types"]
            clade = set(pair)
            kind = "quartet"
        else:
            raise SubsetError(f"{cid}: hypothesis {name} has no recognised form")
        clade -= free | loose
        hypotheses[name] = {"kind": kind, "source": hyp["source"], "clade": sorted(clade),
                            "free": sorted(free | loose),
                            "spec": {k: v for k, v in hyp.items() if k != "source"}}

    rows = []
    for tip in sorted(members):
        extra = extras.get(tip)
        if extra:
            names, deriv, cohort = extra["formal_names"], root / extra["derivative_fasta_path"], False
        else:
            row = manifest[tip]
            names, deriv, cohort = row["formal_names"], root / row["derivative_fasta_path"], True
        if not deriv.is_file():
            raise SubsetError(f"{cid}: derivative for {tip} missing at {deriv}")
        rows.append({
            "tip_id": tip, "in_cohort": str(cohort), "in_panel": str(tip in panel),
            "unit_id": unit_id.get(tip, "outgroup" if tip in outgroups else "extra"),
            "is_type_anchor": str(tip in types), "formal_names": names,
            "gtdb_species": "" if extra else manifest[tip].get("gtdb_species", ""),
            "reasons": ";".join(reasons[tip]),
            "derivative_fasta_path": str(deriv.relative_to(root)),
            "derivative_sha256": sha256_file(deriv)})
    return {"rows": rows, "competing": competing, "units": unit_rows, "hypotheses": hypotheses,
            "outgroups": outgroups, "outgroup_check": outgroup_check, "panel": sorted(panel),
            "members": sorted(members), "free": sorted(free), "cores": cores}


def constraint_newick(hyp: dict, tips: list[str], outgroups: list[str]) -> str:
    """IQ-TREE ``-g`` constraint for one hypothesis on one tip set.

    Tips in ``hyp["free"]`` are omitted, which leaves IQ-TREE free to place them anywhere.
    """
    present = set(tips)
    free = set(hyp["free"])
    clade = [t for t in hyp["clade"] if t in present]
    if hyp["kind"] == "quartet":
        pair = clade
        single = [t for t in hyp["spec"]["quartet_of_types"][1] if t in present]
        outs = [t for t in outgroups if t in present]
        if len(pair) != 2 or len(single) != 1 or not outs:
            raise SubsetError("quartet constraint needs its two types, the third and an outgroup")
        return f"(({pair[0]},{pair[1]}),{single[0]},({','.join(outs)}));"
    rest = sorted(present - set(clade) - free)
    if len(clade) < 2 or len(rest) < 2:
        raise SubsetError("clade constraint needs two tips inside and two outside")
    inner = clade[0] if len(clade) == 1 else "(" + ",".join(clade) + ")"
    return f"({inner},{','.join(rest)});"


def bipartition(hyp: dict, tips: list[str]) -> dict:
    present = set(tips) - set(hyp["free"])
    inside = sorted(t for t in hyp["clade"] if t in present)
    if hyp["kind"] == "quartet":
        single = hyp["spec"]["quartet_of_types"][1]
        return {"side_a": inside, "side_b_named": [t for t in single if t in present],
                "kind": "quartet"}
    return {"side_a": inside, "side_b": sorted(present - set(inside)), "kind": "clade"}


def run_build(args) -> int:
    root = args.root.resolve()
    config_path = args.config
    config = json.loads(config_path.read_text())
    trees = load_trees(config, root)
    manifest_rows = read_tsv(root / config["manifest"])
    manifest = {r["tip_id"]: r for r in manifest_rows}
    cohort_ids = set((root / config["ids"]).read_text().split())
    if cohort_ids != set(manifest):
        raise SubsetError("manifest tip ids differ from the ids file")
    if not trees[METHODS[0]].tips <= cohort_ids:
        raise SubsetError("a tree tip is not in the cohort manifest")
    types = {t for t, r in manifest.items() if "type_anchor" in r["panel_roles"]}

    selected = args.conflict or list(config["conflicts"])
    summary = {}
    for cid in selected:
        spec = config["conflicts"][cid]
        out = args.output / spec["directory"] / "subset"
        if out.exists() and any(out.iterdir()):
            raise SubsetError(f"{out} is not empty; frozen subsets are never overwritten")
        out.mkdir(parents=True, exist_ok=True)
        result = build_conflict(cid, spec, config, trees, manifest, types, root)
        fields = list(result["rows"][0])
        write_tsv(out / "subset.tsv", fields, result["rows"])
        (out / "subset.ids.txt").write_text("\n".join(result["members"]) + "\n")
        (out / "panel.ids.txt").write_text("\n".join(result["panel"]) + "\n")
        write_tsv(out / "units.tsv", list(result["units"][0]), result["units"])
        (out / "competing-clades.json").write_text(json.dumps(result["competing"], indent=1) + "\n")
        hyp_out = {}
        for scope, tips in (("subset", result["members"]), ("panel", result["panel"])):
            cdir = out / "constraints" / scope
            cdir.mkdir(parents=True, exist_ok=True)
            for name, hyp in result["hypotheses"].items():
                (cdir / f"{name}.nwk").write_text(
                    constraint_newick(hyp, tips, result["outgroups"]) + "\n")
                hyp_out.setdefault(name, {"source": hyp["source"], "kind": hyp["kind"],
                                          "free_in_constraint": hyp["free"]})
                hyp_out[name][scope] = bipartition(hyp, tips)
        (out / "hypotheses.json").write_text(json.dumps(hyp_out, indent=1) + "\n")
        (out / "lineage-cores.json").write_text(json.dumps(result["cores"], indent=1) + "\n")
        record = {
            "conflict": cid, "title": spec["title"],
            "config": str(config_path), "config_sha256": sha256_file(config_path),
            "manifest": config["manifest"], "manifest_sha256": sha256_file(root / config["manifest"]),
            "trees": {m: {"path": config["trees"][m],
                          "sha256": sha256_file(root / config["trees"][m])} for m in METHODS},
            "subset_tips": len(result["members"]), "panel_tips": len(result["panel"]),
            "lineage_units": len(result["units"]),
            "shared_units_of_two_or_more": sum(1 for u in result["units"] if u["tips_in_subset"] > 1),
            "outgroups": result["outgroups"], "outgroup_check": result["outgroup_check"],
            "free_in_constraints": result["free"],
            "competing_clade_sizes": {f"{c['method']}:{c['focal_set']}": c["tips"]
                                      for c in result["competing"]},
            "thinning": config["thinning"],
        }
        (out / "build-record.json").write_text(json.dumps(record, indent=1) + "\n")
        summary[cid] = {k: record[k] for k in ("subset_tips", "panel_tips", "lineage_units",
                                                "competing_clade_sizes")}
    print(json.dumps(summary, indent=1))
    return 0


# --------------------------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------------------------

def tips_of(path: Path) -> list[str]:
    text = path.read_text()
    if text.lstrip().startswith(">"):
        return [line[1:].split()[0] for line in text.splitlines() if line.startswith(">")]
    if text.lstrip().startswith("("):
        tree = parse_newick(text.strip().splitlines()[0])
        return [tree.label[n] for n in tree.tips()]
    return [line.strip().split("\t")[0] for line in text.splitlines() if line.strip()]


def run_verify(args) -> int:
    expected = set(args.ids.read_text().split())
    failures = []
    for path in args.files:
        found = tips_of(path)
        dup = len(found) - len(set(found))
        missing, unexpected = expected - set(found), set(found) - expected
        ok = not (dup or missing or unexpected)
        print(f"{'OK ' if ok else 'BAD'} {path}: {len(found)} tips, {dup} duplicated, "
              f"{len(missing)} missing, {len(unexpected)} unexpected")
        if not ok:
            failures.append(str(path))
    return 1 if failures else 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build", help="freeze the focused subsets")
    build.add_argument("--config", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True, help="the Task 17 run directory")
    build.add_argument("--root", type=Path, default=Path("."))
    build.add_argument("--conflict", action="append", help="build only this conflict (repeatable)")
    verify = sub.add_parser("verify", help="check files carry exactly a frozen tip set")
    verify.add_argument("--ids", type=Path, required=True)
    verify.add_argument("files", type=Path, nargs="+")
    args = parser.parse_args(argv)
    try:
        return run_build(args) if args.command == "build" else run_verify(args)
    except SubsetError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
