#!/usr/bin/env python3
"""Compare the three pilot trees on exactly matched tip sets, without merging their support.

Task 10's whole risk is that a comparison manufactures a result.  Three things manufacture one
here, and this tool is built around refusing all three.

**One.  The three trees carry three different things.**  FastTree gives Shimodaira-Hasegawa-like
*local* supports; IQ-TREE gives ultrafast bootstrap and SH-aLRT, which are separate statistics
again; FastME on Verticall distances and Mashtree carry no support at all and none was generated.
So there is no such thing as "the supported clades" across the three.  Each tree declares its own
support convention in the configuration, each is thresholded on its own terms, and nothing in this
file ever adds, averages or compares two support values from different conventions.  The neutral
word used throughout is **retained**: a split that survives its own tree's stated rule.  What
retention *means* differs per tree and is printed beside every column that uses it.

**Two.  A raw Robinson-Foulds distance over near-zero branches counts arbitrary resolutions of
genuine polytomies as conflict.**  Task 07 measured 31% of the ML tree's internal branches as
near-zero-length.  So the primary comparison collapses weak and near-zero branches first, and then
classifies every retained split of A into one of three states against B -- *present*, *unresolved
in B* (absent but contradicted by nothing B retains) and *conflicting* (incompatible with a split
B also retains).  Only the third is disagreement.  A polytomy resolved one way by A and another way
by B lands in "unresolved", which is what it is.  Raw RF is still reported, with its convention
stated, but never on its own.

**Three.  Agreement can be read as replication when it is not.**  All three trees derive from the
same Task 06 `verticall repair` derivatives, so no agreement between them is independent of that
step; and two of the three are distance methods on whole genomes, so their agreement is weaker
evidence than either agreeing with the marker tree.  This tool therefore reports every pair
separately and labels each tree's inference class; it never produces a single "concordance" number
across all three, and it never ranks the methods.

Everything is on exactly matched tip sets, verified against a manifest before anything is
compared.  Patristic comparisons are Spearman rank correlations only: a protein marker distance, a
Verticall vertical-window distance and a Mash distance are three different units and must never be
put on one scale.

The Newick parser, formatter, tree container and distance-preserving pruner are Task 05's
(``scripts/gtdb_context.py``).  There are no third-party imports: ``dendropy``, ``ete3``, ``scipy``
and ``matplotlib`` are all absent from this project's Python, the figures are emitted as SVG text
directly, and the statistics are written out longhand so the numbers survive without a plotting
stack.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gtdb_context import Tree, format_newick, parse_newick, prune_to_tips  # noqa: E402


class ComparisonError(RuntimeError):
    """A configuration or input problem that must stop the run rather than be worked around."""


# --------------------------------------------------------------------------------------
# Support conventions
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SupportConvention:
    """How one tree's internal labels are to be read -- and that they are not another's.

    ``kind`` is one of:

    ``none``
        The method produces no branch support.  Nothing is parsed, nothing is thresholded on
        support, and a tree with this convention can never contribute a "supported" claim.
    ``sh_like_local``
        FastTree's SH-like local support, a bare number on 0-1.  A *local* resampling statistic:
        it is known to be the more liberal of the statistics here and is not a bootstrap.
    ``ufboot``
        A bare ultrafast-bootstrap percentage on 0-100.
    ``alrt_ufboot``
        IQ-TREE's ``SH-aLRT/UFboot`` pair.  Both components are parsed and kept; ``retain_on``
        selects which is thresholded, and the other is still reported.
    ``seed_stability``
        The fraction of independent sketch seeds whose tree contains the split, on 0-1, as
        produced by ``scripts/summarize_seed_stability.py``.  **This is not branch support.**  It
        measures reproducibility under the method's own randomness, not evidence for the split,
        and ``is_branch_support`` is False so that no report can read it as support.
    """

    #: the kinds that are branch support; everything else is reproducibility or nothing at all
    BRANCH_SUPPORT_KINDS = ("sh_like_local", "ufboot", "alrt_ufboot")

    kind: str
    statistic: str
    retain_at: float | None = None
    scale_max: float | None = None
    retain_on: str = "primary"
    note: str = ""

    KINDS = ("none", "sh_like_local", "ufboot", "alrt_ufboot", "seed_stability")

    def __post_init__(self) -> None:
        if self.kind not in self.KINDS:
            raise ComparisonError(f"unknown support kind {self.kind!r}; expected one of {self.KINDS}")
        if self.kind == "none":
            if self.retain_at is not None:
                raise ComparisonError(
                    "a support convention of kind 'none' cannot carry a retention threshold: "
                    "there is no support to threshold")
        elif self.retain_at is None:
            raise ComparisonError(f"support kind {self.kind!r} needs a retain_at threshold")
        if self.kind == "alrt_ufboot" and self.retain_on not in ("alrt", "ufboot"):
            raise ComparisonError("alrt_ufboot needs retain_on set to 'alrt' or 'ufboot'")

    @property
    def carries_support(self) -> bool:
        """Whether the convention carries a per-branch value that can threshold a split.

        True for ``seed_stability`` as well: it is a real per-branch number and it can filter a
        tree.  Whether that number is *support* is a separate question, answered by
        :attr:`is_branch_support`, and only that one licenses the word "supported".
        """
        return self.kind != "none"

    @property
    def is_branch_support(self) -> bool:
        return self.kind in self.BRANCH_SUPPORT_KINDS

    def parse(self, label: str) -> dict[str, float | None]:
        """Read one internal Newick label into named components.

        Returns ``{"primary": value_or_None, ...}``; ``primary`` is the component that
        ``retain_at`` applies to.  An unparseable or empty label gives ``None`` rather than a
        default, so "no value" is never silently read as zero.
        """
        text = (label or "").strip()
        if self.kind == "none":
            return {"primary": None}
        if not text:
            return {"primary": None}
        if self.kind == "alrt_ufboot":
            if "/" not in text:
                return {"primary": None, "alrt": None, "ufboot": None, "unparsed": text}
            left, _, right = text.partition("/")
            try:
                alrt: float | None = float(left)
            except ValueError:
                alrt = None
            try:
                ufboot: float | None = float(right)
            except ValueError:
                ufboot = None
            primary = alrt if self.retain_on == "alrt" else ufboot
            return {"primary": primary, "alrt": alrt, "ufboot": ufboot}
        try:
            value = float(text)
        except ValueError:
            return {"primary": None, "unparsed": text}
        return {"primary": value}

    def is_retained(self, parsed: dict[str, float | None]) -> bool:
        """A split is support-retained when its own statistic reaches its own threshold.

        With ``kind='none'`` this is vacuously true: there is no support evidence either way, so
        support cannot remove a split.  That asymmetry is real and is reported rather than hidden
        -- a tree with no support keeps splits that a supported tree would drop.
        """
        if not self.carries_support:
            return True
        value = parsed.get("primary")
        if value is None:
            return False
        return value >= float(self.retain_at)

    def describe(self) -> dict:
        return {
            "kind": self.kind,
            "statistic": self.statistic,
            "retention_threshold": self.retain_at,
            "scale_max": self.scale_max,
            "thresholded_component": self.retain_on if self.kind == "alrt_ufboot" else "the value",
            "note": self.note,
            "is_branch_support": self.is_branch_support,
            "comparable_with_another_convention": False,
        }

    @classmethod
    def from_config(cls, spec: dict) -> "SupportConvention":
        known = {"kind", "statistic", "retain_at", "scale_max", "retain_on", "note"}
        unknown = set(spec) - known
        if unknown:
            raise ComparisonError(f"unknown support keys {sorted(unknown)}")
        return cls(
            kind=spec.get("kind", "none"),
            statistic=spec.get("statistic", ""),
            retain_at=spec.get("retain_at"),
            scale_max=spec.get("scale_max"),
            retain_on=spec.get("retain_on", "primary" if spec.get("kind") != "alrt_ufboot" else "ufboot"),
            note=spec.get("note", ""),
        )


# --------------------------------------------------------------------------------------
# Trees, tips and splits
# --------------------------------------------------------------------------------------


@dataclass
class LoadedTree:
    """One named comparator tree plus everything needed to read it on its own terms."""

    name: str
    path: Path
    tree: Tree
    tips: list[str]
    method: str
    inference_class: str
    distance_units: str
    support: SupportConvention
    near_zero_branch_length: float
    near_zero_basis: str
    role: str = "primary"
    sha256: str = ""
    provenance: str = ""


def tip_names(tree: Tree) -> list[str]:
    """Tip labels in traversal order, refusing an unnamed or repeated tip.

    A duplicated tip label is not a cosmetic problem: every split is a set of tip names, so two
    tips sharing a name silently merge and the tree quietly loses a taxon.
    """
    names = [tree.label[i] for i in tree.tips()]
    blank = sum(1 for n in names if not n.strip())
    if blank:
        raise ComparisonError(f"{blank} tip(s) have no label")
    seen: dict[str, int] = {}
    for n in names:
        seen[n] = seen.get(n, 0) + 1
    duplicated = sorted(k for k, v in seen.items() if v > 1)
    if duplicated:
        raise ComparisonError(f"duplicate tip labels: {duplicated[:10]}"
                              f"{' ...' if len(duplicated) > 10 else ''}")
    return names


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical(mask: int, full: int) -> int:
    """One of the two complementary sides of a split, chosen deterministically.

    The side *not* containing tip 0 is used.  That is a genuine canonical form for any tip set --
    unlike "the smaller side", which needs a tie-break at even tip counts.
    """
    return mask if not (mask & 1) else (full ^ mask)


def compatible(x: int, y: int, full: int) -> bool:
    """Two bipartitions are compatible when one of their four quadrant intersections is empty.

    This, not set equality, is what distinguishes a genuine conflict from an unresolved branch.
    A split present in one tree and merely absent from the other is compatible with it: the second
    tree left that question open.  Only an incompatible pair asserts contradictory things.
    """
    xc, yc = full ^ x, full ^ y
    return not (x & y) or not (x & yc) or not (xc & y) or not (xc & yc)


@dataclass
class Split:
    """One non-trivial unrooted bipartition of a tree, with what that tree says about it."""

    mask: int
    size: int               # tips on the canonical side
    smaller_side: int       # min(size, n - size); the informative description of the branch
    length: float | None
    label: str
    parsed: dict[str, float | None]
    support_retained: bool
    length_retained: bool
    merged_branches: int = 1

    @property
    def retained(self) -> bool:
        return self.support_retained and self.length_retained

    @property
    def drop_reasons(self) -> list[str]:
        reasons = []
        if not self.support_retained:
            reasons.append("support_below_threshold")
        if not self.length_retained:
            reasons.append("branch_length_at_or_below_near_zero")
        return reasons


def extract_splits(tree: Tree, index: dict[str, int], convention: SupportConvention,
                   near_zero: float) -> dict[int, Split]:
    """Every non-trivial unrooted split of ``tree``, keyed by its canonical bitmask.

    The input Newick strings are written with a trifurcating root, so the root node is an ordinary
    internal node of the unrooted tree and contributes no edge of its own.  A degree-2 root is the
    one case that needs care -- its two edges are a single unrooted edge, and they are merged here
    rather than double-counted or half-dropped.  Pruning can produce one.
    """
    n = len(index)
    full = (1 << n) - 1
    below: dict[int, int] = {}
    count: dict[int, int] = {}
    for node in tree.postorder():
        if not tree.children[node]:
            name = tree.label[node]
            bit = index.get(name)
            below[node] = 0 if bit is None else (1 << bit)
            count[node] = 0 if bit is None else 1
            continue
        mask = 0
        total = 0
        for child in tree.children[node]:
            mask |= below[child]
            total += count[child]
        below[node] = mask
        count[node] = total

    root_children = tree.children[0]
    degree_two_root = len(root_children) == 2

    splits: dict[int, Split] = {}

    def add(mask: int, length: float | None, label: str, branches: int = 1) -> None:
        size = bin(mask).count("1")
        if size < 2 or size > n - 2:
            return  # trivial: a tip, or everything but a tip
        key = canonical(mask, full)
        parsed = convention.parse(label)
        support_ok = convention.is_retained(parsed)
        length_ok = length is None or length > near_zero
        existing = splits.get(key)
        if existing is None:
            splits[key] = Split(
                mask=key,
                size=bin(key).count("1"),
                smaller_side=min(size, n - size),
                length=length,
                label=label,
                parsed=parsed,
                support_retained=support_ok,
                length_retained=length_ok,
                merged_branches=branches,
            )
            return
        # The same bipartition reached twice.  That happens at a degree-2 root, and after a
        # projection onto fewer tips.  Merge rather than discard: sum the lengths, and keep the
        # stronger claim the tree makes about the bipartition.
        existing.merged_branches += branches
        if length is not None:
            existing.length = (existing.length or 0.0) + length
        old = existing.parsed.get("primary")
        new = parsed.get("primary")
        if new is not None and (old is None or new > old):
            existing.label = label
            existing.parsed = parsed
            existing.support_retained = support_ok
        existing.length_retained = (existing.length is None
                                    or existing.length > near_zero)

    for node in range(tree.size):
        if node == 0:
            continue
        if degree_two_root and node == root_children[1]:
            # Its split is the complement of its sibling's; the edges were merged above.
            continue
        length = tree.blen[node]
        label = tree.label[node] if tree.children[node] else ""
        branches = 1
        if degree_two_root and node == root_children[0]:
            # The two root edges are one edge of the unrooted tree: sum their lengths, and take
            # whichever end carries a support label rather than dropping the other's.
            sibling = root_children[1]
            length = (length or 0.0) + (tree.blen[sibling] or 0.0)
            sibling_label = tree.label[sibling] if tree.children[sibling] else ""
            if not label.strip():
                label = sibling_label
            branches = 2
        add(below[node], length, label, branches)
    return splits


def retained_masks(splits: dict[int, Split]) -> set[int]:
    return {k for k, s in splits.items() if s.retained}


# --------------------------------------------------------------------------------------
# Pairwise comparison
# --------------------------------------------------------------------------------------


def conflicting_partners(mask: int, others: list[int], full: int) -> list[int]:
    return [o for o in others if not compatible(mask, o, full)]


def classify_direction(a: dict[int, Split], b: dict[int, Split], full: int) -> dict:
    """Place every retained split of A into exactly one state with respect to B.

    The four states are exhaustive and mutually exclusive:

    ``present``
        B has the same bipartition.  Reported split by whether B retains it too, because a split
        B has but does not retain is agreement of topology without agreement of confidence.
    ``conflicting_with_retained``
        Some split B *retains* is incompatible with it.  This is disagreement, and it is the only
        state that should ever be written up as one.
    ``conflicting_with_dropped_only``
        B resolves that region differently, but only on a branch B's own rule drops.  B is not
        asserting the alternative.
    ``unresolved_in_b``
        B neither has the split nor contradicts it.  A polytomy resolved by A and left open by B
        lands here.  Counting this as conflict is the single commonest way to overstate a
        Robinson-Foulds distance, and it is why it is separated out.
    """
    b_all = list(b)
    b_retained = [k for k, s in b.items() if s.retained]
    present_retained: list[int] = []
    present_dropped: list[int] = []
    conflicting_retained: list[tuple[int, list[int]]] = []
    conflicting_dropped: list[tuple[int, list[int]]] = []
    unresolved: list[int] = []
    for key, split in a.items():
        if not split.retained:
            continue
        if key in b:
            (present_retained if b[key].retained else present_dropped).append(key)
            continue
        against_retained = conflicting_partners(key, b_retained, full)
        if against_retained:
            conflicting_retained.append((key, against_retained))
            continue
        against_any = conflicting_partners(key, b_all, full)
        if against_any:
            conflicting_dropped.append((key, against_any))
            continue
        unresolved.append(key)
    return {
        "present_and_retained_in_b": present_retained,
        "present_but_dropped_in_b": present_dropped,
        "conflicting_with_a_retained_split_of_b": conflicting_retained,
        "conflicting_only_with_a_dropped_split_of_b": conflicting_dropped,
        "unresolved_in_b": unresolved,
    }


def rf_stats(a_keys: set[int], b_keys: set[int]) -> dict:
    only_a, only_b = a_keys - b_keys, b_keys - a_keys
    raw = len(only_a) + len(only_b)
    denominator = len(a_keys) + len(b_keys)
    return {
        "splits_in_a": len(a_keys),
        "splits_in_b": len(b_keys),
        "shared": len(a_keys & b_keys),
        "only_in_a": len(only_a),
        "only_in_b": len(only_b),
        "raw_robinson_foulds": raw,
        "normalised_robinson_foulds": round(raw / denominator, 6) if denominator else 0.0,
        "percent_of_a_splits_absent_from_b": (
            round(100 * len(only_a) / len(a_keys), 1) if a_keys else None),
        "percent_of_b_splits_absent_from_a": (
            round(100 * len(only_b) / len(b_keys), 1) if b_keys else None),
        "identical": not only_a and not only_b,
        "convention": (
            "raw_robinson_foulds is the symmetric difference of the two non-trivial unrooted "
            "split sets, counting each missing split once in each direction; "
            "normalised_robinson_foulds divides it by the total number of splits in the two trees "
            "(splits_in_a + splits_in_b), not by the 2(n-3) maximum. Trees resolved to different "
            "degrees have different denominators, which is why this number is never reported "
            "without the split counts beside it."
        ),
    }


@dataclass
class PanelMetadata:
    """Names and roles for the tips, so a contested split can be described nomenclaturally."""

    type_names: dict[str, list[str]] = field(default_factory=dict)      # tip -> formal names
    type_ranks: dict[str, list[str]] = field(default_factory=dict)
    ingroup_status: dict[str, str] = field(default_factory=dict)
    gtdb_species: dict[str, str] = field(default_factory=dict)
    roles: dict[str, str] = field(default_factory=dict)

    def names_on(self, tips: list[str]) -> list[str]:
        out: list[str] = []
        for tip in tips:
            out.extend(self.type_names.get(tip, []))
        return sorted(out)


def describe_split(mask: int, order: list[str], meta: PanelMetadata, full: int,
                   max_tips: int = 8) -> dict:
    """A split described the way a reader needs it: by size, and by the type names it separates."""
    n = len(order)
    side = [order[i] for i in range(n) if mask >> i & 1]
    other = [order[i] for i in range(n) if not (mask >> i & 1)]
    if len(side) > len(other):
        side, other = other, side
    names_in = meta.names_on(side)
    names_out = meta.names_on(other)
    return {
        "smaller_side_tips": len(side),
        "larger_side_tips": len(other),
        "type_names_on_smaller_side": names_in,
        "type_name_count_on_smaller_side": len(names_in),
        "type_names_on_larger_side_count": len(names_out),
        "separates_named_types": bool(names_in) and bool(names_out),
        "example_tips": sorted(side)[:max_tips],
        "smaller_side_tips_truncated": len(side) > max_tips,
    }


def compare_pair(name_a: str, a: dict[int, Split], name_b: str, b: dict[int, Split],
                 order: list[str], meta: PanelMetadata, backbone_min: int,
                 max_reported: int = 40) -> dict:
    n = len(order)
    full = (1 << n) - 1
    a_ret, b_ret = retained_masks(a), retained_masks(b)

    forward = classify_direction(a, b, full)
    reverse = classify_direction(b, a, full)

    def backbone(keys) -> list[int]:
        return [k for k in keys if min(bin(k).count("1"), n - bin(k).count("1")) >= backbone_min]

    conflict_pairs: set[tuple[int, int]] = set()
    for key, partners in forward["conflicting_with_a_retained_split_of_b"]:
        for p in partners:
            conflict_pairs.add((key, p))
    for key, partners in reverse["conflicting_with_a_retained_split_of_b"]:
        for p in partners:
            conflict_pairs.add((p, key))

    contested = []
    for key, partner in sorted(conflict_pairs,
                               key=lambda kp: -min(bin(kp[0]).count("1"),
                                                   n - bin(kp[0]).count("1"))):
        entry = {
            f"{name_a}_split": describe_split(key, order, meta, full),
            f"{name_b}_split": describe_split(partner, order, meta, full),
            f"{name_a}_support_label": a[key].label or None,
            f"{name_b}_support_label": b[partner].label or None,
            f"{name_a}_branch_length": a[key].length,
            f"{name_b}_branch_length": b[partner].length,
        }
        contested.append(entry)
    type_anchored = [c for c in contested
                     if c[f"{name_a}_split"]["separates_named_types"]
                     or c[f"{name_b}_split"]["separates_named_types"]]

    def direction_summary(d: dict, source: dict[int, Split], other_name: str) -> dict:
        retained_total = sum(1 for s in source.values() if s.retained)
        conflicting = [k for k, _ in d["conflicting_with_a_retained_split_of_b"]]
        conflicting_weak = [k for k, _ in d["conflicting_only_with_a_dropped_split_of_b"]]
        present = d["present_and_retained_in_b"] + d["present_but_dropped_in_b"]
        return {
            "retained_splits": retained_total,
            "present_in_other_tree": len(present),
            "present_and_retained_in_other_tree": len(d["present_and_retained_in_b"]),
            "present_but_dropped_in_other_tree": len(d["present_but_dropped_in_b"]),
            "conflicting_with_a_retained_split": len(conflicting),
            "conflicting_only_with_a_dropped_split": len(conflicting_weak),
            "unresolved_in_other_tree": len(d["unresolved_in_b"]),
            "percent_present_in_other_tree": (
                round(100 * len(present) / retained_total, 1) if retained_total else None),
            "percent_conflicting_with_a_retained_split": (
                round(100 * len(conflicting) / retained_total, 1) if retained_total else None),
            "backbone_retained_splits": len(backbone([k for k, s in source.items() if s.retained])),
            "backbone_conflicting_with_a_retained_split": len(backbone(conflicting)),
            "note": (f"'retained' here means retained under {name_a if source is a else name_b}'s "
                     f"own rule; it is not the same evidence as {other_name}'s."),
        }

    return {
        "trees": [name_a, name_b],
        "tips": n,
        "raw_split_comparison": rf_stats(set(a), set(b)),
        "retained_split_comparison": rf_stats(a_ret, b_ret),
        f"{name_a}_against_{name_b}": direction_summary(forward, a, name_b),
        f"{name_b}_against_{name_a}": direction_summary(reverse, b, name_a),
        "mutually_retained_conflicting_split_pairs": len(conflict_pairs),
        "mutually_retained_conflicting_split_pairs_backbone": sum(
            1 for x, y in conflict_pairs
            if min(bin(x).count("1"), n - bin(x).count("1")) >= backbone_min
            and min(bin(y).count("1"), n - bin(y).count("1")) >= backbone_min),
        "type_anchored_contested_bipartitions": len(type_anchored),
        "contested_bipartitions": contested[:max_reported],
        "contested_bipartitions_truncated": max(0, len(contested) - max_reported),
        "type_anchored_contested_detail": type_anchored[:max_reported],
    }


# --------------------------------------------------------------------------------------
# Patristic ranks -- descriptive only
# --------------------------------------------------------------------------------------


def patristic_distances(tree: Tree, order: list[str]) -> list[float]:
    """All pairwise path lengths, in the fixed order (i, j) for i < j.

    Returned as a flat list in row-major upper-triangular order so two trees on the same tip
    order are directly comparable position by position.  The values are **not** comparable in
    magnitude across methods -- a protein marker branch length, a Verticall vertical-window
    distance and a Mash distance are three different units -- so these are only ever fed to a
    rank correlation.
    """
    position = {name: i for i, name in enumerate(order)}
    node_of: dict[int, int] = {}
    for node in tree.tips():
        name = tree.label[node]
        if name in position:
            node_of[position[name]] = node
    if len(node_of) != len(order):
        missing = sorted(set(range(len(order))) - set(node_of))
        raise ComparisonError(f"{len(missing)} requested tips are not in the tree")

    adjacency: dict[int, list[tuple[int, float]]] = {i: [] for i in range(tree.size)}
    for node in range(tree.size):
        parent = tree.parent[node]
        if parent == -1:
            continue
        length = tree.blen[node] or 0.0
        adjacency[parent].append((node, length))
        adjacency[node].append((parent, length))

    n = len(order)
    out: list[float] = []
    rows: list[list[float]] = []
    for i in range(n):
        start = node_of[i]
        dist = {start: 0.0}
        stack = [start]
        while stack:
            current = stack.pop()
            base = dist[current]
            for nxt, length in adjacency[current]:
                if nxt in dist:
                    continue
                dist[nxt] = base + length
                stack.append(nxt)
        rows.append([dist[node_of[j]] for j in range(n)])
    for i in range(n):
        for j in range(i + 1, n):
            out.append(rows[i][j])
    return out


def average_ranks(values: list[float]) -> list[float]:
    """Ranks with ties averaged, written out longhand so no third-party library is needed."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = mean_rank
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float | None:
    """Spearman's rho via Pearson on average ranks; ties handled, no library needed."""
    if len(x) != len(y):
        raise ComparisonError("rank correlation needs equal-length inputs")
    if len(x) < 3:
        return None
    rx, ry = average_ranks(x), average_ranks(y)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def patristic_report(trees: dict[str, LoadedTree], order: list[str],
                     pruned: dict[str, Tree]) -> dict:
    distances = {name: patristic_distances(pruned[name], order) for name in trees}
    pairs = sorted(trees)
    out: dict = {
        "units_warning": (
            "These are rank correlations only. A BAC120 branch length is expected amino-acid "
            "substitutions per site, a Verticall/FastME branch length is a vertical-window "
            "genomic distance and a Mashtree branch length is a Mash distance from 21-mer "
            "sketches. They are three different units on three different scales; no common "
            "scale, regression or ratio between them is reported, and none should be derived "
            "from these numbers."
        ),
        "tip_pairs": len(next(iter(distances.values()))) if distances else 0,
        "per_tree": {
            name: {
                "units": trees[name].distance_units,
                "min": min(values), "median": statistics.median(values),
                "mean": round(statistics.fmean(values), 8), "max": max(values),
                "negative_values": sum(1 for v in values if v < 0),
            }
            for name, values in distances.items()
        },
        "spearman_rank_correlation": {},
    }
    for i, a in enumerate(pairs):
        for b in pairs[i + 1:]:
            rho = spearman(distances[a], distances[b])
            out["spearman_rank_correlation"][f"{a}_vs_{b}"] = (
                None if rho is None else round(rho, 6))
    return out


# --------------------------------------------------------------------------------------
# Balanced-panel sensitivity
# --------------------------------------------------------------------------------------


def balanced_subsets(order: list[str], groups: dict[str, str], max_per_group: int,
                     replicates: int, seed: int) -> list[list[str]]:
    """Cap every group at ``max_per_group`` tips, so agreement cannot be carried by one big clade.

    The panel is deliberately unbalanced -- the *hormaechei* complex alone holds about 180 of the
    388 tips -- so a headline agreement figure over all tips is partly a statement about how many
    near-identical genomes happen to sit in one group.  Each replicate draws a different balanced
    subset from a seeded generator, and the spread across replicates is reported rather than one
    draw.
    """
    buckets: dict[str, list[str]] = {}
    for tip in order:
        buckets.setdefault(groups.get(tip, "__ungrouped__"), []).append(tip)
    rng = random.Random(seed)
    subsets = []
    for _ in range(replicates):
        chosen: list[str] = []
        for key in sorted(buckets):
            members = sorted(buckets[key])
            if len(members) <= max_per_group:
                chosen.extend(members)
            else:
                chosen.extend(rng.sample(members, max_per_group))
        subsets.append(sorted(chosen))
    return subsets


# --------------------------------------------------------------------------------------
# Figures, written as SVG text
# --------------------------------------------------------------------------------------
#
# ``matplotlib`` is absent from this project's Python and so are ``ete3`` and ``dendropy``.
# Rather than pin a plotting environment for two line drawings, the figures are emitted as SVG
# directly.  That keeps the analysis dependency-free -- the numbers survive without any plotting
# stack, which is the stated requirement -- and it makes the figures deterministic text that
# diffs, rather than binaries that have to be regenerated to be checked.

COLOUR = {
    "retained": "#1f2933",
    "support_dropped": "#b45309",
    "near_zero_dropped": "#9aa5b1",
    "conflict": "#b91c1c",
    "connector": "#7b8794",
    "connector_crossed": "#b45309",
    "text": "#1f2933",
    "faint": "#cbd2d9",
}


def _escape(text: str) -> str:
    return (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def _caption(text: str, x: float, y: float, width_px: float, font: float = 10.5,
             fill: str = "#52606d") -> tuple[list[str], float]:
    """Wrap a caption to the drawing's width and return its lines plus the height used.

    SVG does not wrap text, so a long caption silently runs off the canvas or collides with
    whatever is drawn beside it. These captions carry the support semantics, so losing them is
    not cosmetic.
    """
    import textwrap
    per_line = max(20, int(width_px / (font * 0.52)))
    lines = []
    for paragraph in text.split("\n"):
        lines.extend(textwrap.wrap(paragraph, per_line) or [""])
    step = font * 1.32
    parts = [f'<text x="{x:.0f}" y="{y + i * step:.1f}" font-size="{font}" fill="{fill}">'
             f'{_escape(line)}</text>' for i, line in enumerate(lines)]
    return parts, len(lines) * step


def cladogram_layout(tree: Tree, tip_order: list[str] | None = None) -> tuple[dict, list[int]]:
    """Positions for a rectangular cladogram: x by depth, y by leaf order.

    ``tip_order`` gives a target vertical order; children are sorted by the mean target position
    of their leaves, which is the standard greedy untangling and is what keeps a tanglegram
    readable.  Branch length is deliberately *not* used for x: the three trees are in three
    different units and drawing them to length would invite exactly the cross-scale reading the
    rest of this tool refuses.
    """
    target = {name: i for i, name in enumerate(tip_order)} if tip_order else None
    leaves_below: dict[int, list[int]] = {}
    for node in tree.postorder():
        if not tree.children[node]:
            leaves_below[node] = [node]
        else:
            merged: list[int] = []
            for child in tree.children[node]:
                merged.extend(leaves_below[child])
            leaves_below[node] = merged

    def key(child: int) -> float:
        if target is None:
            return 0.0
        values = [target[tree.label[leaf]] for leaf in leaves_below[child]
                  if tree.label[leaf] in target]
        return statistics.fmean(values) if values else 0.0

    ordered_tips: list[int] = []
    depth: dict[int, int] = {0: 0}
    stack = [0]
    visit_order: list[int] = []
    while stack:
        node = stack.pop()
        visit_order.append(node)
        children = list(tree.children[node])
        if target is not None:
            children.sort(key=key)
        for child in reversed(children):
            depth[child] = depth[node] + 1
            stack.append(child)
        if not children:
            ordered_tips.append(node)
    max_depth = max(depth.values()) or 1
    y: dict[int, float] = {node: float(i) for i, node in enumerate(ordered_tips)}
    for node in tree.postorder():
        if tree.children[node]:
            y[node] = statistics.fmean(y[c] for c in tree.children[node])
    position = {node: (depth[node] / max_depth, y[node]) for node in range(tree.size)}
    return position, ordered_tips


def _tree_paths(tree: Tree, position: dict, x0: float, x1: float, y0: float, row: float,
                colour_of, mirrored: bool) -> list[str]:
    parts = []

    def px(node: int) -> float:
        frac = position[node][0]
        return x1 - frac * (x1 - x0) if mirrored else x0 + frac * (x1 - x0)

    def py(node: int) -> float:
        return y0 + position[node][1] * row

    for node in range(tree.size):
        if node == 0:
            continue
        parent = tree.parent[node]
        colour, width = colour_of(node)
        parts.append(
            f'<path d="M {px(parent):.2f} {py(parent):.2f} L {px(parent):.2f} {py(node):.2f} '
            f'L {px(node):.2f} {py(node):.2f}" fill="none" stroke="{colour}" '
            f'stroke-width="{width}" stroke-linecap="round"/>')
    return parts


def render_tanglegram(left: Tree, right: Tree, left_name: str, right_name: str,
                      left_caption: str, right_caption: str, title: str,
                      row: float = 9.0, label_tips: bool = False,
                      highlight: set[str] | None = None) -> str:
    """Two trees face to face with their shared tips joined, crossings drawn as they fall.

    The tip sets are identical by construction, so a tip is labelled once in the middle rather
    than twice.  The right tree's children are ordered by the mean position of their leaves in the
    left tree -- the standard greedy untangling -- so the crossings that remain are the ones the
    topologies actually force, not an artefact of Newick sibling order.
    """
    left_pos, left_tips = cladogram_layout(left)
    order = [left.label[t] for t in left_tips]
    right_pos, right_tips = cladogram_layout(right, tip_order=order)

    n = len(left_tips)
    panel = 250.0
    label_width = 176.0 if label_tips else 0.0
    connector = 200.0
    margin = 30.0
    left_x0 = margin
    left_x1 = left_x0 + panel
    label_x = left_x1 + 6.0
    connector_x0 = label_x + label_width + 6.0
    connector_x1 = connector_x0 + connector
    right_x0 = connector_x1 + 6.0
    right_x1 = right_x0 + panel
    width = right_x1 + margin
    text_width = width - 2 * margin
    caption_parts, used = _caption(f"LEFT \u2014 {left_caption}\n"
                                   f"RIGHT \u2014 {right_caption}",
                                   margin, 44.0, text_width)
    top = 44.0 + used + 14.0
    height = top + n * row + 46

    plain = lambda node: (COLOUR["retained"], 1.0)  # noqa: E731
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica,Arial,sans-serif">',
        f'<rect width="{width:.0f}" height="{height:.0f}" fill="#ffffff"/>',
        f'<text x="{margin:.0f}" y="28" font-size="15" fill="{COLOUR["text"]}">'
        f'{_escape(title)}</text>',
    ] + caption_parts
    parts += _tree_paths(left, left_pos, left_x0, left_x1, top, row, plain, mirrored=False)
    parts += _tree_paths(right, right_pos, right_x0, right_x1, top, row, plain, mirrored=True)

    left_y = {left.label[t]: top + left_pos[t][1] * row for t in left_tips}
    right_y = {right.label[t]: top + right_pos[t][1] * row for t in right_tips}
    left_rank = {name: i for i, name in enumerate(order)}
    right_rank = {right.label[t]: i for i, t in enumerate(right_tips)}
    crossings = 0
    names = sorted(left_y)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if (left_rank[a] - left_rank[b]) * (right_rank[a] - right_rank[b]) < 0:
                crossings += 1
    for name in names:
        moved = left_rank[name] != right_rank[name]
        colour = COLOUR["connector_crossed"] if moved else COLOUR["connector"]
        opacity = 0.75 if moved else 0.28
        if highlight and name in highlight:
            colour, opacity = COLOUR["conflict"], 0.95
        parts.append(
            f'<line x1="{connector_x0:.1f}" y1="{left_y[name]:.2f}" x2="{connector_x1:.1f}" '
            f'y2="{right_y[name]:.2f}" stroke="{colour}" stroke-width="0.7" '
            f'opacity="{opacity}"/>')
        if label_tips:
            parts.append(
                f'<text x="{label_x:.1f}" y="{left_y[name] + 3:.2f}" font-size="9" '
                f'fill="{COLOUR["text"]}">{_escape(name)}</text>')
    footnote, _ = _caption(
        f"{n} shared tips; {crossings} pairwise tip-order inversions between the two leaf orders "
        f"after greedy untangling. Drawn as cladograms: the three methods\u2019 branch lengths "
        f"are in different units and are deliberately not drawn to scale.",
        margin, height - 30, text_width, font=10.0)
    parts.extend(footnote)
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def render_annotated_tree(tree: Tree, splits: dict[int, Split], index: dict[str, int],
                          conflict_masks: set[int], title: str, caption: str,
                          labels: dict[str, str], row: float = 9.0) -> str:
    """One tree with every branch coloured by why it was kept or collapsed, conflicts marked."""
    position, ordered_tips = cladogram_layout(tree)
    n_tips = len(ordered_tips)
    n = len(index)
    full = (1 << n) - 1
    below: dict[int, int] = {}
    for node in tree.postorder():
        if not tree.children[node]:
            bit = index.get(tree.label[node])
            below[node] = 0 if bit is None else 1 << bit
        else:
            mask = 0
            for child in tree.children[node]:
                mask |= below[child]
            below[node] = mask

    width = 940.0
    caption_parts, used = _caption(caption, 34.0, 44.0, width - 68.0)
    top = 44.0 + used + 14.0
    height = top + n_tips * row + 66
    x0, x1 = 34.0, 560.0

    def colour_of(node: int) -> tuple[str, float]:
        if not tree.children[node]:
            return COLOUR["retained"], 0.9
        key = canonical(below[node], full)
        split = splits.get(key)
        if split is None:
            return COLOUR["faint"], 0.8
        if split.retained:
            return COLOUR["retained"], 1.4
        if not split.support_retained:
            return COLOUR["support_dropped"], 1.0
        return COLOUR["near_zero_dropped"], 1.0

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica,Arial,sans-serif">',
        f'<rect width="{width:.0f}" height="{height:.0f}" fill="#ffffff"/>',
        f'<text x="34" y="28" font-size="15" fill="{COLOUR["text"]}">{_escape(title)}</text>',
    ] + caption_parts
    parts += _tree_paths(tree, position, x0, x1, top, row, colour_of, mirrored=False)

    for node in range(tree.size):
        if not tree.children[node] or node == 0:
            continue
        key = canonical(below[node], full)
        if key not in conflict_masks:
            continue
        x = x0 + position[node][0] * (x1 - x0)
        y = top + position[node][1] * row
        parts.append(f'<rect x="{x - 2.6:.2f}" y="{y - 2.6:.2f}" width="5.2" height="5.2" '
                     f'fill="{COLOUR["conflict"]}"/>')

    for node in ordered_tips:
        name = tree.label[node]
        text = labels.get(name)
        if not text:
            continue
        y = top + position[node][1] * row
        x = x0 + position[node][0] * (x1 - x0)
        parts.append(f'<line x1="{x + 3:.2f}" y1="{y:.2f}" x2="{x1 + 16:.2f}" y2="{y:.2f}" '
                     f'stroke="{COLOUR["faint"]}" stroke-width="0.6"/>')
        parts.append(f'<text x="{x1 + 20:.2f}" y="{y + 3:.2f}" font-size="8.4" '
                     f'fill="{COLOUR["text"]}">{_escape(text)}</text>')

    legend = [
        (COLOUR["retained"], "retained under this tree’s own rule"),
        (COLOUR["support_dropped"], "collapsed: support below this tree’s threshold"),
        (COLOUR["near_zero_dropped"], "collapsed: branch at or below near-zero length"),
        (COLOUR["conflict"], "■ incompatible with a split another method also retains"),
    ]
    for i, (colour, text) in enumerate(legend):
        y = height - 52 + i * 12
        parts.append(f'<rect x="34" y="{y - 6:.0f}" width="16" height="3" fill="{colour}"/>')
        parts.append(f'<text x="56" y="{y:.0f}" font-size="9.5" fill="#52606d">'
                     f'{_escape(text)}</text>')
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


# --------------------------------------------------------------------------------------
# Configuration, metadata and the comparison-set driver
# --------------------------------------------------------------------------------------


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_metadata(panel_path: Path, types_path: Path) -> tuple[PanelMetadata, list[dict]]:
    """Names and roles from the frozen manifests.

    ``manifests/types.tsv`` carries 35 nomenclatural type names on 33 genomes: two accessions each
    carry a species name and its nominotypical subspecies name on separate rows.  Both names are
    kept for every such tip.  Collapsing them to one name per genome would silently drop a name
    from the nomenclatural reading of a contested split, which is the whole point of anchoring.
    """
    panel = read_tsv(panel_path)
    types = read_tsv(types_path)
    meta = PanelMetadata()
    for row in panel:
        tip = row["tip_id"]
        meta.ingroup_status[tip] = row.get("ingroup_status", "")
        # A genome absent from the pinned GTDB release has no species to be balanced within, and
        # the 38 such tips here are diverse rather than a group. Each becomes its own group, so a
        # balanced subset keeps them all instead of sampling three out of a false cluster.
        meta.gtdb_species[tip] = row.get("gtdb_species", "") or f"__unassigned__:{tip}"
        meta.roles[tip] = row.get("panel_roles", "")
    anchors = [r for r in types if r.get("reference_role") == "nomenclatural_type_anchor"]
    for row in anchors:
        accession = row.get("selected_accession_version", "").strip()
        if not accession:
            continue
        meta.type_names.setdefault(accession, []).append(row["taxon_name"])
        meta.type_ranks.setdefault(accession, []).append(row.get("rank", ""))
    for key in meta.type_names:
        meta.type_names[key].sort()
    return meta, panel


def load_trees(config: dict, root: Path) -> dict[str, LoadedTree]:
    trees: dict[str, LoadedTree] = {}
    for name, spec in config["trees"].items():
        path = root / spec["path"]
        if not path.exists():
            raise ComparisonError(f"tree {name!r}: {path} does not exist")
        tree = parse_newick(path.read_text())
        tips = tip_names(tree)
        near_zero = float(spec["near_zero_branch_length"])
        if near_zero < 0:
            raise ComparisonError(f"tree {name!r}: near_zero_branch_length must not be negative")
        trees[name] = LoadedTree(
            name=name,
            path=path,
            tree=tree,
            tips=tips,
            method=spec["method"],
            inference_class=spec["inference_class"],
            distance_units=spec["distance_units"],
            support=SupportConvention.from_config(spec.get("support", {"kind": "none"})),
            near_zero_branch_length=near_zero,
            near_zero_basis=spec.get("near_zero_basis", ""),
            role=spec.get("role", "primary_comparator"),
            sha256=sha256_file(path),
            provenance=spec.get("provenance", ""),
        )
    return trees


def _rule_matcher(rule: dict, meta: PanelMetadata, where: str):
    """A ``drop_where``/``keep_where`` rule as (column lookup, predicate).

    ``values`` matches the whole cell; ``contains`` matches one semicolon-separated item within
    it, which is what ``panel_roles`` needs -- a tip can be ``study_genome;type_anchor`` and a
    rule naming ``type_anchor`` has to find it.
    """
    column = rule["column"]
    lookup = {"ingroup_status": meta.ingroup_status, "gtdb_species": meta.gtdb_species,
              "panel_roles": meta.roles}.get(column)
    if lookup is None:
        raise ComparisonError(f"{where} column {column!r} is not one this tool reads")
    if ("values" in rule) == ("contains" in rule):
        raise ComparisonError(f"{where} on {column!r} needs exactly one of 'values', 'contains'")
    if "values" in rule:
        wanted = set(rule["values"])
        return lookup, lambda cell: cell in wanted
    wanted = set(rule["contains"])
    return lookup, lambda cell: bool(wanted & {p.strip() for p in cell.split(";")})


def resolve_tips(spec: dict, base: list[str], meta: PanelMetadata) -> list[str]:
    keep = set(base)
    for tip in spec.get("drop_tips", []):
        if tip not in keep:
            raise ComparisonError(
                f"comparison set {spec['name']!r} drops {tip!r}, which is not in the base tip set")
        keep.discard(tip)
    for rule in spec.get("drop_where", []):
        lookup, matches = _rule_matcher(rule, meta, "drop_where")
        removed = {t for t in list(keep) if matches(lookup.get(t, ""))}
        if not removed:
            raise ComparisonError(
                f"comparison set {spec['name']!r}: the rule on {rule['column']!r} removed no "
                "tip; a rule that matches nothing is a configuration error, not an empty result")
        keep -= removed
    for rule in spec.get("keep_where", []):
        lookup, matches = _rule_matcher(rule, meta, "keep_where")
        keep &= {t for t in keep if matches(lookup.get(t, ""))}
    for tip in spec.get("keep_tips", []):
        if tip in base:
            keep.add(tip)
    if len(keep) < 4:
        raise ComparisonError(f"comparison set {spec['name']!r} keeps only {len(keep)} tips")
    return sorted(keep)


def project_tree(loaded: LoadedTree, keep: list[str]) -> tuple[Tree, dict]:
    """Prune to ``keep``, preserving patristic distances between the tips that remain.

    The pruner sums a suppressed node's branch length into the surviving branch, so distances
    between retained tips are unchanged, and returns the suppressed internal labels rather than
    discarding them.  The surviving branch keeps its own support value -- the tip-ward one of a
    collapsed chain, whose original bipartition is the closer of the two to the projected one --
    and the count of collapsed labels is reported so the substitution is visible.
    """
    by_name = {loaded.tree.label[i]: i for i in loaded.tree.tips()}
    missing = [t for t in keep if t not in by_name]
    if missing:
        raise ComparisonError(
            f"tree {loaded.name!r} is missing {len(missing)} of the requested tips: "
            f"{sorted(missing)[:8]}")
    if len(keep) == len(by_name):
        return loaded.tree, {"pruned": False, "internal_labels_collapsed": 0}
    result = prune_to_tips(loaded.tree, {by_name[t] for t in keep})
    return result.tree, {
        "pruned": True,
        "dropped_tips": len(by_name) - len(keep),
        "internal_labels_collapsed": sum(len(v) for v in result.collapsed_labels.values()),
    }


def split_table_rows(name: str, splits: dict[int, Split], order: list[str],
                     meta: PanelMetadata) -> list[dict]:
    n = len(order)
    full = (1 << n) - 1
    rows = []
    for key, split in sorted(splits.items(), key=lambda kv: -kv[1].smaller_side):
        described = describe_split(key, order, meta, full, max_tips=4)
        rows.append({
            "tree": name,
            "smaller_side_tips": split.smaller_side,
            "branch_length": split.length,
            "support_label": split.label,
            "support_value": split.parsed.get("primary"),
            "support_retained": split.support_retained,
            "length_retained": split.length_retained,
            "retained": split.retained,
            "drop_reasons": ";".join(split.drop_reasons),
            "merged_branches": split.merged_branches,
            "type_names_on_smaller_side": ";".join(described["type_names_on_smaller_side"]),
            "example_tips": ";".join(described["example_tips"]),
        })
    return rows


def _retention_means(loaded: LoadedTree) -> str:
    """One sentence saying exactly what "retained" means for this tree, in its own vocabulary.

    The word "support" appears only where the convention really is branch support; a
    seed-stability filter is named as the reproducibility filter it is.
    """
    if not loaded.support.carries_support:
        return (f"branch length > {loaded.near_zero_branch_length:g} only; this method produces "
                "no branch support, so no support filter can be applied and this tree "
                "necessarily retains splits a support-bearing tree would drop")
    quantity = ("support" if loaded.support.is_branch_support
                else "a per-branch reproducibility statistic that is NOT branch support")
    return (f"{quantity} >= {loaded.support.retain_at} ({loaded.support.statistic}) "
            f"AND branch length > {loaded.near_zero_branch_length:g}")


def tree_split_summary(loaded: LoadedTree, splits: dict[int, Split], n: int,
                       backbone_min: int) -> dict:
    values = list(splits.values())
    backbone = [s for s in values if s.smaller_side >= backbone_min]
    retained = [s for s in values if s.retained]
    near_zero = [s for s in values if not s.length_retained]
    unsupported = [s for s in values if not s.support_retained]
    lengths = [s.length for s in values if s.length is not None]
    support_values = [s.parsed.get("primary") for s in values
                      if s.parsed.get("primary") is not None]
    return {
        "method": loaded.method,
        "inference_class": loaded.inference_class,
        "role": loaded.role,
        "support_convention": loaded.support.describe(),
        "near_zero_branch_length": loaded.near_zero_branch_length,
        "near_zero_basis": loaded.near_zero_basis,
        "tips": n,
        "non_trivial_splits": len(values),
        "maximum_possible_splits": n - 3,
        "percent_of_maximum_resolution": round(100 * len(values) / (n - 3), 1) if n > 3 else None,
        "retained_splits": len(retained),
        "retained_backbone_splits": sum(1 for s in backbone if s.retained),
        "backbone_splits": len(backbone),
        "splits_at_or_below_near_zero_length": len(near_zero),
        "percent_of_splits_at_or_below_near_zero_length": (
            round(100 * len(near_zero) / len(values), 1) if values else None),
        "splits_below_support_threshold": (
            len(unsupported) if loaded.support.carries_support else None),
        "splits_with_no_parseable_support": (
            sum(1 for s in values if s.parsed.get("primary") is None)
            if loaded.support.carries_support else None),
        "branch_length_median": statistics.median(lengths) if lengths else None,
        "branch_length_negative": sum(1 for v in lengths if v < 0),
        "support_median": statistics.median(support_values) if support_values else None,
        "retention_means": _retention_means(loaded),
    }


def collapse_sensitivity(trees: dict[str, LoadedTree], projected: dict[str, Tree],
                         order: list[str], pairs: list[tuple[str, str]],
                         near_zero_multipliers: list[float],
                         support_thresholds: list[float | None]) -> list[dict]:
    """The headline conflict count across a grid of collapsing rules.

    A number that only holds at one threshold is not a result.  This re-runs the mutually-retained
    conflict count over a grid of near-zero branch-length multipliers and support thresholds, so a
    reader can see whether the conclusion is a property of the trees or of the cut-offs.
    """
    index = {name: i for i, name in enumerate(order)}
    n = len(order)
    full = (1 << n) - 1
    rows = []
    for multiplier in near_zero_multipliers:
        for threshold in support_thresholds:
            splits: dict[str, dict[int, Split]] = {}
            for name, loaded in trees.items():
                convention = loaded.support
                if threshold is not None and convention.carries_support:
                    # Overrides are fractions of each convention's own scale, so one grid can
                    # cover a 0-1 local support and a 0-100 bootstrap without ever putting the
                    # two statistics on one scale.
                    convention = SupportConvention(
                        kind=convention.kind, statistic=convention.statistic,
                        retain_at=threshold * (convention.scale_max or 1.0),
                        scale_max=convention.scale_max,
                        retain_on=convention.retain_on, note=convention.note)
                splits[name] = extract_splits(
                    projected[name], index, convention,
                    loaded.near_zero_branch_length * multiplier)
            row = {
                "near_zero_multiplier": multiplier,
                "support_threshold_override_as_fraction_of_each_scale": threshold,
                "retained_splits": {name: len(retained_masks(s)) for name, s in splits.items()},
            }
            for a, b in pairs:
                a_ret = [k for k, s in splits[a].items() if s.retained]
                b_ret = [k for k, s in splits[b].items() if s.retained]
                conflicts = set()
                for key in a_ret:
                    for partner in conflicting_partners(key, b_ret, full):
                        conflicts.add((key, partner))
                row[f"{a}_vs_{b}_mutually_retained_conflicting_pairs"] = len(conflicts)
                row[f"{a}_vs_{b}_shared_retained_splits"] = len(set(a_ret) & set(b_ret))
            rows.append(row)
    return rows


def run_comparison_set(spec: dict, trees: dict[str, LoadedTree], base_tips: list[str],
                       meta: PanelMetadata, config: dict, outdir: Path) -> dict:
    name = spec["name"]
    participating = spec.get("trees") or [n for n, t in trees.items()
                                          if t.role == "primary_comparator"]
    unknown = [t for t in participating if t not in trees]
    if unknown:
        raise ComparisonError(f"comparison set {name!r} names unknown trees {unknown}")
    keep = resolve_tips(spec, base_tips, meta)
    order = keep
    index = {tip: i for i, tip in enumerate(order)}
    n = len(order)
    backbone_min = int(config.get("backbone_min_smaller_side_tips", 10))

    set_dir = outdir / name
    (set_dir / "trees").mkdir(parents=True, exist_ok=True)
    (set_dir / "tips.txt").write_text("\n".join(order) + "\n")

    projected: dict[str, Tree] = {}
    projection_notes: dict[str, dict] = {}
    for tree_name in participating:
        tree, note = project_tree(trees[tree_name], keep)
        projected[tree_name] = tree
        projection_notes[tree_name] = note
        (set_dir / "trees" / f"{tree_name}.nwk").write_text(format_newick(tree) + "\n")
        got = sorted(tip_names(tree))
        if got != order:
            raise ComparisonError(
                f"comparison set {name!r}: {tree_name} has {len(got)} tips after projection, "
                f"not the {n} requested")

    splits = {t: extract_splits(projected[t], index, trees[t].support,
                                trees[t].near_zero_branch_length)
              for t in participating}

    rows: list[dict] = []
    for tree_name in participating:
        rows.extend(split_table_rows(tree_name, splits[tree_name], order, meta))
    write_tsv(set_dir / "splits.tsv", rows)

    per_tree = {t: tree_split_summary(trees[t], splits[t], n, backbone_min)
                | projection_notes[t]
                for t in participating}

    comparisons: dict[str, dict] = {}
    pair_kinds: dict[str, str] = {}
    all_pairs = [(participating[i], participating[j])
                 for i in range(len(participating)) for j in range(i + 1, len(participating))]
    for a, b in all_pairs:
        key = f"{a}_vs_{b}"
        comparisons[key] = compare_pair(a, splits[a], b, splits[b], order, meta, backbone_min)
        if trees[a].role == "primary_comparator" and trees[b].role == "primary_comparator":
            same_class = trees[a].inference_class == trees[b].inference_class
            pair_kinds[key] = ("method_comparison_between_two_"
                               + ("distance_methods_on_whole_genomes"
                                  if same_class and trees[a].inference_class
                                  == "whole_genome_distance"
                                  else "different_inference_classes"))
        else:
            pair_kinds[key] = f"baseline_or_evidence: {trees[a].role} vs {trees[b].role}"

    # A second, exactly matched rule.  Each method's own rule is the primary one, but it is not
    # the same rule: FastTree can drop a branch for weak support and the two distance methods
    # cannot, because they have no support statistic to drop it with.  So the same comparison is
    # also run with every tree collapsed on branch length alone -- one identical rule for all
    # three -- and both are reported.  Neither is a substitute for the other: the own-rule
    # comparison is the one that respects what each method actually claims, and the matched-rule
    # comparison is the one where the collapsing cannot favour a method.
    length_only = SupportConvention(
        kind="none",
        statistic=("topology only; each method's support is deliberately ignored so that one "
                   "identical rule -- collapse branches at or below that tree's own near-zero "
                   "length -- applies to all three"),
    )
    splits_matched = {t: extract_splits(projected[t], index, length_only,
                                        trees[t].near_zero_branch_length)
                      for t in participating}
    comparisons_matched = {f"{a}_vs_{b}": compare_pair(a, splits_matched[a], b,
                                                       splits_matched[b], order, meta,
                                                       backbone_min, max_reported=15)
                           for a, b in all_pairs}

    result: dict = {
        "name": name,
        "description": spec.get("description", ""),
        "required_by_the_brief": bool(spec.get("required", False)),
        "tips": n,
        "tip_set_sha256": hashlib.sha256(("\n".join(order) + "\n").encode()).hexdigest(),
        "tips_dropped_from_base": len(base_tips) - n,
        "dropped_tips": sorted(set(base_tips) - set(order)),
        "trees_compared": participating,
        "per_tree": per_tree,
        "pair_kind": pair_kinds,
        "comparisons": comparisons,
        "matched_rule_comparisons": {
            "rule": ("support ignored for every tree; a split is retained only when its branch "
                     "length exceeds that tree's own near-zero threshold. The thresholds are "
                     "still per tree because the three branch lengths are in three different "
                     "units, but no tree is filtered by a statistic another tree does not have."),
            "retained_splits": {t: len(retained_masks(splits_matched[t]))
                                for t in participating},
            "comparisons": comparisons_matched,
        },
    }

    method_pairs = [(a, b) for a, b in all_pairs
                    if trees[a].role == "primary_comparator"
                    and trees[b].role == "primary_comparator"]
    sweep = config.get("collapse_sensitivity")
    if sweep and method_pairs:
        result["collapse_sensitivity"] = collapse_sensitivity(
            {t: trees[t] for t in participating}, projected, order, method_pairs,
            sweep.get("near_zero_multipliers", [1.0]),
            sweep.get("support_threshold_overrides", [None]))

    if spec.get("patristic", True):
        result["patristic_ranks"] = patristic_report(
            {t: trees[t] for t in participating}, order, projected)

    balanced = config.get("balanced_panel")
    if balanced and spec.get("balanced_panel_sensitivity", False) and method_pairs:
        result["balanced_panel_sensitivity"] = run_balanced_panel(
            balanced, trees, projected, order, meta, method_pairs, backbone_min)

    rooting_spec = spec.get("root")
    if rooting_spec:
        result["rooting"] = run_rooting(rooting_spec, config, projected, participating, set_dir,
                                        meta, config.get("named_type_groups"))
        result["type_anchor_concordance"] = type_anchor_concordance(result["rooting"])

    if spec.get("figures", False):
        result["figures"] = write_figures(spec, trees, projected, splits, index, order, meta,
                                          comparisons, set_dir, config)
    return result


def run_balanced_panel(balanced: dict, trees: dict[str, LoadedTree], projected: dict[str, Tree],
                       order: list[str], meta: PanelMetadata,
                       method_pairs: list[tuple[str, str]], backbone_min: int) -> dict:
    groups = {"gtdb_species": meta.gtdb_species}[balanced.get("group_column", "gtdb_species")]
    subsets = balanced_subsets(order, groups, int(balanced["max_per_group"]),
                               int(balanced["replicates"]), int(balanced["seed"]))
    per_pair: dict[str, list[dict]] = {f"{a}_vs_{b}": [] for a, b in method_pairs}
    sizes = []
    for subset in subsets:
        sizes.append(len(subset))
        sub_index = {tip: i for i, tip in enumerate(subset)}
        sub_splits = {}
        for tree_name in {t for pair in method_pairs for t in pair}:
            tree, _ = project_tree(
                LoadedTree(name=tree_name, path=Path("."), tree=projected[tree_name],
                           tips=order, method="", inference_class="", distance_units="",
                           support=trees[tree_name].support,
                           near_zero_branch_length=trees[tree_name].near_zero_branch_length,
                           near_zero_basis=""), subset)
            sub_splits[tree_name] = extract_splits(
                tree, sub_index, trees[tree_name].support,
                trees[tree_name].near_zero_branch_length)
        for a, b in method_pairs:
            summary = compare_pair(a, sub_splits[a], b, sub_splits[b], subset, meta, backbone_min,
                                   max_reported=0)
            per_pair[f"{a}_vs_{b}"].append({
                "tips": len(subset),
                "retained_a": summary[f"{a}_against_{b}"]["retained_splits"],
                "retained_b": summary[f"{b}_against_{a}"]["retained_splits"],
                "shared_retained": summary["retained_split_comparison"]["shared"],
                "conflicting_pairs": summary["mutually_retained_conflicting_split_pairs"],
                "percent_of_a_retained_present_in_b":
                    summary[f"{a}_against_{b}"]["percent_present_in_other_tree"],
                "percent_of_b_retained_present_in_a":
                    summary[f"{b}_against_{a}"]["percent_present_in_other_tree"],
            })

    def spread(values: list[float | None]) -> dict | None:
        clean = [v for v in values if v is not None]
        if not clean:
            return None
        return {"n": len(clean), "min": min(clean),
                "median": round(statistics.median(clean), 3),
                "mean": round(statistics.fmean(clean), 3), "max": max(clean)}

    return {
        "purpose": (
            "The panel is deliberately unbalanced, so a whole-panel agreement figure is partly a "
            "statement about how many near-identical genomes sit in one group. Each replicate "
            "caps every GTDB species at the stated number of tips and repeats the comparison."),
        "group_column": balanced.get("group_column", "gtdb_species"),
        "max_per_group": balanced["max_per_group"],
        "replicates": balanced["replicates"],
        "seed": balanced["seed"],
        "subset_size": spread([float(s) for s in sizes]),
        "per_pair": {
            key: {
                "conflicting_pairs": spread([r["conflicting_pairs"] for r in rows]),
                "percent_of_a_retained_present_in_b":
                    spread([r["percent_of_a_retained_present_in_b"] for r in rows]),
                "percent_of_b_retained_present_in_a":
                    spread([r["percent_of_b_retained_present_in_a"] for r in rows]),
                "replicates": rows,
            }
            for key, rows in per_pair.items()
        },
    }


def run_rooting(rooting_spec: dict, config: dict, projected: dict[str, Tree],
                participating: list[str], set_dir: Path, meta: PanelMetadata,
                named_groups: dict[str, list[str]] | None) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from root_tree import root_and_compare  # noqa: E402  (local, tested, Task 07's)

    named = config["outgroups"]
    requested = {k: named[k] for k in rooting_spec["outgroup_sets"]}
    reference_rooting = rooting_spec.get("type_anchor_reading_on",
                                         rooting_spec["outgroup_sets"][0])
    out: dict = {"outgroup_sets": requested, "per_tree": {},
                 "type_anchor_reading_on": reference_rooting}
    rooted_dir = set_dir / "rooted"
    rooted_dir.mkdir(parents=True, exist_ok=True)
    for tree_name in participating:
        newick = format_newick(projected[tree_name]) + "\n"
        report = root_and_compare(newick, requested)
        rooted_trees = report.pop("trees")
        for rooting_name, rooted in rooted_trees.items():
            path = rooted_dir / f"{tree_name}.{rooting_name}.rooted.nwk"
            path.write_text(format_newick(rooted) + "\n")
            report["rootings"][rooting_name]["tree"] = str(path)
        if reference_rooting in rooted_trees:
            report["type_anchor_topology"] = type_anchor_topology(
                rooted_trees[reference_rooting], meta, named_groups)
        out["per_tree"][tree_name] = report
    rooted_ok = {t: all(r["rooted"] for r in out["per_tree"][t]["rootings"].values())
                 for t in participating}
    out["every_rooting_succeeded"] = rooted_ok
    out["all_trees_rootable_on_every_requested_outgroup_set"] = all(rooted_ok.values())
    return out


# --------------------------------------------------------------------------------------
# Type-anchored reading of a rooted tree
# --------------------------------------------------------------------------------------


def type_anchor_topology(rooted: Tree, meta: PanelMetadata,
                         named_groups: dict[str, list[str]] | None = None) -> dict:
    """What a rooted tree says about the nomenclatural type anchors, and nothing about names.

    Two questions, both answerable from topology alone and both directly readable against a
    nomenclatural treatment:

    *Exclusivity.*  For a named set of type anchors, does the smallest clade containing all of
    them contain any *other* named type?  A group whose MRCA clade admits another type's genome
    is not a candidate for a single name at that rank on this tree.

    *Nearest named neighbour.*  For every pair of type anchors, how many other type anchors fall
    inside their MRCA clade?  A pair with none is each other's nearest named neighbour under this
    method.  Comparing that set of pairs across methods is a type-anchored concordance statement
    that does not depend on any threshold or support convention.

    Nothing here selects a name or a rank.  ``manifests/types.tsv`` puts 35 type names on 33
    genomes -- two accessions each carry a species name and its nominotypical subspecies name --
    and both names are carried through, never collapsed to one per genome.
    """
    name_of_node: dict[str, int] = {}
    for node in rooted.tips():
        name_of_node[rooted.label[node]] = node
    anchors = sorted(t for t in meta.type_names if t in name_of_node)
    bit = {tip: i for i, tip in enumerate(anchors)}

    anchors_below: dict[int, int] = {}
    tips_below: dict[int, int] = {}
    for node in rooted.postorder():
        if not rooted.children[node]:
            anchors_below[node] = 1 << bit[rooted.label[node]] if rooted.label[node] in bit else 0
            tips_below[node] = 1
        else:
            mask, count = 0, 0
            for child in rooted.children[node]:
                mask |= anchors_below[child]
                count += tips_below[child]
            anchors_below[node] = mask
            tips_below[node] = count

    def mrca(tips: list[str]) -> int:
        return rooted.mrca([name_of_node[t] for t in tips])

    def describe(tips: list[str]) -> dict:
        node = mrca(tips)
        inside = anchors_below[node]
        wanted = 0
        for t in tips:
            wanted |= 1 << bit[t]
        extra = inside & ~wanted
        extra_tips = [anchors[i] for i in range(len(anchors)) if extra >> i & 1]
        return {
            "type_anchor_tips": sorted(tips),
            "type_names": sorted(n for t in tips for n in meta.type_names.get(t, [])),
            "mrca_clade_tips": tips_below[node],
            "other_type_anchors_inside_the_mrca_clade": extra_tips,
            "other_type_names_inside_the_mrca_clade": sorted(
                n for t in extra_tips for n in meta.type_names.get(t, [])),
            "exclusive": extra == 0,
        }

    pairs: dict[str, int] = {}
    nearest: list[list[str]] = []
    for i, a in enumerate(anchors):
        for b in anchors[i + 1:]:
            node = mrca([a, b])
            intervening = bin(anchors_below[node]).count("1") - 2
            pairs[f"{a}|{b}"] = intervening
            if intervening == 0:
                nearest.append([a, b])

    out: dict = {
        "type_anchor_tips": len(anchors),
        "type_names": sum(len(meta.type_names[t]) for t in anchors),
        "tips_carrying_two_type_names": sorted(
            t for t in anchors if len(meta.type_names[t]) > 1),
        "nearest_named_neighbour_pairs": sorted(nearest),
        "nearest_named_neighbour_pair_count": len(nearest),
        "intervening_type_anchors_per_pair": pairs,
        "named_groups": {},
    }
    for label, names in (named_groups or {}).items():
        tips = sorted({t for t in anchors
                       if set(meta.type_names.get(t, [])) & set(names)})
        missing = [n for n in names
                   if not any(n in meta.type_names.get(t, []) for t in anchors)]
        if len(tips) < 2:
            out["named_groups"][label] = {
                "resolvable": False,
                "reason": "fewer than two of the named types are tips of this tree",
                "type_names_not_found": missing,
            }
            continue
        entry = describe(tips)
        entry["type_names_requested"] = sorted(names)
        entry["type_names_not_found"] = missing
        entry["resolvable"] = True
        out["named_groups"][label] = entry
    return out


def type_anchor_concordance(rooting: dict) -> dict:
    """Where the methods agree and disagree about the type anchors, method by method.

    This is the part of the comparison that is nomenclaturally interpretable, and it is
    deliberately the part with no threshold in it: whether a named group's MRCA clade admits
    another type, and which type anchors are each other's nearest named neighbours, are
    statements about topology alone.  No name and no rank is proposed; agreement and disagreement
    are reported for DEC-004 to decide from.
    """
    readings = {name: report["type_anchor_topology"]
                for name, report in rooting["per_tree"].items()
                if "type_anchor_topology" in report}
    if len(readings) < 2:
        return {"comparable_trees": sorted(readings), "note": "fewer than two rooted readings"}

    groups = sorted({label for r in readings.values() for label in r["named_groups"]})
    per_group = {}
    for label in groups:
        entries = {name: r["named_groups"].get(label) for name, r in readings.items()}
        resolvable = {n: e for n, e in entries.items() if e and e.get("resolvable")}
        exclusive = sorted(n for n, e in resolvable.items() if e["exclusive"])
        not_exclusive = sorted(n for n, e in resolvable.items() if not e["exclusive"])
        per_group[label] = {
            "type_names_requested": (next(iter(resolvable.values()))["type_names_requested"]
                                     if resolvable else []),
            "methods_recovering_the_group_exclusively": exclusive,
            "methods_admitting_another_type_into_the_mrca_clade": not_exclusive,
            "unanimous": len(exclusive) == len(resolvable) and bool(resolvable),
            "per_method": {
                n: {"mrca_clade_tips": e["mrca_clade_tips"],
                    "exclusive": e["exclusive"],
                    "other_type_names_inside_the_mrca_clade":
                        e["other_type_names_inside_the_mrca_clade"]}
                for n, e in sorted(resolvable.items())
            },
        }

    pair_sets = {name: {tuple(p) for p in r["nearest_named_neighbour_pairs"]}
                 for name, r in readings.items()}
    shared = set.intersection(*pair_sets.values()) if pair_sets else set()
    union = set().union(*pair_sets.values()) if pair_sets else set()
    return {
        "comparable_trees": sorted(readings),
        "named_groups": per_group,
        "named_groups_unanimous": sorted(k for k, v in per_group.items() if v["unanimous"]),
        "named_groups_split": sorted(k for k, v in per_group.items() if not v["unanimous"]),
        "nearest_named_neighbour_pairs": {
            "per_method_count": {n: len(v) for n, v in sorted(pair_sets.items())},
            "found_by_every_method": sorted("|".join(p) for p in sorted(shared)),
            "found_by_every_method_count": len(shared),
            "found_by_at_least_one_method_count": len(union),
            "found_by_only_some": {
                "|".join(p): sorted(n for n, v in pair_sets.items() if p in v)
                for p in sorted(union - shared)
            },
        },
    }


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def conflict_masks_for(name: str, splits: dict[str, dict[int, Split]], others: list[str],
                       full: int) -> set[int]:
    """Splits of ``name`` that are incompatible with a split another named tree also retains."""
    mine = [k for k, s in splits[name].items() if s.retained]
    found: set[int] = set()
    for other in others:
        theirs = [k for k, s in splits[other].items() if s.retained]
        for key in mine:
            if key in found:
                continue
            if conflicting_partners(key, theirs, full):
                found.add(key)
    return found


def write_figures(spec: dict, trees: dict[str, LoadedTree], projected: dict[str, Tree],
                  splits: dict[str, dict[int, Split]], index: dict[str, int], order: list[str],
                  meta: PanelMetadata, comparisons: dict, set_dir: Path, config: dict) -> dict:
    figure_spec = spec["figures"]
    figures_dir = set_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    n = len(order)
    full = (1 << n) - 1
    row = float(figure_spec.get("row_height", 9.0))
    label_tips = bool(figure_spec.get("label_tips", False))
    written: dict[str, str] = {}

    for pair in figure_spec.get("tanglegram_pairs", []):
        a, b = pair
        svg = render_tanglegram(
            projected[a], projected[b], a, b,
            f"{a}: {trees[a].method}",
            f"{b}: {trees[b].method}",
            f"Tanglegram, {spec['name']} ({n} tips): {a} against {b}",
            row=row, label_tips=label_tips)
        path = figures_dir / f"tanglegram.{a}-vs-{b}.svg"
        path.write_text(svg)
        written[f"tanglegram_{a}_vs_{b}"] = str(path)

    labels = {tip: " / ".join(meta.type_names[tip]) for tip in order if tip in meta.type_names}
    for tree_name in figure_spec.get("annotated_trees", []):
        others = [t for t in splits if t != tree_name
                  and trees[t].role == "primary_comparator"]
        masks = conflict_masks_for(tree_name, splits, others, full)
        summary = tree_split_summary(trees[tree_name], splits[tree_name], n,
                                     int(config.get("backbone_min_smaller_side_tips", 10)))
        caption = (f"{trees[tree_name].method}. {summary['retained_splits']} of "
                   f"{summary['non_trivial_splits']} non-trivial splits retained; "
                   f"{len(masks)} of those are incompatible with a split another method also "
                   f"retains. Retention here means: {summary['retention_means']}")
        svg = render_annotated_tree(
            projected[tree_name], splits[tree_name], index, masks,
            f"{spec['name']} ({n} tips): {tree_name}, branches coloured by retention",
            caption, labels, row=row)
        path = figures_dir / f"annotated.{tree_name}.svg"
        path.write_text(svg)
        written[f"annotated_{tree_name}"] = str(path)
    return written


# --------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------


def run(config: dict, root: Path, outdir: Path) -> dict:
    manifests = config["manifests"]
    base_tips = [line.strip() for line in
                 (root / manifests["tip_ids"]).read_text().splitlines() if line.strip()]
    if len(set(base_tips)) != len(base_tips):
        raise ComparisonError("the tip manifest contains duplicate ids")
    base_tips = sorted(base_tips)
    meta, panel = load_metadata(root / manifests["panel"], root / manifests["types"])
    trees = load_trees(config, root)

    verification = {}
    for name, loaded in trees.items():
        equal = sorted(loaded.tips) == base_tips
        verification[name] = {
            "path": str(loaded.path),
            "sha256": loaded.sha256,
            "tips": len(loaded.tips),
            "tip_set_equals_manifest": equal,
            "missing_from_tree": sorted(set(base_tips) - set(loaded.tips))[:10],
            "not_in_manifest": sorted(set(loaded.tips) - set(base_tips))[:10],
        }
        if not equal:
            raise ComparisonError(
                f"tree {name!r} does not carry exactly the manifest tip set: "
                f"{len(set(base_tips) - set(loaded.tips))} missing, "
                f"{len(set(loaded.tips) - set(base_tips))} extra. A comparison on silently "
                "different tip sets manufactures conflict; fix the inputs.")

    outdir.mkdir(parents=True, exist_ok=True)
    sets = {}
    for spec in config["comparison_sets"]:
        sets[spec["name"]] = run_comparison_set(spec, trees, base_tips, meta, config, outdir)

    result = {
        "tool": "scripts/compare_phylogenies.py",
        "schema_version": config.get("schema_version"),
        "run_label": config.get("run_label"),
        "base_tip_manifest": str(root / manifests["tip_ids"]),
        "base_tips": len(base_tips),
        "base_tip_set_sha256": hashlib.sha256(
            ("\n".join(base_tips) + "\n").encode()).hexdigest(),
        "tip_verification": verification,
        "trees": {name: {
            "path": str(loaded.path), "sha256": loaded.sha256, "method": loaded.method,
            "role": loaded.role, "inference_class": loaded.inference_class,
            "distance_units": loaded.distance_units,
            "support_convention": loaded.support.describe(),
            "near_zero_branch_length": loaded.near_zero_branch_length,
            "near_zero_basis": loaded.near_zero_basis,
            "provenance": loaded.provenance,
        } for name, loaded in trees.items()},
        "shared_input_caveat": (
            "All comparator trees were built from the same Task 06 `verticall repair` "
            "derivatives, so no agreement between them is independent of that step. Two of the "
            "three primary comparators are whole-genome distance methods, so their agreement is "
            "weaker evidence than either agreeing with the marker ML tree. The 17 September "
            "meeting adopted the three as complementary evidence: complementary is not "
            "independent, and no method is ranked here."),
        "comparison_sets": sets,
    }
    (outdir / "tree-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    write_tsv(outdir / "summary.tsv", summary_rows(result))
    return result


def summary_rows(result: dict) -> list[dict]:
    rows = []
    for set_name, data in result["comparison_sets"].items():
        blocks = [("each_method_own_rule", data["comparisons"]),
                  ("matched_length_only_rule",
                   data["matched_rule_comparisons"]["comparisons"])]
        for rule, block in blocks:
            for key, comparison in block.items():
                a, b = comparison["trees"]
                forward = comparison[f"{a}_against_{b}"]
                reverse = comparison[f"{b}_against_{a}"]
                rows.append({
                    "comparison_set": set_name,
                    "collapsing_rule": rule,
                    "tips": comparison["tips"],
                    "pair": key,
                    "pair_kind": data["pair_kind"][key],
                    "retained_splits_a": forward["retained_splits"],
                    "retained_splits_b": reverse["retained_splits"],
                    "shared_retained_splits":
                        comparison["retained_split_comparison"]["shared"],
                    "percent_a_retained_present_in_b":
                        forward["percent_present_in_other_tree"],
                    "percent_b_retained_present_in_a":
                        reverse["percent_present_in_other_tree"],
                    "a_retained_unresolved_in_b": forward["unresolved_in_other_tree"],
                    "b_retained_unresolved_in_a": reverse["unresolved_in_other_tree"],
                    "a_retained_conflicting_in_b":
                        forward["conflicting_with_a_retained_split"],
                    "b_retained_conflicting_in_a":
                        reverse["conflicting_with_a_retained_split"],
                    "percent_a_retained_conflicting_in_b":
                        forward["percent_conflicting_with_a_retained_split"],
                    "percent_b_retained_conflicting_in_a":
                        reverse["percent_conflicting_with_a_retained_split"],
                    "a_backbone_retained": forward["backbone_retained_splits"],
                    "b_backbone_retained": reverse["backbone_retained_splits"],
                    "a_backbone_conflicting_in_b":
                        forward["backbone_conflicting_with_a_retained_split"],
                    "b_backbone_conflicting_in_a":
                        reverse["backbone_conflicting_with_a_retained_split"],
                    "mutually_retained_conflicting_pairs":
                        comparison["mutually_retained_conflicting_split_pairs"],
                    "mutually_retained_conflicting_pairs_backbone":
                        comparison["mutually_retained_conflicting_split_pairs_backbone"],
                    "type_anchored_contested_bipartitions":
                        comparison["type_anchored_contested_bipartitions"],
                    "raw_splits_a": comparison["raw_split_comparison"]["splits_in_a"],
                    "raw_splits_b": comparison["raw_split_comparison"]["splits_in_b"],
                    "raw_robinson_foulds":
                        comparison["raw_split_comparison"]["raw_robinson_foulds"],
                    "normalised_robinson_foulds":
                        comparison["raw_split_comparison"]["normalised_robinson_foulds"],
                    "percent_of_a_raw_splits_absent_from_b":
                        comparison["raw_split_comparison"]["percent_of_a_splits_absent_from_b"],
                })
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd(),
                        help="directory the configuration's relative paths resolve against")
    parser.add_argument("--allow-existing-output", action="store_true",
                        help="write into an existing output directory; off by default so a run "
                             "cannot silently overwrite another")
    args = parser.parse_args(argv)

    if args.output.exists() and any(args.output.iterdir()) and not args.allow_existing_output:
        raise SystemExit(f"refusing to write into the non-empty {args.output}")
    config = json.loads(args.config.read_text())
    try:
        result = run(config, args.root, args.output)
    except ComparisonError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({
        "run_label": result["run_label"],
        "comparison_sets": {name: {"tips": data["tips"],
                                   "trees": data["trees_compared"]}
                            for name, data in result["comparison_sets"].items()},
        "output": str(args.output),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
