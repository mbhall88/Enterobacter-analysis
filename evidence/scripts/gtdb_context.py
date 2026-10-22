#!/usr/bin/env python3
"""Map project manifests onto a pinned GTDB release and extract a contextual subtree.

Task 05 of the Enterobacter analysis execution plan.

What this does
--------------
1. Reads a pinned GTDB reference tree (Newick, decorated internal labels), the release
   species-cluster table (the accession crosswalk) and one or more project manifests.
2. Assigns every manifest record exactly one **level of evidence** for GTDB membership:

   ``exact_tip``               the record's own accession.version is a tip of the tree
   ``alias_tip``               its GCF/GCA counterpart at the same version is a tip
   ``version_tip``             the same accession numeric at a different version is a tip
   ``clustered_representative`` the genome is in the release but is not a tip; it is
                               clustered into a representative genome that is a tip
   ``absent_or_unknown``       no record of the genome in the supplied release files

   These are five different claims and are never collapsed into one coverage number.
   ``match_basis`` records how the accession was matched (exact / alias / version) so the
   two axes can be cross-tabulated.

3. Prunes the reference tree to a contextual tip set and verifies, programmatically, that
   every tip of the pruned tree is a named tip of the source tree and that patristic
   distances are preserved.
4. Reports how the records distribute across GTDB species clusters, and where candidate
   outgroups sit relative to the clade containing the project's genomes.

This is GTDB context. It is not a substitute for a de novo tree of the project's cohort:
the tree here is GTDB's own bac120 reference tree of species representatives, and the great
majority of project genomes are not tips of it.

The Newick parser is implemented here rather than taken from DendroPy/ETE because neither is
installed in the target interpreter, and because GTDB's decorated internal labels
(``'100.0:f__Enterobacteriaceae'``) are quoted and colon-bearing, which is exactly the case a
naive parser mangles. Parsing is explicit and covered by tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

ACCESSION_RE = re.compile(r"^(?:(GB|RS)_)?(GCA|GCF)_(\d+)\.(\d+)$")
GTDB_PREFIXES = ("GB_", "RS_")
STRUCTURAL = set("(),:;")
WHITESPACE = set(" \t\r\n")


class GtdbContextError(RuntimeError):
    """Raised when an input is malformed or a stated invariant fails."""


# --------------------------------------------------------------------------------------
# Accessions
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Accession:
    """A parsed assembly accession.

    ``db`` is GCA or GCF, ``numeric`` the zero-padded digits, ``version`` the integer suffix.
    The GTDB ``GB_``/``RS_`` prefix is not part of the identity; it is recorded separately.
    """

    db: str
    numeric: str
    version: int
    gtdb_prefix: str = ""

    @property
    def bare(self) -> str:
        return f"{self.db}_{self.numeric}.{self.version}"

    @property
    def unversioned(self) -> str:
        return f"{self.db}_{self.numeric}"

    @property
    def pair_key(self) -> str:
        """Key shared by the GCA and GCF members of one assembly pair, ignoring version."""
        return self.numeric

    @property
    def alias_db(self) -> str:
        return "GCF" if self.db == "GCA" else "GCA"

    def with_db(self, db: str) -> "Accession":
        return Accession(db, self.numeric, self.version)

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.bare


def parse_accession(text: str) -> Accession | None:
    """Parse ``GCA_000016325.1`` or ``RS_GCF_000016325.1``; return None if unparseable."""
    if text is None:
        return None
    match = ACCESSION_RE.match(text.strip())
    if match is None:
        return None
    prefix, db, numeric, version = match.groups()
    return Accession(db, numeric, int(version), f"{prefix}_" if prefix else "")


def strip_gtdb_prefix(label: str) -> str:
    for prefix in GTDB_PREFIXES:
        if label.startswith(prefix):
            return label[len(prefix) :]
    return label


# --------------------------------------------------------------------------------------
# Newick
# --------------------------------------------------------------------------------------


@dataclass
class Tree:
    """A rooted Newick tree stored as parallel arrays; node 0 is the root."""

    parent: list[int]
    children: list[list[int]]
    label: list[str]
    blen: list[float | None]

    @property
    def size(self) -> int:
        return len(self.parent)

    def tips(self) -> list[int]:
        return [i for i in range(self.size) if not self.children[i]]

    def postorder(self) -> list[int]:
        """Iterative post-order over all nodes. Recursion is unusable at GTDB tree depth."""
        order: list[int] = []
        stack: list[tuple[int, bool]] = [(0, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                order.append(node)
                continue
            stack.append((node, True))
            for child in self.children[node]:
                stack.append((child, False))
        return order

    def root_distance(self) -> list[float]:
        """Cumulative branch length from the root to each node (root = 0.0)."""
        dist = [0.0] * self.size
        stack = list(self.children[0])
        for child in stack:
            dist[child] = self.blen[child] or 0.0
        while stack:
            node = stack.pop()
            for child in self.children[node]:
                dist[child] = dist[node] + (self.blen[child] or 0.0)
                stack.append(child)
        return dist

    def ancestors(self, node: int) -> list[int]:
        chain = [node]
        while self.parent[chain[-1]] != -1:
            chain.append(self.parent[chain[-1]])
        return chain

    def mrca(self, nodes: list[int]) -> int:
        if not nodes:
            raise GtdbContextError("mrca() requires at least one node")
        seen = set(self.ancestors(nodes[0]))
        current = nodes[0]
        for node in nodes[1:]:
            chain = self.ancestors(node)
            for candidate in chain:
                if candidate in seen:
                    current = candidate
                    break
            else:  # pragma: no cover - impossible in a connected tree
                raise GtdbContextError("nodes are not in one tree")
            seen = set(self.ancestors(current))
        return current

    def descendant_tips(self, node: int) -> list[int]:
        found: list[int] = []
        stack = [node]
        while stack:
            current = stack.pop()
            if not self.children[current]:
                found.append(current)
            else:
                stack.extend(self.children[current])
        return found


def _read_quoted(text: str, start: int) -> tuple[str, int]:
    """Read a single-quoted Newick label starting at ``text[start] == "'"``.

    Two consecutive quotes inside the label are one literal quote, per the Newick
    convention that GTDB's decorated labels rely on.
    """
    pos = start + 1
    out: list[str] = []
    n = len(text)
    while pos < n:
        char = text[pos]
        if char == "'":
            if pos + 1 < n and text[pos + 1] == "'":
                out.append("'")
                pos += 2
                continue
            return "".join(out), pos + 1
        out.append(char)
        pos += 1
    raise GtdbContextError("unterminated quoted label in Newick input")


def _read_bare(text: str, start: int) -> tuple[str, int]:
    pos = start
    n = len(text)
    out: list[str] = []
    while pos < n:
        char = text[pos]
        if char in STRUCTURAL or char in WHITESPACE:
            break
        # No underscore-to-space conversion: GTDB accession labels use '_' literally.
        out.append(char)
        pos += 1
    return "".join(out), pos


def parse_newick(text: str) -> Tree:
    """Parse one Newick string into a :class:`Tree`.

    Handles quoted labels, labels containing colons (GTDB's ``support:rank__name``
    decoration), branch lengths in scientific notation, and a label on the root.
    """
    parent: list[int] = [-1]
    children: list[list[int]] = [[]]
    label: list[str] = [""]
    blen: list[float | None] = [None]

    def new_child(par: int) -> int:
        node = len(parent)
        parent.append(par)
        children.append([])
        label.append("")
        blen.append(None)
        children[par].append(node)
        return node

    cur = 0
    pos = 0
    n = len(text)
    while pos < n:
        char = text[pos]
        if char in WHITESPACE:
            pos += 1
        elif char == "(":
            cur = new_child(cur)
            pos += 1
        elif char == ",":
            if parent[cur] == -1:
                raise GtdbContextError("comma outside any parenthesis in Newick input")
            cur = new_child(parent[cur])
            pos += 1
        elif char == ")":
            if parent[cur] == -1:
                raise GtdbContextError("unbalanced ')' in Newick input")
            cur = parent[cur]
            pos += 1
        elif char == ";":
            break
        elif char == ":":
            pos += 1
            start = pos
            while pos < n and text[pos] not in STRUCTURAL and text[pos] not in WHITESPACE:
                pos += 1
            raw = text[start:pos]
            try:
                blen[cur] = float(raw)
            except ValueError as exc:
                raise GtdbContextError(f"bad branch length {raw!r}") from exc
        elif char == "'":
            value, pos = _read_quoted(text, pos)
            label[cur] = value
        else:
            value, pos = _read_bare(text, pos)
            label[cur] = value

    if cur != 0:
        raise GtdbContextError("unbalanced parentheses in Newick input")
    return Tree(parent, children, label, blen)


def format_newick(tree: Tree, decimals: int = 6) -> str:
    """Render a :class:`Tree` back to Newick, quoting labels that need it."""

    def render_label(value: str) -> str:
        if value == "":
            return ""
        if any(ch in value for ch in "(),:;' \t") or value != value.strip():
            return "'" + value.replace("'", "''") + "'"
        return value

    def render_length(value: float | None) -> str:
        if value is None:
            return ""
        text = f"{value:.{decimals}f}".rstrip("0").rstrip(".")
        return ":" + (text if text not in ("", "-") else "0")

    parts: dict[int, str] = {}
    for node in tree.postorder():
        if tree.children[node]:
            inner = ",".join(parts.pop(child) for child in tree.children[node])
            text = f"({inner}){render_label(tree.label[node])}"
        else:
            text = render_label(tree.label[node])
        if node != 0:
            text += render_length(tree.blen[node])
        parts[node] = text
    return parts[0] + ";"


@dataclass
class PruneResult:
    tree: Tree
    #: pruned node index -> source node index
    source_of: dict[int, int]
    #: pruned internal node index -> source internal labels collapsed into it, tip-ward first
    collapsed_labels: dict[int, list[str]]


def prune_to_tips(tree: Tree, keep: set[int]) -> PruneResult:
    """Prune ``tree`` to the tips in ``keep``, suppressing unifurcations.

    Branch lengths of suppressed nodes are summed into the surviving branch, so patristic
    distances between retained tips are unchanged. Internal labels of suppressed nodes are
    not silently discarded: they are returned in ``collapsed_labels``.
    """
    missing = [node for node in keep if tree.children[node]]
    if missing:
        raise GtdbContextError(f"prune_to_tips was given {len(missing)} non-tip node(s)")

    retained = [False] * tree.size
    kept_child_count = [0] * tree.size
    for node in tree.postorder():
        if not tree.children[node]:
            retained[node] = node in keep
        else:
            count = sum(1 for child in tree.children[node] if retained[child])
            kept_child_count[node] = count
            retained[node] = count > 0

    if not retained[0]:
        raise GtdbContextError("no requested tip is present in the tree")

    # A source node survives as a node of the pruned tree when it is a kept tip, or an
    # internal node with at least two kept children, or the root of the kept subtree.
    def survives(node: int) -> bool:
        if not retained[node]:
            return False
        if not tree.children[node]:
            return True
        return kept_child_count[node] >= 2

    # The pruned root is the deepest surviving ancestor chain start: walk down from the
    # source root through single kept children until a surviving node is reached.
    source_root = 0
    while not survives(source_root):
        next_nodes = [child for child in tree.children[source_root] if retained[child]]
        if len(next_nodes) != 1:  # pragma: no cover - guarded by survives()
            raise GtdbContextError("inconsistent pruning state at the root")
        source_root = next_nodes[0]

    parent: list[int] = [-1]
    children: list[list[int]] = [[]]
    label: list[str] = [tree.label[source_root]]
    blen: list[float | None] = [None]
    source_of: dict[int, int] = {0: source_root}
    collapsed: dict[int, list[str]] = {}

    # (source node, pruned parent index, accumulated branch length, collapsed labels)
    stack: list[tuple[int, int, float, list[str]]] = []
    for child in tree.children[source_root]:
        if retained[child]:
            stack.append((child, 0, tree.blen[child] or 0.0, []))

    while stack:
        node, pruned_parent, acc_len, dropped = stack.pop()
        if survives(node):
            index = len(parent)
            parent.append(pruned_parent)
            children.append([])
            label.append(tree.label[node])
            blen.append(acc_len)
            children[pruned_parent].append(index)
            source_of[index] = node
            if dropped:
                collapsed[index] = list(dropped)
            for child in tree.children[node]:
                if retained[child]:
                    stack.append((child, index, tree.blen[child] or 0.0, []))
        else:
            # A retained internal node with exactly one kept child: suppress it, carrying
            # its branch length and remembering its decorated label.
            next_nodes = [child for child in tree.children[node] if retained[child]]
            if len(next_nodes) != 1:  # pragma: no cover
                raise GtdbContextError("inconsistent pruning state")
            child = next_nodes[0]
            carried = list(dropped)
            if tree.label[node]:
                carried.append(tree.label[node])
            stack.append((child, pruned_parent, acc_len + (tree.blen[child] or 0.0), carried))

    # The traversal above uses a LIFO stack, so siblings come out reversed. Source node
    # indices are assigned left to right during parsing, so sorting by them restores the
    # source tree's sibling order and makes the emitted Newick deterministic.
    for kids in children:
        kids.sort(key=lambda node: source_of[node])

    return PruneResult(Tree(parent, children, label, blen), source_of, collapsed)


# --------------------------------------------------------------------------------------
# GTDB release index
# --------------------------------------------------------------------------------------


@dataclass
class GtdbEntry:
    accession: Accession
    representative: str
    species: str
    taxonomy: str
    is_representative: bool


@dataclass
class GtdbIndex:
    """Accession-level index of one GTDB release, built from the species-cluster table."""

    by_accession: dict[str, GtdbEntry] = field(default_factory=dict)
    by_pair_key: dict[str, list[GtdbEntry]] = field(default_factory=lambda: defaultdict(list))
    representative_species: dict[str, str] = field(default_factory=dict)
    representative_taxonomy: dict[str, str] = field(default_factory=dict)
    cluster_size: dict[str, int] = field(default_factory=dict)
    qc_failed: dict[str, str] = field(default_factory=dict)
    n_clusters: int = 0
    n_genomes: int = 0

    def lookup_exact(self, accession: Accession) -> GtdbEntry | None:
        return self.by_accession.get(accession.bare)

    def lookup_alias(self, accession: Accession) -> GtdbEntry | None:
        return self.by_accession.get(accession.with_db(accession.alias_db).bare)

    def lookup_other_version(self, accession: Accession) -> GtdbEntry | None:
        """Best different-version match: same db preferred, then the alias db."""
        candidates = [
            entry
            for entry in self.by_pair_key.get(accession.pair_key, [])
            if entry.accession.version != accession.version
        ]
        if not candidates:
            return None
        candidates.sort(
            key=lambda e: (e.accession.db != accession.db, -e.accession.version)
        )
        return candidates[0]


def load_sp_clusters(path: Path, index: GtdbIndex | None = None) -> GtdbIndex:
    """Load ``sp_clusters_r232.tsv``: representatives plus their clustered genomes."""
    index = index or GtdbIndex()
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        needed = [
            "Representative genome",
            "GTDB species",
            "GTDB taxonomy",
            "No. clustered genomes",
            "Clustered genomes",
        ]
        for column in needed:
            if column not in header:
                raise GtdbContextError(f"{path}: missing column {column!r}")
        pos = {name: header.index(name) for name in needed}
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            rep_raw = fields[pos["Representative genome"]]
            species = fields[pos["GTDB species"]]
            taxonomy = fields[pos["GTDB taxonomy"]]
            rep_bare = strip_gtdb_prefix(rep_raw)
            index.n_clusters += 1
            index.representative_species[rep_bare] = species
            index.representative_taxonomy[rep_bare] = taxonomy
            try:
                index.cluster_size[rep_bare] = int(fields[pos["No. clustered genomes"]])
            except (ValueError, IndexError):
                index.cluster_size[rep_bare] = 0

            members = [rep_raw]
            clustered = fields[pos["Clustered genomes"]] if len(fields) > pos["Clustered genomes"] else ""
            if clustered:
                members.extend(clustered.split(","))
            for raw in members:
                raw = raw.strip()
                if not raw:
                    continue
                accession = parse_accession(raw)
                if accession is None:
                    continue
                entry = GtdbEntry(
                    accession=accession,
                    representative=rep_bare,
                    species=species,
                    taxonomy=taxonomy,
                    is_representative=(accession.bare == rep_bare),
                )
                index.by_accession[accession.bare] = entry
                index.by_pair_key[accession.pair_key].append(entry)
                index.n_genomes += 1
    return index


def load_qc_failed(path: Path, index: GtdbIndex) -> GtdbIndex:
    """Load ``qc_failed_r232.tsv`` so that an absence can say *why* where GTDB knows."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if "Accession" not in header:
            raise GtdbContextError(f"{path}: missing column 'Accession'")
        acc_pos = header.index("Accession")
        fail_columns = [(i, name) for i, name in enumerate(header) if name.startswith("Failed ")]
        exclude_pos = header.index("NCBI exclude from RefSeq") if "NCBI exclude from RefSeq" in header else None
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            accession = parse_accession(fields[acc_pos])
            if accession is None:
                continue
            reasons = [
                name[len("Failed ") :]
                for i, name in fail_columns
                if i < len(fields) and fields[i].strip() == "1"
            ]
            if exclude_pos is not None and exclude_pos < len(fields) and fields[exclude_pos].strip():
                reasons.append(f"ncbi_exclude:{fields[exclude_pos].strip()}")
            index.qc_failed[accession.bare] = ";".join(reasons) if reasons else "listed_without_flag"
    return index


# --------------------------------------------------------------------------------------
# Manifest records
# --------------------------------------------------------------------------------------


@dataclass
class Record:
    role: str
    record_id: str
    accession_text: str
    accession: Accession | None
    source: str
    extra: dict[str, str] = field(default_factory=dict)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        text = handle.read()
    # The 11 September membership table has CRLF line endings; normalise on read.
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    lines = [line for line in lines if line != ""]
    if not lines:
        raise GtdbContextError(f"{path}: empty file")
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) < len(header):
            fields += [""] * (len(header) - len(fields))
        rows.append(dict(zip(header, fields)))
    return header, rows


def load_manifest(
    path: Path,
    role: str,
    accession_column: str,
    id_column: str | None,
    keep_columns: tuple[str, ...] = (),
) -> list[Record]:
    header, rows = read_tsv(path)
    if accession_column not in header:
        raise GtdbContextError(
            f"{path}: no column {accession_column!r}; available: {', '.join(header)}"
        )
    records: list[Record] = []
    for number, row in enumerate(rows, start=2):
        raw = (row.get(accession_column) or "").strip()
        if not raw:
            continue
        record_id = (row.get(id_column) or "").strip() if id_column and id_column in header else ""
        if not record_id:
            record_id = raw
        extra = {name: row.get(name, "") for name in keep_columns if name in header}
        records.append(
            Record(
                role=role,
                record_id=record_id,
                accession_text=raw,
                accession=parse_accession(raw),
                source=f"{path}:{number}",
                extra=extra,
            )
        )
    return records


# --------------------------------------------------------------------------------------
# Evidence assignment
# --------------------------------------------------------------------------------------

EVIDENCE_LEVELS = (
    "exact_tip",
    "alias_tip",
    "version_tip",
    "clustered_representative",
    "absent_or_unknown",
)


@dataclass
class Assignment:
    record: Record
    evidence_level: str
    match_basis: str
    matched_accession: str
    representative: str
    representative_is_tip: bool
    gtdb_species: str
    gtdb_taxonomy: str
    absence_basis: str
    note: str


def assign_evidence(
    record: Record,
    index: GtdbIndex,
    tip_accessions: set[str],
) -> Assignment:
    """Assign one record its single level of evidence for GTDB membership."""
    empty = Assignment(
        record=record,
        evidence_level="absent_or_unknown",
        match_basis="none",
        matched_accession="",
        representative="",
        representative_is_tip=False,
        gtdb_species="",
        gtdb_taxonomy="",
        absence_basis="",
        note="",
    )
    if record.accession is None:
        empty.absence_basis = "unparseable_accession"
        empty.note = f"could not parse {record.accession_text!r}"
        return empty

    lookups = (
        ("exact_accession_version", index.lookup_exact(record.accession)),
        ("gca_gcf_alias", index.lookup_alias(record.accession)),
        ("different_version", index.lookup_other_version(record.accession)),
    )
    tip_level = {
        "exact_accession_version": "exact_tip",
        "gca_gcf_alias": "alias_tip",
        "different_version": "version_tip",
    }

    # Identity strength decides first: an exact accession.version match is a stronger claim
    # about *which genome this is* than a same-assembly alias, and an alias is stronger than
    # a different version of the assembly. Only once the identity is fixed does the record
    # split into "this genome is itself a tip" and "this genome is clustered into a tip".
    for basis, entry in lookups:
        if entry is None:
            continue
        is_tip = entry.accession.bare in tip_accessions
        return Assignment(
            record=record,
            evidence_level=tip_level[basis] if is_tip else "clustered_representative",
            match_basis=basis,
            matched_accession=entry.accession.bare,
            representative=entry.representative,
            representative_is_tip=entry.representative in tip_accessions,
            gtdb_species=entry.species,
            gtdb_taxonomy=entry.taxonomy,
            absence_basis="",
            note=(
                ""
                if is_tip or entry.representative in tip_accessions
                else "representative is not a tip of this bac120 tree"
            ),
        )

    for basis, probe in (
        ("exact_accession_version", record.accession.bare),
        ("gca_gcf_alias", record.accession.with_db(record.accession.alias_db).bare),
    ):
        if probe in index.qc_failed:
            empty.absence_basis = "qc_failed_in_release"
            empty.match_basis = basis
            empty.matched_accession = probe
            empty.note = index.qc_failed[probe]
            return empty
    empty.absence_basis = "not_in_supplied_release_files"
    return empty


def genus_of(taxonomy: str) -> str:
    for part in taxonomy.replace("; ", ";").split(";"):
        if part.startswith("g__"):
            return part[3:]
    return ""


def species_of(species_field: str) -> str:
    return species_field[3:] if species_field.startswith("s__") else species_field


# --------------------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------------------


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if value is None else str(value) for value in row) + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------------------


def verify_pruned_tree(
    source: Tree,
    source_tip_of_label: dict[str, int],
    pruned_newick: str,
    requested_labels: set[str],
    sample_pairs: int,
    seed: int,
) -> dict[str, object]:
    """Re-parse the emitted Newick and prove it against the source tree.

    Checks, in order: the emitted text parses; its tip label set equals the requested set;
    every emitted tip label is a named tip of the source tree; labels are unique; and a
    random sample of patristic distances between retained tips is identical in both trees.
    """
    reparsed = parse_newick(pruned_newick)
    emitted = [reparsed.label[node] for node in reparsed.tips()]
    emitted_set = set(emitted)
    result: dict[str, object] = {
        "reparsed_tip_count": len(emitted),
        "requested_tip_count": len(requested_labels),
        "tip_labels_unique": len(emitted_set) == len(emitted),
        "tip_set_equals_requested": emitted_set == requested_labels,
        "unnamed_tips": sum(1 for label in emitted if not label),
    }
    unmapped = sorted(label for label in emitted_set if label not in source_tip_of_label)
    result["tips_not_in_source_tree"] = len(unmapped)
    result["tips_not_in_source_tree_examples"] = unmapped[:5]

    source_dist = source.root_distance()
    pruned_dist = reparsed.root_distance()
    pruned_tip_of_label = {reparsed.label[node]: node for node in reparsed.tips()}

    def patristic(tree: Tree, dist: list[float], a: int, b: int) -> float:
        node = tree.mrca([a, b])
        return dist[a] + dist[b] - 2 * dist[node]

    rng = random.Random(seed)
    # Only labels present in both trees can be compared; an emitted tip missing from the
    # source is already a recorded failure above, and must not crash the check.
    labels = sorted(emitted_set & set(source_tip_of_label))
    checked = 0
    max_delta = 0.0
    worst: tuple[str, str] | None = None
    if len(labels) >= 2:
        for _ in range(sample_pairs):
            a_label, b_label = rng.sample(labels, 2)
            delta = abs(
                patristic(source, source_dist, source_tip_of_label[a_label], source_tip_of_label[b_label])
                - patristic(reparsed, pruned_dist, pruned_tip_of_label[a_label], pruned_tip_of_label[b_label])
            )
            checked += 1
            if delta > max_delta:
                max_delta = delta
                worst = (a_label, b_label)
    result["patristic_pairs_checked"] = checked
    result["patristic_max_abs_difference"] = max_delta
    result["patristic_worst_pair"] = list(worst) if worst else []
    result["patristic_preserved"] = max_delta < 1e-5
    result["passed"] = bool(
        result["tip_labels_unique"]
        and result["tip_set_equals_requested"]
        and result["tips_not_in_source_tree"] == 0
        and result["unnamed_tips"] == 0
        and result["patristic_preserved"]
    )
    return result


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def parse_manifest_spec(spec: str) -> tuple[str, Path, str, str | None]:
    """``role=path#accession_column[,id_column]``."""
    if "=" not in spec:
        raise argparse.ArgumentTypeError(f"--manifest needs role=path#column, got {spec!r}")
    role, rest = spec.split("=", 1)
    if "#" not in rest:
        raise argparse.ArgumentTypeError(f"--manifest needs a #column suffix, got {spec!r}")
    path_text, columns = rest.rsplit("#", 1)
    if "," in columns:
        accession_column, id_column = columns.split(",", 1)
    else:
        accession_column, id_column = columns, None
    return role, Path(path_text), accession_column, id_column


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Map project manifests onto a pinned GTDB release and prune a contextual subtree.",
    )
    parser.add_argument("--study", type=Path, required=True, help="study manifest TSV")
    parser.add_argument("--study-role", default="study", help="role label for study records")
    parser.add_argument("--study-accession-column", default="accession_version")
    parser.add_argument("--study-id-column", default="entity_id")
    parser.add_argument(
        "--manifest",
        action="append",
        default=[],
        metavar="ROLE=PATH#ACCESSION_COL[,ID_COL]",
        help="additional manifest, repeatable (e.g. type_anchor=manifests/types.tsv#selected_accession_version,taxon_name)",
    )
    parser.add_argument("--taxonomy", type=Path, required=True, help="pinned GTDB taxonomy TSV (representatives)")
    parser.add_argument("--tree", type=Path, required=True, help="pinned GTDB reference tree (Newick)")
    parser.add_argument("--tree-md5", default="", help="expected MD5 of --tree; verified when given")
    parser.add_argument(
        "--sp-clusters",
        type=Path,
        required=True,
        help="release species-cluster table; this is the accession crosswalk",
    )
    parser.add_argument("--qc-failed", type=Path, default=None, help="release qc_failed table (optional)")
    parser.add_argument("--release", default="", help="release label recorded in provenance, e.g. R232")
    parser.add_argument(
        "--focus-genus",
        action="append",
        default=[],
        help="GTDB genus to treat as the focus group; repeatable (default: Enterobacter, Enterobacter_B, Enterobacter_D)",
    )
    parser.add_argument(
        "--context-rank-prefix",
        default="f__Enterobacteriaceae",
        help="taxonomy substring defining the contextual subtree scope",
    )
    parser.add_argument(
        "--placement-role",
        default="outgroup_candidate",
        help="manifest role whose records are placed relative to the project clade",
    )
    parser.add_argument(
        "--placement-species",
        action="append",
        default=[],
        metavar="GTDB_SPECIES",
        help=(
            "GTDB species name to place even though no project manifest supplies a genome for it "
            "(e.g. 's__Kosakonia cowanii'); repeatable"
        ),
    )
    parser.add_argument("--sample-pairs", type=int, default=250, help="patristic verification sample size")
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--output", type=Path, required=True, help="run directory for outputs")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    focus_genera = set(args.focus_genus) or {"Enterobacter", "Enterobacter_B", "Enterobacter_D"}
    out = args.output
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    started = time.time()
    log_lines: list[str] = []

    def log(message: str) -> None:
        stamp = f"[{time.time() - started:8.2f}s] {message}"
        log_lines.append(stamp)
        print(stamp, file=sys.stderr, flush=True)

    provenance: dict[str, object] = {
        "release": args.release,
        "focus_genera": sorted(focus_genera),
        "context_rank_prefix": args.context_rank_prefix,
        "inputs": {},
    }

    # --- pinned inputs -----------------------------------------------------------------
    for name, path in (
        ("tree", args.tree),
        ("taxonomy", args.taxonomy),
        ("sp_clusters", args.sp_clusters),
        ("qc_failed", args.qc_failed),
        ("study", args.study),
    ):
        if path is None:
            continue
        provenance["inputs"][name] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    tree_md5 = md5_file(args.tree)
    provenance["inputs"]["tree"]["md5"] = tree_md5
    if args.tree_md5:
        provenance["tree_md5_expected"] = args.tree_md5
        provenance["tree_md5_verified"] = tree_md5 == args.tree_md5
        if tree_md5 != args.tree_md5:
            raise GtdbContextError(
                f"--tree MD5 {tree_md5} does not match the pinned {args.tree_md5}"
            )
        log(f"tree MD5 verified against the refpkg manifest: {tree_md5}")

    # --- release index -----------------------------------------------------------------
    log(f"loading species clusters from {args.sp_clusters}")
    index = load_sp_clusters(args.sp_clusters)
    log(f"  {index.n_clusters:,} clusters, {index.n_genomes:,} genome accessions indexed")
    if args.qc_failed is not None:
        load_qc_failed(args.qc_failed, index)
        log(f"  {len(index.qc_failed):,} QC-failed accessions indexed")

    taxonomy_reps: dict[str, str] = {}
    with args.taxonomy.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            taxonomy_reps[strip_gtdb_prefix(fields[0])] = fields[1]
    log(f"  {len(taxonomy_reps):,} representative taxonomy rows")
    provenance["taxonomy_rows"] = len(taxonomy_reps)
    provenance["sp_cluster_rows"] = index.n_clusters
    provenance["sp_cluster_genomes"] = index.n_genomes

    # --- tree --------------------------------------------------------------------------
    log(f"parsing {args.tree}")
    parse_start = time.time()
    tree = parse_newick(args.tree.read_text(encoding="utf-8"))
    tip_nodes = tree.tips()
    log(f"  parsed {tree.size:,} nodes, {len(tip_nodes):,} tips in {time.time() - parse_start:.1f}s")
    source_tip_of_label = {tree.label[node]: node for node in tip_nodes}
    if len(source_tip_of_label) != len(tip_nodes):
        raise GtdbContextError("source tree has duplicate tip labels")
    tip_node_of_accession: dict[str, int] = {}
    for node in tip_nodes:
        tip_node_of_accession[strip_gtdb_prefix(tree.label[node])] = node
    tip_accessions = set(tip_node_of_accession)
    provenance["tree_nodes"] = tree.size
    provenance["tree_tips"] = len(tip_nodes)

    # --- records -----------------------------------------------------------------------
    records = load_manifest(
        args.study,
        args.study_role,
        args.study_accession_column,
        args.study_id_column,
        keep_columns=("entity_role", "usage_restrictions", "historical_ncbi", "historical_kraken_chromosome"),
    )
    log(f"{args.study}: {len(records):,} {args.study_role} records")
    manifest_counts = {args.study_role: len(records)}
    for spec in args.manifest:
        role, path, accession_column, id_column = parse_manifest_spec(spec)
        extra = load_manifest(
            path,
            role,
            accession_column,
            id_column,
            keep_columns=("taxon_name", "rank", "reference_role", "selection_status", "context_role", "usage_restrictions"),
        )
        records.extend(extra)
        manifest_counts[role] = manifest_counts.get(role, 0) + len(extra)
        log(f"{path}: {len(extra):,} {role} records")
        provenance["inputs"][f"manifest:{role}"] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "accession_column": accession_column,
        }
    provenance["manifest_record_counts"] = manifest_counts

    # --- evidence ----------------------------------------------------------------------
    assignments = [assign_evidence(record, index, tip_accessions) for record in records]
    level_counts: dict[str, Counter] = defaultdict(Counter)
    for assignment in assignments:
        level_counts[assignment.record.role][assignment.evidence_level] += 1
    for role, counts in level_counts.items():
        log(f"{role}: " + ", ".join(f"{level}={counts[level]}" for level in EVIDENCE_LEVELS))

    membership_rows = []
    for assignment in assignments:
        record = assignment.record
        membership_rows.append(
            [
                record.role,
                record.record_id,
                record.accession_text,
                assignment.evidence_level,
                assignment.match_basis,
                assignment.matched_accession,
                assignment.representative,
                "yes" if assignment.representative_is_tip else "no",
                assignment.gtdb_species,
                genus_of(assignment.gtdb_taxonomy),
                assignment.gtdb_taxonomy,
                index.cluster_size.get(assignment.representative, ""),
                assignment.absence_basis,
                assignment.note,
                record.source,
                ";".join(f"{k}={v}" for k, v in record.extra.items() if v),
            ]
        )
    write_tsv(
        out / "tables" / "gtdb-membership.tsv",
        [
            "role",
            "record_id",
            "accession_version",
            "evidence_level",
            "match_basis",
            "matched_gtdb_accession",
            "gtdb_representative",
            "representative_is_tree_tip",
            "gtdb_species",
            "gtdb_genus",
            "gtdb_taxonomy",
            "gtdb_cluster_size",
            "absence_basis",
            "note",
            "source_row",
            "manifest_fields",
        ],
        membership_rows,
    )

    summary_rows = []
    for role in sorted(level_counts):
        for level in EVIDENCE_LEVELS:
            summary_rows.append([role, level, level_counts[role][level]])
        summary_rows.append([role, "TOTAL", sum(level_counts[role].values())])
    write_tsv(out / "tables" / "evidence-level-summary.tsv", ["role", "evidence_level", "records"], summary_rows)

    basis_counter: Counter = Counter()
    for assignment in assignments:
        basis_counter[(assignment.record.role, assignment.evidence_level, assignment.match_basis)] += 1
    write_tsv(
        out / "tables" / "evidence-level-by-match-basis.tsv",
        ["role", "evidence_level", "match_basis", "records"],
        [[role, level, basis, count] for (role, level, basis), count in sorted(basis_counter.items())],
    )

    # --- cluster distribution ----------------------------------------------------------
    cluster_records: dict[str, Counter] = defaultdict(Counter)
    for assignment in assignments:
        if assignment.evidence_level == "absent_or_unknown":
            continue
        cluster_records[assignment.representative][assignment.record.role] += 1
    roles = sorted({assignment.record.role for assignment in assignments})
    distribution_rows = []
    for representative, counts in cluster_records.items():
        species = index.representative_species.get(representative, "")
        taxonomy = index.representative_taxonomy.get(representative, "")
        distribution_rows.append(
            [
                species_of(species),
                genus_of(taxonomy),
                representative,
                "yes" if representative in tip_accessions else "no",
                index.cluster_size.get(representative, ""),
                sum(counts.values()),
                *[counts.get(role, 0) for role in roles],
                taxonomy,
            ]
        )
    distribution_rows.sort(key=lambda row: (-int(row[5]), str(row[0])))
    write_tsv(
        out / "tables" / "gtdb-cluster-distribution.tsv",
        [
            "gtdb_species",
            "gtdb_genus",
            "gtdb_representative",
            "representative_is_tree_tip",
            "gtdb_cluster_size",
            "project_records",
            *[f"records_{role}" for role in roles],
            "gtdb_taxonomy",
        ],
        distribution_rows,
    )

    # --- how GTDB's cluster names relate to the project's nomenclatural anchors ---------
    # GTDB species names come in three shapes, and the difference is the whole DEC-004
    # question: a plain Latin binomial, a suffixed split of one (cloacae_M), and an
    # accession-derived placeholder (sp000016325). Crossing that with "does any of the
    # project's nomenclatural type anchors fall in this cluster" says which project genomes
    # sit in a GTDB cluster that no type strain anchors.
    def name_class(species: str) -> str:
        epithet = species_of(species).split(" ", 1)[-1] if " " in species_of(species) else ""
        if re.fullmatch(r"sp\d+", epithet):
            return "placeholder_sp"
        if re.search(r"_[A-Z]$", epithet):
            return "suffixed_split"
        return "latin_binomial"

    class_rows = []
    class_counter: dict[tuple[str, str, str], Counter] = defaultdict(Counter)
    for representative, counts in cluster_records.items():
        species = index.representative_species.get(representative, "")
        genus = genus_of(index.representative_taxonomy.get(representative, ""))
        anchored = "with_project_type_anchor" if counts.get("type_anchor", 0) else "no_project_type_anchor"
        in_focus = "focus_genus" if genus in focus_genera else "other_genus"
        bucket = class_counter[(in_focus, name_class(species), anchored)]
        bucket["clusters"] += 1
        for role in roles:
            bucket[role] += counts.get(role, 0)
    for (in_focus, shape, anchored), bucket in sorted(class_counter.items()):
        class_rows.append(
            [in_focus, shape, anchored, bucket["clusters"], *[bucket.get(role, 0) for role in roles]]
        )
    write_tsv(
        out / "tables" / "cluster-name-class-summary.tsv",
        [
            "genus_class",
            "gtdb_name_shape",
            "project_type_anchor",
            "gtdb_clusters",
            *[f"records_{role}" for role in roles],
        ],
        class_rows,
    )

    # --- contextual subtrees -----------------------------------------------------------
    context_reps = {
        representative
        for representative, taxonomy in index.representative_taxonomy.items()
        if args.context_rank_prefix in taxonomy and representative in tip_accessions
    }
    project_reps = {
        assignment.representative
        for assignment in assignments
        if assignment.representative and assignment.representative in tip_accessions
    }
    focus_reps = {
        representative
        for representative, taxonomy in index.representative_taxonomy.items()
        if genus_of(taxonomy) in focus_genera and representative in tip_accessions
    }
    log(
        f"context scope {args.context_rank_prefix}: {len(context_reps):,} tips; "
        f"project-linked: {len(project_reps):,}; focus genera: {len(focus_reps):,}"
    )

    subtrees = {
        "context-subtree": sorted(context_reps | project_reps),
        "focus-subtree": sorted(focus_reps | project_reps),
    }
    verification: dict[str, object] = {}
    tip_map_rows: list[list[str]] = []
    records_by_rep: dict[str, list[Assignment]] = defaultdict(list)
    for assignment in assignments:
        if assignment.representative:
            records_by_rep[assignment.representative].append(assignment)

    for name, keep_accessions in subtrees.items():
        keep_nodes = {tip_node_of_accession[accession] for accession in keep_accessions}
        prune_start = time.time()
        pruned = prune_to_tips(tree, keep_nodes)
        newick = format_newick(pruned.tree)
        (out / f"{name}.nwk").write_text(newick + "\n", encoding="utf-8")
        log(f"{name}: {len(keep_nodes):,} tips, {pruned.tree.size:,} nodes, {time.time() - prune_start:.1f}s")

        requested_labels = {tree.label[node] for node in keep_nodes}
        verification[name] = verify_pruned_tree(
            tree, source_tip_of_label, newick, requested_labels, args.sample_pairs, args.seed
        )
        log(f"{name} verification: {'PASS' if verification[name]['passed'] else 'FAIL'}")

        collapsed_rows = []
        for index_pruned, labels in sorted(pruned.collapsed_labels.items()):
            collapsed_rows.append(
                [
                    name,
                    index_pruned,
                    pruned.tree.label[index_pruned],
                    len(labels),
                    " | ".join(labels),
                ]
            )
        write_tsv(
            out / "tables" / f"{name}-collapsed-internal-labels.tsv",
            ["subtree", "pruned_node_index", "pruned_node_label", "collapsed_count", "collapsed_source_labels"],
            collapsed_rows,
        )

        for pruned_node in pruned.tree.tips():
            source_node = pruned.source_of[pruned_node]
            source_label = tree.label[source_node]
            accession = strip_gtdb_prefix(source_label)
            linked = records_by_rep.get(accession, [])
            tip_map_rows.append(
                [
                    name,
                    pruned.tree.label[pruned_node],
                    source_label,
                    source_node,
                    accession,
                    index.representative_species.get(accession, ""),
                    genus_of(index.representative_taxonomy.get(accession, "")),
                    len(linked),
                    ";".join(sorted({a.record.role for a in linked})),
                ]
            )

    write_tsv(
        out / "tables" / "source-tip-map.tsv",
        [
            "subtree",
            "pruned_tip_label",
            "source_tip_label",
            "source_tip_node_index",
            "accession_version",
            "gtdb_species",
            "gtdb_genus",
            "project_records_represented",
            "project_roles",
        ],
        tip_map_rows,
    )

    # --- clade definitions, monophyly, and candidate-outgroup placement -----------------
    # A rooting decision needs more than one clade. "Outside Enterobacter" and "outside
    # every genome in the project's cohort" are different requirements, and GTDB splits the
    # genus name across three genera, so each definition is evaluated separately rather than
    # lumped into one MRCA.
    dist = tree.root_distance()

    def patristic(a: int, b: int) -> float:
        return dist[a] + dist[b] - 2 * dist[tree.mrca([a, b])]

    def labels_on_path(start: int, stop: int) -> list[str]:
        out_labels = []
        node = start
        while node != stop and node != -1:
            node = tree.parent[node]
            if node != -1 and tree.label[node]:
                out_labels.append(tree.label[node])
        return out_labels

    genus_tips: dict[str, list[int]] = defaultdict(list)
    for representative, taxonomy in index.representative_taxonomy.items():
        node = tip_node_of_accession.get(representative)
        if node is not None:
            genus_tips[genus_of(taxonomy)].append(node)

    project_tip_nodes = sorted({tip_node_of_accession[a] for a in project_reps})
    project_focus_tips = sorted(
        node
        for node in project_tip_nodes
        if genus_of(index.representative_taxonomy.get(strip_gtdb_prefix(tree.label[node]), ""))
        in focus_genera
    )

    clade_definitions: list[tuple[str, str, list[int]]] = []
    for genus in sorted(focus_genera):
        if genus_tips.get(genus):
            clade_definitions.append(
                (f"gtdb_genus:{genus}", f"all R232 representatives of g__{genus}", genus_tips[genus])
            )
    if project_focus_tips:
        clade_definitions.append(
            (
                "project_linked_focus_genera",
                "representatives linked to a project record and inside the focus genera",
                project_focus_tips,
            )
        )
    if project_tip_nodes:
        clade_definitions.append(
            (
                "project_linked_all",
                "every representative linked to any project record, wherever GTDB places it",
                project_tip_nodes,
            )
        )

    clade_rows: list[list[str]] = []
    clades: dict[str, dict[str, object]] = {}
    for name, description, defining in clade_definitions:
        mrca_node = tree.mrca(defining)
        under = tree.descendant_tips(mrca_node)
        defining_set = set(defining)
        intruders = [node for node in under if node not in defining_set]
        intruder_genera = Counter(
            genus_of(index.representative_taxonomy.get(strip_gtdb_prefix(tree.label[node]), ""))
            for node in intruders
        )
        # The sister group and the support on the path matter for rooting: a clade that is
        # monophyletic at 1% support is not a boundary a root can be argued from.
        parent_node = tree.parent[mrca_node]
        sister_tips: list[int] = []
        if parent_node != -1:
            for sibling in tree.children[parent_node]:
                if sibling != mrca_node:
                    sister_tips.extend(tree.descendant_tips(sibling))
        sister_genera = Counter(
            genus_of(index.representative_taxonomy.get(strip_gtdb_prefix(tree.label[node]), ""))
            for node in sister_tips
        )
        path_labels = []
        walker = mrca_node
        for _ in range(6):
            if tree.parent[walker] == -1:
                break
            walker = tree.parent[walker]
            path_labels.append(tree.label[walker] or "-")
        clades[name] = {
            "description": description,
            "defining_tips": len(defining),
            "mrca_node": mrca_node,
            "mrca_label": tree.label[mrca_node],
            "tips_under_mrca": len(under),
            "intruder_tips": len(intruders),
            "monophyletic": not intruders,
            "top_intruder_genera": dict(intruder_genera.most_common(10)),
            "parent_label": tree.label[parent_node] if parent_node != -1 else "",
            "sister_tips": len(sister_tips),
            "sister_top_genera": dict(sister_genera.most_common(8)),
            "path_to_root_labels": path_labels,
            "tip_set": set(under),
            "defining_set": defining_set,
            "mrca_of_defining": mrca_node,
        }
        clade_rows.append(
            [
                name,
                description,
                len(defining),
                tree.label[mrca_node],
                len(under),
                len(intruders),
                "yes" if not intruders else "no",
                "; ".join(f"{genus}={count}" for genus, count in intruder_genera.most_common(10)),
                tree.label[parent_node] if parent_node != -1 else "",
                len(sister_tips),
                "; ".join(f"{genus}={count}" for genus, count in sister_genera.most_common(8)),
                " > ".join(path_labels),
            ]
        )
    write_tsv(
        out / "tables" / "clade-summary.tsv",
        [
            "clade",
            "definition",
            "defining_tips",
            "mrca_decorated_label",
            "tips_under_mrca",
            "intruder_tips",
            "monophyletic_in_r232_reference_tree",
            "top_intruder_genera",
            "parent_node_decorated_label",
            "sister_group_tips",
            "sister_group_top_genera",
            "decorated_labels_toward_root",
        ],
        clade_rows,
    )

    placement: dict[str, object] = {
        name: {key: value for key, value in info.items() if key not in ("tip_set", "defining_set")}
        for name, info in clades.items()
    }

    # (display name, manifest accession, evidence level, representative accession)
    candidates: list[tuple[str, str, str, str]] = [
        (
            assignment.record.extra.get("taxon_name", assignment.record.record_id),
            assignment.record.accession_text,
            assignment.evidence_level,
            assignment.representative,
        )
        for assignment in assignments
        if assignment.record.role == args.placement_role and assignment.representative in tip_accessions
    ]
    # Taxa with no verified project genome still have a GTDB representative, and their
    # placement is still evidence. Kosakonia cowanii is the case this exists for.
    species_to_rep = {
        species: representative for representative, species in index.representative_species.items()
    }
    for wanted in args.placement_species:
        representative = species_to_rep.get(wanted)
        if representative is None or representative not in tip_accessions:
            log(f"placement probe {wanted!r}: no representative tip in this release/tree")
            continue
        candidates.append((wanted, "", "gtdb_species_probe_no_project_genome", representative))

    # One row per candidate genome, collapsing the manifest's several accessions per taxon
    # only in the reporting order, never in the evidence.
    placement_rows: list[list[str]] = []
    by_clade_rows: list[list[str]] = []
    seen_reps: set[str] = set()
    for display_name, manifest_accession, level, representative in candidates:
        node = tip_node_of_accession[representative]
        primary = clades.get("project_linked_all")
        reference_mrca = primary["mrca_node"] if primary else 0
        join_point = tree.mrca([node, reference_mrca])
        nearest_node = min(project_tip_nodes, key=lambda other: patristic(node, other))
        placement_rows.append(
            [
                display_name,
                manifest_accession,
                level,
                representative,
                species_of(index.representative_species.get(representative, "")),
                genus_of(index.representative_taxonomy.get(representative, "")),
                *[
                    "inside" if node in clades[name]["tip_set"] else "outside"
                    for name, _, _ in clade_definitions
                ],
                f"{patristic(node, reference_mrca):.5f}",
                f"{patristic(node, nearest_node):.5f}",
                species_of(index.representative_species.get(strip_gtdb_prefix(tree.label[nearest_node]), "")),
                len(tree.ancestors(node)) - len(tree.ancestors(join_point)),
                len(tree.descendant_tips(join_point)),
                tree.label[join_point],
                " > ".join(labels_on_path(node, join_point)[-6:]),
            ]
        )
        if representative in seen_reps:
            continue
        seen_reps.add(representative)
        for name, _, _ in clade_definitions:
            info = clades[name]
            inside = node in info["tip_set"]
            by_clade_rows.append(
                [
                    display_name,
                    representative,
                    species_of(index.representative_species.get(representative, "")),
                    name,
                    "inside" if inside else "outside",
                    "no" if inside else "yes",
                    f"{patristic(node, info['mrca_node']):.5f}",
                    info["mrca_label"],
                ]
            )

    placement_rows.sort(key=lambda row: (str(row[0]), str(row[1])))
    write_tsv(
        out / "tables" / "outgroup-placement.tsv",
        [
            "candidate_taxon",
            "manifest_accession",
            "evidence_level",
            "gtdb_representative_tip",
            "gtdb_species",
            "gtdb_genus",
            *[f"relative_to_{name}" for name, _, _ in clade_definitions],
            "patristic_to_project_linked_all_mrca",
            "patristic_to_nearest_project_tip",
            "nearest_project_tip_species",
            "edges_to_join_point",
            "tips_under_join_point",
            "join_point_decorated_label",
            "decorated_labels_toward_join_point",
        ],
        placement_rows,
    )
    write_tsv(
        out / "tables" / "outgroup-placement-by-clade.tsv",
        [
            "candidate_taxon",
            "gtdb_representative_tip",
            "gtdb_species",
            "clade",
            "position",
            "usable_as_outgroup_for_this_clade",
            "patristic_to_clade_mrca",
            "clade_mrca_decorated_label",
        ],
        by_clade_rows,
    )

    # --- per-record placement of the focus tips ----------------------------------------
    focus_rows = []
    for accession in sorted(focus_reps):
        node = tip_node_of_accession[accession]
        linked = records_by_rep.get(accession, [])
        role_counts = Counter(a.record.role for a in linked)
        focus_rows.append(
            [
                species_of(index.representative_species.get(accession, "")),
                genus_of(index.representative_taxonomy.get(accession, "")),
                accession,
                index.cluster_size.get(accession, ""),
                sum(role_counts.values()),
                *[role_counts.get(role, 0) for role in roles],
            ]
        )
    focus_rows.sort(key=lambda row: (str(row[1]), str(row[0])))
    write_tsv(
        out / "tables" / "focus-genus-clusters.tsv",
        [
            "gtdb_species",
            "gtdb_genus",
            "gtdb_representative",
            "gtdb_cluster_size",
            "project_records",
            *[f"records_{role}" for role in roles],
        ],
        focus_rows,
    )

    # --- provenance and validation ------------------------------------------------------
    provenance["evidence_level_counts"] = {
        role: {level: counts[level] for level in EVIDENCE_LEVELS} for role, counts in level_counts.items()
    }
    provenance["subtree_tip_counts"] = {name: len(keep) for name, keep in subtrees.items()}
    provenance["verification"] = verification
    provenance["placement"] = placement
    provenance["elapsed_seconds"] = round(time.time() - started, 2)
    provenance["all_verifications_passed"] = all(
        bool(result["passed"]) for result in verification.values()
    )
    (out / "validation.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "logs" / "gtdb_context.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    if not provenance["all_verifications_passed"]:
        log("VERIFICATION FAILED; see validation.json")
        return 1
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
