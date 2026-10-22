#!/usr/bin/env python3
"""Audit archived ANI / ANIm / Galah tables for provenance, schema and cohort membership.

Task 12 BUILD FIRST interface (docs/analysis-execution-plan.md).

Contract
--------
* Read a declarative config naming each archived table, the schema and units that the
  declaration claims, and the method evidence that is actually recoverable from the archive.
* Verify each declaration against the file. A schema or unit that cannot be confirmed, or
  that matches more than one registered schema, is REJECTED: the source is inventoried and
  counted but contributes no normalized view.
* Emit a lossless per-source row ledger (every source line is retained verbatim, with a
  disposition and a counted reason).
* Emit explicitly separated normalized views, per source and never merged across sources or
  methods: directed rows, unique unordered pairs, and self comparisons.
* Retain both aligned fractions wherever the source has them; never synthesise one.
* Reject duplicate source keys: an exactly repeated directed key keeps one copy and counts
  the rest; a repeated key with conflicting values excludes every copy.
* Flag self comparisons and accessions outside the Task 01 study manifest.
* Give a counted reason for every excluded or unresolved row.
* Never blend methods. Each normalized file carries exactly one source and one declared
  method; the method index records which is which.

This script measures and records. It makes no taxonomic, threshold or naming decision.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import sys
import time
from collections import Counter, defaultdict

SCRIPT_VERSION = "1.0.0"

# --------------------------------------------------------------------------------------
# Schema registry
# --------------------------------------------------------------------------------------
# role indices: ref, query, ani, af_ref, af_query  (missing roles are absent from the schema)

SCHEMAS = {
    "skani_af7": {
        "header": ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query",
                   "Ref_name", "Query_name"],
        "n_fields": 7,
        "roles": {"ref": 0, "query": 1, "ani": 2, "af_ref": 3, "af_query": 4},
    },
    "skani_af10_species": {
        "header": ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query",
                   "Ref_name", "Query_name", "Ref_Sp", "Query_Sp", "Species_Comparison"],
        "n_fields": 10,
        "roles": {"ref": 0, "query": 1, "ani": 2, "af_ref": 3, "af_query": 4},
        "label_fields": {"ref_label": 7, "query_label": 8, "comparison": 9},
    },
    "skani_af9_newspecies": {
        "header": ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query",
                   "Ref_name", "Query_name", "Ref_NewSpecies", "Query_NewSpecies"],
        "n_fields": 9,
        "roles": {"ref": 0, "query": 1, "ani": 2, "af_ref": 3, "af_query": 4},
        "label_fields": {"ref_label": 7, "query_label": 8},
    },
    "pair_ani3_header": {
        "header": ["Ref_file", "Query_file", "ANI"],
        "n_fields": 3,
        "roles": {"ref": 0, "query": 1, "ani": 2},
    },
    "pair_ani3_noheader": {
        "header": None,
        "n_fields": 3,
        "roles": {"ref": 0, "query": 1, "ani": 2},
    },
    "pair_ani5_species": {
        "header": ["Ref_file", "Query_file", "ANI", "Ref_Sp", "Query_Sp"],
        "n_fields": 5,
        "roles": {"ref": 0, "query": 1, "ani": 2},
        "label_fields": {"ref_label": 3, "query_label": 4},
    },
    "pair_ani6_species_cmp": {
        "header": ["Ref_file", "Query_file", "ANI", "Ref_Sp", "Query_Sp", "Species_Comparison"],
        "n_fields": 6,
        "roles": {"ref": 0, "query": 1, "ani": 2},
        "label_fields": {"ref_label": 3, "query_label": 4, "comparison": 5},
    },
    "pair_ani6_cmp": {
        "header": ["Ref_file", "Query_file", "ANI", "Ref_Sp", "Query_Sp", "Comparison"],
        "n_fields": 6,
        "roles": {"ref": 0, "query": 1, "ani": 2},
        "label_fields": {"ref_label": 3, "query_label": 4, "comparison": 5},
    },
    "pair_ani_rounded6": {
        "header": ["Ref_file", "Query_file", "ANI_Rounded", "Ref", "Query", "Species_Comparison"],
        "n_fields": 6,
        "roles": {"ref": 0, "query": 1, "ani": 2},
        "label_fields": {"ref_label": 3, "query_label": 4, "comparison": 5},
        "ani_is_rounded": True,
    },
    "anim_pair10_trailing": {
        "header": None,
        "n_fields": 10,
        "roles": {"ref": 0, "query": 1, "ani": 2},
        "require_empty": [3, 4, 5, 6, 7, 8, 9],
    },
    "pairlist2": {
        "header": None,
        "n_fields": 2,
        "roles": {"ref": 0, "query": 1},
        "field_pattern": {0: r"\.fna$", 1: r"\.fna$"},
        "structural": "not_cluster_definition",
    },
    "galah_pairmap2": {
        "header": None,
        "n_fields": 2,
        "roles": {"ref": 0, "query": 1},
        "field_pattern": {1: r"\.fna$"},
        "structural": "cluster_definition",
    },
    "galah_assign2": {
        "header": None,
        "n_fields": 2,
        "roles": {"ref": 0},
        "field_pattern": {1: r"^Cluster_\d+$"},
    },
    "label2_header": {"header": ["Genome", "Species"], "n_fields": 2, "roles": {"ref": 0}},
    "label2_noheader": {"header": None, "n_fields": 2, "roles": {"ref": 0},
                        "field_pattern": {0: r"_parsed"},
                        "field_pattern_not": {1: r"\.fna$|^Cluster_\d+$"},
                        "structural": "not_cluster_definition"},
    "label3_header": {"header": ["genome", "revised_species", "revised_subspecies"],
                      "n_fields": 3, "roles": {"ref": 0}},
    "kraken_summary5": {"header": ["filename", "NCBI", "Kraken2_parsed", "Kraken2", "ST"],
                        "n_fields": 5, "roles": {"ref": 0}},
    "genomelist1": {"header": None, "n_fields": 1, "roles": {"ref": 0}},
    "freetext": {"header": None, "n_fields": None, "roles": {}},
}

ANI_BEARING = {"ani_table"}
KEYED_KINDS = {"ani_table", "pair_list", "cluster_table", "label_table", "genome_list"}

ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)")


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------

def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def genome_key(raw):
    """Archive genome key: basename with a single trailing .fna removed. Path prefixes differ
    between archived copies of the same table and are not part of the identity."""
    key = raw.strip().rstrip("\r").split("/")[-1]
    if key.endswith(".fna"):
        key = key[:-4]
    return key


def accession_of(key):
    m = ACCESSION_RE.match(key)
    return m.group(1) if m else ""


def open_text(path):
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def gz(path):
    return gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=6)


class Tsv:
    def __init__(self, path, columns, compress=False):
        self.path = path
        self.columns = columns
        self.fh = gz(path) if compress else open(path, "w", encoding="utf-8", newline="")
        self.w = csv.writer(self.fh, delimiter="\t", lineterminator="\n",
                            quoting=csv.QUOTE_MINIMAL)
        self.w.writerow(columns)
        self.n = 0

    def row(self, values):
        self.w.writerow(values)
        self.n += 1

    def close(self):
        self.fh.close()


def esc(text):
    """Keep a source line verbatim but safe inside one TSV field."""
    return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")


# --------------------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------------------

def load_manifest(path, exclusions_path=None):
    study_acc, legacy_to_acc = set(), {}
    with open_text(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if "accession_version" not in (rdr.fieldnames or []):
            raise SystemExit(f"manifest {path} has no accession_version column")
        for row in rdr:
            acc = (row.get("accession_version") or "").strip()
            if acc:
                study_acc.add(acc)
            gid = (row.get("legacy_metadata_genome_id") or "").strip()
            if gid and acc:
                legacy_to_acc[gid] = acc
    excluded_acc = set()
    if exclusions_path and os.path.exists(exclusions_path):
        with open_text(exclusions_path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                eid = (row.get("entity_id") or "").strip()
                if eid:
                    excluded_acc.add(eid)
    return study_acc, legacy_to_acc, excluded_acc


# --------------------------------------------------------------------------------------
# schema and unit verification
# --------------------------------------------------------------------------------------

def peek(path, has_header, max_rows=4000):
    """Return (header_fields, sample_rows, field_count_histogram, total_line_count)."""
    hist = Counter()
    sample = []
    header = None
    total = 0
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if line.endswith("\r"):
                line = line[:-1]
            if i == 0 and has_header:
                header = line.split("\t")
                continue
            total += 1
            fields = line.split("\t")
            hist[len(fields)] += 1
            if len(sample) < max_rows:
                sample.append(fields)
    return header, sample, hist, total


def looks_like_cluster_definition(rows):
    """A galah --output-cluster-definition table lists every genome exactly once in the
    member column, and every representative is also a member of its own cluster. A
    below-threshold pair list satisfies neither. Two-column tables are otherwise
    indistinguishable from one row's shape, so this content check is what separates them."""
    if not rows:
        return False
    a = [r[0].strip() for r in rows]
    b = [r[1].strip() for r in rows]
    return len(set(b)) == len(b) and set(a) <= set(b)


def schema_candidates(header, sample, hist, sample_complete=True):
    """Every registered schema consistent with what the file actually contains.

    Structural checks (see looks_like_cluster_definition) need the whole file; when the
    sample was truncated, schemas that depend on one are not offered as candidates."""
    out = []
    if not hist:
        return out
    modal_width = hist.most_common(1)[0][0]
    modal_sample = [r for r in sample if len(r) == modal_width]
    for name, spec in SCHEMAS.items():
        if spec["n_fields"] is None:
            continue
        if spec["header"] is not None:
            if header is None:
                continue
            if [h.strip() for h in header] != spec["header"]:
                continue
        else:
            if header is not None:
                continue
        if spec["n_fields"] != modal_width:
            continue
        ok = True
        for idx in spec.get("require_empty", []):
            if any(r[idx].strip() != "" for r in modal_sample):
                ok = False
                break
        if ok:
            for idx, pat in spec.get("field_pattern", {}).items():
                rx = re.compile(pat)
                if not all(rx.search(r[idx].strip()) for r in modal_sample):
                    ok = False
                    break
        if ok:
            for idx, pat in spec.get("field_pattern_not", {}).items():
                rx = re.compile(pat)
                if any(rx.search(r[idx].strip()) for r in modal_sample):
                    ok = False
                    break
        if ok and "structural" in spec:
            if not sample_complete:
                ok = False
            else:
                is_cluster = looks_like_cluster_definition(modal_sample)
                want = spec["structural"] == "cluster_definition"
                ok = (is_cluster == want)
        if ok:
            out.append(name)
    return out


def numeric_summary(sample, idx):
    vals, bad = [], 0
    for r in sample:
        if len(r) <= idx:
            bad += 1
            continue
        t = r[idx].strip()
        if t == "":
            bad += 1
            continue
        try:
            vals.append(float(t))
        except ValueError:
            bad += 1
    return vals, bad


def classify_unit(vals):
    """percent | fraction | ambiguous_fraction_or_percent | out_of_range | no_values"""
    if not vals:
        return "no_values"
    lo, hi = min(vals), max(vals)
    if lo < 0 or hi > 100:
        return "out_of_range"
    if hi > 1.0:
        return "percent"
    return "ambiguous_fraction_or_percent"


# --------------------------------------------------------------------------------------
# per-source audit
# --------------------------------------------------------------------------------------

LEDGER_COLS = [
    "source_id", "source_line_no", "disposition", "reasons",
    "ref_key", "query_key", "ref_accession", "query_accession",
    "ani", "af_ref", "af_query", "raw_line",
]
DIRECTED_COLS = [
    "source_id", "method_declared", "source_line_no",
    "ref_key", "query_key", "ref_accession", "query_accession",
    "ani", "af_ref", "af_query", "cohort_status", "flags",
]
PAIR_COLS = [
    "source_id", "method_declared", "genome_a", "genome_b", "accession_a", "accession_b",
    "n_directed_rows", "reciprocal", "ani_ab", "ani_ba", "ani_abs_delta",
    "af_ref_ab", "af_query_ab", "af_ref_ba", "af_query_ba",
    "cohort_status", "flags",
]
SELF_COLS = ["source_id", "method_declared", "source_line_no", "genome_key", "accession",
             "ani", "af_ref", "af_query", "cohort_status"]


def audit_source(src, out_dir, study_acc, legacy_to_acc, excluded_acc, write_ledger=True):
    path = src["path"]
    kind = src["kind"]
    sid = src["id"]
    counts = Counter()
    inv = {
        "source_id": sid,
        "kind": kind,
        "path": path,
        "exists": os.path.exists(path),
        "bytes": os.path.getsize(path) if os.path.exists(path) else 0,
        "sha256": "",
        "mtime_utc": "",
        "declared_schema": src.get("schema", ""),
        "schema_status": "",
        "schema_candidates": "",
        "declared_method": src.get("method_declared", "unknown"),
        "method_evidence": src.get("method_evidence", ""),
        "software_version_in_file": src.get("software_version_in_file", "unknown"),
        "scope_declared": src.get("scope_declared", "unknown"),
        "filter_declared": src.get("filter_declared", "unknown"),
        "expected_direction": src.get("expected_direction", "unknown"),
        "observed_direction": "",
        "declared_ani_unit": src.get("ani_unit_declared", "absent"),
        "observed_ani_unit": "",
        "declared_af_unit": src.get("af_unit_declared", "absent"),
        "observed_af_unit": "",
        "has_header": src.get("has_header", False),
        "field_count_histogram": "",
        "source_rows": 0,
        "genomes": 0,
        "genomes_resolved_to_accession": 0,
        "self_rows": 0,
        "directed_rows_kept": 0,
        "unordered_pairs": 0,
        "reciprocal_pairs": 0,
        "single_direction_pairs": 0,
        "duplicate_key_rows_identical": 0,
        "duplicate_key_rows_conflicting": 0,
        "rows_missing_ani": 0,
        "rows_unparseable": 0,
        "rows_missing_af": 0,
        "rows_off_modal_field_count": 0,
        "ani_min": "", "ani_max": "",
        "af_min": "",
        "accessions_in_study_manifest": 0,
        "accessions_off_cohort": 0,
        "accessions_in_exclusions": 0,
        "off_cohort_accession_examples": "",
        "rows_both_in_study": 0,
        "rows_any_off_cohort": 0,
        "status": "",
        "reject_reasons": "",
        "ledger_path": "",
        "directed_path": "",
        "pairs_path": "",
        "self_path": "",
    }
    if not inv["exists"]:
        inv["status"] = "rejected"
        inv["reject_reasons"] = "file_missing"
        counts["file_missing"] = 1
        return inv, counts

    inv["sha256"] = sha256_file(path)
    inv["mtime_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(os.path.getmtime(path)))

    if kind == "provenance_log" or src.get("schema") == "freetext":
        with open_text(path) as fh:
            text = fh.read()
        inv["source_rows"] = text.count("\n")
        inv["status"] = "inventoried_only"
        inv["schema_status"] = "free_text_not_parsed"
        counts["free_text_lines"] = inv["source_rows"]
        return inv, counts

    header, sample, hist, total = peek(path, src.get("has_header", False))
    inv["source_rows"] = total
    inv["field_count_histogram"] = ";".join(f"{k}:{v}" for k, v in sorted(hist.items()))

    sample_complete = total <= len(sample) or total == 0
    cands = schema_candidates(header, sample, hist, sample_complete=sample_complete)
    inv["schema_candidates"] = ",".join(cands) if cands else "none"
    declared = src.get("schema", "")
    reject = []
    if declared not in SCHEMAS:
        reject.append("declared_schema_not_registered")
        inv["schema_status"] = "unregistered"
    elif total == 0:
        # header-only file: the header alone can still confirm a headered schema
        if SCHEMAS[declared]["header"] is not None and header is not None and \
                [h.strip() for h in header] == SCHEMAS[declared]["header"]:
            inv["schema_status"] = "confirmed_by_header_no_data_rows"
        else:
            inv["schema_status"] = "no_data_rows_unconfirmed"
            reject.append("no_data_rows")
    elif declared not in cands:
        inv["schema_status"] = "declared_schema_contradicted_by_file"
        reject.append("schema_mismatch")
    elif len(cands) > 1:
        inv["schema_status"] = "ambiguous_multiple_schemas_match"
        reject.append("schema_ambiguous")
    elif SCHEMAS[declared]["header"] is not None:
        inv["schema_status"] = "confirmed_by_header"
    else:
        inv["schema_status"] = "confirmed_by_shape_and_constraints"

    spec = SCHEMAS.get(declared, {"roles": {}, "n_fields": None})
    roles = spec.get("roles", {})

    # unit verification
    modal_width = hist.most_common(1)[0][0] if hist else 0
    modal_sample = [r for r in sample if len(r) == modal_width]
    inv["rows_off_modal_field_count"] = total - hist.get(modal_width, 0)
    if "ani" in roles and total:
        vals, bad = numeric_summary(modal_sample, roles["ani"])
        inv["observed_ani_unit"] = classify_unit(vals)
        declared_unit = src.get("ani_unit_declared", "absent")
        base_declared = "percent" if declared_unit.startswith("percent") else declared_unit
        if inv["observed_ani_unit"] in ("ambiguous_fraction_or_percent", "out_of_range", "no_values"):
            reject.append(f"ani_unit_{inv['observed_ani_unit']}")
        elif inv["observed_ani_unit"] != base_declared:
            reject.append("ani_unit_mismatch")
    else:
        inv["observed_ani_unit"] = "absent"
    if "af_ref" in roles and total:
        vals, _ = numeric_summary(modal_sample, roles["af_ref"])
        inv["observed_af_unit"] = classify_unit(vals)
        if inv["observed_af_unit"] in ("ambiguous_fraction_or_percent", "out_of_range", "no_values"):
            reject.append(f"af_unit_{inv['observed_af_unit']}")
    else:
        inv["observed_af_unit"] = "absent"

    inv["reject_reasons"] = ";".join(sorted(set(reject)))
    if reject:
        inv["status"] = "rejected"
        counts.update({f"rejected:{r}": 1 for r in set(reject)})
        counts["source_rows_not_normalized"] = total
        return inv, counts
    inv["status"] = "accepted"

    if kind not in KEYED_KINDS or "ref" not in roles:
        counts["source_rows_not_normalized"] = total
        return inv, counts

    # ---------------- streaming pass ----------------
    norm = os.path.join(out_dir, "normalized")
    led = os.path.join(out_dir, "ledger")
    os.makedirs(norm, exist_ok=True)
    os.makedirs(led, exist_ok=True)

    ledger = None
    if write_ledger:
        inv["ledger_path"] = os.path.join(led, f"{sid}.ledger.tsv.gz")
        ledger = Tsv(inv["ledger_path"], LEDGER_COLS, compress=True)

    has_ani = "ani" in roles
    has_af = "af_ref" in roles and "af_query" in roles
    method = src.get("method_declared", "unknown")

    directed = None
    if has_ani:
        inv["directed_path"] = os.path.join(norm, f"{sid}.directed.tsv.gz")
        directed = Tsv(inv["directed_path"], DIRECTED_COLS, compress=True)
        inv["self_path"] = os.path.join(norm, f"{sid}.self.tsv")
        selfw = Tsv(inv["self_path"], SELF_COLS)
    else:
        selfw = None

    seen = {}            # directed key -> (value tuple, first line no)
    conflicting = set()
    genomes = set()
    ani_min = ani_max = None
    af_min = None
    pairs = defaultdict(list)  # frozenset -> list of directed records

    with open_text(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if raw.endswith("\r"):
                raw = raw[:-1]
            if lineno == 1 and src.get("has_header", False):
                if ledger:
                    ledger.row([sid, lineno, "header", "header_row", "", "", "", "",
                                "", "", "", esc(raw)])
                counts["header_rows"] += 1
                continue
            f = raw.split("\t")
            reasons = []
            if spec["n_fields"] is not None and len(f) != spec["n_fields"]:
                counts["excluded:field_count_mismatch"] += 1
                inv["rows_unparseable"] += 1
                if ledger:
                    ledger.row([sid, lineno, "excluded", "field_count_mismatch",
                                "", "", "", "", "", "", "", esc(raw)])
                continue

            rkey = genome_key(f[roles["ref"]])
            qkey = genome_key(f[roles["query"]]) if "query" in roles else ""
            racc = accession_of(rkey)
            qacc = accession_of(qkey) if qkey else ""
            if not racc and rkey in legacy_to_acc:
                racc = legacy_to_acc[rkey]
            if not qacc and qkey in legacy_to_acc:
                qacc = legacy_to_acc[qkey]
            genomes.add(rkey)
            if qkey:
                genomes.add(qkey)

            ani = af_r = af_q = ""
            if has_ani:
                t = f[roles["ani"]].strip()
                if t == "":
                    inv["rows_missing_ani"] += 1
                    counts["excluded:missing_ani"] += 1
                    if ledger:
                        ledger.row([sid, lineno, "excluded", "missing_ani", rkey, qkey,
                                    racc, qacc, "", "", "", esc(raw)])
                    continue
                try:
                    ani = float(t)
                except ValueError:
                    inv["rows_unparseable"] += 1
                    counts["excluded:ani_not_numeric"] += 1
                    if ledger:
                        ledger.row([sid, lineno, "excluded", "ani_not_numeric", rkey, qkey,
                                    racc, qacc, "", "", "", esc(raw)])
                    continue
                ani_min = ani if ani_min is None else min(ani_min, ani)
                ani_max = ani if ani_max is None else max(ani_max, ani)
            if has_af:
                for idx, slot in ((roles["af_ref"], "r"), (roles["af_query"], "q")):
                    t = f[idx].strip()
                    v = ""
                    if t != "":
                        try:
                            v = float(t)
                        except ValueError:
                            v = ""
                    if slot == "r":
                        af_r = v
                    else:
                        af_q = v
                if af_r == "" or af_q == "":
                    inv["rows_missing_af"] += 1
                    reasons.append("missing_aligned_fraction")
                else:
                    m = min(af_r, af_q)
                    af_min = m if af_min is None else min(af_min, m)

            # cohort membership
            rin = racc in study_acc if racc else False
            qin = qacc in study_acc if qacc else (not qkey)
            if qkey:
                cohort = "both_in_study" if (rin and qin) else (
                    "both_off_cohort" if (not rin and not qin) else "one_off_cohort")
            else:
                cohort = "in_study" if rin else "off_cohort"
            if not rin and racc:
                reasons.append("ref_off_study_manifest")
            if qkey and qacc and not qin:
                reasons.append("query_off_study_manifest")
            if racc == "":
                reasons.append("ref_accession_unresolved")
            if qkey and qacc == "":
                reasons.append("query_accession_unresolved")

            if not has_ani:
                if ledger:
                    ledger.row([sid, lineno, "keyed_no_ani", ";".join(reasons) or "ok",
                                rkey, qkey, racc, qacc, "", "", "", esc(raw)])
                counts["keyed_rows_without_ani"] += 1
                continue

            if qkey and rkey == qkey:
                inv["self_rows"] += 1
                counts["self_comparison_rows"] += 1
                if selfw:
                    selfw.row([sid, method, lineno, rkey, racc, ani, af_r, af_q, cohort])
                if ledger:
                    ledger.row([sid, lineno, "self_pair",
                                ";".join(["self_comparison"] + reasons), rkey, qkey, racc, qacc,
                                ani, af_r, af_q, esc(raw)])
                continue

            dkey = (rkey, qkey)
            vals = (ani, af_r, af_q)
            if dkey in seen:
                prev_vals, prev_line = seen[dkey]
                if prev_vals == vals:
                    inv["duplicate_key_rows_identical"] += 1
                    counts["excluded:duplicate_source_key_identical"] += 1
                    if ledger:
                        ledger.row([sid, lineno, "excluded",
                                    f"duplicate_source_key_identical;first_at_line_{prev_line}",
                                    rkey, qkey, racc, qacc, ani, af_r, af_q, esc(raw)])
                else:
                    conflicting.add(dkey)
                    inv["duplicate_key_rows_conflicting"] += 1
                    counts["excluded:duplicate_source_key_conflicting"] += 1
                    if ledger:
                        ledger.row([sid, lineno, "excluded",
                                    f"duplicate_source_key_conflicting;first_at_line_{prev_line}",
                                    rkey, qkey, racc, qacc, ani, af_r, af_q, esc(raw)])
                continue
            seen[dkey] = (vals, lineno)

            directed.row([sid, method, lineno, rkey, qkey, racc, qacc, ani, af_r, af_q,
                          cohort, ";".join(reasons)])
            pairs[frozenset((rkey, qkey))].append((rkey, qkey, ani, af_r, af_q, cohort, reasons))
            inv["directed_rows_kept"] += 1
            counts["normalized_directed_rows"] += 1
            if cohort == "both_in_study":
                inv["rows_both_in_study"] += 1
            else:
                inv["rows_any_off_cohort"] += 1
            if ledger:
                ledger.row([sid, lineno, "normalized_directed", ";".join(reasons) or "ok",
                            rkey, qkey, racc, qacc, ani, af_r, af_q, esc(raw)])

    if ledger:
        ledger.close()

    # retract directed rows whose key later turned out to be conflicting
    if conflicting and has_ani:
        directed.close()
        tmp = inv["directed_path"] + ".tmp"
        removed = 0
        with gzip.open(inv["directed_path"], "rt", encoding="utf-8", newline="") as src_fh, \
                gz(tmp) as dst_fh:
            rdr = csv.reader(src_fh, delimiter="\t")
            wtr = csv.writer(dst_fh, delimiter="\t", lineterminator="\n")
            wtr.writerow(next(rdr))
            for row in rdr:
                if (row[3], row[4]) in conflicting:
                    removed += 1
                    continue
                wtr.writerow(row)
        os.replace(tmp, inv["directed_path"])
        inv["directed_rows_kept"] -= removed
        inv["duplicate_key_rows_conflicting"] += removed
        counts["excluded:duplicate_source_key_conflicting"] += removed
        counts["normalized_directed_rows"] -= removed
        for k in list(pairs):
            pairs[k] = [r for r in pairs[k] if (r[0], r[1]) not in conflicting]
            if not pairs[k]:
                del pairs[k]
    elif has_ani:
        directed.close()

    if has_ani:
        selfw.close()
        inv["pairs_path"] = os.path.join(norm, f"{sid}.pairs.tsv.gz")
        pw = Tsv(inv["pairs_path"], PAIR_COLS, compress=True)
        recip = single = 0
        for key, recs in pairs.items():
            a, b = sorted(key)
            fwd = next((r for r in recs if r[0] == a), None)
            rev = next((r for r in recs if r[0] == b), None)
            ani_ab = fwd[2] if fwd else ""
            ani_ba = rev[2] if rev else ""
            delta = abs(ani_ab - ani_ba) if (fwd and rev) else ""
            reciprocal = bool(fwd and rev)
            if reciprocal:
                recip += 1
            else:
                single += 1
            cohort = recs[0][5]
            flags = sorted({x for r in recs for x in r[6]})
            if reciprocal and delta != "" and delta > 0:
                flags.append("reciprocal_ani_disagreement")
            pw.row([sid, method, a, b, accession_of(a) or legacy_to_acc.get(a, ""),
                    accession_of(b) or legacy_to_acc.get(b, ""),
                    len(recs), reciprocal, ani_ab, ani_ba, delta,
                    fwd[3] if fwd else "", fwd[4] if fwd else "",
                    rev[3] if rev else "", rev[4] if rev else "",
                    cohort, ";".join(flags)])
        pw.close()
        inv["unordered_pairs"] = len(pairs)
        inv["reciprocal_pairs"] = recip
        inv["single_direction_pairs"] = single
        counts["unordered_pairs"] = len(pairs)
        if recip and not single:
            inv["observed_direction"] = "directed_both"
        elif single and not recip:
            inv["observed_direction"] = "undirected_single"
        elif recip and single:
            inv["observed_direction"] = "mixed"
        else:
            inv["observed_direction"] = "none"
        inv["ani_min"] = ani_min if ani_min is not None else ""
        inv["ani_max"] = ani_max if ani_max is not None else ""
        inv["af_min"] = af_min if af_min is not None else ""

    accs = {accession_of(g) or legacy_to_acc.get(g, "") for g in genomes}
    accs.discard("")
    inv["genomes"] = len(genomes)
    inv["genomes_resolved_to_accession"] = len(accs)
    inv["accessions_in_study_manifest"] = len(accs & study_acc)
    off = sorted(accs - study_acc)
    inv["accessions_off_cohort"] = len(off)
    inv["accessions_in_exclusions"] = len(accs & excluded_acc)
    inv["off_cohort_accession_examples"] = ",".join(off[:10])
    counts["genomes"] = len(genomes)
    return inv, counts


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--manifest", required=True,
                    help="Task 01 study manifest (accession_version column)")
    ap.add_argument("--exclusions", default="manifests/exclusions.tsv")
    ap.add_argument("--output", required=True, help="fresh run directory; must not be reused")
    ap.add_argument("--only", action="append", default=None,
                    help="restrict to these source ids (repeatable)")
    ap.add_argument("--no-ledger", action="store_true",
                    help="skip the lossless row ledger (inventory and views only)")
    ap.add_argument("--fail-on-reject", action="store_true",
                    help="exit non-zero if any source is rejected")
    ap.add_argument("--allow-existing", action="store_true")
    args = ap.parse_args(argv)

    started = time.time()
    out = args.output
    if os.path.isdir(out) and os.listdir(out) and not args.allow_existing:
        if any(f.endswith(".tsv") for f in os.listdir(out)):
            raise SystemExit(f"refusing to overwrite a populated run directory: {out}")
    os.makedirs(out, exist_ok=True)

    with open_text(args.config) as fh:
        cfg = json.load(fh)
    sources = cfg["sources"]
    if args.only:
        wanted = set(args.only)
        sources = [s for s in sources if s["id"] in wanted]
        missing = wanted - {s["id"] for s in sources}
        if missing:
            raise SystemExit(f"unknown source ids: {sorted(missing)}")

    ids = [s["id"] for s in sources]
    dupe = [k for k, v in Counter(ids).items() if v > 1]
    if dupe:
        raise SystemExit(f"duplicate source ids in config: {dupe}")

    study_acc, legacy_to_acc, excluded_acc = load_manifest(args.manifest, args.exclusions)

    inv_rows, all_counts = [], []
    for src in sources:
        t0 = time.time()
        inv, counts = audit_source(src, out, study_acc, legacy_to_acc, excluded_acc,
                                   write_ledger=not args.no_ledger)
        inv["elapsed_seconds"] = round(time.time() - t0, 2)
        inv_rows.append(inv)
        all_counts.append((src["id"], counts))
        print(f"[{inv['status']:>16}] {src['id']:<52} rows={inv['source_rows']:>9} "
              f"pairs={inv['unordered_pairs']:>9} {inv['elapsed_seconds']:>7.2f}s",
              file=sys.stderr, flush=True)

    inv_cols = list(inv_rows[0].keys()) if inv_rows else []
    w = Tsv(os.path.join(out, "source-inventory.tsv"), inv_cols)
    for r in inv_rows:
        w.row([r.get(c, "") for c in inv_cols])
    w.close()

    ct = Tsv(os.path.join(out, "count-transitions.tsv"),
             ["source_id", "category", "count"])
    for sid, counts in all_counts:
        for k in sorted(counts):
            ct.row([sid, k, counts[k]])
    ct.close()

    mi = Tsv(os.path.join(out, "method-view-index.tsv"),
             ["source_id", "declared_method", "status", "schema_status", "observed_direction",
              "observed_ani_unit", "directed_view", "pairs_view", "self_view", "ledger"])
    for r in inv_rows:
        mi.row([r["source_id"], r["declared_method"], r["status"], r["schema_status"],
                r["observed_direction"], r["observed_ani_unit"],
                os.path.relpath(r["directed_path"], out) if r["directed_path"] else "",
                os.path.relpath(r["pairs_path"], out) if r["pairs_path"] else "",
                os.path.relpath(r["self_path"], out) if r["self_path"] else "",
                os.path.relpath(r["ledger_path"], out) if r["ledger_path"] else ""])
    mi.close()

    # identical file content groups
    by_hash = defaultdict(list)
    for r in inv_rows:
        if r["sha256"]:
            by_hash[r["sha256"]].append(r["source_id"])
    dg = Tsv(os.path.join(out, "duplicate-content-groups.tsv"),
             ["sha256", "n_sources", "source_ids"])
    for h, sids in sorted(by_hash.items()):
        if len(sids) > 1:
            dg.row([h, len(sids), ",".join(sorted(sids))])
    dg.close()

    rejected = [r for r in inv_rows if r["status"] == "rejected"]
    summary = {
        "script": os.path.basename(__file__),
        "script_version": SCRIPT_VERSION,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
        "elapsed_seconds": round(time.time() - started, 2),
        "config": args.config,
        "config_sha256": sha256_file(args.config),
        "manifest": args.manifest,
        "manifest_sha256": sha256_file(args.manifest),
        "manifest_accessions": len(study_acc),
        "exclusions": args.exclusions if os.path.exists(args.exclusions) else "",
        "exclusions_sha256": sha256_file(args.exclusions) if os.path.exists(args.exclusions) else "",
        "exclusion_entities": len(excluded_acc),
        "sources_total": len(inv_rows),
        "sources_accepted": sum(1 for r in inv_rows if r["status"] == "accepted"),
        "sources_rejected": len(rejected),
        "sources_inventoried_only": sum(1 for r in inv_rows if r["status"] == "inventoried_only"),
        "rejected_sources": {r["source_id"]: r["reject_reasons"] for r in rejected},
        "total_source_rows": sum(r["source_rows"] for r in inv_rows),
        "total_directed_rows_kept": sum(r["directed_rows_kept"] for r in inv_rows),
        "total_unordered_pairs": sum(r["unordered_pairs"] for r in inv_rows),
        "total_self_rows": sum(r["self_rows"] for r in inv_rows),
        "ledger_written": not args.no_ledger,
    }
    with open(os.path.join(out, "audit-summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.fail_on_reject and rejected:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
