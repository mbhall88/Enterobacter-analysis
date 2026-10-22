#!/usr/bin/env python3
"""Snapshot LPSN records for every name relevant to the P6a nomenclature register.

Read-only HTTP GETs of public LPSN pages. Raw HTML is kept unchanged with a
SHA-256 and UTC retrieval time; parsed fields are written to a TSV. The genus
page is used to enumerate all Enterobacter names (including synonyms and names
now placed in other genera), and each species page is followed to its
subspecies children. Extra paths can be added with --extra.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://lpsn.dsmz.de"
UA = "Mozilla/5.0 (P6a nomenclature register; read-only snapshot)"

FIELDS = [
    "Name",
    "Category, proposed as",
    "Category",
    "Type strain",
    "Nomenclatural type",
    "Valid publication",
    "Effective publication",
    "Original publication",
    "IJSEM list",
    "Nomenclatural status",
    "Taxonomic status",
    "Correct name",
    "Parent taxon",
    "Record number",
]


def to_text(raw: str) -> str:
    t = re.sub(r"<script.*?</script>", "", raw, flags=re.S)
    t = re.sub(r"<style.*?</style>", "", t, flags=re.S)
    t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    t = re.sub(r"<br\s*/?>", "\n", t)
    t = re.sub(r"</(p|div|tr|li|h\d|td|th)>", "\n", t)
    t = re.sub(r"<[^>]+>", "", t)
    t = html.unescape(t)
    t = re.sub(r"[ \t]+", " ", t)
    t = re.sub(r"\n\s*\n+", "\n", t)
    return t


def record_block(text: str) -> str:
    i = text.find("Name:")
    j = text.find("Not found what you were looking for")
    return text[i:j] if i >= 0 else ""


def parse_fields(block: str) -> dict[str, str]:
    lines = [ln.strip() for ln in block.split("\n")]
    out: dict[str, str] = {}
    labels = {f + ":" for f in FIELDS} | {
        "Etymology:", "Pronunciation (Latin), gender:", "Gender:", "See detailed strain information at",
        "Conduct genome-based taxonomy at", "16S rRNA gene:", "Analyse", "FASTA", "Risk group:",
        "Synonyms:", "Notes:", "Assigned by:", "Linking:", "Previous versions:", "Citation: If you use this resource in a publication, please cite us.",
        "Children:", "Parent taxon:",
    }
    k = 0
    while k < len(lines):
        ln = lines[k]
        for f in FIELDS:
            lab = f + ":"
            if ln.startswith(lab):
                val = ln[len(lab):].strip()
                k2 = k + 1
                # value continues on following lines until the next label
                while not val and k2 < len(lines) and lines[k2] not in labels:
                    val = lines[k2]
                    k2 += 1
                if f not in out:
                    out[f] = val
                break
        k += 1
    # synonyms table and notes
    syn = []
    if "Synonyms:" in lines:
        s = lines.index("Synonyms:") + 1
        if s + 1 < len(lines) and lines[s] == "Name" and lines[s + 1] == "Kind":
            s += 2
        while s + 1 < len(lines) and lines[s] not in labels:
            syn.append(f"{lines[s]} [{lines[s + 1]}]")
            s += 2
    out["Synonyms"] = " | ".join(syn)
    notes = []
    if "Notes:" in lines:
        s = lines.index("Notes:") + 1
        while s < len(lines) and lines[s] not in labels:
            if not lines[s].startswith("😷"):
                notes.append(lines[s])
            s += 1
    out["Notes"] = " | ".join(notes)
    return out


def fetch(url: str) -> tuple[int, bytes]:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.status, r.read()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--extra", nargs="*", default=[])
    ap.add_argument("--sleep", type=float, default=0.6)
    a = ap.parse_args()
    raw_dir = a.out / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    retrieval = []
    seen: set[str] = set()
    queue = ["/genus/enterobacter"] + list(a.extra)
    records = []
    while queue:
        path = queue.pop(0)
        if path in seen:
            continue
        seen.add(path)
        url = BASE + path
        when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        try:
            status, body = fetch(url)
        except Exception as e:  # keep going; record the failure
            retrieval.append({"path": path, "url": url, "retrieved_utc": when, "http_status": "error",
                              "bytes": 0, "sha256": "", "raw_file": "", "error": repr(e)})
            continue
        fn = raw_dir / (path.strip("/").replace("/", "__") + ".html")
        fn.write_bytes(body)
        retrieval.append({"path": path, "url": url, "retrieved_utc": when, "http_status": status,
                          "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(),
                          "raw_file": str(fn), "error": ""})
        raw = body.decode("utf-8", "replace")
        if path == "/genus/enterobacter":
            for m in sorted(set(re.findall(r'href="(/(?:species|subspecies)/enterobacter-[a-z0-9-]+)"', raw))):
                queue.append(m)
        elif path.startswith("/species/enterobacter-"):
            stem = path.split("/species/")[1]
            for m in sorted(set(re.findall(r'href="(/subspecies/' + re.escape(stem) + r'-[a-z0-9-]+)"', raw))):
                queue.append(m)
        if path.startswith(("/species/", "/subspecies/")):
            f = parse_fields(record_block(to_text(raw)))
            f["path"] = path
            f["url"] = url
            f["retrieved_utc"] = when
            f["raw_sha256"] = retrieval[-1]["sha256"]
            records.append(f)
        time.sleep(a.sleep)
    with open(a.out / "retrieval.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(retrieval[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(retrieval)
    cols = ["path", "url", "retrieved_utc", "raw_sha256"] + FIELDS + ["Synonyms", "Notes"]
    with open(a.out / "lpsn-records.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({c: r.get(c, "") for c in cols})
    json.dump({"n_pages": len(retrieval), "n_records": len(records),
               "errors": [r for r in retrieval if r["http_status"] != 200]},
              open(a.out / "fetch-summary.json", "w"), indent=2)
    print(f"pages={len(retrieval)} records={len(records)} errors={sum(r['http_status'] != 200 for r in retrieval)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
