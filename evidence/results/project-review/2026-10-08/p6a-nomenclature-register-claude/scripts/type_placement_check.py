#!/usr/bin/env python3
"""Metadata competitor screen for the merged treatments (revised 9 Oct 2026).

For every validly published name that LPSN lists under Enterobacter, select a
type assembly (project-verified type panel first; otherwise an NCBI
type-material assembly whose strain identifiers match an LPSN type-strain
designation), then read its GTDB species from the official, MD5-verified
GTDB R232 bacterial taxonomy (bac120_taxonomy_r232.tsv.gz, all genomes, not
only representatives). The screen reports whether any selected type falls in
a GTDB R232 species occupied by a merged group's study members in the frozen
accession overlay. It is a metadata screen: names whose type has no selected
assembly are listed with the reason, never treated as screened.

Revision against the 8 Oct snapshot (kept unchanged in checks/type-placement/):
- NCBI Datasets queries use filters.is_type_material=true and follow
  next_page_token to the end; every raw page is saved (gzip) with its SHA-256,
  total_count and a completeness flag. A failed query aborts the run.
- Every NCBI candidate and the selection decision are written out.
- GTDB species come only from the hashed R232 file. A row that cannot be
  resolved there is labelled as such; if the optional --api-fallback is used,
  the raw card is saved and the row is labelled unversioned API metadata with
  its retrieval date, never as R232.
Read-only public GET requests; no genome download.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UA = "Mozilla/5.0 (P6a nomenclature register; read-only)"
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/"
R232_FILE = "bac120_taxonomy_r232.tsv.gz"


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def slug(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "-", x).strip("-")


class Recorder:
    """Saves every raw response (gzip) and one retrieval row per request."""

    def __init__(self, out: Path):
        self.out = out
        self.rows: list[dict] = []

    def get(self, url: str, kind: str, stem: str, attempts: int = 4) -> dict:
        err = None
        for i in range(attempts):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    raw, status = r.read(), r.status
                break
            except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
                err = e
                time.sleep(2 * (i + 1))
        else:
            self.rows.append({"kind": kind, "url": url, "retrieved_utc": now(), "http_status": "",
                              "error": repr(err), "raw_file": "", "raw_sha256": "", "bytes": ""})
            raise RuntimeError(f"request failed after {attempts} attempts: {url}: {err!r}")
        f = self.out / "raw" / kind / f"{stem}.json.gz"
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(gzip.compress(raw, mtime=0))
        self.rows.append({"kind": kind, "url": url, "retrieved_utc": now(), "http_status": status, "error": "",
                          "raw_file": str(f.relative_to(self.out)), "raw_sha256": hashlib.sha256(raw).hexdigest(),
                          "bytes": len(raw)})
        return json.loads(raw.decode())


def binomial(name: str) -> str:
    n = name.replace('"', "").replace(" corrig.", "")
    m = re.match(r"([A-Z][a-z]+) ([a-z]+)( subsp\. ([a-z]+))?", n)
    return f"{m.group(1)} {m.group(2)}" if m else n


def norm(x: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", x.upper())


def alias_match(aliases: list[str], idents: list[str]) -> bool:
    """Exact designation match after normalisation; whole identifier or any
    whitespace token, allowing a trailing type-strain 'T' superscript."""
    ids = set()
    for i in idents:
        if not i:
            continue
        ids.add(norm(i))
        ids.update(norm(t) for t in i.split())
    ids.discard("")
    for al in aliases:
        a = norm(al)
        if a and (a in ids or a + "T" in ids):
            return True
    return False


def ncbi_type_reports(taxon: str, rec: Recorder, pages: list) -> tuple[list[dict], dict]:
    """All type-material dataset reports for a taxon name, following next_page_token."""
    reports: list[dict] = []
    token, page, total = None, 1, None
    while True:
        url = (NCBI + urllib.parse.quote(taxon) + "/dataset_report?filters.is_type_material=true&page_size=1000"
               + (f"&page_token={urllib.parse.quote(token)}" if token else ""))
        d = rec.get(url, "ncbi", f"{slug(taxon)}.p{page:02d}")
        if page == 1:
            total = d.get("total_count")
        got = d.get("reports", []) or []
        reports += got
        token = d.get("next_page_token")
        pages.append({"query": taxon, "page": page, "url": url, "total_count": total if total is not None else "",
                      "n_reports": len(got), "next_page_token_present": bool(token),
                      "raw_file": rec.rows[-1]["raw_file"], "raw_sha256": rec.rows[-1]["raw_sha256"],
                      "retrieved_utc": rec.rows[-1]["retrieved_utc"]})
        time.sleep(0.4)
        if not token:
            break
        page += 1
    complete = len(reports) == (total or 0)
    return reports, {"total_count": total or 0, "n_reports": len(reports), "pages": page, "complete": complete}


def candidates(taxon: str, aliases: list[str], reports: list[dict]) -> list[dict]:
    out = []
    for r in reports:
        tm = r.get("type_material") or {}
        b = r.get("assembly_info", {}).get("biosample", {}) or {}
        attrs = {a.get("name"): a.get("value", "") for a in b.get("attributes", [])}
        inf = r.get("organism", {}).get("infraspecific_names") or {}
        idents = [inf.get("strain", ""), inf.get("isolate", ""), attrs.get("strain", ""), attrs.get("isolate", ""),
                  attrs.get("culture_collection", ""), attrs.get("type-material", "")]
        idents += [x.get("value", "") for x in b.get("sample_ids", [])]
        idents = [y.strip() for x in idents for y in re.split(r"[;,]| = ", x or "") if y.strip()]
        out.append({"query": taxon, "accession": r["accession"], "type_label": tm.get("type_label", ""),
                    "type_display_text": tm.get("type_display_text", ""),
                    "organism_name": r.get("organism", {}).get("organism_name", ""),
                    "strain": inf.get("strain", "") or inf.get("isolate", ""),
                    "culture_collection": attrs.get("culture_collection", ""),
                    "identifiers_compared": " | ".join(dict.fromkeys(idents)),
                    "alias_match": alias_match(aliases, idents)})
    return out


def rank_key(x: dict):
    return (not x["alias_match"], not x["accession"].startswith("GCA_"), x["type_label"] != "TYPE_MATERIAL", x["accession"])


def load_r232(gz: Path, md5sum: Path) -> tuple[dict, dict, str]:
    want = None
    for line in open(md5sum):
        h, f = line.split()
        if f.lstrip("./") == R232_FILE:
            want = h
    data = gz.read_bytes()
    got = hashlib.md5(data).hexdigest()
    if want is None or got != want:
        raise SystemExit(f"R232 taxonomy MD5 mismatch or missing: expected {want}, got {got}")
    exact, stems = {}, {}
    for line in gzip.decompress(data).decode().splitlines():
        key, tax = line.split("\t")
        acc = key.split("_", 1)[1]
        sp = tax.split(";")[-1]
        exact[acc] = (key, sp)
        stems.setdefault(acc.split("_", 1)[1].split(".")[0], []).append((key, sp))
    return exact, stems, got


def r232_lookup(acc: str, exact: dict, stems: dict) -> tuple[str, str, str]:
    stem, ver = acc.split("_", 1)[1].split(".")
    for a in (acc, ("GCF_" if acc.startswith("GCA_") else "GCA_") + f"{stem}.{ver}"):
        if a in exact:
            key, sp = exact[a]
            return sp, "gtdb_r232_full_taxonomy", key
    if stem in stems:
        hits = stems[stem]
        sp = ";".join(sorted({s for _, s in hits}))
        return sp, "gtdb_r232_full_taxonomy_other_assembly_version", ";".join(k for k, _ in hits)
    return "", "not_in_gtdb_r232", ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lpsn", required=True, type=Path)
    ap.add_argument("--overlay", required=True, type=Path)
    ap.add_argument("--r232", required=True, type=Path, help="official bac120_taxonomy_r232.tsv.gz")
    ap.add_argument("--r232-md5", required=True, type=Path, help="official MD5SUM.txt of the release")
    ap.add_argument("--reps", type=Path, help="local R232 representative taxonomy (cross-check only)")
    ap.add_argument("--panel", required=True, type=Path, help="Task02 approved types.tsv")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--api-fallback", action="store_true",
                    help="query the unversioned GTDB API card for rows absent from R232 (labelled, raw saved)")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    exact, stems, md5 = load_r232(a.r232, a.r232_md5)
    reps = {}
    if a.reps:
        for line in open(a.reps):
            k, tax = line.rstrip("\n").split("\t")
            reps[k.split("_", 1)[1]] = tax.split(";")[-1]
    panel = {}
    for x in csv.DictReader(open(a.panel), delimiter="\t"):
        if x["selected_accession_version"]:
            al = [y.strip() for y in re.split(r"[;|]", x["type_strain_aliases"]) if y.strip()]
            panel[x["taxon_name"]] = (x["selected_accession_version"], al)
    rec = Recorder(a.out)
    pages: list[dict] = []
    cand_rows: list[dict] = []
    rows = []
    for r in csv.DictReader(open(a.lpsn), delimiter="\t"):
        if not r["path"].startswith("/species/enterobacter-"):
            continue
        if "validly published under the ICNP" not in r["Nomenclatural status"]:
            continue
        listed = binomial(r["Name"])
        correct = binomial(r["Correct name"]) if r["Correct name"] else listed
        aliases = [x.strip() for x in r["Type strain"].split(";") if x.strip()]
        acc, sel, strain, label, n_cand, completeness = "", "", "", "", 0, ""
        unmatched = []
        for key, (pacc, palias) in panel.items():
            if aliases and alias_match(aliases, palias):
                acc, sel, strain, label = pacc, "project_verified_type_panel:" + key, "; ".join(palias[:3]), "project_panel"
                break
        if not acc:
            hits, comp = [], []
            for q in dict.fromkeys([listed, correct]):
                reps_q, info = ncbi_type_reports(q, rec, pages)
                comp.append(f"{q}: {info['n_reports']}/{info['total_count']} in {info['pages']} page(s)"
                            + ("" if info["complete"] else " INCOMPLETE"))
                hits += candidates(q, aliases, reps_q)
            completeness = "; ".join(comp)
            hits.sort(key=rank_key)
            n_cand = len({h["accession"] for h in hits})
            chosen = hits[0]["accession"] if hits and hits[0]["alias_match"] else ""
            for i, h in enumerate(hits, 1):
                cand_rows.append({"lpsn_name": r["Name"], "rank": i, "selected": bool(chosen) and i == 1, **h})
            if chosen:
                h = hits[0]
                acc, sel, strain, label = h["accession"], "ncbi_type_material_alias_match", (h["strain"] + " " + h["culture_collection"]).strip(), h["type_label"]
            else:
                sel = "no_alias_matched_type_assembly"
                strain = "; ".join(sorted({h["accession"] + ":" + h["strain"] for h in hits}))
                unmatched = sorted({h["accession"] for h in hits})
        if acc:
            sp, src, key = r232_lookup(acc, exact, stems)
            release = "GTDB R232 (232.0; bac120_taxonomy_r232.tsv.gz MD5 " + md5 + ")" if src.startswith("gtdb_r232") else ""
            if src == "not_in_gtdb_r232" and a.api_fallback:
                d = rec.get(f"https://gtdb-api.ecogenomic.org/genome/{acc}/card", "gtdb-api", acc)
                sp = (d.get("metadataTaxonomy") or {}).get("gtdbSpecies") or ""
                src, key = "gtdb_api_card_unversioned", rec.rows[-1]["raw_file"]
                release = "unversioned GTDB API metadata retrieved " + rec.rows[-1]["retrieved_utc"][:10]
            coverage = "type_assembly_selected"
        else:
            sp, src, key, release = "", "no_type_assembly", "", ""
            coverage = "ncbi_type_candidates_failed_alias_match" if n_cand else "no_ncbi_type_material_candidate"
        cur_genus = correct.split()[0]
        rows.append({"lpsn_name": r["Name"], "lpsn_correct_name": r["Correct name"],
                     "name_category": "current_enterobacter_name" if cur_genus == "Enterobacter" else f"historical_combination_now_{cur_genus}",
                     "type_strain_lpsn": r["Type strain"], "type_assembly": acc, "selection": sel,
                     "matched_identifiers_or_candidates": strain, "type_label": label,
                     "n_ncbi_type_candidates": n_cand, "ncbi_query_completeness": completeness,
                     "coverage_status": coverage, "gtdb_species": sp, "gtdb_release": release, "gtdb_source": src,
                     "gtdb_record": key, "local_r232_reps_species": reps.get(acc, "") if acc else "",
                     "unmatched_candidates_gtdb_r232_species": "; ".join(
                         f"{c}:{r232_lookup(c, exact, stems)[0] or 'not_in_gtdb_r232'}" for c in unmatched)})
    groups = {
        "xiangfangensis_A+B": "Enterobacter xiangfangensis", "roggenkampii": "Enterobacter roggenkampii",
        "chengduensis": "Enterobacter chengduensis", "mori": "Enterobacter mori", "asburiae": "Enterobacter asburiae",
        "cloacae": "Enterobacter cloacae", "hormaechei_narrow": "Enterobacter hormaechei",
    }
    occ: dict[str, set] = {g: set() for g in groups}
    for x in csv.DictReader(open(a.overlay), delimiter="\t"):
        for g, sp in groups.items():
            if x["project_species"] == sp and x["gtdb_r232_species"]:
                occ[g].add(x["gtdb_r232_species"])
    for row in rows:
        sp = row["gtdb_species"]
        row["falls_in_merged_group_gtdb_species"] = ";".join(sorted(g for g, s in occ.items() if sp and sp in s))

    def write(name: str, data: list[dict]):
        with open(a.out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)

    write("type-placement.tsv", rows)
    write("ncbi-type-candidates.tsv", cand_rows)
    write("ncbi-pages.tsv", pages)
    write("retrieval-log.tsv", rec.rows)
    json.dump({g: sorted(s) for g, s in occ.items()}, open(a.out / "merged-group-gtdb-species.json", "w"), indent=2)
    json.dump({"gtdb_taxonomy_file": str(a.r232), "gtdb_taxonomy_md5_verified": md5,
               "gtdb_taxonomy_sha256": hashlib.sha256(a.r232.read_bytes()).hexdigest(),
               "n_r232_genomes": len(exact), "run_utc": now()},
              open(a.out / "gtdb-source.json", "w"), indent=2)
    n = len(rows)
    sel_n = sum(bool(r["type_assembly"]) for r in rows)
    print(f"names={n} with_type_assembly={sel_n} without={n - sel_n} "
          f"gtdb_sources={dict((s, sum(r['gtdb_source'] == s for r in rows)) for s in sorted({r['gtdb_source'] for r in rows}))} "
          f"ncbi_queries={len({p['query'] for p in pages})} ncbi_pages={len(pages)} "
          f"incomplete={sum('INCOMPLETE' in r['ncbi_query_completeness'] for r in rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
