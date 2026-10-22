#!/usr/bin/env python3
"""Snapshot public culture-collection catalogue entries for act-relevant type strains.

Read-only GETs of public catalogue pages (no ordering, login or contact). Only
collections whose public pages are server-rendered and parseable are queried
(DSMZ, CRBIP/CIP, JCM, CCUG); other designations are recorded as not checked.
A listing is evidence of a current catalogue entry on the retrieval date, not a
certificate of deposit or proof of viability at the future publication date.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
import re
import sys
import time
import urllib.request
from pathlib import Path

UA = "Mozilla/5.0 (P6a nomenclature register; read-only catalogue snapshot)"
COUNTRY = {"DSM": "Germany", "CIP": "France", "LMG": "Belgium", "CCM": "Czech Republic", "CCUG": "Sweden",
           "NCIMB": "United Kingdom", "JCM": "Japan", "KCTC": "Republic of Korea", "GDMCC": "China",
           "CCTCC": "China", "CGMCC": "China", "NCTC": "United Kingdom", "ATCC": "USA"}

# Type strains of the names involved in proposed acts (LPSN designations).
STRAINS = {
    "E. xiangfangensis (type subspecies of the proposed arrangement)": "10-17; CCUG 62994; LMG 27195; NCIMB 14836",
    "E. hormaechei subsp. oharae (basonym of proposed comb. nov.)": "CIP 108490; DSM 16687; EN-314",
    "E. hormaechei subsp. steigerwaltii (basonym of proposed comb. nov.)": "CIP 108489; DSM 16691; EN-562",
    "E. hormaechei subsp. hoffmannii (basonym of proposed comb. nov.)": "CCM 9140; DSM 14563; EN-114; LMG 30171",
    "E. roggenkampii (type subspecies of the proposed arrangement)": "DSM 16690; EN-117; LMG 30172",
    "E. quasiroggenkampii (basonym of proposed comb. nov.)": "GDMCC 1.1742; KCTC 52992; WCHECL1060",
    "E. chengduensis (type subspecies of the proposed arrangement)": "CCTCC AB 2017105; CCUG 76941; GDMCC 1.1216; KCTC 52993; WCHECl-C4",
    "E. rongchengensis (basonym of proposed comb. nov.)": "170250; GDMCC 1.3670; JCM 36189",
}


def url_for(coll: str, num: str) -> str | None:
    if coll == "DSM":
        return f"https://www.dsmz.de/collection/catalogue/details/culture/DSM-{num}"
    if coll == "CIP":
        return f"https://catalogue-crbip.pasteur.fr/fiche_catalogue.xhtml?crbip=CIP%20{num}"
    if coll == "JCM":
        return f"https://www.jcm.riken.jp/cgi-bin/jcm/jcm_number?JCM={num}"
    if coll == "CCUG":
        return f"https://www.ccug.se/strain?id={num}"
    return None


def to_text(raw: str) -> str:
    t = re.sub(r"<script.*?</script>|<style.*?</style>|<!--.*?-->", "", raw, flags=re.S)
    t = re.sub(r"<[^>]+>", "\n", t)
    t = html.unescape(t)
    return re.sub(r"\s+", " ", t)


def assess(coll: str, num: str, t: str) -> tuple[str, str]:
    n = re.escape(num)
    if coll == "DSM":
        ok = re.search(r"DSM No\.:\s*" + n + r", Type strain", t)
        order = "Add to Cart" in t
        return ("listed_type_strain" if ok else "not_found"), ("orderable_marker_present" if order else "no_order_marker")
    if coll == "CIP":
        ok = re.search(r"CIP " + n + r"T", t) and "Type strain" in t
        return ("listed_type_strain" if ok else "not_found"), ("order_marker_present" if "Order" in t else "no_order_marker")
    if coll == "JCM":
        ok = re.search(n + r"T\s*<--", t) or (re.search(n, t) and "Type strain" in t)
        avail = re.search(r"A \(Freeze-dried or L-dried culture\)", t)
        return ("listed_type_strain" if ok else "not_found"), ("supply_form_listed" if avail else "no_supply_form")
    if coll == "CCUG":
        ok = re.search(r"CCUG " + n + r"T", t)
        return ("listed_type_strain" if ok else "not_found"), ("type_strain_collection_page" if "Type Strain" in t else "")
    return "not_checked", ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for taxon, desig in STRAINS.items():
        for d in [x.strip() for x in desig.split(";")]:
            m = re.match(r"([A-Z]+)\s+(.+)$", d)
            coll, num = (m.group(1), m.group(2)) if m else ("", d)
            if coll not in COUNTRY:
                rows.append({"taxon": taxon, "designation": d, "collection": coll or "strain designation",
                             "country": "", "url": "", "retrieved_utc": "", "http_status": "", "sha256": "",
                             "catalogue_result": "not_a_collection_number", "supply_marker": ""})
                continue
            url = url_for(coll, num.replace(" ", "%20") if coll != "CIP" else num)
            if not url:
                rows.append({"taxon": taxon, "designation": d, "collection": coll, "country": COUNTRY[coll], "url": "",
                             "retrieved_utc": "", "http_status": "", "sha256": "",
                             "catalogue_result": "not_checked_no_parseable_public_page", "supply_marker": ""})
                continue
            when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
            try:
                req = urllib.request.Request(url, headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=60) as r:
                    body, status = r.read(), r.status
            except Exception as e:
                rows.append({"taxon": taxon, "designation": d, "collection": coll, "country": COUNTRY[coll], "url": url,
                             "retrieved_utc": when, "http_status": "error", "sha256": "",
                             "catalogue_result": "fetch_error:" + repr(e)[:80], "supply_marker": ""})
                continue
            fn = raw / (re.sub(r"[^A-Za-z0-9.]+", "_", d) + ".html")
            fn.write_bytes(body)
            res, mark = assess(coll, num, to_text(body.decode("utf-8", "replace")))
            rows.append({"taxon": taxon, "designation": d, "collection": coll, "country": COUNTRY[coll], "url": url,
                         "retrieved_utc": when, "http_status": status, "sha256": hashlib.sha256(body).hexdigest(),
                         "catalogue_result": res, "supply_marker": mark})
            time.sleep(0.6)
    with open(a.out / "collection-checks.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # per-taxon summary: countries with a verified listing vs designated countries
    summ = []
    for taxon in STRAINS:
        rs = [r for r in rows if r["taxon"] == taxon]
        verified = sorted({r["country"] for r in rs if r["catalogue_result"] == "listed_type_strain"})
        designated = sorted({r["country"] for r in rs if r["country"]})
        summ.append({"taxon": taxon, "designated_collection_countries": "; ".join(designated),
                     "n_designated_countries": len(designated),
                     "catalogue_verified_countries": "; ".join(verified), "n_verified_countries": len(verified),
                     "two_country_listing_verified_here": "yes" if len(verified) >= 2 else "no"})
    with open(a.out / "collection-summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    print("\n".join(f"{s['taxon']}: verified {s['catalogue_verified_countries'] or '-'} of {s['designated_collection_countries']}" for s in summ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
