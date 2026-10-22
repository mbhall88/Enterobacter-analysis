#!/usr/bin/env python3
"""Current NCBI records for type assemblies named in the register.

Read-only NCBI Datasets v2 accession reports. Values are NCBI-reported
(assembly statistics, type-material flag, version status); nothing is
recomputed from sequence. Raw JSON responses are kept with SHA-256.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

UA = "Mozilla/5.0 (P6a nomenclature register; read-only)"
TYPES = {
    "Enterobacter xiangfangensis / E. hormaechei subsp. xiangfangensis": "GCA_001729785.1",
    "Enterobacter hormaechei subsp. oharae": "GCA_001729705.1",
    "Enterobacter hormaechei subsp. steigerwaltii": "GCA_001729725.1",
    "Enterobacter hormaechei subsp. hoffmannii / E. hoffmannii": "GCA_001729745.1",
    "Enterobacter intestinihominis": "GCA_040096145.1",
    "Enterobacter hoffmannii subsp. nototheniae": "GCA_052922205.1",
    "Enterobacter hormaechei": "GCA_019048245.1",
    "Enterobacter roggenkampii": "GCA_001729805.1",
    "Enterobacter quasiroggenkampii": "GCA_001039365.2",
    "Enterobacter chengduensis": "GCA_001984825.2",
    "Enterobacter rongchengensis": "GCA_034008275.1",
    "Enterobacter mori": "GCA_000211415.1",
    "Enterobacter quasimori": "GCA_003964905.1",
    "Enterobacter asburiae": "GCA_001521715.1",
    "Enterobacter dykesii": "GCA_008364625.2",
    "Enterobacter cloacae": "GCA_000025565.1",
    "Enterobacter cloacae subsp. dissolvens": "GCA_009660455.1",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for taxon, acc in TYPES.items():
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report"
        when = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=90) as r:
            body = r.read()
        (raw / f"{acc}.json").write_bytes(body)
        d = json.loads(body)
        rep = (d.get("reports") or [{}])[0]
        st = rep.get("assembly_stats", {})
        ai = rep.get("assembly_info", {})
        bs = ai.get("biosample", {}) or {}
        attrs = {x.get("name"): x.get("value", "") for x in bs.get("attributes", [])}
        rows.append({
            "taxon": taxon, "accession": acc, "retrieved_utc": when, "sha256": hashlib.sha256(body).hexdigest(),
            "ncbi_organism": rep.get("organism", {}).get("organism_name", ""),
            "strain": (rep.get("organism", {}).get("infraspecific_names") or {}).get("strain", ""),
            "culture_collection": attrs.get("culture_collection", ""),
            "type_material": (rep.get("type_material") or {}).get("type_display_text", ""),
            "assembly_level": ai.get("assembly_level", ""), "assembly_status": ai.get("assembly_status", ""),
            "paired_accession": rep.get("paired_accession", ""),
            "total_sequence_length": st.get("total_sequence_length", ""),
            "gc_percent": st.get("gc_percent", ""), "number_of_contigs": st.get("number_of_contigs", ""),
            "biosample": bs.get("accession", ""), "bioproject": (ai.get("bioproject_lineage") or [{}])[0].get("bioprojects", [{}])[0].get("accession", ""),
            "ani_check_status": (rep.get("average_nucleotide_identity") or {}).get("taxonomy_check_status", ""),
            "ani_match_status": (rep.get("average_nucleotide_identity") or {}).get("match_status", ""),
        })
        time.sleep(0.4)
    with open(a.out / "type-assemblies-ncbi.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"records={len(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
