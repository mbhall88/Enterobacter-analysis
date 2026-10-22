#!/usr/bin/env python3
"""Validate the P6a formal-act register against its inputs and the completion criteria.

Each check is recorded with pass/fail and detail in checks/validation.json; the
script exits non-zero if any check fails. Read-only except for that JSON.

Checks r1*-r3* were added on 9 Oct 2026 for the coordinator review
(https://github.com/mbhall88/Enterobacter/issues/27#issuecomment-6059320669):
R1 exact G+C, R2 screen coverage and wording, R3 GTDB/NCBI provenance.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_register as br  # noqa: E402

PLACEMENT = "checks/type-placement-20261009"  # revised screen; checks/type-placement/ is the 8 Oct snapshot
ASTRA_EXACT = "results/project-review/2026-10-08/p6a-coordinator-review-astra/outputs/exact-type-base-counts.tsv"
TEXT_DOCS = ["manuscript-contribution.md", "outputs/nomenclature-register.md", "outputs/author-action-list.md", "README.md",
             "logs/revision-handoff-comment.md"]
# Wording the review (R2) asked to remove: each phrase overstated the screen's coverage.
OVERCLAIMS = ["type genome of every validly published name", "No type genome exists", "No type genome is available",
              "has no type genome", "all 53 validly published", "No earlier-dated type falls inside any merged group",
              "No earlier-dated type falls in any group"]
NUMWORD = {"two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9}

FROZEN = [
    "results/ani-cases/2026-10-03/task12-astra-cp102-approved-final-reconciliation/accession-assignments.tsv",
    "manifests/study-candidates.tsv",
    "results/reference-panel/2026-09-17/task02-approved-manifests-astra/types.tsv",
    "docs/review-decisions-manuscript-map.md",
    "docs/meeting-decisions-2026-10-01-session2.md",
    "results/project-review/2026-10-01/task16-astra-item6-comment/naming-rules.md",
]


def tsv(p: Path) -> list[dict]:
    return list(csv.DictReader(open(p), delimiter="\t"))


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--run", required=True, type=Path)
    a = ap.parse_args()
    run, repo = a.run, a.repo
    out = run / "outputs"
    reg = {r["register_id"]: r for r in tsv(out / "formal-act-register.tsv")}
    tmap = tsv(out / "treatment-map.tsv")
    tdoc = {r["register_id"]: r for r in tsv(out / "type-documentation.tsv")}
    prio = tsv(out / "priority-check.tsv")
    epi = tsv(out / "epithet-uniqueness.tsv")
    names = {r["lpsn_path"]: r for r in tsv(out / "name-status.tsv")}
    overlay = tsv(repo / FROZEN[0])
    lpsn_ret = tsv(run / "sources/lpsn/retrieval.tsv")
    checks = []

    def check(cid, criterion, ok, detail):
        checks.append({"id": cid, "criterion": criterion, "pass": bool(ok), "detail": detail})

    # (a) every adopted treatment mapped
    n_map = sum(int(r["n_records"]) for r in tmap)
    unknown = sorted({r["register_id"] for r in tmap} - set(reg))
    remap = collections.Counter(br.map_row(r) for r in overlay)
    check("a1", "a", n_map == len(overlay) and not unknown,
          f"overlay rows={len(overlay)}; mapped via treatment-map={n_map}; unknown register ids={unknown}")
    labels = {(r["project_species"], r["project_subspecies"]) for r in overlay if r["project_species"] or r["rule_id"]}
    check("a2", "a", all(sum(1 for t in tmap if (t["project_species"], t["project_subspecies"]) == lab) >= 1 for lab in labels),
          f"distinct adopted label pairs={len(labels)} all present in treatment-map")
    cls = collections.Counter(r["entry_class"].split(":")[0] for r in reg.values())
    unnamed = {k: v["entry_class"] for k, v in reg.items() if v["entry_class"].startswith(("unnamed", "unassigned", "infrasubspecific"))}
    check("a3", "a", any(v == "unnamed:included_no_naming" for v in unnamed.values()) and any(v == "unnamed:deferred_naming" for v in unnamed.values()),
          f"entry classes={dict(cls)}; included unnamed={sorted(k for k, v in unnamed.items() if v == 'unnamed:included_no_naming')}; "
          f"deferred={sorted(k for k, v in unnamed.items() if v == 'unnamed:deferred_naming')}")
    noact_missing = [k for k, v in reg.items() if not v["act_required"]]
    check("a4", "a", not noact_missing, f"entries without explicit act_required: {noact_missing}")
    unmapped_entries = [k for k, v in reg.items() if int(v["n_overlay_records"]) == 0 and v["entry_class"] not in
                        ("no_act:existing_name_reinstated", "no_act:circumscription_broadened", "no_act:existing_names_match", "infrasubspecific:informal")]
    check("a5", "a", not unmapped_entries, f"register entries with no overlay rows (excluding parent/informal entries): {unmapped_entries}")

    # (b) proposed acts carry rule/name/type evidence and a gap disposition; no unsupported formal status
    acts = {k: v for k, v in reg.items() if v["act_required"].startswith("yes")}
    need = ["proposed_name_and_citation", "basonym", "homotypic_synonyms", "heterotypic_synonyms", "priority_basis",
            "icnp_rule_basis", "type_strain", "type_genome", "original_description_and_emendations",
            "documentary_gaps", "gap_disposition", "author_action_ids", "evidence"]
    missing = {k: [c for c in need if not v[c].strip()] for k, v in acts.items()}
    missing = {k: m for k, m in missing.items() if m}
    check("b1", "b", len(acts) == 8 and not missing, f"acts={sorted(acts)}; missing fields={missing}")
    bad_status = [k for k, v in acts.items() if v["formal_status_statement"] != br.PENDING]
    check("b2", "b", not bad_status, f"every proposed act states it is not validly published until IJSEM publication; violations={bad_status}")
    claims = [k for k, v in reg.items() if not v["act_required"].startswith("yes") and "comb. nov." in v["proposed_name_and_citation"]]
    check("b3", "b", not claims, f"no-act entries styled as nomenclatural acts: {claims}")
    ret_bad = []
    for k, v in reg.items():
        if v["entry_class"] == "no_act:name_unchanged":
            p = "/species/enterobacter-" + k.split("P6A-RET-")[1]
            n = names[p]
            if "validly published under the ICNP" not in n["nomenclatural_status"] or not n["taxonomic_status"].startswith("correct name"):
                ret_bad.append(k)
    check("b4", "b", not ret_bad, f"retained names validly published and correct in LPSN snapshot; failures={ret_bad}")
    nvp = names["/species/enterobacter-adelaidei"]
    check("b5", "b", nvp["nomenclatural_status"] == "not validly published" and "not validly published" in reg["P6A-NVP-01"]["proposed_name_and_citation"],
          f"'E. adelaidei' LPSN status={nvp['nomenclatural_status']!r}; register keeps quotation marks/qualification")
    epi_bad = [r["register_id"] for r in epi if r["distinct_type_strain_lists"] != "1"]
    check("b6", "b", not epi_bad, f"Rule 12b/13c: each proposed epithet used in Enterobacter only with one type; violations={epi_bad}")
    sub_paths = [r["path"] for r in lpsn_ret if r["path"].startswith(("/subspecies/enterobacter-xiangfangensis-",
                                                                         "/subspecies/enterobacter-roggenkampii-",
                                                                         "/subspecies/enterobacter-chengduensis-"))]
    crawled = {r["path"] for r in lpsn_ret}
    parents_crawled = all(p in crawled for p in ("/species/enterobacter-xiangfangensis", "/species/enterobacter-roggenkampii", "/species/enterobacter-chengduensis"))
    check("b7", "b", parents_crawled and not sub_paths,
          f"no existing LPSN subspecies under the three target species (species pages crawled for children={parents_crawled}; found={sub_paths})")
    place_bad = {k: d["overlay_register_id"] for k, d in tdoc.items() if d["overlay_register_id"] != k}
    check("b8", "b", len(tdoc) == 8 and not place_bad, f"each act's type genome sits in that act's overlay group; mismatches={place_bad}")
    by_acc = {r["accession"]: br.map_row(r) for r in overlay}
    syn_bad = {acc: by_acc.get(acc) for rid, lst in br.SYN_TYPES.items() for acc, _ in lst if by_acc.get(acc) != rid}
    check("b9", "b", not syn_bad, f"heterotypic-synonym types (intestinihominis, nototheniae) inside P6A-ACT-03; mismatches={syn_bad}")
    nontype = [k for k, d in tdoc.items() if d["ncbi_type_material"] != "assembly from type material"]
    check("b10", "b", not nontype, f"NCBI flags every act's type genome as type material; failures={nontype}")
    gap_bad = [k for k, d in tdoc.items() if int(d["n_listed_countries"]) < 2 and "resolvable_before_submission" not in reg[k]["gap_disposition"]]
    check("b11", "b", not gap_bad, f"acts with <2 catalogue-verified countries carry a pre-submission disposition; failures={gap_bad}; "
                                   f"verified-two-country={sorted(k for k, d in tdoc.items() if int(d['n_listed_countries']) >= 2)}")

    # (c) decisions retained; priority consistent; contradictions
    ga = [r for r in overlay if r["rule_id"] == "group-a"]
    check("c1", "c", len(ga) == 173 and all(r["project_species"] == "Enterobacter xiangfangensis" and r["project_subspecies"] == "hoffmannii" for r in ga)
          and "subsp. hoffmannii (Sutton et al. 2018)" in reg["P6A-ACT-03"]["proposed_name_and_citation"] and "DSM 14563" in reg["P6A-ACT-03"]["type_strain"],
          f"Group A rows={len(ga)}, all one hoffmannii subspecies; ACT-03 retains type DSM 14563")
    grp = collections.defaultdict(list)
    for r in prio:
        grp[r["adopted_taxon_rank"]].append(r)
    prio_bad = {}
    for g, rows in grp.items():
        adopted = [r for r in rows if r["note"].startswith("adopted")]
        others = [r for r in rows if not r["note"].startswith("adopted")]
        if len(adopted) != 1 or any(o["order_key"] <= adopted[0]["order_key"] for o in others):
            prio_bad[g] = [(r["name"], r["order_key"]) for r in rows]
    check("c2", "c", not prio_bad, f"adopted name/epithet is earliest in every merged group ({len(grp)} groups); failures={prio_bad}")
    tp = tsv(run / PLACEMENT / "type-placement.tsv")
    earlier = []
    adopted_dates = {"xiangfangensis_A+B": "2014", "roggenkampii": "2018", "chengduensis": "2019", "mori": "2011",
                     "asburiae": "1988", "cloacae": "1960", "hormaechei_narrow": "1990"}
    for r in tp:
        for g in filter(None, r["falls_in_merged_group_gtdb_species"].split(";")):
            yr = re.findall(r"\b(1[89]\d\d|20\d\d)\b", r["lpsn_name"].replace("(Approved Lists 1980)", ""))[-1]
            if yr < adopted_dates[g]:
                earlier.append((g, r["lpsn_name"]))
    check("c3", "c", not earlier, f"no validly published Enterobacter species type with an earlier date falls in a merged group's GTDB species; hits={earlier}; "
                                 f"names without a selected type assembly (not screened): {[r['lpsn_name'] for r in tp if not r['type_assembly']]}")
    spec = {"group-b": "Enterobacter xiangfangensis", "group-a": "Enterobacter xiangfangensis", "roggenkampii-quasi": "Enterobacter roggenkampii",
            "chengduensis-rongchengensis": "Enterobacter chengduensis", "mori-quasimori": "Enterobacter mori", "asburiae-dykesii": "Enterobacter asburiae",
            "cloacae-dissolvens": "Enterobacter cloacae", "narrow-hormaechei": "Enterobacter hormaechei"}
    dec_bad = [r["accession"] for r in overlay if r["rule_id"] in spec and r["project_species"] != spec[r["rule_id"]]]
    check("c4", "c", not dec_bad, f"species decisions in the overlay match the register's adopted species; mismatches={dec_bad[:5]}")

    # (d) deliverables present; posting/no-contact constraints
    need_files = ["outputs/nomenclature-register.md", "outputs/author-action-list.md", "manuscript-contribution.md", "README.md"]
    absent = [f for f in need_files if not (run / f).exists()]
    check("d1", "d", not absent, f"deliverables present; missing={absent}")
    posts = sorted((run / "logs").glob("*comment*.md"))
    bad_post = [p.name for p in posts if not open(p).read().startswith("**Posted by:** Claude Opus 5.5 (`claude-opus-5-5[1m]`)")]
    check("d2", "d", posts and not bad_post, f"GitHub post drafts start with the required attribution: {[p.name for p in posts]}; failures={bad_post}")

    # frozen inputs unchanged
    rec = dict(line.split()[::-1] for line in open(run / "inputs.sha256") if line.strip())
    fro = []
    for f in FROZEN:
        clean = subprocess.run(["git", "-C", str(repo), "diff", "--quiet", "HEAD", "--", f]).returncode == 0
        fro.append((f, clean, sha(repo / f) == rec.get(f)))
    check("f1", "frozen", all(c and h for _, c, h in fro), f"frozen inputs clean vs HEAD and equal to recorded sha256: {fro}")
    circ = tsv(run / "checks/circumscription/circumscription-stats.tsv")
    check("f2", "evidence", all(r["n_excluded"] == "0" and r["n_measured"] == r["n_members_overlay"] for r in circ),
          f"circumscription stats measured every member with verified hashes ({len(circ)} taxa)")

    # (R1) exact G+C from integer base counts
    tg = tsv(run / "checks/type-assemblies/type-genome-gc.tsv")
    r1a = []
    for r in tg:
        rows = tsv(repo / r["per_contig_counts"])  # repo-relative path recorded by type_genome_gc.sh
        tot = [sum(int(x[k]) for x in rows) for k in ("length", "A", "C", "G", "T")]
        ok = (tot == [int(r[k]) for k in ("total_length", "A", "C", "G", "T")] and len(rows) == int(r["n_sequences"])
              and int(r["non_acgt_bases"]) == tot[0] - sum(tot[1:])
              and f"{100 * (tot[2] + tot[3]) / sum(tot[1:]):.2f}" == r["gc_mol_percent_acgt"])
        if not ok:
            r1a.append(r["accession"])
    check("r1a", "R1", len(tg) == 17 and not r1a,
          f"17 type genomes: per-contig integer counts (seqkit -C) sum to the summary and G+C = 100(G+C)/(A+C+G+T) from integers; failures={r1a}")
    ex = repo / ASTRA_EXACT
    if ex.exists():
        ast = {r["accession"]: r for r in tsv(ex)}
        r1b = [r["accession"] for r in tg if [int(r[k]) for k in "ACGT"] != [int(ast[r["accession"]][k]) for k in "ACGT"]
               or r["gc_mol_percent_acgt"] != f"{float(ast[r['accession']]['gc_exact']):.2f}"]
        check("r1b", "R1", not r1b, f"counts and G+C equal the coordinator's independent exact counts ({ASTRA_EXACT}); mismatches={r1b}")
    t = subprocess.run([sys.executable, str(HERE / "test_gc_exact.py")], capture_output=True, text=True)
    check("r1c", "R1", t.returncode == 0, "regression test scripts/test_gc_exact.py: " + t.stdout.strip().splitlines()[-1])
    mem = tsv(run / "checks/circumscription/member-genome-stats.tsv")
    r1d = []
    for r in mem:
        with gzip.open(run / "checks/circumscription/base-counts" / f"{r['accession']}.base-counts.tsv.gz", "rt") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        tot = [sum(int(x[k]) for x in rows) for k in ("length", "A", "C", "G", "T")]
        if tot != [int(r[k]) for k in ("total_length", "A", "C", "G", "T")] or f"{100 * (tot[2] + tot[3]) / sum(tot[1:]):.2f}" != r["gc_mol_percent_acgt"]:
            r1d.append(r["accession"])
    check("r1d", "R1", len(mem) == 824 and not r1d, f"{len(mem)} circumscription members: per-contig integer counts reproduce size and G+C; failures={r1d[:5]}")
    gcmap = {r["accession"]: r for r in tg}
    circ_by_n = collections.defaultdict(list)
    for r in circ:
        circ_by_n[int(r["n_measured"])].append(r)
    txt = open(run / "manuscript-contribution.md").read()
    r1e, n_type, n_range = [], 0, 0
    for acc, bp, gcv in re.findall(r"(GCA_\d+\.\d) \(([\d,]+) bp; (?:G\+C )?([\d.]+) mol%\)", txt):
        n_type += 1
        if acc not in gcmap or int(bp.replace(",", "")) != int(gcmap[acc]["total_length"]) or gcv != gcmap[acc]["gc_mol_percent_acgt"]:
            r1e.append(acc)
    pat = (r"(\d+|two|three|four|five|six|seven|eight|nine) genomes[.;] (?:Genome sizes range from |Genome sizes are |size )"
           r"([\d.]+)(?: to |–)([\d.]+) Mb[,;] (?:and )?(?:genomic )?G\+C(?: content)?(?: from)? ([\d.]+)(?: to |–)([\d.]+) mol%")
    for n, s0, s1, g0, g1 in re.findall(pat, txt):
        n_range += 1
        n = NUMWORD.get(n) or int(n)
        if not any(f"{int(r['size_min_bp']) / 1e6:.2f}" == s0 and f"{int(r['size_max_bp']) / 1e6:.2f}" == s1
                   and r["gc_min"] == g0 and r["gc_max"] == g1 for r in circ_by_n[n]):
            r1e.append(f"range n={n} {g0}-{g1}")
    check("r1e", "R1", not r1e and n_type >= 14 and n_range == 10,
          f"manuscript type-genome values ({n_type}) and taxon ranges ({n_range}) equal the recomputed tables; mismatches={r1e}")

    # (R2) metadata competitor screen: coverage, missingness and wording
    cov = tsv(out / "priority-screen-coverage.tsv")
    lp = [r for r in tsv(run / "sources/lpsn/lpsn-records.tsv") if r["path"].startswith("/species/enterobacter-")
          and "validly published under the ICNP" in r["Nomenclatural status"]]
    st = collections.Counter(r["coverage_status"] for r in cov)
    miss = {r["lpsn_name"].split(" (")[0].split(" ")[1]: (r["name_category"], r["coverage_status"]) for r in cov
            if r["coverage_status"] != "type_assembly_selected"}
    want = {"nimipressuralis": ("historical_combination_now_Lelliottia", "ncbi_type_candidates_failed_alias_match"),
            "pulveris": ("historical_combination_now_Franconibacter", "no_ncbi_type_material_candidate"),
            "pyrinus": ("historical_combination_now_Pluralibacter", "no_ncbi_type_material_candidate"),
            "siamensis": ("current_enterobacter_name", "no_ncbi_type_material_candidate")}
    check("r2a", "R2", len(cov) == len(lp) == 53 and st["type_assembly_selected"] == 49 and miss == want
          and all(r["type_assembly"] for r in cov if r["coverage_status"] == "type_assembly_selected"),
          f"audited={len(cov)} (LPSN validly published Enterobacter species records={len(lp)}); status={dict(st)}; unscreened={miss}")
    pages = tsv(run / PLACEMENT / "ncbi-pages.tsv")
    byq = collections.defaultdict(list)
    for p in pages:
        byq[p["query"]].append(p)
    inc = [q for q, ps in byq.items() if ps[-1]["next_page_token_present"] != "False"
           or sum(int(p["n_reports"]) for p in ps) != int(ps[0]["total_count"] or 0)]
    check("r2b", "R2", not inc and all("INCOMPLETE" not in r["ncbi_query_completeness"] for r in cov),
          f"every NCBI type-material query followed next_page_token to the end and returned total_count reports ({len(byq)} queries, {len(pages)} pages); incomplete={inc}")
    over = {d: [p for p in OVERCLAIMS if p.lower() in open(run / d).read().lower()] for d in TEXT_DOCS if (run / d).exists()}
    over = {d: v for d, v in over.items() if v}
    named = all("metadata competitor screen" in open(run / d).read() for d in ("manuscript-contribution.md", "outputs/nomenclature-register.md"))
    check("r2c", "R2", not over and named,
          f"coverage overclaims removed from {TEXT_DOCS}; 'metadata competitor screen' named in Methods and register; hits={over}")

    # (R3) GTDB/NCBI provenance
    md5_official = next(line.split()[0] for line in open(run / "sources/gtdb/MD5SUM.txt") if line.split()[1].lstrip("./") == "bac120_taxonomy_r232.tsv.gz")
    md5_file = hashlib.md5(open(run / "sources/gtdb/bac120_taxonomy_r232.tsv.gz", "rb").read()).hexdigest()
    gsrc = json.load(open(run / PLACEMENT / "gtdb-source.json"))
    check("r3a", "R3", md5_official == md5_file == gsrc["gtdb_taxonomy_md5_verified"] and "v232" in open(run / "sources/gtdb/VERSION.txt").read(),
          f"GTDB R232 bac120 taxonomy MD5 {md5_file} equals the release MD5SUM.txt entry and the screen's record; VERSION.txt v232")
    sel = [r for r in tp if r["type_assembly"]]
    mislabel = [r["lpsn_name"] for r in tp if r["gtdb_species"] and not (
        (r["gtdb_source"].startswith("gtdb_r232_full_taxonomy") and r["gtdb_release"].startswith("GTDB R232") and md5_file in r["gtdb_release"])
        or (r["gtdb_source"] == "gtdb_api_card_unversioned" and r["gtdb_release"].startswith("unversioned") and "R232" not in r["gtdb_release"]))]
    hdr = open(run / PLACEMENT / "type-placement.tsv").readline()
    check("r3b", "R3", not mislabel and "gtdb_r232_species" not in hdr.rstrip("\n").split("\t") and len(sel) == 49 and all(r["gtdb_species"] for r in sel),
          f"every selected type assembly's GTDB species carries its release; R232 only from the MD5-verified file; "
          f"sources={dict(collections.Counter(r['gtdb_source'] for r in tp))}; mislabelled={mislabel}")
    rl = tsv(run / PLACEMENT / "retrieval-log.tsv")
    raw_bad = [r["raw_file"] for r in rl if r["error"] or not r["raw_file"]
               or hashlib.sha256(gzip.decompress((run / PLACEMENT / r["raw_file"]).read_bytes())).hexdigest() != r["raw_sha256"]]
    cand = tsv(run / PLACEMENT / "ncbi-type-candidates.tsv")
    ncbi_sel = [r for r in tp if r["selection"] == "ncbi_type_material_alias_match"]
    sel_ok = all(any(c["lpsn_name"] == r["lpsn_name"] and c["accession"] == r["type_assembly"] and c["selected"] == "True" for c in cand) for r in ncbi_sel)
    check("r3c", "R3", rl and not raw_bad and sel_ok,
          f"{len(rl)} raw responses saved with matching SHA-256; every NCBI-selected type ({len(ncbi_sel)}) traced to a saved candidate row; failures={raw_bad}")
    snap = run / "checks/type-placement"
    before = dict(line.split()[::-1] for line in open(run / "superseded/2026-10-08-pre-review.sha256"))
    kept = all(sha(snap / f) == before.get(f"./checks/type-placement/{f}") for f in ("type-placement.tsv", "api-log.json", "merged-group-gtdb-species.json"))
    check("r3d", "R3", kept, "the 8 Oct snapshot (checks/type-placement/) is unchanged; the refreshed screen is written separately to " + PLACEMENT)

    res = {"n_checks": len(checks), "n_pass": sum(c["pass"] for c in checks), "checks": checks}
    (run / "checks").mkdir(exist_ok=True)
    json.dump(res, open(run / "checks/validation.json", "w"), indent=1, ensure_ascii=False)
    for c in checks:
        print(("PASS" if c["pass"] else "FAIL"), c["id"], c["detail"][:220])
    print(f"{res['n_pass']}/{res['n_checks']} checks passed")
    return 0 if res["n_pass"] == res["n_checks"] else 1


if __name__ == "__main__":
    sys.exit(main())
