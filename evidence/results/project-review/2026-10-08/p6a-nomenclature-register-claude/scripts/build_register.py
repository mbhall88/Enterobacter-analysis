#!/usr/bin/env python3
"""Build the P6a formal-act register from fixed inputs and saved primary sources.

Inputs are read-only: the frozen final accession overlay (adopted treatments),
the LPSN snapshot taken in this run, the type-genome and catalogue checks of
this run. Curated register entries below carry their evidence pointers; every
statement of nomenclatural status is traceable to a saved source (LPSN record,
Validation List, protologue or ICNP 2025 text). Nothing here changes a
taxonomic decision: adopted groupings are inputs, and proposed acts remain
proposals until validly published in IJSEM.

Outputs (in --out):
  formal-act-register.tsv   one row per adopted treatment / register entry
  treatment-map.tsv         every overlay label combination -> register entry
  name-status.tsv           LPSN status of every name the register cites
  type-documentation.tsv    type strain, deposits, catalogue check, type genome
  priority-check.tsv        competing names per adopted taxon, ordered by date
  epithet-uniqueness.tsv    Rule 12b/13c check for each proposed epithet
  priority-screen-coverage.tsv  every audited name: screened or not, and why (added 9 Oct 2026, review R2)
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

DEC = "https://github.com/mbhall88/Enterobacter/issues/2#issuecomment-"
D_S2 = DEC + "5925503549"       # DEC-004 session 2 (A+B, roggenkampii, chengduensis, mori, adelaidei, isolateB)
D_S1 = DEC + "5923562552"       # DEC-004 session 1 (A/B one species, cloacae/dissolvens, asburiae/dykesii, bugandensis/nematophilus)
D_GA = DEC + "5947129754"       # Group A one hoffmannii subspecies, GV numbering
D_CD = DEC + "5950936288"       # cloacae/dissolvens reaffirmed; C3 uncertainty accepted
D_RES = DEC + "5925858585"      # cloacae_P, nimipressuralis, residual genus
D_RHB = DEC + "5942807222"      # RHBSTW-00994 unnamed lineage
D_SEQ = DEC + "5881570584"      # SeqCode/Candidatus deferred exploration
D_PSE = DEC + "5950340532"
D_VON = DEC + "5950557263"
D_KOB = DEC + "5960193786"
D_LUD = DEC + "5961101406"
D_QH = DEC + "5963112449"
D_SOLI = DEC + "5963404731"
MAP = "docs/review-decisions-manuscript-map.md"

ICNP = "ICNP 2025 Revision (Oren et al. 2026, IJSEM 76:006979; sources/icnp)"
VL230 = "VL230 (Oren & Göker, IJSEM 2026;76(7):007224; sources/ijsem/validation-list-230-PMC13426511.html)"
VL225 = "VL225 (IJSEM 2025; sources/ijsem/validation-list-225-PMC12483762.html)"
VL221 = "VL221 (IJSEM 2025;75(1):006562; sources/ijsem/validation-list-221-PMC12281875.html)"
DEP_REQ = ("Rules 27(3), 30(3)(b), 30(4) and 32a(4): type strain of a species/subspecies, including new combinations, "
           "deposited in >=2 publicly accessible collections in different countries; evidence that cultures are present, "
           "viable and available at IJSEM publication; no deposit requiring third-party/national-authority permission")
PENDING = "Proposed only; not validly published until published in IJSEM (Rule 27(1), Note 1)"

# ----------------------------------------------------------------------------
# Register entries. Keys used for every row; '' where not applicable.
# ----------------------------------------------------------------------------
COLS = ["register_id", "entry_class", "adopted_treatment", "decision_source", "parent_entry",
        "act_required", "proposed_name_and_citation", "rank", "basonym", "homotypic_synonyms",
        "heterotypic_synonyms", "priority_basis", "icnp_rule_basis", "type_strain", "type_genome",
        "original_description_and_emendations", "documentary_gaps", "gap_disposition",
        "author_action_ids", "formal_status_statement", "evidence"]

E: list[dict] = []


def add(**kw):
    row = {c: "" for c in COLS}
    row.update(kw)
    E.append(row)


# ---- Proposed new combinations --------------------------------------------
add(register_id="P6A-ACT-01", entry_class="proposed_act:comb_nov",
    adopted_treatment="Group B subdivision oharae retained as a subspecies of E. xiangfangensis (A+B species)",
    decision_source=D_S2, parent_entry="P6A-SPN-01", act_required="yes",
    proposed_name_and_citation="Enterobacter xiangfangensis subsp. oharae (Hoffmann et al. 2016) [authors of this paper] comb. nov.",
    rank="subspecies",
    basonym="Enterobacter hormaechei subsp. oharae Hoffmann et al. 2016 (VL172, IJSEM 66:4299-4305; effective: J Clin Microbiol 2005;43:3297-3303)",
    homotypic_synonyms="Enterobacter hormaechei subsp. oharae Hoffmann et al. 2016 (basonym)",
    heterotypic_synonyms="none identified (no other validly published type falls in this subspecies)",
    priority_basis="Only epithet with a type in the adopted subspecies; 12b/13c: no other Enterobacter epithet 'oharae'",
    icnp_rule_basis="34a (comb. nov. on transfer to another species); 23a Note 2 (subspecific epithet retained on transfer); 27(2)(a)-(d) incl. 2023 synonym listing; 27(3)/30(3)(b)/30(4)/32a(4) deposits; 33a/34b citation",
    type_strain="EN-314 = DSM 16687 = CIP 108490 (LPSN; protologue: 'DSMZ 16687' misprint noted by LPSN)",
    type_genome="GCA_001729705.1 (NCBI: assembly from type material; DSM 16687)",
    original_description_and_emendations="Hoffmann et al. 2005 J Clin Microbiol 43:3297-3303 (sources/literature/hoffmann2005-jcm-PMC1169129.html); no emendation recorded by LPSN",
    documentary_gaps="Deposit availability evidence at publication. Catalogue listing verified 8 Oct 2026 for DSM 16687 (Germany) and CIP 108490 (France): two countries listed; not a viability/availability certificate.",
    gap_disposition="resolvable_before_submission: authors obtain current availability evidence from DSMZ and CRBIP (P6A-AA-01)",
    author_action_ids="P6A-AA-01;P6A-AA-05;P6A-Q02", formal_status_statement=PENDING,
    evidence="sources/lpsn/raw/subspecies__enterobacter-hormaechei-oharae.html; checks/collections/collection-checks.tsv; checks/type-assemblies/type-assemblies-ncbi.tsv")

add(register_id="P6A-ACT-02", entry_class="proposed_act:comb_nov",
    adopted_treatment="Group B subdivision steigerwaltii retained as a subspecies of E. xiangfangensis",
    decision_source=D_S2, parent_entry="P6A-SPN-01", act_required="yes",
    proposed_name_and_citation="Enterobacter xiangfangensis subsp. steigerwaltii (Hoffmann et al. 2016) [authors of this paper] comb. nov.",
    rank="subspecies",
    basonym="Enterobacter hormaechei subsp. steigerwaltii Hoffmann et al. 2016 (VL172, IJSEM 66:4299-4305; effective: J Clin Microbiol 2005;43:3297-3303)",
    homotypic_synonyms="Enterobacter hormaechei subsp. steigerwaltii Hoffmann et al. 2016 (basonym)",
    heterotypic_synonyms="none identified",
    priority_basis="Only epithet with a type in the adopted subspecies; no other Enterobacter epithet 'steigerwaltii'",
    icnp_rule_basis="34a; 23a Note 2; 27(2)(a)-(d); 27(3)/30(3)(b)/30(4)/32a(4); 33a/34b",
    type_strain="EN-562 = DSM 16691 = CIP 108489 (LPSN; protologue 'DSMZ 16691' misprint noted by LPSN)",
    type_genome="GCA_001729725.1 (NCBI: assembly from type material; DSM 16691)",
    original_description_and_emendations="Hoffmann et al. 2005 (as for P6A-ACT-01); no emendation recorded",
    documentary_gaps="Deposit availability evidence at publication. Catalogue listing verified 8 Oct 2026 for DSM 16691 (Germany) and CIP 108489 (France).",
    gap_disposition="resolvable_before_submission (P6A-AA-01)",
    author_action_ids="P6A-AA-01;P6A-AA-05;P6A-Q02", formal_status_statement=PENDING,
    evidence="sources/lpsn/raw/subspecies__enterobacter-hormaechei-steigerwaltii.html; checks/collections/collection-checks.tsv")

add(register_id="P6A-ACT-03", entry_class="proposed_act:comb_nov",
    adopted_treatment="All Group A (173 genomes incl. hoffmannii, intestinihominis and nototheniae types) = one subspecies of E. xiangfangensis named hoffmannii, type DSM 14563",
    decision_source=D_GA + " ; " + D_CD, parent_entry="P6A-SPN-01", act_required="yes",
    proposed_name_and_citation="Enterobacter xiangfangensis subsp. hoffmannii (Sutton et al. 2018) [authors of this paper] comb. nov.",
    rank="subspecies",
    basonym="Enterobacter hormaechei subsp. hoffmannii Sutton et al. 2018 (VL184, IJSEM 68:3379-3393; effective: F1000Res 2018;7:521). Also the basonym of E. hoffmannii (Sutton et al. 2018) Sedláček et al. 2026 (VL230)",
    homotypic_synonyms=("Enterobacter hormaechei subsp. hoffmannii Sutton et al. 2018 (basonym); Enterobacter hoffmannii (Sutton et al. 2018) Sedláček et al. 2026 "
                        "(VL230, priority no. 26); Enterobacter hoffmannii subsp. hoffmannii (Sutton et al. 2018) Sedláček et al. 2026 (VL230); "
                        "'Enterobacter hoffmannii' (Sutton et al. 2018) Wu et al. 2020 (not validly published)"),
    heterotypic_synonyms=("Enterobacter intestinihominis Hitch et al. 2025 (VL225, priority no. 28; type CLA-AC-H004 = DSM 118557 = LMG 33610); "
                          "Enterobacter hoffmannii subsp. nototheniae Sedláček et al. 2026 (VL230; type P5473 = CCM 8629 = LMG 34031)"),
    priority_basis=("Rule 42/23a Note 2: earliest legitimate epithet among types in Group A is hoffmannii (2018, VL184); intestinihominis 2025 (VL225), "
                    "hoffmannii as species 2026 (VL230, homotypic), nototheniae 2026 (VL230). Epithet 'hoffmannii' as E. hoffmannii is homotypic, so Rule 12b is not contravened"),
    icnp_rule_basis="34a; 23a Note 2; 42 (union at subspecies rank, earliest epithet); 12b (same type); 27(2)(b) synonym list incl. heterotypic synonyms; 27(3)/30(3)(b)/30(4)/32a(4); 33a/34b",
    type_strain="EN-114 = DSM 14563 = CCM 9140 (= LMG 30171 per protologue; VL230 footnote 22: LMG 30171 'no documentation was supplied')",
    type_genome="GCA_001729745.1 (NCBI: assembly from type material; DSM 14563)",
    original_description_and_emendations=("Sutton et al. 2018 F1000Res 7:521 (sources/literature/sutton2018-f1000res-PMC6097438.html); "
                                          "explicit description of E. hoffmannii subsp. hoffmannii in Sedláček et al. 2025 Czech Polar Rep 15:224-248 (LPSN note); no emendation recorded"),
    documentary_gaps=("Deposit availability evidence at publication. VL230 (Jul 2026) accepted documentation for DSM 14563 (Germany) and CCM 9140 (Czech Republic); "
                      "LMG 30171 undocumented. Catalogue listing verified here only for DSM 14563 (CCM, LMG pages not machine-checkable)."),
    gap_disposition="resolvable_before_submission (P6A-AA-01): cite DSM 14563 and CCM 9140; do not rely on LMG 30171 unless documented",
    author_action_ids="P6A-AA-01;P6A-AA-02;P6A-AA-05;P6A-Q02", formal_status_statement=PENDING,
    evidence="sources/lpsn/raw/subspecies__enterobacter-hormaechei-hoffmannii.html; sources/lpsn/raw/species__enterobacter-hoffmannii.html; " + VL230 + "; " + VL225)

add(register_id="P6A-ACT-04", entry_class="proposed_act:comb_nov",
    adopted_treatment="E. quasiroggenkampii included in E. roggenkampii as a subspecies",
    decision_source=D_S2, parent_entry="P6A-SPN-03", act_required="yes",
    proposed_name_and_citation="Enterobacter roggenkampii subsp. quasiroggenkampii (Wu et al. 2020) [authors of this paper] comb. nov.",
    rank="subspecies",
    basonym="Enterobacter quasiroggenkampii Wu et al. 2020 (VL196, IJSEM 70:5596-5601; effective: mSystems 2020;5:e00527-20)",
    homotypic_synonyms="Enterobacter quasiroggenkampii Wu et al. 2020 (basonym; species rank)",
    heterotypic_synonyms="none identified",
    priority_basis="Species rank: roggenkampii (VL184, 2018) precedes quasiroggenkampii (VL196, 2020), so the parent species is E. roggenkampii (Rule 42); subspecific epithet retained (Rule 50b)",
    icnp_rule_basis="34a (2025: species lowered to subspecies with epithet retained is a comb. nov.); 50b; 45; 42; 27(2)(a)-(d); 27(3)/30(3)(b)/30(4)/32a(4); 33a/34b",
    type_strain="WCHECL1060 = GDMCC 1.1742 = KCTC 52992",
    type_genome="GCA_001039365.2 (NCBI: assembly from type material; WCHECL1060)",
    original_description_and_emendations="Wu et al. 2020 mSystems 5:e00527-20 (sources/literature/wu2020-msystems-PMC7406230.html); protologue type-pair values with E. roggenkampii DSM 16690T: ANI 95.37%, isDDH 64.8% (Table 3)",
    documentary_gaps=("Neither designated collection (GDMCC, China; KCTC, Republic of Korea) has a machine-checkable public catalogue page; no listing verified. "
                      "Rule 30(4) access condition must be confirmed for both deposits."),
    gap_disposition="resolvable_before_submission but highest documentary risk (P6A-AA-01, P6A-AA-03); fallback is a new deposit by the authors if either is unavailable",
    author_action_ids="P6A-AA-01;P6A-AA-03;P6A-AA-05;P6A-AA-07;P6A-Q02", formal_status_statement=PENDING,
    evidence="sources/lpsn/raw/species__enterobacter-quasiroggenkampii.html; sources/literature/wu2020-msystems-PMC7406230.html")

add(register_id="P6A-ACT-05", entry_class="proposed_act:comb_nov",
    adopted_treatment="E. chengduensis and E. rongchengensis as one species with subspecies; parent chengduensis (supersedes 28 Sep separate-species decision)",
    decision_source=D_S2, parent_entry="P6A-SPN-04", act_required="yes",
    proposed_name_and_citation="Enterobacter chengduensis subsp. rongchengensis (He et al. 2025) [authors of this paper] comb. nov.",
    rank="subspecies",
    basonym="Enterobacter rongchengensis He et al. 2025 (VL221, IJSEM 75(1):006562, priority no. 12; effective: Microbiol Spectr 2024;12:e00292-24)",
    homotypic_synonyms="Enterobacter rongchengensis He et al. 2025 (basonym; species rank)",
    heterotypic_synonyms="none identified",
    priority_basis="Species rank: chengduensis (VL188, 2019) precedes rongchengensis (VL221, 2025) (Rule 42)",
    icnp_rule_basis="34a; 50b; 45; 42; 27(2)(a)-(d); 27(3)/30(3)(b)/30(4)/32a(4); 33a/34b",
    type_strain="170250 = GDMCC 1.3670 = JCM 36189",
    type_genome="GCA_034008275.1 (NCBI: assembly from type material; 170250)",
    original_description_and_emendations=("He et al. 2024 Microbiol Spectr 12:e00292-24 (sources/literature/he2024-spectrum-PMC11302248.html); "
                                          "protologue Table 1 type-pair values with E. chengduensis GCA_001984825: ANI 95.58%, isDDH 64.4% (the narrative text swaps the two strain IDs; the table is internally consistent)"),
    documentary_gaps="JCM 36189 (Japan) catalogue listing verified 8 Oct 2026; GDMCC 1.3670 (China) not machine-checkable; Rule 30(4) access condition to confirm.",
    gap_disposition="resolvable_before_submission (P6A-AA-01, P6A-AA-03)",
    author_action_ids="P6A-AA-01;P6A-AA-03;P6A-AA-05;P6A-AA-07;P6A-Q02", formal_status_statement=PENDING,
    evidence="sources/lpsn/raw/species__enterobacter-rongchengensis.html; " + VL221 + "; sources/literature/he2024-spectrum-PMC11302248.html")

# ---- Type subspecies --------------------------------------------------------
TSP_RULE = ("40d (type subspecies created automatically by valid publication of a subspecies excluding the species type, unless explicitly proposed in the same "
            "effective publication; Note 1: must conform to 28b and 32a); 13d/45/40b (same epithet as species); 34a and 40d example (explicit proposal styled comb. nov.); 40c Note (no 12b conflict)")
add(register_id="P6A-TSP-01", entry_class="type_subspecies:rule_40d",
    adopted_treatment="Group B subdivision containing the E. xiangfangensis type (accepted B subdivisions retained)",
    decision_source=D_S2, parent_entry="P6A-SPN-01", act_required="yes (automatic under Rule 40d, or explicit; author choice P6A-Q01)",
    proposed_name_and_citation="Enterobacter xiangfangensis subsp. xiangfangensis (Gu et al. 2014) [authors of this paper]; if explicitly proposed: comb. nov.",
    rank="subspecies",
    basonym="Enterobacter xiangfangensis Gu et al. 2014 (IJSEM 64:2650-2656)",
    homotypic_synonyms="Enterobacter xiangfangensis Gu et al. 2014 (species rank); Enterobacter hormaechei subsp. xiangfangensis (Gu et al. 2014) Sutton et al. 2018 (VL184)",
    heterotypic_synonyms="none identified",
    priority_basis="Name fixed by the species type (Rules 13d, 40b, 45)",
    icnp_rule_basis=TSP_RULE,
    type_strain="10-17 = CCUG 62994 = LMG 27195 = NCIMB 14836",
    type_genome="GCA_001729785.1 (NCBI: assembly from type material; LMG 27195)",
    original_description_and_emendations="Gu et al. 2014 IJSEM 64:2650-2656 (not open access; LPSN record saved); Sutton et al. 2018 comb. nov. as E. hormaechei subsp.; no emendation recorded",
    documentary_gaps="CCUG 62994 (Sweden) catalogue listing verified; LMG 27195 (Belgium) and NCIMB 14836 (UK) not machine-checkable. Rule 40d Note 1 applies 32a (incl. 30(3)/(4)) to the type subspecies name.",
    gap_disposition="resolvable_before_submission (P6A-AA-01)",
    author_action_ids="P6A-AA-01;P6A-Q01;P6A-Q02", formal_status_statement=PENDING,
    evidence="sources/lpsn/raw/species__enterobacter-xiangfangensis.html; sources/lpsn/raw/subspecies__enterobacter-hormaechei-xiangfangensis.html")

add(register_id="P6A-TSP-02", entry_class="type_subspecies:rule_40d",
    adopted_treatment="E. roggenkampii subspecies containing the E. roggenkampii type (81 genomes)",
    decision_source=D_S2, parent_entry="P6A-SPN-03", act_required="yes (automatic under Rule 40d, or explicit; P6A-Q01)",
    proposed_name_and_citation="Enterobacter roggenkampii subsp. roggenkampii (Sutton et al. 2018) [authors of this paper]; if explicitly proposed: comb. nov.",
    rank="subspecies", basonym="Enterobacter roggenkampii Sutton et al. 2018 (VL184, IJSEM 68:3379-3393)",
    homotypic_synonyms="Enterobacter roggenkampii Sutton et al. 2018 (species rank)", heterotypic_synonyms="none identified",
    priority_basis="Name fixed by the species type", icnp_rule_basis=TSP_RULE,
    type_strain="EN-117 = DSM 16690 = LMG 30172", type_genome="GCA_001729805.1 (NCBI: assembly from type material; DSM 16690)",
    original_description_and_emendations="Sutton et al. 2018 F1000Res 7:521; no emendation recorded",
    documentary_gaps="DSM 16690 (Germany) catalogue listing verified; LMG 30172 (Belgium) not machine-checkable; protologue states the LMG deposit.",
    gap_disposition="resolvable_before_submission (P6A-AA-01)", author_action_ids="P6A-AA-01;P6A-Q01;P6A-Q02",
    formal_status_statement=PENDING, evidence="sources/lpsn/raw/species__enterobacter-roggenkampii.html; sources/literature/sutton2018-f1000res-PMC6097438.html")

add(register_id="P6A-TSP-03", entry_class="type_subspecies:rule_40d",
    adopted_treatment="E. chengduensis subspecies containing the E. chengduensis type (3 genomes)",
    decision_source=D_S2, parent_entry="P6A-SPN-04", act_required="yes (automatic under Rule 40d, or explicit; P6A-Q01)",
    proposed_name_and_citation="Enterobacter chengduensis subsp. chengduensis (Wu et al. 2019) [authors of this paper]; if explicitly proposed: comb. nov.",
    rank="subspecies", basonym="Enterobacter chengduensis Wu et al. 2019 (VL188, IJSEM 69:1844-1846)",
    homotypic_synonyms="Enterobacter chengduensis Wu et al. 2019 (species rank)", heterotypic_synonyms="none identified",
    priority_basis="Name fixed by the species type", icnp_rule_basis=TSP_RULE,
    type_strain="WCHECl-C4 = CCTCC AB 2017105 = CCUG 76941 = GDMCC 1.1216 = KCTC 52993",
    type_genome="GCA_001984825.2 (NCBI: assembly from type material; WCHECl-C4)",
    original_description_and_emendations="Wu et al. 2019 Antonie van Leeuwenhoek 112:491-500 (not open access; LPSN record saved); no emendation recorded",
    documentary_gaps="CCUG 76941 (Sweden) catalogue listing verified; CCTCC/GDMCC (China) and KCTC (Republic of Korea) not machine-checkable.",
    gap_disposition="resolvable_before_submission (P6A-AA-01)", author_action_ids="P6A-AA-01;P6A-Q01;P6A-Q02",
    formal_status_statement=PENDING, evidence="sources/lpsn/raw/species__enterobacter-chengduensis.html")

# ---- Species-level treatments without a required nomenclatural act ---------
add(register_id="P6A-SPN-01", entry_class="no_act:existing_name_reinstated",
    adopted_treatment="Groups A+B form one species named E. xiangfangensis; B subdivisions retained; all A one hoffmannii subspecies",
    decision_source=D_S1 + " ; " + D_S2 + " ; " + D_GA, act_required="no (taxonomic opinion; validly published name with same type already exists)",
    proposed_name_and_citation="Enterobacter xiangfangensis Gu et al. 2014 (optional: emend. [authors of this paper], Rule 35; P6A-Q02)",
    rank="species", basonym="none (original name)",
    homotypic_synonyms="Enterobacter hormaechei subsp. xiangfangensis (Gu et al. 2014) Sutton et al. 2018 (distinct rank)",
    heterotypic_synonyms=("At species rank: Enterobacter intestinihominis Hitch et al. 2025; Enterobacter hoffmannii (Sutton et al. 2018) Sedláček et al. 2026. "
                          "Subspecies names transferred by P6A-ACT-01..03"),
    priority_basis=("Rule 42/24b: xiangfangensis validly published IJSEM 64:2650 (Aug 2014) precedes every other species name whose type falls in A+B "
                    "(intestinihominis 2025, hoffmannii 2026), among the 49 names whose type assembly the metadata competitor screen selected (53 audited; 4 without a selected assembly). "
                    "E. hormaechei type is outside A+B (narrow hormaechei). E. siamensis (VL159, IJSEM 64:2927; no type assembly selected) is on a later page of the same volume, "
                    "so it cannot displace xiangfangensis at species rank (Rule 24b(4)); that argument covers only xiangfangensis; "
                    "LPSN note: 'E. xiangfangensis is the correct name instead if this species is regarded as a separate species'"),
    icnp_rule_basis="23a; 24b(3),(4); 37b (circumscription change does not change a name); 38; 42; 35 (emend. optional)",
    type_strain="10-17 = CCUG 62994 = LMG 27195 = NCIMB 14836", type_genome="GCA_001729785.1",
    original_description_and_emendations="Gu et al. 2014 IJSEM 64:2650-2656; no emendation recorded",
    documentary_gaps="None for the name. The species description must cover A+B (emendation optional).",
    gap_disposition="no_gap (name available); emendation is author choice P6A-Q02", author_action_ids="P6A-Q02;P6A-AA-06",
    formal_status_statement="Validly published name; using it at species rank is a taxonomic opinion and needs no validation",
    evidence="sources/lpsn/raw/species__enterobacter-xiangfangensis.html; checks/type-placement-20261009/type-placement.tsv; outputs/priority-screen-coverage.tsv")

add(register_id="P6A-SPN-02", entry_class="no_act:circumscription_narrowed",
    adopted_treatment="Narrow E. hormaechei (13 study genomes + type); frozen hormaechei reference retained",
    decision_source=D_S2, act_required="no",
    proposed_name_and_citation="Enterobacter hormaechei O'Hara et al. 1990 (optional emend.; P6A-Q02)",
    rank="species", basonym="none",
    homotypic_synonyms="Enterobacter hormaechei subsp. hormaechei (O'Hara et al. 1990) Hoffmann et al. 2016 (type subspecies, Rule 40d, VL172) remains the only subspecies once ACT-01..03 are published",
    heterotypic_synonyms="none",
    priority_basis="Type of E. hormaechei retained in the narrowed species (Rules 40a/40b)",
    icnp_rule_basis="37b; 40a/40b; 35",
    type_strain="0992-77 = ATCC 49162 = CCUG 27126 = CIP 103441 = DSM 12409", type_genome="GCA_019048245.1",
    original_description_and_emendations="O'Hara et al. 1989 J Clin Microbiol 27:2046-2049 (VL32, IJSB 40:105-106); Hoffmann et al. 2005 subspecies; no emendation recorded",
    documentary_gaps="none for the name", gap_disposition="no_gap; author choice whether to emend and whether to cite subsp. hormaechei (P6A-Q02, P6A-Q05)",
    author_action_ids="P6A-Q02;P6A-Q05", formal_status_statement="Validly published name unchanged",
    evidence="sources/lpsn/raw/species__enterobacter-hormaechei.html; sources/lpsn/raw/subspecies__enterobacter-hormaechei-hormaechei.html")

add(register_id="P6A-SPN-03", entry_class="no_act:circumscription_broadened",
    adopted_treatment="One E. roggenkampii species (87 genomes) with quasiroggenkampii as a subspecies; ECC-249 outside",
    decision_source=D_S2, act_required="no (species name); subspecies acts P6A-ACT-04/P6A-TSP-02",
    proposed_name_and_citation="Enterobacter roggenkampii Sutton et al. 2018 (optional emend.; P6A-Q02)", rank="species", basonym="none",
    heterotypic_synonyms="Enterobacter quasiroggenkampii Wu et al. 2020 (as basonym of P6A-ACT-04)",
    priority_basis="roggenkampii VL184 (2018) precedes quasiroggenkampii VL196 (2020)",
    icnp_rule_basis="37b; 42; 45; 35", type_strain="EN-117 = DSM 16690 = LMG 30172", type_genome="GCA_001729805.1",
    original_description_and_emendations="Sutton et al. 2018; no emendation recorded",
    documentary_gaps="Published type-pair isDDH 64.8% (<70%) must be addressed in the justification (not a Code requirement)",
    gap_disposition="manuscript_justification (P6A-AA-06)", author_action_ids="P6A-AA-06;P6A-Q02",
    formal_status_statement="Validly published name unchanged", evidence="sources/literature/wu2020-msystems-PMC7406230.html")

add(register_id="P6A-SPN-04", entry_class="no_act:circumscription_broadened",
    adopted_treatment="One E. chengduensis species (5 genomes) with rongchengensis as a subspecies",
    decision_source=D_S2, act_required="no (species name); subspecies acts P6A-ACT-05/P6A-TSP-03",
    proposed_name_and_citation="Enterobacter chengduensis Wu et al. 2019 (optional emend.; P6A-Q02)", rank="species", basonym="none",
    heterotypic_synonyms="Enterobacter rongchengensis He et al. 2025 (as basonym of P6A-ACT-05)",
    priority_basis="chengduensis VL188 (2019) precedes rongchengensis VL221 (2025)",
    icnp_rule_basis="37b; 42; 45; 35", type_strain="WCHECl-C4 = CCTCC AB 2017105 = CCUG 76941 = GDMCC 1.1216 = KCTC 52993", type_genome="GCA_001984825.2",
    original_description_and_emendations="Wu et al. 2019; no emendation recorded",
    documentary_gaps="Published type-pair isDDH 64.4% (<70%) must be addressed in the justification; only 5 genomes (3 + 2)",
    gap_disposition="manuscript_justification (P6A-AA-06)", author_action_ids="P6A-AA-06;P6A-Q02",
    formal_status_statement="Validly published name unchanged", evidence="sources/literature/he2024-spectrum-PMC11302248.html")

add(register_id="P6A-SPN-05", entry_class="no_act:heterotypic_synonymy",
    adopted_treatment="One E. mori species incl. E. quasimori and boundary genome GCA_025244905.1",
    decision_source=D_S2, act_required="no (synonymy is a taxonomic opinion)",
    proposed_name_and_citation="Enterobacter mori Zhu et al. 2011 emend. Wu and Zong 2020 (optional further emend.; P6A-Q02)", rank="species", basonym="none",
    heterotypic_synonyms="Enterobacter quasimori Wu et al. 2021 (VL198) — proposed later heterotypic synonym; Enterobacter tabaci Duan et al. 2016 (already a synonym, Wu and Zong 2020)",
    priority_basis="mori IJSEM 61:2769-2774 (2011) precedes tabaci (VL171, 2016) and quasimori (VL198, 2021)",
    icnp_rule_basis="23a; 24b(3); 42; 37b; 35", type_strain="R18-2 = CGMCC 1.10322 = DSM 26271 = LMG 25706", type_genome="GCA_000211415.1",
    original_description_and_emendations="Zhu et al. 2011 IJSEM 61:2769-2774; emend. Wu and Zong 2020 IJSEM 70:1055-1058",
    documentary_gaps="Published type-pair ANI 95.32% (Wu et al. 2020 Table 3, quasimori vs mori LMG 25706T) to be cited in justification",
    gap_disposition="manuscript_justification (P6A-AA-06)", author_action_ids="P6A-AA-06;P6A-Q02",
    formal_status_statement="Validly published name unchanged; synonymy is a proposal, published as a taxonomic opinion",
    evidence="sources/lpsn/raw/species__enterobacter-mori.html; sources/lpsn/raw/species__enterobacter-quasimori.html")

add(register_id="P6A-SPN-06", entry_class="no_act:heterotypic_synonymy",
    adopted_treatment="One broad E. asburiae species incl. E. dykesii; internal lineages informal; outside pair species-unassigned",
    decision_source=D_S1, act_required="no (synonymy is a taxonomic opinion)",
    proposed_name_and_citation="Enterobacter asburiae Brenner et al. 1988 emend. Hoffmann et al. 2005 (optional further emend.; P6A-Q02)", rank="species", basonym="none",
    heterotypic_synonyms="Enterobacter dykesii Cho et al. 2021 (VL200) — proposed later heterotypic synonym; Enterobacter muelleri Kämpfer et al. 2015 (already a synonym, Sutton et al. 2018; type GCA_900180435.1 in GTDB s__Enterobacter cloacae_M within this group)",
    priority_basis="asburiae VL25 (IJSB 38:220, 1988) precedes muelleri (2015) and dykesii (VL200, 2021)",
    icnp_rule_basis="23a; 24b(3); 42; 37b; 35", type_strain="1497-78 = ATCC 35953 = CIP 103358 = DSM 17506 = JCM 6051 = NCTC 12123", type_genome="GCA_001521715.1",
    original_description_and_emendations="Brenner et al. 1986 J Clin Microbiol (VL25); emend. Hoffmann et al. 2005",
    documentary_gaps="none for the name", gap_disposition="no_gap; emendation author choice", author_action_ids="P6A-AA-06;P6A-Q02",
    formal_status_statement="Validly published name unchanged; synonymy is a proposal", evidence="sources/lpsn/raw/species__enterobacter-asburiae.html; sources/lpsn/raw/species__enterobacter-dykesii.html")

add(register_id="P6A-SPN-07", entry_class="no_act:existing_names_match",
    adopted_treatment="E. cloacae one species with subspecies cloacae and dissolvens (supersedes earlier disagreement; ANIb follow-up complete)",
    decision_source=D_S1 + " ; " + D_CD, act_required="no (adopted treatment equals current LPSN correct names)",
    proposed_name_and_citation="Enterobacter cloacae (Jordan 1890) Hormaeche and Edwards 1960 (Approved Lists 1980)", rank="species", basonym="'Bacillus cloacae' Jordan 1890",
    homotypic_synonyms="", heterotypic_synonyms="",
    priority_basis="Approved Lists 1980", icnp_rule_basis="24a; 40d (subsp. cloacae created 2005)",
    type_strain="ATCC 13047 = DSM 30054 = …", type_genome="GCA_000025565.1",
    original_description_and_emendations="Hoffmann et al. 2005 (subspecies; VL104)", documentary_gaps="none", gap_disposition="no_gap",
    formal_status_statement="Validly published names unchanged", evidence="sources/lpsn/raw/species__enterobacter-cloacae.html")
add(register_id="P6A-SPN-07a", entry_class="no_act:existing_names_match", adopted_treatment="E. cloacae subsp. cloacae (41 study genomes)",
    decision_source=D_S1 + " ; " + D_CD, parent_entry="P6A-SPN-07", act_required="no",
    proposed_name_and_citation="Enterobacter cloacae subsp. cloacae (Jordan 1890) Hoffmann et al. 2005", rank="subspecies",
    icnp_rule_basis="40d", type_genome="GCA_000025565.1", documentary_gaps="none", gap_disposition="no_gap",
    formal_status_statement="Validly published (VL104)", evidence="sources/lpsn/raw/subspecies__enterobacter-cloacae-cloacae.html")
add(register_id="P6A-SPN-07b", entry_class="no_act:existing_names_match", adopted_treatment="E. cloacae subsp. dissolvens (4 study genomes + type)",
    decision_source=D_S1 + " ; " + D_CD, parent_entry="P6A-SPN-07", act_required="no",
    proposed_name_and_citation="Enterobacter cloacae subsp. dissolvens (Rosen 1922) Hoffmann et al. 2005", rank="subspecies",
    homotypic_synonyms="Enterobacter dissolvens (Rosen 1922) Brenner et al. 1988; Erwinia dissolvens (Rosen 1922) Burkholder 1948",
    icnp_rule_basis="—", type_strain="ATCC 23373 = DSM 16657 = …", type_genome="GCA_009660455.1",
    documentary_gaps="none", gap_disposition="no_gap", formal_status_statement="Validly published (VL104); current correct name in LPSN",
    evidence="sources/lpsn/raw/subspecies__enterobacter-cloacae-dissolvens.html")

# ---- Retained species: names unchanged ------------------------------------
RETAINED = {
    "bugandensis": (D_S1, "species name retained (final; masking qualifications)"),
    "nematophilus": (D_S1, "species name retained (final; masking qualifications)"),
    "pseudoroggenkampii": (D_S1 + " ; " + D_PSE, "separate from quasiroggenkampii; GCA_009755685.1 assigned"),
    "vonholyi": (D_VON, "species level only"),
    "sichuanensis": (D_S2, "separate species; two label corrections"),
    "chinensis": (D_S2, "separate species (type anchor only)"),
    "kobei": (D_KOB, "unchanged binomial; ten routine assignments approved"),
    "ludwigii": (D_LUD, "unchanged binomial; eight routine assignments approved"),
    "soli": (D_S2 + " ; " + D_SOLI, "retained; CP102 relatively divergent member"),
    "quasihormaechei": (D_QH, "unchanged binomial; three routine assignments approved"),
    "cancerogenus": ("unchanged current NCBI binomial (overlay basis)", "no author decision needed"),
    "chuandaensis": ("unchanged current NCBI binomial (overlay basis)", "no author decision needed"),
    "huaxiensis": ("frozen type identity (overlay basis)", "no author decision needed"),
    "wuhouensis": ("unchanged current NCBI binomial (overlay basis)", "no author decision needed"),
    "oligotrophicus": ("unchanged current NCBI binomial (overlay basis)", "no author decision needed"),
    "pasteurii": ("unchanged current NCBI binomial (overlay basis)", "no author decision needed"),
}
for ep, (src, note) in RETAINED.items():
    add(register_id=f"P6A-RET-{ep}", entry_class="no_act:name_unchanged", adopted_treatment=f"E. {ep}: {note}",
        decision_source=src, act_required="no", proposed_name_and_citation=f"LPSN:/species/enterobacter-{ep}", rank="species",
        documentary_gaps="none", gap_disposition="no_gap",
        formal_status_statement="Validly published, correct name in LPSN snapshot (see name-status.tsv)",
        evidence=f"sources/lpsn/raw/species__enterobacter-{ep}.html")

add(register_id="P6A-NVP-01", entry_class="no_act:name_not_validly_published",
    adopted_treatment="'Enterobacter adelaidei' group (4 study genomes; CZ-1 inside with boundary flag); unvalidated-name qualification retained",
    decision_source=D_S2, act_required="no (no act approved; decision retains the qualification)",
    proposed_name_and_citation="'Enterobacter adelaidei' Siderius et al. 2024 (in quotation marks; not validly published)", rank="species",
    type_strain="ECC3473 (nomenclatural type per LPSN; LPSN lists no culture-collection deposits)",
    original_description_and_emendations="Siderius et al. 2024 Microbiol Res 288:127867 (effective publication only)",
    documentary_gaps="Name has no standing; validation would require a Rule 27/30 compliant type deposit, which LPSN does not record",
    gap_disposition="no_act_by_decision: cite in quotation marks with 'not validly published'; any validation is outside this paper's approved acts (P6A-Q04)",
    author_action_ids="P6A-Q04", formal_status_statement="Not validly published (LPSN); must not be presented as a correct name",
    evidence="sources/lpsn/raw/species__enterobacter-adelaidei.html")

add(register_id="P6A-CTX-01", entry_class="no_act:context_taxon", adopted_treatment="RHB15-C17 as Lelliottia nimipressuralis context; GTDB crosswalk retained",
    decision_source=D_RES, act_required="no", proposed_name_and_citation="Lelliottia nimipressuralis (Carter 1945) Brady et al. 2013", rank="species",
    documentary_gaps="none", gap_disposition="no_gap", formal_status_statement="Validly published, correct name (LPSN)",
    evidence="sources/lpsn/raw/species__lelliottia-nimipressuralis.html")

# ---- Unnamed groups ----------------------------------------------------------
for rid, tr, src in [
        ("P6A-UNI-01", "cloacae_P pair (GCA_014169315.1, GCA_015708675.1): one unnamed lineage, genus qualified", D_RES),
        ("P6A-UNI-02", "RHBSTW-00994 (GCA_013782625.1): separate unnamed species-level Enterobacter lineage; genus-evidence flag", D_RHB),
        ("P6A-UNI-03", "ECC-175: separate unnamed Enterobacter lineage (Site decision 2 Oct)", MAP + " (residual genus section)")]:
    add(register_id=rid, entry_class="unnamed:included_no_naming", adopted_treatment=tr, decision_source=src, act_required="no (no name proposed)",
        proposed_name_and_citation="Enterobacter sp. / unnamed lineage", gap_disposition="no_act: report as unnamed; no type designated",
        formal_status_statement="No name; no Candidatus/SeqCode proposal")
add(register_id="P6A-UND-01", entry_class="unnamed:deferred_naming",
    adopted_treatment="Soli-related isolateB: separate unnamed lineage with QC/AF qualifications; formal naming deferred",
    decision_source=D_S2 + " ; " + D_SEQ, act_required="no (deferred)", proposed_name_and_citation="unnamed",
    documentary_gaps="Any later ICNP name needs a deposited type strain; SeqCode/Candidatus routes are a deferred exploration only",
    gap_disposition="deferred_by_decision: exclude from this paper's acts; may be listed as future work", author_action_ids="P6A-Q06",
    formal_status_statement="No name; naming deferred")
for rid, tr, src in [
        ("P6A-UNA-01", "ECC-249: species-unassigned (outside roggenkampii/pseudoroggenkampii)", D_S2 + " ; " + D_S1),
        ("P6A-UNA-02", "ECC-219 + strain 93: species-unassigned outside bugandensis and nematophilus", D_S1),
        ("P6A-UNA-03", "asburiae outside pair (GCA_019968765.1, GCA_025502525.1): species-unassigned", D_S1),
        ("P6A-UNA-04", "Peripheral Group B genome GCA_039734725.2: E. xiangfangensis, subspecies unassigned", D_S2 + " ; " + MAP)]:
    add(register_id=rid, entry_class="unassigned:not_a_taxon", adopted_treatment=tr, decision_source=src, act_required="no",
        proposed_name_and_citation="(no taxon)", gap_disposition="no_act", formal_status_statement="Genome-level assignment status only")
add(register_id="P6A-INF-01", entry_class="infrasubspecific:informal", parent_entry="P6A-ACT-03",
    adopted_treatment="GV01-GV25 numbered informal genomovars within subsp. hoffmannii (99.5% ANI, complete linkage)",
    decision_source=D_GA + " ; " + MAP, act_required="no (Rule 14a: below subspecies not regulated)",
    proposed_name_and_citation="GV01-GV25 (informal; not Latin, not names)", icnp_rule_basis="5d; 14a; Appendix 10",
    gap_disposition="no_act: must not be styled as names or 'genomovar nov.'", formal_status_statement="Unregulated informal labels")
add(register_id="P6A-INF-02", entry_class="infrasubspecific:informal", parent_entry="P6A-SPN-06",
    adopted_treatment="E. asburiae internal lineages (descriptive)", decision_source=D_S1, act_required="no (Rule 14a)",
    icnp_rule_basis="14a", gap_disposition="no_act", formal_status_statement="Unregulated descriptive labels")
add(register_id="P6A-OOS-01", entry_class="not_a_treatment",
    adopted_treatment="Rows without an adopted taxon: analysis-scope restrictions (study), declared out-of-genus context, outgroups, Task01 public exclusions",
    decision_source="overlay assignment_status/record_role", act_required="no", gap_disposition="outside_register_scope",
    formal_status_statement="No naming disposition implied")


def map_row(r: dict) -> str:
    sp, ssp, rule = r["project_species"], r["project_subspecies"], r["rule_id"]
    if sp == "Enterobacter xiangfangensis":
        return {"xiangfangensis": "P6A-TSP-01", "oharae": "P6A-ACT-01", "steigerwaltii": "P6A-ACT-02",
                "hoffmannii": "P6A-ACT-03", "": "P6A-UNA-04"}[ssp]
    if sp == "Enterobacter roggenkampii":
        return {"roggenkampii": "P6A-TSP-02", "quasiroggenkampii": "P6A-ACT-04"}[ssp]
    if sp == "Enterobacter chengduensis":
        return {"chengduensis": "P6A-TSP-03", "rongchengensis": "P6A-ACT-05"}[ssp]
    if sp == "Enterobacter cloacae":
        return {"cloacae": "P6A-SPN-07a", "dissolvens": "P6A-SPN-07b"}[ssp]
    fixed = {"Enterobacter hormaechei": "P6A-SPN-02", "Enterobacter mori": "P6A-SPN-05",
             "Enterobacter asburiae": "P6A-SPN-06", "Enterobacter adelaidei": "P6A-NVP-01",
             "Lelliottia nimipressuralis": "P6A-CTX-01"}
    if sp in fixed and not ssp:
        return fixed[sp]
    if sp.startswith("Enterobacter ") and sp.split()[1] in RETAINED and not ssp:
        return "P6A-RET-" + sp.split()[1]
    if sp == "Enterobacter sp.":
        return {"RHBSTW00994": "P6A-UNI-02", "ECC-175": "P6A-UNI-03", "isolateB": "P6A-UND-01", "ECC-249": "P6A-UNA-01",
                "bug-nemat-boundary": "P6A-UNA-02", "asburiae-outside-pair": "P6A-UNA-03"}[rule]
    if not sp and rule == "cloacae_P":
        return "P6A-UNI-01"
    if not sp and not ssp and r["assignment_status"] in ("excluded_or_restricted_analysis_scope", "retained_context_restriction"):
        return "P6A-OOS-01"
    raise KeyError(f"unmapped overlay row {r['accession']}: {sp!r} {ssp!r} {rule!r}")


TYPE_ACCS = {  # register entry -> (type genome accession, LPSN path)
    "P6A-ACT-01": ("GCA_001729705.1", "/subspecies/enterobacter-hormaechei-oharae"),
    "P6A-ACT-02": ("GCA_001729725.1", "/subspecies/enterobacter-hormaechei-steigerwaltii"),
    "P6A-ACT-03": ("GCA_001729745.1", "/subspecies/enterobacter-hormaechei-hoffmannii"),
    "P6A-ACT-04": ("GCA_001039365.2", "/species/enterobacter-quasiroggenkampii"),
    "P6A-ACT-05": ("GCA_034008275.1", "/species/enterobacter-rongchengensis"),
    "P6A-TSP-01": ("GCA_001729785.1", "/species/enterobacter-xiangfangensis"),
    "P6A-TSP-02": ("GCA_001729805.1", "/species/enterobacter-roggenkampii"),
    "P6A-TSP-03": ("GCA_001984825.2", "/species/enterobacter-chengduensis"),
}
SYN_TYPES = {  # heterotypic-synonym types expected inside the act's taxon
    "P6A-ACT-03": [("GCA_040096145.1", "/species/enterobacter-intestinihominis"),
                   ("GCA_052922205.1", "/subspecies/enterobacter-hoffmannii-nototheniae")],
}
COLL_TAXON = {  # check_collections.py taxon key -> register entry
    "E. xiangfangensis (type subspecies of the proposed arrangement)": "P6A-TSP-01",
    "E. hormaechei subsp. oharae (basonym of proposed comb. nov.)": "P6A-ACT-01",
    "E. hormaechei subsp. steigerwaltii (basonym of proposed comb. nov.)": "P6A-ACT-02",
    "E. hormaechei subsp. hoffmannii (basonym of proposed comb. nov.)": "P6A-ACT-03",
    "E. roggenkampii (type subspecies of the proposed arrangement)": "P6A-TSP-02",
    "E. quasiroggenkampii (basonym of proposed comb. nov.)": "P6A-ACT-04",
    "E. chengduensis (type subspecies of the proposed arrangement)": "P6A-TSP-03",
    "E. rongchengensis (basonym of proposed comb. nov.)": "P6A-ACT-05",
}
# Priority ordering keys: (year, volume, page or list number) from the LPSN/VL records.
PRIORITY = [
    ("xiangfangensis_A+B (species)", "Enterobacter xiangfangensis Gu et al. 2014", "/species/enterobacter-xiangfangensis", (2014, 64, 2650), "IJSEM 64:2650-2656 (original article)", "adopted"),
    ("xiangfangensis_A+B (species)", "Enterobacter siamensis Khunthongpan et al. 2014", "/species/enterobacter-siamensis", (2014, 64, 2927), "VL159, IJSEM 64:2927-2929", "no type assembly selected (no NCBI type-material candidate, 9 Oct screen; type-strain identity unresolved, Opinion 114), so its group is unknown; later page than xiangfangensis in the same volume (Rule 24b(4)), which settles competition with xiangfangensis only"),
    ("xiangfangensis_A+B (species)", "Enterobacter intestinihominis Hitch et al. 2025", "/species/enterobacter-intestinihominis", (2025, 75, 0), "VL225 (no. 28)", "type in A+B; later"),
    ("xiangfangensis_A+B (species)", "Enterobacter hoffmannii (Sutton et al. 2018) Sedláček et al. 2026", "/species/enterobacter-hoffmannii", (2026, 76, 0), "VL230 (no. 26)", "type in A+B; later"),
    ("Group A (subspecies)", "Enterobacter hormaechei subsp. hoffmannii Sutton et al. 2018", "/subspecies/enterobacter-hormaechei-hoffmannii", (2018, 68, 3379), "VL184, IJSEM 68:3379-3393", "adopted epithet"),
    ("Group A (subspecies)", "Enterobacter intestinihominis Hitch et al. 2025", "/species/enterobacter-intestinihominis", (2025, 75, 0), "VL225", "later (species epithet)"),
    ("Group A (subspecies)", "Enterobacter hoffmannii subsp. nototheniae Sedláček et al. 2026", "/subspecies/enterobacter-hoffmannii-nototheniae", (2026, 76, 0), "VL230", "later"),
    ("roggenkampii (species)", "Enterobacter roggenkampii Sutton et al. 2018", "/species/enterobacter-roggenkampii", (2018, 68, 3379), "VL184", "adopted"),
    ("roggenkampii (species)", "Enterobacter quasiroggenkampii Wu et al. 2020", "/species/enterobacter-quasiroggenkampii", (2020, 70, 5596), "VL196", "later"),
    ("chengduensis (species)", "Enterobacter chengduensis Wu et al. 2019", "/species/enterobacter-chengduensis", (2019, 69, 1844), "VL188", "adopted"),
    ("chengduensis (species)", "Enterobacter rongchengensis He et al. 2025", "/species/enterobacter-rongchengensis", (2025, 75, 0), "VL221 (no. 12)", "later"),
    ("mori (species)", "Enterobacter mori Zhu et al. 2011", "/species/enterobacter-mori", (2011, 61, 2769), "IJSEM 61:2769-2774 (original article)", "adopted"),
    ("mori (species)", "Enterobacter tabaci Duan et al. 2016", "/species/enterobacter-tabaci", (2016, 66, 3761), "VL171", "later; existing synonym"),
    ("mori (species)", "Enterobacter quasimori Wu et al. 2021", "/species/enterobacter-quasimori", (2021, 71, 0), "VL198", "later"),
    ("asburiae (species)", "Enterobacter asburiae Brenner et al. 1988", "/species/enterobacter-asburiae", (1988, 38, 220), "VL25, IJSB 38:220-222", "adopted"),
    ("asburiae (species)", "Enterobacter muelleri Kämpfer et al. 2015", "/species/enterobacter-muelleri", (2015, 65, 0), "IJSEM 65 (original article)", "later; existing synonym"),
    ("asburiae (species)", "Enterobacter dykesii Cho et al. 2021", "/species/enterobacter-dykesii", (2021, 71, 0), "VL200", "later"),
    ("cloacae (species)", "Enterobacter cloacae (Jordan 1890) Hormaeche and Edwards 1960", "/species/enterobacter-cloacae", (1980, 30, 225), "Approved Lists 1980", "adopted"),
    ("cloacae (species)", "Enterobacter dissolvens (Rosen 1922) Brenner et al. 1988", "/species/enterobacter-dissolvens", (1988, 38, 220), "VL25", "later at species rank; subspecies name current"),
]
EPITHETS = {"P6A-ACT-01": "oharae", "P6A-ACT-02": "steigerwaltii", "P6A-ACT-03": "hoffmannii",
            "P6A-ACT-04": "quasiroggenkampii", "P6A-ACT-05": "rongchengensis", "P6A-TSP-01": "xiangfangensis",
            "P6A-TSP-02": "roggenkampii", "P6A-TSP-03": "chengduensis"}


def write(path: Path, rows: list[dict], cols: list[str]):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlay", required=True, type=Path)
    ap.add_argument("--lpsn", required=True, type=Path)
    ap.add_argument("--collections", required=True, type=Path)
    ap.add_argument("--type-ncbi", required=True, type=Path)
    ap.add_argument("--type-gc", required=True, type=Path)
    ap.add_argument("--circ", required=True, type=Path)
    ap.add_argument("--placement", required=True, type=Path, help="checks/type-placement-<date>/type-placement.tsv")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    overlay = list(csv.DictReader(open(a.overlay), delimiter="\t"))
    lpsn = {r["path"]: r for r in csv.DictReader(open(a.lpsn), delimiter="\t")}
    ids = {e["register_id"] for e in E}
    assert len(ids) == len(E), "duplicate register ids"

    # treatment map
    groups = collections.OrderedDict()
    by_acc = {}
    for r in overlay:
        rid = map_row(r)
        assert rid in ids, rid
        by_acc[r["accession"]] = (rid, r)
        k = (r["project_species"], r["project_subspecies"], r["rule_id"], r["record_role"], rid)
        groups[k] = groups.get(k, 0) + 1
    tm = [{"project_species": k[0], "project_subspecies": k[1], "rule_id": k[2], "record_role": k[3],
           "register_id": k[4], "n_records": n} for k, n in sorted(groups.items())]
    write(a.out / "treatment-map.tsv", tm, list(tm[0].keys()))
    counts = collections.Counter()
    study = collections.Counter()
    for rid, r in by_acc.values():
        counts[rid] += 1
        study[rid] += r["record_role"] == "study"
    for e in E:
        e["n_overlay_records"] = counts.get(e["register_id"], 0)
        e["n_study_records"] = study.get(e["register_id"], 0)
    cols = COLS[:5] + ["n_overlay_records", "n_study_records"] + COLS[5:]
    # retained species: fill names from the LPSN snapshot
    for e in E:
        p = e["proposed_name_and_citation"]
        if p.startswith("LPSN:"):
            rec = lpsn[p[5:]]
            e["proposed_name_and_citation"] = rec["Name"]
            e["type_strain"] = rec["Type strain"]
    write(a.out / "formal-act-register.tsv", E, cols)

    # name status for every name the register cites
    paths = sorted({p for _, _, p, *_ in PRIORITY} | {v[1] for v in TYPE_ACCS.values()} |
                   {p for v in SYN_TYPES.values() for _, p in v} |
                   {f"/species/enterobacter-{ep}" for ep in RETAINED} |
                   {"/species/enterobacter-hormaechei", "/subspecies/enterobacter-hormaechei-hormaechei",
                    "/subspecies/enterobacter-hormaechei-xiangfangensis", "/subspecies/enterobacter-hoffmannii-hoffmannii",
                    "/species/enterobacter-hoffmannii-1", "/species/enterobacter-adelaidei",
                    "/subspecies/enterobacter-cloacae-cloacae", "/subspecies/enterobacter-cloacae-dissolvens",
                    "/species/lelliottia-nimipressuralis"})
    ns = []
    for p in paths:
        r = lpsn[p]
        ns.append({"lpsn_path": p, "name": r["Name"], "category": r["Category, proposed as"] or r["Category"],
                   "nomenclatural_status": r["Nomenclatural status"], "taxonomic_status": r["Taxonomic status"],
                   "correct_name_if_synonym": r["Correct name"], "type_strain": r["Type strain"] or r["Nomenclatural type"],
                   "valid_publication": r["IJSEM list"] or r["Valid publication"], "effective_publication": r["Effective publication"] or r["Original publication"],
                   "synonyms": r["Synonyms"], "lpsn_retrieved_utc": r["retrieved_utc"], "lpsn_raw_sha256": r["raw_sha256"]})
    write(a.out / "name-status.tsv", ns, list(ns[0].keys()))

    # type documentation
    coll = collections.defaultdict(list)
    for r in csv.DictReader(open(a.collections), delimiter="\t"):
        coll[COLL_TAXON[r["taxon"]]].append(r)
    ncbi = {r["accession"]: r for r in csv.DictReader(open(a.type_ncbi), delimiter="\t")}
    gc = {r["accession"]: r for r in csv.DictReader(open(a.type_gc), delimiter="\t")}
    td = []
    for rid, (acc, path) in TYPE_ACCS.items():
        cs = coll[rid]
        listed = sorted({f"{c['designation']} ({c['country']})" for c in cs if c["catalogue_result"] == "listed_type_strain"})
        unchecked = sorted({f"{c['designation']} ({c['country']})" for c in cs if c["catalogue_result"].startswith("not_checked")})
        countries = sorted({c["country"] for c in cs if c["country"]})
        orow = by_acc.get(acc, (None, {}))[1]
        td.append({"register_id": rid, "lpsn_name": lpsn[path]["Name"], "type_strain_lpsn": lpsn[path]["Type strain"],
                   "designated_collection_countries": "; ".join(countries), "n_designated_countries": len(countries),
                   "catalogue_listed_8oct2026": "; ".join(listed), "n_listed_countries": len({c['country'] for c in cs if c['catalogue_result'] == 'listed_type_strain'}),
                   "not_machine_checkable": "; ".join(unchecked),
                   "type_genome": acc, "ncbi_type_material": ncbi[acc]["type_material"], "ncbi_strain": ncbi[acc]["strain"],
                   "ncbi_culture_collection": ncbi[acc]["culture_collection"], "assembly_level": ncbi[acc]["assembly_level"],
                   "genome_size_bp": gc[acc]["total_length"], "gc_mol_percent": gc[acc]["gc_mol_percent_acgt"],
                   "fasta_sha256": gc[acc]["fasta_sha256"],
                   "overlay_placement": f"{orow.get('project_species', '')} | {orow.get('project_subspecies', '')} | {orow.get('record_role', '')}",
                   "overlay_register_id": by_acc.get(acc, ("",))[0]})
    write(a.out / "type-documentation.tsv", td, list(td[0].keys()))

    # priority
    pr = []
    for grp, name, path, key, cite, note in PRIORITY:
        pr.append({"adopted_taxon_rank": grp, "name": name, "lpsn_valid_publication": lpsn[path]["IJSEM list"] or lpsn[path]["Valid publication"],
                   "valid_publication_short": cite, "order_key": "%d-%03d-%05d" % key, "note": note})
    write(a.out / "priority-check.tsv", pr, list(pr[0].keys()))

    # metadata competitor screen: coverage and missingness of every audited name
    reason = {
        "type_assembly_selected": "screened",
        "no_ncbi_type_material_candidate": "not screened: no NCBI type-material assembly returned for the listed or current name",
        "ncbi_type_candidates_failed_alias_match": ("not screened: NCBI type-material candidates exist, but none carries an LPSN "
                                                    "type-strain designation; no type genome selected (absence not established)"),
    }
    cov = []
    for r in csv.DictReader(open(a.placement), delimiter="\t"):
        cov.append({"lpsn_name": r["lpsn_name"], "lpsn_correct_name": r["lpsn_correct_name"], "name_category": r["name_category"],
                    "coverage_status": r["coverage_status"], "screen_result": reason[r["coverage_status"]],
                    "type_assembly": r["type_assembly"], "selection": r["selection"],
                    "n_ncbi_type_candidates": r["n_ncbi_type_candidates"], "ncbi_query_completeness": r["ncbi_query_completeness"],
                    "gtdb_species": r["gtdb_species"], "gtdb_release": r["gtdb_release"], "gtdb_source": r["gtdb_source"],
                    "falls_in_merged_group_gtdb_species": r["falls_in_merged_group_gtdb_species"],
                    "unselected_candidates_gtdb_r232_species": r["unmatched_candidates_gtdb_r232_species"]})
    write(a.out / "priority-screen-coverage.tsv", cov, list(cov[0].keys()))

    # epithet uniqueness within Enterobacter (Rule 12b/13c)
    def own_epithet(p: str) -> str:
        # species page: its specific epithet; subspecies page: its subspecific epithet (the species
        # epithet inside a subspecies name belongs to the species' type, not the subspecies')
        stem = p.split("/")[-1]
        if stem.endswith("-1"):
            stem = stem[:-2]
        return stem.split("-")[-1]

    eu = []
    for rid, ep in EPITHETS.items():
        hits = [r for p, r in lpsn.items() if p.startswith(("/species/enterobacter-", "/subspecies/enterobacter-"))
                and own_epithet(p) == ep]
        types = sorted({r["Type strain"] or r["Nomenclatural type"] for r in hits})
        eu.append({"register_id": rid, "epithet": ep, "n_enterobacter_names_with_epithet": len(hits),
                   "names": " | ".join(sorted(r["Name"] for r in hits)), "distinct_type_strain_lists": len(types),
                   "type_strain_lists": " || ".join(types)})
    write(a.out / "epithet-uniqueness.tsv", eu, list(eu[0].keys()))
    print(f"register_entries={len(E)} overlay_rows={len(overlay)} mapped={len(by_acc)} treatment_groups={len(tm)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
