"""Synthetic tests for label_transfer.py: every clade reading, trigger and sealing check (design section 7)."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import label_transfer as lt  # noqa: E402

K, C = lt.OUTGROUPS["klebsiella_aerogenes"], lt.OUTGROUPS["cronobacter_sakazakii"]
A, B = "Enterobacter alpha", "Enterobacter beta"
LABELS = {K: None, C: None, "a1": (A, "x"), "a2": (A, "x"), "a3": (A, "y"), "b1": (B, ""), "b2": (B, ""), "u1": None}
# v1 sits with a3 (A subsp. y); v2 with the unassigned u1; v3 with b1+b2; v4 joins A and B; v5 joins a2+a3 (mixed subspecies).
TREE = (f"({K}:1,{C}:1,(((a1:1,(a2:1,v5:1):1):1,(a3:1,v1:1):1):1,v4:1,((b1:1,b2:1,v3:1):1,(u1:1,v2:1):1):1):1);")


def readings(tree=TREE, rooting=K):
    return lt.read_clades(lt.rooted(tree, rooting), LABELS, ["v1", "v2", "v3", "v4", "v5"])


def test_clade_readings():
    r = readings()
    assert (r["v1"]["species"], r["v1"]["subspecies"], r["v1"]["clade_size"]) == (A, "y", 1)
    assert r["v2"]["species"] == lt.UNASSIGNED
    assert (r["v3"]["species"], r["v3"]["subspecies"]) == (B, "")
    assert r["v4"]["species"] == lt.UNASSIGNED  # its clade reaches u1
    assert (r["v5"]["species"], r["v5"]["subspecies"]) == (A, "x")


def test_mixed_and_ambiguous_readings():
    tree = f"({K}:1,{C}:1,((a1:1,a3:1,v1:1):1,(a2:1,b1:1,v2:1):1,(b2:1,u1:1):1):1);"
    r = lt.read_clades(lt.rooted(tree, K), LABELS, ["v1", "v2"])
    assert (r["v1"]["species"], r["v1"]["subspecies"]) == (A, lt.AMBIGUOUS)
    assert r["v2"]["species"] == lt.MIXED


def test_other_validation_genomes_are_ignored():
    tree = f"({K}:1,{C}:1,((v1:1,v2:1):1,a3:1):1,(a1:1,a2:1,b1:1,b2:1,u1:1):1);"
    r = lt.read_clades(lt.rooted(tree, K), LABELS, ["v1", "v2"])
    assert r["v1"]["species"] == A and r["v1"]["clade_size"] == 1


def test_both_rootings_agree_inside_the_ingroup():
    assert readings(rooting=K) == readings(rooting=C)


def test_changed_labels():
    paper = lt.rooted(f"({K}:1,{C}:1,(((a1:1,a2:1):1,a3:1):1,(b1:1,b2:1):1,u1:1):1);", K)
    same = lt.rooted(TREE, K)
    assert lt.changed_labels(paper, same, LABELS) == set()
    moved = lt.rooted(f"({K}:1,{C}:1,((a1:1,a2:1):1,(a3:1,(b1:1,b2:1):1):1,u1:1,v1:1):1);", K)
    changed = lt.changed_labels(paper, moved, LABELS)
    assert ("species", A) in changed and ("subspecies", A, "y") not in changed  # single-tip labels never change
    assert ("species", B) not in changed


def comparators():
    rows = [dict(accession="r1", label_kind="named", project_species=A, project_subspecies="x", is_package_reference="True",
                 whole_fasta="/c/r1.fna"),
            dict(accession="c2", label_kind="named", project_species=A, project_subspecies="y", is_package_reference="False",
                 whole_fasta="/c/c2.fna"),
            dict(accession="c3", label_kind="named", project_species=B, project_subspecies="", is_package_reference="False",
                 whole_fasta="/c/c3.fna"),
            dict(accession="u9", label_kind="unassigned_or_context", project_species="", project_subspecies="",
                 is_package_reference="False", whole_fasta="/c/u9.fna")]
    return {r["whole_fasta"]: r for r in rows}


def hit(ref, ani, afr=90, afq=90, query="/q/g.fna"):
    return {"Ref_file": f"/c/{ref}.fna", "Query_file": query, "ANI": str(ani), "Align_fraction_ref": str(afr),
            "Align_fraction_query": str(afq)}


def test_ani_evidence():
    ev = lt.ani_evidence([hit("r1", 99.95, 95, 92), hit("c2", 97), hit("c3", 95.0, 80, 79.9), hit("u9", 95.5)],
                         comparators(), "/q/g.fna")
    assert ev["S"] == [A] and ev["U"] == ["u9"] and ev["near_clone"] is True
    assert [h["accession"] for h in ev["hits"]] == ["r1", "c2", "u9", "c3"]
    top = lt.top_within_species(ev, A)
    assert (top["subspecies"], top["is_reference"], top["tie"]) == ("x", True, False)
    assert lt.top_within_species(ev, "Enterobacter gamma")["subspecies"] is None


def test_top_comparator_tie_and_non_reference():
    ev = lt.ani_evidence([hit("r1", 97.5), hit("c2", 97.5)], comparators(), "/q/g.fna")
    top = lt.top_within_species(ev, A)
    assert top["tie"] and top["subspecies"] is None and not top["is_reference"]
    ev = lt.ani_evidence([hit("r1", 97.0), hit("c2", 97.5)], comparators(), "/q/g.fna")
    assert lt.top_within_species(ev, A)["is_reference"] is False


@pytest.mark.parametrize("hits, message", [
    ([hit("zz", 99)], "outside the comparator manifest"),
    ([hit("r1", 99, query="/q/other.fna")], "is not the sampled FASTA"),
    ([hit("r1", 99), hit("r1", 98)], "duplicate comparator"),
    ([hit("r1", 101)], "outside 0-100"),
])
def test_ani_evidence_refuses(hits, message):
    with pytest.raises(ValueError, match=message):
        lt.ani_evidence(hits, comparators(), "/q/g.fna")


def ev(S=(A,), U=(), bac=(A, "x"), vert=(A, "x"), alt=None, markers=None, vok=True, changed=(), top=None):
    def reading(sp, ss):
        return dict(species=sp, subspecies=ss, clade_size=3)
    alt = alt or {}
    rd = {m: {"klebsiella_aerogenes": reading(*(bac if m == "bac120" else vert)),
              "cronobacter_sakazakii": reading(*alt.get(m, bac if m == "bac120" else vert))} for m in lt.METHODS}
    top = top or dict(subspecies="x", tie=False, is_reference=True, accession="r1", ani=99.0)
    return dict(S=list(S), U=list(U), readings=rd, verticall_ok=vok, changed=set(changed),
                markers=markers or dict(unique_markers=118, translation_table="11", aligned=True), top=lambda sp: top)


def test_automatic_single_and_tree_resolved():
    d = lt.decide(ev())
    assert d["species_auto"] and d["species_label_source"] == "auto_single_species"
    assert (d["subspecies_label"], d["subspecies_label_source"], d["top_comparator_is_reference"]) == ("x", "automatic", True)
    d = lt.decide(ev(S=(A, B)))
    assert d["species_label"] == A and d["species_label_source"] == "auto_tree_resolved"
    d = lt.decide(ev(S=(B,), bac=(B, ""), vert=(B, "")))
    assert d["species_label"] == B and d["subspecies_label"] == "" and d["subspecies_reasons"] == []


@pytest.mark.parametrize("kwargs, trigger", [
    (dict(S=()), "S_empty"),
    (dict(U=("u9",)), "unassigned_comparator_passes"),
    (dict(bac=(lt.MIXED, "")), "bac120_clade_mixed"),
    (dict(vert=(lt.UNASSIGNED, "")), "verticall_clade_unassigned"),
    (dict(S=(A, B), vert=(B, "")), "trees_disagree"),
    (dict(S=(B,)), "T_not_in_S"),
    (dict(markers=dict(unique_markers=0, translation_table="", aligned=False)), "gtdbtk_not_identified_or_aligned"),
    (dict(markers=dict(unique_markers=99, translation_table="11", aligned=True)), "low_markers"),
    (dict(markers=dict(unique_markers=118, translation_table="4", aligned=True)), "translation_table_not_11"),
    (dict(vok=False), "verticall_pair_missing"),
    (dict(changed=[("species", A)]), "changed_species_clade"),
    (dict(S=(A, B), changed=[("species", B)]), "changed_species_clade"),
    (dict(alt={"bac120": (B, "")}), "bac120_rooting_sensitive"),
    (dict(alt={"verticall": (lt.MIXED, "")}), "verticall_rooting_sensitive"),
])
def test_each_trigger_sends_the_genome_to_michael(kwargs, trigger):
    d = lt.decide(ev(**kwargs))
    assert trigger in d["triggers"] and not d["species_auto"] and d["species_label"] == ""


@pytest.mark.parametrize("kwargs, reason", [
    (dict(vert=(A, "y")), "trees_subspecies_disagree_or_ambiguous"),
    (dict(bac=(A, lt.AMBIGUOUS)), "trees_subspecies_disagree_or_ambiguous"),
    (dict(alt={"verticall": (A, "y")}), "verticall_subspecies_rooting_sensitive"),
    (dict(changed=[("subspecies", A, "x")]), "changed_subspecies_clade"),
    (dict(top=dict(subspecies=None, tie=True, is_reference=False, accession="r1;c2", ani=97.5)), "top_comparator_subspecies_tie"),
    (dict(top=dict(subspecies="", tie=False, is_reference=False, accession="c5", ani=97.5)), "top_comparator_without_subspecies"),
    (dict(top=dict(subspecies="y", tie=False, is_reference=False, accession="c2", ani=97.5)), "top_comparator_subspecies_differs"),
])
def test_subspecies_only_edges_keep_the_automatic_species(kwargs, reason):
    d = lt.decide(ev(**kwargs))
    assert d["species_auto"] and d["species_label"] == A
    assert reason in d["subspecies_reasons"] and d["subspecies_label"] == ""


SPECIES_IDS = {A: "SP-A", B: "SP-B"}
SUBSPECIES_IDS = {f"{A} subsp. x": "SS-AX", f"{A} subsp. y": "SS-AY"}


def preseal_row(code, acc, to_michael="", species=A, sub="x", sp_source="auto_single_species", ss_source="automatic",
                S=A, top_ref="True", top=None):
    top = top if top is not None else {A: dict(subspecies="x", tie=False, is_reference=top_ref == "True", accession="r1", ani=98.0)}
    return dict(coded_id=code, accession=acc, query_sha256="h" + acc, S=S, near_clone="False", to_michael=to_michael,
                species_label=species, subspecies_label=sub, species_label_source=sp_source,
                subspecies_label_source=ss_source, top_comparator_is_reference=top_ref, top_by_species=json.dumps(top))


def sample_rows(*accs):
    return {a: dict(assembly_unit=a, component_id=f"PRV-{a}", sample_group="Enterobacter alpha", qc_stratum="high_quality")
            for a in accs}


def michael_label(code, decision="species", species=A, sub="", confidence="high", rationale="clear clade"):
    return {code: dict(coded_id=code, decision=decision, species=species, subspecies=sub, confidence=confidence,
                       rationale=rationale)}


def test_seal_merges_and_writes_truth():
    edge = dict(to_michael="species", species="", sub="", sp_source="", ss_source="")
    pre = [preseal_row("E1", "g1"),
           preseal_row("E2", "g2", S=f"{A};{B}", **edge),
           preseal_row("E3", "g3", to_michael="subspecies", sub="", ss_source=""),
           preseal_row("E4", "g4", S="", **edge),
           preseal_row("E5", "g5", S=B, **edge),
           preseal_row("E6", "g6", S="", **edge),
           preseal_row("E7", "g7", **edge),
           preseal_row("E8", "g8", **edge),
           preseal_row("E9", "g9", to_michael="subspecies", sub="", ss_source="")]
    labels = (michael_label("E2", species=B) | michael_label("E3", sub="y") | michael_label("E4", decision="outside_genus", species="")
              | michael_label("E5", decision="unresolved", species="")
              | michael_label("E6", decision="species_without_reference", species="Enterobacter novus")
              | michael_label("E7", sub="x") | michael_label("E8", sub="undetermined") | michael_label("E9", sub="none"))
    sealed, truth = lt.seal_rows(pre, labels, sample_rows(*[f"g{i}" for i in range(1, 10)]), SPECIES_IDS, SUBSPECIES_IDS)
    t = {r["genome_id"]: r for r in truth}
    assert (t["g1"]["truth_kind"], t["g1"]["label_source"], t["g1"]["truth_subspecies_id"], t["g1"]["subspecies_structural"]) == \
        ("named", "auto_single_species", "SS-AX", "True")
    assert (t["g2"]["truth_species_id"], t["g2"]["label_source"], t["g2"]["ani_species_passing"]) == ("SP-B", "michael", "SP-A;SP-B")
    assert (t["g3"]["truth_subspecies_id"], t["g3"]["subspecies_structural"]) == ("SS-AY", "False")
    assert (t["g4"]["truth_kind"], t["g4"]["label_source"]) == ("outside_genus", "michael")
    assert t["g5"]["truth_kind"] == "unresolved"
    assert (t["g6"]["truth_kind"], t["g6"]["truth_species_id"]) == ("no_reference", "")
    # Michael's subspecies equal to an untied reference top comparator cannot be wrong either (R4-m2).
    assert (t["g7"]["truth_subspecies_id"], t["g7"]["subspecies_structural"], t["g7"]["label_source"]) == ("SS-AX", "True", "michael")
    assert (t["g8"]["truth_subspecies_id"], t["g8"]["subspecies_structural"]) == ("undetermined", "False")
    assert (t["g9"]["truth_subspecies_id"], t["g9"]["subspecies_structural"]) == ("", "False")
    assert all(r["query_sha256"] == "h" + r["genome_id"] and r["block"] == f"PRV-{r['genome_id']}" for r in truth)
    s = {r["coded_id"]: r for r in sealed}
    assert s["E3"]["subspecies_label_source"] == "michael" and s["E2"]["michael_confidence"] == "high"
    assert (s["E8"]["subspecies_label"], s["E9"]["subspecies_label"]) == ("undetermined", "none")


@pytest.mark.parametrize("top, structural", [
    ({A: dict(subspecies="x", tie=False, is_reference=True, accession="r1", ani=98.0)}, "True"),
    ({A: dict(subspecies="x", tie=False, is_reference=False, accession="c1", ani=98.0)}, "False"),
    ({A: dict(subspecies=None, tie=True, is_reference=False, accession="r1;c2", ani=98.0)}, "False"),
    ({A: dict(subspecies="y", tie=False, is_reference=True, accession="r2", ani=98.0)}, "False"),
    ({}, "False"),
])
def test_subspecies_structural_whoever_labels(top, structural):
    pre = [preseal_row("E7", "g7", to_michael="species", species="", sub="", sp_source="", ss_source="", top=top)]
    _, truth = lt.seal_rows(pre, michael_label("E7", sub="x"), sample_rows("g7"), SPECIES_IDS, SUBSPECIES_IDS)
    assert truth[0]["subspecies_structural"] == structural


@pytest.mark.parametrize("labels, message", [
    ({}, "cover exactly the edge cases"),
    (michael_label("E2") | michael_label("E1"), "cover exactly the edge cases"),
    (michael_label("E2") | michael_label("E9"), "unknown coded IDs"),
    (michael_label("E2", rationale=" "), "incomplete or invalid"),
    (michael_label("E2", confidence="sure"), "incomplete or invalid"),
    (michael_label("E2", decision="maybe"), "incomplete or invalid"),
    (michael_label("E2", species=""), "needs a species"),
    (michael_label("E2", sub="zz"), "unknown subspecies"),
    (michael_label("E2", sub=""), "unknown subspecies"),
    (michael_label("E2", species=B, sub="x"), "has no subspecies in the package"),
    (michael_label("E2", species="Enterobacter novus"), "not a package species"),
    (michael_label("E2", species="Enterobacter alpah", sub="x"), "not a package species"),
    (michael_label("E2", decision="species_without_reference", species=A), "a species name the package lacks"),
    (michael_label("E2", decision="species_without_reference", species=""), "a species name the package lacks"),
    (michael_label("E2", decision="species_without_reference", species="Enterobacter novus", sub="x"), "takes no subspecies"),
    (michael_label("E2", decision="species_without_reference", species=A + " "), "within 2 edits"),
    (michael_label("E2", decision="species_without_reference", species="enterobacter ALPHA"), "within 2 edits"),
    (michael_label("E2", decision="species_without_reference", species="Enterobacter alpah"), "within 2 edits"),
    (michael_label("E2", decision="species_without_reference", species="Enterobacter alph"), "within 2 edits"),
    (michael_label("E2", decision="unresolved", species="", sub="x"), "takes no species or subspecies"),
    (michael_label("E2", decision="outside_genus"), "takes no species or subspecies"),
])
def test_seal_refuses(labels, message):
    pre = [preseal_row("E1", "g1"),
           preseal_row("E2", "g2", to_michael="species", species="", sub="", sp_source="", ss_source="")]
    with pytest.raises(ValueError, match=message):
        lt.seal_rows(pre, labels, sample_rows("g1", "g2"), SPECIES_IDS, SUBSPECIES_IDS)


def test_species_without_reference_name_is_normalised():
    pre = [preseal_row("E6", "g6", to_michael="species", species="", sub="", sp_source="", ss_source="")]
    sealed, truth = lt.seal_rows(pre, michael_label("E6", decision="species_without_reference", species=" Enterobacter  novus "),
                                 sample_rows("g6"), SPECIES_IDS, SUBSPECIES_IDS)
    assert sealed[0]["species_label"] == "Enterobacter novus" and truth[0]["truth_kind"] == "no_reference"


def test_sealed_truth_scores_with_the_frozen_scorer(tmp_path):
    """seal -> score_validation.main, including Michael's undetermined subspecies on an automatic species (R5-M1)."""
    import score_validation as sv
    edge = dict(to_michael="species", species="", sub="", sp_source="", ss_source="")
    pre = [preseal_row("E1", "g1"),
           preseal_row("E3", "g3", to_michael="subspecies", sub="", ss_source=""),
           preseal_row("E6", "g6", S="", **edge),
           preseal_row("E7", "g7", **edge),
           preseal_row("E8", "g8", **edge)]
    labels = (michael_label("E3", sub="undetermined") | michael_label("E6", decision="species_without_reference",
                                                                      species="Enterobacter novus")
              | michael_label("E7", sub="y") | michael_label("E8", sub="undetermined"))
    _, truth = lt.seal_rows(pre, labels, sample_rows("g1", "g3", "g6", "g7", "g8"), SPECIES_IDS, SUBSPECIES_IDS)
    write_tsv(tmp_path / "truth.tsv", truth)
    manifest = []
    for t in truth:
        report = tmp_path / f"{t['genome_id']}.report.tsv"
        write_tsv(report, [dict(experimental_eligibility_status="completed", experimental_eligibility_reasons="",
                                hypothetical_species_id="SP-A", hypothetical_subspecies_id="SS-AX", near_boundary_species_id="",
                                query_sha256=t["query_sha256"])])
        manifest.append(dict(genome_id=t["genome_id"], report_tsv=str(report),
                             report_sha256=hashlib.sha256(report.read_bytes()).hexdigest()))
    write_tsv(tmp_path / "results.tsv", manifest)
    sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results.tsv"), "--group", "label_source",
             "--output", str(tmp_path / "scored")])
    scored = {r["genome_id"]: r for r in csv.DictReader((tmp_path / "scored" / "scored-genomes.tsv").open(), delimiter="\t")}
    assert (scored["g1"]["species_outcome"], scored["g1"]["subspecies_outcome"], scored["g1"]["subspecies_structural"]) == \
        ("correct", "correct", "True")
    assert (scored["g3"]["label_source"], scored["g3"]["subspecies_outcome"]) == ("auto_single_species", "truth_undetermined")
    assert scored["g6"]["species_outcome"] == "wrong"
    assert scored["g7"]["subspecies_outcome"] == "wrong_subspecies"
    assert scored["g8"]["subspecies_outcome"] == "truth_undetermined"


def test_seal_refuses_an_automatic_species_outside_the_package():
    with pytest.raises(ValueError, match="automatic species"):
        lt.seal_rows([preseal_row("E1", "g1", species="Enterobacter novus")], {}, sample_rows("g1"), SPECIES_IDS, SUBSPECIES_IDS)


def test_verticall_ok(tmp_path):
    def row(g, q="12", t="12", zero="0", dup="0", expected="12"):
        return dict(assembly_id=g, directed_rows_as_query=q, directed_rows_as_target=t, expected_directed_pairs_each_direction=expected,
                    zero_result_pairs=zero, pairs_with_secondary_results="0", pairs_with_duplicate_rows=dup)
    write_tsv(tmp_path / "review.tsv", [row("ok"), row("secondary", q="15", t="14"), row("short", t="11"),
                                        row("missing", q="11", zero="1"), row("duplicate", q="13", dup="1")])
    ok = lt._verticall_ok(tmp_path / "review.tsv", ["ok", "secondary", "short", "missing", "duplicate", "absent"], 12)
    assert ok == dict(ok=True, secondary=True, short=False, missing=False, duplicate=False, absent=False)
    with pytest.raises(ValueError, match="expects 12 pairs per genome, not 13"):
        lt._verticall_ok(tmp_path / "review.tsv", ["ok"], 13)


def test_subspecies_only_edge_cannot_change_species():
    pre = [preseal_row("E3", "g3", to_michael="subspecies", sub="", ss_source="")]
    with pytest.raises(ValueError, match="keeps the automatic species"):
        lt.seal_rows(pre, michael_label("E3", species=B), sample_rows("g3"), SPECIES_IDS, SUBSPECIES_IDS)


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_prepare_end_to_end(tmp_path, monkeypatch):
    comp = [dict(accession=a, label_kind="named" if lab else "unassigned_or_context", project_species=lab[0] if lab else "",
                 project_subspecies=lab[1] if lab else "", is_package_reference=str(a == "a3"), whole_fasta=f"/c/{a}.fna",
                 whole_sha256="x") for a, lab in LABELS.items() if a not in (K, C)]
    write_tsv(tmp_path / "comp.tsv", comp)
    write_tsv(tmp_path / "tips.tsv", [dict(accession=a, tip_label_kind="named" if lab else "unassigned_or_context",
                                           project_species=lab[0] if lab else "", project_subspecies=lab[1] if lab else "")
                                      for a, lab in LABELS.items()])
    vs = ["v1", "v2", "v3", "v4", "v5"]
    write_tsv(tmp_path / "sample.tsv", [dict(assembly_unit=v, component_id=f"PRV-{v}", sample_group="x", qc_stratum="high_quality")
                                        for v in vs])
    monkeypatch.setattr(lt, "COMPARATORS", tmp_path / "comp.tsv")
    monkeypatch.setattr(lt, "TREE_TIPS", tmp_path / "tips.tsv")
    monkeypatch.setattr(lt, "SAMPLE", tmp_path / "sample.tsv")
    manifest = []
    for v, ref in zip(vs, ["a3", "u1", "b1", "a1", "a2"]):
        fasta = tmp_path / f"{v}.fna"
        fasta.write_text(f">{v}\nACGT\n")
        hits = tmp_path / f"{v}.hits.tsv"
        write_tsv(hits, [hit(ref, 98.5, query=str(fasta)), hit("b2", 85, 40, 40, query=str(fasta))])
        manifest.append(dict(accession=v, query_fasta=str(fasta), query_sha256=hashlib.sha256(fasta.read_bytes()).hexdigest(),
                             hits_tsv=str(hits), hits_sha256=hashlib.sha256(hits.read_bytes()).hexdigest()))
    write_tsv(tmp_path / "manifest.tsv", manifest)
    (tmp_path / "rebuilt.nwk").write_text(TREE)
    paper = f"({K}:1,{C}:1,((a1:1,a2:1):1,a3:1):1,((b1:1,b2:1):1,u1:1):1);"
    (tmp_path / "paper.nwk").write_text(paper)
    write_tsv(tmp_path / "markers.tsv", [dict(name=f"full_{v}", number_unique_genes="118" if v != "v3" else "60") for v in vs])
    (tmp_path / "tt.tsv").write_text("".join(f"full_{v}\t11\n" for v in vs))
    (tmp_path / "aligned.txt").write_text("".join(f"full_{v}\n" for v in vs))
    write_tsv(tmp_path / "review.tsv", [dict(assembly_id=v, zero_result_pairs="0", pairs_with_duplicate_rows="0",
                                             directed_rows_as_query="14" if v == "v1" else "12",
                                             directed_rows_as_target={"v1": "13", "v5": "11"}.get(v, "12"),
                                             expected_directed_pairs_each_direction="12") for v in vs])
    out = tmp_path / "out"
    lt.main(["prepare", "--skani-manifest", str(tmp_path / "manifest.tsv"), "--bac120-tree", str(tmp_path / "rebuilt.nwk"),
             "--verticall-tree", str(tmp_path / "rebuilt.nwk"), "--paper-bac120-tree", str(tmp_path / "paper.nwk"),
             "--paper-verticall-tree", str(tmp_path / "paper.nwk"), "--marker-summary", str(tmp_path / "markers.tsv"),
             "--translation-table", str(tmp_path / "tt.tsv"), "--aligned-ids", str(tmp_path / "aligned.txt"),
             "--verticall-assembly-review", str(tmp_path / "review.tsv"), "--gtdbtk-prefix", "full_", "--output", str(out)])
    rows = {r["accession"]: r for r in csv.DictReader((out / "preseal-labels.tsv").open(), delimiter="\t")}
    assert rows["v1"]["species_label_source"] == "auto_single_species" and rows["v1"]["subspecies_label"] == "y"
    assert rows["v1"]["top_comparator_is_reference"] == "True" and rows["v1"]["to_michael"] == ""
    assert rows["v1"]["verticall_ok"] == "True"  # secondary-result rows above the expected count are complete
    top = json.loads(rows["v1"]["top_by_species"])
    assert set(top) == {A, B} and top[A]["is_reference"] and top[A]["subspecies"] == "y" and not top[B]["is_reference"]
    assert "unassigned_comparator_passes" in rows["v2"]["triggers"] and rows["v2"]["to_michael"] == "species"
    assert "low_markers" in rows["v3"]["triggers"]
    assert "verticall_clade_unassigned" in rows["v4"]["triggers"]
    assert "verticall_pair_missing" in rows["v5"]["triggers"]
    packets = sorted(p.stem for p in (out / "packets").iterdir())
    assert packets == sorted(r["coded_id"] for r in rows.values() if r["to_michael"])
    assert lt.coded_id("v1") not in packets and "v2" not in (out / "packets" / f"{lt.coded_id('v2')}.md").read_text()
    coded = (out / "trees" / "bac120.klebsiella_aerogenes.coded.nwk").read_text()
    assert lt.coded_id("v1") in coded and "v1:" not in coded
    summary = json.loads((out / "summary.json").read_text())
    assert summary["genomes"] == 5 and summary["changed_labels"] == []
    with pytest.raises(SystemExit):
        lt.main(["prepare", "--skani-manifest", str(tmp_path / "manifest.tsv"), "--bac120-tree", str(tmp_path / "rebuilt.nwk"),
                 "--verticall-tree", str(tmp_path / "rebuilt.nwk"), "--paper-bac120-tree", str(tmp_path / "paper.nwk"),
                 "--paper-verticall-tree", str(tmp_path / "paper.nwk"), "--marker-summary", str(tmp_path / "markers.tsv"),
                 "--translation-table", str(tmp_path / "tt.tsv"), "--aligned-ids", str(tmp_path / "aligned.txt"),
                 "--verticall-assembly-review", str(tmp_path / "review.tsv"), "--output", str(out)])
