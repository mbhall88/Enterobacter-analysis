"""Synthetic-record tests for score_validation.py; every outcome category in design section 3."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_validation as sv  # noqa: E402

SP, SP2, SS, SS2 = "ECTX-SP008", "ECTX-SP009", "ECTX-SS003", "ECTX-SS004"


def result(sp="", ss="", reasons="", note="", status="completed"):
    return dict(experimental_eligibility_status=status, experimental_eligibility_reasons=reasons,
                hypothetical_species_id=sp, hypothetical_subspecies_id=ss, near_boundary_species_id=note, query_sha256="x")


def truth(kind="named", sp=SP, ss="", expected="", expected_ss="", gid="g1", layer="L1", block="", component="", source=None,
          passing="", sub_structural="", query="x"):
    if source is None:
        source = "project_label" if kind in sv.TRUTH_LAYERS else "none"
    return dict(genome_id=gid, layer=layer, truth_kind=kind, label_source=source, truth_species_id=sp, truth_subspecies_id=ss,
                ani_species_passing=passing, subspecies_structural=sub_structural, expected_species_ids=expected,
                expected_subspecies_id=expected_ss, block=block, component_id=component, query_sha256=query)


@pytest.mark.parametrize("res, species, subspecies, note", [
    (result(SP, SS), "correct", "correct", ""),
    (result(SP, SS2), "correct", "wrong_subspecies", ""),
    (result(SP), "correct", "species_only", ""),
    (result(reasons="nearest_below_species_ANI"), "withheld", "withheld_at_species", "absent"),
    (result(reasons="nearest_below_species_ANI", note=SP), "withheld", "withheld_at_species", "names_truth_species"),
    (result(reasons="nearest_below_species_ANI", note=SP2), "withheld", "withheld_at_species", "names_other_species"),
    (result(SP2, "ECTX-SS007"), "wrong", "wrong_species", ""),
])
def test_named_with_subspecies(res, species, subspecies, note):
    out = sv.score(truth(ss=SS), res)
    assert (out["species_outcome"], out["subspecies_outcome"], out["note_outcome"]) == (species, subspecies, note)


def test_named_without_subspecies_counts_overcall_separately():
    assert sv.score(truth(), result(SP, SS))["subspecies_outcome"] == "overcall_signal"
    assert sv.score(truth(), result(SP))["subspecies_outcome"] == "not_applicable"


@pytest.mark.parametrize("res", [result(SP, SS), result(SP), result(SP2, SS), result(reasons="nearest_below_species_ANI")])
def test_undetermined_subspecies_is_not_scored(res):
    for source in ("michael", "auto_single_species"):  # Michael's subspecies on his own or an automatic species
        t = dict(truth(ss=sv.UNDETERMINED, source=source, layer="expert", passing=SP), subspecies_label_source="michael")
        out = sv.score(t, res)
        assert out["subspecies_outcome"] == "truth_undetermined"
        assert out["species_outcome"] == sv.score(truth(ss=SS, source=source, layer="expert", passing=SP), res)["species_outcome"]
    for ss_source in ("automatic", ""):
        with pytest.raises(ValueError, match="only Michael's subspecies labels"):
            sv.score(dict(truth(ss=sv.UNDETERMINED, source="michael", layer="expert", passing=SP),
                          subspecies_label_source=ss_source), res)


@pytest.mark.parametrize("kind", ["outside_genus", "no_reference"])
def test_special_labels(kind):
    assert sv.score(truth(kind, sp=""), result(SP))["species_outcome"] == "wrong"
    withheld = sv.score(truth(kind, sp=""), result(reasons="species_coverage_not_met", note=SP))
    assert (withheld["species_outcome"], withheld["note_outcome"]) == ("correct_withheld", "names_a_species")


@pytest.mark.parametrize("kind", ["unresolved", "context", "crosswalk_none"])
def test_descriptive_kinds(kind):
    assert sv.score(truth(kind, sp=""), result(SP))["species_outcome"] == "descriptive_result"
    assert sv.score(truth(kind, sp=""), result(reasons="no_estimate"))["species_outcome"] == "descriptive_withheld"


@pytest.mark.parametrize("res, species, subspecies, note", [
    (result(SP, SS), "consistent", "consistent", ""),
    (result(SP, SS2), "consistent", "discordant", ""),
    (result(SP), "consistent", "consistent_at_species_no_subspecies", ""),
    (result("ECTX-SP025", "ECTX-SS010"), "discordant", "", ""),
    (result(reasons="two_parent_support_alert"), "withheld", "", "absent"),
    (result(reasons="nearest_below_species_ANI", note=SP), "withheld", "", "names_expected_species"),
    (result(reasons="nearest_below_species_ANI", note="ECTX-SP025"), "withheld", "", "names_other_species"),
])
def test_crosswalk(res, species, subspecies, note):
    out = sv.score(truth("crosswalk", sp="", expected=f"{SP};{SP2}", expected_ss=SS), res)
    assert (out["species_outcome"], out["subspecies_outcome"], out["note_outcome"]) == (species, subspecies, note)


def test_broader_crosswalk_without_subspecies():
    out = sv.score(truth("crosswalk", sp="", expected=f"{SP2};ECTX-SP025"), result("ECTX-SP025", "ECTX-SS009"))
    assert (out["species_outcome"], out["subspecies_outcome"]) == ("consistent", "")


def test_failed_runs_stay_counted():
    assert sv.score(truth(), None)["species_outcome"] == "run_failed"
    assert sv.score(truth(), result(reasons="native_search_failed", status="failed"))["species_outcome"] == "run_failed"


@pytest.mark.parametrize("t, res, message", [
    (truth("mystery"), result(SP), "unknown truth_kind"),
    (truth(), result(ss=SS, reasons="x"), "subspecies result without species"),
    (truth(), result(SP, note=SP), "near-boundary note alongside"),
    (truth(), result(), "withheld record without a reason"),
    (truth(sp=""), result(SP), "named truth without a species"),
    (truth("crosswalk", sp=""), result(SP), "crosswalk truth without expected"),
    (truth(source="guess"), result(SP), "unknown label_source"),
    (truth(source="none"), result(SP), "need a label_source"),
    (truth("crosswalk", sp="", expected=SP, source="michael"), result(SP), "crosswalk rows must use label_source none"),
    (truth(query="y"), result(SP), "query_sha256 does not match"),
    (truth(query=""), result(SP), "query_sha256 does not match"),
    (truth(sub_structural="yes"), result(SP), "subspecies_structural must be"),
    (truth(passing=SP), result(SP), "expert-stratum rows only"),
    (truth(source="auto_single_species"), result(SP), "exactly the label's species"),
    (truth(source="auto_single_species", passing=f"{SP};{SP2}"), result(SP), "exactly the label's species"),
    (truth(source="auto_tree_resolved", passing=SP), result(SP), "several passing species"),
    (truth(source="auto_tree_resolved", passing=f"{SP2};ECTX-SP025"), result(SP), "several passing species"),
])
def test_contract_violations_raise(t, res, message):
    with pytest.raises(ValueError, match=message):
        sv.score(t, res)


@pytest.mark.parametrize("n", [32, 39, 40, 45, 47, 48, 77])
def test_clopper_pearson_zero_matches_closed_form(n):
    # For k = 0 the exact two-sided upper limit is 1 - (alpha/2)^(1/n).
    assert sv.clopper_pearson(0, n) == pytest.approx([0.0, 1 - 0.025 ** (1 / n)], abs=2e-6)


@pytest.mark.parametrize("k, n, expected", [(5, 10, [0.187086, 0.812914]), (10, 10, [0.025 ** 0.1, 1.0])])
def test_clopper_pearson(k, n, expected):
    assert sv.clopper_pearson(k, n) == pytest.approx(expected, abs=2e-6)


def test_clopper_pearson_empty():
    assert sv.clopper_pearson(0, 0) is None


def test_summary_counts_and_denominators():
    rows = [sv.score(truth(gid=f"g{i}"), r) for i, r in enumerate(
        [result(SP), result(SP), result(SP2), result(reasons="species_coverage_not_met;two_parent_support_alert")])]
    rows.append(sv.score(truth("outside_genus", sp="", gid="o"), result(reasons="no_estimate")))
    rows.append(sv.score(truth("unresolved", sp="", gid="u"), result(SP)))
    rows.append(sv.score(truth(gid="f"), None))
    s = sv.summarise(rows)
    assert s["accuracy_denominator"] == 6 and s["species_results"] == 3 and s["wrong_results"] == 1
    assert s["named_coverage"] == [3, 5] and s["run_failed"] == 1
    assert s["withheld_reasons"] == {"no_estimate": 1, "species_coverage_not_met": 1, "two_parent_support_alert": 1}


def test_block_sensitivity_keeps_one_per_block():
    rows = [dict(sv.score(truth(gid=f"g{i}"), result(SP)), block=f"b{i % 3}") for i in range(9)]
    assert sv.block_sensitivity(rows)["accuracy_denominator"] == 3
    layer3 = [dict(sv.score(truth("crosswalk", sp="", expected=SP, gid=f"p{i}"), result(SP)), block=f"c{i % 3}") for i in range(9)]
    assert sv.block_sensitivity(layer3)["crosswalk_species"] == {"consistent": 3}
    with pytest.raises(ValueError, match="no block"):
        sv.block_sensitivity(rows + [dict(sv.score(truth(gid="z"), result(SP)), block="")])


def test_descriptive_rows_keep_their_label_source():
    assert sv.score(truth("context", sp="", source="project_label"), result(SP))["species_outcome"] == "descriptive_result"
    assert sv.score(truth("unresolved", sp="", source="michael"), result(reasons="no_estimate"))["species_outcome"] == \
        "descriptive_withheld"


def test_tool_species_outside_the_passing_set_is_flagged_for_any_label_source():
    rows = [sv.score(truth(source="michael", sp=SP2, passing=SP2), result(SP)),
            sv.score(truth(source="auto_tree_resolved", passing=f"{SP};{SP2}"), result(SP)),
            sv.score(truth("outside_genus", sp="", source="michael"), result(SP)),
            sv.score(truth(), result(SP2))]  # Layer 1 rows carry no passing set
    assert [r["result_outside_passing"] for r in rows] == [True, False, True, False]
    assert sv.summarise(rows)["tool_species_outside_passing"] == 2


def test_structural_subspecies_rows_are_counted_apart():
    rows = [sv.score(truth(ss=SS, source="auto_single_species", passing=SP, sub_structural="True", gid="s"), result(SP, SS)),
            sv.score(truth(ss=SS, source="auto_single_species", passing=SP, sub_structural="False", gid="o"), result(SP, SS2)),
            sv.score(truth(ss=SS, source="michael", passing=SP, gid="m"), result(SP, SS)),
            sv.score(truth(ss=SS, source="michael", passing=SP, gid="w"), result(SP))]
    s = sv.summarise(rows)
    assert s["subspecies_structural_outcomes"] == {"correct": 1}
    assert s["wrong_subspecies"] == [1, 2] and s["wrong_subspecies_ci95"] is not None


def test_component_weighting_gives_each_component_weight_one():
    rows = [dict(sv.score(truth("crosswalk", sp="", expected=SP, gid=f"g{i}"), result(SP if i else SP2)),
                 component_id="big" if i < 4 else f"c{i}") for i in range(6)]
    weights = sv.component_weighted(rows)
    assert sum(weights.values()) == pytest.approx(3) and weights["discordant"] == pytest.approx(0.25)


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_crossed_groups_and_grouped_component_weighting(tmp_path):
    rows, reports = [], []
    for i, (kind_, level) in enumerate([("broader", "Contig"), ("broader", "Contig"), ("broader", "Scaffold"), ("equivalent", "Contig")]):
        path = tmp_path / f"r{i}.tsv"
        write_tsv(path, [result(SP if i != 1 else SP2)])
        reports.append(dict(genome_id=f"p{i}", report_tsv=str(path), report_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        rows.append(dict(truth("crosswalk", sp="", expected=SP, gid=f"p{i}", layer="L3", block="c1" if i < 2 else f"c{i}",
                               component="c1" if i < 2 else f"c{i}"), crosswalk_kind=kind_, assembly_level=level))
    write_tsv(tmp_path / "truth.tsv", rows)
    write_tsv(tmp_path / "results.tsv", reports)
    sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results.tsv"),
             "--group", "crosswalk_kind,assembly_level", "--output", str(tmp_path / "out")])
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    key = "L3|crosswalk_kind=broader,assembly_level=Contig"
    assert summary["groups"][key]["crosswalk_species"] == {"consistent": 1, "discordant": 1}
    assert summary["component_weighted"][key] == {"consistent": 0.5, "discordant": 0.5}
    assert summary["block_sensitivity"]["L3"]["genomes"] == 3
    write_tsv(tmp_path / "truth2.tsv", [dict(rows[0], block="")])
    write_tsv(tmp_path / "results2.tsv", reports[:1])
    with pytest.raises(ValueError, match="no block"):
        sv.main(["--truth", str(tmp_path / "truth2.tsv"), "--results", str(tmp_path / "results2.tsv"),
                 "--output", str(tmp_path / "out2")])


def test_end_to_end(tmp_path):
    reports = []
    for gid, res in (("a", result(SP, SS)), ("b", result(reasons="nearest_below_species_ANI", note=SP))):
        path = tmp_path / f"{gid}.tsv"
        write_tsv(path, [res])
        reports.append(dict(genome_id=gid, report_tsv=str(path), report_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    write_tsv(tmp_path / "results.tsv", reports)
    write_tsv(tmp_path / "truth.tsv", [dict(truth(gid=g, ss=SS, block=g), subset="reserve") for g in ("a", "b", "c")])
    sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results.tsv"),
             "--group", "subset", "--output", str(tmp_path / "out")])
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    layer = summary["groups"]["L1"]
    assert layer["species_outcomes"] == {"correct": 1, "withheld": 1, "run_failed": 1}
    assert layer["note_outcomes"] == {"names_truth_species": 1}
    assert "L1|subset=reserve" in summary["groups"] and summary["block_sensitivity"]["L1"]["accuracy_denominator"] == 3
    with pytest.raises(SystemExit):
        sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results.tsv"),
                 "--output", str(tmp_path / "out")])
    with pytest.raises(ValueError, match=r"--group columns not in the truth table: \['subst'\]"):
        sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results.tsv"),
                 "--group", "subst", "--output", str(tmp_path / "out3")])


def test_tampered_report_and_unknown_genomes_rejected(tmp_path):
    path = tmp_path / "a.tsv"
    write_tsv(path, [result(SP)])
    write_tsv(tmp_path / "results.tsv", [dict(genome_id="a", report_tsv=str(path), report_sha256="0" * 64)])
    write_tsv(tmp_path / "truth.tsv", [truth(gid="a", block="b")])
    with pytest.raises(ValueError, match="hash mismatch"):
        sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results.tsv"),
                 "--output", str(tmp_path / "o1")])
    write_tsv(tmp_path / "results2.tsv", [dict(genome_id="zz", report_tsv=str(path), report_sha256="0" * 64)])
    with pytest.raises(ValueError, match="not in truth"):
        sv.main(["--truth", str(tmp_path / "truth.tsv"), "--results", str(tmp_path / "results2.tsv"),
                 "--output", str(tmp_path / "o2")])
    write_tsv(tmp_path / "truth2.tsv", [truth(gid="a", block="b"), truth(gid="a", block="b")])
    with pytest.raises(ValueError, match="duplicate"):
        sv.main(["--truth", str(tmp_path / "truth2.tsv"), "--results", str(tmp_path / "results.tsv"),
                 "--output", str(tmp_path / "o3")])


def test_structural_labels_are_excluded_from_the_wrong_result_rate():
    rows = [sv.score(truth(gid="s1", source="auto_single_species", passing=SP), result(SP)),
            sv.score(truth(gid="s2", source="auto_single_species", passing=SP), result(reasons="species_coverage_not_met")),
            sv.score(truth(gid="t1", source="auto_tree_resolved", sp=SP2, passing=f"{SP};{SP2}"), result(SP)),
            sv.score(truth(gid="m1", source="michael", passing=SP), result(SP)),  # Michael chose the only passing species
            sv.score(truth(gid="m2", source="michael", sp=SP2, passing=SP), result(SP)),  # trees overruled ANI
            sv.score(truth("outside_genus", sp="", gid="m3", source="michael"), result(reasons="no_estimate"))]
    assert [r["structural"] for r in rows] == [True, True, False, True, False, False]
    s = sv.summarise(rows)
    assert s["structural_rows"] == 3 and s["structural_rows_wrong"] == 0
    assert s["species_results"] == 2 and s["wrong_results"] == 2
    assert s["named_coverage"] == [4, 5]  # coverage still counts every named row
    assert s["label_sources"] == {"auto_single_species": 2, "auto_tree_resolved": 1, "michael": 3}
    # Layer 1 project labels are never structural.
    assert sv.score(truth(), result(SP))["structural"] is False
    # A structural row scored wrong would expose a scope difference between the labelling and tool skani runs.
    odd = sv.score(truth(source="michael", sp=SP2, passing=SP2), result(SP))
    assert sv.summarise([odd])["structural_rows_wrong"] == 1 and sv.summarise([odd])["wrong_results"] == 0
