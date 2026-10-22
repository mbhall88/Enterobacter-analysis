"""Tests for draw_disagreement_review.py."""
import csv
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import draw_disagreement_review as dd  # noqa: E402


def test_all_reviewed_at_or_below_the_cap():
    assert dd.allocate({"a": 12, "b": 8}) == {"a": 12, "b": 8}


def test_one_slot_per_form_then_proportional():
    counts = {"a": 40, "b": 10, "c": 1, "d": 1}
    slots = dd.allocate(counts)
    assert sum(slots.values()) == 20 and slots["c"] == 1 and slots["d"] == 1
    # 16 slots left over 39 and 9 remaining discordants: 13.0 and 3.0.
    assert (slots["a"], slots["b"]) == (14, 4)


def test_more_forms_than_slots():
    counts = {f"f{i}": 2 for i in range(25)}
    slots = dd.allocate(counts)
    assert sum(slots.values()) == 20 and max(slots.values()) == 1
    assert {f for f, s in slots.items() if s} == set(sorted(counts, key=dd.order)[:20])


@pytest.mark.parametrize("counts", [{"a": 21}, {"a": 30, "b": 2, "c": 3}, {"a": 7, "b": 7, "c": 7, "d": 7}])
def test_caps_and_totals(counts):
    slots = dd.allocate(counts)
    assert sum(slots.values()) == 20 and all(slots[f] <= counts[f] for f in counts)


def test_draw_eligibility_and_order():
    scored, truth = [], []
    for i in range(30):
        gid = f"g{i:02d}"
        kind = "crosswalk" if i < 28 else "crosswalk_none"
        species = "discordant" if i % 2 == 0 else "consistent"
        sub = "discordant" if i in (1, 3) else ""
        scored.append(dict(genome_id=gid, truth_kind=kind, species_outcome=species, subspecies_outcome=sub))
        truth.append(dict(genome_id=gid, ncbi_name_normalised="A" if i < 15 else "B"))
    chosen, record = dd.draw(scored, truth, expert={"g00"})
    ids = {c["genome_id"] for c in chosen}
    assert "g00" not in ids and "g28" not in ids and {"g01", "g03"} <= ids
    assert record["eligible"] == 15 and len(chosen) == 15  # 13 species-discordant (g02..g26 even) + g01 + g03
    again, _ = dd.draw(list(reversed(scored)), truth, expert={"g00"})
    assert again == chosen


def test_main(tmp_path):
    def write(path, rows):
        with path.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    write(tmp_path / "scored.tsv", [dict(genome_id=f"g{i}", truth_kind="crosswalk", species_outcome="discordant",
                                         subspecies_outcome="") for i in range(25)])
    write(tmp_path / "truth.tsv", [dict(genome_id=f"g{i}", ncbi_name_normalised=f"N{i % 3}") for i in range(25)])
    write(tmp_path / "expert.tsv", [dict(assembly_unit="g0")])
    dd.main(["--scored", str(tmp_path / "scored.tsv"), "--truth", str(tmp_path / "truth.tsv"),
             "--expert-sample", str(tmp_path / "expert.tsv"), "--output", str(tmp_path / "out")])
    record = json.loads((tmp_path / "out" / "draw.json").read_text())
    assert record["eligible"] == 24 and record["chosen"] == 20
    with pytest.raises(SystemExit):
        dd.main(["--scored", str(tmp_path / "scored.tsv"), "--truth", str(tmp_path / "truth.tsv"),
                 "--expert-sample", str(tmp_path / "expert.tsv"), "--output", str(tmp_path / "out")])
