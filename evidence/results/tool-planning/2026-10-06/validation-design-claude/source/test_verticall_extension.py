"""Synthetic tests for verticall_extension.py."""
import hashlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verticall_extension as ve  # noqa: E402

HEADER = "assembly_a\tassembly_b\talignment_count\tmedian_vertical_window_distance\n"


def paper(tmp_path, rows=(("c1", "c2"), ("c2", "c1"), ("c1", "c3"), ("c3", "c1"), ("c2", "c3"), ("c3", "c2"), ("c2", "c3"))):
    table = tmp_path / "paper.tsv"
    table.write_text(HEADER + "".join(f"{a}\t{b}\t5\t0.01\n" for a, b in rows))
    sha = tmp_path / "paper.sha256"
    sha.write_text(hashlib.sha256(table.read_bytes()).hexdigest() + "  paper.tsv\n")
    return table, sha


def part(tmp_path, name, rows, header=False):
    path = tmp_path / name
    path.write_text((HEADER if header else "") + "".join(f"{a}\t{b}\t5\t0.02\n" for a, b in rows))
    return path


def new_ids(tmp_path, ids=("n1",)):
    path = tmp_path / "new.txt"
    path.write_text("".join(i + "\n" for i in ids))
    return path


def test_skip_list(tmp_path):
    table, sha = paper(tmp_path)
    out = tmp_path / "skip.tsv"
    record = ve.skip_list(table, sha, out, expected_pairs=6)  # the secondary row repeats a pair
    lines = out.read_text().splitlines()
    assert lines[0] == "assembly_a\tassembly_b" and len(lines) == 7 and record["directed_pairs"] == 6
    with pytest.raises(ValueError, match="expected 5"):
        ve.skip_list(table, sha, tmp_path / "skip2.tsv", expected_pairs=5)
    sha.write_text("0" * 64 + "\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        ve.skip_list(table, sha, tmp_path / "skip3.tsv")


def test_merge(tmp_path):
    table, sha = paper(tmp_path)
    p1 = part(tmp_path, "p1.tsv", [("n1", "c1"), ("c1", "n1")], header=True)
    p2 = part(tmp_path, "p2.tsv", [("n1", "c2"), ("c2", "n1"), ("n1", "c3")])
    p3 = part(tmp_path, "p3.tsv", [])
    out = tmp_path / "merged.tsv"
    record = ve.merge(table, sha, [p1, p2, p3], new_ids(tmp_path), out, expected_parts=3)
    lines = out.read_text().splitlines()
    assert lines[0] + "\n" == HEADER and sum(line.startswith("assembly_a") for line in lines) == 1
    assert (record["paper_rows"], record["new_rows"], record["merged_rows"]) == (7, 5, 12)
    with pytest.raises(ValueError, match="2 parts given, expected 3"):
        ve.merge(table, sha, [p1, p2], new_ids(tmp_path), tmp_path / "merged2.tsv", expected_parts=3)


@pytest.mark.parametrize("parts, message", [
    ([("p1", [("n1", "c1")], False)], "header differs"),
    ([("p1", [("n1", "c1")], True), ("p2", [("c1", "n1")], True)], "only part 1 may"),
    ([("p1", [("c1", "c2")], True)], "repeats a paper pair"),
    ([("p1", [("c3", "c4")], True)], "names no new genome"),
])
def test_merge_refuses(tmp_path, parts, message):
    table, sha = paper(tmp_path)
    paths = [part(tmp_path, f"{n}.tsv", rows, header) for n, rows, header in parts]
    with pytest.raises(ValueError, match=message):
        ve.merge(table, sha, paths, new_ids(tmp_path), tmp_path / "merged.tsv")


def test_merge_refuses_a_paper_table_naming_a_new_genome(tmp_path):
    table, sha = paper(tmp_path, rows=(("c1", "n1"),))
    with pytest.raises(ValueError, match="already names a new genome"):
        ve.merge(table, sha, [part(tmp_path, "p1.tsv", [], True)], new_ids(tmp_path), tmp_path / "m.tsv")


def write(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return path


def test_stage_map(tmp_path):
    freeze = write(tmp_path / "freeze.tsv", ["tip_id", "derivative_fasta_path"], [["c1", "d/c1.fna"], ["c2", "d/c2.fna"]])
    repair = write(tmp_path / "repair.tsv", ["tip_id", "derivative_fasta_path", "derivative_sha256"],
                   [["c1", "d/c1.fna", "h1"], ["c2", "d/c2.fna", "h2"], ["x9", "d/x9.fna", "h9"]])
    cols = ["tip_id", "source_fasta_sha256", "derivative_fasta_path", "derivative_sha256"]
    new = write(tmp_path / "new.tsv", cols, [["n1", "q1", "e/n1.fna", "hn"]])
    query = write(tmp_path / "query.tsv", ["accession", "query_sha256"], [["n1", "q1"]])
    record = ve.stage_map(repair, freeze, new, query, tmp_path / "map.tsv", tmp_path / "ids.txt")
    assert record == dict(cohort=2, new=1, total=3)
    assert (tmp_path / "ids.txt").read_text().split() == ["c1", "c2", "n1"]
    assert "x9" not in (tmp_path / "map.tsv").read_text()
    clash = write(tmp_path / "clash.tsv", cols, [["c1", "q1", "e/c1.fna", "h"]])
    clash_query = write(tmp_path / "clash-query.tsv", ["accession", "query_sha256"], [["c1", "q1"]])
    with pytest.raises(ValueError, match="clash"):
        ve.stage_map(repair, freeze, clash, clash_query, tmp_path / "m2.tsv", tmp_path / "i2.txt")
    moved = write(tmp_path / "freeze2.tsv", ["tip_id", "derivative_fasta_path"], [["c1", "other/c1.fna"]])
    with pytest.raises(ValueError, match="derivative path differs"):
        ve.stage_map(repair, moved, new, query, tmp_path / "m3.tsv", tmp_path / "i3.txt")
    nohash = write(tmp_path / "nohash.tsv", cols, [["n1", "q1", "e/n1.fna", ""]])
    with pytest.raises(ValueError, match="missing hashes"):
        ve.stage_map(repair, freeze, nohash, query, tmp_path / "m4.tsv", tmp_path / "i4.txt")
    other = write(tmp_path / "other-query.tsv", ["accession", "query_sha256"], [["n1", "q2"]])
    with pytest.raises(ValueError, match="not the query manifest's FASTA"):
        ve.stage_map(repair, freeze, new, other, tmp_path / "m5.tsv", tmp_path / "i5.txt")
    extra = write(tmp_path / "extra-query.tsv", ["accession", "query_sha256"], [["n1", "q1"], ["n2", "q2"]])
    with pytest.raises(ValueError, match="differ from the query manifest"):
        ve.stage_map(repair, freeze, new, extra, tmp_path / "m6.tsv", tmp_path / "i6.txt")
