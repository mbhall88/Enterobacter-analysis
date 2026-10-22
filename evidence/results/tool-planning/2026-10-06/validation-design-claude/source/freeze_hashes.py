"""Hash every item validation-design.md section 8 freezes, for the issue 3 freeze post.

Each file is read and hashed now; where the design records a hash, the recomputed value must match it.
The Verticall table is also checked against its own .sha256 file. Package source files are hashed one
by one and summarised as the SHA256 of the sorted "sha256  path" lines.

Usage: freeze_hashes.py OUTPUT_DIR   (refuses an existing directory)
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
RUN = Path(__file__).resolve().parents[1]
T11 = "results/full-cohort/2026-09-22/task11-opus-full"
A08 = "outputs/design-artifacts-attempt08"
TP = "results/tool-planning/2026-10-06"

# (group, role, path relative to ROOT or RUN, expected hash prefix…suffix or None)
ITEMS = [
    ("design", "validation design (revision 6)", ("run", "validation-design.md"), "6e09d49f…7a0d"),
    ("design", "frozen dispositions (attempt03)", ("run", "outputs/exposure-audit-attempt03/reserve-exposure-dispositions.tsv"), "71d0d86b…da81"),
    ("artifacts", "NCBI-name crosswalk", ("run", f"{A08}/ncbi-name-crosswalk.tsv"), "cedeefaa"),
    ("artifacts", "stratified expert sample (n=100)", ("run", f"{A08}/expert-sample-option-stratified-n100.tsv"), "f2787df9"),
    ("artifacts", "Layer 1 manifest", ("run", f"{A08}/layer1-manifest.tsv"), "91e798a8"),
    ("artifacts", "labelling comparators", ("run", f"{A08}/labelling-comparators.tsv"), "0298309c"),
    ("artifacts", "cohort tree tips", ("run", f"{A08}/cohort-tree-tips.tsv"), "424fb3d1"),
    ("artifacts", "Layer 3 pool strata", ("run", f"{A08}/layer3-pool-strata.tsv"), "f99c844e"),
    ("artifacts", "artifact build summary", ("run", f"{A08}/summary.json"), "875d298f"),
    ("code", "scorer", ("run", "source/score_validation.py"), "2afa6d76…a19f"),
    ("code", "scorer tests", ("run", "source/test_score_validation.py"), "90f9206d…ec8ec"),
    ("code", "label transfer", ("run", "source/label_transfer.py"), "4fc90dab…71a8"),
    ("code", "label transfer tests", ("run", "source/test_label_transfer.py"), "2702cae4…80b4"),
    ("code", "Verticall extension", ("run", "source/verticall_extension.py"), "c7a46181…cf6"),
    ("code", "Verticall extension tests", ("run", "source/test_verticall_extension.py"), "330d038a…a33a7d"),
    ("code", "disagreement draw", ("run", "source/draw_disagreement_review.py"), "11b6b408"),
    ("code", "disagreement draw tests", ("run", "source/test_draw_disagreement_review.py"), "2b19c01f…fe1c"),
    ("code", "artifact builder", ("run", "source/build_design_artifacts.py"), "0d140649…0aa57"),
    ("labelling inputs", "derivative map (repair)", ("root", f"{T11}/repair/source-derivative-map.tsv"), "d073eb87…eaef"),
    ("labelling inputs", "derivative map (freeze-1089)", ("root", f"{T11}/freeze-1089/source-derivative-map.tsv"), None),
    ("labelling inputs", "BAC120 user MSA", ("root", f"{T11}/marker/alignments/full.user_msa.fasta.gz"), None),
    ("labelling inputs", "BAC120 filtered alignment", ("root", f"{T11}/marker/filtered/alignment.faa"), None),
    ("labelling inputs", "BAC120 column map", ("root", f"{T11}/marker/filtered/column-map.tsv"), None),
    ("labelling inputs", "BAC120 partitions", ("root", f"{T11}/marker/filtered/partitions.nex"), None),
    ("labelling inputs", "Verticall paper table", ("root", f"{T11}/verticall-1089/full-1089-pairwise.tsv"), "963576bd…1ca0"),
    ("labelling inputs", "Verticall skip list", ("run", "outputs/verticall-skip-list-attempt01/paper-directed-pairs.tsv"), "5a5209ef…a733"),
    ("labelling inputs", "paper BAC120 tree", ("root", f"{T11}/pruned-1089/marker/fasttree-wag.nwk"), None),
    ("labelling inputs", "paper Verticall tree", ("root", f"{T11}/pruned-1089/verticall/median-first.nwk"), None),
    ("wrappers", "repair_pilot_inputs.py", ("root", "scripts/repair_pilot_inputs.py"), None),
    ("wrappers", "run_marker_cohort.sh", ("root", "scripts/run_marker_cohort.sh"), None),
    ("wrappers", "subset_marker_alignment.py", ("root", "scripts/subset_marker_alignment.py"), None),
    ("wrappers", "run_verticall_cohort.sh", ("root", "scripts/run_verticall_cohort.sh"), None),
    ("wrappers", "qc_verticall.py", ("root", "scripts/qc_verticall.py"), None),
    ("wrappers", "gtdb_context.py (imported by label_transfer.py)", ("root", "scripts/gtdb_context.py"), None),
    ("wrappers", "root_tree.py (imported by label_transfer.py)", ("root", "scripts/root_tree.py"), None),
    ("pins", "reference package", ("root", f"{TP}/portable-reference-package-astra/relocated-package/portable.json"),
     "12aa0705b774de8ea262869fdc4c85d513c3aa97d39c674f08e85d1c9055b88f"),
    ("pins", "rank rule", ("root", f"{TP}/frozen-eligibility-claude/inputs/rank-rule.json"),
     "5849c5db6649678963215ca9363a5dfae2ef2b0282d96946aa843814fafa4188"),
    ("pins", "eligibility configuration 2.1", ("root", f"{TP}/near-boundary-claude/inputs/frozen-eligibility-config-2.1.json"),
     "3188518eb7bdaad287e7f335ffb18cadfbb41b46308f018f5b69ff33cab7c626"),
    ("pins", "MLST configuration", ("root", f"{TP}/mlst-runner-integration-astra/inputs/mlst-config.json"),
     "1fce9a33cf423967e0c089830ee3e6549930b5d45761ab95ac662ccd1430a937"),
    ("pins", "skani 0.3.2 binary", ("abs", "/sw/local/rocky9/noarch/rcc/software/gtdb-tk/gtdbtk_2.7.2/bin/skani"),
     "b3000043feb4ce8a5200b3716bae155f34f73e2a123b8b055d574a7897a14151"),
]
RANK_RULE_CANONICAL = "de3ec713ce3b902127c4402e893d84bb74ea38bb6faf5f82ae65fb03b9be4018"


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def matches(digest, expected):
    if expected is None:
        return True
    if "…" in expected:
        head, tail = expected.split("…")
        return digest.startswith(head) and digest.endswith(tail)
    return digest.startswith(expected)


def resolve(where):
    kind, rel = where
    return {"run": RUN / rel, "root": ROOT / rel, "abs": Path(rel)}[kind], (rel if kind != "run" else str((RUN / rel).relative_to(ROOT)))


def main(argv):
    out = Path(argv[0])
    if out.exists():
        sys.exit("refusing to overwrite output")
    rows, failures = [], []
    for group, role, where, expected in ITEMS:
        path, shown = resolve(where)
        digest = sha256(path)
        ok = matches(digest, expected)
        if not ok:
            failures.append(f"{shown}: {digest} does not match recorded {expected}")
        rows.append(dict(group=group, role=role, path=shown, sha256=digest, bytes=path.stat().st_size,
                         recorded=expected or "", matches_record=ok if expected else ""))
    table = ROOT / T11 / "verticall-1089/full-1089-pairwise.tsv"
    recorded_table = (ROOT / T11 / "verticall-1089/full-1089-pairwise.sha256").read_text().split()[0]
    if next(r["sha256"] for r in rows if r["path"].endswith("full-1089-pairwise.tsv")) != recorded_table:
        failures.append("Verticall paper table differs from its .sha256 file")
    sys.path.insert(0, str(ROOT / "src"))
    from enterobacter_classifier.rank import _canonical_digest  # noqa: E402
    canonical = _canonical_digest(json.loads((ROOT / f"{TP}/frozen-eligibility-claude/inputs/rank-rule.json").read_text()))
    if canonical != RANK_RULE_CANONICAL:
        failures.append(f"rank rule canonical digest {canonical} != {RANK_RULE_CANONICAL}")
    package = []
    for path in sorted((ROOT / "src/enterobacter_classifier").rglob("*")):
        if path.is_file() and "__pycache__" not in path.parts:
            package.append(dict(path=str(path.relative_to(ROOT)), sha256=sha256(path), bytes=path.stat().st_size))
    package_digest = hashlib.sha256("".join(f"{p['sha256']}  {p['path']}\n" for p in package).encode()).hexdigest()
    if failures:
        sys.exit("\n".join(failures))
    out.mkdir(parents=True)
    for name, data in (("freeze-hashes.tsv", rows), ("package-source.tsv", package)):
        with (out / name).open("x") as handle:
            handle.write("\t".join(data[0]) + "\n")
            handle.writelines("\t".join(str(v) for v in r.values()) + "\n" for r in data)
    record = dict(items=len(rows), package_files=len(package), package_source_digest=package_digest,
                  rank_rule_canonical_digest=canonical, verticall_table_matches_sha256_file=True,
                  source_sha256=sha256(__file__))
    (out / "summary.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main(sys.argv[1:])
