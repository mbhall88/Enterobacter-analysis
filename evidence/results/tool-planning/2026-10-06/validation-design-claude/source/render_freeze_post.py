"""Render the issue 3 freeze post from a freeze_hashes.py output directory.

Usage: render_freeze_post.py FREEZE_DIR TARGETED_CHECK_README OUTPUT_MD   (refuses an existing output)
"""
import csv
import json
import sys
from pathlib import Path

HEAD = """**Posted by:** Claude Opus 5.5 (`claude-opus-5-5[1m]`)

**Validation design frozen: revision 6 (DEC-005).**

- **Review.** The third independent review passed with 8 minor findings. Three follow-up checks by the same reviewer:
  - **revision 4:** fixes accepted; 3 minor findings;
  - **revision 5:** one major defect, a scoring stop if Michael left a subspecies undetermined on an automatic species;
  - **revision 6:** both fixes confirmed, nothing left open (`{check}`).

  None of these findings needed a decision from Michael.
- **Sequences.** No public validation genome sequence was accessed before this post. Download (design section 8, step 2) starts after it.
- **Tests.** {tests}
- **How the hashes were made.** Every file below was hashed again for this post by `source/freeze_hashes.py` (SHA256 `{source}`). Each hash the design records was checked and matches. The Verticall table also matches its own `.sha256` file.
- **Paths.** Paths are repository-relative. Files in the design's run directory, `{run}`, are shown relative to it.
"""

GROUPS = [
    ("design", "Design and frozen dispositions"),
    ("artifacts", "Design artifacts (`outputs/design-artifacts-attempt08/`)"),
    ("code", "Frozen code, with tests"),
    ("labelling inputs", "Labelling inputs (section 7.1)"),
    ("wrappers", "Project wrappers and the modules `label_transfer.py` imports"),
    ("pins", "Object under test"),
]

TAIL = """
**Package source.** The package source is pinned by these hashes, not by a commit. `reporting.py` and `workflow.py` include the near-boundary integration (configuration 2.1), which is not yet committed. It has been unchanged since 6 October, 23:02 AEST.
- **Files:** {n} under `src/enterobacter_classifier/`.
- **Digest:** `{digest}`, the SHA256 of the sorted `sha256  path` lines in `package-source.tsv`.

<details><summary>Per-file package hashes</summary>

| File | SHA256 |
|---|---|
{package_rows}
</details>

**The rank rule's canonical digest** is `{canonical}`, the value configuration 2.1 binds.

**Restated at freeze** (third review, n5):
- **Michael's own labels.** A genome Michael labels whose label is the only species passing 95/80 is also excluded from the wrong-result rate. Such a label cannot be wrong. The posted Q-AL wording named only automatic labels.
- **ANI scope.** Labelling ANI is whole-assembly, a disclosed departure from DEC-007's chromosome scope.

**Next steps, under separate claims:**
1. Download the 4,857 public assemblies with `datasets` 18.25.1 (dehydrated, rehydrate, `md5sum -c`). No classifier run.
2. Label the expert sample: skani, the two tree rebuilds, then `label_transfer.py prepare`. Michael labels only the edge cases, then `seal`, and the sealed hash is posted here.
3. Run and score: Layer 1 first, then Layer 3 and Layer 2.

Positive-assignment activation stays Michael's decision.

Manifest: `{run}{freeze_dir}/` (`freeze-hashes.tsv`, `package-source.tsv`, `summary.json`).
"""


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(argv):
    freeze, check, output = Path(argv[0]), argv[1], Path(argv[2])
    tests = argv[3] if len(argv) > 3 else ""
    if output.exists():
        sys.exit("refusing to overwrite output")
    rows = read_tsv(freeze / "freeze-hashes.tsv")
    package = read_tsv(freeze / "package-source.tsv")
    summary = json.loads((freeze / "summary.json").read_text())
    run = "results/tool-planning/2026-10-06/validation-design-claude/"
    parts = [HEAD.format(check=check, tests=tests, source=summary["source_sha256"], run=run)]
    for key, title in GROUPS:
        parts.append(f"\n**{title}**\n\n| Item | Path | SHA256 |\n|---|---|---|\n")
        parts.append("".join(f"| {r['role']} | `{r['path'].removeprefix(run)}` | `{r['sha256']}` |\n" for r in rows if r["group"] == key))
    freeze_rel = str(freeze.resolve()).split(run, 1)[-1]
    parts.append(TAIL.format(n=len(package), digest=summary["package_source_digest"], canonical=summary["rank_rule_canonical_digest"],
                             package_rows="\n".join(f"| `{p['path'].removeprefix('src/enterobacter_classifier/')}` | `{p['sha256']}` |"
                                                    for p in package),
                             freeze_dir=freeze_rel, run=run))
    with output.open("x") as handle:
        handle.write("".join(parts))


if __name__ == "__main__":
    main(sys.argv[1:])
