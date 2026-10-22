#!/usr/bin/env python3
"""List files changed or added by the 9 Oct 2026 revision against the before-copies.

Compares the current run files (and the two repo files outside the run) with
superseded/2026-10-08-pre-review.sha256. Writes TSV to stdout:
path, status (modified|new|unchanged), sha256_before, sha256_after.
Raw API pages and per-contig base-count files are summarised by directory.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
from pathlib import Path

REPO_FILES = {"repo-files/docs-task-status-12-nomenclature-register-claude.md": "docs/task-status/12-nomenclature-register-claude.md",
              "repo-files/work-claims-issue-27-owner.json": ".work/claims/issue-27/owner.json"}
SKIP_PARTS = {"superseded", "sources", "inputs", "__pycache__"}
REVISION_START = dt.datetime.fromisoformat("2026-10-09T04:52:00+10:00").timestamp()  # before-copies taken


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, type=Path)
    a = ap.parse_args()
    before = {}
    for line in open(a.run / "superseded/2026-10-08-pre-review.sha256"):
        h, f = line.split(maxsplit=1)
        before[f.strip()[2:]] = h
    rows = []
    bulk: dict[str, list] = {}
    for p in sorted(x for x in a.run.rglob("*") if x.is_file()):
        rel = p.relative_to(a.run).as_posix()
        if SKIP_PARTS & set(rel.split("/")[:1]) or "__pycache__" in rel:
            continue
        if "/raw/" in rel or "/base-counts/" in rel:
            d = rel.split("/raw/")[0] + "/raw/" if "/raw/" in rel else rel.split("/base-counts/")[0] + "/base-counts/"
            bulk.setdefault(d, []).append(p.stat().st_mtime)
            continue
        h = sha(p)
        st = "new" if rel not in before else ("unchanged" if before[rel] == h else "modified")
        rows.append((rel, st, before.get(rel, ""), h))
    for k, v in REPO_FILES.items():
        h = sha(Path(v))
        rows.append((v, "unchanged" if before[k] == h else "modified", before[k], h))
    for f in sorted(a.run.glob("sources/gtdb/*")):
        rows.append((f.relative_to(a.run).as_posix(), "new", "", sha(f)))
    for d, mt in sorted(bulk.items()):
        if max(mt) < REVISION_START:
            st = f"unchanged ({len(mt)} files, all last modified before the revision; not copied)"
        elif min(mt) >= REVISION_START:
            st = f"new ({len(mt)} files; hashes in checks-raw.sha256)"
        else:
            st = f"MIXED ({len(mt)} files; check)"
        rows.append((d, st, "", ""))
    print("path\tstatus\tsha256_before\tsha256_after")
    for r in rows:
        print("\t".join(r))
    return 0


if __name__ == "__main__":
    sys.exit(main())
