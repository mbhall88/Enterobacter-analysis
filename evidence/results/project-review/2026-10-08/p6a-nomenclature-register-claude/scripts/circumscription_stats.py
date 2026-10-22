#!/usr/bin/env python3
"""Genome size and G+C ranges for adopted taxa (emended-description drafting aid).

Descriptive only: no labels, rules or denominators change. Members come from the
frozen final accession overlay (study records plus reference-only type anchors);
sequence scope is the whole deposited assembly listed in the frozen study
manifest (or, for reference-only types, the verified type FASTA). Each file's
SHA-256 is checked against the manifest before use; mismatches are excluded and
reported. G+C is computed over A/C/G/T only, from integer per-contig counts
(`seqkit fx2tab --base-count`) summed before dividing; the per-contig counts
are kept (gzip) under OUT/base-counts/. Revised 9 Oct 2026 (coordinator review
R1): the earlier version rebuilt counts from rounded percentages (-B).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

TAXA = {  # (project_species, project_subspecies or None for all members)
    "E. xiangfangensis (A+B)": ("Enterobacter xiangfangensis", None),
    "E. xiangfangensis subsp. xiangfangensis (B)": ("Enterobacter xiangfangensis", "xiangfangensis"),
    "E. xiangfangensis subsp. oharae (B)": ("Enterobacter xiangfangensis", "oharae"),
    "E. xiangfangensis subsp. steigerwaltii (B)": ("Enterobacter xiangfangensis", "steigerwaltii"),
    "E. xiangfangensis subsp. hoffmannii (A)": ("Enterobacter xiangfangensis", "hoffmannii"),
    "E. hormaechei (narrow)": ("Enterobacter hormaechei", None),
    "E. roggenkampii": ("Enterobacter roggenkampii", None),
    "E. roggenkampii subsp. roggenkampii": ("Enterobacter roggenkampii", "roggenkampii"),
    "E. roggenkampii subsp. quasiroggenkampii": ("Enterobacter roggenkampii", "quasiroggenkampii"),
    "E. chengduensis": ("Enterobacter chengduensis", None),
    "E. chengduensis subsp. chengduensis": ("Enterobacter chengduensis", "chengduensis"),
    "E. chengduensis subsp. rongchengensis": ("Enterobacter chengduensis", "rongchengensis"),
    "E. mori": ("Enterobacter mori", None),
    "E. asburiae": ("Enterobacter asburiae", None),
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def base_counts(p: Path) -> str:
    """Per-contig integer A/C/G/T counts (seqkit fx2tab -C; case-insensitive)."""
    return subprocess.run(["seqkit", "fx2tab", "-n", "-l", "-C", "A", "-C", "C", "-C", "G", "-C", "T", str(p)],
                          check=True, capture_output=True, text=True).stdout


def sum_counts(table: str) -> tuple[int, int, int, int, int]:
    L = a = c = g = t = 0
    for line in table.splitlines():
        f = line.split("\t")
        L += int(f[1])
        a += int(f[2])
        c += int(f[3])
        g += int(f[4])
        t += int(f[5])
    return L, a, c, g, t


def measure(item, pc_dir: Path | None = None):
    acc, path, expected = item
    p = Path(path)
    if not p.exists():
        return acc, path, "", "missing_file", 0, 0.0, (0, 0, 0, 0)
    got = sha256(p)
    if expected and got != expected:
        return acc, path, got, "sha256_mismatch", 0, 0.0, (0, 0, 0, 0)
    table = base_counts(p)
    if pc_dir is not None:
        (pc_dir / f"{acc}.base-counts.tsv.gz").write_bytes(
            gzip.compress(("name\tlength\tA\tC\tG\tT\n" + table).encode(), mtime=0))
    L, a, c, g, t = sum_counts(table)
    return acc, path, got, "ok" if expected else "ok_no_expected_hash", L, 100 * (g + c) / (a + c + g + t), (a, c, g, t)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlay", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--type-paths", required=True, type=Path, help="ACC=PATH lines for reference-only types")
    ap.add_argument("--type-gc", required=True, type=Path, help="type-genome-gc.tsv with verified type file hashes")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    man = {}
    for r in csv.DictReader(open(a.manifest), delimiter="\t"):
        path = r["resolved_whole_path"] or r["whole_path"]
        sha = r["resolved_whole_sha256"] or r["whole_sha256"]
        man[r["accession_version"]] = (path, sha)
    tpath = dict(line.strip().split("=", 1) for line in open(a.type_paths) if "=" in line)
    tsha = {r["accession"]: r["fasta_sha256"] for r in csv.DictReader(open(a.type_gc), delimiter="\t")}
    members: dict[str, list[str]] = {k: [] for k in TAXA}
    roles = {}
    for r in csv.DictReader(open(a.overlay), delimiter="\t"):
        for k, (sp, ssp) in TAXA.items():
            if r["project_species"] == sp and (ssp is None or r["project_subspecies"] == ssp):
                members[k].append(r["accession"])
                roles[r["accession"]] = r["record_role"]
    accs = sorted({x for v in members.values() for x in v})
    items = []
    for acc in accs:
        if acc in man:
            items.append((acc,) + man[acc])
        elif acc in tpath:
            items.append((acc, tpath[acc], tsha.get(acc, "")))
        else:
            items.append((acc, "", ""))
    pc = a.out / "base-counts"
    pc.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(a.threads) as ex:
        res = {r[0]: r for r in ex.map(lambda it: measure(it, pc), items)}
    with open(a.out / "member-genome-stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "record_role", "fasta_path", "fasta_sha256", "status", "total_length", "gc_mol_percent_acgt",
                    "A", "C", "G", "T"])
        for acc in accs:
            _, p, s, st, L, gc, cnt = res[acc]
            w.writerow([acc, roles.get(acc, ""), p, s, st, L, f"{gc:.2f}" if L else "", *cnt])
    with open(a.out / "circumscription-stats.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["taxon", "n_members_overlay", "n_measured", "n_excluded", "size_min_bp", "size_median_bp",
                    "size_max_bp", "gc_min", "gc_median", "gc_max"])
        for k in TAXA:
            ok = [res[x] for x in members[k] if res[x][3].startswith("ok")]
            sizes = [x[4] for x in ok]
            gcs = [x[5] for x in ok]
            w.writerow([k, len(members[k]), len(ok), len(members[k]) - len(ok),
                        min(sizes), int(statistics.median(sizes)), max(sizes),
                        f"{min(gcs):.2f}", f"{statistics.median(gcs):.2f}", f"{max(gcs):.2f}"])
    bad = [r for r in res.values() if not r[3].startswith("ok")]
    print(f"genomes={len(accs)} measured={len(accs) - len(bad)} excluded={len(bad)}")
    for r in bad:
        print("EXCLUDED", r[0], r[3], r[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
