"""R4-m1 check: label_transfer._verticall_ok on the complete paper QC gate (1,089 genomes, 1,088 pairs each)."""
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import label_transfer as lt  # noqa: E402

REVIEW = lt.ROOT / "results/full-cohort/2026-09-22/task11-opus-full/verticall-1089/qc/gate-before-any-matrix/assembly-review.tsv"
rows = list(csv.DictReader(REVIEW.open(), delimiter="\t"))
ok = lt._verticall_ok(REVIEW, [r["assembly_id"] for r in rows], 1088)
print(json.dumps(dict(review=str(REVIEW.relative_to(lt.ROOT)), review_sha256=hashlib.sha256(REVIEW.read_bytes()).hexdigest(),
                      genomes=len(rows), complete=sum(ok.values()), incomplete=sorted(g for g, v in ok.items() if not v),
                      with_secondary_rows=sum(int(r["pairs_with_secondary_results"]) > 0 for r in rows),
                      label_transfer_sha256=lt.sha256(Path(lt.__file__))), indent=2))
