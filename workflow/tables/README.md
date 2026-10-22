# Reconstruct the manuscript tables

This package reconstructs the 23 numbered supplementary TSVs from saved measurement products and selected metadata. It uses Python's standard library; it does not infer trees, measure ANI, rerun classifiers or refresh taxonomy. The expected hashes identify the author-approved table versions.

```sh
python3 rebuild_tables.py --output rebuilt-tables
python3 -m unittest -v test_rebuild.py
```

The output directory must be new. Each compressed input and its decoded content are hash-checked before any output is written. Successful reconstruction means all 23 TSVs have the exact expected bytes, including column order, formatting, empty fields and line endings. The result is recorded in `reconstruction-checks.json`.

The original table routines are identified by source path and hash in `routine-provenance.json`. Dataset counts use the original calculation with only the private documentary receipt-writing removed. Assignment/marker joins, type-display IDs, reference-column selection, coded JSONL views, safety-case views and timing-column views follow their original manuscript transformations. Twelve tables were originally copied from earlier saved measurements/registers; this package repeats that copy from those pre-manuscript products. Their inclusion is not a claim to rerun the underlying scientific analyses. See the analysis repository's workflow and source manifests for those upstream steps.

`input-inventory.json` is the explicit input list. It distinguishes the original source SHA256 from the exported compressed and decoded hashes. Selected metadata inputs omit unused fields; count-only views contain one ordinal per source row. Rescue JSONL preserves the nested scientific fields used by the original renderer while omitting unused private receipts and aliases. Table S12/S14's `source_record_sha256` therefore remains the original full-file digest, not the hash of the public projection. Cost inputs replace only private archive/member fields that the renderer discards. No hidden query identity mapping is needed or supplied.

Some source labels and discussion URLs remain in the frozen scientific outputs as dated provenance. The supplementary dataset's source index, interpretation note and public evidence mappings explain them without requiring access to private discussions. Table S20 is the historical nomenclatural register; its draft authority placeholders and older questions are not the final manuscript's formal naming statements.

The original command constructing the twelve-case safety catalogue was not recovered during this packaging review. Its exact retained input is hashed in the inventory, and this package reproduces the manuscript view of those records. This limitation does not prevent exact table reconstruction; it limits the claim about reproducing the earlier case-selection step.
