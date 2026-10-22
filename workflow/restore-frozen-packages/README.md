# Restore the packages used in the classifier evaluations

This helper restores the exact original reference, MLST and ANIm package bytes used by the frozen Python and Rust evaluations. It verifies the archives before writing, restores only named files, and checks every restored file and original bound manifest. It does not run a classifier or change any classification rule.

Place these four archives together in a directory:

- `enterobacter-frozen-evaluation-packages-noticed-2026-10-22.tar.gz`, from the analysis evidence archive;
- `enterobacter-reference-2026-10-15-reference.tar.gz`, from the published reference database;
- `enterobacter-reference-2026-10-15-primary-fasta.tar.gz`, from the same database;
- `enterobacter-reference-2026-10-15-mlst-historical-provenance-trimmed.tar.gz`, from the same database.

`restoration-manifest.json` pins their exact sizes and SHA256 values. The reference database's later sketch files are deliberately unused: they are not byte-identical to the evaluated sketches. Its six matching small metadata assets, 140 exact ANIm FASTAs, and MLST manifest and allele/profile files can be reused. The supplemental archive supplies the original sketches, wrappers and remaining metadata, plus the frozen Python MLST BLAST database.

Run with Python 3.10 or later and about 1.1 GB of free output space:

```sh
python3 restore_packages.py --archive-dir /path/to/archives --output restored-packages
```

The output directory must not already exist. The helper rejects path traversal, symbolic/hard links in archives, duplicate members/targets and checksum differences. It never uses `extractall`. If restoration fails after output creation, an `INCOMPLETE` marker remains; preserve it for diagnosis or remove that failed output directory before retrying. A successful run writes `restoration-receipt.json` and removes the marker.

The 249 restored files (1,054,707,423 bytes) are organised as:

- `python-reference/`: 42 original files, including the manifest with SHA256 `12aa0705b774de8ea262869fdc4c85d513c3aa97d39c674f08e85d1c9055b88f`;
- `rust-reference/`: 44 files, including manifest SHA256 `5743c25bcd706242ff8a28863fd6af2bf60e76f4c6af8b50d3f77deaf2e01dae`;
- `mlst/`: the original manifest, attribution and eight allele/profile files;
- `rust-anim-reference/`: 140 whole-assembly reference FASTAs and their two original wrappers;
- `python-mlst-blast/`: 11 original BLAST database inputs/index files.

Original bound files retain their historical strings, including local path names and documentary decision references. Editing those files changes their hashes. The restoration manifest itself contains relative destinations and archive names only.

## What this establishes, and what it does not

The restoration test verifies exact package bytes. The frozen Python source consists of the 19 files supplied separately with their original freeze hashes. The Rust evaluation used commit `834f4c288e9c163f316f3668f0219004d0adc6da`, not the later public release or macOS-parity commit. Both must be kept distinct when comparing results.

A fully portable end-to-end rerun still needs:

1. The original Python native runtimes (skani 0.3.2, mlst 2.32.2, BLAST 2.16.0+ and any2fasta) or a declared, validated replacement. Executables are not redistributed by this helper. Their original hashes and configuration are supplied elsewhere in the evidence package.
2. A documented path adapter or the original runtime layout. The Python MLST configuration binds executable/database paths and hashes; restoring its data does not rewrite that configuration or silently approve a new configuration hash.
3. The exact execution manifests, masks, frozen truth tables and query-derivative recipes. Downloading the 4,857 public assemblies alone does not specify the 5,803 native and 9,306 masked evaluations. Reference-only chromosome queries and synthetic scope/fragmentation conditions must use their recorded recipes and input hashes.
4. The original Rust source and compatible build/runtime environment. Recompilation is not expected to reproduce executable bytes across platforms; the scientific parity checks and saved evidence remain relevant.
5. Independent execution of that assembled workflow. No biological inference, classification or performance test is part of package restoration.

Separately, the phylogenetic reproducibility packet records missing historical chromosome-extraction commands, the unavailable original Mashtree wrapper bytes, the lack of an exact Task11 environment lock, and unrecorded verbatim final pruning invocations. Those gaps are not repaired by classifier-package restoration. Exact analysed sequences, preserved alignments, original trees and final trees support a narrower, explicit reproducibility claim.

## Helper checks

```sh
python3 -B -m unittest discover -s . -p test_restore_packages.py -v
```

Tests cover successful restoration, existing-output protection, traversal and absolute paths, link headers, corrupt archive/member hashes, missing members, duplicate targets and original manifest identity. The preparation worker also tested all 249 real targets against the supplied archives and reread every resulting file; see the retained restoration receipt.
