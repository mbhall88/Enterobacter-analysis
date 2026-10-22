# Reproducing the analyses

## Three different operations

1. **Redraw the published figures.** `workflow/figures/reproduce.py` reads included measurements, counts and trees. It recomputes the displayed aggregations and draws the figures. This path has been tested in a separate directory with access to the original research workspace and network denied.
2. **Recalculate summaries from saved evidence.** The selected scripts and per-genome/per-pair data in `evidence/` record the computations behind the manuscript tables and interpretations. The figure runner tests the aggregations that it uses. The table runner reconstructs all 23 supplementary TSVs exactly from upstream products, including the 12 original copy operations; neither claims to rerun every upstream scientific summariser.
3. **Repeat the biological computations.** ANI, genome alignment, marker detection and tree inference require the exact sequence scope, program versions and parameters. The historical command records describe those runs. A complete rerun has not been performed for this export.

## Paths and provenance

Files under `evidence/` retain their original relative organisation so that a saved result can be traced to its scripts and parameters. Except for the three hash-bound configurations listed in [content and rights](content-and-rights.md), historical absolute workspace prefixes have been replaced by `${PROJECT_ROOT}`, `${EVALUATION_ROOT}` and `${USER_HOME}`. These substitutions are documented in `provenance/analysis-export-map.json`; the file records both original and exported SHA256 values. They do not alter measured ANI, assignments, topology or sequence identifiers.

The archived commands and scripts are historical evidence. They are not all portable executable entry points: some refer to former job folders, environment locations or intermediate files. Review the command's inputs before adapting it for a rerun. The separate figure workflow supplies its own portable path adapter; its source changes are recorded in `provenance/figure-adaptations.patch`.

The original student archive, clinical validation cohort, author drafts and meeting discussions are outside this export. A historical path or discussion URL in a manifest is a provenance pointer, not a promise that the private working directory or conversation is publicly available. Companion scientific summaries and source mappings supply the public route to the relevant evidence.

## Sequence scope matters

Downloading the same assembly accession does not necessarily recreate the analysed input. Some analyses used an extracted chromosome; others used all assembly records. Exact accession.version, input scope and sequence hashes must be matched before a result is compared. In particular, the original extraction commands for the student's 1,063 chromosomes are not available. The analysed bytes must be preserved rather than substituting a newly downloaded whole assembly under the same identifier.

The main global inference initially used 1,090 genomes. The later 1,089-tip comparison and 1,088-tip sensitivity are distinct outputs. Five changed initial inputs were recovered at their original hashes; one originally pinned Mashtree wrapper was not found byte for byte in the available history. Recorded commands, seeds and final outputs remain available. Historical pruning commands were not all retained, so an independently implemented reconstruction must be identified as such.

An exact environment export is available for the Task17 Panaroo/core-gene work. The original Task11 environment directory is empty; surviving tool/version and command records are retained instead. Do not describe a fresh environment resolved from those versions as the original lockfile.

## Validation boundaries

The saved figure check verifies 18 PNGs, the two dataset/reference accounting tables and the type-tree display geometry. It does not establish that alternative tool versions, fresh database queries or new sequence extractions would give identical biological results. The exact input/output manifests and explicit historical gaps are part of the reproducibility record.

## Exact archived inputs

The [archive catalogue](../provenance/archive-catalogue.json) identifies the four final sequence, alignment/tree and frozen-package archives. Its completed readback receipts supersede the build-time `readback_pending` values in the original inventories. The [analysis-data record](https://doi.org/10.6084/m9.figshare.34742013) supplies these four archives. The frozen-package restoration workflow checks all 249 restored files and four original binding manifests; it does not rerun classification.
