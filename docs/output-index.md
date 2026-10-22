# Find analysis outputs

Start with the full-cohort trees below. Smaller trees used for figures and targeted follow-up analyses are listed separately. Files retain their original names and locations so their provenance remains traceable.

## Main full-cohort trees

These three files contain the matched 1,089-tip comparison. They include reference/context genomes and outgroups; this is not the 1,070-assembly reporting denominator. They are Newick text files: open a link, then use **Raw** or **Download raw file** to save the tree for a tree viewer.

| Method | Tree file | Interpretation |
| --- | --- | --- |
| BAC120 marker tree | [bac120_fasttree.nwk](../evidence/results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/full-with-20-seed-consensus-1089/trees/bac120_fasttree.nwk) | FastTree marker-protein tree; internal support values are SH-like local support. |
| Verticall distance tree | [verticall_fastme.nwk](../evidence/results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/full-with-20-seed-consensus-1089/trees/verticall_fastme.nwk) | FastME tree from Verticall distances; no bootstrap support was generated. |
| Mashtree consensus | [mashtree_consensus_20.nwk](../evidence/results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/full-with-20-seed-consensus-1089/trees/mashtree_consensus_20.nwk) | Consensus across 20 sketch seeds; seed agreement is not a sequence-site bootstrap. |

[Comparison settings and method definitions](../evidence/config/tree-comparison-full.json) describe the analysis. Other seed draws, rooting choices and sensitivity trees are retained in the catalogue below; they are not interchangeable with this comparison.

## Trees shown in the figures

| Output | Location |
| --- | --- |
| Figure 2: 35 type-anchor projection | [Three projected trees](../evidence/results/full-cohort/2026-09-22/task11-opus-full/comparison-1089/type-anchor-projection-35/trees); [Displayed strain labels](../data/figure-seeds/source-data/figure2-strain-labels.tsv) |
| Figures S1–S7: group-focused tree panels | [Newick files and figure source data](../expected/figure-data) — look for `figureS1-*.nwk` through `figureS7-*.nwk` |
| All four main and fourteen supplementary figures | [PNG, PDF and SVG figures](../expected/figures) |

The figure trees are projections/subsets for display, not additional full-genome tree inferences. Use the full-cohort files above when you need all sampled tips.

## Targeted trees and follow-up evidence

These folders contain retained analyses and sensitivities, including alternatives that were not adopted as a primary result. An alternative tree or topology constraint is not itself a taxonomic conclusion.

| Question | Evidence folder |
| --- | --- |
| bugandensis / nematophilus local trees | [Browse evidence](../evidence/results/ani-cases/2026-09-27/task12-astra-nematophilus-expanded) |
| asburiae / dykesii local trees | [Browse evidence](../evidence/results/ani-cases/2026-09-28/task12-astra-asburiae-local) |
| Neighbouring genera and genus context | [Browse evidence](../evidence/results/ani-cases/2026-10-01/task12-astra-genus-tree) |
| Group A: hoffmannii / intestinihominis / nototheniae populations | [Browse evidence](../evidence/results/ani-cases/2026-10-02/task12-astra-group-a-populations) |
| Conflicting marker/core-genome placements | [Browse evidence](../evidence/results/conflict-resolution/2026-09-23/task17-opus55-resolve) |

**[Search all 418 saved tree files](tree-files.tsv).** This tab-separated catalogue lists the filename, category, full repository path, clickable GitHub URL, size and hash. Search it for a group, method, accession or filename. Constraint files are marked separately; some test files contain multiple candidate topologies. The file count is not a count of independent analyses.

## Tables, assignments and measurements

| Material | Location |
| --- | --- |
| Tables S01–S23 (individual TSVs) | [Open](../expected/tables) |
| Accession-level reviewed assignments | [Open](../metadata/assignments.tsv) |
| Figure source tables and measurements | [Open](../expected/figure-data) |
| Supplementary-table reconstruction and source mappings | [Open](../workflow/tables/README.md) |
| Indexed Excel workbook and deposited tables | [Supplementary-table record](https://doi.org/10.6084/m9.figshare.34742010.v1) |

The assignment catalogue has 1,098 entries, including reference/context records, rather than the 1,070-genome descriptive denominator. Proposed names reflect the study's treatment; they do not establish valid publication of nomenclatural acts.

## Alignments, sequences and frozen reference packages

Large files are in the **[versioned analysis-data record on Figshare](https://doi.org/10.6084/m9.figshare.34742013.v1)**. Download only the archive you need:

| Archive | Compressed size | Contents |
| --- | --- | --- |
| `enterobacter-analysed-sequences-2026-10-22.tar.gz` | 1.717 GB | Exact analysed sequence files |
| `enterobacter-analysis-alignments-trees-2026-10-22.tar.gz` | 1.243 GB | Alignments, masks, trees and related saved analysis evidence |
| `enterobacter-context-sequences-2026-10-22.tar.gz` | 0.011 GB | Additional context sequences |
| `enterobacter-frozen-evaluation-packages-noticed-2026-10-22.tar.gz` | 0.124 GB | Frozen evaluation reference packages and required third-party notices |

**[Search the archive contents before downloading](archive-files.tsv).** The catalogue gives each member's exact path, size and SHA256. Archive-level hashes are in the [archive catalogue](../provenance/archive-catalogue.json).

For example, the full-cohort filtered BAC120 alignment is in the alignments archive at:

```text
results/full-cohort/2026-09-22/task11-opus-full/marker/filtered/alignment.faa
```

After downloading the archive, extract just that file:

```sh
tar -xzf enterobacter-analysis-alignments-trees-2026-10-22.tar.gz \
  results/full-cohort/2026-09-22/task11-opus-full/marker/filtered/alignment.faa
```

The alignment retains its original inference-stage membership. The matched 1,089-tip comparison trees were subsequently pruned; consult the [comparison settings](../evidence/config/tree-comparison-full.json) before assuming alignment rows and displayed tree tips are identical.

## Reproduce or trace an output

- [Rebuild the figures](../README.md#rebuild-the-figures).
- [Rebuild the supplementary tables](../README.md#rebuild-the-supplementary-tables).
- [Reconstruction scope and historical gaps](reproduction.md).
- [Complete evidence/provenance index](../provenance/evidence-index.json).

The public release of Pterobacter and its user database are separate from the frozen materials used for the paper's evaluation. See [software and reference versions](../README.md#software-and-reference-versions).
