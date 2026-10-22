# Enterobacter genomic reassessment: analysis materials

Analysis materials for *Genomic reassessment of Enterobacter species and subspecies, with proposals for new subspecies combinations*, by Monique Clifton, Jocelyn Hume, Michael B. Hall and Leah W. Roberts.

This repository contains the saved measurements, trees, sequence manifests and analysis scripts behind the paper. To classify a new assembly, use [Pterobacter](https://github.com/mbhall88/pterobacter) and its [user guide](https://mbhall88.github.io/pterobacter/).

## Rebuild the figures

The four main figures and fourteen supplementary figures can be regenerated from the included measurements. This does not require the original HPC workspace, a database download or a new phylogenetic analysis.

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r environments/figures-requirements.txt
python workflow/figures/reproduce.py --output rebuilt-figures
```

The output directory must be new. Figures are written to `rebuilt-figures/outputs/figures/`, with commands' output in `logs/` and validation in `checks/`. The tested environment is recorded in `environments/`.

All 18 PNGs reproduced byte for byte in the recorded check. The tree display also checks all 1,785 pairwise path lengths and the original unrooted splits, branch lengths and support values. These checks establish display reconstruction, not a repeat of sequence alignment or tree inference. PDF metadata and SVG identifiers can differ between runs even when the graphics are unchanged.

## Rebuild the supplementary tables

```sh
python3 workflow/tables/rebuild_tables.py --output rebuilt-tables
```

This standard-library workflow reconstructs all 23 numbered TSVs exactly from included upstream measurements and metadata. Eleven tables repeat joins, counts or selected-field views; twelve repeat the original copy of an earlier saved measurement/register table. It does not rerun their biological analyses. Inputs and the original transformations are documented in [the table workflow](workflow/tables/README.md).

## Find the evidence

| Location | Contents |
| --- | --- |
| `workflow/figures/` | Portable, tested figure reconstruction |
| `data/figure-seeds/` | Saved display labels and summary inputs |
| `evidence/` | Selected analysis scripts, commands, manifests, measurements and final trees, retaining their original relative organisation |
| `metadata/assignments.tsv` | Accession-level reviewed assignment catalogue; its 1,098 entries are not the 1,070-genome reporting denominator |
| `expected/` | Approved figure files and selected reference summaries for comparison |
| `provenance/` | Original and exported hashes, documented adaptations and reproduction checks |
| `docs/reproduction.md` | What can be reconstructed and the remaining historical limitations |

The numbered supplementary tables and indexed workbook form a separate dataset. Their exact scientific files are preserved; companion source mappings explain legacy private discussion links and locate the public analysis evidence. Table S20 is a dated nomenclatural working register, not the final authority for the paper's taxonomic descriptions.

## Software and reference versions

The released application is [Pterobacter 0.1.0](https://github.com/mbhall88/pterobacter/releases/tag/0.1.0), source commit `6e4376e5d417439ffd593aa1f8d94b0aff9c65d8`. The frozen Rust evaluation used commit `834f4c288e9c163f316f3668f0219004d0adc6da`; these versions must not be conflated. Figure 4 reports the saved Python-baseline comparison described in the paper, not a fresh validation of the public software release.

The Pterobacter reference database has its own versioned record: [2026-10-15, DOI 10.6084/m9.figshare.34640883.v1](https://doi.org/10.6084/m9.figshare.34640883.v1). That DOI identifies the software reference package, not all the analysis data in this repository.

## Interpreting the names

Assignments record this study's reviewed taxonomic treatment. Proposed combinations and synonyms do not become validly published names merely by appearing in this repository. See the paper for the evidence, qualifications and nomenclatural proposals.

## Verify the downloaded repository

Run `python3 verify_export.py` before creating environments or reconstruction outputs inside this directory. It checks every listed file and rejects extra files except Git metadata. The [evidence index](provenance/evidence-index.json) supplies exported paths and hashes; the [archive catalogue](provenance/archive-catalogue.json) identifies the separate large inputs. See [content and rights](docs/content-and-rights.md) for provenance exceptions and approved terms.

## Data records

The [supplementary tables and workbook](https://doi.org/10.6084/m9.figshare.34742010) and the [exact analysis inputs and frozen evidence](https://doi.org/10.6084/m9.figshare.34742013) have separate Figshare records. The archive catalogue records the exact files and hashes. The repository source is available at [Enterobacter-analysis](https://github.com/mbhall88/Enterobacter-analysis).

The authors' code is MIT licensed; their original data, figures and text are CC BY 4.0. Third-party material retains its own terms; see [content and rights](docs/content-and-rights.md).
