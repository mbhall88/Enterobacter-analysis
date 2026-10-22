# Interpreting the supplementary-table source records

This note explains the three Table S20 records that cite the internal taxonomy-review map. It draws on the final author-reviewed supplement, sections S06–S07, and the unchanged table entries. It does not reproduce private discussions or assign new labels.

## ECC-175 and the residual lineages

Table S20 record `P6A-UNI-03` (TSV row 39) retains ECC-175 as a separate unnamed Enterobacter lineage. Supplement S06, “Other boundaries and unnamed lineages”, also leaves the two GTDB cloacae_P genomes unnamed and retains RHB15-C17 as Lelliottia nimipressuralis context. RHBSTW-00994 is treated as an Enterobacter lineage with uncertainty about its genus placement in the broader trees. These treatments apply to the specified accessions; GTDB suffixes are not formal taxonomic names. Table S01 contains the accession-level assignments and qualifications.

## The species-only xiangfangensis assignment

Table S20 record `P6A-UNA-04` (TSV row 44) assigns GCA_039734725.2 to Enterobacter xiangfangensis without a subspecies. Supplement S06 distinguishes that reviewed treatment from a separate nearest-type sensitivity analysis, which places the genome in a descriptive steigerwaltii bin. That bin does not supply a subspecies assignment. The within/between-subspecies comparison in Figure S3 therefore uses the 437 assigned genomes and excludes this species-only genome. Figure S2 retains the broader population comparison.

## Informal genomovars within the hoffmannii subspecies

Table S20 record `P6A-INF-01` (TSV row 45) records 25 numbered genomovars within the proposed Enterobacter xiangfangensis subsp. hoffmannii. Supplement S07 defines these by complete-linkage clustering at inclusive ANI 99.5%, requiring every within-group pair to meet the threshold. The 173-genome definition set includes one type genome absent from the study collection. Ten groups contain multiple sampled members and fifteen are singletons. Eleven pairs crossing groups also reach the threshold, so the groups are an explicit clustering convention rather than perfectly separated ANI bands. Memberships and cross-group pairs are in Tables S09–S10; Figure S10 displays the saved ANI matrix. These informal groups are below the ranks regulated by the ICNP and are not proposed as formally named taxa.

## The candidate topology files

The 41 small tree files accompanying the Table S05 source map preserve the topologies tested in focused marker/core-gene analyses. A file's inclusion does not mean its topology was accepted. In particular, the expanded ludwigii constraint labelled as the marker hypothesis placed ludwigii outside all sampled in-genus lineages, unlike its actual marker-tree position. Its rejection therefore cannot be used as evidence against that marker placement. Supplement S04 and Table S05 identify the tested hypotheses, conditional likelihood comparisons and excluded interpretations. The trio's non-converged searches and uncertain internal branching remain qualifications; no additional test was run for this export.

## Source versions

The source map records exact saved-file hashes and locations within the proposed public repository. LPSN links come from the canonical-link metadata in the saved catalogue pages. Those URLs identify the public service; the current page can differ from the recorded snapshot. The saved HTML pages are not redistributed in this export. Full likelihood reconstruction additionally needs the alignment, software, parameters and commands listed in the analysis workflow; a topology file alone is not a complete rerun package.
