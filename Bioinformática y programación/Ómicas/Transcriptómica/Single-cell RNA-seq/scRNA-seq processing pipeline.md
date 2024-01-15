---
tags:
  - pipeline
nf-core: https://nf-co.re/scrnaseq/2.5.0
---

A comparative analysis of scRNA-seq mapping tools can be found at: [Comparative analysis of common alignment tools for single-cell RNA sequencing | GigaScience | Oxford Academic](https://doi.org/10.1093%2Fgigascience%2Fgiac001).

## Sequencing QC

- RSeQC (RNA-seq specific)
- FastQC

Outputs can be summarized with MultiQC.

## Sequence preprocessing

## Mapping

- CellRanger (for Chromium 10x data)
- Standard aligners
	- Bowtie2
	- Tophat2
	- STAR

## Mapping QC

## Computing read counts

- Salmon
- Kallisto
- CellRanger
- HTSeq
- Cufflinks
- featureCounts

The output of this processing workflow should be a UMI count matrix. 

- Rows: transcripts.
- Columns: cells.

The `mtx` and `h5ad` formats are used to store gene expression matrices, with `h5ad` offering additional capabilities for storing metadata and annotations. The matrix is loaded in R, and subject to [[scRNA-seq analysis]].