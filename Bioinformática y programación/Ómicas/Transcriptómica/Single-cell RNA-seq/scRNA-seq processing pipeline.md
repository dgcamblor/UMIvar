---
tags:
  - pipeline
nf-core: https://nf-co.re/scrnaseq/2.5.0
---
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

- HTSeq
- Cufflinks
- featureCounts

The output of this processing workflow should be a count matrix.

- Rows: transcripts.
- Columns: cells.