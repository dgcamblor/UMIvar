---
tags:
  - pipeline
  - phd
description: Particularies of ctDNA analysis that a standard methylation pipeline (WGBS) should account for.
---

The steps of ctDNA WGBS analysis pipeline mostly follow the basics of the [[Bisulfite conversion processing pipeline]], while needing to account for the particularities of cfDNA and ctDNA detection. There are no apparent differences in analysis strategies between ctDNA methylation and conventional data analysis [@luoLiquidBiopsyMethylation2021; @huangCellFreeDNAMethylation2019] A good review on ctDNA methylation methods is: [@huangCellFreeDNAMethylation2019].

There are three essential approximations:

- Studying cfDNA methylation patterns as a whole, and establishing differential methylation via [[DM analysis]]. For example: [[ctDNA methylation studies — methodologies and results#^af677d]]. This approximation analyzes data as if it were from tissue.
- Deconvolution of cell types.
- Detecting ctDNA-specific patterns via deconvolution algorithms.

## Simple analysis

This approach studies cfDNA as a whole, treating the analysis as if it were from tissue. Thus, the standard methods apply.

To reduce the background signal, one commonly used strategy is to determine the methylation profile of tumor-free PBMCs as a negative control. By comparing DMRs between **cancer cfDNA** and **healthy cfDNA** to those between cancer cfDNA and PBMC genomic DNA, the shared regions are considered to be tumor-specific DMRs.

## Cell type deconvolution

The deconvolution into cell types can provide valuable information to assess the origin of the cfDNA (vs. ctDNA).

![[Cell type deconvolution methods#Epigenomic deconvolution]]

