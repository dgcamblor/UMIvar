---
githubUrl: https://github.com/BiodataAnalysisGroup/UMIc
paperUrl: https://pubmed.ncbi.nlm.nih.gov/34122513/
dedupStrategy: Consensus
---

## Installation

```r
# Dependencies
install.packages(c("tidyverse", "data.table", "stringdist", "pryr"))
BiocManager::install(c("Biostrings", "ShortRead"))
```

```bash
git clone https://github.com/BiodataAnalysisGroup/UMIc
```

Using UMIc requires providing the main script (`UMIsProject.R`) the corresponding input parameters.

## Input

Input data must be provided in FASTQ files, with the UMI assumed to be placed at the beginning of each sequence. Works mainly with **folders**. The input folder must contain the paired-end reads like `_R1.fastq.gz`.

## Parameters

Parameters must be selected in `UMIsProject.R`. For the [[Comparative analysis of UMI deduplication software in variant calling error correction (planning)]]:

```
pairedData <- TRUE
UMIlocation <- "R1 & R2"
UMIlength <- 6
sequenceLength <- 150
countsCutoff <- 1
UMIdistance <- 4
sequenceDistance <- 100
inputsFolder <- "/home/dgonzalez/dedup_bm/data/SRR13200987/UMIc"
outputsFolder <- "/home/dgonzalez/dedup_bm/data/SRR13200987/UMIc"
```

## Method

> - paired end libraries and UMI in both Read1 and Read2, each read of both files is separated in two parts, UMI and sequence. A new combined UMI12 is constructed by the union of the two UMIs for each sequence. The sequences of the Read1 file are grouped by the UMI12 and their IDs are used to find their corresponding Read2 sequences.

The same as the [[UMI-tools]] command `umi_tools extract` with different UMIs in R1 and R2.