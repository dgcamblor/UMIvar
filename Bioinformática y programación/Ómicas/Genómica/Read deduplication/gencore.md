---
githubUrl: https://github.com/OpenGene/gencore
paperUrl: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3280-9
dedupStrategy: Consensus
---

## Installation

```bash
conda install -c bioconda gencore
```

## Input

The input of gencore are paired-end sequencing data. Gencore accepts a sorted BAM, with its corresponding reference fasta. UMIs must be inside the read query names, as per the result of [fastp](https://github.com/OpenGene/fastp).

```
@NS500713:64:HFKJJBGXY:1:11101:1675:1101:AAAAAAAA
```

UMI appended at the end, separated by `:`, 8th field.

Note that gencore performs consensus with or without UMIs, although they are preferred. 
## Basic usage

```bash
gencore -i input.sorted.bam -o output.bam -r hg19.fasta -b test.bed
```