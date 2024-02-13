---
aliases:
  - Adapter removal
---

## Software for adapter trimming

The main adapter trimming software is:

- cutadapt
- Trimmomatic
- Trim Galore! -> A wrapper of [[cutadapt]] and [[FastQC]], that automatically performes adapter trimming and quality control.  It has additional functionality to removed biased methylation in [[RRBS]].

## Particularities

### Adapter trimming in Illumina

>Illumina FASTQ file generation pipelines include an adapter trimming option for the removal of adapter sequences from the 3’ ends of reads. Libraries prepared with Illumina library prep kits require adapter trimming only on the 3’ ends of reads, because adapter sequences are not found on the 5’ ends.
>To understand why adapter sequences are found only on the 3’ ends of the reads, it helps first to understand where the sequencing primers anneal to the library template on a flow cell.

![[Pasted image 20240207133553.png]]