---
tags:
  - pipeline
  - phd
---

## FASTQ Quality Control and preprocessing

### Initial quality control

![[FASTQ quality control]]

### Preprocessing

1. **Extraction of UMIs (if they are available)**.
	- [[UMI-tools]] (extract)

2. **Adapter trimming.** Trimming can be performed using a known adapter, or an adapter list.
	- [[cutadapt]] 

3. **Filtering out reads with low quality.** A commonly used threshold is `20` (Q20), corresponding to a call accuracy of 99%. Aiming for a quality of `30` (Q30) is ideal for variant calling in critical settings such as clinical research.
	- [[prinseq]] 

## Mapping and BAM file processing

- [[Base Quality Score Recalibration]].

## Variant calling

## Variant annotation