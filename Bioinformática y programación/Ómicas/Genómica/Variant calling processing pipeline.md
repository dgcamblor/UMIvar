---
tags:
  - pipeline
  - phd
---

## FASTQ Quality Control and preprocessing

### Initial quality control


### Preprocessing

1. Extraction of UMIs (if they are available)
	- [[UMI-tools]] (extract)
	
2. Adapter trimming. Trimming can be performed using a known adapter, or an adapter list.
	- [[cutadapt]] 

3. Removing reads with low quality. A good minimum mean quality cutoff is `30`. 

## Mapping and BAM file processing

- [[Base Quality Score Recalibration]].

## Variant calling

## Variant annotation