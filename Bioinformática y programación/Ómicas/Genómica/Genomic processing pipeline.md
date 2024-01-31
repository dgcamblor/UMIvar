---
tags:
  - pipeline
  - phd
---
Some important reviews that assess good practices in variant calling:

- Best practices for variant calling in clinical sequencing [@koboldtBestPracticesVariant2020].

## FASTQ Quality Control and preprocessing

### Initial quality control

![[FASTQ quality control#Software for FASTQ QC]]

### Preprocessing

1. **Extraction of UMIs (if they are available)**.
	- [[UMI-tools]] (extract) -> UMI in the header, separated by "\_"

2. **Adapter (and primer) trimming.** Trimming can be performed using a known adapter, or an adapter list. Alternatively, the primers can be clipped after mapping (see [[Primer removal alternatives]])
	- [[cutadapt]] 

4. **Filtering out reads with low quality.** A commonly used threshold is `20` (Q20), corresponding to a call accuracy of 99%. Aiming for a quality of `30` (Q30) is ideal for variant calling in critical settings such as clinical research.
	- [[prinseq]] 

### Quality control after preprocessing

The same QC software used in [[#Initial quality control]] is applied to check that the FASTQ preprocessing steps have been performed correctly.

## Mapping and BAM file processing

### Mapping to the reference genome

![[Genome mapping#Standard aligners]]

### BAM file preprocessing

The output of the mapping is a SAM file.

1. **Converting the SAM file to a BAM file.** This is done with [[Samtools]]: `samtools view -S -h -b file.sam > file.bam`

2. **Filtering low mapping quality reads.** The `view` command can also be used to filter reads that do not meet a certain minimum mapping quality with the `-q` parameter. A common practice is to filter reads with a MAPQ value below `30`, which is used as a cutoff for retaining high quality mappings.

3. **Sorting and indexing.**

#### Deduplication (or not)

![[Read deduplication]]

#### Posterior BAM processing

1. **Cleaning the BAM file.** This can be done with [[Picard#CleanSam]].

2. [[Base Quality Score Recalibration]].

#### Specific processing

- Lofreq requires indel qualities. They are not added with [[Base Quality Score Recalibration]], so `lofreq indelqual` should be used (`--dindel` for Illumina data, `--uniform` for non Illumina).

## Variant calling

### Variant calling proper

![[Variant calling software]]

### Filtering false positives

- Mutect2 calls can be filtered with [[GATK#FilterMutectCalls]].

- In the final VCF, variants can be selected that meet a minimum criteria.
	- **Minimum AF**: 

### Normalization

Using [[bcftools]]: `bcftools norm`.

### Variant calling analysis

![[Variant calling analysis]]
## CNV determination

## Variant annotation