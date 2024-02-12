---
tags:
  - pipeline
  - phd
nf-core: https://nf-co.re/methylseq
---
> To assign:
> - [GitHub - zhou-lab/biscuit: BISulfite-seq CUI Toolkit](https://github.com/zhou-lab/biscuit). DNA methylation and mutation calling.

The ENCODE Project gathers a series of WGBS data and processing standards at: [Whole-Genome Bisulfite Sequencing Data Standards and Processing Pipeline – ENCODE](https://www.encodeproject.org/data-standards/wgbs/). Some important standards are:

- Each replicate should have 30X [coverage](https://www.encodeproject.org/data-standards/terms/#read-depth).
- The C to T conversion rate should be ≥98%.
- The CpG quantification should have a [Pearson correlation](https://www.encodeproject.org/data-standards/terms/#concordance) of ≥0.8 for sites with ≥10X coverage.

## FASTQ preprocessing and quality control

The quality of the FASTQ reads must be initially assessed using the habitual software.

![[FASTQ quality control#Software for FASTQ QC]]



- Adapter trimming: Trimmomatic, Trim Galore, cutadapt.
- Data visualization: IGV, Methylation plotter, WBSA.

## Alignment and quality control

### Alignment

![[Alignment in epigenomics processing#Bisulfite conversion-based methods]]

### BAM quality control

Once aligned, control quality is very important in order to avoid miscalled C-T conversions [@huangCellFreeDNAMethylation2019].

> [!danger]
> Incomplete bisulfite conversion is a problem, because it causes false positives. 

To perform quality control on bisulfite conversion rate, there are two main approaches:

- The use of **spike-in control** sequences. These are sequences that are either completely unmethylated (lambda phage DNA) or completely methylated (pUC19 plasmid). Rates higher than 98.5% ensure the absence of bias [@gongAnalysisPerformanceAssessment2022].
- Studying the post-bisulfite methylation patterns both in CpG (normal) and in CHG and CHH contexts (where methylation is not encountered in mammals).
	- CpG methylation levels must fall in the expected 45-65% range [@fooxSEQC2EpigenomicsQuality2021].
	- CHG and CHH methylation rates should be close to the expected 0% range (which would mean a 100% conversion efficiency). 

Some additional quality control steps are:

- Filtering out C/T SNPs is highly recommended.
- Removal of PCR duplicates (by genome position) -> [[GATK#Picard#MarkDuplicates]]

## Methylation calling and analysis

At each CpG site, the methylation levels are calculated by looking at all the reads that cover that position.

$$
Methylation\ levels\ at\ CpG\ site = C / C+T
$$

Software choices: 

- [[Bismark]].
- [[MethylDackel]]. Groups cytosine into CpG, CHG and CHH contexts. Outputs a [[bedGraph]].
- [[methylKit]]. This R package includes a myriad of functions for clustering, sample quality visualization, [[Differential methylation (DM) analysis]]

If using a [[Directional vs. non-directional bisulfite-converted libraries#Directional bisulfite sequencing|Directional bisulfite sequencing]] method,`MethylDacke mergeContext` can be used to produce one value per CpG dinucleotide.

Additional reference: [@kimMsPIPEPipelineAnalysis2022].

Finally, the called methylation is subject to [[Differential methylation (DM) analysis]].

![[Differential methylation (DM) analysis]]

## Annotation

The genomic sites/regions of interest, such as individual CpGs, differentially methylated CpGs or regions, etc., must be annotated in order to gain functional insights. For this purpose, one can use:

- [[annotatr]] (R package)
- [[methylKit]] 
- methylGSA. Approaches: ORA, GSEA.

## Data visualization

- [[methylKit]]

## References

- Huang, Jinyong, y Liang Wang. «Cell-Free DNA Methylation Profiling Analysis—Technologies and Bioinformatics». _Cancers_ 11, n.o 11 (6 de noviembre de 2019): 1741. [https://doi.org/10.3390/cancers11111741](https://doi.org/10.3390/cancers11111741).
- Gong, Ting, Heather Borgard, Zao Zhang, Shaoqiu Chen, Zitong Gao, y Youping Deng. «Analysis and performance assessment of the whole genome bisulfite sequencing data workflow: currently available tools and a practical guide to advance DNA methylation studies». _Small methods_ 6, n.o 3 (marzo de 2022): e2101251. [https://doi.org/10.1002/smtd.202101251](https://doi.org/10.1002/smtd.202101251).