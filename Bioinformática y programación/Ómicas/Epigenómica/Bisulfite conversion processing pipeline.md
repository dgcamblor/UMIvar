---
tags:
  - pipeline
  - phd
  - write
nf-core: https://nf-co.re/methylseq/2.6.0
---
## Alignment

Alignment to the reference genome depends on the strategy used.

- **Bisulfite conversion-based methods.** Requires special aligners, as the sequence itself has been altered by the bisulfite conversion and does not align directly to the reference genome.
- **Enrichment-based methods.** The standard aligners can be used.

## FASTQ preprocessing and quality control

The quality of the reads must be assessed using the habitual software.

- Quality control proper: FastQC
- Adapter trimming: Trimmomatic, Trim Galore, cutadapt.
- Data visualization: IGV, Methylation plotter, WBSA.

## Alignment and quality control

### Alignment

For bisulfite sequencing, bisulfite converted DNA does not align directly to the reference genome. Two algorithms are available: wild card and three-letter.

- Wild card algorithm. Both Cs and Ts map into Cs in the reference genome.
- Three-letter algorithm. Converst all Cs in the reference genome and the reads into Ts, so that standard aligners can be used.

Three-letter aligners seem to outperform wildcard aligners in running time and peak memory usage [@gongAnalysisPerformanceAssessment2022].

Software choices:

- Bismark
- BS Seeker 2
- BWA-Meth

The SEQC2 study and additional studies seem to suggest that BWA-METH has the highest uniquely mapped read rates and the lowest unmapped reads [@fooxSEQC2EpigenomicsQuality2021; @gongAnalysisPerformanceAssessment2022].

### Quality control

Once aligned, control quality is very important in order to avoid miscalled C-T conversions [@huangCellFreeDNAMethylation2019].

> [!danger]
> Incomplete bisulfite conversion is a problem, because it causes false positives. 
> - Spìke-in sequences with unmethylated Cs are generally used and analyzed to measure the bisulfite conversion rate. For example, one alternative is to use unmethylated lambda phage DNA.
> - Another approximation in whole genome sequencing would be to study the post-bisulfite methylation patterns in CHG and CHH contexts. Here, methylation rates should be close to the expected 0% range (which would mean a 100% conversion efficiency). [[Metilación del ADN#^dd9ff0]].
> Rates higher than 98.5% ensure the absence of bias [@gongAnalysisPerformanceAssessment2022].

Some quality control steps:

- Filtering out C/T SNPs is highly recommended.
- Removal of PCR duplicates (by genome position) -> [[GATK#Picard#MarkDuplicates]]

## Methylation calling and analysis

$$
Methylation\ levels\ at\ CpG\ site = C / C+T
$$
Software choices: 

- Bismark.
- [[MethylDackel]]. Groups cytosine into CpG, CHG and CHH contexts. Outputs a [[bedGraph]].

Additional reference: [@kimMsPIPEPipelineAnalysis2022].

Finally, the called methylation is subject to [[DM analysis]].

![[DM analysis]]

## Annotation

The genomic sites/regions of interest, such as individual CpGs, differentially methylated CpGs or regions, etc., must be annotated in order to gain functional insights. For this purpose, one can use:

- [[annotatr]] (R package)
- methylKit
- methylGSA. Approaches: ORA, GSEA.

## Data visualization

- methylKit

## References

- Huang, Jinyong, y Liang Wang. «Cell-Free DNA Methylation Profiling Analysis—Technologies and Bioinformatics». _Cancers_ 11, n.o 11 (6 de noviembre de 2019): 1741. [https://doi.org/10.3390/cancers11111741](https://doi.org/10.3390/cancers11111741).
- Gong, Ting, Heather Borgard, Zao Zhang, Shaoqiu Chen, Zitong Gao, y Youping Deng. «Analysis and performance assessment of the whole genome bisulfite sequencing data workflow: currently available tools and a practical guide to advance DNA methylation studies». _Small methods_ 6, n.o 3 (marzo de 2022): e2101251. [https://doi.org/10.1002/smtd.202101251](https://doi.org/10.1002/smtd.202101251).