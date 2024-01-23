---
paperStatus: planning
---

NOTE: Search for variant caller benchmarking papers. See what pipelines do they use.
## Approach

- Comparing the effect of a series of open source UMI deduplication softwares in the endstream variant calling process:
	- [[UMI-tools]]
	- [[fgbio]]
	- [[UMICollapse]]
	- [[UMIc]] (only UMIvar?)
	- [[gencore]]
- Tools such as zUMIs and umis (UMICollapse paper) are more RNA-seq oriented and, thus, were discarded. [[AmpUMI]] was considered but finally discarded due to it only working in settings where R1 and R2 overlap. [[UMIc]] does not handle correctly data with deep coverage.

- Comparing the effect against a reference of deduplication: [[GATK]] MarkDuplicates ([MarkDuplicates (Picard) – GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)). 

- The comparison will be in terms of endstream **variant calling efficiency** (decrease in FP, increase in TP).

- The comparison will be done in two datasets:
	- SEQC2
	- UMIvar



## Deduplication software

``` dataview
table dedupStrategy, githubUrl, paperUrl
from "Bioinformática/Genómica/Read deduplication"
```


## Datasets

### SEQC2 Liquid Biopsy

[Ultra-deep sequencing data from a liquid biopsy proficiency study demonstrating analytic validity | Scientific Data](https://www.nature.com/articles/s41597-022-01276-8). BRP dataset was chosen because it was sequenced with UMI tags, and it proved to have the best results).

#### General characteristics
![[Genomic benchmarking resources#SEQC2 Liquid Biopsy]]

#### Specific data


### UMIvar
## Benchmarking references

- [Evaluating assembly and variant calling software for strain-resolved analysis of large DNA viruses](https://academic.oup.com/bib/article/22/3/bbaa123/5868070)
- [Best practices for benchmarking germline small-variant calls in human genomes](https://www.nature.com/articles/s41587-019-0054-x)
- [Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery - BMC Genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08365-3)
- [A review of somatic single nucleotide variant calling algorithms for next-generation sequencing data](https://www.sciencedirect.com/science/article/pii/S2001037017300946)
- [Gencore: an efficient tool to generate consensus reads for error suppressing and duplicate removing of NGS data - BMC Bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3280-9?ref=https://githubhelp.com)
- [biorxiv.org/content/10.1101/2022.06.03.494742v2.full.pdf](https://www.biorxiv.org/content/10.1101/2022.06.03.494742v2.full.pdf)
