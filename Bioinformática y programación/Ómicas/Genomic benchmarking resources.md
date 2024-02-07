## Genome in a Bottle

[Genome in a Bottle | NIST](https://www.nist.gov/programs-projects/genome-bottle)

## SEQC2 Liquid Biopsy

The papers can be found at:

- [Ultra-deep sequencing data from a liquid biopsy proficiency study demonstrating analytic validity | Scientific Data](https://www.nature.com/articles/s41597-022-01276-8)
- [Evaluating the analytical validity of circulating tumor DNA sequencing assays for precision oncology | Nature Biotechnology](https://www.nature.com/articles/s41587-021-00857-z#code-availability)

360 datasets generated from two reference samples. Samples were:

- **Sample A:** cancer cell lines, genotyped for high confidence variants.
- **Sample B:** normal male cell line.
- **Sample D:** mix of A and B (1:4). Also termed "LBx-high". Median VAF: 1%. Majority above 0.5%
- **Sample E:** mix of A and B (1:24). Also termed "LBx-low". Median VAF: 0.2%. Majority above 0.1%
- **Sample F:** mix of A and B (1:124). 

### Burning Rock Project data

> BRP (Burning Rock Dx) assay provided the highest accuracy in the project, and we therefore chose the BRP datasets for our benchmark. 
> [Processing UMI Datasets at High Accuracy and Efficiency with the Sentieon ctDNA Analysis Pipeline | bioRxiv](https://www.biorxiv.org/content/10.1101/2022.06.03.494742v1)

The library preparation was performed with the Burning Rock HS UMI library preparation kit, which means it is a **hybridization-based assay**.

Sequencing data can be found at:

- LBx-high (Sample Df): [PRJNA677999 BURNING ROCK SampleDf - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA677999+BURNING+ROCK+SampleDf)
- LBx-low (Sample Ef): 
	- 25 ng input: [PRJNA677999 BURNING ROCK SampleEf 25ng - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA677999+BURNING+ROCK+SampleEf+25ng)
	- 10 ng input: [PRJNA677999 BURNING ROCK SampleEf 10ng - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA677999+BURNING+ROCK+SampleEf+10ng)

## SEQC2 Epigenomics (EpiQC)

Epigenomic analysis of seven well-characterized human cell lines used for genomic benchmarking by the GIAB Consortium.

[The SEQC2 epigenomics quality control (EpiQC) study | Genome Biology | Full Text](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02529-2#availability-of-data-and-materials)

All data sequenced is stored under the accession numbers [SRR13050956–SRR13051274](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA200694). They also point to the PRJNA646948 BioProject. Methylation [[bedGraph]] files are available in GEO: [GSE186383](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173688).

They also perform a comparison of alignment and methylation calling pipelines.

## Synthetic datasets


## Variant simulation

- BamSurgeon
- UMIvar