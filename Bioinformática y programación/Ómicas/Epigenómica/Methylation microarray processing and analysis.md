Although Illumina provides the GenomeStudio software to analyze Infinium methylation microarrays, the most versatility can be achieved using R packages from Bioconductor:

- SeSAMe
- RnBeads
- wateRmelon
- minfi

## Processing

From now on, we will be working with `minfi`, although similar steps are to be followed in other  packages. The input for `minfi` are [[IDAT format|IDAT]] files.

![[Pasted image 20240325195052.png]]

### Reading the data

The [[IDAT format|IDAT]] files are 

### Filtering problematic probes

- C/T SNPs
- Probes located in sex chromosomes
- Cross-reactive
- Multimapping probes
- Probes with detection p-values > 0.05 (0.01)

### Normalization

Normalization procedures currently available for correcting output include:

- Quantile normalization (QN)
- Beta-mixture quantile normalization (BMIQ)
- Subset-quantiles within microarray normalization (SWAN)

## References

- [Analysis of 450k data using minfi](https://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html)