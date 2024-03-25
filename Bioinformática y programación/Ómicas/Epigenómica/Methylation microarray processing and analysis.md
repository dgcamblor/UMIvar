Although Illumina provides the GenomeStudio software to analyze Infinium methylation microarrays, the most versatility can be achieved using R packages from Bioconductor:

- SeSAMe
- RnBeads
- wateRmelon
- minfi

Some microarray concepts for their bioinformatic analysis are:

- Array -> One sample.
- Slide -> Physical slide containing 12 arrays (6 x 2 grid).
- Plate -> Physical plate containing at most 8 slides (96 arrays). The plate determines the batch

## Processing

From now on, we will be working with `minfi`, although similar steps are to be followed in other  packages. The input for `minfi` are [[IDAT format|IDAT]] files.

![[Pasted image 20240325195052.png]]

### Reading the data

The [[IDAT format|IDAT]] files are read by `minfi` using their filenames or the directory path. The user can also read from a sample sheet.

The result is a `RGChannelSet` object, the initial `minfi` object that contains the raw intensities in the green and red channels. This object also contains the intensities of the internal control probes.

```
pheno_data <- pData(RGSet) 
pheno_data[,1:6]
```

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