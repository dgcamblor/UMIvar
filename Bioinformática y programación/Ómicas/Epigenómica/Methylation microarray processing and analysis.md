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

From now on, we will be working with `minfi`, although similar steps are to be followed in other packages. This package is highly compatible with other packages. The input for `minfi` are [[IDAT format|IDAT]] files.

![[Pasted image 20240325195052.png]]

### Reading the data

The [[IDAT format|IDAT]] files are read by `minfi` using their filenames or the directory path. The user can also read from a sample sheet.

The result is a `RGChannelSet` object, the initial `minfi` object that contains the raw intensities in the green and red channels. This object also contains the intensities of the internal control probes.

```
pheno_data <- pData(RGSet) 
pheno_data[,1:6]
```

### Filtering problematic probes

Several probes are known to be problematic due to a myriad of reasons, and it is advisable to remove them from the analysis. 

- C/T SNPs -> Probes with a common SNP at the CpG site, which can lead to a false positive signal. They can be identified with the `getSNPs()` function.

- Probes located in sex chromosomes -> They can trigger false positives in the differential methylation analysis. They can be identified with the `getSex()` function.

- Cross-reactive/multimapping probes -> Some of the probes can map to multiple locations in the genome, leading to confounding signals. Cross reactive probes have been documented in the literature ([Chen et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23314698/); [Pidsley et al., 2016](https://pubmed.ncbi.nlm.nih.gov/23314698/)). Some packages like `maxprobes` allow for easy removal in `minfi` objects.

- Probes with detection p-values > 0.01. The detection p-value is a measure of the confidence in the signal of a probe, comparing the total signal at the probe (Methylated + Unmethylated) with the background signal (estimated from the negative control probes). Probes with detection p-values > 0.01 (or 0.05) are considered unreliable.

### Normalization

Normalization is an important step in the analysis of the methylation microarray data to correct for technical biases:

- Background signal (within-array): Background signal in the arrays
- Color bias (within-array): Differences  in the intensity measurement fidelity between the two dyes.
- Probe-type bias (within-array) (most important within-array bias): Corrects for the differences in the signal between the two probe types (Infinium I and II). Type I probes are known to have a higher dynamic range than type II probes. Infinium II assays use the same bead to measure both the methylated and unmethylated signal, so the measurement of one signal can affect the measurement of the other signal.
- Batch/array-specific effects (between-arrays): Corrects for the differences in the signal between the arrays. These differences can be due to the batch in which the arrays were processed, the slide, or the position of the array in the slide.

Normalization procedures currently available for correcting output include:

- Illumina normalization (GenomeStudio). Background correction and control normalization. Implemented in the `preprocessIllumina` function (reverse-engineered from GenomeStudio).

- Quantile normalization (QN)

- Beta-mixture quantile normalization (BMIQ)

- Subset-quantiles within microarray normalization (SWAN). Within-array normalization method that corrects for probe-type bias (Infinium I and II).

- Peak-based correction (PBC)

- Functional normalization (Funnorm)

- Normal-exponential convolution using out-of-band probes (Noob)

- Single-sample noov (SSnoob)

For 450K data BMIQ appears to be the most effective normalization method for probe-type bias.

Some bibliography:

- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10008016/

## References

- https://academic.oup.com/bib/article/15/6/929/179607 -> Filtering and normalization
- [Analysis of 450k data using minfi](https://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html)
- [Maternal obesity and gestational diabetes reprogram the methylome of offspring beyond birth by inducing epigenetic signatures in metabolic and developmental pathways - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9985842/)
- https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html