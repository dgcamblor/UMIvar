
Two commonly used measures are used to report the methylation levels: Beta values and M values

## Beta values

Beta values are computed using the methylation (M) and unmethylation (U) signals:

$$
\beta = \frac{M}{M+U+100}
$$

DNA methylation beta values are continuous variables measured between 0 and 1.

- A number close to 0 indicates non methylation of the region.
- A number close to 1 indicates complete methylation of the region.

```text
probe  sample_A  Sample_B
cg_1      0.6       0.7
cg_2      0.2       0.3
cg_3      0.8       0.9
cg_4      0.2       0.9
cg_5      0.3       0.6
cg_6      0.1       0.4
```

## M values

M values are also computed using the methylation (M) and unmethylation (U) signals, and use the logit transformation:

$$
M = \log_2\left(\frac{M}{U}\right)
$$

M values are continuous variables that can take any value between $-\infty$ and $\infty$.

## Comparison

The paper [**Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587) compares the two methods and concludes that:

- Beta values are more intuitive and easier to interpret in a biological context.
- M values are more statistically valid and are more appropriate for statistical analysis. This is because the distribution of M values is more symmetric and homoscedastic than that of Beta values.

```text
For **microarray data**, the average methylation level of all CpG sites in a cluster is used to represent methylation level of that cluster. A cluster’s methylation level is marked as “not available” (NA) if less than half of its CpG sites have methylation measurements.

For **WGBS data**, the methylation level of a CpG cluster is calculated as the ratio between the number of methylated cytosines and the total number of cytosines within the cluster. However, if the total number of cytosines in the reads aligned to a CpG cluster is less than a given threshold (30 as used in the paper), the methylation level of this cluster is considered as NA.
```