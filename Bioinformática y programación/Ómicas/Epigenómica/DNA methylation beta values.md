
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

For **microarray data**, the average methylation level of all CpG sites in a cluster is used to represent methylation level of that cluster. A cluster’s methylation level is marked as “not available” (NA) if less than half of its CpG sites have methylation measurements.

For **WGBS data**, the methylation level of a CpG cluster is calculated as the ratio between the number of methylated cytosines and the total number of cytosines within the cluster. However, if the total number of cytosines in the reads aligned to a CpG cluster is less than a given threshold (30 as used in the paper), the methylation level of this cluster is considered as NA.