---
tags:
  - pipeline
---

## Oncogenic-related analysis
### Oncoplot

A **Oncoplot** or **Waterfall plot** is a plot to represent a set of mutations in a set of samples. A color code can be used to designate the type of mutation, which needs to establish a mutation type hierarchy (to deal with multiple mutations in the same gene). 

Generally, the x-axis represents the samples and the y-axis represents the genes. The height of the bar represents the number of mutations in the gene in the sample. The color of the bar represents the type of mutation.

### dN/dS calculation

The dN/dS ratio compares the rate of non-synonymous substitutions (dN) to the rate of synonymous substitutions (dS) in protein-coding regions.

-  A dN/dS ratio greater than 1 suggests positive selection promoting change.
- A ratio of 1 indicates neutrality.
- A ratio less than 1 indicates purifying selection suppressing protein change

### Tumor Mutational Burden (TMB)

The **Tumor Mutational Burden (TMB)** is a measure of the total number of mutations found in the tumor cells. It is calculated by counting the number of somatic mutations (SNVs and indels) per megabase of sequenced DNA. Two alternative approaches can be followed:

- Non-synonymous somatic mutations/Mb
- All somatic mutations/Mb

The key challenges of TMB adoption are the inconsistency of tumor mutational burden measurement among assays and the lack of a meaningful threshold for TMB classification.