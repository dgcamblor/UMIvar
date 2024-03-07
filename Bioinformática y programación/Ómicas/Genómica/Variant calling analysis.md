---
tags:
  - pipeline
---

## Oncoprint

An **Oncoprint** or **Waterfall plot** is a plot to represent a set of mutations in a set of samples. A color code can be used to designate the type of mutation, which needs to establish a mutation type hierarchy (to deal with multiple mutations in the same gene). 

Generally, the x-axis represents the samples and the y-axis represents the genes. The height of the bar represents the number of mutations in the gene in the sample. The color of the bar represents the type of mutation.

## dN/dS calculation

The dN/dS ratio compares the rate of non-synonymous substitutions (dN) to the rate of synonymous substitutions (dS) in protein-coding regions.

-  A dN/dS ratio greater than 1 suggests positive selection promoting change.
- A ratio of 1 indicates neutrality.
- A ratio less than 1 indicates purifying selection suppressing protein change

Available packages:

- [dndscv](https://github.com/im3sanger/dndscv)

## Tumor Mutational Burden (TMB)

The **Tumor Mutational Burden (TMB)** is a measure of the total number of mutations found in the tumor cells. It is calculated by counting the number of somatic mutations (SNVs and indels) per megabase of sequenced DNA. Two alternative approaches can be followed:

- Non-synonymous somatic mutations/Mb
- All somatic mutations/Mb

The key challenges of TMB adoption are the inconsistency of tumor mutational burden measurement among assays and the lack of a meaningful threshold for TMB classification.

## Mutation enrichment analysis

Mutation enrichment analysis consists on determining the top genes that contain somatic mutations in a certain cohort.

- [MutEnricher](https://github.com/asoltis/MutEnricher)

## SNV clonality

[Clonal evolution of chemotherapy-resistant urothelial carcinoma | Nature Genetics](http://dx.doi.org/10.1038/ng.3692)
## SNV nucleotide variant signatures

SNVs are partitioned into six mutation classes corresponding to six types of base pair substitution: C>A, C>G, C>T, T >A, T>C, T>G.