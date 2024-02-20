DM analysis can be done at different levels:

- DMR  -> Differentially methylated regions 
- DML -> Differentially methylated loci

The differential methylation between two conditions (case vs. control) is performed using standard statistical tests. As the distribution of the methylation level among the study population is unknown, a [[nonparametric test]] is preferentially adopted in methylation studies. For this purpose, [[limma]] can be used [@huangCellFreeDNAMethylation2019]. Additionally, for multiple-group comparison, ANOVA should be performed.

Prediction using CpG sites can be performed using [[Lasso regression]], which penalizes and eliminates CpG sites that are not relevant to the model.

## Binning strategy

To increase the statistical power, bins of a certain number of nucleotides can be considered for the differential methylation analysis. CpGs next to each other are most frequently correlated, therefore it is sensible to group them together in bins, thereby reducing the number of hypothesis tests.

## Bisulfite sequencing

Differential methylation analysis for bisulfite sequencing can be performed with:

- [[methylKit]] -> Designed for RRBS, but also works with WGBS.