The **mutual information (MI)** metric about two quantities is a measurement of how much the knowledge about one variable reduces the uncertainty about the other. The advantage to [[Correlation]] is that it detects all kinds of relationships, not just linear ones.

When MI is zero, the quantities are independent. There's no upper bound.

Requirements:

- **Discrete variables.** Mutual information can be performed directly.
- **Continuous variables.** May need discretization (for faster computing times or algorithm requirements).

## Python

Scikit-learn has two mutual information metrics in its `feature_selection` module: one for real-valued targets (`mutual_info_regression`) and one for categorical targets (`mutual_info_classif`).

## R

MI can be computed with the `infotheo` package.

## References

- [Mutual Information based Feature Selection Based for Ml | Medium](https://guhanesvar.medium.com/feature-selection-based-on-mutual-information-gain-for-classification-and-regression-d0f86ea5262a)