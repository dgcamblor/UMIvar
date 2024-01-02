Data imputation is not always needed. Some algorithms handle missing data.

## Mean


## Most frequent

Works well with categorical features.

## k-NN

The [[K-Nearest Neighbors]] algorithm allows to predict new data points based on similarity to existing observations in the data. The new value is the weighted mean of the valus in the closest neighbors.

## Multiple Imputations (MI)

Multiple Imputation works in two steps:
1. Generating replacement values and repeating the procedure many times, ending up with many datasets that have been imputated. 
2. Analyzing the datasets and combining the results.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4638176/#:~:text=Multiple%20imputation%20entails%20two%20stages,sets%20and%20combining%20the%20results
## References

- https://towardsdatascience.com/6-different-ways-to-compensate-for-missing-values-data-imputation-with-examples-6022d9ca0779
- 