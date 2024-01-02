Feature engineering involves the transformation of the variables in a dataset to produce more meaningful new variables for a ML model to be trained on. This is specially important when considering that models have difficulty learning complex relationships between variables.

Feature engineering requires a certain amount of knowledge about the data at hand. #insight

The metric of [[Mutual information]] can be useful in establishing relationships between variables: by looking into features with low MI, you can develop them further to increase their relation to the target variable. A [[PCA (Principal Component Analysis)]] can also be useful, finding insights among the loadings.

Some common transformations are:

- Mathematical transformations. Normalization is useful for [[Simple linear regression]].
- Counts (adding up categorical variables). Helpful for tree models.
- Grouped transforms. Adding information to rows based on group information (e.g., mean of a variable in a group). The information should be that of the [[Training dataset]].
- Adding groups formed with [[K-means clustering]] or [[PCA (Principal Component Analysis)]].
- [[Target encoding]]. 