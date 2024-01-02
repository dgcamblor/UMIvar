
## Accuracy

Accuracy measures the percentage of correct predictions made by the model

- Can be misleading if the dataset is imbalanced, meaning one class has significantly more samples than the other
- May not be the best metric to use if the cost of false positives and false negatives is different

## ROC-AUC

ROC-AUC measures the ability of the model to distinguish between positive and negative classes. Takes into account the trade-off between true positive rate (TPR) and false positive rate (FPR).

- Useful when the dataset is imbalanced
- Can be more informative than accuracy when the cost of false positives and false negatives is different