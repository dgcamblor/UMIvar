---
field: statistics
---

The confusion matrix is a two-dimensional contigency table that contains the number of predicted and actual labels, allowing to visualize the performance of a test or algorithm.

| | Predicted: 0 | Predicted: 1 |
| --- | --- | --- |
| Actual: 0 | True Negative (TN) | False Positive (FP) |
| Actual: 1 | False Negative (FN) | True Positive (TP) |

## Metrics

A description of the direct metrics is:

- True Positive (TP): The number of positive cases correctly classified as positive.
- True Negative (TN): The number of negative cases correctly classified as negative.
- False Positive (FP) ([[Error tipo I (FP)]]): The number of negative cases incorrectly classified as positive.
- False Negative (FN) ([[Error tipo II (FN)]]): The number of positive cases incorrectly classified as negative.

The metrics derived from the confusion matrix are:

- Sensitivity or recall (sensiblidad): The proportion of positive cases correctly classified as positive. It is the probability that a positive case is correctly classified as positive. It is calculated as $TP/(TP+FN)$.
- Specificity (especificidad): The proportion of negative cases correctly classified as negative. It is the probability that a negative case is correctly classified as negative. It is calculated as $TN/(TN+FP)$.
- Precision (precisión): The proportion of positive cases correctly classified as positive. It is the probability that a positive case is correctly classified as positive. It is calculated as $TP/(TP+FP)$.
## Other metrics

[[ROC Curve]]