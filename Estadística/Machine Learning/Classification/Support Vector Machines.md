---
machineLearningProblem:
  - Classification
---
Maximal Margin Classifier -> Sensitive to outliers. Allow misclasifications -> Soft margin. Soft Margin Classifier = Support Vector Classifier.

Maximal Margin Classifiers and Support Vector Classifiers are not good when the separation needs to be doubled in a line, etc. For other dimensions.

SVM add an extra dimension.
1.       Start with data in a relatively low dimension.
2.       Move the data to a higher dimension
3.       Find a support vector classifier that separates the higher dimensional data into two groups.

SVMs use Kernel functions. The polynomial kernel adds dimensions. We can find the d parameter with cross validation.

Kernel Trick -> Data is not actually transformed to a higher dimension.

R and d are computed using cross-validation

## Kernel SVM

