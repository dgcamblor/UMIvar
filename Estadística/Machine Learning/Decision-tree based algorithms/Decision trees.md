---
machineLearningProblem:
  - Classification
  - Regression
needsDummies: false
needsImputation: false
---

## Complexity parameter (cp)

The **complexity parameter (cp)** in rpart is a stopping parameter that controls the size of the decision tree and helps to select the optimal tree size. It is the minimum improvement in the model needed at each node, based on the cost complexity of the model defined as the sum of the misclassification rate and a complexity penalty. The complexity penalty is proportional to the number of terminal nodes in the tree, so increasing the value of cp results in a smaller tree with fewer terminal nodes.