---
machineLearningProblem:
  - Regression
---

Lasso is a regression analysis that performs both variable selection and [[Regularization]].

Minimizes: the sum of the squared residuals + labda x |slope|

Lasso regression can make coefficients go to 0 as we increase the lambda parameter. It can exclude useless variables from equations. [[Ridge regression]] is better when most variables are useful, Lasso is better to exclude useless variables.