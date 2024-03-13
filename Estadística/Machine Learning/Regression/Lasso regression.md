---
machineLearningProblem:
  - Regression
---

Lasso is a regression analysis that performs both variable selection and [[Regularization]].

Minimizes: the sum of the squared residuals + labda x |slope|

Lasso regression can make coefficients go to 0 as we increase the lambda parameter. It can exclude useless variables from equations. [[Ridge regression]] is better when most variables are useful, Lasso is better to exclude useless variables.

## Lasso regression in R

Lasso regression in R can be performed using the `glmnet` package. Documentation can be found at: [CRAN - Package glmnet](https://cran.r-project.org/web/packages/glmnet/).

It admits fitting a Generalized Linear Model (GLM) with regularization via [[Lasso regression]] or [[Elastic net]]. Admits models such as [[Regresión de Cox]].

```
# Use alpha = 1 for LASSO
cl_fit <- glmnet(x = predictor_matrix, y = surv_obj, family = "cox", alpha = 1)
```