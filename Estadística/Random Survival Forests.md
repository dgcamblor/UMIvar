The main difference with [[Random forests]] is the handling of censored data.

> Comparing RSF and [[Regresión de Cox]] models with traditional scores, it was observed that the AUC of RSF model was higher than that of traditional scores, suggesting that RSF had a better predictive performance and application value in disease prognosis of HS patients. https://bmcmedinformdecismak.biomedcentral.com/articles/10.1186/s12911-023-02293-2#:~:text=Comparing%20RSF%20and%20Cox%20models,disease%20prognosis%20of%20HS%20patients.

## RSF in R

```r
install.packages("randomForestSRC")
library(randomForestSRC)
```

```r
fit_rsf <- rfsrc(formula = Surv(time, status) ~ ., data = data, ntree = 500)
imp_rsf <- vimp(fit_rsf)
```