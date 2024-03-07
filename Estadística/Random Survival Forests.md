The main difference with [[Random forests]] is the handling of censored data.
## RSF in R

```r
install.packages("randomForestSRC")
library(randomForestSRC)
```

```r
fit_rsf <- rfsrc(formula = Surv(time, status) ~ ., data = data, ntree = 500)
imp_rsf <- vimp(fit_rsf)
```