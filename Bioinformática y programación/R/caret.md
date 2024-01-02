---
tags:
  - package
---

The `caret` package serves as a unified interface to multiple machine learning algorithms. 
## Data partition

By default, `createDataPartition()` creates a stratified random split of the data.

```r
createDataPartition()
```

## Preprocessing

### Dummy variables (one-hot encoding)

```r
dummies <- dummyVars(survived ~ ., data = etitanic)
```

## Model selection and training

![[TrainAlgo.png]]

### Resampling for parameter tuning

Options for resampling are:

- *k*-fold cross-validation (once or repeated)
- Leave-one-out cross-validation
- Bootstrap (performed by default using the `trainControl()` function in `train()` internally).

```r
ctrl <- trainControl(
	CV method = "repeatedcv", 
	savePredictions = TRUE,  # Store predicted values in train object
	number = 10,  # 10-fold
	repeats = 10  # Repeated 10 times
	)
```

### Tuning grid

By default, caret will estimate a [[Tuning grid]] for each method. However, sometimes the defaults are not the most sensible given the nature of the data. The `tuneGrid` argument allows the user to specify a custom grid of tuning parameters as opposed to simply using what exists implicitly.

### Training the model

Training the model of interest is performed using the `train()` function.

```r
xgb_fit <- train(x = X_train,
				 y = train$Survived,
				 method = "xgbTree",
				 trControl = ctrl,
				 metric = "Accuracy",
				 verbose = "False")
```

The `train()` function can perform data preprocessing prior to model fitting. The `preProcess()` function is used for that purpose.

### Comparing models

The easiest way to compare trained models is using the `MLeval` package ([CRAN - Package MLeval](https://cran.r-project.org/web/packages/MLeval/index.html)).

## References

- [The caret Package](https://topepo.github.io/caret/)
- [Machine Learning con R y caret](https://cienciadedatos.net/documentos/41_machine_learning_con_r_y_caret#Introducci%C3%B3n)

