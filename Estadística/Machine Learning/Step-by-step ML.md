---
tags:
  - how-to
---

Prior to initiating the ML process, you must be sure ML is the correct application.

## Defining the objective and getting the data

Clearly defining the problem at hand. Getting data relevant to the problem.

## Exploratory data analysis (EDA)

The **whole data** must be explored in order to identify key inputs that might be useful for the ML process. Some key elements are:

- Summary statistics
- Data visualization: histograms, boxplots, etc.
- Correlation analysis. Examining the correlations between features in the data. Using [[Correlation]] or [[Mutual information]].

Along with the EDA some processing of the dataset can be done: codification of the variables as factors, etc.

## Splitting the data

Splitting the data ([[Train-Test Split]]) into two datasets: [[Training dataset]] and [[Testing dataset]]. The [[Training dataset]] is often split into the real [[Training dataset]] and a [[Validation dataset]].

More on the three types of datasets: [About Train, Validation and Test Sets in Machine Learning | by Tarang Shah | Towards Data Science](https://towardsdatascience.com/train-validation-and-test-sets-72cb40cba9e7)

## Data preprocessing

The data preprocessing step is focused on modifying the dataset so it can serve as input for the different machine learning algorithms. These steps should be done in the [[Training dataset]] separately, so that the [[Testing dataset]] is not biased.

- **Data cleaning.**
	- Handling missing values. **Are the values missing because it weren't recorded or because they don't exist?** If they weren't recorded, one of several [[Data imputation methods]] can be used.
	- [[One-Hot encoding]] the categorical variables.
	- [[Scaling]] and/or [[Normalizing]]. Important in distance-based algorithms and algorithms that assume normality.
	- Dealing with data inconsistencies and outliers.
- **[[Feature engineering]].** Creating features to improve the model.

## Training the models

A selection of [[Machine learning algorithms]] is trained using the [[Training dataset]]. This involves choosing the most suitable algorithms and [[Hyperparameters]]. 

- The algorithms are chosen based on the type of problem (classification, regression, clustering, etc.) and the type of data (numerical, categorical, etc.).
- The [[Hyperparameters]] are chosen based on the type of algorithm and the type of data.

Models are implemented in Python using [[Scikit-learn]], or in R using [[caret]].

## Validation and tuning

The models are validated using a [[Validation dataset]]. The best model is selected and its [[Hyperparameters]] are [[Hyperparameter tuning|tuned]]. The model is retrained using the [[Training dataset]] and the [[Validation dataset]].

A technique to create splits and avoid overfitting is [[Cross-validation]]. This technique allows to tune the [[Hyperparameters]].

## Final evaluation

The [[Testing dataset]] has been kept separate to test the final model. This step provides unbiased assessment of the model's performance.

The [[Testing dataset]] must be used sparingly and only at the end of the process. If the model is tuned using the [[Testing dataset]], it will not be able to provide an unbiased assessment of the model's performance.

## Production steps

If the model has passed the statistical testing, you must make it ready for production, including:

- **Live experiment.** Deploy the model in a live experiment to test it in real-world conditions.
- **Monitor and mantain.** Monitor the model's performance and retrain it if necessary.

## References

- [Step-by-step guide to AI projects - YouTube](https://www.youtube.com/watch?v=2caALBeiMAo)