Cross-validation is a technique that allows to evaluate a ML model by splitting a [[Training dataset]] into different [[Training dataset|Training datasets]] and [[Validation dataset|Validation datasets]]. The model is trained on each of the [[Training dataset|Training datasets]] and evaluated on the corresponding [[Validation dataset|Validation datasets]]. The final performance of the model is the average of the performance on each of the [[Training dataset|Training datasets]] and [[Validation dataset|Validation datasets]].

Note that cross-validation is performed on the [[Training dataset]]. As a general rule, during training and model tuning, the model should not be exposed to the [[Testing dataset]] until the very end, where the model is finally evaluated.

There are several types of cross-validation:

- K-fold cross-validation
- Leave-one-out cross-validation

## K-fold cross-validation

In K-fold cross-validation, the [[Training dataset]] is split into K different [[Training dataset|Training datasets]] and [[Validation dataset|Validation datasets]].

Sampling can be stratified, meaning that the proportion of samples of each class is preserved in each [[Training dataset]] and [[Validation dataset]].

## Leave-one-out cross-validation

In leave-one-out cross-validation, the [[Training dataset]] is split into N different [[Training dataset|Training datasets]] and [[Validation dataset|Validation datasets]], where N is the number of samples in the [[Training dataset]]. Each [[Training dataset|Training dataset]] contains all but one sample of the [[Training dataset]], and the corresponding [[Validation dataset|Validation dataset]] contains the remaining sample.