# BayesNormTest
R package for Bayesian normality test using local priors and non-local priors

To install the package use

```
library(devtools)
install_github("JingPan1003/BayesNormTest")
library(BayesNormTest)

# In order to use the nested functionality, please install the following packages
install_github("FJRubio67/twopiece")
library(twopiece)
library(numDeriv)
```

To apply the Bayesian normality test using local and non-local priors, please see

```
?BayesNormTestLocal
?BayesNormTestNonLocal
```
