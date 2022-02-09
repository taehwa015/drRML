# `drRML`: R package for estimating differences in restricted mean lifetimes

## Overview

`drRML` is the R package to estimate differences in restricted mean lifetimes (RML) of two treatments using pseudo-observations.

## Installation
```r
library(devtools)
devtools::install_github(repo='taehwa015/drRML')
```

## Status of development

The code provided here is only being provided for research purposes.

## Documentation

Vignette is available at [here](http://htmlpreview.github.io/?https://github.com/taehwa015/drRML/blob/master/vignettes/drRML.html)

## Usage

The `drRML` package provides a double-robust estimator for differences in restricted mean lifetimes using pseudo-observations.
See Choi et al. (2022+) for more detailed description of the method.

Below example is the GSE6532 data application in the article.
```r
library(drRML)
data(gse)
Y = gse$Y
Censor = gse$Censor
A = gse$trt
X = gse[, 3:6]
L = 365 * 5
dr_rml(Y = Y, Censor = Censor, A = A, Xps = X, Xreg = X,
        L = L, PS = "logit", Reg = "lm", nboot = 10)

library(tmle)
library(pseudo)
prmst = pseudomean(Y, Censor, L)
Xt = cova = as.data.frame(model.matrix( ~ -1 + Age + Size + Grade + Er, data = dat))
fit = tmle(
  Y = prmst, A = A, W = cova,
  Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
  g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step")
)
# Other methods
tmsl = fit$estimates$ATE$psi
library(grf)
cffit = causal_forest(Xt, prmst, A)
cf = mean(cffit$predictions)
}
```



