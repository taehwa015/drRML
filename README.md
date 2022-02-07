# Double-robust methods for estimating differences in restricted mean lifetimes using pseudo-observations

## Installation
```r
library(devtools)
devtools::install_github(repo='taehwa015/drRML')
```

## Status of developement

The code provided here is only being provided for research purposes.

## Usage

The `drRML` package provides a double-robust estimator for differences in restricted mean lifetimes using pseudo-observations.
See Choi et al. (2022+) for Detailed description of the method.

Below example is the GSE6532 data application in the article.
```r
data(gse)
 Y = dat$Y
 Censor = dat$Censor
 A = dat$trt
 X = dat[, 3:6]
 L = 365 * 5
 dr_rml(Y = Y, Censor = Censor, A = A, X = X, L = L, 
        PS = "logit", Reg = "lm", nboot = 10)

 library(tmle)
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



