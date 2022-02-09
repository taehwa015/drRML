# `drRML`: R package for estimating differences in restricted mean lifetimes

## Overview

`drRML` is the R package to estimate differences in restricted mean lifetimes (RML) of two treatments using pseudo-observations.
Define $T$ and $C$ as the event time and right-censoring times, respectively.
Then we observe $\{Y_i = T_i\wedge C_i,\Delta_i = I(T_i\le C_i),Z_i\}_{i=1}^n$, 
where $Z_i$ is $p-$dimensional baseline covariate vector.
We adopt the pseudo-observation (available in `pseudo` R package)
to replace unknown event times,
$$
\hat \theta^*_{i}
	= n\cdot\int_{0}^{L}\hat{S}(t) dt - (n-1)\cdot\int_{0}^{L}\hat{S}^{-i}(t)dt,
$$
where $L$ is fixed time point for RML and $\hat S(t)$ is the empirical survival estimates, such as the Kaplan-Meier estimator.
Then, the double robust estimator for average causal effect is defined as
$$
\hat\delta_{DR} = \hat\mu_1^{DR}-\hat\mu_0^{DR},
$$
where $\hat\mu_a^{DR}~(a=0,1)$ can be expressed as
$$
\hat\mu_a^{DR} = n^{-1}\sum_{i=1}^n
	\frac{I(A_i=a)}{\tilde \pi_a(Z_i;\hat\gamma)}\hat\theta^*_i
	+ n^{-1}\sum_{i=1}^n \left\{1- \frac{I(A_i=a)}{\tilde\pi_a(Z_i;\hat\gamma)} \right\} m(a,Z_i;\hat\beta).
$$
Here, $\tilde \pi_a(Z_i;\gamma) = a\pi(Z;\gamma) + (1-a) \{1-\pi(Z;\gamma) \}$ 
is the propensity score, estimated by logistic regression or machine learning approaches (e.g., `SuperLearner`, `gbm` packages in R).
In addition, $m(a,Z_i;\beta)$ is the outcome regression model, estimated by
linear regression or machine learning approaches (e.g., `SuperLearner` in R).
See Choi et al. (2022+) for more detailed description of the method.



## Installation
```r
library(devtools)
devtools::install_github(repo='taehwa015/drRML')
```

## Status of development

The code provided here is only being provided for research purposes.

## Usage

The `drRML` package provides a double-robust estimator for differences in restricted mean lifetimes using pseudo-observations.

Below example is the GSE6532 data application in the article.
```r
library(drRML)
data(gse)
 Y = gse$Y
 Censor = gse$Censor
 A = gse$trt
 X = gse[, 3:6]
 L = 365 * 5
 dr_rml(Y = Y, Censor = Censor, A = A, X = X, L = L, 
        PS = "logit", Reg = "lm", nboot = 10)

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



