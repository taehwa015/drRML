#' @importFrom stats as.formula binomial lm predict sd
NULL
#' Fit the DR-RML model
#'
#' Estimate double-robust average causal estimator for restricted mean lifetime.
#'
#'
#' @param Y: event time.
#' @param Censor: censoring indicator, 1: observed; 0: right-censored.
#' @param A: binary treatment indicator 1: treated; 0: untreated.
#' @param X: baseline covariate.
#' @param L: a scalar value to truncate the timepoint for RML, \code{L} < max(\code{Y}).
#' @param PS: specify the propensity score model one of following options:
#' \itemize{
#'   \item \code{logit}: parametric logistic regression.
#'   \item \code{logit2}: improved parametric logistic regression (Cao et al., 2009).
#'   \item \code{SL}: super learner (package: \pkg{SuperLearner}) with libraries \code{SL.glm, SL.glm.interaction, SL.step}.
#'   \item \code{GBM}: Generalized boosted model (package: \pkg{gbm}).
#' }
#' @param Reg: specify the outcome regression model one of follwoing options:
#' \itemize{
#'   \item \code{lm}: parametric gaussian regression.
#'   \item \code{lm2}: improved parametric gaussian regression (Cao et al., 2009).
#'   \item \code{SL}: super learner (package: \pkg{SuperLearner}) with libraries \code{SL.glm, SL.glm.interaction, SL.step}.
#' }
#' @param nboot: a numeric value to specify the number of bootstrap for standard error of ACE. If standard error estiamtes is unnecessary, use \code{nboot = 0}.
#' 
#' 
#' @return \code{dr_RML} returns a list with the following components:
#' \itemize{
#'   \item mu1: average effect for treated group.
#'   \item mu0: average effect for untreated group.
#'   \item ace: average causal effect mu1 - mu0.
#'   \item se: standard error estimates for ace.
#' }
#' 
#' @details 
#' see Choi et al., (2022+) for detailed method explanation.
#' 
#' @references 
#' Cao, W., Tsiatis, A. A., & Davidian, M. (2009). Improving efficiency and robustness of the doubly robust estimator for a population mean with incomplete data. \emph{Biometrika}, \bold{96}(3), 723--734.
#' 
#' Choi, S., Choi, T., Lee, H. Y., Han, S. W., Bandyopadhyay, D. (2022+). Double-robust methods for estimating differences in restricted mean lifetimes using pseudo-observations. \emph{In revision}.
#' 
#' 
#' @examples
#' \dontrun{
#' simdata = function(n, tau) {
#'   require(mvtnorm)
#'   expit = function(x)
#'     exp(x) / (1 + exp(x))
#'   Sigma = diag(rep(1, 3))
#'   Sigma[1, 3] = Sigma[3, 1] = 0.2
#'   Z = rmvnorm(n,
#'               mean = rep(0, 3),
#'               sigma = Sigma,
#'               method = "chol")
#'   Z1 = Z[, 1]
#'   Z2 = Z[, 2]
#'   Z3 = Z[, 3]
#'   # trt indicator (Bernoulli)
#'   A = rbinom(n, 1, expit(-0.5 * Z1 - 0.5 * Z2))# Z1 and Z2 are involved in trt assignment
#'   par.t = exp((-2.5 - 1.5 * Z1 - 1 * Z2 - 0.7 * Z3) * (A == 0) +
#'                 (-3 - 1 * Z1 - 0.9 * Z2 - 1 * Z3) * (A == 1))
#'   T = rexp(n, par.t)                   # true surv time
#'   C = rexp(n, tau)                     # independent censoring
#'   # tau as a controller of censoring rate
#'   Y = pmin(T, C)                       # observed time
#'   delta = ifelse(T <= C, 1, 0)         # censoring indicator
#'   simdata = data.frame(Y = round(Y, 5), delta, A, Z1, Z2, Z3, T)
#'   simdata[order(Y), ]
#' }
#' L = 10       
#' n = 600      
#' tau = exp(-4.5)
#' truepar = ifelse(L == 10, 0.871, 1.682)
#'
#' set.seed(123)
#' dt = simdata(n, tau)
#' Y = dt$Y
#' Censor = dt$delta
#' A = dt$A
#'
#' X = dt[, 4:6]
#' dr_rml(Y = Y, Censor = Censor, A = A, X = X, L = L, 
#'        PS = "logit", Reg = "lm", nboot = 10)$ace
#'
#' # Data Application
#' data(gse)
#' Y = gse$Y
#' Censor = gse$Censor
#' A = gse$trt
#' X = gse[, 3:6]
#' L = 365 * 5
#' dr_rml(Y = Y, Censor = Censor, A = A, X = X, L = L, 
#'        PS = "logit", Reg = "lm", nboot = 10)
#'
#' library(tmle)
#' library(pseudo)
#' prmst = pseudomean(Y, Censor, L)
#' Xt = cova = as.data.frame(model.matrix( ~ -1 + Age + Size + Grade + Er, data = gse))
#' fit = tmle(
#'   Y = prmst, A = A, W = cova,
#'   Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
#'   g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step")
#' )
#' # Other methods
#' tmsl = fit$estimates$ATE$psi
#' library(grf)
#' cffit = causal_forest(Xt, prmst, A)
#' cf = mean(cffit$predictions)
#' }
#' @export
#' 
dr_rml = function(Y, Censor, A, X, L, PS = c("logit", "logit2", "SL", "GBM"), Reg = c("lm", "lm2", "SL"), nboot) {
  est_func = function(Y, Censor, A, X, L, PS = c("logit", "logit2", "SL", "GBM"), Reg = c("lm", "lm2", "SL")) {
    n = nrow(X)
    p = ncol(X)
    colnames(X) = paste0("X", 1:p)
    prmst = pseudo::pseudomean(time = Y,
                               event = Censor,
                               tmax = L)
    dat = data.frame(prmst, A, X)
    formul = as.formula(paste0("A~", paste0("X", 1:p, collapse = "+")))
    if (PS == "logit") {
      mod.p = stats::glm(formul, family = stats::binomial(), data = dat)
      p.1 = stats::predict.glm(mod.p, type = "response")  # p.1 = Prob(A=1|Z)
    } else if (PS == "logit2") {
      formul = stats::as.formula(paste0("A~", paste0("X", 1:p, collapse = "+")))
      mod.p = stats::glm(formul, family = stats::binomial(), data = dat)
      p.10 = stats::predict.glm(mod.p, type = "response")  # p.1 = Prob(A=1|Z)
      idx1 = which(A == 1)
      idx0 = which(A == 0)
      d = log(mean(A * (1 + exp(p.10)) / exp(p.10)))
      p.1 = as.numeric(exp(d + p.10) / (1 + exp(p.10)))
    } else if (PS == "SL") {
      options(warn = -1)
      Xt = with(dat, model.matrix(as.formula(paste0(
        "~-1+", paste0("X", 1:p, collapse = "+")
      ))))
      Xt = as.data.frame(Xt)
      slpi = SuperLearner::SuperLearner(
        Y = A,
        X = Xt,
        SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
        family = binomial()
      )
      p.1 = slpi$SL.predict
    } else if (PS == "GBM") {
      gfit = gbm::gbm(formul, data = dat, distribution = "bernoulli")
      p.1 = gbm::predict.gbm(gfit, dat, n.trees = 100, type = "response")
    }
    formul = as.formula(paste0("prmst~", paste0("X", 1:p, collapse = "+")))
    if (Reg == "lm") {
      mod.m0 = predict(lm(formul, data = subset(dat, A == 0)), newdata = dat)
      mod.m1 = predict(lm(formul, data = subset(dat, A == 1)), newdata = dat)
    } else if (Reg == "lm2") {
      wt1 = ((1 - p.1) / p.1 ^ 2)[idx1]
      wt0 = ((p.1) / (1 - p.1) ^ 2)[idx0]
      formul = as.formula(paste0("prmst~", paste0("X", 1:p, collapse = "+")))
      mod.m1 = predict(lm(formul, data = subset(dat, A == 1), weights = wt1), dat)
      mod.m0 = predict(lm(formul, data = subset(dat, A == 0), weights = wt0), dat)
    } else if (Reg == "SL") {
      mod.m0 = with(dat,
                    SuperLearner::SuperLearner(
                      Y = prmst[A == 1],
                      X = Xt[A == 1, ],
                      newX = Xt,
                      SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
                      family = gaussian()
                    ))$SL.predict
      mod.m1 = with(dat,
                    SuperLearner::SuperLearner(
                      Y = prmst[A == 0],
                      X = Xt[A == 0, ],
                      newX = Xt,
                      SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
                      family = gaussian()
                    ))$SL.predict
    }
    mu1 = (prmst * A) / p.1 - ((A - p.1) * mod.m1) / p.1
    mu0 = (prmst * (1 - A)) / (1 - p.1) -
      (((1 - A) - (1 - p.1)) * mod.m0) / (1 - p.1)
    res = list(mu1 = mean(mu1),
               mu0 = mean(mu0),
               ace = mean(mu1 - mu0))
    
    return(res)
  }
  
  res = est = est_func(Y, Censor, A, X, L, PS, Reg)
  if (nboot > 1) {
    res$se = sd(replicate(nboot, {
      n = length(Y)
      ind = sample(n, n, replace = TRUE)
      est_func(Y[ind], Censor[ind], A[ind], X[ind, ], L, PS, Reg)$ace
    }))
  }
  return(res)
}
