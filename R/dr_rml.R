#' Fit the DR-RML model
#'
#' Estimate double-robust average causal estimator for restricted mean lifetime
#'
#'
#' @param Y: outcome of interest
#' @param Censor: censoring indicator, 1: observed; 0: right-censored
#' @param A: binary treatment indicator 1: treated; 0: untreated
#' @param X: baseline covariate
#' @param L: restrict timepoint to for RML
#' @param PS: propensity score model
#' @param Reg: outcome regression mo  del
#' @param nboot: number of bootstrap for standard error
#' @return Returns of list with the following components:
#' \itemize{
#'   \item mu1: average effect for treated group
#'   \item mu0: average effect for untreated group
#'   \item ace: average causal effect mu1 - mu0
#'   \item se: standard error estimates for ace with bootstrap
#' }
#' @export
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
#'   Y = pmin(T, C)                        # observed time
#'   delta = ifelse(T <= C, 1, 0)             # censoring indicator
#'   simdata = data.frame(Y = round(Y, 5), delta, A, Z1, Z2, Z3, T)
#'   simdata[order(Y), ]
#' }
#' L = 10       #tmax (setting | L=10; L=20)
#' n = 600      #sample size
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
#' dr_rml(Y = Y,
#'        Censor = Censor,
#'        A = A,
#'        X = X,
#'        L = L,
#'        PS = "logis",
#'        Reg = "lm",
#'        nboot = 10)$ace
#'
#' # Data Application
#' data(gse)
#' Y = dat$Y
#' Censor = dat$Censor
#' A = dat$trt
#' X = dat[, 3:6]
#' L = 365 * 5
#' dr_rml(Y = Y, Censor = Censor, A = A, X = X,
#'        L = L, PS = "logis", Reg = "lm", nboot = 10)
#'
#' library(tmle)
#' prmst = pseudomean(Y, Censor, L)
#' Xt = cova = as.data.frame(model.matrix( ~ -1 + Age + Size + Grade + Er, data = dat))
#' fit = tmle(Y = prmst, A = A, W = cova,
#'   Q.SL.library = c("SL.glm", "SL.glm.interaction", "SL.step"),
#'   g.SL.library =  c("SL.glm", "SL.glm.interaction", "SL.step"))
#' # Other methods
#' tmsl = fit$estimates$ATE$psi
#' library(grf)
#' cffit = causal_forest(Xt, prmst, A)
#' cf = mean(cffit$predictions)
#' }

dr_rml = function(Y,
                  Censor,
                  A,
                  X,
                  L,
                  PS = c("logis", "logis2", "SL", "GBM"),
                  Reg = c("lm", "lm2", "SL"),
                  nboot) {
  res = est = est_func(Y, Censor, A, X, L, PS, Reg)
  if (nboot > 1) res$se = se_func(Y, Censor, A, X, L, PS, Reg, nboot)

  return(res)
}
