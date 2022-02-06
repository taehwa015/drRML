se_func = function(Y,
                   Censor,
                   A,
                   X,
                   L,
                   PS = c("logis", "logis2", "SL", "GBM"),
                   Reg = c("lm", "lm2", "SL"),
                   nboot) {
  se = sd(replicate(nboot, {
    n = length(Y)
    ind = sample(n, n, replace = TRUE)
    est_func(Y[ind], Censor[ind], A[ind], X[ind, ], L, PS, Reg)$ace
  }))

  return(se)
}
