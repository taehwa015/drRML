est_func = function(Y,
                    Censor,
                    A,
                    X,
                    L,
                    PS = c("logis", "logis2", "SL", "GBM"),
                    Reg = c("lm", "lm2", "SL")) {
  n = nrow(X)
  p = ncol(X)
  colnames(X) = paste0("X", 1:p)
  prmst = pseudo::pseudomean(time = Y,
                             event = Censor,
                             tmax = L)
  dat = data.frame(prmst, A, X)
  formul = as.formula(paste0("A~", paste0("X", 1:p, collapse = "+")))
  if (PS == "logis") {
    mod.p = glm(formul, family = binomial, data = dat)
    p.1 = predict(mod.p, type = "response")  # p.1 = Prob(A=1|Z)
  } else if (PS == "logis2") {
    formul = as.formula(paste0("A~", paste0("X", 1:p, collapse = "+")))
    mod.p = glm(formul, family = binomial, data = dat)
    p.10 = predict(mod.p, type = "response")  # p.1 = Prob(A=1|Z)
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
