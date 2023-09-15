"DCCt_fit" <- function(theta.0 = 0.9, innov, m = 5) {
  
  Y <- innov
  n <- dim(Y)[1]
  d <- dim(Y)[2]
  
  D <- matrix(0, n, d)
  E <- matrix(0, n, d)
  
  garch.omega.alpha.beta <- matrix(0, d, 3)
  
  for (i in 1:d) {
    fit.temp <- garchFit(data = Y[, i], include.mean = FALSE, trace = FALSE)
    garch.omega.alpha.beta[i, ] <- fit.temp@fit$matcoef[, 1]
    D[, i] <- sqrt(fit.temp@h.t)
    E[, i] <- Y[, i] / sqrt(fit.temp@h.t)
  }
  
  DCCt.params <- DCCt_est(theta = theta.0, innov = E, m = m)
  print(DCCt.params)
  DCCt.Sigma  <- DCCt_sigma(theta = DCCt.params, innov = Y, m = m)
  
  return(list(Sigma.t = DCCt.Sigma$Sigma.t, R.t = DCCt.Sigma$R.t, coef = garch.omega.alpha.beta, lambda = DCCt.params))
}

# END