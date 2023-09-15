"DCCe.fit" <- function(theta.0 = 0.9, innov) {
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
  
  DCCe.params <- est.DCCe(theta = theta.0, innov = E)
  DCCe.params
  DCCe.Sigma <- sigma.DCCe(theta = DCCe.params, innov = Y)
  
  return(list(Sigma.t = DCCe.Sigma$Sigma.t, R.t = DCCe.Sigma$R.t, coef = garch.omega.alpha.beta, lambda = DCCe.params))
}

# END