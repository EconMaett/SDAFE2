"DCCt_sigma" <- function(theta, innov, m = 5) {
  Y <- innov
  n <- dim(Y)[1]
  d <- dim(Y)[2]
  
  D <- matrix(0, n, d)
  E <- matrix(0, n, d)
  
  for (i in 1:d) {
    fit.temp <- garchFit(data = Y[, i], include.mean = FALSE, trace = FALSE)
    D[, i] <- sqrt(fit.temp@h.t)
    E[, i] <- Y[, i] / sqrt(fit.temp@h.t)
  }
  
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta12 <- (1 - theta1 - theta2)
  
  Gamma <- cor(D)
  Gamma.t <- Gamma
  
  Sigma.t <- array(0, c(d, d, n))
  R.t <- array(0, c(d, d, n))
  Sigma.t[, , 1:m] <- var(Y)
  R.t[, , 1:m] <- Gamma
  
  for (t in (m + 1):n) {
    Gamma.t <- theta12 * Gamma + theta1 * Gamma.t + theta2 * cor(D[(t - 1):(t - m), ])
    R.t[, , t] <- Gamma.t
    Sigma.t[, , t] <- diag(D[t, ]) %*% R.t[, , t] %*% diag(D[t, ])
  }
  
  return(list(Sigma.t = Sigma.t, R.t = R.t))
}

# END