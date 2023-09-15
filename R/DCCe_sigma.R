"DCCe_sigma" <- function(theta, innov) {
  
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
  
  lambda <- theta[1]
  # alpha = theta[1]; beta = theta[2];	omega = (1-alpha-beta)
  
  S <- cor(E)
  Q <- S
  Q.temp <- diag(d)
  
  Sigma.t <- array(0, c(d, d, n))
  R.t <- array(0, c(d, d, n))
  Sigma.t[, , 1] <- var(Y)
  R.t[, , 1] <- S
  
  for (t in 2:n) {
    Q <- (1 - lambda) * E[t - 1, ] %*% t(E[t - 1, ]) + lambda * Q
    # Q = S * omega + alpha * E[t-1,] %*% t(E[t-1,]) + beta * Q
    diag(Q.temp) <- 1 / sqrt(diag(Q))
    R.t[, , t] <- (Q.temp) %*% Q %*% (Q.temp)
    Sigma.t[, , t] <- diag(D[t, ]) %*% R.t[, , t] %*% diag(D[t, ])
  }
  
  return(list(Sigma.t = Sigma.t, R.t = R.t))
}

# END