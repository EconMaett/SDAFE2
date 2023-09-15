"EWMA_sigma" <- function(lambda, innov) {
  clambda <- 1 - lambda
  Sigma.hat <- var(innov)
  n <- dim(innov)[1]
  d <- dim(innov)[2]
  Sigma.t <- array(0, c(d, d, n))
  Sigma.t[, , 1:2] <- Sigma.hat
  for (i in 3:n) {
    atm1 <- innov[(i - 1), ]
    at <- innov[i, ]
    denom <- 1 - lambda^(i - 1)
    # Sigma.hat = clambda * tcrossprod(atm1) + lambda * Sigma.hat # approx
    Sigma.t[, , i] <- (clambda / denom) * tcrossprod(atm1) + (lambda * (1 - lambda^(i - 2)) / denom) * Sigma.t[, , (i - 1)] # exact
  }
  return(Sigma.t)
}

# END