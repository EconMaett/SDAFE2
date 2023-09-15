"EWMA_nllik" <- function(lambda, innov) {
  clambda <- 1 - lambda
  Sigma.hat <- var(innov)
  invSigma.hat <- chol2inv(chol(Sigma.hat)) 
  # invSigma.hat = solve(Sigma.hat)
  detSigma.hat <- det(Sigma.hat)
  llik <- -0.5 * log(detSigma.hat) - 0.5 * crossprod(innov[1, ], invSigma.hat) %*% innov[1, ]
  llik <- llik - 0.5 * log(detSigma.hat) - 0.5 * crossprod(innov[2, ], invSigma.hat) %*% innov[2, ]
  n <- dim(innov)[1]
  for (i in 3:n) {
    atm1 <- innov[(i - 1), ]
    at <- innov[i, ]
    denom <- 1 - lambda^(i - 1)
    # Sigma.hat = clambda * tcrossprod(atm1) + lambda * Sigma.hat #approx
    Sigma.hat <- (clambda / denom) * tcrossprod(atm1) + (lambda * (1 - lambda^(i - 2)) / denom) * Sigma.hat # exact
    invSigma.hat <- chol2inv(chol(Sigma.hat)) 
    # invSigma.hat = solve(Sigma.hat)
    detSigma.hat <- det(Sigma.hat)
    llik <- llik - 0.5 * (log(detSigma.hat) + crossprod(at, invSigma.hat) %*% at)
  }
  nllik <- -llik
  return(nllik)
}

# END