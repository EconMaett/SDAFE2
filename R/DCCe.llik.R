"llik.DCCe" <- function(theta, innov) {
  # llik for the correlation component
  E <- innov
  n <- dim(E)[1]
  d <- dim(E)[2]
  
  lambda <- theta[1]
  # alpha = theta[1]; beta = theta[2];	omega = (1-alpha-beta)
  
  S <- cor(E)
  Q <- S
  Q.temp <- diag(d)
  
  llik <- 0
  for (t in 2:n) {
    Q <- (1 - lambda) * E[t - 1, ] %*% t(E[t - 1, ]) + lambda * Q
    # Q = S * omega + alpha * E[t-1,] %*% t(E[t-1,]) + beta * Q
    diag(Q.temp) <- 1 / sqrt(diag(Q))
    R <- (Q.temp) %*% Q %*% (Q.temp)
    llik <- llik + log(det(R)) + t(E[t, ]) %*% solve(R) %*% E[t, ] - t(E[t, ]) %*% E[t, ]
  }
  return(llik / 2)
}

# END