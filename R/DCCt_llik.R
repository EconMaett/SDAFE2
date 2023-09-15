"DCCt.llik" <- function(theta, innov, m = 5) {
  # llik for the correlation component
  E <- innov
  n <- dim(E)[1]
  d <- dim(E)[2]
  
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta12 <- (1 - theta1 - theta2)
  
  Gamma <- cor(E)
  Gamma.t <- Gamma
  
  llik <- 0
  for (t in (m + 1):n) {
    Gamma.t <- theta12 * Gamma + theta1 * Gamma.t + theta2 * cor(E[(t - 1):(t - m), ])
    llik <- llik + log(det(Gamma.t)) + t(E[t, ]) %*% solve(Gamma.t) %*% E[t, ] - t(E[t, ]) %*% E[t, ]
  }
  return(llik / 2)
}

# END