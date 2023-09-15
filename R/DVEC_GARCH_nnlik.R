# **********************************
# DVEC GARCH(1,1) estimation ----
# only for d = 2 components
# uses fGarch and mnormt packages 
# **********************************
"DVEC_GARCH_nllik" <- function(param, innov) {
  omega <- param[1:3]
  alpha <- param[4:6]
  beta  <- param[7:9]
  Y  <- as.matrix(innov)
  d  <- ncol(Y)
  n  <- nrow(Y)
  mu <- numeric(d)
  Sigma.dvec <- array(numeric(n * d * d), c(d, d, n))
  Sigma.dvec[, , 1] <- var(Y)
  V1   <- Sigma.dvec[1, 1, 1]
  V2   <- Sigma.dvec[2, 2, 1]
  V12  <- Sigma.dvec[1, 2, 1]
  llik <- 0
  for (t in 2:n) {
    V1  <- omega[1] + alpha[1] * Y[(t - 1), 1]^2 + beta[1] * V1
    V2  <- omega[2] + alpha[2] * Y[(t - 1), 2]^2 + beta[2] * V2
    V12 <- omega[3] + alpha[3] * Y[(t - 1), 1] * Y[(t - 1), 2] + beta[3] * V12
    Sigma.dvec[1, 1, t] <- V1
    Sigma.dvec[2, 2, t] <- V2
    Sigma.dvec[1, 2, t] <- V12
    Sigma.dvec[2, 1, t] <- V12
    llik <- llik + dmnorm(Y[t, ], mean = mu, varcov = Sigma.dvec[, , t], log = TRUE)
  }
  return(-llik)
}

# END