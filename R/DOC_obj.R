# Compute the unconstrained objective function for the DOC in volatility algorithm.
# For use in estimating mixing parameter, theta (px1).
# Time series vector Z (nxd) should be standardized.
# Truncation level for computing Huber's function c (1x1) (non-negative).
# Max_lag L (1x1) for lags 0,1,...,L.
# Diagonal weighting matrix.
# See Matteson and Tsay (2010)
"DOC_obj" <- function(theta, Z, c, L) {
  p <- length(theta)
  p2 <- 2 * p
  q <- L * p2 + p
  N <- L + 1
  sizeZ <- dim(Z)
  n <- sizeZ[1]
  d <- sizeZ[2]
  
  # Compute W from theta
  W <- theta2W(theta)
  
  # Separate signals with estimate of W: S (nxd)
  S <- Z %*% t(W)
  
  # Transforming signals with Hubers function, saving column means
  SH <- myHuber(S, c)
  SH.bar <- colMeans(SH)
  
  # Computing objective function
  f.bar <- numeric(q)
  iter <- 0:(2 * L) * p + 1
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      SH.barXSH.bar <- SH.bar[i] * SH.bar[j]
      f.bar[iter[1]] <- crossprod(SH[, i], SH[, j]) / n - SH.barXSH.bar
      for (ell in 1:L) {
        f.bar[iter[2 * ell]] <- crossprod(SH[-(1:ell), i], SH[-((n - ell + 1):n), j]) / n - SH.barXSH.bar
        f.bar[iter[2 * ell + 1]] <- crossprod(SH[-(1:ell), j], SH[-((n - ell + 1):n), i]) / n - SH.barXSH.bar
      }
      iter <- iter + 1
    }
  }
  phi <- 1 - (0:L) / N
  phi.total <- sum(phi)
  phi <- phi / phi.total
  # PHI = c( rep(phi[1],p), rep(phi[-1],each = p2) )
  PHI <- c(rep(phi[1] / p, p), rep(phi[-1] / p2, each = p2))
  PHI <- diag(PHI)
  return(crossprod(f.bar, PHI) %*% f.bar)
}

# END