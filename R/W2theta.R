# Decompose a d by d orthogonal matrix W into the product of
# d(d-1)/2 Givens rotation matrices. Returns theta, the d(d-1)/2 by 1
# vector of angles, theta.
# W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2
"W2theta" <- function(W) {
  if (dim(W)[1] != dim(W)[2]) {
    stop("W must be a square matrix")
  }
  d <- dim(W)[1]
  if(sum(abs(t(W)%*%W  - diag(d))) > 1e-10){stop("W must be an orthogonal matrix")}
  theta <- numeric(d * (d - 1) / 2)
  index <- 1
  for (j in d:2) {
    for (i in 1:(j - 1)) {
      x <- W[j, j]
      y <- W[i, j]
      theta.temp <- atan2(y, x)
      Q.ij <- givens_rotation(theta.temp, d, c(i, j))
      W.temp <- Q.ij %*% W
      W <- W.temp
      theta[index] <- theta.temp
      index <- index + 1
    }
  }
  return(-1 * rev(theta))
}

# END