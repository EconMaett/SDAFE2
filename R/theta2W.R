# For a vector of angles theta, returns W, a d x d Givens rotation matrix:
# W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2
# if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
"theta2W" <- function(theta) {
  d <- (sqrt(8 * length(theta) + 1) + 1) / 2
  if (d - floor(d) != 0) {
    stop("theta must have length: d(d-1)/2")
  }
  W <- diag(d)
  index <- 1
  for (j in 2:d) {
    for (i in (j - 1):1) {
      Q.ij <- givens.rotation(theta[index], d, c(i, j))
      W <- Q.ij %*% W
      index <- index + 1
    }
  }
  return(W)
}

# END