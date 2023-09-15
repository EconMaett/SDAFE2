"orthogonalize" <- function(M) {
  return(solve(matrix.sqrt(M %*% t(M))) %*% M)
}

# END