"orthogonalize" <- function(M) {
  return(solve(matrix_sqrt(M %*% t(M))) %*% M)
}

# END