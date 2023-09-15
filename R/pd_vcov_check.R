pd.vcov.check <- function(Sigma) {
  dex <- NULL
  if (is.na(dim(Sigma)[3])) {
    sva <- svd(Sigma)
    temp <- min(sva$d)
    if (min(sva$d) <= 0) {
      dex <- c(dex, t)
    }
  } else {
    N <- dim(Sigma)[3]
    temp <- numeric(N)
    for (t in 1:N) {
      sva <- svd(Sigma[, , t])
      temp[t] <- min(sva$d)
      if (min(sva$d) <= 0) {
        dex <- c(dex, t)
      }
    }
  }
  return(list(index = dex, min.sv = temp))
}

# END