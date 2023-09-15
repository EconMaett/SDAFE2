"DOC_test" <- function(A, m) {
  N <- dim(A)[1]
  k <- dim(A)[2]
  temp <- numeric(m)
  pval <- numeric(m)
  out <- as.data.frame(matrix(0, m + 1, 4))
  names(out) <- c("m", "Q(m)", "d.f.", "p-value")
  out[, 1] <- 0:m
  
  Q.temp <- N * sum(cor(A)[lower.tri(cor(A), diag = FALSE)]^2)
  out[1, 2] <- Q.temp
  df.temp <- k * (k - 1) / 2
  out[1, 3] <- df.temp
  out[1, 4] <- 1 - pchisq(Q.temp, df.temp)
  
  for (j in 1:m) {
    ccf <- cor(A[-(1:j), ], A[-((N - j + 1):N), ])
    Q.temp <- Q.temp + N * (N + 2) * sum(ccf[lower.tri(ccf, diag = FALSE)]^2) / (N - j)
    Q.temp <- Q.temp + N * (N + 2) * sum(ccf[upper.tri(ccf, diag = FALSE)]^2) / (N - j)
    out[(j + 1), 2] <- Q.temp
    df.temp <- df.temp + k * (k - 1)
    out[(j + 1), 3] <- df.temp
    out[(j + 1), 4] <- 1 - pchisq(Q.temp, df.temp)
  }
  return(round(out, 3))
}

# END