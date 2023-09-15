# multivariate Ljung--Box tests ----
# Compute multivariate Ljung-Box test statistics
"mLjungBox" <- function(x, lag = 1, df.adj = 0) {
  x <- as.matrix(x)
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  g0 <- var(x)
  ginv <- solve(g0)
  qm <- 0.0
  df <- 0
  out <- matrix(0, lag, 4)
  for (i in 1:lag) {
    x1 <- x[(i + 1):nr, ]
    x2 <- x[1:(nr - i), ]
    g <- cov(x1, x2)
    g <- g * (nr - i - 1) / (nr - 1)
    h <- t(g) %*% ginv %*% g %*% ginv
    qm <- qm + nr * nr * sum(diag(h)) / (nr - i)
    df <- df + nc * nc
    pv <- 1 - pchisq(qm, df - df.adj)
    print(c(i, qm, pv))
    out[i, ] <- c(i, round(qm, 2), round(df - df.adj), round(pv, 3))
  }
  output <- as.data.frame(matrix(out[lag, ], 1, 4))
  names(output) <- c("K", "Q(K)", "d.f.", "p-value")
  return(output)
}

# END