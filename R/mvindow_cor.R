# A,B: return series
# win: window size.
"mvwindow_cor" <- function(A, B, win = 30) {
  T <- length(A)
  ctemp <- cor(A, B)
  cor <- rep(ctemp, T)
  if (win < (T + 1)) {
    for (i in win:T) {
      ist <- i - win + 1
      x <- A[ist:i]
      y <- B[ist:i]
      v <- cor(x, y)
      cor[i] <- v
    }
  }
  return(list(correlation = cor))
}

# END