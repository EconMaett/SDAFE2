# rt: return series
# win: window size.
"mvwindow_var" <- function(rt, win = 30) {
  T <- length(rt)
  vtemp <- var(rt)
  vol <- rep(vtemp, T)
  if (win < (T + 1)) {
    for (i in win:T) {
      ist <- i - win + 1
      x <- rt[ist:i]
      v <- var(x)
      vol[i] <- v
    }
  }
  return(list(volatility = vol))
}

# END