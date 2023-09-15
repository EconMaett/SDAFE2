# Compute Huber's function with truncation level 'c'.
"myHuber" <- function(S, c) {
  if (c < 0) {
    stop("truncation level 'c' must be non-negative")
  }
  S <- abs(S)
  SH <- S^2 * ifelse(S <= c, 1, 0) + (2 * c * S - c^2) * ifelse(S > c, 1, 0)
  return(SH)
}

# END
