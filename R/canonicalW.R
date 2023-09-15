"canonicalW" <- function(W) {
  
  if (dim(W)[1] != dim(W)[2]) {
    stop("W must be a square matrix")
  }
  
  d <- dim(W)[1]
  
  if(sum(abs(t(W)%*%W  - diag(d))) > 1e-10){ stop("W must be an orthogonal matrix") }
  
  W.temp <- W
  W.new <- matrix(0, d, d)
  
  for (i in 1:d) {
    index <- which.max(abs(W))
    row.index <- index %% d # row
    row.index <- ifelse(row.index == 0, d, row.index)
    col.index <- ceiling(index / d) # col
    w.i <- W.temp[row.index, ]
    W.new[col.index, ] <- w.i * ifelse(w.i[col.index] < 0, -1, 1)
    W[row.index, ] <- 0
    W[, col.index] <- 0
  }
  
  if (det(W.new) < 0) {
    W.new[d, ] <- -W.new[d, ]
  }
  
  return(W.new)
}

# END