# For a given angle theta, returns a d x d Givens rotation matrix
# Ex: for i < j , d = 2:  (c -s)
#                         (s  c)
"givens_rotation" <- function(theta = 0, d = 2, which = c(1, 2)) {
  c <- cos(theta)
  s <- sin(theta)
  M <- diag(d)
  a <- which[1]
  b <- which[2]
  M[a, a] <- c
  M[b, b] <- c
  M[a, b] <- -s
  M[b, a] <- s
  return(M)
}

# END