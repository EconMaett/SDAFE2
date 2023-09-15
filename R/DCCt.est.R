"DCCt.est" <- function(theta, innov, m = 5) {
  out <- optim(
    theta, llik.DCCt, 
    lower = c(0.0001, 0.0001), 
    upper = c(0.9999, 0.999), 
    innov = innov, m = m, 
    method = "L-BFGS-B", 
    hessian = FALSE, 
    control = list(trace = 6)
    )
  # out = optim(theta, llik.DCCt, innov = innov, m = m, hessian = FALSE, control=list(trace=10))
  return(out$par)
}

# END