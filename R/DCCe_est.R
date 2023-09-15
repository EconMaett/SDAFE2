"DCCe_est" <- function(theta, innov) {
  out <- optim(theta, llik.DCCe, lower = 0.001, upper = 0.999, innov = innov, method = "L-BFGS-B", hessian = FALSE)
  return(out$par)
}

# END