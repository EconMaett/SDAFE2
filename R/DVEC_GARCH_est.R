"DVEC_GARCH_est" <- function(param, innov) {
  out <- optim(par = param, fn = DVEC_GARCH_nllik, lower = 0.0001, upper = c(rep(Inf, 3), rep(0.9999, 6)), innov = innov, method = "L-BFGS-B", hessian = TRUE, control = list(trace = 6))
  theta.hat <- out$par
  theta.hat.se <- sqrt(diag(solve(out$hessian)))
  return(list(theta.hat = theta.hat, theta.hat.se = theta.hat.se))
}

# END