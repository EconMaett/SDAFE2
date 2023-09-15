"ewma.est" <- function(lambda.0, innov) {
  out <- optim(
    par = lambda.0, 
    fn = nllik.ewma, 
    lower = 0.001, 
    upper = 0.999, 
    innov = innov, 
    method = "L-BFGS-B", 
    hessian = TRUE
    )
  lambda.hat <- out$par
  lambda.hat.se <- 1 / sqrt(out$hessian)
  return(list(lambda.hat = lambda.hat, lambda.hat.se = lambda.hat.se))
}

# END