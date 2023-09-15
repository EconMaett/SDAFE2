"DOC_garch" <- function(E, L = 4., c = 2.25, theta.ini = NULL,
                        n.ahead = 10, common.tail.index = FALSE) {
  E <- as.matrix(E)
  d <- ncol(E)
  n <- nrow(E)
  p <- d * (d - 1) / 2
  L <- L
  c <- c
  
  U <- t(matrix.sqrt.inv(var(E)))
  U.inv <- t(matrix.sqrt(var(E)))
  Z <- E %*% U
  
  if (is.null(theta.ini)) {
    theta.ini <- as.matrix(rep(0, p))
  } else {
    theta.ini <- as.matrix(theta.ini)
  }
  
  # Perform DOC in Volatility Estimation
  out <- optim(
    par = theta.ini, 
    fn = DOC.obj, 
    gr = NULL,
    Z = Z, 
    c = c, 
    L = L, 
    method = "L-BFGS-B",
    lower = -pi, 
    upper = pi,
    control = list(), 
    hessian = FALSE
  )
  W <- canonicalW(theta2W(out$par))
  W.inv <- t(W)
  theta.hat <- W2theta(W)
  S <- Z %*% t(W)
  M.inv <- U %*% t(W)
  M <- W %*% U.inv
  
  # Esitmate sigma squared for DOCs and PCs:
  V.doc <- array(numeric(n * d * d), c(d, d, n))
  V.pca <- array(numeric(n * d * d), c(d, d, n))
  H.doc <- array(numeric(n * d * d), c(d, d, n))
  H.pca <- array(numeric(n * d * d), c(d, d, n))
  V.doc.forecast <- array(numeric(n.ahead * d * d), c(d, d, n.ahead))
  V.pca.forecast <- array(numeric(n.ahead * d * d), c(d, d, n.ahead))
  H.doc.forecast <- array(numeric(n.ahead * d * d), c(d, d, n.ahead))
  H.pca.forecast <- array(numeric(n.ahead * d * d), c(d, d, n.ahead))
  
  DOC.resid <- matrix(0, n, d)
  PC.resid <- matrix(0, n, d)
  
  DOC.est <- matrix(0, d, 4, dimnames = list(NULL, c("omega", "alpha1", "beta1", "shape")))
  PC.est <- matrix(0, d, 3, dimnames = list(
    NULL,
    c("omega", "alpha1", "beta1")
  ))
  
  if (common.tail.index == TRUE) {
    tail.index.grid <- seq(4.1, 10, 0.1)
    n.grid <- length(tail.index.grid)
    tail.index.llh <- numeric(n.grid)
    for (g in 1:n.grid) {
      for (i in 1:d) {
        fitDOC <- garchFit(~ garch(1, 1),
                           data = S[, i],
                           include.mean = FALSE,
                           cond.dist = "std", shape = tail.index.grid[g], include.shape = F,
                           trace = FALSE
        ) # tail.index.grid[g], trace = F)
        tail.index.llh[g] <- tail.index.llh[g] - fitDOC@fit$llh
      }
    }
    
    est.common.tail.index <- tail.index.grid[which.max(tail.index.llh)]
    
    for (i in 1:d) {
      fitDOC <- garchFit(~ garch(1, 1),
                         data = S[, i],
                         include.mean = FALSE,
                         cond.dist = "std", shape = est.common.tail.index,
                         include.shape = F, trace = F
      )
      
      DOC.resid[, i] <- fitDOC@residuals / fitDOC@sigma.t
      V.doc[i, i, ] <- fitDOC@h.t
      V.doc.forecast[i, i, ] <- predict(fitDOC, n.ahead = n.ahead)$standardDeviation^2
      DOC.est[i, ] <- c(fitDOC@fit$par, est.common.tail.index)
      
      fitPC <- garchFit(~ garch(1, 1),
                        data = Z[, i],
                        include.mean = FALSE,
                        cond.dist = "norm", trace = F
      )
      PC.resid[, i] <- fitPC@residuals / fitPC@sigma.t
      V.pca[i, i, ] <- fitPC@h.t
      V.pca.forecast[i, i, ] <- predict(fitPC, n.ahead = n.ahead)$standardDeviation^2
      PC.est[i, ] <- fitPC@fit$par
    }
  }
  
  
  if (common.tail.index == FALSE) {
    for (i in 1:d) {
      fitDOC <- garchFit(~ garch(1, 1),
                         data = S[, i],
                         include.mean = FALSE,
                         cond.dist = "std",
                         include.shape = TRUE, trace = F
      )
      
      DOC.resid[, i] <- fitDOC@residuals / fitDOC@sigma.t
      V.doc[i, i, ] <- fitDOC@h.t
      V.doc.forecast[i, i, ] <- predict(fitDOC, n.ahead = n.ahead)$standardDeviation^2
      DOC.est[i, ] <- fitDOC@fit$par
      
      fitPC <- garchFit(~ garch(1, 1),
                        data = Z[, i],
                        include.mean = FALSE,
                        cond.dist = "norm", trace = F
      )
      PC.resid[, i] <- fitPC@residuals / fitPC@sigma.t
      V.pca[i, i, ] <- fitPC@h.t
      V.pca.forecast[i, i, ] <- predict(fitPC, n.ahead = n.ahead)$standardDeviation^2
      PC.est[i, ] <- fitPC@fit$par
    }
  }
  
  tM <- t(M)
  tU.inv <- t(U.inv)
  for (t in 1:n) {
    H.doc[, , t] <- tM %*% V.doc[, , t] %*% M
    H.pca[, , t] <- tU.inv %*% V.pca[, , t] %*% U.inv
  }
  for (t in 1:n.ahead) {
    H.doc.forecast[, , t] <- tM %*% V.doc.forecast[, , t] %*% M
    H.pca.forecast[, , t] <- tU.inv %*% V.pca.forecast[, , t] %*% U.inv
  }
  # matrices W, U, M returned for dx1 (column) components s_t, z_t
  return(
    list(
      Sigma.doc = H.doc, Sigma.pca = H.pca, Z.hat = Z, S.hat = S,
      coef.doc = DOC.est, coef.pca = PC.est, theta.hat = theta.hat,
      W.hat = W, U.hat = t(U), M.hat = tM, DOC.resid = DOC.resid, PC.resid = PC.resid,
      Sigma.doc.forecast = H.doc.forecast, Sigma.pca.forecast = H.pca.forecast
    )
  )
}

# END